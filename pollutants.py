import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# ------------------------
# DATI INPUT
# ------------------------
# Formula elementare del combustibile secco (base molare equivalente)
# C_a H_b O_c N_d
aC = 1.0
bH = 1.67
cO = 0.51
dN = 0.048

# Umidità del combustibile (wet basis): frazione massica di H2O nel combustibile tal quale
moisture = 0.112  # 11.2%

# Portate annue (kg)
mass_secco = 800
mass_umido = 259
mass_fanghi = 70
mass_tot = mass_secco + mass_umido + mass_fanghi

# LHV su base secca (J/kg_dry)
LHV_dry = 19.04e6

# Range equivalence ratio
phi_list = np.linspace(0.5, 1.2, 100)

# Condizioni ambiente
Tamb = 298.15
Pamb = ct.one_atm

# Portata di riferimento (serve solo per scalare poi a kg/anno; i risultati specifici non dipendono da questo valore)
m_b_wet = 2.0  # kg/h (tal quale)

# ------------------------
# PROPRIETÀ ELEMENTARI E CONVERSIONI DI BASE
# ------------------------
AW = {"C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007}

# Massa molare del “pseudo-combustibile” (kg/kmol, numericamente uguale a g/mol)
MW_fuel = aC * AW["C"] + bH * AW["H"] + cO * AW["O"] + dN * AW["N"]

# Frazioni massiche elementari nel secco
wC_dry = (aC * AW["C"]) / MW_fuel
wH_dry = (bH * AW["H"]) / MW_fuel
wO_dry = (cO * AW["O"]) / MW_fuel
wN_dry = (dN * AW["N"]) / MW_fuel

# Portate (kg/h): secco e acqua
m_b_dry = m_b_wet * (1.0 - moisture)  # kg_dry/h
m_H2O = m_b_wet * moisture  # kg/h

# Portate elementari nel secco (kg/h)
m_C = wC_dry * m_b_dry
m_H = wH_dry * m_b_dry
m_O = wO_dry * m_b_dry
m_N = wN_dry * m_b_dry

# ------------------------
# STECHIOMETRIA DELL’OSSIGENO E COMPOSIZIONE ARIA
# ------------------------
# O2 stechiometrico (moli O2 per mole di combustibile secco):
# O2,st = a + b/4 - c/2   per C_a H_b O_c
O2_st_mol_per_mol_fuel = aC + bH / 4.0 - cO / 2.0

# Conversione in rapporto massico O2/kg_dry (kg_O2 per kg_dry)
O2_per_kg_dry = (O2_st_mol_per_mol_fuel * 32.0) / MW_fuel

# Portata O2 stechiometrica (kg/h)
m_O2_st = O2_per_kg_dry * m_b_dry

# Aria: N2/O2 = 3.76 in moli; convertiamo in rapporto di massa per lavorare in kg/h
N2_per_O2_molar = 3.76
N2_per_O2_mass = (N2_per_O2_molar * 28.0134) / 31.998  # kg_N2/kg_O2

# ------------------------
# Cantera mechanism (deve includere specie atomiche C,H,O,N e specie principali di combustione)
# ------------------------
MECH = "gri30_gasifier.yaml"


# ------------------------
# FUNZIONE: LHV residua nei prodotti (J/kg miscela)
# ------------------------
def LHV_mass(g: ct.Solution) -> float:
    """
    Stima della LHV residua del gas: quanta energia chimica rimane nei prodotti
    (per kg di miscela), calcolata aggiungendo O2 stechiometrico per ossidare
    completamente la miscela a CO2 e H2O a condizioni di riferimento.

    Serve per gestire correttamente le miscele ricche: se nei prodotti restano CO/H2,
    quella è energia non convertita in calore e quindi riduce la temperatura adiabatica.
    """
    Tin, Pin = g.T, g.P
    X_in = g.X

    Tref = 300.0
    Pref = ct.one_atm

    X = X_in.copy()
    X[X < 1e-12] = 0.0
    s = X.sum()
    if s <= 0:
        return 0.0
    X /= s

    g.TPX = Tref, Pref, X

    iO2 = g.species_index("O2")
    iH2O = g.species_index("H2O")
    iCO2 = g.species_index("CO2")

    # Conta atomi “combustibili” in specie non già completamente ossidate
    N_H = N_C = 0.0
    for k in range(g.n_species):
        if X[k] == 0.0:
            continue
        # specie già ossidate non contribuiscono
        if k in (iH2O, iCO2):
            continue
        for e in range(g.n_elements):
            n = g.n_atoms(k, e)
            if n == 0:
                continue
            en = g.element_name(e)
            if en == "H":
                N_H += X[k] * n
            elif en == "C":
                N_C += X[k] * n

    if (N_H + N_C) == 0.0:
        g.TPX = Tin, Pin, X_in
        return 0.0

    # O atoms needed: H -> H2O needs 1 O per 2 H; C -> CO2 needs 2 O per C
    O_atoms_needed = (N_H / 2.0) + (2.0 * N_C)
    O2_needed = O_atoms_needed / 2.0

    Xmix = X.copy()
    Xmix[iO2] += O2_needed
    Xmix /= Xmix.sum()

    # Enthalpy of reactants (with added O2) at reference conditions
    g.TPX = Tref, Pref, Xmix
    h_react = g.enthalpy_mass

    # Enthalpy of equilibrium products at reference conditions
    g.equilibrate("TP")
    h_prod = g.enthalpy_mass

    # Restore original state
    g.TPX = Tin, Pin, X_in

    return max(h_react - h_prod, 0.0)


# ------------------------
# FUNZIONE: calcolo Tad con iterazione TP -> LHV residua -> HP
# ------------------------
def compute_Tad(phi: float, debug: bool = False) -> ct.Solution:
    """
    Calcolo della temperatura adiabatica in equilibrio a pressione costante.

    Procedura:
    1) si costruiscono i reagenti (elementi + aria + umidità) su base massica
    2) equilibrio TP a Tamb per ottenere una miscela molecolare consistente
    3) si impone un target di entalpia (HP) aggiungendo:
       DH = LHV_dry * (m_b_dry / m_f)
       e sottraendo l’energia residua nei prodotti (LHV_mass)
    4) si itera su una temperatura T1 finché la soluzione si stabilizza
    """
    gas = ct.Solution(MECH)
    idx = gas.species_index

    # Indici specie
    iO2 = idx("O2")
    iN2 = idx("N2")
    iH2O = idx("H2O")
    iC = idx("C")
    iH = idx("H")
    iO = idx("O")
    iN = idx("N")

    lam = 1.0 / phi  # eccesso d’aria (lambda)

    # Portate aria (kg/h)
    m_O2 = m_O2_st * lam
    m_N2 = m_O2 * N2_per_O2_mass
    m_air = m_O2 + m_N2

    # Massa totale reagenti (kg/h)
    m_f = m_b_wet + m_air

    # Vettore portate specie (kg/h) -> frazioni massiche
    ms = np.zeros(gas.n_species)
    ms[iO2] = m_O2
    ms[iN2] = m_N2
    ms[iH2O] = m_H2O

    # Combustibile come elementi (kg/h)
    ms[iC] = m_C
    ms[iH] = m_H
    ms[iO] = m_O
    ms[iN] = m_N

    Y_in = ms / ms.sum()

    # Passo 1: equilibrio TP a Tamb per ottenere specie molecolari
    gas.TPY = Tamb, Pamb, Y_in
    gas.equilibrate("TP")
    Y_eq = gas.Y.copy()

    # Energia chimica disponibile per kg di miscela (base secca)
    DH = LHV_dry * (m_b_dry / m_f)  # J/kg_miscela

    # Iterazione su T
    T1 = 600.0
    tol = 1.0
    it = 0

    while tol > 1e-2 and it < 60:
        it += 1

        gas.TPY = T1, Pamb, Y_eq
        gas.equilibrate("TP")

        # Energia residua nei prodotti (per miscele ricche diventa significativa)
        LHV_g = LHV_mass(gas)

        # Chiusura adiabatica (HP) con entalpia target
        gas.TP = Tamb, Pamb
        h_amb = gas.enthalpy_mass
        h_target = h_amb + (DH - LHV_g)

        gas.HP = h_target, Pamb
        gas.equilibrate("HP")

        Tnew = gas.T
        tol = abs(Tnew - T1)

        if debug:
            released_frac = 0.0 if DH <= 0 else max(0.0, 1.0 - LHV_g / DH)
            print("\n" + "=" * 70)
            print(f"phi = {phi:.3f} | lambda = {lam:.3f}")
            print(f"m_b_wet = {m_b_wet:.3f} kg/h, m_b_dry = {m_b_dry:.3f} kg_dry/h")
            print(
                f"m_O2 = {m_O2:.3f} kg/h, m_air = {m_air:.3f} kg/h, m_f = {m_f:.3f} kg/h"
            )
            print(f"DH (energia chimica) = {DH/1e6:.3f} MJ/kg_miscela")
            print(
                f"iter {it:02d}: T1={T1:8.2f} K -> Tnew={Tnew:8.2f} K | "
                f"LHV_g={LHV_g/1e6:6.3f} MJ/kg | frazione_rilasciata≈{released_frac:5.3f} | err={tol:7.4f}"
            )

        T1 = Tnew

    return gas


# ------------------------
# CICLO SUL PHI
# ------------------------
T_flame = []
NO_ppm = []
CO_ppm = []
O2_ppm = []
CO2_ppm = []

for phi in phi_list:
    gas = compute_Tad(phi, debug=False)

    T_flame.append(gas.T)
    NO_ppm.append(gas["NO"].X[0] * 1e6)
    CO_ppm.append(gas["CO"].X[0] * 1e6)
    O2_ppm.append(gas["O2"].X[0] * 1e6)
    CO2_ppm.append(gas["CO2"].X[0] * 1e6)

# ------------------------
# GRAFICI
# ------------------------
plt.figure(figsize=(10, 6))
plt.plot(phi_list, T_flame, "r-", label="T adiabatica [K]")
plt.xlabel("Equivalence ratio φ")
plt.ylabel("Temperatura [K]")
plt.grid(True)
plt.title("Temperatura adiabatica (equilibrio) vs equivalence ratio")
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(phi_list, NO_ppm, label="NO [ppm]")
plt.plot(phi_list, CO_ppm, label="CO [ppm]")
plt.plot(phi_list, O2_ppm, label="O2 [ppm]")
plt.plot(phi_list, CO2_ppm, label="CO2 [ppm]")
plt.xlabel("Equivalence ratio φ")
plt.ylabel("Concentrazione [ppm]")
plt.grid(True)
plt.title("Andamento dei prodotti di combustione (equilibrio) vs equivalence ratio")
plt.legend()
plt.show()

# ------------------------
# TROVA PHI OTTIMALE (minimizza NO + CO)
# ------------------------
pollutants = np.array(NO_ppm) + np.array(CO_ppm)
idx_opt = int(np.argmin(pollutants))
phi_opt = float(phi_list[idx_opt])
print(f"Equivalence ratio ottimale (min NO+CO): φ = {phi_opt:.2f}")

# ------------------------
# SIMULAZIONE CON φ OTTIMALE (con stampe di controllo)
# ------------------------
gas = compute_Tad(phi_opt, debug=True)

# ------------------------
# CONVERSIONE IN kg/anno e g/kWh
# ------------------------
# Ore/anno coerenti con massa secca annua e portata secca (kg_dry/h)
hours_per_year = mass_secco / m_b_dry

# Portata totale prodotti (kg/h) ~ portata totale reagenti (kg/h)
lam_opt = 1.0 / phi_opt
m_O2_opt = m_O2_st * lam_opt
m_N2_opt = m_O2_opt * N2_per_O2_mass
m_air_opt = m_O2_opt + m_N2_opt
m_f_opt = m_b_wet + m_air_opt

# kg/h delle specie = Y * m_f
m_NO_h = gas["NO"].Y[0] * m_f_opt
m_CO_h = gas["CO"].Y[0] * m_f_opt
m_CO2_h = gas["CO2"].Y[0] * m_f_opt

# kg/anno
m_NO_anno = m_NO_h * hours_per_year
m_CO_anno = m_CO_h * hours_per_year
m_CO2_anno = m_CO2_h * hours_per_year

# g/kWh (energia termica su base secca annua)
E_tot_J = mass_secco * LHV_dry
E_tot_kWh = E_tot_J / 3.6e6

g_NO_kWh = (m_NO_anno * 1e3) / E_tot_kWh
g_CO_kWh = (m_CO_anno * 1e3) / E_tot_kWh
g_CO2_kWh = (m_CO2_anno * 1e3) / E_tot_kWh

print("\n--- RISULTATI ANNUALI ---")
print(f"NO:  {m_NO_anno:.2f} kg/anno, {g_NO_kWh:.2f} g/kWh")
print(f"CO:  {m_CO_anno:.2f} kg/anno, {g_CO_kWh:.2f} g/kWh")
print(f"CO2: {m_CO2_anno:.2f} kg/anno, {g_CO2_kWh:.2f} g/kWh")
print(f"Temperatura adiabatica φ ottimale: {gas.T:.1f} K")
