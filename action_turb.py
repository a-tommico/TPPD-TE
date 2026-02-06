from pyfluids import Fluid, FluidsList, Input
import matplotlib.pyplot as plt
import numpy as np

# --- 1. System Setup ---
fluid_name = FluidsList.MDM
fluid = Fluid(fluid_name)
mass_flow_rate = 2.0  # kg/s

# --- 2. Cycle Calculation ---

# State 1 (Pump Inlet): Saturated liquid at T = 95°C
T_cond = 95
state1 = fluid.with_state(Input.temperature(T_cond), Input.quality(0))
p_cond = state1.pressure

# State 2 (Pump Outlet): Compressed to P = 10 bar
p_high = 10e5  # 10 bar
pump_eff = 0.70

# Isentropic step
state2_isen = state1.isentropic_compression_to_pressure(p_high)
h2_isen = state2_isen.enthalpy
# Actual step
w_pump_ideal = h2_isen - state1.enthalpy
w_pump_actual = w_pump_ideal / pump_eff
h2 = state1.enthalpy + w_pump_actual
state2 = fluid.with_state(Input.pressure(p_high), Input.enthalpy(h2))

# State 3 (Turbine Inlet): Saturated Vapor at P = 10 bar
state3 = fluid.with_state(Input.pressure(p_high), Input.quality(100))

# State 4 (Turbine Outlet): Expansion to P_cond
turb_eff = 0.85
state4_isen = state3.isentropic_expansion_to_pressure(p_cond)
h4_isen = state4_isen.enthalpy
w_turb_ideal = state3.enthalpy - h4_isen
w_turb_actual = w_turb_ideal * turb_eff
h4 = state3.enthalpy - w_turb_actual
state4 = fluid.with_state(Input.pressure(p_cond), Input.enthalpy(h4))

# --- 3. Regenerator Logic ---
effectiveness = 0.80
cp_hot = state4.specific_heat
cp_cold = state2.specific_heat

C_hot = mass_flow_rate * cp_hot
C_cold = mass_flow_rate * cp_cold
C_min = min(C_hot, C_cold)

Q_max = C_min * (state4.temperature - state2.temperature)
Q_regen = effectiveness * Q_max

h5 = state4.enthalpy - (Q_regen / mass_flow_rate)
state5 = fluid.with_state(Input.pressure(p_cond), Input.enthalpy(h5))

h2_prime = state2.enthalpy + (Q_regen / mass_flow_rate)
state2_prime = fluid.with_state(Input.pressure(p_high), Input.enthalpy(h2_prime))

# --- 4. Performance Metrics ---
W_turb = mass_flow_rate * (state3.enthalpy - state4.enthalpy)
W_pump = mass_flow_rate * (state2.enthalpy - state1.enthalpy)
W_net = W_turb - W_pump
Q_in = mass_flow_rate * (state3.enthalpy - state2_prime.enthalpy)
Q_out = mass_flow_rate * (state5.enthalpy - state1.enthalpy)
efficiency = (W_net / Q_in) * 100 if Q_in > 0 else 0

# --- 5. Output Data ---
print(f"--- Cycle Results (MDM) ---")
print(f"{'State':<15} {'T (°C)':<10} {'P (bar)':<10} {'h (kJ/kg)':<10} {'s (kJ/kgK)':<10}")
print("-" * 60)

states = {'1 (Pump In)': state1, '2 (Pump Out)': state2, "2' (Regen Out)": state2_prime, 
          '3 (Turb In)': state3, '4 (Turb Out)': state4, '5 (Cond In)': state5}

for name, st in states.items():
    print(f"{name:<15} {st.temperature:<10.2f} {st.pressure/1e5:<10.2f} {st.enthalpy/1e3:<10.2f} {st.entropy/1e3:<10.3f}")

print("-" * 60)
print(f"Net Power (per kg/s): {W_net/1000:.2f} kW")
print(f"Thermal Efficiency:   {efficiency:.2f} %")
print(f"Condenser Heat:       {Q_out/1000:.2f} kW")
print(f"Input thermal power:  {Q_in/1000:.2f} kW")

# ------------------------------------
# DIMENSIONAMENTO TURBINA
# ------------------------------------

# Dati di input
delta_h_turb = state3.enthalpy - state4.enthalpy  # J/kg
N_stages = 4 # numero di stadi di espansione
beta = p_high / p_cond  # rapporto di pressione totale
alfa = 15 # angolo di uscita dallo statore in gradi (flow deflection angle)

print(f"\n--- Dimensionamento Turbina ---")
print(f"Salto entalpico totale turbina: {delta_h_turb/1000:.2f} kJ/kg")
print(f"Rapporto di espansione totale turbina: {beta:.2f}")

# TRIANGOLI DI VELOCITÀ
# 0: ingresso statore, 1: uscita statore / ingresso rotore, 2: uscita rotore
delta_h = delta_h_turb / N_stages  # J/kg salto entalpico per ogni stadio
beta_s = beta ** (1 / N_stages)  # rapporto di pressione per stadio
phi = 0.97 # speed reduction coefficient at the stator (phi=c1/c1_is), valore tabulato per alfa=30° e beta_st=4
psi = phi # speed reduction coefficient at the rotor (psi=w2/w2_is), assunto uguale a phi (symmetrical design)

print(f"Salto entalpico per stadio: {delta_h/1000:.2f} kJ/kg")
print(f"Rapporto di pressione per stadio: {beta_s:.2f}")

# Triangoli di velocità allo statore

c0 = 0 # velocità di ingresso (assunta nulla, dato che c_in << c_out)
c1_is = np.sqrt(2 * delta_h + c0**2)  # velocità ideale in uscita dallo statore
c1 = phi * c1_is  # velocità reale in uscita dallo statore
c_s = state3.sound_speed  # velocità del suono (Hp: gas ideale)

print(f"\n--- Triangoli di Velocità allo Statore ---")

print(f"Velocità ideale in uscita dallo statore (c1_is): {c1_is:.2f} m/s")
print(f"Velocità reale in uscita dallo statore (c1): {c1:.2f} m/s")
print(f"Velocità del suono alla temperatura di ingresso (c_s): {c_s:.2f} m/s")
print(f"Numero di Mach alla uscita dallo statore (Ma): {c1/c_s:.2f}")

u_opt = c1 * np.cos(np.radians(alfa)) / 2  # velocità di rotazione ottimale

c1_t = u_opt # componente tangenziale di c1
c1_a = c1 * np.sin(np.radians(alfa))  # componente assiale di c1
c1_t = c1 * np.cos(np.radians(alfa))  # componente tangenziale di c1

w1_a = c1_a  # componente assiale di w1
w1_t = c1_t - u_opt  # componente tangenziale di w1
w1 = np.sqrt(w1_t**2 + w1_a**2)  # velocità relativa in ingresso al rotore

beta_1 = np.arctan(w1_a / w1_t)  # angolo di ingresso al rotore

print(f"Velocità di rotazione ottimale (u_opt): {u_opt:.2f} m/s")
print(f"Velocità relativa in ingresso al rotore (w1): {w1:.2f} m/s")
print(f"Angolo di uscita dallo statore (alfa): {alfa:.2f}°")
print(f"Angolo di ingresso al rotore (beta_1): {np.degrees(beta_1):.2f}°")

# Triangoli di velocità al rotore

w2_id = w1
beta_2 = np.pi - beta_1  # angolo di uscita dal rotore (assunto simmetrico all'ingresso)

w2 = psi * w2_id  # velocità relativa reale in uscita dal rotore
w2_a = w2 * np.sin(beta_2)  # componente assiale di w2
w2_t = w2 * np.cos(beta_2)  # componente tangenziale di w2

c2_a = w2_a  # componente assiale di c2
c2_t = u_opt + w2_t  # componente tangenziale di c2
c2 = np.sqrt(c2_t**2 + c2_a**2)  # velocità assoluta in uscita dal rotore
alfa_2 = np.arctan(c2_a / c2_t)  # angolo di uscita dal rotore rispetto alla velocità assoluta

print(f"\n--- Triangoli di Velocità al Rotore ---")

print(f"Velocità relativa reale in uscita dal rotore (w2): {w2:.2f} m/s")
print(f"Angolo di uscita dal rotore (beta_2): {np.degrees(beta_2):.2f}°")
print(f"Velocità assoluta in uscita dal rotore (c2): {c2:.2f} m/s")
print(f"Componente tangenziale di c2 (c2_t): {c2_t:.2f} m/s")
print(f"Componente assiale di c2 (c2_a): {c2_a:.2f} m/s")
print(f"Angolo di uscita dal rotore rispetto alla velocità assoluta (alfa_2): {np.degrees(alfa_2):.2f}°")

# DIMENSIONAMENTO GEOMETRICO
p = 1 # coppie polari generatore
f = 50 # frequenza rete elettrica (Hz)
K = 2 # rapporto di trasmissione
n_rot = K * (f * 60) / p  # velocità di rotazione (giri/min)
omega_rot = n_rot * (2 * np.pi) / 60  # velocità angolare (rad/s)

p_3stator = state3.pressure / beta_s  # pressione dopo statore
h_3stator_is = state3.enthalpy - delta_h  # entalpia ideale dopo statore

state3_stator = fluid.with_state(Input.pressure(p_3stator), Input.enthalpy(h_3stator_is))

D_m = 2 * u_opt / omega_rot  # diametro medio turbina
rho = state3_stator.density  # densità fluido alla temperatura di ingresso (kg/m³) A INGRESSO DEL ROTORE
parz = 0.15 # grado di parzializzazione
b = mass_flow_rate / (rho * 0.95 * c1_a * np.pi * D_m * parz)  # altezza del canale di flusso (m), con coefficiente di riempimento 0.95

print(f"\n--- Dimensionamento Geometrico ---")
print(f"Velocità di rotazione (n_rot): {n_rot:.2f} giri/min")
print(f"Pressione dopo statore (p_3stator): {p_3stator/1e5:.2f} bar")
print(f"Densità fluido a ingresso rotore (rho): {rho:.2f} kg/m³")
print(f"Diametro medio turbina (D_m): {D_m:.2f} m")
print(f"Altezza del canale di flusso (b): {b:.10f} m")

# RAPPRESENTAZIONE GRAFICA TRIANGOLI DI VELOCITÀ

# Parametri di visualizzazione
scale = 1.0
y0, y01, y1, y2 =-55, 45, 50, 150  # Offset verticale per separare i due stadi
offset = 5       # Distanza del testo dalla punta/centro

plt.figure(figsize=(10, 8))

# --- STADIO 0 (Ingresso statore) ---
plt.arrow(-100, y0, 0, 2, head_width=2, head_length=4, fc='blue', ec='blue')
plt.text(-100, y0 + 2 + offset, r"$\vec{c}_0 \approx 0$", color='blue', fontsize=12, ha='center')

# --- STADIO 1 (Ingresso rotore/uscita statore) ---
# c1
plt.arrow(0, y1, c1_t*scale, c1_a*scale, head_width=2, head_length=4, fc='blue', ec='blue')
plt.text(c1_t*scale/2, y1 + c1_a*scale/2 + offset, r"$\vec{c}_1$", color='blue', fontsize=12, ha='center')

# w1 (parte dalla fine di u se u parte dall'origine, o viceversa. 
# Qui lo facciamo partire dall'origine per confrontare gli angoli)
plt.arrow(0, y1, w1_t*scale, w1_a*scale, head_width=2, head_length=4, fc='green', ec='green')
plt.text(w1_t*scale/2 - offset, y1 + w1_a*scale/2, r"$\vec{w}_1$", color='green', fontsize=12, ha='right')

# u (congiunge le punte o chiude il triangolo)
plt.arrow(w1_t*scale, y1 + w1_a*scale, (c1_t-w1_t)*scale, 0, head_width=2, head_length=4, fc='red', ec='red')
plt.text((w1_t + c1_t)*scale/2, y1 + w1_a*scale + offset, r"$\vec{u}$", color='red', fontsize=12, ha='center')


# --- STADIO 2 (Uscita rotore) ---
# c2
plt.arrow(0, y2, c2_t*scale, c2_a*scale, head_width=2, head_length=4, fc='blue', ec='blue')
plt.text(c2_t*scale/2 + offset, y2 + c2_a*scale/2, r"$\vec{c}_2$", color='blue', fontsize=12)

# w2
plt.arrow(0, y2, w2_t*scale, w2_a*scale, head_width=2, head_length=4, fc='green', ec='green')
plt.text(w2_t*scale/2 - offset, y2 + w2_a*scale/2, r"$\vec{w}_2$", color='green', fontsize=12, ha='right')

# u
plt.arrow(w2_t*scale, y2 + w2_a*scale, (c2_t-w2_t)*scale, 0, head_width=2, head_length=4, fc='red', ec='red')
plt.text((w2_t + c2_t)*scale/2, y2 + w2_a*scale - offset*2, r"$\vec{u}$", color='red', fontsize=12, ha='center')

# Linea Ingresso Rotore (corrispondente a y_in)
plt.hlines(y=y1, xmin=-150, xmax=100, colors='gray', linestyles='--', linewidth=1, alpha=0.7)
plt.text(100, y1, " Ingresso Rotore", va='center', fontsize=10, color='gray', fontstyle='italic')

# Linea Uscita Rotore (corrispondente a y_out)
plt.hlines(y=y2, xmin=-150, xmax=100, colors='gray', linestyles='--', linewidth=1, alpha=0.7)
plt.text(100, y2, " Uscita Rotore", va='center', fontsize=10, color='gray', fontstyle='italic')

#Linea Ingresso Statore
plt.hlines(y=y0, xmin=-150, xmax=100, colors='gray', linestyles='--', linewidth=1, alpha=0.7)
plt.text(100, y0, " Ingresso Statore", va='center', fontsize=10, color='gray', fontstyle='italic')

#Linea Uscita Statore
plt.hlines(y=y01, xmin=-150, xmax=100, colors='gray', linestyles='--', linewidth=1, alpha=0.7)
plt.text(100, y01, " Uscita Statore", va='center', fontsize=10, color='gray', fontstyle='italic')

# -------------------------------------------------
# LINEA MEDIA DELLA PALA DEL ROTORE (Bezier cubica)
# -------------------------------------------------

# Lunghezza "virtuale" della pala nel piano t–a
L = 0.8 * (y2 - y1)

# Versori delle velocità relative
e_w1 = np.array([w1_t, w1_a]) / np.sqrt(w1_t**2 + w1_a**2)
e_w2 = np.array([w2_t, w2_a]) / np.sqrt(w2_t**2 + w2_a**2)

# Punti di controllo Bezier
P0 = np.array([0.0, y1])                 # bordo d'attacco
P3 = np.array([0.0, y2])                 # bordo d'uscita

P1 = P0 + L * e_w1                       # direzione w1
P2 = P3 - L * e_w2                       # direzione w2 (tangente in uscita)

# Parametro curva
t = np.linspace(0, 1, 200)

# Bezier cubica
B = ((1 - t)**3)[:, None] * P0 \
  + (3 * (1 - t)**2 * t)[:, None] * P1 \
  + (3 * (1 - t) * t**2)[:, None] * P2 \
  + (t**3)[:, None] * P3

# Plot linea media
plt.plot(B[:, 0], B[:, 1], 'k', linewidth=1, linestyle='--', label="Linea media rotore")

# DA SISTEMARE PLOT LINEA STATORE

# -------------------------------------------------
# LINEA MEDIA DELLA PALA DELO STATORE (Bezier cubica)
# -------------------------------------------------

# Lunghezza "virtuale" della pala nel piano t–a
LL = 0.8 * (y01 - y0)

# Versori delle velocità relative
e_c0 = np.array([0, 2]) / np.sqrt(2**2)
e_w1 = np.array([w1_t, w1_a]) / np.sqrt(w1_t**2 + w1_a**2)

# Punti di controllo Bezier
P0 = np.array([-100, y0])                 # bordo d'attacco
P3 = np.array([0.0, y01])                 # bordo d'uscita

P1 = P0 + LL * e_c0                       # direzione c0
P2 = P3 - LL * e_w1                       # direzione w1 (tangente in uscita)

# Parametro curva
tt = np.linspace(0, 1, 200)

# Bezier cubica
BB = ((1 - tt)**3)[:, None] * P0 \
  + (3 * (1 - tt)**2 * tt)[:, None] * P1 \
  + (3 * (1 - tt) * tt**2)[:, None] * P2 \
  + (tt**3)[:, None] * P3

# Plot linea media
plt.plot(BB[:, 0], BB[:, 1], 'k', linewidth=1, linestyle='--', label="Linea media statore")

# Evidenzia bordo d'attacco e uscita
# plt.scatter([PP1[0], PP0[0]], [PP1[1], PP0[1]], color='k', zorder=2)

plt.legend()

# Formattazione grafico
plt.axis('equal')
plt.title("Triangoli di Velocità - ORC Turbine Stage")
plt.xlabel("Direzione tangenziale →")
plt.ylabel("Direzione assiale →")
plt.xticks([])
plt.yticks([])
plt.show()