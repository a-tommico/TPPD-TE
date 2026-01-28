from pyfluids import Fluid, FluidsList, Input
import matplotlib.pyplot as plt
import numpy as np

# --- 1. System Setup ---
fluid_name = FluidsList.MDM
fluid = Fluid(fluid_name)
mass_flow_rate = 1.0  # kg/s

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

print(f"Velocità ideale in uscita dallo statore (c1_is): {c1_is:.2f} m/s")
print(f"Velocità reale in uscita dallo statore (c1): {c1:.2f} m/s")
print(f"Velocità del suono alla temperatura di ingresso (c_s): {c_s:.2f} m/s")

u_opt = c1 * np.cos(np.radians(alfa))  # velocità di rotazione ottimale
c1_t = u_opt # componente tangenziale di c1
w1 = c1 * np.sin(np.radians(alfa))  # velocità relativa in ingresso al rotore
c1_a = w1  # componente assiale di c1
beta_1 = np.arctan(w1 / (c1_t - u_opt))  # angolo di ingresso al rotore

print(f"Velocità di rotazione ottimale (u_opt): {u_opt:.2f} m/s")
print(f"Angolo di uscita dallo statore (alfa): {alfa:.2f}°")
print(f"Angolo di ingresso al rotore (beta_1): {np.degrees(beta_1):.2f}°")

# Triangoli di velocità al rotore

w2_is = c1 # velocità relativa ideale in uscita dal rotore
w2 = psi * w2_is  # velocità relativa reale in uscita dal rotore
beta_2 = 180 - alfa # angolo di uscita dal rotore (assunto uguale ad alfa per simmetria)

w2_t = w2 * np.cos(np.radians(beta_2))  # componente tangenziale di w2
w2_a = w2 * np.sin(np.radians(beta_2))  # componente assiale di w2

print(f"Velocità relativa reale in uscita dal rotore (w2): {w2:.2f} m/s")
print(f"Angolo di uscita dal rotore (beta_2): {beta_2:.2f}°")
print(f"Componente tangenziale di w2 (w2_t): {np.abs(w2_t):.2f} m/s")
print(f"Componente assiale di w2 (w2_a): {w2_a:.2f} m/s")

c2_t = u_opt + w2_t  # componente tangenziale di c2
alfa_2 = np.arctan(w2_a / c2_t)  # angolo di uscita dal rotore rispetto alla velocità assoluta
print(f"Componente tangenziale di c2 (c2_t): {c2_t:.2f} m/s")
print(f"Angolo di uscita dal rotore rispetto alla velocità assoluta (alfa_2): {np.degrees(alfa_2):.2f}°")

# DIMENSIONAMENTO GEOMETRICO
p = 1 # coppie polari generatore
f = 50 # frequenza rete elettrica (Hz)
n_rot = (f * 60) / p  # velocità di rotazione (giri/min)
omega_rot = n_rot * (2 * np.pi) / 60  # velocità angolare (rad/s)

D_m = 2 * u_opt / omega_rot  # diametro medio turbina
rho = state3.density  # densità fluido alla temperatura di ingresso (kg/m³)
b = mass_flow_rate / (rho * 0.95 * c1_a * np.pi * D_m)  # altezza del canale di flusso (m), con coefficiente di riempimento 0.95

print(f"\n--- Dimensionamento Geometrico ---")
print(f"Velocità di rotazione (n_rot): {n_rot:.2f} giri/min")
print(f"Diametro medio turbina (D_m): {D_m:.2f} m")
print(f"Altezza del canale di flusso (b): {b:.10f} m")

plt.figure(figsize=(8, 10))

# Offset assiale delle sezioni
y0 = 0      # ingresso statore
y1 = 60     # uscita statore / ingresso rotore
y2 = 120    # uscita rotore

scale = 1.0  # fattore di scala grafica

# ----------------------
# INGRESSO STATORE
# ----------------------
plt.arrow(0, y0, 0, 0, head_width=2)
plt.text(5, y0 + 5, r"$\vec{c}_0$")
plt.text(-40, y0, "Ingresso statore", fontsize=10, va="center")

# ----------------------
# USCITA STATORE / INGRESSO ROTORE
# ----------------------
# c1
plt.arrow(0, y1, c1_t*scale, c1_a*scale,
          head_width=3, length_includes_head=True)
plt.text(c1_t*1.05, y1 + c1_a*1.05, r"$\vec{c}_1$")

# U
plt.arrow(0, y1, u_opt*scale, 0,
          head_width=3, length_includes_head=True)
plt.text(u_opt*1.05, y1, r"$\vec{u}$")

# w1 = c1 - U
plt.arrow(u_opt*scale, y1, (c1_t - u_opt)*scale, c1_a*scale,
          head_width=3, length_includes_head=True)
plt.text(u_opt + (c1_t - u_opt)*1.05,
         y1 + c1_a*1.05, r"$\vec{w}_1$")

plt.text(-40, y1, "Uscita statore\nIngresso rotore", fontsize=10, va="center")

# ----------------------
# USCITA ROTORE
# ----------------------
# c2
plt.arrow(0, y2, c2_t*scale, w2_a*scale,
          head_width=3, length_includes_head=True)
plt.text(c2_t*1.05, y2 + w2_a*1.05, r"$\vec{c}_2$")

# U
plt.arrow(0, y2, u_opt*scale, 0,
          head_width=3, length_includes_head=True)
plt.text(u_opt*1.05, y2, r"$\vec{u}$")

# w2
plt.arrow(u_opt*scale, y2, w2_t*scale, w2_a*scale,
          head_width=3, length_includes_head=True)
plt.text(u_opt + w2_t*1.05,
         y2 + w2_a*1.05, r"$\vec{w}_2$")

plt.text(-40, y2, "Uscita rotore", fontsize=10, va="center")

# ----------------------
# IMPOSTAZIONI GRAFICHE
# ----------------------
plt.axvline(0, linewidth=0.8)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel("Componente tangenziale [m/s]")
plt.ylabel("Direzione assiale →")
plt.title("Triangoli di velocità – stadio turbina assiale ORC (MDM)")
plt.grid(True)

plt.show()
