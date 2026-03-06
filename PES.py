from pyfluids import Fluid, FluidsList, Input
import matplotlib.pyplot as plt
import numpy as np

fluid = Fluid(FluidsList.MDM)
mass_flow_rate = 2.11
T_cond = 100
T_evap = 270
pump_eff = 0.70
turb_eff = 0.85

# State 1 (Pump Inlet): Saturated liquid
state1 = fluid.with_state(Input.temperature(T_cond), Input.quality(0))
p_cond = state1.pressure

# State 2 (Pump Outlet): Compressed to P_evap
p_evap = fluid.with_state(Input.temperature(T_evap), Input.quality(100)).pressure
state2 = state1.compression_to_pressure(p_evap, pump_eff*100)
w_pump = state2.enthalpy - state1.enthalpy

# State 3 (Turbine Inlet): Saturated Vapor
state3 = fluid.with_state(Input.pressure(p_evap), Input.quality(100))

# State 4 (Turbine Outlet): Expansion to P_cond
state4 = state3.expansion_to_pressure(p_cond, turb_eff*100)
w_turb = state3.enthalpy - state4.enthalpy

# Regenerator Logic
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
state2_prime = fluid.with_state(Input.pressure(p_evap), Input.enthalpy(h2_prime))

# Performance Metrics
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

# -----------------------------------------
# CALCOLO PES
# -----------------------------------------

eta_el_rif = 0.3
eta_th_rif = 0.8
eta_el_best = 0.42
eta_th_best = 0.9

LHV_mean = 17.09 # MJ/kg
yearly_fuel_consumption = 1124.96 # tons/year
m_dot_fuel = yearly_fuel_consumption * 1e3 / (365*24*3600) # kg/s
Q_in_fuel = m_dot_fuel * LHV_mean * 1e6 # W

P_el_ORC = W_net # W
P_th_ORC = Q_out # W

eta_el_ORC = P_el_ORC / Q_in_fuel
eta_th_ORC = P_th_ORC / Q_in_fuel

print(f"Electrical Efficiency of ORC: {eta_el_ORC*100:.2f} %")
print(f"Thermal Efficiency of ORC: {eta_th_ORC*100:.2f} %")

PES = 1 - 1 / (eta_el_ORC / eta_el_rif + eta_th_ORC / eta_th_rif)

print(f"Primary Energy Saving (PES) of ORC: {PES*100:.2f} %")

# retta di riferimento
x_ref = np.linspace(0, eta_el_rif, 100)
y_ref = eta_th_rif * (1 - x_ref/eta_el_rif)

# retta best (opzionale)
x_best = np.linspace(0, eta_el_best, 100)
y_best = eta_th_best * (1 - x_best/eta_el_best)

plt.figure(figsize=(6,6))

# rette
plt.plot(x_ref, y_ref, label='Reference PES = 0')
plt.plot(x_best, y_best, color='green', label='Best technology')

# punto ORC
plt.scatter(eta_el_ORC, eta_th_ORC, color='black', label='ORC')


plt.xlabel(r'Electrical efficiency $\eta_{el}$')
plt.ylabel(r'Thermal efficiency $\eta_{th}$')

plt.xlim(0,1)
plt.ylim(0,1)

plt.grid(True)
plt.legend()

plt.show()