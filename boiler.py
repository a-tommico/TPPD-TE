import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from pyfluids import Fluid, FluidsList, Input

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

# State 2'': Saturated Liquid at P = 10 bar
state2_second = fluid.with_state(Input.pressure(p_high), Input.quality(0))

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

# ---------------------
# BOILER SIZING
# ---------------------

# Input data
Q_ORC = Q_in # W
m_MDM = mass_flow_rate # kg/s
cp_to = 2.48e3 # J/kgK
eta_boiler = 0.90
h_evap = 2.26e6 # J/kg
e = 0.99 # combustion efficiency
air_excess = 0.30 # 30% excess air
th_losses = 0.02 # 2% thermal losses
cp_flue = 1.1e3 # J/kgK

# Fuel
m_dry = 800e3 # kg
m_wetdried = 259e3 # kg
m_sludgedried = 70e3 # kg
m_tot = m_dry + m_wetdried + m_sludgedried # kg

LHV_dry = 19.01e6 # J/kg
LHV_wetdry = 19.15e6 # J/kg
LHV_sludgedry = 18.72e6 # J/kg

LHV_dry_real = (1-0.0663)*LHV_dry - 0.0663*h_evap # J/kg
LHV_wet_real = (1-0.15)*LHV_wetdry - 0.15*h_evap # J/kg
LHV_sludge_real = (1-0.15)*LHV_sludgedry - 0.15*h_evap # J/kg

hum = (0.0663 * m_dry + 0.15 * m_wetdried + 0.15 * m_sludgedried) / m_tot

print(f"LHV dry = {LHV_dry_real/1e6:.2f} MJ/kg")
print(f"LHV wet = {LHV_wet_real/1e6:.2f} MJ/kg")
print(f"LHV sludge = {LHV_sludge_real/1e6:.2f} MJ/kg")

LHV_mean = (m_dry*LHV_dry_real + m_wetdried*LHV_wet_real + m_sludgedried*LHV_sludge_real) / m_tot # J/kg

# Combustion
fuel = {
    'C': 1.0,
    'H': 1.67,
    'O': 0.51,
    'N': 0.048
}
fuel = ct.Species(name="fuel", composition=fuel)
MM_fuel = fuel.molecular_weight # g/mol
print(f"Molecular weight fuel: {MM_fuel:.2f} g/mol")

MM_air = 28.96 # g/mol
print(f"Molecular weight air: {MM_air:.2f} g/mol")

st_O2 = (1 + 1.67/4 - 0.51/2) # mol O2 per mol fuel
st_air = st_O2 * (1 + 3.76) # mol air per mol fuel
print(f"Stoichiometric air-fuel ratio (molar): {st_air:.2f} mol air/mol fuel")

fuel_to_air_ratio = (MM_fuel) / (st_air * MM_air) # kg fuel per kg air
print(f"Stoichiometric air-fuel ratio (mass): {1/fuel_to_air_ratio:.2f} kg air/kg fuel")

# Boiler calculations
Q_boiler = Q_ORC / eta_boiler # W

m_dot_fuel = Q_boiler / (LHV_mean * e) # kg/s
m_dot_air = m_dot_fuel / fuel_to_air_ratio * (1 + air_excess) # kg/s
m_dot_flue = m_dot_fuel + m_dot_air # kg/s

P_th_flue = (e - eta_boiler - th_losses) * LHV_mean * m_dot_fuel # W
T_flue = (298.15 + P_th_flue / (m_dot_flue * cp_flue)) - 273.15 # °C
print(f"Flue gas temperature at the chimney: {T_flue:.2f} °C")

# Maximum thermal power from fuel
m_yearly = m_dot_fuel * 3600 * 24 * 365 # kg/year

m_dot = m_tot / (365*24*3600) # kg/s
Q_th = m_dot * LHV_mean # W

print(f"Q_ORC: {Q_ORC/1e3:.2f} kW")
print(f"Q_th: {Q_th/1e3:.2f} kW")
print(f"Q_ORC_theoretical: {Q_th*eta_boiler/1e3:.2f} kW")
print(f"LHV_mean: {LHV_mean/1e6:.2f} MJ/kg")
print(f"yearly fuel consumption: {m_yearly/1e3:.2f} tons/year")
print(f"available fuel mass: {m_tot/1e3:.2f} tons")

