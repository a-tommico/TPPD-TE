from pyfluids import Fluid, FluidsList, Input
import matplotlib.pyplot as plt
import numpy as np

fluid = Fluid(FluidsList.MDM)
mass_flow_rate = 3.0
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

# Output Data
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

# Fixed Plotting Logic
plt.figure(figsize=(10, 7))

t_crit = fluid.critical_temperature
T_range = np.linspace(0, t_crit, 300)

s_liq, T_liq_plot, h_liq, p_liq_plot = [], [], [], []
s_vap, T_vap_plot, h_vap, p_vap_plot = [], [], [], []

for T in T_range:
    try:
        # Calculate State (Liquid)
        f = fluid.with_state(Input.temperature(T), Input.quality(0))
        s_liq.append(f.entropy / 1000)      # kJ/kgK
        T_liq_plot.append(f.temperature)    # °C
        h_liq.append(f.enthalpy / 1000)     # kJ/kg
        p_liq_plot.append(f.pressure / 1e5) # bar
    except:
        pass
    try:
        # Calculate State (Vapor)
        f = fluid.with_state(Input.temperature(T), Input.quality(100))
        s_vap.append(f.entropy / 1000)      # kJ/kgK
        T_vap_plot.append(f.temperature)    # °C
        h_vap.append(f.enthalpy / 1000)     # kJ/kg
        p_vap_plot.append(f.pressure / 1e5) # bar
    except:
        pass

# Plot the Curves
plt.plot(s_liq, T_liq_plot, 'k-', linewidth=1.5, label='Saturated Liquid (Q=0)')
plt.plot(s_vap, T_vap_plot, 'k--', linewidth=1.5, label='Saturated Vapor (Q=1)')

# Add saturation points to the plot for clarity
state_sat_liq = fluid.with_state(Input.pressure(p_evap), Input.quality(0))
state_sat_vap = fluid.with_state(Input.pressure(p_cond), Input.quality(100))

cycle_states = [state1, state2, state2_prime, state_sat_liq, state3, state4, state5, state_sat_vap, state1]
s_cycle = [st.entropy / 1000 for st in cycle_states]
T_cycle = [st.temperature for st in cycle_states]
plt.plot(s_cycle, T_cycle, 'b-o', linewidth=2, label='ORC Cycle')
labels = ['1', '2', "2'", '', '3', '4', '5', '']
offsets = [(-15, -15), (-15, 10), (-15, 10), (5, 5), (0, 10), (10, 10), (10, -15), (5, -10)]

for i, txt in enumerate(labels):
    plt.annotate(txt, (s_cycle[i], T_cycle[i]), xytext=offsets[i], 
                 textcoords='offset points', fontweight='bold', fontsize=11, color='blue')

# Condenser Pressure (Low)
s_range_low = np.linspace(min(s_liq), max(s_vap), 50)
T_iso_low = []
for s in s_range_low:
    try:
        f = fluid.with_state(Input.pressure(p_cond), Input.entropy(s * 1000))
        T_iso_low.append(f.temperature)
    except: T_iso_low.append(None)
plt.plot(s_range_low, T_iso_low, 'g:', linewidth=2, alpha=0.7, label=f'P_cond ({p_cond/1e5:.2f} bar)')

# Evaporator Pressure (High)
s_range_high = np.linspace(min(s_liq), max(s_vap), 50)
T_iso_high = []
for s in s_range_high:
    try:
        f = fluid.with_state(Input.pressure(p_evap), Input.entropy(s * 1000))
        T_iso_high.append(f.temperature)
    except: T_iso_high.append(None)
plt.plot(s_range_high, T_iso_high, 'r:', linewidth=2, alpha=0.7, label=f'P_evap ({p_evap/1e5:.2f} bar)')

# 6. Formatting
plt.xlabel('Entropy (kJ/kg·K)')
plt.ylabel('Temperature (°C)')
plt.title('T-s Diagram: Regenerative ORC (MDM)')
plt.legend(loc='upper left')
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.ylim(bottom=0, top=t_crit + 20) # Start Y-axis at 0°C
plt.tight_layout()

# --- 7. p-h Diagram (New) ---
plt.figure(figsize=(10, 7))

# Plot Saturation Curves
plt.plot(h_liq, p_liq_plot, 'k-', linewidth=1.5, label='Saturated Liquid')
plt.plot(h_vap, p_vap_plot, 'k--', linewidth=1.5, label='Saturated Vapor')

# Cycle Points for p-h
h_cycle = [st.enthalpy / 1000 for st in cycle_states]
p_cycle = [st.pressure / 1e5 for st in cycle_states]

plt.plot(h_cycle, p_cycle, 'b-o', linewidth=2, label='ORC Cycle')

# Annotate points on p-h
for i, txt in enumerate(labels):
    plt.annotate(txt, (h_cycle[i], p_cycle[i]), xytext=offsets[i], 
                 textcoords='offset points', fontweight='bold', fontsize=11, color='blue')

plt.xlabel('Enthalpy (kJ/kg)')
plt.ylabel('Pressure (bar)')
plt.title('p-h Diagram: Regenerative ORC (MDM)')
plt.legend(loc='upper left')
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.yscale('log')
plt.tight_layout()
plt.show()