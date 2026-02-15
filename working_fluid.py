from pyfluids import Fluid, FluidsList, Input
import pandas as pd

# Boundary Conditions
T_oil_supply = 310  # °C (Limit of Thermal Oil)
T_pinch_evap = 30   # °C (Min delta T for Heat Exchanger)
T_max_cycle_limit = T_oil_supply - T_pinch_evap # 280°C

T_water_req = 90    # °C (Required for Dryer)
T_pinch_cond = 10   # °C
T_cond = T_water_req + T_pinch_cond # 100°C

# Component Efficiencies
eta_turb = 0.85
eta_pump = 0.70
regen_eff = 0.80

# --- 2. Define Candidates ---
candidates = [
    {"name": "MDM", "fluid": FluidsList.MDM, "type": "Siloxane"},
    {"name": "MM", "fluid": FluidsList.MM, "type": "Siloxane"},
    {"name": "D4", "fluid": FluidsList.D4, "type": "Siloxane (Heavy)"},
    {"name": "Toluene", "fluid": FluidsList.Toluene, "type": "Hydrocarbon"},
    {"name": "Cyclopentane", "fluid": FluidsList.CycloPentane, "type": "Hydrocarbon"},
]

results = []

print(f"--- Simulation Boundary Conditions ---")
print(f"Source Limit (Oil): {T_oil_supply}°C -> Max Cycle T: {T_max_cycle_limit}°C")
print(f"Sink Target (Dryer): {T_water_req}°C -> Condensing T: {T_cond}°C")
print("-" * 60)

for item in candidates:
    try:
        fluid = Fluid(item["fluid"])
        name = item["name"]
        
        t_crit = fluid.critical_temperature
        t_evap_target = min(T_max_cycle_limit, t_crit - 20)
        
        if t_evap_target <= T_cond + 20:
            print(f"Skipping {name}: Critical temp too low for 100°C condensation.")
            continue

        # State 1: Condenser Outlet (Sat Liquid)
        st1 = fluid.with_state(Input.temperature(T_cond), Input.quality(0))
        p_cond = st1.pressure
        
        # State 3: Turbine Inlet (Sat Vapor at T_evap)
        st3 = fluid.with_state(Input.temperature(t_evap_target), Input.quality(1))
        p_evap = st3.pressure
        
        # State 2: Pump Outlet
        st2 = st1.compression_to_pressure(p_evap, eta_pump*100)
        w_pump = st2.enthalpy - st1.enthalpy
        
        # State 4: Turbine Outlet
        st4 = st3.expansion_to_pressure(p_cond, eta_turb*100)
        w_turb = st3.enthalpy - st4.enthalpy
        
        # Regenerator
        # Max heat recoverable: from T4 down to T2
        q_regen_max = st4.enthalpy - st4.with_state(Input.temperature(st2.temperature), Input.pressure(p_cond)).enthalpy
        # Actually we simplify: Q_regen = eff * Cp_min * dT. 
        # Liquid heats up from h2 to h2'
        # Vapor cools from h4 to h5
        # Energy balance: h2' - h2 = h4 - h5.
        # Max T limit is T4.
        
        # Simplified Check for Dry/Wet fluid
        expansion_type = "Dry" if st4.entropy >= st3.entropy else "Wet" 
        # Note: Ideally compare to Sat Vapor entropy at P_cond, but for ORC usually Dry.
        
        # Efficiency Calculation
        q_in = st3.enthalpy - st2.enthalpy # Without regen for raw comparison first
        # With Regen (approximate):
        # We assume we can heat liquid up to T4 - 10°C (approach) if regen is used
        t_regen_liquid_out = min(st4.temperature - 10, t_evap_target - 5)
        h2_regen = fluid.with_state(Input.temperature(t_regen_liquid_out), Input.pressure(p_evap)).enthalpy
        
        # Check if Regen is actually possible (T4 > T2)
        if st4.temperature > st2.temperature + 10:
             q_in_regen = st3.enthalpy - h2_regen
             efficiency = (w_turb - w_pump) / q_in_regen
             regen_status = "YES"
        else:
             efficiency = (w_turb - w_pump) / q_in
             regen_status = "NO"

        # C. Store Results
        results.append({
            "Fluid": name,
            "Type": item["type"],
            "T_evap (°C)": round(t_evap_target, 1),
            "P_evap (bar)": round(p_evap/1e5, 2),
            "P_cond (bar)": round(p_cond/1e5, 2),
            "Press. Ratio": round(p_evap/p_cond, 1),
            "Efficiency (%)": round(efficiency * 100, 2),
            "T_crit (°C)": round(t_crit, 1),
            "Turbine Out T (°C)": round(st4.temperature, 1),
            "Expansion Type": expansion_type,
        })

    except Exception as e:
        print(f"Error calculating {item['name']}: {e}")

df = pd.DataFrame(results)
print("\n--- Fluid Comparison Matrix (Condensing at 100°C) ---")
print(df)