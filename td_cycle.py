from pyfluids import Fluid, FluidsList, Input
import matplotlib.pyplot as plt
import numpy as np
from sim_const import FLUID_NAME, T_COND, T_EVAP, MDM_MFR, ETA_PUMP, ETA_TURB, ETA_REGEN

def calculate_cycle(fluid_name=FLUID_NAME, T_cond=T_COND, T_evap=T_EVAP, mass_flow_rate=MDM_MFR, 
                    pump_eff=ETA_PUMP, turb_eff=ETA_TURB, regen_eff=ETA_REGEN):
    fluid = Fluid(fluid_name)

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
    # Q_max = min(C_hot, C_cold) * (T4 - T2)
    C_hot = mass_flow_rate * state4.specific_heat
    C_cold = mass_flow_rate * state2.specific_heat
    C_min = min(C_hot, C_cold)

    Q_max = C_min * (state4.temperature - state2.temperature)
    Q_regen = regen_eff * Q_max

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

    return {
        'states': {'1': state1, '2': state2, "2'": state2_prime, '3': state3, '4': state4, '5': state5},
        'W_turb': W_turb,
        'W_pump': W_pump,
        'W_net': W_net,
        'Q_in': Q_in,
        'Q_out': Q_out,
        'Q_regen': Q_regen,
        'efficiency': efficiency,
        'fluid': fluid,
        'p_cond': p_cond,
        'p_evap': p_evap
    }

def plot_cycle(results):
    fluid = results['fluid']
    states_dict = results['states']
    p_evap = results['p_evap']
    p_cond = results['p_cond']
    
    # Extract states for plotting (ordered for the cycle loop)
    # 1 -> 2 -> 2' -> sat_liq_evap -> 3 -> 4 -> 5 -> sat_vap_cond -> 1
    state_sat_liq_evap = fluid.with_state(Input.pressure(p_evap), Input.quality(0))
    state_sat_vap_cond = fluid.with_state(Input.pressure(p_cond), Input.quality(100))
    
    cycle_points = [
        states_dict['1'], states_dict['2'], states_dict["2'"], 
        state_sat_liq_evap, states_dict['3'], states_dict['4'], 
        states_dict['5'], state_sat_vap_cond, states_dict['1']
    ]
    
    s_cycle = [st.entropy / 1000 for st in cycle_points]
    T_cycle = [st.temperature for st in cycle_points]
    h_cycle = [st.enthalpy / 1000 for st in cycle_points]
    p_cycle = [st.pressure / 1e5 for st in cycle_points]

    # --- T-s Diagram ---
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    
    t_crit = fluid.critical_temperature
    T_range = np.linspace(0, t_crit, 200)
    s_liq, s_vap = [], []
    for T in T_range:
        try:
            s_liq.append(fluid.with_state(Input.temperature(T), Input.quality(0)).entropy / 1000)
            s_vap.append(fluid.with_state(Input.temperature(T), Input.quality(100)).entropy / 1000)
        except: break
    
    plt.plot(s_liq, T_range[:len(s_liq)], 'k-')
    plt.plot(s_vap, T_range[:len(s_vap)], 'k--')
    plt.plot(s_cycle, T_cycle, 'b-o', label='ORC Cycle')

    # per-point offsets to avoid overlap with curve and each other
    ts_offsets = {'1': (-18, -12), '2': (5, -12), "2'": (5, 5), '3': (5, 5), '4': (5, -12), '5': (-18, 5)}
    for key, st in states_dict.items():
        ox, oy = ts_offsets[key]
        plt.annotate(key, (st.entropy/1000, st.temperature),
                     textcoords="offset points", xytext=(ox, oy), fontsize=9,
                     arrowprops=dict(arrowstyle="-", color='gray', lw=0.8))

    plt.xlabel('Entropy (kJ/kg·K)')
    plt.ylabel('Temperature (°C)')
    plt.title('T-s Diagram')
    plt.grid(True, alpha=0.3)

    # --- p-h Diagram ---
    plt.subplot(1, 2, 2)
    h_liq, h_vap = [], []
    for T in T_range:
        try:
            h_liq.append(fluid.with_state(Input.temperature(T), Input.quality(0)).enthalpy / 1000)
            h_vap.append(fluid.with_state(Input.temperature(T), Input.quality(100)).enthalpy / 1000)
        except: break
        
    plt.plot(h_liq, [fluid.with_state(Input.temperature(T), Input.quality(0)).pressure/1e5 for T in T_range[:len(h_liq)]], 'k-')
    plt.plot(h_vap, [fluid.with_state(Input.temperature(T), Input.quality(100)).pressure/1e5 for T in T_range[:len(h_vap)]], 'k--')
    plt.plot(h_cycle, p_cycle, 'r-o', label='ORC Cycle')

    ph_offsets = {'1': (-18, -12), '2': (5, 5), "2'": (5, 5), '3': (5, 5), '4': (5, -12), '5': (-18, -12)}
    for key, st in states_dict.items():
        ox, oy = ph_offsets[key]
        plt.annotate(key, (st.enthalpy/1000, st.pressure/1e5),
                     textcoords="offset points", xytext=(ox, oy), fontsize=9,
                     arrowprops=dict(arrowstyle="-", color='gray', lw=0.8))

    plt.yscale('log')
    plt.xlabel('Enthalpy (kJ/kg)')
    plt.ylabel('Pressure (bar)')
    plt.title('p-h Diagram')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    res = calculate_cycle()
    st = res['states']

    print("=" * 55)
    print(f"{'STATE TABLE':^55}")
    print("=" * 55)
    print(f"{'State':<8} {'T (°C)':<10} {'P (bar)':<10} {'h (kJ/kg)':<12} {'s (kJ/kgK)':<12}")
    print("-" * 55)
    for name, s in st.items():
        print(f"{name:<8} {s.temperature:<10.2f} {s.pressure/1e5:<10.3f} {s.enthalpy/1e3:<12.2f} {s.entropy/1e3:<12.4f}")

    print("=" * 55)
    print(f"{'COMPONENT  Δh  (kJ/kg)':<55}")
    print("-" * 55)
    print(f"  Pump       Δh = {(st['2'].enthalpy  - st['1'].enthalpy)/1e3:+.2f} kJ/kg")
    print(f"  Regen (c)  Δh = {(st['2\''].enthalpy - st['2'].enthalpy)/1e3:+.2f} kJ/kg")
    print(f"  Evaporator Δh = {(st['3'].enthalpy  - st['2\''].enthalpy)/1e3:+.2f} kJ/kg")
    print(f"  Turbine    Δh = {(st['3'].enthalpy  - st['4'].enthalpy)/1e3:+.2f} kJ/kg")
    print(f"  Regen (h)  Δh = {(st['4'].enthalpy  - st['5'].enthalpy)/1e3:+.2f} kJ/kg")
    print(f"  Condenser  Δh = {(st['5'].enthalpy  - st['1'].enthalpy)/1e3:+.2f} kJ/kg")

    print("=" * 55)
    print(f"{'POWER & HEAT  (kW)':<55}")
    print("-" * 55)
    print(f"  W_turb  = {res['W_turb']/1e3:.2f} kW")
    print(f"  W_pump  = {res['W_pump']/1e3:.2f} kW")
    print(f"  W_net   = {res['W_net']/1e3:.2f} kW")
    print(f"  Q_in    = {res['Q_in']/1e3:.2f} kW")
    print(f"  Q_regen = {res['Q_regen']/1e3:.2f} kW")
    print(f"  Q_out   = {res['Q_out']/1e3:.2f} kW")
    print(f"  P ratio = {res['p_evap']/res['p_cond']:.2f}")
    print(f"  η_cycle = {res['efficiency']:.2f} %")
    print("=" * 55)

    plot_cycle(res)
