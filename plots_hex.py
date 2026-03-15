import matplotlib.pyplot as plt
import numpy as np
from pyfluids import Fluid, Input
from sim_const import *
from td_cycle import calculate_cycle
from boiler import solve_boiler, calculate_fuel_properties

def plot_hex_profiles():
    # 1. Get Simulation Data
    orc = calculate_cycle()
    boiler = solve_boiler(orc['Q_in'])
    fuel_props = calculate_fuel_properties()
    fluid = orc['fluid']
    m_mdm = MDM_MFR
    
    # 2. Boiler / Oil Heater (Flue Gas -> Thermal Oil)
    T_flame = T_AMB + (boiler['m_dot_fuel'] * fuel_props['lhv_mean'] * ETA_COMB * (1-TH_LOSSES)) / (boiler['m_dot_flue'] * CP_FLUE)
    q_boiler = boiler['Q_to_oil']
    # Counter-current: Fumes In (T_flame) meets Oil Out (T_OIL_SUPPLY)
    fumes_T_b = [T_flame, boiler['T_flue_post_oil']]
    oil_T_b = [T_OIL_SUPPLY, T_OIL_RETURN]
    
    # 3. Evaporator (Oil -> MDM)
    p_evap = orc['p_evap']
    st2_prime = orc['states']["2'"]
    st3 = orc['states']['3']
    st_sat_liq_evap = fluid.with_state(Input.pressure(p_evap), Input.quality(0))
    
    Q_preheat = m_mdm * (st_sat_liq_evap.enthalpy - st2_prime.enthalpy)
    Q_evap = m_mdm * (st3.enthalpy - st_sat_liq_evap.enthalpy)
    Q_total_evap = Q_preheat + Q_evap
    
    # Oil T profile
    T_oil_mid = T_OIL_SUPPLY - (Q_evap / Q_total_evap) * (T_OIL_SUPPLY - T_OIL_RETURN)
    # Counter-current: Oil In (SUPPLY) meets MDM Out (st3)
    mdm_T_ev = [st3.temperature, st_sat_liq_evap.temperature, st2_prime.temperature]
    oil_T_ev = [T_OIL_SUPPLY, T_oil_mid, T_OIL_RETURN]
    
    # 4. Regenerator (MDM Gas -> MDM Liquid)
    st4 = orc['states']['4']
    st5 = orc['states']['5']
    st2 = orc['states']['2']
    Q_regen = m_mdm * (st4.enthalpy - st5.enthalpy)
    # Counter-current: Hot In (st4) meets Cold Out (st2_prime)
    hot_T_reg = [st4.temperature, st5.temperature]
    cold_T_reg = [st2_prime.temperature, st2.temperature]

    # 5. Condenser (MDM -> Water)
    p_cond = orc['p_cond']
    st1 = orc['states']['1']
    st_sat_vap_cond = fluid.with_state(Input.pressure(p_cond), Input.quality(100))
    
    Q_desuper = m_mdm * (st5.enthalpy - st_sat_vap_cond.enthalpy)
    Q_cond = m_mdm * (st_sat_vap_cond.enthalpy - st1.enthalpy)
    Q_total_cond = Q_desuper + Q_cond
    
    T_w_in = 60
    T_w_out = 90
    T_w_mid = T_w_out - (Q_desuper / Q_total_cond) * (T_w_out - T_w_in)
    # Counter-current: MDM In (st5) meets Water Out (90)
    mdm_T_co = [st5.temperature, st_sat_vap_cond.temperature, st1.temperature]
    water_T_co = [T_w_out, T_w_mid, T_w_in]

    # 6. Air Preheater (APH: Fumes -> Air)
    Q_aph = boiler['m_dot_flue'] * CP_FLUE * (boiler['T_flue_post_oil'] - T_STACK_TARGET)
    # Counter-current: Fumes In (T_flue_post_oil) meets Air Out (T_air_preheated)
    fumes_T_aph = [boiler['T_flue_post_oil'], T_STACK_TARGET]
    air_T_aph = [boiler['T_air_preheated'], T_AMB]

    # --- PLOTTING ---
    fig, axs = plt.subplots(3, 2, figsize=(14, 15))
    
    # 1. Boiler
    axs[0, 0].plot([0, q_boiler/1e3], fumes_T_b, 'r-o', label='Flue Gas')
    axs[0, 0].plot([0, q_boiler/1e3], oil_T_b, 'b-o', label='Thermal Oil')
    axs[0, 0].set_title("Boiler: Fumes -> Oil")
    
    # 2. Evaporator
    q_coords_evap = [0, Q_evap/1e3, Q_total_evap/1e3]
    axs[0, 1].plot(q_coords_evap, oil_T_ev, 'r-o', label='Thermal Oil')
    axs[0, 1].plot(q_coords_evap, mdm_T_ev, 'b-o', label='MDM')
    axs[0, 1].set_title("Evaporator: Oil -> MDM")
    
    # 3. Regenerator
    axs[1, 0].plot([0, Q_regen/1e3], hot_T_reg, 'r-o', label='Hot MDM (Gas)')
    axs[1, 0].plot([0, Q_regen/1e3], cold_T_reg, 'b-o', label='Cold MDM (Liq)')
    axs[1, 0].set_title("Regenerator: MDM -> MDM")

    # 4. Condenser
    q_coords_cond = [0, Q_desuper/1e3, Q_total_cond/1e3]
    axs[1, 1].plot(q_coords_cond, mdm_T_co, 'r-o', label='MDM')
    axs[1, 1].plot(q_coords_cond, water_T_co, 'b-o', label='Water')
    axs[1, 1].set_title("Condenser: MDM -> Water")

    # 5. Air Preheater
    axs[2, 0].plot([0, Q_aph/1e3], fumes_T_aph, 'r-o', label='Flue Gas')
    axs[2, 0].plot([0, Q_aph/1e3], air_T_aph, 'b-o', label='Combustion Air')
    axs[2, 0].set_title("APH: Fumes -> Air")
    
    # Hide the empty 6th subplot
    axs[2, 1].axis('off')

    for ax in axs.flat:
        if ax.axison:
            ax.set_xlabel("Cumulative Heat Exchanged (kW)")
            ax.set_ylabel("Temperature (°C)")
            ax.legend()
            ax.grid(True, alpha=0.3)
        
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_hex_profiles()
