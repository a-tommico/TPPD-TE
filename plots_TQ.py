import matplotlib.pyplot as plt
from pyfluids import Input
from sim_const import (
    MDM_MFR,
    T_OIL_SUPPLY, T_OIL_RETURN,
    T_W_IN, T_W_OUT,
    PINCH_EVAP, PINCH_COND, PINCH_REGEN,
)
from td_cycle import calculate_cycle

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"


def _annotate_pinch(ax, q_kw, T_hot, T_cold, label="ΔT"):
    """Draw a vertical double-headed arrow at the pinch point."""
    ax.annotate(
        "", xy=(q_kw, T_cold), xytext=(q_kw, T_hot),
        arrowprops=dict(arrowstyle="<->", color="green", lw=1.5),
    )
    ax.text(q_kw, (T_hot + T_cold) / 2, f"  {label}\n  {T_hot - T_cold:.1f} K",
            color="green", fontsize=8, va="center")


def plot_hex_TQ():
    orc = calculate_cycle()
    fluid  = orc['fluid']
    m      = MDM_MFR
    p_evap = orc['p_evap']
    p_cond = orc['p_cond']

    st1  = orc['states']['1']
    st2  = orc['states']['2']
    st2p = orc['states']["2'"]
    st3  = orc['states']['3']
    st4  = orc['states']['4']
    st5  = orc['states']['5']

    _, axs = plt.subplots(1, 3, figsize=(18, 6))

    # ── 1. EVAPORATOR  (Oil → MDM, counter-current) ──────────────────────────
    # MDM path: 2' (subcooled) → sat-liq → sat-vap (state 3)
    st_sat_liq_evap = fluid.with_state(Input.pressure(p_evap), Input.quality(0))

    Q_evap_lat  = m * (st3.enthalpy              - st_sat_liq_evap.enthalpy)  # latent
    Q_evap_sens = m * (st_sat_liq_evap.enthalpy  - st2p.enthalpy)             # preheat
    Q_evap_tot  = Q_evap_lat + Q_evap_sens

    # x=0: hot end (oil in / MDM vapor out)  |  x=Q_evap_tot: cold end
    # Pinch is at the bubble-point boundary (x = Q_evap_lat)
    T_oil_mid = T_OIL_SUPPLY - (Q_evap_lat / Q_evap_tot) * (T_OIL_SUPPLY - T_OIL_RETURN)

    Q_ev     = [0,            Q_evap_lat/1e3,              Q_evap_tot/1e3]
    T_oil_ev = [T_OIL_SUPPLY, T_oil_mid,                   T_OIL_RETURN]
    T_mdm_ev = [st3.temperature, st_sat_liq_evap.temperature, st2p.temperature]

    axs[0].plot(Q_ev, T_oil_ev, 'r-o', label='Thermal Oil (hot)')
    axs[0].plot(Q_ev, T_mdm_ev, 'b-o', label='MDM (cold)')
    # Pinch annotation at bubble-point boundary
    _annotate_pinch(axs[0], Q_evap_lat/1e3,
                    T_oil_mid, st_sat_liq_evap.temperature,
                    label=f"PP={PINCH_EVAP:.0f} K req.")
    axs[0].set_title("Evaporator  (Oil → MDM)")
    axs[0].set_xlabel("Cumulative heat  Q  (kW)")
    axs[0].set_ylabel("Temperature  (°C)")
    axs[0].legend()
    axs[0].grid(True, alpha=0.3)

    # ── 2. REGENERATOR  (MDM vapor → MDM liquid, counter-current) ────────────
    # Hot side: st4 → st5  |  Cold side (counter-current): st2' → st2
    # Pinch at the end where streams are closest (check both ends)
    Q_regen = m * (st4.enthalpy - st5.enthalpy)

    dT_hot_end  = st4.temperature  - st2p.temperature   # x = 0
    dT_cold_end = st5.temperature  - st2.temperature    # x = Q_regen

    Q_rg      = [0,                Q_regen/1e3]
    T_hot_rg  = [st4.temperature,  st5.temperature]
    T_cold_rg = [st2p.temperature, st2.temperature]

    axs[1].plot(Q_rg, T_hot_rg,  'r-o', label='Hot MDM — vapor')
    axs[1].plot(Q_rg, T_cold_rg, 'b-o', label='Cold MDM — liquid')
    # Annotate the tighter end
    if dT_cold_end <= dT_hot_end:
        _annotate_pinch(axs[1], Q_regen/1e3,
                        st5.temperature, st2.temperature,
                        label=f"PP={PINCH_REGEN:.0f} K req.")
    else:
        _annotate_pinch(axs[1], 0,
                        st4.temperature, st2p.temperature,
                        label=f"PP={PINCH_REGEN:.0f} K req.")
    axs[1].set_title("Regenerator  (MDM → MDM)")
    axs[1].set_xlabel("Cumulative heat  Q  (kW)")
    axs[1].set_ylabel("Temperature  (°C)")
    axs[1].legend()
    axs[1].grid(True, alpha=0.3)

    # ── 3. CONDENSER  (MDM → Water, counter-current) ─────────────────────────
    # MDM path: st5 (superheated) → sat-vap → sat-liq (state 1)
    st_sat_vap_cond = fluid.with_state(Input.pressure(p_cond), Input.quality(100))

    Q_cond_desuper = m * (st5.enthalpy              - st_sat_vap_cond.enthalpy)
    Q_cond_lat     = m * (st_sat_vap_cond.enthalpy  - st1.enthalpy)
    Q_cond_tot     = Q_cond_desuper + Q_cond_lat

    # x=0: hot end (MDM in / water out)  |  x=Q_cond_tot: cold end
    # Pinch is at the dew-point boundary (x = Q_cond_desuper)
    T_w_mid = T_W_OUT - (Q_cond_desuper / Q_cond_tot) * (T_W_OUT - T_W_IN)

    Q_co       = [0,               Q_cond_desuper/1e3,          Q_cond_tot/1e3]
    T_mdm_co   = [st5.temperature, st_sat_vap_cond.temperature, st1.temperature]
    T_water_co = [T_W_OUT,         T_w_mid,                     T_W_IN]

    axs[2].plot(Q_co, T_mdm_co,   'r-o', label='MDM (hot)')
    axs[2].plot(Q_co, T_water_co, 'b-o', label='Water (cold)')
    # Pinch annotation at dew-point boundary
    _annotate_pinch(axs[2], Q_cond_desuper/1e3,
                    st_sat_vap_cond.temperature, T_w_mid,
                    label=f"PP={PINCH_COND:.0f} K req.")
    axs[2].set_title("Condenser  (MDM → Water)")
    axs[2].set_xlabel("Cumulative heat  Q  (kW)")
    axs[2].set_ylabel("Temperature  (°C)")
    axs[2].legend()
    axs[2].grid(True, alpha=0.3)

    for ax in axs:
        ax.set_box_aspect(1)

    plt.tight_layout()
    plt.savefig("/home/red/Desktop/utils/TPPD-TE/plots/hex_TQ.png", dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    plot_hex_TQ()