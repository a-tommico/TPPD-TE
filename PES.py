import matplotlib.pyplot as plt
import numpy as np
from sim_const import LHV_MEAN, M_TOT
from td_cycle import calculate_cycle

LHV_mean = LHV_MEAN

# -----------------------------------------
# ORC CYCLE  (from td_cycle)
# -----------------------------------------
orc = calculate_cycle()
W_net = orc['W_net']   # W  — net electrical power
Q_out = orc['Q_out']   # W  — condenser thermal output

# -----------------------------------------
# FUEL CONSUMPTION  (all available fuel over 1 year)
# -----------------------------------------
m_dot_fuel = M_TOT / (365 * 24 * 3600)   # kg/s
Q_in_fuel  = m_dot_fuel * LHV_mean        # W

# -----------------------------------------
# PES  (EU Directive 2012/27/EU definition)
# -----------------------------------------
eta_el_rif  = 0.30
eta_th_rif  = 0.80
eta_el_best = 0.42
eta_th_best = 0.90

eta_el_ORC = W_net / Q_in_fuel
eta_th_ORC = Q_out / Q_in_fuel

PES_ref  = 1 - 1 / (eta_el_ORC / eta_el_rif  + eta_th_ORC / eta_th_rif)
PES_best = 1 - 1 / (eta_el_ORC / eta_el_best + eta_th_ORC / eta_th_best)

print(f"LHV_mean:                  {LHV_mean/1e6:.3f} MJ/kg")
print(f"Yearly fuel consumption:   {M_TOT/1e3:.2f} tons/year")
print(f"W_net:                     {W_net/1e3:.2f} kW")
print(f"Q_out (thermal):           {Q_out/1e3:.2f} kW")
print(f"Q_in_fuel:                 {Q_in_fuel/1e3:.2f} kW")
print(f"η_el ORC:                  {eta_el_ORC*100:.2f} %")
print(f"η_th ORC:                  {eta_th_ORC*100:.2f} %")
print(f"PES (vs reference):        {PES_ref*100:.2f} %")
print(f"PES (vs best technology):  {PES_best*100:.2f} %")

# -----------------------------------------
# PLOT
# -----------------------------------------
x_ref  = np.linspace(0, eta_el_rif,  100)
x_best = np.linspace(0, eta_el_best, 100)
x_1    = np.linspace(0, 1, 100)

plt.figure(figsize=(6, 6))
plt.plot(x_ref,  eta_th_rif  * (1 - x_ref  / eta_el_rif),  label='Reference  PES = 0')
plt.plot(x_best, eta_th_best * (1 - x_best / eta_el_best), color='green', label='Best technology')
plt.plot(x_1, 1 - x_1, color='black', linestyle='--', label=r'$\eta_{tot} = 1$')

plt.scatter(eta_el_ORC, eta_th_ORC, color='red', zorder=5, label='ORC')
plt.hlines(eta_th_ORC, 0, eta_el_ORC, color='red', linestyle='--', linewidth=0.8)
plt.vlines(eta_el_ORC, 0, eta_th_ORC, color='red', linestyle='--', linewidth=0.8)

plt.xlabel(r'Electrical efficiency $\eta_{el}$')
plt.ylabel(r'Thermal efficiency $\eta_{th}$')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
