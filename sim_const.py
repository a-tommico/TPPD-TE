from pyfluids import FluidsList

# --- 1. Working Fluid Selection ---
FLUID_NAME = FluidsList.MDM

# --- 2. ORC Cycle Boundary Conditions ---
T_COND = 100.0  # °C
T_EVAP = 270.0  # °C
MDM_MFR = 2.11  # kg/s (Reference mass flow)

# --- 3. Component Efficiencies ---
ETA_TURB = 0.85
ETA_PUMP = 0.70
ETA_REGEN = 0.80
ETA_BOILER = 0.90 # Global boiler efficiency (approx)
ETA_COMB = 0.99   # Combustion efficiency (e)
TH_LOSSES = 0.02  # Thermal losses in boiler

# --- 4. Thermal Oil Circuit ---
T_OIL_SUPPLY = 310.0 # °C
T_OIL_RETURN = 280.0 # °C
CP_OIL = 2.48e3      # J/kgK (Generic thermal oil)

# --- 5. Boiler & Flue Gas Parameters ---
T_AMB = 25.0         # °C
T_STACK_TARGET = 150.0 # °C
CP_FLUE = 1.1e3      # J/kgK
CP_AIR = 1.005e3     # J/kgK
AIR_EXCESS = 0.30    # 30% excess air
PINCH_POINT_BOILER = 20.0 # K

# --- 6. Fuel Properties ---
# Mass (kg)
M_DRY = 800e3
M_WETDRIED = 259e3
M_SLUDGEDRIED = 70e3
M_TOT = M_DRY + M_WETDRIED + M_SLUDGEDRIED

# LHV (J/kg)
LHV_DRY_REF = 19.01e6
LHV_WET_REF = 19.15e6
LHV_SLUDGE_REF = 18.72e6

H_EVAP = 2.26e6 # J/kg (Water enthalpy of vaporization)

# Moisture content
HUM_DRY = 0.0663
HUM_WET = 0.15
HUM_SLUDGE = 0.15

# Fuel composition (Molar)
FUEL_COMP = {
    'C': 1.0,
    'H': 1.67,
    'O': 0.51,
    'N': 0.048
}
