"""
All biological and numerical parameters for the dFBA simulation are
centralised here. Adjust these values to explore different batch culture
scenarios without touching any other module.

Biological basis:
  - Kinetic constants (Vmax, Km) follow Michaelis-Menten uptake kinetics
    commonly reported for E. coli core metabolism.
  - Initial concentrations mimic a typical M9 minimal-medium batch culture.
  - ATP yield estimates are derived from stoichiometric coefficients of the
    E. coli core metabolic model (aerobic ~29-32 ATP/glucose;
    anaerobic ~2 ATP/glucose via glycolysis alone).
"""

# Time span (hours)
T_START: float = 0.0  # Start time
T_END: float = 20.0  # End time (hours)
N_EVAL: int = 500  # Number of evaluation points for dense output

# Initial conditions  [biomass, glucose, oxygen, acetate]
# Units:  biomass → g_DW/L   |  metabolites → mmol/L
BIOMASS_0: float = 0.01  # g dry-weight / L
GLUCOSE_0: float = 10.0  # mmol / L
OXYGEN_0: float = 0.21  # mmol / L  (air-saturated aqueous ~0.21 mmol/L)
ACETATE_0: float = 0.0  # mmol / L  (no acetate at start)

# Convenience vector (order must match dynamic_system state vector)
Y0 = [BIOMASS_0, GLUCOSE_0, OXYGEN_0, ACETATE_0]

# Michaelis-Menten kinetic parameters for exchange reactions
# Glucose uptake  (EX_glc__D_e)
VMAX_GLC: float = 10.0  # mmol / (g_DW · h)  — maximum glucose uptake rate
KM_GLC: float = 0.015  # mmol / L            — half-saturation constant

# Oxygen uptake   (EX_o2_e)
VMAX_O2: float = 15.0  # mmol / (g_DW · h)  — maximum oxygen uptake rate
KM_O2: float = 0.001  # mmol / L            — half-saturation constant

# Acetate secretion / uptake exchange (EX_ac_e) — allowed to vary freely
#   COBRApy default bounds are used; secretion is driven by model optimisation.

# Numerical solver settings
METHOD: str = "BDF"  # Backward Differentiation Formula (stiff ODEs)
RTOL: float = 1e-6  # Relative tolerance
ATOL: float = 1e-8  # Absolute tolerance
EPSILON: float = 1e-8  # Substrate floor — triggers infeasible event

# ATP yield reference values (for annotation / annotation overlay on plots)
ATP_AEROBIC: float = 30.0  # Approximate ATP per glucose (aerobic)
ATP_ANAEROBIC: float = 2.0  # Approximate ATP per glucose (anaerobic / fermentation)

# Output directory for saved plots
RESULTS_DIR: str = "results"

# Exchange reaction IDs in the E. coli core / textbook model
REACTION_GLUCOSE: str = "EX_glc__D_e"
REACTION_OXYGEN: str = "EX_o2_e"
REACTION_ACETATE: str = "EX_ac_e"
REACTION_BIOMASS: str = "BIOMASS_Ecoli_core_w_GAM"  # Objective reaction
