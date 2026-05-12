"""
Translates the current external metabolite concentrations into
flux bounds for the COBRApy metabolic model at each integration step.

Biological rationale:
  Real bacterial cells exhibit saturable nutrient uptake described by
  Michaelis-Menten (MM) kinetics.  As substrate concentration [S] falls,
  the maximum achievable uptake rate (v) declines smoothly:

          v = Vmax * [S] / (Km + [S])

  This prevents numerical artefacts where the model would allow infinite
  consumption of a nearly-exhausted substrate.

  Exchange reactions in COBRApy use a sign convention where IMPORT
  (consumption) is represented by a NEGATIVE flux.  Therefore the computed
  kinetic rate is negated before being assigned as the lower_bound.
"""

import logging
import cobra
from src.config import (  # type: ignore
    VMAX_GLC,
    KM_GLC,
    VMAX_O2,
    KM_O2,
    EPSILON,
    REACTION_GLUCOSE,
    REACTION_OXYGEN,
)

logger = logging.getLogger(__name__)


def _michaelis_menten(vmax: float, km: float, concentration: float) -> float:
    """
    Compute Michaelis-Menten uptake rate.

    Parameters
    ----------
    vmax : float
        Maximum uptake rate  [mmol / (g_DW · h)]
    km : float
        Half-saturation constant  [mmol / L]
    concentration : float
        Current substrate concentration  [mmol / L]

    Returns
    -------
    float
        Uptake rate ≥ 0.  Returns 0 when concentration is at or below
        the numerical epsilon floor.
    """
    if concentration <= EPSILON:
        return 0.0
    return vmax * concentration / (km + concentration)


def add_dynamic_bounds(model: cobra.Model, y: list) -> None:
    """
    Update the exchange-reaction lower bounds of *model* using the current
    external metabolite concentrations in the state vector *y*.

    State vector layout
    -------------------
    y[0] : biomass   [g_DW / L]   (not used for bound calculation)
    y[1] : glucose   [mmol / L]
    y[2] : oxygen    [mmol / L]
    y[3] : acetate   [mmol / L]   (not used for bound calculation here;
                                   COBRApy chooses secretion optimally)

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to modify in-place.
    y : list or array-like
        Current ODE state vector.
    """
    _, glucose, oxygen, _ = y

    # --- Glucose ---
    glc_rate = _michaelis_menten(VMAX_GLC, KM_GLC, glucose)
    model.reactions.get_by_id(REACTION_GLUCOSE).lower_bound = -glc_rate
    logger.debug(
        "  [Glc=%.4f mmol/L]  uptake bound = %.4f mmol/(gDW·h)", glucose, glc_rate
    )

    # --- Oxygen ---
    o2_rate = _michaelis_menten(VMAX_O2, KM_O2, oxygen)
    model.reactions.get_by_id(REACTION_OXYGEN).lower_bound = -o2_rate
    logger.debug(
        "  [O2 =%.4f mmol/L]  uptake bound = %.4f mmol/(gDW·h)", oxygen, o2_rate
    )
