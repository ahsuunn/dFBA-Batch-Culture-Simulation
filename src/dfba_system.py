"""
ODE State vector y = [X, Glc, O2, Ace]
  X   : biomass   [g_DW / L]
  Glc : glucose   [mmol / L]
  O2  : oxygen    [mmol / L]
  Ace : acetate   [mmol / L]

The infeasibility event fires when BOTH glucose AND oxygen drop below EPSILON,
signalling complete substrate exhaustion (end of batch culture).
"""

import logging
import numpy as np
import cobra
from src.config import (  # type: ignore
    EPSILON,
    REACTION_GLUCOSE,
    REACTION_OXYGEN,
    REACTION_ACETATE,
)
from src.dynamic_bounds import add_dynamic_bounds  # type: ignore

logger = logging.getLogger(__name__)

# Module-level cache for the last FBA result
_last_solution: dict = {
    "mu": 0.0,
    "v_glc": 0.0,
    "v_o2": 0.0,
    "v_ace": 0.0,
    "feasible": True,
}


def _run_fba(model: cobra.Model) -> dict:
    """Solve the FBA LP and return key fluxes."""
    solution = model.optimize()
    if solution.status != "optimal":
        return {"mu": 0.0, "v_glc": 0.0, "v_o2": 0.0, "v_ace": 0.0, "feasible": False}

    fluxes = solution.fluxes
    return {
        "mu": max(0.0, solution.objective_value),
        "v_glc": max(0.0, -fluxes.get(REACTION_GLUCOSE, 0.0)),
        "v_o2": max(0.0, -fluxes.get(REACTION_OXYGEN, 0.0)),
        "v_ace": fluxes.get(REACTION_ACETATE, 0.0),
        "feasible": True,
    }


def dynamic_system(t: float, y: list, model: cobra.Model) -> list:
    """
    ODE right-hand side for the dFBA system.

    dX/dt   =  mu  * X
    dGlc/dt = -v_glc * X
    dO2/dt  = -v_o2  * X
    dAce/dt =  v_ace * X
    """
    global _last_solution

    X, Glc, O2, Ace = y
    X = max(X, 0.0)
    Glc = max(Glc, 0.0)
    O2 = max(O2, 0.0)
    Ace = max(Ace, 0.0)

    with model:
        add_dynamic_bounds(model, [X, Glc, O2, Ace])
        result = _run_fba(model)

    _last_solution = result

    if not result["feasible"] or X < EPSILON:
        return [0.0, 0.0, 0.0, 0.0]

    mu = result["mu"]
    v_glc = result["v_glc"]
    v_o2 = result["v_o2"]
    v_ace = result["v_ace"]

    return [mu * X, -v_glc * X, -v_o2 * X, v_ace * X]


def glucose_event(t: float, y: list, model: cobra.Model) -> float:
    """
    Terminal event: fires when glucose drops below EPSILON.
    Returns glucose concentration (positive initially, crosses zero on depletion).
    """
    return y[1] - EPSILON


glucose_event.terminal = True
glucose_event.direction = -1


def oxygen_event(t: float, y: list, model: cobra.Model) -> float:
    """
    Non-terminal event: fires when O2 drops below EPSILON.
    Used to mark the aerobic->anaerobic transition time.
    """
    return y[2] - EPSILON


oxygen_event.terminal = False
oxygen_event.direction = -1
