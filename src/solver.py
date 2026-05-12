"""
Orchestrates the scipy.integrate.solve_ivp integration using the
Backward Differentiation Formula (BDF) method.

Events tracked:
  - glucose_event (terminal=True):  stops when glucose < epsilon
  - oxygen_event  (terminal=False): records when O2 < epsilon (transition marker)
"""

import logging
import numpy as np
from scipy.integrate import solve_ivp
import cobra

from src.config import Y0, T_START, T_END, N_EVAL, METHOD, RTOL, ATOL  # type: ignore
from src.dfba_system import dynamic_system, glucose_event, oxygen_event  # type: ignore

logger = logging.getLogger(__name__)


def run_dfba(model, y0=None, t_span=None, n_eval=None):
    """
    Run the dFBA simulation.

    Parameters
    ----------
    model : cobra.Model
    y0 : list, optional  [X0, Glc0, O2_0, Ace0]
    t_span : tuple, optional  (t_start, t_end) hours
    n_eval : int, optional  number of output time points

    Returns
    -------
    OdeResult
      .t          time array
      .y          state matrix  rows: [X, Glc, O2, Ace]
      .t_events   [0] = glucose depletion time, [1] = O2 depletion time
    """
    y0 = y0 or Y0
    t_span = t_span or (T_START, T_END)
    n_eval = n_eval or N_EVAL

    t_eval = np.linspace(t_span[0], t_span[1], n_eval)

    logger.info("=" * 60)
    logger.info("Starting dFBA simulation")
    logger.info("  Method : %s | t_span : %s h", METHOD, t_span)
    logger.info("  y0 : X=%.4f | Glc=%.2f | O2=%.4f | Ace=%.2f", *y0)
    logger.info("=" * 60)

    def _system(t, y):
        return dynamic_system(t, y, model)

    def _glc_event(t, y):
        return glucose_event(t, y, model)

    def _o2_event(t, y):
        return oxygen_event(t, y, model)

    _glc_event.terminal = True
    _glc_event.direction = -1
    _o2_event.terminal = False
    _o2_event.direction = -1

    solution = solve_ivp(
        fun=_system,
        t_span=t_span,
        y0=y0,
        method=METHOD,
        t_eval=t_eval,
        events=[_glc_event, _o2_event],
        rtol=RTOL,
        atol=ATOL,
        dense_output=True,
    )

    # Log summary
    if solution.t_events[0].size > 0:
        logger.info("Glucose exhausted at t = %.4f h", solution.t_events[0][0])
    if solution.t_events[1].size > 0:
        logger.info("Oxygen transition at t = %.4f h", solution.t_events[1][0])
    if not solution.t_events[0].size and not solution.t_events[1].size:
        logger.info("Integration completed to t = %.4f h", solution.t[-1])

    logger.info("Status: %s | Steps: %d", solution.message, len(solution.t))
    logger.info("=" * 60)

    return solution
