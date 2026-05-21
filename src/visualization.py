"""
Generates a publication-quality multi-panel figure from dFBA results,
visualizing the aerobic-to-anaerobic respiratory transition in E. coli
batch culture.

Panel layout (2 rows × 3 columns):
  [0,0] Biomass growth (g_DW/L)
  [0,1] Glucose depletion (mmol/L)
  [0,2] Oxygen depletion (mmol/L)
  [1,0] Acetate secretion (mmol/L)  — fermentation byproduct
  [1,1] Estimated ATP yield per glucose (aerobic→anaerobic transition)
  [1,2] Specific growth rate µ (h⁻¹)
"""

import os
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

from src.config import (
    RESULTS_DIR,
    OXYGEN_0,
    ATP_AEROBIC,
    ATP_ANAEROBIC,
)

logger = logging.getLogger(__name__)

# Colour palette───
DARK_BG = "#0d1117"
PANEL_BG = "#161b22"
GRID_CLR = "#30363d"
TEXT_CLR = "#e6edf3"
ACCENT = "#58a6ff"

COLOURS = {
    "biomass": "#39d353",  # green
    "glucose": "#f78166",  # orange-red
    "oxygen": "#79c0ff",  # light blue
    "acetate": "#d2a8ff",  # purple
    "atp": "#ffa657",  # gold
    "mu": "#56d364",  # bright green
    "vline": "#ff7b72",  # red (event line)
}


def _apply_dark_style(fig, axes_flat):
    """Apply consistent dark theme to all axes."""
    fig.patch.set_facecolor(DARK_BG)
    for ax in axes_flat:
        ax.set_facecolor(PANEL_BG)
        ax.tick_params(colors=TEXT_CLR, labelsize=9)
        ax.xaxis.label.set_color(TEXT_CLR)
        ax.yaxis.label.set_color(TEXT_CLR)
        ax.title.set_color(TEXT_CLR)
        for spine in ax.spines.values():
            spine.set_edgecolor(GRID_CLR)
        ax.grid(True, color=GRID_CLR, linewidth=0.5, linestyle="--", alpha=0.6)


def _draw_event_line(ax, t_event, label="O₂ exhausted"):
    """Draw a dashed vertical line at the infeasibility / depletion event."""
    if t_event is not None:
        ax.axvline(
            t_event,
            color=COLOURS["vline"],
            linewidth=1.2,
            linestyle="--",
            alpha=0.85,
            label=label,
        )


def _estimate_atp_yield(
    o2_arr, o2_0=OXYGEN_0, atp_aero=ATP_AEROBIC, atp_ana=ATP_ANAEROBIC
):
    """
    Linearly interpolate ATP yield between aerobic and anaerobic extremes
    based on fractional oxygen availability.

    ATP_yield = ATP_anaerobic + (ATP_aerobic - ATP_anaerobic) * (O2 / O2_0)

    This is a simplified heuristic for visualisation purposes.
    The real yield is computed internally by the FBA optimiser.
    """
    frac = np.clip(o2_arr / o2_0, 0.0, 1.0)
    return atp_ana + (atp_aero - atp_ana) * frac


def _compute_mu(t_arr, x_arr):
    """
    Estimate specific growth rate µ from the biomass time-series using
    a numerical finite-difference approximation:
        µ(t) ≈ (1/X) · dX/dt
    """
    dt = np.diff(t_arr)
    dX = np.diff(x_arr)
    X_m = 0.5 * (x_arr[:-1] + x_arr[1:])  # mid-point biomass
    t_m = 0.5 * (t_arr[:-1] + t_arr[1:])

    with np.errstate(divide="ignore", invalid="ignore"):
        mu = np.where(X_m > 1e-10, dX / (X_m * dt), 0.0)

    return t_m, mu


def plot_results(solution, save: bool = True, show: bool = True) -> str:
    """
    Generate the 2x3 multi-panel dFBA results figure.

    Parameters
    ----------
    solution : OdeResult
        Output from run_dfba().
    save : bool
        If True, save PNG to RESULTS_DIR.
    show : bool
        If True, display the figure interactively.

    Returns
    -------
    str
        Absolute path to the saved PNG (or empty string if save=False).
    """
    t = solution.t
    X = solution.y[0]
    Glc = solution.y[1]
    O2 = solution.y[2]
    Ace = solution.y[3]

    # Detect O2 transition event (index 1) and glucose exhaustion (index 0)
    t_o2_event = None
    t_glc_event = None
    if solution.t_events and len(solution.t_events) > 1:
        if solution.t_events[1].size > 0:
            t_o2_event = solution.t_events[1][0]
        if solution.t_events[0].size > 0:
            t_glc_event = solution.t_events[0][0]
    # Use O2 event as the primary visual transition line
    t_event = t_o2_event or t_glc_event

    # Derived quantities
    atp_yield = _estimate_atp_yield(O2)
    t_mu, mu = _compute_mu(t, X)

    # Figure layout
    fig = plt.figure(figsize=(16, 9), dpi=130)
    fig.patch.set_facecolor(DARK_BG)

    gs = GridSpec(
        2,
        3,
        figure=fig,
        hspace=0.45,
        wspace=0.35,
        left=0.07,
        right=0.97,
        top=0.88,
        bottom=0.1,
    )

    ax_bio = fig.add_subplot(gs[0, 0])
    ax_glc = fig.add_subplot(gs[0, 1])
    ax_o2 = fig.add_subplot(gs[0, 2])
    ax_ace = fig.add_subplot(gs[1, 0])
    ax_atp = fig.add_subplot(gs[1, 1])
    ax_mu = fig.add_subplot(gs[1, 2])

    axes_flat = [ax_bio, ax_glc, ax_o2, ax_ace, ax_atp, ax_mu]
    _apply_dark_style(fig, axes_flat)

    # Panel 1: Biomass
    ax_bio.plot(t, X, color=COLOURS["biomass"], linewidth=2.2, label="Biomass X")
    ax_bio.fill_between(t, 0, X, color=COLOURS["biomass"], alpha=0.12)
    _draw_event_line(ax_bio, t_event)
    ax_bio.set_title("① Biomass Growth", fontweight="bold", fontsize=10)
    ax_bio.set_xlabel("Time (h)")
    ax_bio.set_ylabel("g DW / L")
    ax_bio.legend(fontsize=8, facecolor=PANEL_BG, labelcolor=TEXT_CLR)

    # Panel 2: Glucose
    ax_glc.plot(t, Glc, color=COLOURS["glucose"], linewidth=2.2, label="Glucose")
    ax_glc.fill_between(t, 0, Glc, color=COLOURS["glucose"], alpha=0.12)
    _draw_event_line(ax_glc, t_event)
    ax_glc.set_title("② Glucose Depletion", fontweight="bold", fontsize=10)
    ax_glc.set_xlabel("Time (h)")
    ax_glc.set_ylabel("mmol / L")
    ax_glc.legend(fontsize=8, facecolor=PANEL_BG, labelcolor=TEXT_CLR)

    # Panel 3: Oxygen
    ax_o2.plot(t, O2, color=COLOURS["oxygen"], linewidth=2.2, label="O₂")
    ax_o2.fill_between(t, 0, O2, color=COLOURS["oxygen"], alpha=0.12)
    _draw_event_line(ax_o2, t_event)
    ax_o2.set_title("③ Oxygen Depletion", fontweight="bold", fontsize=10)
    ax_o2.set_xlabel("Time (h)")
    ax_o2.set_ylabel("mmol / L")
    ax_o2.legend(fontsize=8, facecolor=PANEL_BG, labelcolor=TEXT_CLR)

    # Panel 4: Acetate
    ax_ace.plot(
        t, Ace, color=COLOURS["acetate"], linewidth=2.2, label="Acetate (ferment.)"
    )
    ax_ace.fill_between(t, 0, Ace, color=COLOURS["acetate"], alpha=0.12)
    _draw_event_line(ax_ace, t_event)
    ax_ace.set_title("④ Acetate Secretion", fontweight="bold", fontsize=10)
    ax_ace.set_xlabel("Time (h)")
    ax_ace.set_ylabel("mmol / L")
    ax_ace.legend(fontsize=8, facecolor=PANEL_BG, labelcolor=TEXT_CLR)

    # Panel 5: ATP yield estimate
    ax_atp.plot(
        t, atp_yield, color=COLOURS["atp"], linewidth=2.2, label="Est. ATP yield"
    )
    ax_atp.axhline(
        ATP_AEROBIC,
        color=COLOURS["oxygen"],
        linewidth=1,
        linestyle=":",
        alpha=0.7,
        label=f"Aerobic max (~{int(ATP_AEROBIC)})",
    )
    ax_atp.axhline(
        ATP_ANAEROBIC,
        color=COLOURS["glucose"],
        linewidth=1,
        linestyle=":",
        alpha=0.7,
        label=f"Anaerobic min (~{int(ATP_ANAEROBIC)})",
    )
    _draw_event_line(ax_atp, t_event)
    ax_atp.set_title("⑤ Est. ATP Yield / Glucose", fontweight="bold", fontsize=10)
    ax_atp.set_xlabel("Time (h)")
    ax_atp.set_ylabel("ATP molecules / glucose")
    ax_atp.set_ylim(-1, ATP_AEROBIC + 3)
    ax_atp.legend(fontsize=7, facecolor=PANEL_BG, labelcolor=TEXT_CLR)

    # Panel 6: Specific growth rate
    ax_mu.plot(t_mu, mu, color=COLOURS["mu"], linewidth=2.0, label="µ (h⁻¹)")
    ax_mu.fill_between(t_mu, 0, np.clip(mu, 0, None), color=COLOURS["mu"], alpha=0.12)
    _draw_event_line(ax_mu, t_event)
    ax_mu.set_title("⑥ Specific Growth Rate µ", fontweight="bold", fontsize=10)
    ax_mu.set_xlabel("Time (h)")
    ax_mu.set_ylabel("h⁻¹")
    ax_mu.legend(fontsize=8, facecolor=PANEL_BG, labelcolor=TEXT_CLR)

    # Main title
    fig.suptitle(
        "dFBA Simulation: Aerobic→Anaerobic Respiratory Transition in E. coli Batch Culture",
        color=TEXT_CLR,
        fontsize=13,
        fontweight="bold",
        y=0.96,
    )

    # Annotation for transition zone
    if t_event is not None:
        fig.text(
            0.5,
            0.005,
            f"⚠  Substrate exhaustion at t ≈ {t_event:.2f} h  |  "
            "Dashed red line marks aerobic→anaerobic transition",
            ha="center",
            fontsize=8.5,
            color=COLOURS["vline"],
        )

    # Save
    out_path = ""
    if save:
        os.makedirs(RESULTS_DIR, exist_ok=True)
        out_path = os.path.join(RESULTS_DIR, "dfba_results.png")
        fig.savefig(
            out_path, dpi=130, bbox_inches="tight", facecolor=DARK_BG, edgecolor="none"
        )
        logger.info("Plot saved -> %s", os.path.abspath(out_path))

    if show:
        plt.show()

    plt.close(fig)
    return out_path
