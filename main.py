"""
Aerobic-to-Anaerobic Respiratory Transition in E. coli
using Dynamic Flux Balance Analysis (dFBA)

Course   : IF3211 Domain-Specific Computation
Topic    : Computational Modelling of Cellular Respiration & ATP Dynamics
Method   : dFBA (Static Optimisation Approach) — COBRApy + SciPy BDF

Usage
-----
    python main.py [--no-show] [--no-save] [--t-end HOURS] [--glc MMOL]

Quick start:
    python main.py
"""

import argparse
import logging
import sys
import time

# Logging setup
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

# Imports
from src.model_loader import load_ecoli_model  # type: ignore
from src.solver import run_dfba  # type: ignore
from src.visualization import plot_results  # type: ignore
from src.config import Y0, T_START, T_END  # type: ignore


def parse_args():
    p = argparse.ArgumentParser(
        description="dFBA Batch Culture Simulation — E. coli Aerobic/Anaerobic Transition"
    )
    p.add_argument(
        "--no-show", action="store_true", help="Do not display the plot window"
    )
    p.add_argument("--no-save", action="store_true", help="Do not save the plot PNG")
    p.add_argument(
        "--t-end",
        type=float,
        default=T_END,
        metavar="HOURS",
        help=f"Simulation end time in hours (default: {T_END})",
    )
    p.add_argument(
        "--glc",
        type=float,
        default=Y0[1],
        metavar="MMOL",
        help=f"Initial glucose concentration mmol/L (default: {Y0[1]})",
    )
    p.add_argument(
        "--o2",
        type=float,
        default=Y0[2],
        metavar="MMOL",
        help=f"Initial oxygen concentration mmol/L (default: {Y0[2]})",
    )
    p.add_argument(
        "--biomass",
        type=float,
        default=Y0[0],
        metavar="G_DW_L",
        help=f"Initial biomass g_DW/L (default: {Y0[0]})",
    )
    return p.parse_args()


def main():
    args = parse_args()

    print()
    print("=" * 65)
    print("  dFBA Batch Culture Simulation")
    print("  E. coli Aerobic -> Anaerobic Respiratory Transition")
    print("=" * 65)
    print()

    # 1. Load metabolic model
    logger.info("Step 1/3 — Loading E. coli core metabolic model …")
    t0 = time.time()
    model = load_ecoli_model()
    logger.info("Model ready in %.2f s", time.time() - t0)

    # 2. Run dFBA
    logger.info("Step 2/3 — Running dFBA integration (BDF solver) …")
    y0_custom = [args.biomass, args.glc, args.o2, 0.0]
    t1 = time.time()
    solution = run_dfba(
        model=model,
        y0=y0_custom,
        t_span=(T_START, args.t_end),
    )
    logger.info("Integration finished in %.2f s", time.time() - t1)

    # 3. Visualise
    logger.info("Step 3/3 — Generating plots …")
    out_path = plot_results(
        solution,
        save=not args.no_save,
        show=not args.no_show,
    )

    print()
    print("=" * 65)
    print("  Simulation complete.")
    if out_path:
        print(f"  Plot saved -> {out_path}")
    print("=" * 65)
    print()

    # Quick console summary
    t_arr = solution.t
    X_arr = solution.y[0]
    O2_arr = solution.y[2]

    print(f"  Final biomass    : {X_arr[-1]:.4f} g_DW/L")
    print(f"  Final O2         : {O2_arr[-1]:.6f} mmol/L")
    print(f"  Time steps       : {len(t_arr)}")
    if solution.t_events and len(solution.t_events[0]) > 0:
        print(f"  Substrate exhaustion at t = {solution.t_events[0][0]:.3f} h")
    print()


if __name__ == "__main__":
    main()
