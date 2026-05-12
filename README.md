# dFBA Batch Culture Simulation

> **IF3211 Domain-Specific Computation**
> Computational Modelling of Aerobic-to-Anaerobic Respiratory Transition in _E. coli_ using Dynamic Flux Balance Analysis (dFBA)

---

## Project Overview

This project simulates how _Escherichia coli_ shifts from efficient **aerobic respiration** (~30 ATP/glucose) to low-yield **anaerobic fermentation** (~2 ATP/glucose) as oxygen is progressively depleted in a sealed batch culture.

The simulation uses a **hybrid dFBA approach**:

- **COBRApy** — loads the _E. coli_ core genome-scale metabolic model and solves the FBA linear programme at each time step
- **SciPy `solve_ivp` (BDF)** — drives the outer ODE loop that tracks biomass, glucose, oxygen, and acetate concentrations over time

---

## Project Structure

```
dFBA-Batch-Culture-Simulation/
├── main.py                  # CLI entry point
├── pyproject.toml           # Dependencies
├── src/
│   ├── config.py            # All simulation parameters
│   ├── model_loader.py      # COBRApy textbook model loader
│   ├── dynamic_bounds.py    # Michaelis-Menten kinetic bounds
│   ├── dfba_system.py       # ODE RHS + infeasibility event
│   ├── solver.py            # solve_ivp BDF wrapper
│   └── visualization.py    # Multi-panel Matplotlib plots
├── notebooks/
│   └── dfba_simulation.ipynb  # Annotated Jupyter Notebook
├── results/                 # Auto-created; PNG plots saved here
└── docs/
    └── Proyek Respirasi Sel & ATP_ Judul & Rencana.txt
```

---

## Installation

```bash
# 1. Create and activate virtual environment
python -m venv .venv
.venv\Scripts\activate        # Windows
# source .venv/bin/activate   # Linux/macOS

# 2. Install dependencies
pip install cobra scipy numpy matplotlib

# Optional: for the Jupyter Notebook
pip install jupyter ipykernel
```

---

## Running the Simulation

### CLI (Recommended)

```bash
python main.py
```

### With custom parameters

```bash
python main.py --t-end 15 --glc 20 --o2 0.3

# Options:
#   --t-end HOURS     Simulation end time (default: 20.0 h)
#   --glc   MMOL      Initial glucose mmol/L (default: 10.0)
#   --o2    MMOL      Initial oxygen mmol/L  (default: 0.21)
#   --biomass G_DW_L  Initial biomass g_DW/L (default: 0.01)
#   --no-show         Skip interactive plot window
#   --no-save         Skip saving PNG to results/
```

### Jupyter Notebook

```bash
jupyter notebook notebooks/dfba_simulation.ipynb
```

---

## Output

Running `main.py` produces a **6-panel figure** saved to `results/dfba_results.png`:

| Panel | Variable            | Biological Significance               |
| ----- | ------------------- | ------------------------------------- |
| 1     | Biomass (g_DW/L)    | Logistic-like growth curve            |
| 2     | Glucose (mmol/L)    | Primary carbon source depletion       |
| 3     | Oxygen (mmol/L)     | Terminal electron acceptor exhaustion |
| 4     | Acetate (mmol/L)    | Fermentation overflow metabolite      |
| 5     | Est. ATP yield      | Aerobic→anaerobic efficiency drop     |
| 6     | Growth rate µ (h⁻¹) | Real-time metabolic activity          |

A **red dashed line** marks the moment of substrate exhaustion.

---

## Algorithm: dFBA Static Optimisation Approach (SOA)

```
for each time step t:
  1. Read current [X, Glc, O2, Ace] from ODE state
  2. Compute Michaelis-Menten uptake rates → set exchange bounds
  3. Solve FBA linear programme (COBRApy) → get µ, v_glc, v_o2, v_ace
  4. Compute dX/dt, dGlc/dt, dO2/dt, dAce/dt
  5. If FBA infeasible → trigger terminal event → stop
```

---

## Key Parameters (`src/config.py`)

| Parameter  | Value           | Description                      |
| ---------- | --------------- | -------------------------------- |
| `VMAX_GLC` | 10.0 mmol/(g·h) | Max glucose uptake rate          |
| `KM_GLC`   | 0.015 mmol/L    | Glucose half-saturation constant |
| `VMAX_O2`  | 15.0 mmol/(g·h) | Max oxygen uptake rate           |
| `KM_O2`    | 0.001 mmol/L    | Oxygen half-saturation constant  |
| `METHOD`   | BDF             | ODE integration method (stiff)   |
| `T_END`    | 20.0 h          | Default simulation duration      |

---

## Research Question

> _"How does progressive oxygen depletion in a batch culture dynamically regulate the stoichiometric flux boundaries of E. coli metabolism, and what is the implication on the total ATP acquisition efficiency curve versus fermentation byproduct excretion?"_

---

## Team Roles

| Role                          | Responsibility                                                              |
| ----------------------------- | --------------------------------------------------------------------------- |
| **Biology Researcher**        | Literature review, parameter validation, Introduction & Conclusion chapters |
| **Python Architect**          | `dfba_system.py`, `solver.py`, `dynamic_bounds.py`, repository management   |
| **Data Analyst / Visualiser** | `visualization.py`, Results & Discussion chapter, plot interpretation       |
| **Multimedia / Presentation** | Report formatting, slide deck, video recording & editing                    |

---

## References

1. COBRApy dFBA Tutorial — https://cobrapy.readthedocs.io/en/latest/dfba.html
2. Orth, J.D., Thiele, I., Palsson, B.Ø. (2010). _What is flux balance analysis?_ Nature Biotechnology, 28, 245–248.
3. Mahadevan, R., Edwards, J.S., Doyle, F.J. (2002). _Dynamic flux balance analysis of diauxic growth in Escherichia coli._ Biophysical Journal, 83(3), 1331–1340.
