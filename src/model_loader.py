"""
Loads the E. coli 'textbook' (core) metabolic model and prepares it
for dynamic FBA by adjusting default exchange-reaction bounds.

The textbook model bundled with COBRApy is derived from the E. coli core
metabolic reconstruction (Orth et al., 2010) and contains 72 metabolites
and 95 reactions, including all central carbon metabolic pathways:
  - Glycolysis / Gluconeogenesis
  - Pentose Phosphate Pathway
  - TCA Cycle
  - Oxidative Phosphorylation
  - Fermentation branch points (acetate, ethanol, formate secretion)
"""

import logging
import cobra
from cobra.io import load_model

logger = logging.getLogger(__name__)


def load_ecoli_model() -> cobra.Model:
    """
    Load the E. coli core 'textbook' metabolic model from COBRApy's
    built-in data store.

    The model's default glucose and oxygen exchange bounds are relaxed so
    that dynamic bounds can be imposed entirely by add_dynamic_bounds().

    Returns
    -------
    cobra.Model
        Ready-to-use metabolic model with objective set to biomass
        production (BIOMASS_Ecoli_core_w_GAM).
    """
    logger.info("Loading E. coli core ('textbook') model via COBRApy …")
    model = load_model("textbook")

    # Relax glucose uptake bound — the dynamic system will cap it via
    # Michaelis-Menten kinetics at each integration step.
    model.reactions.get_by_id("EX_glc__D_e").lower_bound = -1000.0

    # Allow free oxygen uptake initially — will be capped dynamically.
    model.reactions.get_by_id("EX_o2_e").lower_bound = -1000.0

    # Allow acetate secretion (positive) and uptake (negative).
    model.reactions.get_by_id("EX_ac_e").lower_bound = -1000.0

    logger.info(
        "Model loaded: %d reactions, %d metabolites",
        len(model.reactions),
        len(model.metabolites),
    )
    return model
