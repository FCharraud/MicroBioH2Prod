#!/usr/bin/env python3
"""
Flux Analysis with FDH Constraints

This script performs Flux Balance Analysis (FBA) on the iIA409 model under
three different scenarios:
1. Wild-type (with formate dehydrogenase)
2. FDH knockout (without formate dehydrogenase)
3. Double constraint (FDH + HPS limitations)

The script compares metabolic flux distributions and exports results to Excel.

Author: Sylvain Davidson
Created: 2024-12-19
Last Modified: 2024-12-19
Version: 1.0.0
License: MIT
"""

import cobra
import pandas as pd
from pathlib import Path


# Constants - File paths
MODEL_INPUT_PATH = Path("../data/iIA409_Ca_W_Cu_formate_ex.xml")
OUTPUT_PATH = Path("../output/ilA409_W_Ca_Cu_flux_ngam7.2.xlsx")

# Constants - Metabolic parameters
NGAM_MAINTENANCE = 7.2  # mmol ATP/gDW/h - Non-growth associated maintenance
O2_UPTAKE = 25.0  # mmol/gDW/h
NO3_UPTAKE = 2.0  # mmol/gDW/h
FORMATE_UPTAKE = 0.0  # mmol/gDW/h
METHANOL_UPTAKE = 0.0  # mmol/gDW/h

# Constants - Reaction IDs
BIOMASS_ID = "BOF"
NGAM_ID = "NGAM"
FDH_ID = "FDHTungs"  # Formate dehydrogenase (tungsten-dependent)
HPS_ID = "HPS"  # Hexulose-6-phosphate synthase
FAE_ID = "FAE"  # Formaldehyde-activating enzyme

# Constants - Exchange reaction IDs
EX_O2_ID = "EX_o2_e"
EX_NO3_ID = "EX_no3_e"
EX_FORMATE_ID = "EX_for_e"
EX_METHANOL_ID = "EX_meoh_e"
EX_CH4_ID = "EX_ch4_e"

# Constants - Constraint values
HPS_LOWER_BOUND_RELAXED = -0.01
FAE_CONSTRAINT_FACTOR = 1.0  # Used in double constraint scenario


def load_model(model_path):
    """
    Load SBML metabolic model from file.

    Parameters
    ----------
    model_path : Path or str
        Path to the SBML model file

    Returns
    -------
    cobra.Model
        Loaded metabolic model

    Raises
    ------
    FileNotFoundError
        If the model file does not exist
    """
    if not Path(model_path).exists():
        raise FileNotFoundError(
            f"Model file not found: {model_path}\n"
            f"Ensure you have run model_formate_modification.py first."
        )
    
    print(f"Loading model: {model_path}")
    model = cobra.io.read_sbml_model(str(model_path))
    print(f"Model loaded: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites\n")
    return model


def setup_medium(model, ch4=None, o2=O2_UPTAKE, no3=NO3_UPTAKE, 
                 formate=FORMATE_UPTAKE, methanol=METHANOL_UPTAKE):
    """
    Configure growth medium composition.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model
    ch4 : float, optional
        Methane uptake rate (if None, uses model default)
    o2 : float
        Oxygen uptake rate
    no3 : float
        Nitrate uptake rate
    formate : float
        Formate uptake rate
    methanol : float
        Methanol uptake rate

    Returns
    -------
    dict
        Medium composition
    """
    medium = model.medium
    
    if ch4 is not None:
        medium[EX_CH4_ID] = ch4
    medium[EX_O2_ID] = o2
    medium[EX_NO3_ID] = no3
    medium[EX_FORMATE_ID] = formate
    medium[EX_METHANOL_ID] = methanol
    
    model.medium = medium
    return medium


def extract_fluxes(model):
    """
    Extract all reaction fluxes from optimized model.

    Parameters
    ----------
    model : cobra.Model
        Optimized metabolic model

    Returns
    -------
    pd.DataFrame
        DataFrame with reaction ID, definition, flux, and absolute flux,
        sorted by absolute flux (descending)
    """
    flux_data = []
    
    for reaction in model.reactions:
        # Call summary() to ensure consistent behavior with original script
        summary = reaction.summary()
        flux = reaction.flux
        flux_data.append({
            "Reaction": reaction.id,
            "Definition": reaction.reaction,
            "Flux": flux,
            "Absolute Flux": abs(flux)
        })
    
    flux_df = pd.DataFrame(flux_data)
    flux_df_sorted = flux_df.sort_values(by="Absolute Flux", ascending=False)
    
    return flux_df_sorted


def run_wildtype_analysis(model):
    """
    Run FBA analysis on wild-type model (no constraints).

    Parameters
    ----------
    model : cobra.Model
        Metabolic model

    Returns
    -------
    tuple
        (solution, flux_dataframe, hps_max, fae_max, ch4_consumption, fae_ch4_ratio)
    
    Raises
    ------
    ValueError
        If optimization fails
    """
    print("=" * 70)
    print("WILD-TYPE ANALYSIS (with FDH)")
    print("=" * 70)
    
    # Setup
    model.reactions.get_by_id(NGAM_ID).lower_bound = NGAM_MAINTENANCE
    model.objective = model.reactions.get_by_id(BIOMASS_ID)
    
    medium = setup_medium(model)
    
    # Optimize
    solution = model.optimize()
    
    if solution.status != "optimal":
        raise ValueError(f"Optimization failed. Status: {solution.status}")
    
    # Extract key fluxes
    hps_max = solution.fluxes[HPS_ID]
    fae_max = solution.fluxes[FAE_ID]
    ch4_consumption = solution.fluxes[EX_CH4_ID]
    fae_ch4_ratio = abs(solution.fluxes[FAE_ID] / solution.fluxes[EX_CH4_ID])
    
    print(f"Optimization result: {solution.objective_value:.4f} h⁻¹")
    print(f"Medium: {medium}")
    print(model.summary())
    print(f"\nKey fluxes:")
    print(f"  HPS: {hps_max:.4f}")
    print(f"  FAE: {fae_max:.4f}")
    print(f"  CH4 consumption: {ch4_consumption:.4f}")
    print(f"  FAE/CH4 ratio: {fae_ch4_ratio:.4f}\n")
    
    flux_df = extract_fluxes(model)
    
    return solution, flux_df, hps_max, fae_max, ch4_consumption, fae_ch4_ratio


def run_fdh_knockout_analysis(model):
    """
    Run FBA analysis with FDH knockout.

    Simulates the effect of removing formate dehydrogenase, preventing
    formate oxidation to CO2.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model (will be copied)

    Returns
    -------
    tuple
        (solution, flux_dataframe)
    """
    print("=" * 70)
    print("FDH KNOCKOUT ANALYSIS (without FDH)")
    print("=" * 70)
    
    # Create copy to avoid modifying original
    fdh_model = model.copy()
    
    # Setup
    fdh_model.reactions.get_by_id(NGAM_ID).lower_bound = NGAM_MAINTENANCE
    
    # Knockout FDH
    fdh_model.reactions.get_by_id(FDH_ID).knock_out()
    print(f"Knocked out reaction: {FDH_ID}")
    
    medium = setup_medium(fdh_model)
    
    # Optimize
    solution = fdh_model.optimize()
    
    print(f"Optimization result: {solution.objective_value:.4f} h⁻¹")
    print(f"Medium: {medium}\n")
    
    flux_df = extract_fluxes(fdh_model)
    
    return solution, flux_df


def run_double_constraint_analysis(model, hps_max, fae_max, ch4_consumption, fae_ch4_ratio):
    """
    Run FBA analysis with double constraints (FDH + HPS).

    Applies constraints on both FDH and HPS to study metabolic redistribution.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model (will be copied)
    hps_max : float
        Maximum HPS flux from wild-type
    fae_max : float
        Maximum FAE flux from wild-type
    ch4_consumption : float
        CH4 consumption rate from wild-type
    fae_ch4_ratio : float
        FAE/CH4 ratio from wild-type

    Returns
    -------
    tuple
        (solution, flux_dataframe)
    """
    print("=" * 70)
    print("DOUBLE CONSTRAINT ANALYSIS (FDH + HPS)")
    print("=" * 70)
    
    # Create copy
    constrained_model = model.copy()
    
    # Setup maintenance
    constrained_model.reactions.get_by_id(NGAM_ID).lower_bound = NGAM_MAINTENANCE
    constrained_model.reactions.get_by_id(NGAM_ID).upper_bound = NGAM_MAINTENANCE
    
    # Apply HPS constraint based on FAE/CH4 ratio
    hps_upper = (1 - fae_ch4_ratio) * abs(ch4_consumption)
    constrained_model.reactions.get_by_id(HPS_ID).upper_bound = hps_upper
    constrained_model.reactions.get_by_id(HPS_ID).lower_bound = HPS_LOWER_BOUND_RELAXED
    
    # Apply FAE constraint
    constrained_model.reactions.get_by_id(FAE_ID).lower_bound = fae_max
    
    # Fix CH4 consumption
    constrained_model.reactions.get_by_id(EX_CH4_ID).lower_bound = ch4_consumption
    
    print(f"Applied constraints:")
    print(f"  NGAM: {NGAM_MAINTENANCE} (fixed)")
    print(f"  HPS upper bound: {hps_upper:.4f}")
    print(f"  HPS lower bound: {HPS_LOWER_BOUND_RELAXED}")
    print(f"  FAE lower bound: {fae_max:.4f}")
    print(f"  CH4 consumption: {ch4_consumption:.4f} (fixed)")
    
    medium = setup_medium(constrained_model)
    
    # Optimize
    solution = constrained_model.optimize()
    
    print(f"\nOptimization result: {solution.objective_value:.4f} h⁻¹")
    print(f"Medium: {medium}")
    print(constrained_model.summary())
    print()
    
    flux_df = extract_fluxes(constrained_model)
    
    return solution, flux_df


def save_results(wild_fluxes, fdh_fluxes, double_fluxes, output_path):
    """
    Save flux analysis results to Excel file.

    Creates an Excel file with three sheets:
    - 'with Fdh': Wild-type flux distribution
    - 'without Fdh': FDH knockout flux distribution
    - 'dble constrain': Double constraint flux distribution

    Parameters
    ----------
    wild_fluxes : pd.DataFrame
        Wild-type flux data
    fdh_fluxes : pd.DataFrame
        FDH knockout flux data
    double_fluxes : pd.DataFrame
        Double constraint flux data
    output_path : Path or str
        Output Excel file path
    """
    # Ensure output directory exists
    output_dir = Path(output_path).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print(f"Saving results to: {output_path}")
    
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        wild_fluxes.to_excel(writer, sheet_name="with Fdh", index=False)
        fdh_fluxes.to_excel(writer, sheet_name="without Fdh", index=False)
        double_fluxes.to_excel(writer, sheet_name="dble constrain", index=False)
    
    print("Results saved successfully!")
    print("=" * 70)


def main():
    """
    Main execution function.

    Workflow:
    1. Load model with formate exchange
    2. Run wild-type FBA analysis
    3. Run FDH knockout analysis
    4. Run double constraint analysis
    5. Save all results to Excel
    """
    print("\n" + "=" * 70)
    print("FLUX ANALYSIS WITH FDH CONSTRAINTS")
    print("=" * 70 + "\n")
    
    try:
        # Load model
        model = load_model(MODEL_INPUT_PATH)
        
        # Run analyses
        (wild_solution, wild_fluxes, hps_max, fae_max, ch4_cons, fae_ch4_ratio) = run_wildtype_analysis(model)
        
        fdh_solution, fdh_fluxes = run_fdh_knockout_analysis(model)
        
        double_solution, double_fluxes = run_double_constraint_analysis(
            model, hps_max, fae_max, ch4_cons, fae_ch4_ratio
        )
        
        # Save results
        save_results(wild_fluxes, fdh_fluxes, double_fluxes, OUTPUT_PATH)
        
        print("\nAnalysis completed successfully!\n")
        
    except Exception as e:
        print(f"\nError: {e}")
        raise


if __name__ == "__main__":
    main()
