#!/usr/bin/env python3
"""
Dynamic FBA Simulation - With Tungsten

This script performs dynamic Flux Balance Analysis (dFBA) to simulate bacterial
growth with tungsten-containing enzymes. The tungstate-dependent formate
dehydrogenase (FDHTungs) is active, and tungsten is supplied in the growth medium.

The simulation includes:
- Tungsten-dependent formate dehydrogenase (FDHTungs) active pathway
- Formate inhibition on FAE
- Dynamic metabolic constraints based on Monod kinetics
- Maintenance energy requirements from methane
- Volumetric oxygen transfer limitation

Author: Sylvain Davidson
Created: 2025-10-03
Last Modified: 2025-10-03
Version: 1.0.0 First release
License: MIT
"""

import cobra
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
import pandas as pd


# ==============================================================================
# FILE PATHS
# ==============================================================================
MODEL_INPUT_PATH = "../data/iIA409_Ca_W_Cu_formate_ex_modif_sucrose.xml"
OUTPUT_PATH = "../output/dFBA_iLA409_CaWCu.xlsx"


# ==============================================================================
# SIMULATION PARAMETERS
# ==============================================================================
DELTA_T = 0.01       # Time step (hours)
SIMULATION_TIME = 70  # Total simulation time (hours)
OUTPUT_INTERVAL = 20  # Data output every N iterations


# ==============================================================================
# METABOLIC PARAMETERS
# ==============================================================================
NGAM = 7.2  # Non-growth associated ATP maintenance (mmol ATP/gDW/h)


# ==============================================================================
# MASS TRANSFER COEFFICIENTS (h⁻¹)
# ==============================================================================
KLA_CH4 = 45.3  # Methane mass transfer coefficient
KLA_O2 = 53.7   # Oxygen mass transfer coefficient
KLA_CO2 = 51.3  # CO2 mass transfer coefficient


# ==============================================================================
# EQUILIBRIUM CONCENTRATIONS (mM)
# ==============================================================================
EQ_O2 = 0.025  # Dissolved oxygen equilibrium
EQ_CH4 = 0.63  # Dissolved methane equilibrium


# ==============================================================================
# AFFINITY CONSTANTS - Monod kinetics (mM)
# ==============================================================================
KS_CH4 = 0.002  # 2 µM, from Hanson & Hanson 1996
KS_O2 = 0.002   # 2 µM, from van Bodegom 2001
KS_NO3 = 0.05   # organisms with nrtABC enzyme: 30 µM cyanobacteria Anabaena (Swapnil 2021),
                # 1-10 µM for Synechocystis (Omata & al 1993),
                # 10-90 µM for Aspergillus (Akhtar 2015)
KS_PI = 0       # mM
KS_SO4 = 0      # mM
KS_MG2 = 0      # mM
KS_CA2 = 0      # mM
KS_CU2 = 0      # mM
KS_FE2 = 0      # mM
KS_W = 0        # mM


# ==============================================================================
# INHIBITION CONSTANTS
# ==============================================================================
KI_FOR = 10  # Formate inhibition constant (mM)


# ==============================================================================
# MAXIMUM VOLUMETRIC RATES
# ==============================================================================
FVOL_O2 = 12.5  # Maximum volumetric O2 consumption observed experimentally (mM/h)


# ==============================================================================
# INITIAL CONDITIONS
# ==============================================================================
BIOMASS = 0.1   # g/L
C_CA2 = 0.088   # mM
C_MG2 = 4.15    # mM
C_PI = 27       # mM
C_SO4 = 4.16    # mM
C_FE2 = 0.08    # mM
C_CU2 = 0.09    # mM
C_NO3 = 10      # mM
C_FOR = 0       # mM
C_W = 0.1       # mM
C_CH4 = 0.63    # mM
C_O2 = 0.025    # mM
C_CO2 = 0       # mM


def main():
    """
    Main execution function.

    Steps:
        1. Load model and run initial optimization to determine maximum fluxes.
        2. Reload and prepare model with simulation constraints.
        3. Initialize concentrations, medium and run first optimization (t=0).
        4. Initialize simulation variables and data arrays.
        5. Run the dFBA simulation loop over SIMULATION_TIME.
        6. Save results to Excel file.
    """
    print("=" * 80)
    print("DYNAMIC FBA SIMULATION - WITH TUNGSTEN")
    print("=" * 80)
    print(f"Simulation time: {SIMULATION_TIME} hours")
    print(f"Time step: {DELTA_T} hours")
    print(f"Output interval: every {OUTPUT_INTERVAL} iterations")
    print("=" * 80)
    print()

    # =========================================================================
    # STEP 1: Load model and run initial optimization to determine maximum fluxes
    # =========================================================================
    print(f"Loading model: {MODEL_INPUT_PATH}")
    model_init = cobra.io.read_sbml_model(MODEL_INPUT_PATH)
    print(f"Model loaded successfully: {len(model_init.reactions)} reactions")
    print()

    # Model objective and NGAM constraint
    model_init.reactions.get_by_id('NGAM').lower_bound = NGAM
    model_init.reactions.get_by_id('NGAM').upper_bound = NGAM
    model_init.objective = {model_init.reactions.get_by_id('BOF'): 1}

    # Medium for initial test run
    print("Running initial optimization to determine maximum fluxes...")
    medium_t0 = model_init.medium
    medium_t0["EX_o2_e"] = 25
    medium_t0["EX_for_e"] = 0.000
    model_init.medium = medium_t0

    solution_t0 = model_init.optimize()
    fae_max = solution_t0.fluxes["FAE"]
    hps_max = solution_t0.fluxes["HPS"]
    r_fae_ch4 = abs(solution_t0.fluxes["FAE"] / solution_t0.fluxes["EX_ch4_e"])

    print("\n--- Initial Optimization Results ---")
    print(f"R_fae_ch4 (FAE/CH4 ratio): {r_fae_ch4}")
    print(f"\nInitial medium composition (mmol/gDW/h):")
    print(medium_t0)
    print(f"\nOptimization solution (growth rate: {solution_t0.objective_value:.3f} h⁻¹):")
    print(solution_t0)
    print(f"\nKey fluxes (mmol/gDW/h):")
    print(f"  EX_ch4_e: {solution_t0.fluxes['EX_ch4_e']}")
    print(f"  BOF (growth): {solution_t0.fluxes['BOF']}")
    print(f"  EX_o2_e: {solution_t0.fluxes['EX_o2_e']}")
    print(f"  EX_no3_e: {solution_t0.fluxes['EX_no3_e']}")

    # Maximum uptake fluxes from the first run (mmol/gDW/h)
    fmax_o2 = -solution_t0.fluxes["EX_o2_e"]
    fmax_ch4 = -solution_t0.fluxes["EX_ch4_e"]
    fmax_no3 = -solution_t0.fluxes["EX_no3_e"]
    fmax_pi = -solution_t0.fluxes["EX_pi_e"]
    fmax_so4 = -solution_t0.fluxes["EX_so4_e"]
    fmax_mg2 = -solution_t0.fluxes["EX_mg2_e"]
    fmax_ca2 = -solution_t0.fluxes["EX_ca2_e"]
    fmax_cu2 = -solution_t0.fluxes["EX_cu2_e"]
    fmax_fe2 = -solution_t0.fluxes["EX_fe2_e"]
    fmax_w = -solution_t0.fluxes["EX_tungs_e"]

    print(f"\nMaximum uptake fluxes (mmol/gDW/h):")
    print(f"  Fmax_O2:  {fmax_o2:.4f}")
    print(f"  Fmax_CH4: {fmax_ch4:.4f}")
    print(f"  Fmax_NO3: {fmax_no3:.4f}")
    print()

    # =========================================================================
    # STEP 2: Prepare model for simulation loop
    # =========================================================================
    print("Preparing model for simulation loop...")
    # Re-import model for loop (copy() is not reliable with cobrapy)
    model = cobra.io.read_sbml_model(MODEL_INPUT_PATH)
    print(f"Model reloaded: {len(model.reactions)} reactions")

    # Apply simulation-specific constraints
    model.reactions.get_by_id('NGAM').lower_bound = NGAM
    model.reactions.get_by_id('NGAM').upper_bound = NGAM
    model.objective = {model_init.reactions.get_by_id('BOF'): 1}

    # =========================================================================
    # STEP 3: Initialize concentrations, medium and run first optimization (t=0)
    # =========================================================================
    # Mutable copies of initial conditions for the simulation loop
    biomass = BIOMASS
    c_ca2 = C_CA2
    c_mg2 = C_MG2
    c_pi = C_PI
    c_so4 = C_SO4
    c_fe2 = C_FE2
    c_cu2 = C_CU2
    c_no3 = C_NO3
    c_for = C_FOR
    c_w = C_W
    c_ch4 = C_CH4
    c_o2 = C_O2
    c_co2 = C_CO2

    # Medium at first iteration: O2 and CH4 unconstrained, nutrients Monod-limited
    medium = model.medium
    medium['EX_o2_e']    = fmax_o2
    medium['EX_ch4_e']   = fmax_ch4
    medium['EX_no3_e']   = fmax_no3 * c_no3 / (KS_NO3 + c_no3)
    medium['EX_pi_e']    = fmax_pi  * c_pi  / (KS_PI  + c_pi)
    medium['EX_so4_e']   = fmax_so4 * c_so4 / (KS_SO4 + c_so4)
    medium['EX_for_e']   = c_for
    medium['EX_ca2_e']   = fmax_ca2 * c_ca2 / (KS_CA2 + c_ca2)
    medium['EX_mg2_e']   = fmax_mg2 * c_mg2 / (KS_MG2 + c_mg2)
    medium['EX_fe2_e']   = fmax_fe2 * c_fe2 / (KS_FE2 + c_fe2)
    medium['EX_cu2_e']   = fmax_cu2 * c_cu2 / (KS_CU2 + c_cu2)
    medium['EX_tungs_e'] = fmax_w   * c_w   / (KS_W   + c_w)
    model.medium = medium

    solution = model.optimize()

    # =========================================================================
    # STEP 4: Initialize simulation variables and data arrays
    # =========================================================================
    total_iterations = int(SIMULATION_TIME / DELTA_T)
    cumul_ch4 = 0
    cumul_o2 = 0

    tab_biomass = []
    tab_time = []
    tab_co2 = []
    tab_ex_ch4_e = []
    tab_ex_o2_e = []
    tab_ex_ca2_e = []
    tab_ex_mg2_e = []
    tab_ex_pi_e = []
    tab_ex_so4_e = []
    tab_ex_fe2_e = []
    tab_ex_cu2_e = []
    tab_ex_no3_e = []
    tab_ex_w_e = []
    tab_ex_for = []
    tab_ex_biomass = []
    tab_cumul_ch4 = []
    tab_cumul_o2 = []
    tab_fae = []

    # First line of the growth data (t=0)
    tab_biomass.append(biomass)
    tab_time.append(0)
    tab_co2.append(c_co2)
    tab_ex_ch4_e.append(c_ch4)
    tab_ex_o2_e.append(c_o2)
    tab_ex_ca2_e.append(c_ca2)
    tab_ex_mg2_e.append(c_mg2)
    tab_ex_pi_e.append(c_pi)
    tab_ex_so4_e.append(c_so4)
    tab_ex_fe2_e.append(c_fe2)
    tab_ex_cu2_e.append(c_cu2)
    tab_ex_no3_e.append(c_no3)
    tab_ex_w_e.append(c_w)
    tab_ex_for.append(c_for)
    tab_ex_biomass.append(solution.fluxes['BOF'])
    tab_cumul_ch4.append(cumul_ch4)
    tab_cumul_o2.append(cumul_o2)
    tab_fae.append(solution.fluxes['FAE'])

    # =========================================================================
    # STEP 5: Main simulation loop
    # =========================================================================
    print(f"Starting simulation loop ({total_iterations} iterations)...")
    print("\n--- Simulation Progress ---")
    print("Column headers:")
    print("  temps: time (h)")
    print("  μ: growth rate (h⁻¹)")
    print("  o2: O2 consumption flux (mmol/gDW/h)")
    print("  ch4: CH4 consumption flux (mmol/gDW/h)")
    print("  no3: NO3 consumption flux (mmol/gDW/h) | NO3 concentration (mM)")
    print("=" * 80)
    print()

    for step in range(total_iterations):
        conso_o2 = solution.fluxes["EX_o2_e"]
        conso_ch4 = solution.fluxes["EX_ch4_e"]
        conso_no3 = solution.fluxes["EX_no3_e"]

        # Concentration updates: consumption and production
        biomass += biomass * solution.fluxes['BOF'] * DELTA_T
        c_o2 = EQ_O2  # Dissolved oxygen is regulated
        c_ch4 += (KLA_CH4 * (EQ_CH4 - c_ch4) - abs(conso_ch4) * biomass) * DELTA_T
        if c_ch4 < KS_CH4 / 1000:
            c_ch4 = KS_CH4 / 1000
        cumul_o2  += abs(conso_o2) * biomass * DELTA_T
        cumul_ch4 += biomass * abs(conso_ch4) * DELTA_T
        c_no3 += solution.fluxes['EX_no3_e'] * biomass * DELTA_T
        c_ca2 += solution.fluxes['EX_ca2_e'] * biomass * DELTA_T
        c_mg2 += solution.fluxes['EX_mg2_e'] * biomass * DELTA_T
        c_pi  += solution.fluxes['EX_pi_e']  * biomass * DELTA_T
        c_so4 += solution.fluxes['EX_so4_e'] * biomass * DELTA_T
        c_fe2 += solution.fluxes['EX_fe2_e'] * biomass * DELTA_T
        c_cu2 += solution.fluxes['EX_cu2_e'] * biomass * DELTA_T
        c_w   += solution.fluxes['EX_tungs_e'] * biomass * DELTA_T
        c_co2 += solution.fluxes['EX_CO2']   * biomass * DELTA_T
        c_for += solution.fluxes['EX_for_e'] * biomass * DELTA_T

        # Formate inhibition on FAE
        inhibition_formate = 1 / (1 + c_for / KI_FOR)

        # New boundaries: oxygen
        ex_o2 = model.reactions.get_by_id('EX_o2_e')
        if abs(conso_o2) * biomass < FVOL_O2:
            f_o2 = (fmax_o2 * c_o2 / (KS_O2 + c_o2))
            ex_o2.lower_bond = -f_o2
        else:  # Volumetric oxygen rate constraint
            ex_o2.lower_bound = -abs(FVOL_O2 / biomass)
            f_o2 = FVOL_O2 / biomass

        # New boundaries: nutrients (Monod kinetics) and new medium
        medium = model.medium
        medium['EX_no3_e']   = fmax_no3 * c_no3 / (KS_NO3 + c_no3)
        medium['EX_pi_e']    = fmax_pi  * c_pi  / (KS_PI  + c_pi)
        medium['EX_so4_e']   = fmax_so4 * c_so4 / (KS_SO4 + c_so4)
        medium['EX_for_e']   = c_for
        medium['EX_ca2_e']   = fmax_ca2 * c_ca2 / (KS_CA2 + c_ca2)
        medium['EX_mg2_e']   = fmax_mg2 * c_mg2 / (KS_MG2 + c_mg2)
        medium['EX_fe2_e']   = fmax_fe2 * c_fe2 / (KS_FE2 + c_fe2)
        medium['EX_cu2_e']   = fmax_cu2 * c_cu2 / (KS_CU2 + c_cu2)
        medium['EX_tungs_e'] = fmax_w   * c_w   / (KS_W   + c_w)
        model.medium = medium

        # Formaldehyde ventilation constraint (HPS upper bound linked to CH4 consumption)
        model.reactions.get_by_id('HPS').upper_bound = (1 - r_fae_ch4) * abs(conso_ch4)
        model.reactions.get_by_id('HPS').lower_bound = -0.01

        # New optimization
        solution = model.optimize()

        # Execute model.summary() for better precision, catching data
        if (step + 1) % OUTPUT_INTERVAL == 0:
            model.summary()
            print("temps:", step * DELTA_T, "μ:", solution.objective_value, "o2:", conso_o2,
                  "ch4:", conso_ch4, "no3:", conso_no3, c_no3)

            tab_biomass.append(biomass)
            tab_time.append(step * DELTA_T)
            tab_co2.append(c_co2)
            tab_ex_ch4_e.append(c_ch4)
            tab_ex_o2_e.append(c_o2)
            tab_ex_ca2_e.append(c_ca2)
            tab_ex_mg2_e.append(c_mg2)
            tab_ex_pi_e.append(c_pi)
            tab_ex_so4_e.append(c_so4)
            tab_ex_fe2_e.append(c_fe2)
            tab_ex_cu2_e.append(c_cu2)
            tab_ex_no3_e.append(c_no3)
            tab_ex_w_e.append(c_w)
            tab_ex_for.append(c_for)
            tab_ex_biomass.append(solution.fluxes['BOF'])
            tab_cumul_ch4.append(cumul_ch4)
            tab_cumul_o2.append(cumul_o2)
            tab_fae.append(solution.fluxes['FAE'])

    # Creation of Excel workbook
    wb = Workbook()
    ws = wb.active
    ws.title = "dFBA Results"

    # Update Excel file with all accumulated data
    df = pd.DataFrame({
        "Time (h)": tab_time,
        "Biomass (g)": tab_biomass,
        "Formate (mmol)": tab_ex_for,
        "CO2 (mmol)": tab_co2,
        "Methane (mmol)": tab_cumul_ch4,
        "Oxygen (mmol)": tab_cumul_o2,
        "Ca2+ (mmol)": tab_ex_ca2_e,
        "Mg (mmol)": tab_ex_mg2_e,
        "Phosphate (mmol)": tab_ex_pi_e,
        "Sulfate (mmol)": tab_ex_so4_e,
        "Fe3+ (mmol)": tab_ex_fe2_e,
        "Cu2+ (mmol)": tab_ex_cu2_e,
        "Nitrate (mmol)": tab_ex_no3_e,
        "Tungsten (mmol)": tab_ex_w_e,
        "µ (h-1)": tab_ex_biomass,
        "C_ch4 (mM)": tab_ex_ch4_e,
        "fae (mmole/h/g)": tab_fae
    })
    for row in dataframe_to_rows(df, index=False, header=True):
        ws.append(row)
    wb.save(OUTPUT_PATH)

    # =========================================================================
    # STEP 6: Simulation complete
    # =========================================================================
    print()
    print("=" * 80)
    print(f"Results saved to: {OUTPUT_PATH}")
    print("Simulation completed!")
    print("=" * 80)


if __name__ == "__main__":
    main()