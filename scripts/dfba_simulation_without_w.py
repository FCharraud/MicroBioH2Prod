#!/usr/bin/env python3
"""
Dynamic FBA Simulation - Without Tungsten

This script performs dynamic Flux Balance Analysis (dFBA) to simulate bacterial
growth without tungsten-containing enzymes. The tungstate-dependent formate
dehydrogenase (FDHTungs) is knocked out, and the tungsten coefficient in the
biomass objective function (BOF) is reduced to NEW_COEFF = 0.01, requiring
external formate as an electron donor.

The simulation includes:
- FDHTungs knockout and biomass tungsten content reduction (NEW_COEFF = 0.01)
- Formate-dependent growth
- Formate inhibition on FAE
- Dynamic metabolic constraints
- Maintenance energy requirements from methane

Author: Sylvain Davidson
Created: 2025-10-06
Last Modified: 2025-10-06
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
MODEL_INPUT_PATH = "../data/iIA409_Ca_W_Cu_formate_ex.xml"
OUTPUT_PATH = "../output/dFBA_iIA409_fdh_CaWCu.xlsx"


# ==============================================================================
# SIMULATION PARAMETERS
# ==============================================================================
DELTA_T = 0.01  # Time step (hours)
SIMULATION_TIME = 200  # Total simulation time (hours)
OUTPUT_INTERVAL = 20  # Data output every N iterations


# ==============================================================================
# METABOLIC PARAMETERS
# ==============================================================================
NGAM = 7.2  # ATP maintenance (mmol ATP/gDW/h)
ATP_YIELD_EFFICIENCY = 0.925  # ATP yield efficiency factor
NEW_COEFF = 0.01  # Tungsten biomass coefficient


# ==============================================================================
# MASS TRANSFER COEFFICIENTS (h⁻¹)
# ==============================================================================
KLA_CH4 = 45.3  # Methane mass transfer coefficient


# ==============================================================================
# EQUILIBRIUM CONCENTRATIONS (mM)
# ==============================================================================
EQ_O2 = 0.025  # Dissolved oxygen equilibrium
EQ_CH4 = 0.63  # Dissolved methane equilibrium


# ==============================================================================
# AFFINITY CONSTANTS - Monod kinetics (mM)
# ==============================================================================
KS_CH4 = 0.002  # 2 µM, from Hanson & Hanson 1996
KS_O2 = 0.002  # 2 µM, from van Bodegom 2001
KS_NO3 = 0.05  # 50 µM, organisms with nrtABC enzyme
KS_PI = 0  # mM
KS_SO4 = 0  # mM
KS_MG2 = 0  # mM
KS_CA2 = 0  # mM
KS_CU2 = 0  # mM
KS_FE2 = 0  # mM


# ==============================================================================
# INHIBITION CONSTANTS
# ==============================================================================
KI_FOR = 10  # Formate inhibition constant (mM)


# ==============================================================================
# MAXIMUM VOLUMETRIC RATES
# ==============================================================================
FVOL_O2 = 12  # Maximum volumetric O2 consumption (mM/h)


# ==============================================================================
# INITIAL CONDITIONS
# ==============================================================================
BIOMASS = 0.1  # g/L
C_CA2 = 0.088  # mM
C_MG2 = 4.15  # mM
C_PI = 27  # mM
C_SO4 = 4.16  # mM
C_FE2 = 0.08  # mM
C_CU2 = 0.09  # mM
C_NO3 = 10  # mM
C_FOR = 4  # mM
C_W = 0  # mM
C_CH4 = 0.63  # mM
C_O2 = 0.025  # mM
C_DIC = 0  # mM


def main():
    """
    Main execution function.

    Steps:
        1. Load model and run initial optimization to determine maximum fluxes.
        2. Reload and prepare model with FDHTungs knockout and simulation constraints.
        3. Initialize concentrations, medium and run first optimization (t=0).
        4. Initialize simulation variables and data arrays.
        5. Run the dFBA simulation loop over SIMULATION_TIME.
        6. Save results to Excel file.
    """
    print("=" * 80)
    print("DYNAMIC FBA SIMULATION - WITHOUT TUNGSTEN")
    print("=" * 80)
    print(f"Simulation time: {SIMULATION_TIME} hours")
    print(f"Time step: {DELTA_T} hours")
    print(f"Output interval: every {OUTPUT_INTERVAL} iterations")
    print("=" * 80)
    print()

    # =========================================================================
    # STEP 1: Load model and run initial optimization
    # =========================================================================
    print(f"Loading model: {MODEL_INPUT_PATH}")
    model_init = cobra.io.read_sbml_model(MODEL_INPUT_PATH)
    print(f"Model loaded successfully: {len(model_init.reactions)} reactions")
    print()

    # Remove tungsten from biomass equation
    print("Modifying biomass equation...")
    bof = model_init.reactions.get_by_id('BOF')
    w_metabolite = model_init.metabolites.get_by_id('tungs_c')
    bof.add_metabolites({w_metabolite: NEW_COEFF})
    print(f"Tungsten biomass coefficient set to: {NEW_COEFF}")
    print()

    # First run for max fluxes data
    print("Running initial optimization to determine maximum fluxes...")
    model_init.reactions.get_by_id('NGAM').lower_bound = NGAM
    model_init.objective = {model_init.reactions.get_by_id('BOF'): 1}

    # Medium definition
    medium_t0 = model_init.medium
    medium_t0["EX_o2_e"] = 25
    medium_t0["EX_for_e"] = 0
    model_init.medium = medium_t0

    solution_t0 = model_init.optimize()
    hps_max = solution_t0.fluxes["HPS"]
    fae_max = solution_t0.fluxes["FAE"]
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
    print(f"  FAE: {solution_t0.fluxes['FAE']}")

    # Max flux obtain from the first run
    fmax_o2 = -solution_t0.fluxes["EX_o2_e"]  # mmole/h/g
    fmax_ch4 = -solution_t0.fluxes["EX_ch4_e"]
    fmax_no3 = -solution_t0.fluxes["EX_no3_e"]
    fmax_pi = -solution_t0.fluxes["EX_pi_e"]
    fmax_so4 = -solution_t0.fluxes["EX_so4_e"]
    fmax_mg2 = -solution_t0.fluxes["EX_mg2_e"]
    fmax_ca2 = -solution_t0.fluxes["EX_ca2_e"]
    fmax_cu2 = -solution_t0.fluxes["EX_cu2_e"]
    fmax_fe2 = -solution_t0.fluxes["EX_fe2_e"]

    print(f"\nMaximum uptake fluxes (mmol/gDW/h):")
    print(f"  Fmax_O2: {fmax_o2:.4f}")
    print(f"  Fmax_CH4: {fmax_ch4:.4f}")
    print(f"  Fmax_NO3: {fmax_no3:.4f}")
    print()

    # =========================================================================
    # STEP 2: Prepare model for simulation loop
    # =========================================================================
    print("Preparing model for simulation loop...")
    # Re-import model for loop (like original - copy() doesn't work reliably)
    model = cobra.io.read_sbml_model(MODEL_INPUT_PATH)
    print(f"Model reloaded: {len(model.reactions)} reactions")

    # Remove tungsten from biomass
    bof = model.reactions.get_by_id('BOF')
    w_metabolite = model.metabolites.get_by_id('tungs_c')
    bof.add_metabolites({w_metabolite: NEW_COEFF})

    # Apply simulation-specific constraints
    model.reactions.get_by_id('FDHTungs').upper_bound = 0.00
    model.reactions.get_by_id('NGAM').lower_bound = NGAM
    model.reactions.get_by_id('NGAM').upper_bound = NGAM
    model.reactions.get_by_id('HPS').upper_bound = hps_max
    # Use model_init reference like in original
    model.objective = {model_init.reactions.get_by_id('BOF'): 1}
    print("FDHTungs knocked out, constraints applied")
    print()

    # =========================================================================
    # STEP 3: Initialize concentrations and medium
    # =========================================================================
    # Initial biomass and concentrations (using global constants)
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
    c_dic = C_DIC

    # Medium before first iteration
    medium = model.medium
    medium['EX_o2_e'] = fmax_o2
    medium['EX_ch4_e'] = fmax_ch4
    medium['EX_no3_e'] = fmax_no3 * C_NO3 / (KS_NO3 + C_NO3)
    medium['EX_pi_e'] = fmax_pi * C_PI / (KS_PI + C_PI)
    medium['EX_so4_e'] = fmax_so4 * C_SO4 / (KS_SO4 + C_SO4)
    medium['EX_for_e'] = C_FOR
    medium['EX_ca2_e'] = fmax_ca2 * C_CA2 / (KS_CA2 + C_CA2)
    medium['EX_mg2_e'] = fmax_mg2 * C_MG2 / (KS_MG2 + C_MG2)
    medium['EX_fe2_e'] = fmax_fe2 * C_FE2 / (KS_FE2 + C_FE2)
    medium['EX_cu2_e'] = fmax_cu2 * C_CU2 / (KS_CU2 + C_CU2)
    medium['EX_tungs_e'] = C_W
    model.medium = medium

    solution = model.optimize()

    # =========================================================================
    # STEP 4: Initialize simulation variables
    # =========================================================================
    total_iterations = int(SIMULATION_TIME / DELTA_T)
    cumul_ch4 = 0
    cumul_o2 = 0
    cumul_co2 = 0
    f_ch4 = fmax_ch4

    # Data catching frame
    tab_biomass = []
    tab_time = []
    tab_cumul_co2 = []
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
    tab_dic = []
    tab_fae = []
    tab_c_ch4 = []

    # First line of the growth data (t0)
    tab_biomass.append(biomass)
    tab_time.append(0)
    tab_cumul_co2.append(cumul_co2)
    tab_ex_ch4_e.append(abs(solution.fluxes['EX_ch4_e']))
    tab_c_ch4.append(c_ch4)
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
    tab_dic.append(c_dic)
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
    print("  formate: formate concentration (mM)")
    print("  ATP/CH4: ATP yield per CH4 consumed")
    print("=" * 80)
    print()

    for step in range(total_iterations):
        conso_o2 = solution.fluxes["EX_o2_e"]
        conso_ch4 = solution.fluxes["EX_ch4_e"]
        conso_no3 = solution.fluxes["EX_no3_e"]
        fae = solution.fluxes["FAE"]

        # Calcul of concentrations, consumption and production
        biomass += biomass * solution.fluxes['BOF'] * DELTA_T
        c_o2 = EQ_O2  # Disolved oxygen is regulated
        c_ch4 += (KLA_CH4 * (EQ_CH4 - c_ch4) - abs(conso_ch4) * biomass) * DELTA_T
        if c_ch4 < KS_CH4 / 1000:
            c_ch4 = KS_CH4 / 1000
        cumul_o2 += abs(conso_o2) * biomass * DELTA_T
        cumul_ch4 += biomass * abs(conso_ch4) * DELTA_T
        cumul_co2 += solution.fluxes['EX_CO2'] * biomass * DELTA_T
        c_no3 += conso_no3 * biomass * DELTA_T
        c_ca2 += solution.fluxes['EX_ca2_e'] * biomass * DELTA_T
        c_mg2 += solution.fluxes['EX_mg2_e'] * biomass * DELTA_T
        c_pi += solution.fluxes['EX_pi_e'] * biomass * DELTA_T
        c_so4 += solution.fluxes['EX_so4_e'] * biomass * DELTA_T
        c_fe2 += solution.fluxes['EX_fe2_e'] * biomass * DELTA_T
        c_cu2 += solution.fluxes['EX_cu2_e'] * biomass * DELTA_T
        c_w += solution.fluxes['EX_tungs_e'] * biomass * DELTA_T
        c_for += solution.fluxes['EX_for_e'] * biomass * DELTA_T

        # Formate retro-inhibition
        inhibition_formate = 1 / (1 + c_for / KI_FOR)

        # New boundaries
        ex_o2 = model.reactions.get_by_id('EX_o2_e')
        if abs(conso_o2) * biomass < FVOL_O2:
            f_o2 = (fmax_o2 * c_o2 / (KS_O2 + c_o2))
            ex_o2.lower_bound = -f_o2
        else:  # Volumic oxygen rate constrain
            ex_o2.lower_bound = -abs(FVOL_O2 / biomass)

        f_no3 = (fmax_no3 * c_no3 / (KS_NO3 + c_no3))
        f_ca2 = (fmax_ca2 * c_ca2) / (KS_CA2 + c_ca2)
        f_mg2 = (fmax_mg2 * c_mg2) / (KS_MG2 + c_mg2)
        f_pi = (fmax_pi * c_pi) / (KS_PI + c_pi)
        f_so4 = (fmax_so4 * c_so4) / (KS_SO4 + c_so4)
        f_fe2 = (fmax_fe2 * c_fe2) / (KS_FE2 + c_fe2)
        f_cu2 = (fmax_cu2 * c_cu2) / (KS_CU2 + c_cu2)

        ex_ch4 = model.reactions.get_by_id('EX_ch4_e')
        y_atp_ch4 = 3 + 3 * (1 - r_fae_ch4)  # ATP/CH4 yield for maintenance
        if abs(conso_ch4) < NGAM / (y_atp_ch4 * ATP_YIELD_EFFICIENCY ):
            # Giving enough methane to sustain maintenance
            ex_ch4.lower_bound = -NGAM / (y_atp_ch4 * ATP_YIELD_EFFICIENCY )
        else:
            # Keeping the same formaldehyde ventilation in absence of formate oxidation
            f_ch4 = fae / r_fae_ch4
            ex_ch4.lower_bound = -f_ch4
            model.reactions.get_by_id('FAE').upper_bound = fae_max * inhibition_formate
            model.reactions.get_by_id('HPS').lower_bound = -1
            model.reactions.get_by_id('HPS').upper_bound = f_ch4 * (1 - r_fae_ch4)

        # New medium
        medium = model.medium
        medium['EX_ch4_e'] = f_ch4
        medium['EX_no3_e'] = f_no3
        medium['EX_pi_e'] = f_pi
        medium['EX_so4_e'] = f_so4
        medium['EX_for_e'] = c_for
        medium['EX_ca2_e'] = f_ca2
        medium['EX_mg2_e'] = f_mg2
        medium['EX_fe2_e'] = f_fe2
        medium['EX_cu2_e'] = f_cu2
        medium['EX_tungs_e'] = c_w
        model.medium = medium

        # New optimization
        solution = model.optimize()

        # Execute model.summary() for better precision, catching data
        if (step + 1) % OUTPUT_INTERVAL == 0:
            model.summary()
            print("temps:", (step + 1) * DELTA_T, "μ:", solution.objective_value, "o2:", conso_o2,
                  "ch4:", conso_ch4, "no3:", conso_no3, c_no3,
                  "formate:", c_for, "ATP/CH4:", y_atp_ch4)

            tab_biomass.append(biomass)
            tab_time.append((step + 1) * DELTA_T)
            tab_cumul_co2.append(cumul_co2)
            tab_ex_ch4_e.append(abs(solution.fluxes['EX_ch4_e']))
            tab_c_ch4.append(c_ch4)
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
            tab_dic.append(c_dic)
            tab_fae.append(solution.fluxes['FAE'])

    # Creation of Excel file
    wb = Workbook()
    ws = wb.active
    ws.title = "dFBA Results"

    # Create a DataFrame for each data series
    data = {
        "Time (h)": tab_time,
        "Biomass (g)": tab_biomass,
        "Formate (mmol)": tab_ex_for,
        "CO2 (mmol)": tab_cumul_co2,
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
        "C_ch4 (mM)": tab_c_ch4,
        "fae (mmole/h/g)": tab_fae,
        "ch4 (mmole/h/g)": tab_ex_ch4_e
    }

    df = pd.DataFrame(data)

    # Adding the DataFrame to the Excel file
    for row in dataframe_to_rows(df, index=False, header=True):
        ws.append(row)

    # =========================================================================
    # STEP 6: Save final Excel file
    # =========================================================================
    print()
    print("=" * 80)
    print(f"Saving results to: {OUTPUT_PATH}")
    wb.save(OUTPUT_PATH)
    print("Excel file saved successfully!")
    print("=" * 80)
    print("Simulation completed!")


if __name__ == "__main__":
    main()