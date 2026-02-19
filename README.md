# MicroBioH2Prod

**Studies on the metabolism of methanotrophs for efficient methane capture and utilisation**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue?logo=python&logoColor=white)](https://www.python.org/)
[![COBRApy](https://img.shields.io/badge/COBRApy-0.30.0-green)](https://cobrapy.readthedocs.io/)
[![Pandas](https://img.shields.io/badge/pandas-2.2.3-150458?logo=pandas&logoColor=white)](https://pandas.pydata.org/)
[![Organism](https://img.shields.io/badge/Organism-M._alcaliphilum_20Z-brightgreen)](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/568706/)
[![Model](https://img.shields.io/badge/GSM_Model-iIA409-orange)](https://gitlab.sirius-web.org/RSF/20ZR_CS_GSM_model)

## Overview

This repository contains metabolic modeling scripts for analyzing methanotrophic bacteria, specifically focusing on *Methylomicrobium alcaliphilum* 20Z. The project uses constraint-based metabolic modeling (Flux Balance Analysis - FBA) and dynamic FBA (dFBA) to study:

- Formate metabolism and transport
- Effects of formate dehydrogenase (FDH) knockout
- Growth kinetics with and without tungsten
- Metabolic flux distributions under various constraints

The models are based on the genome-scale metabolic model iIA409, extended with calcium, tungsten, and copper metabolism.

## Project Structure

```
MicroBioH2Prod/
├── README.md                                      # This file
├── LICENSE                                        # MIT License
├── requirements.txt                               # Python dependencies
├── data/                                          # SBML model files
│   ├── iIA409_Ca_W_Cu.xml                         # Base metabolic model downloaded by model_formate_modification.py
│   ├── iIA409_Ca_W_Cu_formate_ex.xml              # model created from iIA409_Ca_W_Cu.xml by model_formate_modification.py
│   └── iIA409_Ca_W_Cu_formate_ex_modif_sucrose.xml
├── scripts/                                       # Analysis scripts
│   ├── model_formate_modification.py              # Add formate exchange to model
│   ├── flux_analysis_fdh_constraint.py            # FBA with FDH constraints
│   ├── dfba_simulation_with_w.py                  # Dynamic FBA with tungsten
│   └── dfba_simulation_without_w.py               # Dynamic FBA without tungsten
└── output/                                        # Generated results (Excel files)
    ├── ilA409_W_Ca_Cu_flux_ngam7.2.xlsx           # created by flux_analysis_fdh_constraint.py
    ├── dFBA_iLA409_CaWCu.xlsx                     # created by dfba_simulation_with_w.py
    └── dFBA_iIA409_fdh_CaWCu.xlsx                 # created by dfba_simulation_without_w.py
```

## Installation

### Prerequisites

- Python 3.9
- pip package manager

### Dependencies

Install required packages:

```bash
pip install -r requirements.txt
```

Main dependencies:
- COBRApy (Constraint-Based Reconstruction and Analysis)
- pandas
- openpyxl

### Model Download

The base metabolic model `iIA409_Ca_W_Cu.xml` is **automatically downloaded** by the `model_formate_modification.py` script from:
https://gitlab.sirius-web.org/RSF/20ZR_CS_GSM_model

The script will:
- Create the `data/` directory if it doesn't exist
- Download the latest version of the model from the master branch
- Overwrite any existing model file

**Manual download (optional):** If you prefer to download manually, place the model file in the `data/` directory before running the scripts.

## Usage

### 1. Model Modification - Add Formate Exchange

This script downloads the base model and modifies it to include formate import/export capabilities.

```bash
python scripts/model_formate_modification.py
```

**Actions performed:**
- Downloads `iIA409_Ca_W_Cu.xml` from GitLab (creates `data/` directory if needed)
- Adds extracellular formate metabolite (`for_e`)
- Creates formate transport reaction
- Adds formate exchange boundary
- Adjusts biomass formate content (coefficient: 0.69)
- Sets default bounds for NGAM, CH4 exchange, and formate exchange

**Inputs:**
- None (downloads automatically from GitLab)

**Outputs:**
- `data/iIA409_Ca_W_Cu.xml` - Downloaded base model
- `data/iIA409_Ca_W_Cu_formate_ex.xml` - Modified model with formate exchange

### 2. Flux Analysis with FDH Constraints

Performs FBA under three scenarios:
1. Wild-type (with formate dehydrogenase)
2. FDH knockout (without formate dehydrogenase)
3. Double constraint (FDH + HPS limitations)

```bash
python scripts/flux_analysis_fdh_constraint.py
```

**Inputs:**
- `data/iIA409_Ca_W_Cu_formate_ex.xml` - Model with formate exchange

**Outputs:**
- `output/ilA409_W_Ca_Cu_flux_ngam7.2.xlsx` - Excel file with three sheets:
  - `with Fdh`: Wild-type flux distribution
  - `without Fdh`: FDH knockout flux distribution
  - `dble constrain`: Double constraint flux distribution

**Parameters:**
- NGAM (maintenance): 7.2 mmol ATP/gDW/h
- O2 uptake: 25 mmol/gDW/h
- NO3 uptake: 2 mmol/gDW/h
- Formate: 0 mmol/gDW/h

### 3. Dynamic FBA - With Tungsten

Simulates bacterial growth over time with tungsten-containing enzymes (FDH active).

```bash
python scripts/dfba_simulation_with_w.py
```

**Inputs:**
- `data/iIA409_Ca_W_Cu_formate_ex_modif_sucrose.xml` - Modified model

**Outputs:**
- `output/dFBA_iLA409_CaWCu.xlsx` - Time-series data including:
  - Biomass concentration
  - Substrate consumption (CH4, O2, NO3, etc.)
  - Product formation (CO2, formate)
  - Growth rate (μ)
  - Metabolic fluxes (FAE)

**Simulation parameters:**
- Time: 70 hours
- Time step (ΔT): 0.01 hours
- Initial biomass: 0.1 g/L
- NGAM: 7.2 mmol ATP/gDW/h
- Mass transfer coefficients (Kla):
  - CH4: 45.3 h⁻¹
  - O2: 53.7 h⁻¹
  - CO2: 51.3 h⁻¹
- Affinity constants (Ks):
  - CH4: 0.002 mM
  - O2: 0.002 mM
  - NO3: 0.05 mM
- Formate inhibition constant (Ki): 10 mM

### 4. Dynamic FBA - Without Tungsten

Simulates growth without tungsten, requiring formate as external electron donor.

```bash
python scripts/dfba_simulation_without_w.py
```

**Inputs:**
- `data/iIA409_Ca_W_Cu_formate_ex.xml` - Model with formate exchange

**Outputs:**
- `output/dFBA_iIA409_fdh_CaWCu.xlsx` - Time-series growth data

**Key differences from tungsten simulation:**
- FDH is knocked out (upper bound = 0)
- Initial formate: 4 mM
- Tungsten biomass coefficient reduced to 0.01
- Simulation time: 200 hours
- Formate inhibition on FAE included

## Metabolic Model Details

### Base Model
- **Organism**: *Methylomicrobium alcaliphilum* 20Z
- **Model ID**: iIA409
- **Extensions**: Calcium, Tungsten, Copper metabolism
- **Reactions**: ~409 reactions
- **Metabolites**: Includes formate, methane, oxygen, nitrate, and trace elements

### Key Reactions
- **MMO** (Methane Monooxygenase): CH4 → CH3OH
- **FAE** (Formaldehyde-activating enzyme): CH2O → CH2=THF
- **HPS** (Hexulose-6-phosphate synthase): CH2O + Ru5P → H6P
- **FDH** (Formate dehydrogenase, tungsten-dependent): Formate → CO2
- **NGAM** (Non-growth associated maintenance): ATP hydrolysis

## Output Files

All simulation results are saved as Excel files in the `output/` directory with:
- Time-series data for all metabolites
- Flux distributions
- Biomass growth curves
- Substrate/product concentrations

## Scientific Background

Methanotrophic bacteria convert methane (a potent greenhouse gas) into useful products and biomass. This work investigates:

1. **Formate metabolism**: Role of formate as electron donor and carbon source
2. **Tungsten dependency**: Impact of tungsten-containing formate dehydrogenase
3. **Metabolic engineering**: Constraining key pathways (HPS, FAE) to redirect carbon flux
4. **Growth kinetics**: Dynamic modeling of substrate-limited growth

## References

Original model repository:
- https://gitlab.sirius-web.org/RSF/20ZR_CS_GSM_model

COBRApy documentation:
- https://cobrapy.readthedocs.io/

## Author

**Sylvain Davidson**

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions, issues, and feature requests are welcome!

## Version History

- **v1.0.0** (2025-02-16): Initial release
  - Model formate modification script
  - Flux analysis with FDH constraints
  - Dynamic FBA simulations (with/without tungsten)

---

**Note**: This is a research project. Results should be validated experimentally before drawing biological conclusions.
