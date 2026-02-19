#!/usr/bin/env python3
"""
Model Formate Modification Script

This script modifies the iIA409 metabolic model to include formate import/export
capabilities. It adds extracellular formate metabolite, transport reaction, and
exchange boundary, along with adjustments to biomass composition.

Author: Sylvain Davidson
Created: 2025-10-07
Last Modified: 2025-10-07
Version: 1.0.0
License: MIT
"""

import cobra
from pathlib import Path
import urllib.request
import urllib.error


# Constants
MODEL_URL = "https://gitlab.sirius-web.org/RSF/20ZR_CS_GSM_model/-/raw/master/Data/CS_GSM_model/file_collection.files/models/iIA409_Ca_W_Cu.xml"
MODEL_INPUT_PATH = Path("../data/iIA409_Ca_W_Cu.xml")
MODEL_OUTPUT_PATH = Path("../data/iIA409_Ca_W_Cu_formate_ex.xml")

# Metabolite IDs
FORMATE_EXTRACELLULAR = "for_e"
FORMATE_CYTOPLASMIC = "for_c"
METHANE_EXCHANGE = "ch4_e"

# Reaction IDs
FORMATE_TRANSPORT_ID = "for_transport"
FORMATE_EXCHANGE_ID = "EX_FORMATE"
NGAM_ID = "NGAM"
METHANE_EXCHANGE_ID = "EX_ch4_e"
BIOMASS_ID = "BOF"

# Default bounds
TRANSPORT_LOWER_BOUND = -1000.0
TRANSPORT_UPPER_BOUND = 1000.0
EXCHANGE_LOWER_BOUND = 0.0
EXCHANGE_UPPER_BOUND = 0.0
NGAM_UPPER_BOUND = 1000.0
METHANE_LOWER_BOUND = -1000.0
METHANE_UPPER_BOUND = 1000.0

# Biomass coefficients
FORMATE_BIOMASS_COEFFICIENT = 0.69


def download_model(url, output_path):
    """
    Download SBML model from GitLab repository.

    Parameters
    ----------
    url : str
        URL to download the model from
    output_path : Path or str
        Path where to save the downloaded model

    Raises
    ------
    urllib.error.URLError
        If download fails
    """
    output_path = Path(output_path)
    output_dir = output_path.parent
    
    # Create data directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Downloading model from GitLab...")
    print(f"URL: {url}")
    
    try:
        # Add User-Agent header to avoid 403 Forbidden errors from GitLab
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response, open(output_path, 'wb') as out_file:
            out_file.write(response.read())
        print(f"Model downloaded successfully to: {output_path}")
    except urllib.error.URLError as e:
        raise urllib.error.URLError(
            f"Failed to download model from {url}\n"
            f"Error: {e}\n"
            f"Please check your internet connection and the URL."
        )


def load_model(model_path):
    """
    Load SBML metabolic model from file.
    
    Downloads the model from GitLab first, then loads it.

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
        If the model file does not exist after download attempt
    """
    # Always download the model first
    download_model(MODEL_URL, model_path)
    
    print(f"Loading model from: {model_path}")
    model = cobra.io.read_sbml_model(str(model_path))
    print(f"Model loaded successfully: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")
    return model


def add_formate_metabolite(model):
    """
    Add extracellular formate metabolite to the model.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model to modify

    Returns
    -------
    cobra.Metabolite
        Created formate metabolite
    """
    formate_ext = cobra.Metabolite(
        id=FORMATE_EXTRACELLULAR,
        name="formate",
        compartment="e"
    )
    
    model.add_metabolites([formate_ext])
    print(f"Added metabolite: {formate_ext.id} - {formate_ext.name}")
    
    return formate_ext


def add_formate_transport(model):
    """
    Add formate transport reaction between cytoplasm and extracellular space.

    The reaction transports formate from cytoplasm to extracellular:
    for_c <=> for_e

    Parameters
    ----------
    model : cobra.Model
        Metabolic model to modify

    Returns
    -------
    cobra.Reaction
        Created transport reaction
    """
    transport = cobra.Reaction(FORMATE_TRANSPORT_ID)
    transport.name = "formate transport"
    transport.subsystem = "Exchange"
    transport.lower_bound = TRANSPORT_LOWER_BOUND
    transport.upper_bound = TRANSPORT_UPPER_BOUND
    
    # Define stoichiometry: for_c -> for_e
    transport.add_metabolites({
        model.metabolites.get_by_id(FORMATE_CYTOPLASMIC): -1.0,
        model.metabolites.get_by_id(FORMATE_EXTRACELLULAR): 1.0
    })
    
    model.add_reactions([transport])
    print(f"Added reaction: {transport.id} - {transport.name}")
    print(f"  Reaction: {transport.reaction}")
    
    return transport


def add_formate_exchange(model):
    """
    Add exchange boundary for extracellular formate.

    This allows the model to import/export formate from/to the environment.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model to modify

    Returns
    -------
    cobra.Reaction
        Created exchange reaction
    """
    exchange = model.add_boundary(
        model.metabolites.get_by_id(FORMATE_EXTRACELLULAR),
        type="exchange"
    )
    
    print(f"Added exchange: {exchange.id} - {exchange.name}")
    
    return exchange


def set_default_bounds(model):
    """
    Set default flux bounds for key reactions.

    Sets bounds for:
    - Formate exchange (closed by default)
    - NGAM (maintenance)
    - Methane exchange

    Parameters
    ----------
    model : cobra.Model
        Metabolic model to modify
    """
    # Formate exchange (closed by default)
    model.reactions.get_by_id(FORMATE_EXCHANGE_ID).lower_bound = EXCHANGE_LOWER_BOUND
    model.reactions.get_by_id(FORMATE_EXCHANGE_ID).upper_bound = EXCHANGE_UPPER_BOUND
    print(f"Set {FORMATE_EXCHANGE_ID} bounds: "
          f"[{EXCHANGE_LOWER_BOUND}, {EXCHANGE_UPPER_BOUND}]")
    
    # NGAM (maintenance energy)
    model.reactions.get_by_id(NGAM_ID).upper_bound = NGAM_UPPER_BOUND
    print(f"Set {NGAM_ID} upper bound: {NGAM_UPPER_BOUND}")
    
    # Methane exchange
    model.reactions.get_by_id(METHANE_EXCHANGE_ID).lower_bound = METHANE_LOWER_BOUND
    model.reactions.get_by_id(METHANE_EXCHANGE_ID).upper_bound = METHANE_UPPER_BOUND
    print(f"Set {METHANE_EXCHANGE_ID} bounds: "
          f"[{METHANE_LOWER_BOUND}, {METHANE_UPPER_BOUND}]")


def modify_biomass_formate_content(model):
    """
    Modify biomass equation to include formate with specified coefficient.

    Adjusts the formate requirement in the biomass objective function (BOF).

    Parameters
    ----------
    model : cobra.Model
        Metabolic model to modify
    """
    bof = model.reactions.get_by_id(BIOMASS_ID)
    formate = model.metabolites.get_by_id(FORMATE_CYTOPLASMIC)
    
    # Add or update formate coefficient in biomass
    bof.add_metabolites({formate: FORMATE_BIOMASS_COEFFICIENT})
    
    # Get the actual coefficient (may differ if metabolite already existed)
    actual_coeff = bof.metabolites[formate]
    print(f"Biomass formate content set to: {actual_coeff}")


def save_model(model, output_path):
    """
    Save modified model to SBML file.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model to save
    output_path : Path or str
        Output file path
    """
    # Ensure output directory exists
    output_dir = Path(output_path).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Saving modified model to: {output_path}")
    cobra.io.write_sbml_model(model, str(output_path))
    print("Model saved successfully")


def main():
    """
    Main execution function.

    Workflow:
    1. Load base model
    2. Add formate metabolite
    3. Add formate transport reaction
    4. Add formate exchange boundary
    5. Set default bounds
    6. Modify biomass formate content
    7. Save modified model
    """
    print("=" * 70)
    print("Model Formate Modification")
    print("=" * 70)
    print()
    
    try:
        # Load model
        model = load_model(MODEL_INPUT_PATH)
        print()
        
        # Add formate metabolism
        print("Adding formate metabolism...")
        add_formate_metabolite(model)
        add_formate_transport(model)
        add_formate_exchange(model)
        print()
        
        # Set constraints
        print("Setting default bounds...")
        set_default_bounds(model)
        print()
        
        # Modify biomass
        print("Modifying biomass composition...")
        modify_biomass_formate_content(model)
        print()
        
        # Save model
        save_model(model, MODEL_OUTPUT_PATH)
        print()
        
        print("=" * 70)
        print("Modification completed successfully!")
        print("=" * 70)
        
    except Exception as e:
        print(f"\nError: {e}")
        raise


if __name__ == "__main__":
    main()
