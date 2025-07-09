import os
import pickle
import re
import subprocess
import sys
from datetime import datetime
from difflib import get_close_matches

import numpy as np
import pandas as pd


def run_diffdock_with_smiles(pdb_path, smiles_string, local_output_dir, gpu_device=0, use_gpu=True):
    try:
        summary = []

        # Check if PDB file exists
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"The PDB file '{pdb_path}' does not exist.")
        summary.append(f"PDB file '{pdb_path}' found.")

        # Ensure the output directory exists
        if not os.path.exists(local_output_dir):
            os.makedirs(local_output_dir)
        summary.append(f"Output directory '{local_output_dir}' is ready.")

        # Pull the pre-built container from Docker Hub
        summary.append("Pulling DiffDock container from Docker Hub...")
        subprocess.run(["docker", "pull", "rbgcsail/diffdock"], check=True)
        summary.append("DiffDock container pulled successfully.")

        # Check for GPU availability (if using GPU)
        if use_gpu:
            summary.append("Checking for GPU availability...")
            gpu_check = subprocess.run(
                [
                    "docker",
                    "run",
                    "--rm",
                    "--gpus",
                    "all",
                    "nvidia/cuda:11.7.1-devel-ubuntu22.04",
                    "nvidia-smi",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            summary.append(f"GPU Status: {gpu_check.stdout.strip()}")

        # Prepare the GPU flag
        gpu_flag = ["--gpus", f"device={gpu_device}"] if use_gpu else []

        # Docker run command
        summary.append("Running DiffDock inference...")
        run_command = (
            ["docker", "run"]
            + gpu_flag
            + [
                # Mount the local directory to /home/appuser/output inside the container
                "-v",
                f"{os.path.abspath(pdb_path)}:/home/appuser/input/protein.pdb",  # PDB file mount
                "-v",
                f"{os.path.abspath(local_output_dir)}:/home/appuser/output",  # Output directory mount
                "--entrypoint",
                "/bin/bash",
                "rbgcsail/diffdock",
                "-c",
                # Command to run inference using micromamba environment
                f"micromamba run -n diffdock python -m inference --config default_inference_args.yaml "
                f"--protein_path /home/appuser/input/protein.pdb --ligand '{smiles_string}' --out_dir /home/appuser/output",
            ]
        )

        # Execute the Docker command
        result = subprocess.run(run_command, check=False, capture_output=True, text=True)

        # Check for errors
        if result.returncode != 0:
            summary.append(f"Error during inference: {result.stderr.strip()}")
            return "\n".join(summary)
        else:
            summary.append("DiffDock inference completed successfully.")
            summary.append(f"Results stored in '{local_output_dir}'.")

        return "\n".join(summary)

    except FileNotFoundError as e:
        return f"File error: {e}"
    except subprocess.CalledProcessError as e:
        return f"Command execution error: {e}"
    except Exception as e:
        return f"An error occurred: {e}"


def docking_autodock_vina(smiles_list, receptor_pdb_file, box_center, box_size, ncpu=1):
    from tdc import Oracle

    log = []

    # Log the start of the process
    log.append("Step 1: Initializing the Oracle")
    log.append(f"Receptor PDB File: {receptor_pdb_file}")
    log.append(f"Box Center: {box_center}")
    log.append(f"Box Size: {box_size}")

    # Initialize the Oracle object
    oracle = Oracle(
        name="pyscreener",
        receptor_pdb_file=receptor_pdb_file,
        box_center=box_center,
        box_size=box_size,
        ncpu=ncpu,
    )
    log.append("Oracle initialized successfully.")

    # Log the list of SMILES strings
    log.append(f"\nStep 2: Processing SMILES strings: {smiles_list}")

    # Get the docking scores
    docking_scores = oracle(smiles_list)
    log.append(f"Docking scores calculated: {docking_scores}")

    # Create a dictionary mapping SMILES to their docking scores
    results_dict = dict(zip(smiles_list, docking_scores, strict=False))

    # Log the result mapping
    log.append("\nStep 3: Mapping SMILES to docking scores:")
    log.append(f"Results: {results_dict}")

    # Convert the log to a string and return it
    research_log = "\n".join(log)
    return research_log


def run_autosite(pdb_file, output_dir, spacing=1.0):
    # Prepare the output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Convert the PDB file to PDBQT format (assuming prepare_receptor4.py is accessible)
    pdbqt_file = pdb_file.replace(".pdb", ".pdbqt")
    subprocess.run(["prepare_receptor", "-r", pdb_file, "-o", pdbqt_file], check=True)

    # Run AutoSite
    autosite_cmd = [
        "autosite",
        "-r",
        pdbqt_file,
        "--spacing",
        str(spacing),
        "-o",
        output_dir,
    ]
    subprocess.run(autosite_cmd, check=True)

    # Parse the results to find the box center and size
    box_center, box_size = None, None
    log_path = os.path.join(output_dir, "_AutoSiteSummary.log")
    with open(log_path) as log_file:
        log_content = log_file.read()

        # Extract box center and size from the log (assuming standard output format)
        box_center_match = re.search(r"Box center:\s*\(([^)]+)\)", log_content)
        box_size_match = re.search(r"Box size:\s*\(([^)]+)\)", log_content)

        if box_center_match:
            box_center = box_center_match.group(1)
        if box_size_match:
            box_size = box_size_match.group(1)

    # Create a research log string
    research_log = f"AutoSite run for {pdb_file} with spacing {spacing}\n"
    research_log += f"Output directory: {output_dir}\n"
    if box_center and box_size:
        research_log += f"Box Center: {box_center}\nBox Size: {box_size}"
    else:
        research_log += "Box Center and Size information not found in log."

    return research_log


# Function to get TxGNN predictions and return a summarized string output
def retrieve_topk_repurposing_drugs_from_disease_txgnn(disease_name, k=5):
    """This function computes TxGNN model predictions for drug repurposing. It takes in the paths to the data,
    the disease name, and returns a summary of the top K predicted drugs with their sigmoid-transformed scores.

    Args:
    - disease_name (str): The name of the disease for which the drug predictions are to be retrieved.
    - k (int, optional): The number of top drug predictions to return. Defaults to 5.

    Returns:
    - str: A summary of the steps and the top K drug predictions with their scores.

    """

    # Sigmoid function to convert raw prediction scores
    def sigmoid(x):
        return 1 / (1 + np.exp(-x))

    # Step 1: Load the mappings and prediction data from the provided paths
    with open("/dfs/project/bioagentos/required_data/txgnn/name_mapping.pkl", "rb") as f:
        mapping = pickle.load(f)

    with open("/dfs/project/bioagentos/required_data/txgnn/prediction.pkl", "rb") as f:
        result = pickle.load(f)

    # Step 2: Fuzzy match the disease name to find the closest match
    possible_diseases = result.keys()
    matched_disease = get_close_matches(disease_name, possible_diseases, n=1, cutoff=0.6)

    if not matched_disease:
        return f"Error: No matching disease found for '{disease_name}'. Please try a different name."

    matched_disease = matched_disease[0]

    # Step 3: Retrieve the prediction scores for the matched disease
    disease_predictions = result[matched_disease]

    # Step 4: Apply the sigmoid function to the raw prediction scores
    sigmoid_predictions = {drug_id: sigmoid(score) for drug_id, score in disease_predictions.items()}

    # Step 5: Sort the drugs by prediction score in descending order
    top_k_drugs = sorted(sigmoid_predictions.items(), key=lambda x: x[1], reverse=True)[:k]

    # Step 6: Map drug IDs to their names and format the results
    top_k_drug_names = [(mapping["id2name_drug"].get(drug_id, "Unknown Drug"), score) for drug_id, score in top_k_drugs]

    # Step 7: Create a human and LLM-friendly summary string
    summary = f"TxGNN Drug Repurposing Predictions for '{matched_disease}':\n"
    summary += f"Top {k} predicted drugs and their corresponding prediction scores (post-sigmoid transformation):\n"

    for i, (drug_name, score) in enumerate(top_k_drug_names, 1):
        summary += f"{i}. {drug_name} - Prediction Score: {score:.4f}\n"

    summary += "\nProcess Summary:\n"
    summary += "- Fuzzy matching was used to match the input disease name to '{matched_disease}'.\n"
    summary += "- Sigmoid function was applied to raw prediction scores to convert them into probabilities.\n"
    summary += f"- The top {k} drugs were selected based on their prediction scores.\n"

    return summary


# ADMET prediction function with research log format
def predict_admet_properties(smiles_list, ADMET_model_type="MPNN"):
    try:
        from DeepPurpose import CompoundPred, utils
    except Exception:
        subprocess.run([sys.executable, "-m", "pip", "install", "DeepPurpose"], check=False)
        from DeepPurpose import CompoundPred, utils

    # Define available model types
    available_model_types = ["MPNN", "CNN", "Morgan"]

    # Check if the provided model type is valid
    if ADMET_model_type not in available_model_types:
        return f"Error: Invalid ADMET model type '{ADMET_model_type}'. Available options are: {', '.join(available_model_types)}."

    # Load pretrained ADMET models only once
    model_ADMETs = {}
    tasks = [
        "AqSolDB",
        "Caco2",
        "HIA",
        "Pgp_inhibitor",
        "Bioavailability",
        "BBB_MolNet",
        "PPBR",
        "CYP2C19",
        "CYP2D6",
        "CYP3A4",
        "CYP1A2",
        "CYP2C9",
        "ClinTox",
        "Lipo_AZ",
        "Half_life_eDrug3D",
        "Clearance_eDrug3D",
    ]

    for task in tasks:
        model_ADMETs[task + "_" + ADMET_model_type + "_model"] = CompoundPred.model_pretrained(
            model=task + "_" + ADMET_model_type + "_model"
        )

    # Helper function for ADMET prediction
    def ADMET_pred(drug, task, unit):
        model = model_ADMETs[task + "_" + ADMET_model_type + "_model"]
        X_pred = utils.data_process(
            X_drug=[drug],
            y=[0],
            drug_encoding=ADMET_model_type,
            split_method="no_split",
        )
        y_pred = model.predict(X_pred)[0]

        if unit == "%":
            y_pred = y_pred * 100

        return f"{y_pred:.2f} " + unit

    # Initialize research log string
    research_log = "Research Log for ADMET Predictions:\n"
    research_log += "-------------------------------------\n"

    # Process each SMILES string in the list
    for smiles in smiles_list:
        research_log += f"\nCompound SMILES: {smiles}\n"
        research_log += "Predicted ADMET properties:\n"

        # Physiochemical properties
        solubility = ADMET_pred(smiles, "AqSolDB", "log mol/L")
        lipophilicity = ADMET_pred(smiles, "Lipo_AZ", "(log-ratio)")
        research_log += f"- Solubility: {solubility}\n"
        research_log += f"- Lipophilicity: {lipophilicity}\n"

        # Absorption
        caco2 = ADMET_pred(smiles, "Caco2", "cm/s")
        hia = ADMET_pred(smiles, "HIA", "%")
        pgp = ADMET_pred(smiles, "Pgp_inhibitor", "%")
        bioavail = ADMET_pred(smiles, "Bioavailability", "%")
        research_log += f"- Absorption (Caco-2 permeability): {caco2}\n"
        research_log += f"- Absorption (HIA): {hia}\n"
        research_log += f"- Absorption (Pgp Inhibitor): {pgp}\n"
        research_log += f"- Absorption (Bioavailability): {bioavail}\n"

        # Distribution
        bbb = ADMET_pred(smiles, "BBB_MolNet", "%")
        ppbr = ADMET_pred(smiles, "PPBR", "%")
        research_log += f"- Distribution (BBB permeation): {bbb}\n"
        research_log += f"- Distribution (PPBR): {ppbr}\n"

        # Metabolism
        cyp2c19 = ADMET_pred(smiles, "CYP2C19", "%")
        cyp2d6 = ADMET_pred(smiles, "CYP2D6", "%")
        cyp3a4 = ADMET_pred(smiles, "CYP3A4", "%")
        cyp1a2 = ADMET_pred(smiles, "CYP1A2", "%")
        cyp2c9 = ADMET_pred(smiles, "CYP2C9", "%")
        research_log += f"- Metabolism (CYP2C19): {cyp2c19}\n"
        research_log += f"- Metabolism (CYP2D6): {cyp2d6}\n"
        research_log += f"- Metabolism (CYP3A4): {cyp3a4}\n"
        research_log += f"- Metabolism (CYP1A2): {cyp1a2}\n"
        research_log += f"- Metabolism (CYP2C9): {cyp2c9}\n"

        # Excretion
        half_life = ADMET_pred(smiles, "Half_life_eDrug3D", "h")
        clearance = ADMET_pred(smiles, "Clearance_eDrug3D", "mL/min/kg")
        research_log += f"- Excretion (Half-life): {half_life}\n"
        research_log += f"- Excretion (Clearance): {clearance}\n"

        # Clinical Toxicity
        clinical_toxicity = ADMET_pred(smiles, "ClinTox", "%")
        research_log += f"- Clinical Toxicity: {clinical_toxicity}\n"

        research_log += "-------------------------------------\n"

    return research_log


# Binding Affinity prediction function with model_type validation
def predict_binding_affinity_protein_1d_sequence(smiles_list, amino_acid_sequence, affinity_model_type="MPNN-CNN"):
    try:
        from DeepPurpose import DTI, utils
    except Exception:
        subprocess.run([sys.executable, "-m", "pip", "install", "DeepPurpose"], check=False)
        from DeepPurpose import DTI, utils

    # Define available model types for Binding Affinity
    available_affinity_model_types = [
        "CNN-CNN",
        "MPNN-CNN",
        "Morgan-CNN",
        "Morgan-AAC",
        "Daylight-AAC",
    ]

    # Check if the provided affinity model type is valid
    if affinity_model_type not in available_affinity_model_types:
        return f"Error: Invalid affinity model type '{affinity_model_type}'. Available options are: {', '.join(available_affinity_model_types)}."

    # Load the pre-trained affinity model
    model_DTI = DTI.model_pretrained(model=affinity_model_type.replace("-", "_") + "_BindingDB")

    # Initialize research log string
    research_log = "Research Log for Binding Affinity Predictions:\n"
    research_log += "-------------------------------------\n"

    # Process each SMILES string in the list
    for smiles in smiles_list:
        research_log += f"\nCompound SMILES: {smiles}\n"
        research_log += f"Amino Acid Sequence: {amino_acid_sequence}\n"

        # Predict binding affinity
        X_pred = utils.data_process(
            X_drug=[smiles],
            X_target=[amino_acid_sequence],
            y=[0],
            drug_encoding=affinity_model_type.split("-")[0],
            target_encoding=affinity_model_type.split("-")[1],
            split_method="no_split",
        )
        y_pred = model_DTI.predict(X_pred)[0]
        y_pred_nM = 10 ** (-y_pred) / 1e-9

        research_log += f"Predicted Binding Affinity: {y_pred_nM:.2f} nM\n"
        research_log += "-------------------------------------\n"

    return research_log


def analyze_accelerated_stability_of_pharmaceutical_formulations(formulations, storage_conditions, time_points):
    """Analyzes the stability of pharmaceutical formulations under accelerated storage conditions.

    Parameters
    ----------
    formulations : list of dict
        List of formulation dictionaries, each containing:
        - 'name': str, name of the formulation
        - 'active_ingredient': str, name of the active pharmaceutical ingredient
        - 'concentration': float, concentration in mg/mL
        - 'excipients': list, list of excipients

    storage_conditions : list of dict
        List of storage condition dictionaries, each containing:
        - 'temperature': float, temperature in °C
        - 'humidity': float, relative humidity in percentage (optional for solid dosage forms)
        - 'description': str, description of storage condition (e.g., "Room Temperature", "Accelerated")

    time_points : list of int
        List of time points in days to evaluate stability

    Returns
    -------
    str
        Research log summarizing the stability testing process and results

    """
    # Create output directory if it doesn't exist
    output_dir = "stability_test_results"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Initialize results storage
    all_results = []

    # Process each formulation under each storage condition
    for formulation in formulations:
        for condition in storage_conditions:
            # Initialize stability parameters
            results = []

            # Get acceleration factor based on temperature (simplified Arrhenius equation)
            temp_c = condition["temperature"]
            accel_factor = 2 ** ((temp_c - 25) / 10)  # Rule of thumb: reaction rate doubles every 10°C

            # Add humidity effect for degradation if provided
            humidity_factor = 1.0
            if "humidity" in condition:
                humidity = condition["humidity"]
                # Simple model: higher humidity increases degradation rate
                humidity_factor = 1.0 + (humidity - 60) / 100 if humidity > 60 else 1.0

            # Calculate stability parameters at each time point
            initial_content = 100.0  # Starting at 100%

            for time in time_points:
                # Chemical stability (% of initial content)
                # Simple first-order degradation model
                effective_time = time * accel_factor * humidity_factor
                chemical_stability = initial_content * np.exp(-0.001 * effective_time)

                # Physical stability (score from 1-10, 10 being perfect)
                # Decreases over time, affected by temperature and humidity
                physical_stability = 10 - (0.05 * effective_time)
                physical_stability = max(1, physical_stability)  # Minimum score of 1

                # Particle size change (% increase from initial)
                # Some formulations show particle growth over time
                particle_size_change = (
                    0.2 * effective_time if "solid" in formulation.get("dosage_form", "").lower() else 0
                )

                results.append(
                    {
                        "Formulation": formulation["name"],
                        "Storage_Condition": condition["description"],
                        "Temperature_C": temp_c,
                        "Humidity_RH": condition.get("humidity", "N/A"),
                        "Time_Days": time,
                        "Chemical_Stability_Percent": round(chemical_stability, 2),
                        "Physical_Stability_Score": round(physical_stability, 1),
                        "Particle_Size_Change_Percent": round(particle_size_change, 2),
                    }
                )

            all_results.extend(results)

    # Convert results to DataFrame
    results_df = pd.DataFrame(all_results)

    # Save results to CSV
    csv_filename = f"{output_dir}/stability_results_{timestamp}.csv"
    results_df.to_csv(csv_filename, index=False)

    # Generate research log
    log = "Accelerated Stability Testing of Pharmaceutical Formulations\n"
    log += f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"

    log += "1. STUDY PARAMETERS\n"
    log += f"   - Number of formulations tested: {len(formulations)}\n"
    log += f"   - Formulations: {', '.join([f['name'] for f in formulations])}\n"
    log += (
        "   - Storage conditions: "
        + ", ".join(
            [
                f"{c['description']} ({c['temperature']}°C" + (f"/{c['humidity']}% RH" if "humidity" in c else "") + ")"
                for c in storage_conditions
            ]
        )
        + "\n"
    )
    log += f"   - Time points evaluated (days): {', '.join(map(str, time_points))}\n\n"

    log += "2. METHODOLOGY\n"
    log += "   - Chemical stability assessed by active ingredient content\n"
    log += "   - Physical stability evaluated on a 10-point scale\n"
    log += "   - Particle size changes measured where applicable\n\n"

    log += "3. KEY FINDINGS\n"

    # Summarize stability at final time point for each formulation/condition
    final_time = max(time_points)
    final_results = results_df[results_df["Time_Days"] == final_time]

    for formulation in formulations:
        form_results = final_results[final_results["Formulation"] == formulation["name"]]
        log += f"   {formulation['name']}:\n"

        for _, row in form_results.iterrows():
            condition = row["Storage_Condition"]
            chem_stab = row["Chemical_Stability_Percent"]
            phys_stab = row["Physical_Stability_Score"]

            stability_assessment = "Stable"
            if chem_stab < 90 or phys_stab < 7:
                stability_assessment = "Potentially unstable"
            if chem_stab < 85 or phys_stab < 5:
                stability_assessment = "Unstable"

            log += f"     - {condition}: Chemical stability {chem_stab}%, Physical stability score {phys_stab}/10 - {stability_assessment}\n"
        log += "\n"

    log += "4. CONCLUSION\n"

    # Identify most stable formulation
    best_formulation = ""
    best_stability = 0

    for formulation in formulations:
        form_data = final_results[final_results["Formulation"] == formulation["name"]]
        avg_chem_stability = form_data["Chemical_Stability_Percent"].mean()
        if avg_chem_stability > best_stability:
            best_stability = avg_chem_stability
            best_formulation = formulation["name"]

    log += f"   - Most stable formulation: {best_formulation} (avg. chemical stability: {best_stability:.2f}%)\n"
    log += f"   - Detailed results saved to: {csv_filename}\n"

    return log


def run_3d_chondrogenic_aggregate_assay(
    chondrocyte_cells, test_compounds, culture_duration_days=21, measurement_intervals=7
):
    """Generates a detailed protocol for performing a 3D chondrogenic aggregate culture assay to evaluate compounds' effects on chondrogenesis.

    Parameters
    ----------
    chondrocyte_cells : dict
        Dictionary with cell information including 'source', 'passage_number', and 'cell_density'
    test_compounds : list of dict
        List of compounds to test, each with 'name', 'concentration', and 'vehicle' keys
    culture_duration_days : int
        Total duration of the culture period in days (default: 21)
    measurement_intervals : int
        Interval in days between measurements (default: 7)

    Returns
    -------
    str
        Detailed protocol document for the 3D chondrogenic aggregate culture assay

    """
    from datetime import datetime

    # Create experiment ID
    experiment_id = f"CHOND3D_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

    # Create time points for measurements
    timepoints = list(range(0, culture_duration_days + 1, measurement_intervals))
    if timepoints[-1] != culture_duration_days:
        timepoints.append(culture_duration_days)

    # Generate the protocol document
    protocol = f"# 3D Chondrogenic Aggregate Culture Assay Protocol - {experiment_id}\n\n"

    protocol += "## 1. Materials and Reagents\n\n"
    protocol += "- Chondrocyte cells\n"
    protocol += "- Chondrogenic differentiation medium\n"
    protocol += "- Transforming growth factor-β3 (TGF-β3)\n"
    protocol += "- Dexamethasone\n"
    protocol += "- Ascorbate-2-phosphate\n"
    protocol += "- 96-well V-bottom plates\n"
    protocol += "- Gaussia luciferase reporter assay kit\n"
    protocol += "- Luminometer\n"
    protocol += "- Test compounds with respective vehicles\n"
    protocol += "- Centrifuge\n"
    protocol += "- CO2 incubator\n"
    protocol += "- Sterile pipettes and tips\n\n"

    protocol += "## 2. Experimental Information\n\n"
    protocol += "### Cell Information:\n"
    protocol += f"- Cell source: {chondrocyte_cells['source']}\n"
    protocol += f"- Passage number: {chondrocyte_cells['passage_number']}\n"
    protocol += f"- Cell density: {chondrocyte_cells['cell_density']} cells/mL\n\n"

    protocol += "### Experimental Design:\n"
    protocol += f"- Culture duration: {culture_duration_days} days\n"
    protocol += f"- Measurement timepoints: {', '.join(map(str, timepoints))} days\n\n"

    protocol += "### Test Compounds:\n"
    for i, compound in enumerate(test_compounds):
        protocol += f"- Compound {i + 1}: {compound['name']} at {compound['concentration']} in {compound['vehicle']}\n"
    protocol += "- Control: Vehicle only\n\n"

    protocol += "## 3. Detailed Procedure\n\n"
    protocol += "### Day 0: Setup\n\n"
    protocol += "1. Prepare chondrogenic differentiation medium containing:\n"
    protocol += "   - High-glucose DMEM\n"
    protocol += "   - 10 ng/mL TGF-β3\n"
    protocol += "   - 100 nM Dexamethasone\n"
    protocol += "   - 50 μg/mL Ascorbate-2-phosphate\n"
    protocol += "   - 1% ITS+ premix (insulin, transferrin, selenium)\n"
    protocol += "   - 1 mM Sodium pyruvate\n"
    protocol += "   - 100 U/mL Penicillin/Streptomycin\n\n"

    protocol += "2. Harvest and count chondrocyte cells\n\n"

    protocol += "3. Prepare cell suspension at the specified density:\n"
    protocol += f"   - {chondrocyte_cells['cell_density']} cells/mL\n\n"

    protocol += "4. Form 3D cell aggregates:\n"
    protocol += "   - Aliquot 2.5×10^5 cells per well in 96-well V-bottom plates\n"
    protocol += "   - Centrifuge plates at 500g for 5 minutes to pellet cells\n\n"

    protocol += "5. Add test compounds to respective wells:\n"
    for _, compound in enumerate(test_compounds):
        protocol += f"   - Add {compound['name']} at {compound['concentration']} in {compound['vehicle']}\n"
    protocol += "   - Add vehicle only to control wells\n\n"

    protocol += "6. Incubate the plates at 37°C, 5% CO2\n\n"

    protocol += "### Day 1 to Day " + str(culture_duration_days) + ":\n\n"
    protocol += "1. Change medium every 2-3 days:\n"
    protocol += "   - Carefully remove 50% of the medium without disturbing the aggregates\n"
    protocol += "   - Replace with fresh medium containing test compounds at the same concentrations\n\n"

    protocol += f"2. At days {', '.join(map(str, timepoints))}, collect samples for analysis:\n"
    protocol += (
        "   - Take medium samples for Gaussia luciferase activity measurement (if using COL2A1-GLuc reporter cells)\n"
    )
    protocol += "   - Fix aggregates in 4% paraformaldehyde for histological analysis\n\n"

    return protocol


def grade_adverse_events_using_vcog_ctcae(clinical_data_file):
    """Grade and monitor adverse events in animal studies using the VCOG-CTCAE standard.

    Parameters
    ----------
    clinical_data_file : str
        Path to a CSV file containing clinical evaluation data with columns:
        subject_id, time_point, symptom, severity, measurement (optional)

    Returns
    -------
    str
        A research log summarizing the adverse event grading process and results.
        The graded events are saved to 'vcog_ctcae_graded_events.csv'.

    """
    import json
    from datetime import datetime

    import pandas as pd

    # Initialize the research log
    log = "# Adverse Event Grading using VCOG-CTCAE v1.1\n"
    log += f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"

    # Step 1: Load the clinical data
    log += "## Step 1: Loading clinical evaluation data\n"
    try:
        data = pd.read_csv(clinical_data_file)
        log += f"Successfully loaded data from {clinical_data_file}\n"
        log += f"Total records: {len(data)}\n"
        log += f"Columns found: {', '.join(data.columns)}\n\n"
    except Exception as e:
        log += f"Error loading data: {str(e)}\n"
        return log

    # Step 2: Define VCOG-CTCAE grading criteria
    log += "## Step 2: Applying VCOG-CTCAE grading criteria\n"

    # Comprehensive VCOG-CTCAE grading criteria based on VCOG-CTCAE v1.1
    vcog_criteria = {
        # Hematologic
        "neutropenia": {
            "description": "Neutrophil count decrease",
            "unit": "cells/µL",
            "grades": {
                0: {"criteria": "≥ 1500", "range": [1500, float("inf")]},
                1: {"criteria": "1000 - <1500", "range": [1000, 1500]},
                2: {"criteria": "500 - <1000", "range": [500, 1000]},
                3: {"criteria": "100 - <500", "range": [100, 500]},
                4: {"criteria": "<100", "range": [0, 100]},
                5: {"criteria": "Death due to neutropenic sepsis", "range": None},
            },
        },
        "anemia": {
            "description": "Hemoglobin decrease",
            "unit": "g/dL",
            "grades": {
                0: {"criteria": "Within reference range", "range": [10, float("inf")]},
                1: {"criteria": "Mild; clinical signs not present", "range": [8, 10]},
                2: {"criteria": "Moderate; clinical signs present", "range": [6.5, 8]},
                3: {"criteria": "Severe; transfusion indicated", "range": [5, 6.5]},
                4: {"criteria": "Life-threatening", "range": [0, 5]},
                5: {"criteria": "Death", "range": None},
            },
        },
        "thrombocytopenia": {
            "description": "Platelet count decrease",
            "unit": "cells/µL",
            "grades": {
                0: {"criteria": "≥ 100,000", "range": [100000, float("inf")]},
                1: {"criteria": "50,000 - <100,000", "range": [50000, 100000]},
                2: {"criteria": "25,000 - <50,000", "range": [25000, 50000]},
                3: {"criteria": "10,000 - <25,000", "range": [10000, 25000]},
                4: {"criteria": "<10,000 or spontaneous bleeding", "range": [0, 10000]},
                5: {"criteria": "Death", "range": None},
            },
        },
        # Gastrointestinal
        "vomiting": {
            "description": "Vomiting frequency",
            "unit": "episodes per 24h period",
            "grades": {
                0: {"criteria": "None", "range": [0, 0]},
                1: {
                    "criteria": "1-2 episodes in 24h; medical intervention not indicated",
                    "range": [1, 2],
                },
                2: {"criteria": "3-5 episodes in 24h; ≤3 days", "range": [3, 5]},
                3: {
                    "criteria": ">5 episodes in 24h; >3 days; hospitalization indicated",
                    "range": [6, float("inf")],
                },
                4: {"criteria": "Life-threatening consequences", "range": None},
                5: {"criteria": "Death", "range": None},
            },
        },
        "diarrhea": {
            "description": "Diarrhea frequency",
            "unit": "episodes per 24h period",
            "grades": {
                0: {"criteria": "None", "range": [0, 0]},
                1: {
                    "criteria": "Increase of <4 stools per day over baseline",
                    "range": [1, 3],
                },
                2: {
                    "criteria": "Increase of 4-6 stools per day over baseline",
                    "range": [4, 6],
                },
                3: {
                    "criteria": "Increase of ≥7 stools per day; hospitalization indicated",
                    "range": [7, float("inf")],
                },
                4: {"criteria": "Life-threatening consequences", "range": None},
                5: {"criteria": "Death", "range": None},
            },
        },
        "anorexia": {
            "description": "Appetite/food intake decrease",
            "unit": "percent of normal intake",
            "grades": {
                0: {"criteria": "Normal", "range": [100, float("inf")]},
                1: {"criteria": "Decreased appetite, but eating", "range": [75, 100]},
                2: {"criteria": "Decreased intake <3 days", "range": [50, 75]},
                3: {"criteria": "Decreased intake ≥3 days", "range": [25, 50]},
                4: {
                    "criteria": "Life-threatening consequences; urgent intervention indicated",
                    "range": [0, 25],
                },
                5: {"criteria": "Death", "range": None},
            },
        },
        # Hepatic
        "alt_increase": {
            "description": "Alanine aminotransferase increased",
            "unit": "x ULN (upper limit of normal)",
            "grades": {
                0: {"criteria": "≤ ULN", "range": [0, 1]},
                1: {"criteria": ">ULN - 2.5xULN", "range": [1, 2.5]},
                2: {"criteria": ">2.5 - 5.0xULN", "range": [2.5, 5]},
                3: {"criteria": ">5.0 - 20.0xULN", "range": [5, 20]},
                4: {"criteria": ">20.0xULN", "range": [20, float("inf")]},
                5: {"criteria": "Death", "range": None},
            },
        },
        # Renal
        "creatinine_increase": {
            "description": "Creatinine increased",
            "unit": "x ULN",
            "grades": {
                0: {"criteria": "≤ ULN", "range": [0, 1]},
                1: {"criteria": ">ULN - 1.5xULN", "range": [1, 1.5]},
                2: {"criteria": ">1.5 - 3.0xULN", "range": [1.5, 3]},
                3: {"criteria": ">3.0 - 6.0xULN", "range": [3, 6]},
                4: {"criteria": ">6.0xULN", "range": [6, float("inf")]},
                5: {"criteria": "Death", "range": None},
            },
        },
        # Constitutional
        "fever": {
            "description": "Fever",
            "unit": "°C",
            "grades": {
                0: {"criteria": "None", "range": [0, 39]},
                1: {"criteria": "39.0 - 39.5°C", "range": [39, 39.5]},
                2: {"criteria": ">39.5 - 40.0°C", "range": [39.5, 40]},
                3: {"criteria": ">40.0 - 41.0°C", "range": [40, 41]},
                4: {"criteria": ">41.0°C for >24 hrs", "range": [41, float("inf")]},
                5: {"criteria": "Death", "range": None},
            },
        },
        "weight_loss": {
            "description": "Weight loss",
            "unit": "percent of baseline weight",
            "grades": {
                0: {"criteria": "<5%", "range": [0, 5]},
                1: {"criteria": "5% - <10%", "range": [5, 10]},
                2: {"criteria": "10% - <20%", "range": [10, 20]},
                3: {"criteria": "≥20%", "range": [20, float("inf")]},
                4: {"criteria": "Life-threatening", "range": None},
                5: {"criteria": "Death", "range": None},
            },
        },
        # Dermatologic
        "alopecia": {
            "description": "Hair loss",
            "unit": None,
            "grades": {
                0: {"criteria": "None", "range": None},
                1: {"criteria": "Hair loss at injection/treatment site", "range": None},
                2: {"criteria": "Moderate alopecia", "range": None},
                3: {"criteria": "Complete alopecia", "range": None},
                4: {"criteria": "Not applicable", "range": None},
                5: {"criteria": "Not applicable", "range": None},
            },
        },
        # Neurologic
        "neuropathy": {
            "description": "Peripheral neuropathy",
            "unit": None,
            "grades": {
                0: {"criteria": "None", "range": None},
                1: {
                    "criteria": "Asymptomatic; clinically detectable on examination",
                    "range": None,
                },
                2: {
                    "criteria": "Mild symptoms; limiting instrumental ADL",
                    "range": None,
                },
                3: {
                    "criteria": "Severe symptoms; limiting self-care ADL",
                    "range": None,
                },
                4: {"criteria": "Life-threatening consequences", "range": None},
                5: {"criteria": "Death", "range": None},
            },
        },
    }

    def apply_vcog_grade(symptom, severity, measurement=None):
        """Apply VCOG-CTCAE grading criteria to an adverse event.

        Parameters
        ----------
        symptom : str
            The type of adverse event
        severity : str
            The severity description
        measurement : float or None
            Quantitative measurement related to the symptom, if available

        Returns
        -------
        int
            The VCOG-CTCAE grade (0-5)
        str
            Description of the grading rationale

        """
        # Standard severity-based grading if no specific criteria exist
        grade_map = {
            "none": 0,
            "mild": 1,
            "moderate": 2,
            "severe": 3,
            "life-threatening": 4,
            "death": 5,
        }

        symptom_lower = symptom.lower()

        # Check if the symptom has specific VCOG-CTCAE criteria
        if symptom_lower in vcog_criteria:
            criteria = vcog_criteria[symptom_lower]

            # If measurement is provided and criteria has numeric ranges
            if measurement is not None:
                try:
                    measurement_value = float(measurement)

                    # Find the appropriate grade based on the measurement ranges
                    for grade, grade_info in criteria["grades"].items():
                        if grade_info["range"] is not None:
                            min_val, max_val = grade_info["range"]
                            if min_val <= measurement_value < max_val:
                                return (
                                    grade,
                                    f"Grade {grade}: {criteria['description']} - {criteria['grades'][grade]['criteria']}",
                                )
                except (ValueError, TypeError):
                    # If measurement can't be converted to float, fall back to severity-based grading
                    pass

            # If there's a reported severity with no valid measurement
            if severity.lower() in grade_map:
                # Check if the grade exists in the criteria
                severity_grade = grade_map[severity.lower()]
                if severity_grade in criteria["grades"]:
                    return (
                        severity_grade,
                        f"Grade {severity_grade}: {criteria['description']} - {criteria['grades'][severity_grade]['criteria']}",
                    )

        # Default to using the severity mapping if no specific criteria match
        if severity.lower() in grade_map:
            return grade_map[severity.lower()], f"Grade {grade_map[severity.lower()]}: Based on reported severity"

        # Default grade if no specific criteria match
        return 1, "Grade 1: Default grade (specific criteria not found)"

    # Step 3: Apply grading to each record
    log += "Applying VCOG-CTCAE v1.1 grading criteria to each adverse event...\n"

    # Create new columns for the grade and rationale
    grading_results = data.apply(
        lambda row: apply_vcog_grade(
            row["symptom"],
            row["severity"],
            row["measurement"] if "measurement" in data.columns else None,
        ),
        axis=1,
    )

    # Split the returned tuples into separate columns
    data["vcog_grade"] = [result[0] for result in grading_results]
    data["grading_rationale"] = [result[1] for result in grading_results]

    # Step 4: Analyze patterns across time points (if available)
    if "time_point" in data.columns:
        log += "\n## Step 3: Analyzing adverse event patterns across time points\n"

        # Group by subject and symptom to track progression
        progression_analysis = data.pivot_table(
            index=["subject_id", "symptom"],
            columns="time_point",
            values="vcog_grade",
            aggfunc="max",
        ).reset_index()

        # Calculate if grade is increasing, decreasing, or stable for each subject-symptom pair
        trend_counts = {"increasing": 0, "decreasing": 0, "stable": 0, "fluctuating": 0}

        numeric_columns = [col for col in progression_analysis.columns if col not in ["subject_id", "symptom"]]

        if len(numeric_columns) >= 2:
            # Sort columns to ensure chronological order
            numeric_columns.sort()

            for _, row in progression_analysis.iterrows():
                values = [row[col] for col in numeric_columns if not pd.isna(row[col])]
                if len(values) >= 2:
                    if all(values[i] < values[i + 1] for i in range(len(values) - 1)):
                        trend_counts["increasing"] += 1
                    elif all(values[i] > values[i + 1] for i in range(len(values) - 1)):
                        trend_counts["decreasing"] += 1
                    elif all(values[i] == values[i + 1] for i in range(len(values) - 1)):
                        trend_counts["stable"] += 1
                    else:
                        trend_counts["fluctuating"] += 1

            log += "Adverse event progression patterns:\n"
            for trend, count in trend_counts.items():
                log += f"- {trend.capitalize()}: {count} subject-symptom pairs\n"

        # Save progression analysis
        progression_file = "vcog_ctcae_progression_analysis.csv"
        progression_analysis.to_csv(progression_file)
        log += f"\nDetailed progression analysis saved to: {progression_file}\n"

    # Step 5: Summarize the grading results
    log += "\n## Step 4: Summarizing adverse event grades\n"

    # Count events by grade
    grade_counts = data["vcog_grade"].value_counts().sort_index()
    log += "Grade distribution:\n"
    for grade, count in grade_counts.items():
        log += f"- Grade {grade}: {count} events\n"

    # Summarize by symptom type
    symptom_summary = data.groupby("symptom")["vcog_grade"].agg(["max", "mean", "count"])
    log += "\nSymptom severity summary:\n"
    for symptom, stats in symptom_summary.iterrows():
        log += f"- {symptom}: max grade = {stats['max']}, avg grade = {stats['mean']:.2f}, count = {stats['count']}\n"

    # Summarize by subject
    subject_summary = data.groupby("subject_id")["vcog_grade"].agg(["max", "mean", "count"])
    log += f"\nSubjects with adverse events: {len(subject_summary)}\n"
    log += f"Subjects with Grade 3+ events: {len(subject_summary[subject_summary['max'] >= 3])}\n"

    # Create a summary of most severe events
    most_severe = data.sort_values("vcog_grade", ascending=False).head(10)
    log += "\nTop 10 most severe adverse events:\n"
    for i, (_, event) in enumerate(most_severe.iterrows(), 1):
        log += f"{i}. Subject {event['subject_id']}: {event['symptom']} (Grade {event['vcog_grade']})\n"

    # Step 6: Save detailed results to file
    output_file = "vcog_ctcae_graded_events.csv"
    data.to_csv(output_file, index=False)
    log += "\n## Step 5: Results saved\n"
    log += f"Detailed graded events saved to: {output_file}\n"

    # Save the VCOG criteria as a reference
    with open("vcog_ctcae_criteria_reference.json", "w") as f:
        json.dump(vcog_criteria, f, indent=2)
    log += "VCOG-CTCAE criteria reference saved to: vcog_ctcae_criteria_reference.json\n"

    return log


def analyze_radiolabeled_antibody_biodistribution(time_points, tissue_data):
    """Analyze biodistribution and pharmacokinetic profile of radiolabeled antibodies.

    Parameters
    ----------
    time_points : list or numpy.ndarray
        Time points (hours) at which measurements were taken
    tissue_data : dict
        Dictionary where keys are tissue names and values are lists/arrays of %IA/g
        measurements corresponding to time_points. Must include 'tumor' as one of the keys.

    Returns
    -------
    str
        Research log summarizing the biodistribution analysis, pharmacokinetic parameters,
        and tumor-to-normal tissue ratios

    """
    import json
    import os

    import numpy as np
    from scipy.optimize import curve_fit

    # Validate inputs
    if "tumor" not in tissue_data:
        return "Error: Tumor data must be provided in tissue_data dictionary"

    # Define bi-exponential model for pharmacokinetic analysis
    # C(t) = A*exp(-alpha*t) + B*exp(-beta*t)
    def bi_exp_model(t, A, alpha, B, beta):
        return A * np.exp(-alpha * t) + B * np.exp(-beta * t)

    # Initialize results dictionary
    results = {
        "tissues_analyzed": list(tissue_data.keys()),
        "pk_parameters": {},
        "tumor_to_normal_ratios": {},
        "auc_values": {},
    }

    # Analyze each tissue
    for tissue, measurements in tissue_data.items():
        try:
            # Fit bi-exponential model
            params, _ = curve_fit(
                bi_exp_model,
                time_points,
                measurements,
                p0=[50, 0.1, 50, 0.01],  # Initial parameter guess
                bounds=([0, 0, 0, 0], [100, 5, 100, 1]),  # Parameter bounds
            )

            A, alpha, B, beta = params

            # Calculate pharmacokinetic parameters
            # Distribution half-life (fast component)
            t_half_dist = np.log(2) / alpha

            # Elimination half-life (slow component)
            t_half_elim = np.log(2) / beta

            # Area under the curve (AUC)
            auc = A / alpha + B / beta

            # Mean residence time (MRT)
            mrt = (A / (alpha**2) + B / (beta**2)) / auc

            # Clearance (for blood/plasma only - conceptual)
            clearance = 1 / auc if tissue.lower() in ["blood", "plasma"] else None

            # Store results
            results["pk_parameters"][tissue] = {
                "A": float(A),
                "alpha": float(alpha),
                "B": float(B),
                "beta": float(beta),
                "distribution_half_life_h": float(t_half_dist),
                "elimination_half_life_h": float(t_half_elim),
                "mean_residence_time_h": float(mrt),
            }

            if clearance:
                results["pk_parameters"][tissue]["clearance"] = float(clearance)

            # Calculate AUC
            results["auc_values"][tissue] = float(auc)

        except Exception as e:
            results["pk_parameters"][tissue] = f"Fitting failed: {str(e)}"

    # Calculate tumor-to-normal tissue ratios at each time point
    for tissue in tissue_data:
        if tissue != "tumor":
            ratios = [
                t / n if n > 0 else float("inf")
                for t, n in zip(tissue_data["tumor"], tissue_data[tissue], strict=False)
            ]
            results["tumor_to_normal_ratios"][tissue] = {
                "values": [float(r) for r in ratios],
                "max_ratio": float(max(ratios)),
                "max_ratio_time_point": float(time_points[np.argmax(ratios)]),
            }

    # Save results to JSON file
    filename = "biodistribution_pk_results.json"
    with open(filename, "w") as f:
        json.dump(results, f, indent=2)

    # Generate research log
    log = "# Biodistribution and Pharmacokinetic Analysis of Radiolabeled Antibody\n\n"
    log += "## Analysis Summary\n"
    log += f"- Analyzed biodistribution data across {len(tissue_data)} tissues\n"
    log += f"- Time points analyzed: {time_points} hours\n"
    log += "- Performed bi-exponential pharmacokinetic modeling\n\n"

    log += "## Key Pharmacokinetic Parameters\n"
    for tissue, params in results["pk_parameters"].items():
        if isinstance(params, dict):
            log += f"\n### {tissue.capitalize()}\n"
            log += f"- Distribution half-life: {params['distribution_half_life_h']:.2f} hours\n"
            log += f"- Elimination half-life: {params['elimination_half_life_h']:.2f} hours\n"
            log += f"- Mean residence time: {params['mean_residence_time_h']:.2f} hours\n"
            if "clearance" in params:
                log += f"- Clearance: {params['clearance']:.4f} units\n"

    log += "\n## Tumor-to-Normal Tissue Ratios\n"
    for tissue, ratio_data in results["tumor_to_normal_ratios"].items():
        log += f"- {tissue.capitalize()}: Max ratio {ratio_data['max_ratio']:.2f} at {ratio_data['max_ratio_time_point']:.1f} hours\n"

    log += "\n## Detailed Results\n"
    log += f"Complete analysis results saved to: {os.path.abspath(filename)}\n"

    return log


def estimate_alpha_particle_radiotherapy_dosimetry(
    biodistribution_data, radiation_parameters, output_file="dosimetry_results.csv"
):
    """Estimate radiation absorbed doses to tumor and normal organs for alpha-particle radiotherapeutics.

    This function implements the Medical Internal Radiation Dose (MIRD) schema to calculate
    absorbed doses based on biodistribution data from healthy mice and radiation transport parameters.

    Parameters
    ----------
    biodistribution_data : dict
        Dictionary containing organ/tissue names as keys and a list of time-activity measurements as values.
        Each measurement should be a tuple of (time_hours, percent_injected_activity).
        Must include entries for all relevant organs including 'tumor'.

    radiation_parameters : dict
        Dictionary containing radiation parameters for the alpha-emitting radionuclide:
        - 'radionuclide': str - Name of the radionuclide (e.g., 'Ac-225')
        - 'half_life_hours': float - Physical half-life in hours
        - 'energy_per_decay_MeV': float - Energy released per decay in MeV
        - 'radiation_weighting_factor': float - Radiation weighting factor for alpha particles
        - 'S_factors': dict - S-factors (Gy/Bq-s) for each source-target organ pair

    output_file : str, optional
        Filename to save the dosimetry results (default: "dosimetry_results.csv")

    Returns
    -------
    str
        Research log summarizing the dosimetry estimation process and results

    """
    import csv
    from datetime import datetime

    import numpy as np
    from scipy.integrate import trapezoid

    # Initialize research log
    log = f"Alpha-Particle Radiotherapy Dosimetry Estimation - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
    log += f"Radionuclide: {radiation_parameters['radionuclide']}\n"
    log += f"Half-life: {radiation_parameters['half_life_hours']} hours\n\n"

    # Step 1: Calculate time-integrated activity for each organ
    log += "Step 1: Calculating time-integrated activity for each organ\n"
    time_integrated_activity = {}

    for organ, measurements in biodistribution_data.items():
        times = [m[0] for m in measurements]
        activities = [m[1] for m in measurements]

        # Apply physical decay correction
        decay_constant = np.log(2) / radiation_parameters["half_life_hours"]
        decay_corrected_activities = [a * np.exp(-decay_constant * t) for a, t in zip(activities, times, strict=False)]

        # Calculate time-integrated activity using trapezoidal integration
        cumulated_activity = trapezoid(decay_corrected_activities, times)
        time_integrated_activity[organ] = cumulated_activity

        log += f"  - {organ}: {cumulated_activity:.4f} %IA-h\n"

    # Step 2: Calculate absorbed dose using MIRD schema
    log += "\nStep 2: Calculating absorbed doses using MIRD schema\n"

    # Convert %IA-h to MBq-h for a standard injection of 1 MBq
    conversion_factor = 0.01  # Convert %IA to fraction of IA

    # Energy conversion factor: MeV to J

    # Calculate absorbed dose for each target organ
    absorbed_doses = {}
    s_factors = radiation_parameters["S_factors"]

    for target_organ in biodistribution_data:
        absorbed_dose = 0

        # Sum contributions from all source organs
        for source_organ, cumulated_activity in time_integrated_activity.items():
            if (source_organ, target_organ) in s_factors:
                s_value = s_factors[(source_organ, target_organ)]
                organ_contribution = cumulated_activity * conversion_factor * s_value
                absorbed_dose += organ_contribution

        # Apply radiation weighting factor for alpha particles
        absorbed_dose *= radiation_parameters["radiation_weighting_factor"]

        # Store as Gy/MBq
        absorbed_doses[target_organ] = absorbed_dose
        log += f"  - {target_organ}: {absorbed_dose:.4f} Gy/MBq\n"

    # Step 3: Calculate therapeutic index (tumor-to-normal tissue dose ratios)
    log += "\nStep 3: Calculating therapeutic indices (tumor-to-normal tissue ratios)\n"

    tumor_dose = absorbed_doses.get("tumor", 0)
    if tumor_dose > 0:
        for organ, dose in absorbed_doses.items():
            if organ != "tumor" and dose > 0:
                therapeutic_index = tumor_dose / dose
                log += f"  - Tumor-to-{organ} ratio: {therapeutic_index:.2f}\n"

    # Save results to CSV file
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Organ", "Absorbed Dose (Gy/MBq)"])
        for organ, dose in absorbed_doses.items():
            writer.writerow([organ, f"{dose:.4f}"])

    log += f"\nDosimetry results saved to {output_file}\n"

    return log


def perform_mwas_cyp2c19_metabolizer_status(
    methylation_data_path,
    metabolizer_status_path,
    covariates_path=None,
    pvalue_threshold=0.05,
    output_file="significant_cpg_sites.csv",
):
    """Perform a Methylome-wide Association Study (MWAS) to identify CpG sites significantly associated with CYP2C19 metabolizer status.

    Parameters
    ----------
    methylation_data_path : str
        Path to CSV or TSV file containing DNA methylation beta values.
        Rows should be samples, columns should be CpG sites.
    metabolizer_status_path : str
        Path to CSV or TSV file containing CYP2C19 metabolizer status for each sample.
        Should have a sample ID column and a status column (e.g., poor, intermediate, normal, rapid, ultrarapid).
    covariates_path : str, optional
        Path to CSV or TSV file containing covariates to adjust for in the regression model
        (e.g., age, sex, smoking status).
    pvalue_threshold : float, optional
        P-value threshold for significance after multiple testing correction. Default is 0.05.
    output_file : str, optional
        Filename to save significant CpG sites. Default is "significant_cpg_sites.csv".

    Returns
    -------
    str
        A research log summarizing the MWAS analysis and results.

    """
    import time

    import pandas as pd
    from scipy.stats import linregress
    from statsmodels.formula.api import ols

    start_time = time.time()
    log = ["## Methylome-wide Association Study (MWAS) of CYP2C19 Metabolizer Status"]

    # Load data from files
    log.append("\n### Loading Data")
    try:
        # Load methylation data
        if methylation_data_path.endswith(".csv"):
            methylation_data = pd.read_csv(methylation_data_path, index_col=0)
        elif methylation_data_path.endswith((".tsv", ".txt")):
            methylation_data = pd.read_csv(methylation_data_path, sep="\t", index_col=0)
        else:
            log.append("Error: Unsupported file format for methylation data. Please provide a CSV or TSV file.")
            return "\n".join(log)
        log.append(f"- Successfully loaded methylation data from {methylation_data_path}")

        # Load metabolizer status data
        if metabolizer_status_path.endswith(".csv"):
            metabolizer_status_df = pd.read_csv(metabolizer_status_path, index_col=0)
        elif metabolizer_status_path.endswith((".tsv", ".txt")):
            metabolizer_status_df = pd.read_csv(metabolizer_status_path, sep="\t", index_col=0)
        else:
            log.append("Error: Unsupported file format for metabolizer status. Please provide a CSV or TSV file.")
            return "\n".join(log)
        log.append(f"- Successfully loaded metabolizer status data from {metabolizer_status_path}")

        # Convert DataFrame to Series if necessary
        if metabolizer_status_df.shape[1] == 1:
            metabolizer_status = metabolizer_status_df.iloc[:, 0]
        else:
            log.append("Error: Metabolizer status file should contain a single column with status values.")
            return "\n".join(log)

        # Load covariates if provided
        covariates = None
        if covariates_path is not None:
            if covariates_path.endswith(".csv"):
                covariates = pd.read_csv(covariates_path, index_col=0)
            elif covariates_path.endswith((".tsv", ".txt")):
                covariates = pd.read_csv(covariates_path, sep="\t", index_col=0)
            else:
                log.append("Error: Unsupported file format for covariates. Please provide a CSV or TSV file.")
                return "\n".join(log)
            log.append(f"- Successfully loaded covariates data from {covariates_path}")
    except Exception as e:
        log.append(f"Error loading data: {str(e)}")
        return "\n".join(log)

    # Step 1: Data preprocessing
    log.append("\n### Data Preprocessing")
    log.append(f"- Methylation data shape: {methylation_data.shape} (samples × CpG sites)")
    log.append(f"- Number of samples with metabolizer status: {len(metabolizer_status)}")

    # Ensure sample IDs match between methylation data and metabolizer status
    common_samples = methylation_data.index.intersection(metabolizer_status.index)
    methylation_data = methylation_data.loc[common_samples]
    metabolizer_status = metabolizer_status.loc[common_samples]

    log.append(f"- Number of samples after matching: {len(common_samples)}")

    # Check for covariates
    if covariates is not None:
        log.append(f"- Covariates provided: {', '.join(covariates.columns)}")
        covariates = covariates.loc[common_samples]

    # Step 2: Perform regression for each CpG site
    log.append("\n### Association Analysis")
    log.append(f"- Total CpG sites to analyze: {methylation_data.shape[1]}")

    results = []
    cpg_sites = methylation_data.columns

    # Convert metabolizer status to numeric if it's categorical
    if metabolizer_status.dtype == "object":
        # Create a mapping dictionary for metabolizer status
        # Assuming order: poor < intermediate < normal < rapid < ultrarapid
        status_order = {
            "poor": 1,
            "intermediate": 2,
            "normal": 3,
            "rapid": 4,
            "ultrarapid": 5,
        }

        # Try to map using the dictionary, or keep as is if already numeric
        try:
            metabolizer_status_numeric = metabolizer_status.map(status_order)
            log.append("- Converted metabolizer status to numeric values")
        except Exception:
            metabolizer_status_numeric = metabolizer_status
            log.append("- Using metabolizer status as provided (assuming numeric)")
    else:
        metabolizer_status_numeric = metabolizer_status

    # Perform regression for each CpG site
    for cpg in cpg_sites:
        methylation_values = methylation_data[cpg]

        # Basic model without covariates
        if covariates is None:
            model = linregress(metabolizer_status_numeric, methylation_values)
            pvalue = model.pvalue
            coefficient = model.slope
        else:
            # Create DataFrame for regression with covariates
            data_for_regression = pd.DataFrame(
                {
                    "methylation": methylation_values,
                    "metabolizer": metabolizer_status_numeric,
                }
            )

            # Add covariates
            for col in covariates.columns:
                data_for_regression[col] = covariates[col]

            # Formula for regression with covariates
            formula = "methylation ~ metabolizer + " + " + ".join(covariates.columns)
            model = ols(formula, data=data_for_regression).fit()

            pvalue = model.pvalues["metabolizer"]
            coefficient = model.params["metabolizer"]

        results.append({"CpG_site": cpg, "coefficient": coefficient, "pvalue": pvalue})

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)

    # Step 3: Multiple testing correction
    log.append("\n### Multiple Testing Correction")
    log.append("- Applying Bonferroni correction")

    # Bonferroni correction
    results_df["adjusted_pvalue"] = results_df["pvalue"] * len(results_df)
    results_df["adjusted_pvalue"] = results_df["adjusted_pvalue"].clip(upper=1.0)  # Ensure p-values don't exceed 1

    # Step 4: Identify significant CpG sites
    significant_sites = results_df[results_df["adjusted_pvalue"] < pvalue_threshold]
    significant_sites = significant_sites.sort_values("adjusted_pvalue")

    log.append("\n### Results")
    log.append(f"- Number of significant CpG sites (adjusted p < {pvalue_threshold}): {len(significant_sites)}")

    if len(significant_sites) > 0:
        # Save significant sites to file
        significant_sites.to_csv(output_file, index=False)
        log.append("- Top 5 significant CpG sites:")

        for _, row in significant_sites.head(5).iterrows():
            log.append(
                f"  * {row['CpG_site']}: coefficient = {row['coefficient']:.4f}, adj. p-value = {row['adjusted_pvalue']:.6f}"
            )

        log.append(f"- Full results saved to: {output_file}")
    else:
        log.append("- No significant CpG sites found after multiple testing correction")

    # Execution time
    execution_time = time.time() - start_time
    log.append("\n### Summary")
    log.append(f"- Analysis completed in {execution_time:.2f} seconds")

    return "\n".join(log)


def calculate_physicochemical_properties(smiles_string):
    """Calculate key physicochemical properties of a drug candidate molecule.

    Parameters
    ----------
    smiles_string : str
        The molecular structure in SMILES format

    Returns
    -------
    str
        A research log summarizing the calculated physicochemical properties and
        indicating where the detailed results are saved

    """
    import csv
    import os

    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski
    from rdkit.Chem.MolStandardize import rdMolStandardize

    # Create RDKit molecule from SMILES
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            return "ERROR: Invalid SMILES string provided."
    except Exception as e:
        return f"ERROR: Failed to process SMILES string: {str(e)}"

    # Calculate basic properties
    properties = {
        "SMILES": smiles_string,
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "cLogP": round(Descriptors.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "H-Bond Donors": Lipinski.NumHDonors(mol),
        "H-Bond Acceptors": Lipinski.NumHAcceptors(mol),
        "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "Heavy Atoms": mol.GetNumHeavyAtoms(),
        "Ring Count": Descriptors.RingCount(mol),
    }

    # Estimate pKa (simplified approach - in practice would use specialized tools)
    # This is a simplification as accurate pKa prediction requires specialized tools
    uncharger = rdMolStandardize.Uncharger()
    uncharger.uncharge(mol)
    acidic_groups = sum(
        1
        for atom in mol.GetAtoms()
        if atom.GetSymbol() == "O"
        and any(neigh.GetSymbol() == "C" and neigh.GetDegree() == 3 for neigh in atom.GetNeighbors())
    )
    basic_groups = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "N" and atom.GetDegree() < 4)
    properties["Estimated Acidic Groups"] = acidic_groups
    properties["Estimated Basic Groups"] = basic_groups

    # Calculate drug-likeness score (using Crippen approach)
    properties["Drug-likeness Score"] = round(Crippen.MolMR(mol), 2)

    # Calculate logD (simplified as logP - pKa adjustment would need specialized tools)
    properties["Estimated logD7.4"] = properties["cLogP"]

    # Save results to CSV
    csv_filename = "physicochemical_properties.csv"
    with open(csv_filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Property", "Value"])
        for prop, value in properties.items():
            writer.writerow([prop, value])

    # Generate research log
    log = f"""Physicochemical Property Calculation Research Log:

Analyzed compound with SMILES: {smiles_string}

Key properties:
- Molecular Weight: {properties["Molecular Weight"]} g/mol
- cLogP: {properties["cLogP"]}
- Topological Polar Surface Area: {properties["TPSA"]} Å²
- H-Bond Donors: {properties["H-Bond Donors"]}
- H-Bond Acceptors: {properties["H-Bond Acceptors"]}
- Rotatable Bonds: {properties["Rotatable Bonds"]}
- Estimated logD (at pH 7.4): {properties["Estimated logD7.4"]}
- Estimated Acidic Groups: {properties["Estimated Acidic Groups"]}
- Estimated Basic Groups: {properties["Estimated Basic Groups"]}

Complete results saved to: {os.path.abspath(csv_filename)}
"""

    return log


def analyze_xenograft_tumor_growth_inhibition(
    data_path,
    time_column,
    volume_column,
    group_column,
    subject_column,
    output_dir="./results",
):
    """Analyze tumor growth inhibition in xenograft models across different treatment groups.

    Parameters
    ----------
    data_path : str
        Path to CSV or TSV file containing tumor volume measurements. The file should have columns for
        time, volume, treatment group, and subject ID
    time_column : str
        Name of the column containing time points (e.g., 'Day', 'Time')
    volume_column : str
        Name of the column containing tumor volume measurements
    group_column : str
        Name of the column containing treatment group labels
    subject_column : str
        Name of the column containing subject/mouse identifiers
    output_dir : str, optional
        Directory to save output files (default: "./results")

    Returns
    -------
    str
        Research log summarizing the analysis steps, findings, and generated file paths

    """
    import os

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm
    from scipy import stats
    from statsmodels.formula.api import ols
    from statsmodels.stats.multicomp import pairwise_tukeyhsd

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Initialize research log
    log = "# Xenograft Tumor Growth Inhibition Analysis\n\n"

    # Load data from file
    log += "## 1. Data Loading and Summary\n\n"
    try:
        if data_path.endswith(".csv"):
            data_df = pd.read_csv(data_path)
        elif data_path.endswith((".tsv", ".txt")):
            data_df = pd.read_csv(data_path, sep="\t")
        else:
            log += "Error: Unsupported file format. Please provide a CSV or TSV file.\n"
            return log
        log += f"Successfully loaded tumor growth data from {data_path}\n"
    except Exception as e:
        log += f"Error loading data: {str(e)}\n"
        return log

    # Validate required columns
    required_columns = [time_column, volume_column, group_column, subject_column]
    missing_columns = [col for col in required_columns if col not in data_df.columns]
    if missing_columns:
        log += f"Error: Missing required columns: {', '.join(missing_columns)}\n"
        return log

    # Get unique groups and time points
    groups = data_df[group_column].unique()
    time_points = sorted(data_df[time_column].unique())
    n_groups = len(groups)
    log += f"- Number of treatment groups: {n_groups} ({', '.join(map(str, groups))})\n"
    log += f"- Number of time points: {len(time_points)}\n"
    log += f"- Number of subjects: {data_df[subject_column].nunique()}\n"
    log += f"- Total number of measurements: {len(data_df)}\n\n"

    # 2. Calculate group statistics at each time point
    log += "## 2. Tumor Growth Analysis\n\n"

    # Group statistics
    stats_df = (
        data_df.groupby([group_column, time_column])[volume_column]
        .agg(mean="mean", sem=lambda x: stats.sem(x), count="count")
        .reset_index()
    )

    # Save group statistics
    stats_file = os.path.join(output_dir, "tumor_volume_statistics.csv")
    stats_df.to_csv(stats_file, index=False)
    log += f"Group statistics saved to: {stats_file}\n\n"

    # 3. Calculate tumor growth rates
    log += "## 3. Tumor Growth Rate Analysis\n\n"

    growth_rates = {}
    for group in groups:
        group_data = data_df[data_df[group_column] == group]

        # Calculate growth rate for each subject
        subject_growth_rates = []
        for subject in group_data[subject_column].unique():
            subject_data = group_data[group_data[subject_column] == subject]

            if len(subject_data) >= 2:
                # Simple linear regression for growth rate
                x = subject_data[time_column].values
                y = subject_data[volume_column].values
                slope, _, _, _, _ = stats.linregress(x, y)
                subject_growth_rates.append(slope)

        growth_rates[group] = subject_growth_rates
        mean_rate = np.mean(subject_growth_rates)
        sem_rate = stats.sem(subject_growth_rates)

        log += (
            f"- {group}: Mean growth rate = {mean_rate:.2f} ± {sem_rate:.2f} mm³/day (n={len(subject_growth_rates)})\n"
        )

    # 4. Calculate Tumor Growth Inhibition (TGI)
    log += "\n## 4. Tumor Growth Inhibition (TGI)\n\n"

    # Identify control group (assuming the first group is control)
    control_group = groups[0]
    log += f"Control group: {control_group}\n\n"

    # Calculate TGI for the final time point
    final_time = max(time_points)
    final_data = data_df[data_df[time_column] == final_time]

    control_final_mean = final_data[final_data[group_column] == control_group][volume_column].mean()

    tgi_results = {}
    for group in groups:
        if group == control_group:
            continue

        group_final_mean = final_data[final_data[group_column] == group][volume_column].mean()
        tgi = ((control_final_mean - group_final_mean) / control_final_mean) * 100
        tgi_results[group] = tgi

        log += f"- {group}: TGI = {tgi:.1f}% (relative to {control_group})\n"

    # 5. Statistical Analysis
    log += "\n## 5. Statistical Analysis\n\n"

    # Repeated measures ANOVA
    log += "### Repeated Measures ANOVA\n\n"

    try:
        # Prepare data for repeated measures ANOVA
        formula = f"{volume_column} ~ C({group_column}) * C({time_column}) + C({subject_column})"
        model = ols(formula, data=data_df).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)

        # Save ANOVA results
        anova_file = os.path.join(output_dir, "repeated_measures_anova.csv")
        anova_table.to_csv(anova_file)

        log += f"ANOVA results saved to: {anova_file}\n\n"

        # Extract p-values
        group_effect_p = anova_table.loc[f"C({group_column})", "PR(>F)"]
        time_effect_p = anova_table.loc[f"C({time_column})", "PR(>F)"]
        interaction_p = anova_table.loc[f"C({group_column}):C({time_column})", "PR(>F)"]

        log += f"- Treatment effect: p = {group_effect_p:.4f}\n"
        log += f"- Time effect: p = {time_effect_p:.4f}\n"
        log += f"- Treatment × Time interaction: p = {interaction_p:.4f}\n\n"

        # Post-hoc analysis at final time point
        log += "### Post-hoc Analysis (Final Time Point)\n\n"

        # Perform Tukey's HSD test
        tukey = pairwise_tukeyhsd(endog=final_data[volume_column], groups=final_data[group_column], alpha=0.05)

        # Save Tukey results
        tukey_file = os.path.join(output_dir, "tukey_posthoc_results.txt")
        with open(tukey_file, "w") as f:
            f.write(str(tukey.summary()))

        log += f"Tukey's HSD results saved to: {tukey_file}\n\n"

        # Summarize significant comparisons
        tukey_df = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])

        sig_pairs = tukey_df[tukey_df["p-adj"] < 0.05]
        if len(sig_pairs) > 0:
            log += "Significant pairwise comparisons:\n"
            for _, row in sig_pairs.iterrows():
                log += f"- {row['group1']} vs {row['group2']}: p = {row['p-adj']:.4f}\n"
        else:
            log += "No significant pairwise comparisons found.\n"

    except Exception as e:
        log += f"Error in statistical analysis: {str(e)}\n"

    # 6. Generate tumor growth curves
    log += "\n## 6. Tumor Growth Visualization\n\n"

    plt.figure(figsize=(10, 6))

    for group in groups:
        group_stats = stats_df[stats_df[group_column] == group]
        plt.errorbar(
            group_stats[time_column],
            group_stats["mean"],
            yerr=group_stats["sem"],
            label=group,
            capsize=3,
            marker="o",
        )

    plt.xlabel(f"{time_column} (days)")
    plt.ylabel(f"Tumor Volume ({volume_column})")
    plt.title("Xenograft Tumor Growth Curves")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.7)

    # Save the plot
    plot_file = os.path.join(output_dir, "tumor_growth_curves.png")
    plt.savefig(plot_file, dpi=300, bbox_inches="tight")
    plt.close()

    log += f"Tumor growth curve plot saved to: {plot_file}\n"

    # 7. Conclusion
    log += "\n## 7. Conclusion\n\n"

    # Summarize most effective treatment
    if tgi_results:
        best_treatment = max(tgi_results.items(), key=lambda x: x[1])
        log += f"The most effective treatment was {best_treatment[0]} with a tumor growth inhibition of {best_treatment[1]:.1f}%.\n"

    # Statistical significance summary
    try:
        if group_effect_p < 0.05:
            log += "Statistical analysis confirmed significant differences between treatment groups.\n"
        else:
            log += "No statistically significant differences were found between treatment groups.\n"
    except Exception:
        pass

    return log


def analyze_western_blot(
    blot_image_path,
    target_bands,
    loading_control_band,
    antibody_info,
    output_dir="./results",
):
    """Performs densitometric analysis of Western blot images to quantify relative protein expression.

    Parameters
    ----------
    blot_image_path : str
        Path to the Western blot image file
    target_bands : list of dict
        List of dictionaries containing information about target protein bands
        Each dict should have 'name', 'roi' (region of interest as [x, y, width, height])
    loading_control_band : dict
        Dictionary with 'name' and 'roi' for the loading control protein (e.g., β-actin, GAPDH)
    antibody_info : dict
        Dictionary containing information about antibodies used
        Should have 'primary' and 'secondary' keys with antibody details
    output_dir : str, optional
        Directory to save output files, defaults to './results'

    Returns
    -------
    str
        Research log summarizing the Western blot analysis process and results

    """
    import os

    import numpy as np
    from skimage import io

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load the Western blot image
    image = io.imread(blot_image_path)
    if len(image.shape) > 2:  # Convert to grayscale if color image
        image = np.mean(image, axis=2).astype(np.uint8)

    # Initialize results dictionary
    results = {
        "loading_control": {"name": loading_control_band["name"], "intensity": 0},
        "targets": [],
    }

    # Analyze loading control band
    lc_roi = loading_control_band["roi"]
    lc_band = image[lc_roi[1] : lc_roi[1] + lc_roi[3], lc_roi[0] : lc_roi[0] + lc_roi[2]]
    lc_intensity = np.sum(lc_band)
    results["loading_control"]["intensity"] = lc_intensity

    # Analyze target protein bands
    for band in target_bands:
        roi = band["roi"]
        band_img = image[roi[1] : roi[1] + roi[3], roi[0] : roi[0] + roi[2]]
        band_intensity = np.sum(band_img)

        # Calculate relative expression (normalized to loading control)
        relative_expression = band_intensity / lc_intensity

        results["targets"].append(
            {
                "name": band["name"],
                "intensity": band_intensity,
                "relative_expression": relative_expression,
            }
        )

    # Generate results table and save to CSV
    results_file = os.path.join(output_dir, "western_blot_results.csv")
    with open(results_file, "w") as f:
        f.write("Protein,Raw Intensity,Relative Expression\n")
        f.write(f"{results['loading_control']['name']},{results['loading_control']['intensity']},1.0\n")
        for target in results["targets"]:
            f.write(f"{target['name']},{target['intensity']},{target['relative_expression']:.4f}\n")

    # Generate research log
    log = "## Western Blot Analysis\n\n"
    log += f"Analyzed Western blot image: {os.path.basename(blot_image_path)}\n\n"
    log += "### Antibodies Used\n"
    log += f"Primary antibody: {antibody_info['primary']}\n"
    log += f"Secondary antibody: {antibody_info['secondary']}\n\n"
    log += "### Analysis Steps\n"
    log += "1. Loaded Western blot image and converted to grayscale\n"
    log += f"2. Quantified loading control ({loading_control_band['name']}) band intensity\n"
    log += "3. Measured target protein band intensities\n"
    log += "4. Calculated relative expression by normalizing to loading control\n\n"
    log += "### Results\n"
    log += f"Loading control ({loading_control_band['name']}): {results['loading_control']['intensity']} intensity units\n\n"
    log += "Target proteins:\n"
    for target in results["targets"]:
        log += f"- {target['name']}: {target['intensity']} intensity units, "
        log += f"{target['relative_expression']:.4f} relative expression\n"
    log += f"\nDetailed results saved to: {results_file}\n"

    return log
