#!/usr/bin/env python3
###############################################################################
# Script: Generate annotation files for phylogenetic tree visualization
#
# This script:
# 1. Reads workflow configuration from config.yaml
# 2. Fetches organism and strain metadata from NCBI (if needed)
# 3. Generates annotation CSV files for tree visualization
# 4. Combines:
#    - Ingroup genomes
#    - One optional outgroup genome
#    - Remaining reference genomes (if present)
#
# Output annotation files are written to:
#   <results_dir>/Annotation/<genus>_annotation.csv
###############################################################################

import os
import sys
import subprocess
import yaml
import pandas as pd

# ---------------------------
# Ensure dependencies
# ---------------------------
# Check whether Biopython is installed.
# If not present, automatically install required Python packages.
# This allows the script to run in clean environments.
try:
    import Bio
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython", "pandas", "pyyaml"])
    import Bio

from Bio import Entrez

# ---------------------------
# Load config
# ---------------------------
# Read configuration values from config.yaml.
# This file must define paths such as:
# - results_dir
# - genome_summary
# - (optionally) ncbi_api_key
with open("config.yaml") as f:
    config = yaml.safe_load(f)

results_dir = config["results_dir"]                 # Base workflow results directory
genome_summary_file = config["genome_summary"]      # Genome summary CSV file
ncbi_api_key = config.get("ncbi_api_key", None)     # Optional NCBI API key

# Configure NCBI Entrez access
# An email address is REQUIRED by NCBI for all Entrez queries
Entrez.api_key = ncbi_api_key
Entrez.email = "your_email@example.com"  # REQUIRED by NCBI

# ---------------------------
# Fetch organism + strain from NCBI
# ---------------------------
def fetch_ncbi_metadata(accession):
    """
    Fetch organism name and strain information from NCBI
    for a given genome accession (GCF or GCA).

    Returns:
        organism (str): Organism name
        strain (str): Strain name (if available)
    """
    try:
        # Search for the assembly record using the accession
        handle = Entrez.esearch(db="assembly", term=accession, retmode="xml")
        record = Entrez.read(handle)
        handle.close()

        # If no record is found, return accession as fallback
        if not record["IdList"]:
            print(f"No NCBI record found for {accession}")
            return accession, ""

        # Retrieve the full assembly summary
        uid = record["IdList"][0]
        handle = Entrez.esummary(db="assembly", id=uid, report="full")
        summary = Entrez.read(handle)
        handle.close()

        # Extract organism and strain information
        docsum = summary["DocumentSummarySet"]["DocumentSummary"][0]

        organism = docsum.get("Organism", accession)
        strain = docsum.get("InfraspecificName", "") or ""
        strain = strain.replace("strain=", "").strip()

        return organism, strain

    except Exception as e:
        # Fail gracefully if any NCBI query errors occur
        print(f"Error fetching {accession}: {e}")
        return accession, ""

# ---------------------------
# Generate ingroup + outgroup + remaining genomes annotation
# ---------------------------
def generate_annotation(input_tsv, output_tsv, outgroup_id=None):
    """
    Generate a tree annotation CSV file for a single genus.

    Inputs:
    - input_tsv: NCBI_genomes_summary_info.tsv from GToTree
    - output_tsv: Output annotation CSV path
    - outgroup_id: Optional outgroup accession (GCF/GCA)

    The output file contains formatting columns used by tree visualization tools.
    """

    # Read ingroup genome summary produced by GToTree
    df = pd.read_csv(input_tsv, sep="\t")

    # Keep only required columns and rename them for clarity
    df = df[["input_accession", "organism_name", "infraspecific_name"]]
    df = df.rename(columns={
        "input_accession": "NCBI Accession",
        "organism_name": "genome",
        "infraspecific_name": "Strain"
    })

    # Clean strain names
    df["Strain"] = df["Strain"].str.replace("strain=", "", regex=False).str.strip()

    # Construct a combined Genome identifier column
    df["Genome"] = (
        df["NCBI Accession"].astype(str) + "   " +
        df["genome"].astype(str) + "   " +
        df["Strain"].astype(str)
    ).str.replace(r"\s+", "_", regex=True)

    # Convert GCF accession to GCA for display labels
    df["Accession_GCA"] = df["NCBI Accession"].str.replace("GCF_", "GCA_", regex=False)

    # Create formatted label strings for tree visualization
    df["Label"] = (
        "<i>" + df["genome"] + "</i><n> " + df["Strain"] + "</n>"
        + "<sup size='0.5'>T</sup><n> (" + df["Accession_GCA"] + ")<n>"
    )

    # Remove helper column
    df = df.drop(columns=["Accession_GCA"])

    # Default visualization styling for ingroup genomes
    df["Pos"] = -1
    df["Color"] = "#000000"
    df["Style"] = "normal"
    df["Size"] = 1
    df["Rotation"] = 0

    # ------------------------
    # Add outgroup genome
    # ------------------------
    if outgroup_id:
        # Fetch metadata for outgroup from NCBI
        organism, strain = fetch_ncbi_metadata(outgroup_id)
        gca_id = outgroup_id.replace("GCF_", "GCA_")

        # Construct outgroup annotation row
        out_row = {
            "NCBI Accession": outgroup_id,
            "genome": organism,
            "Strain": strain,
            "Genome": f"{outgroup_id}   {organism}   {strain}".replace(" ", "_"),
            "Label": f"<i>{organism}</i><n> {strain}</n><sup size='0.5'>T</sup><n> ({gca_id})<n>",
            "Pos": -1,
            "Color": "#0000ff",   # Blue color for outgroup
            "Style": "bold",
            "Size": 1,
            "Rotation": 0
        }

        # Append outgroup row
        df = pd.concat([df, pd.DataFrame([out_row])], ignore_index=True)

    # ------------------------
    # Add remaining reference genomes (if present)
    # ------------------------
    # These come from Fasta_genomes_summary_info.tsv
    # Genome and strain names are left blank intentionally
    fasta_summary_path = os.path.join(
        os.path.dirname(input_tsv),
        "Fasta_genomes_summary_info.tsv"
    )

    if os.path.exists(fasta_summary_path):
        fasta_df = pd.read_csv(fasta_summary_path, sep="\t")

        # Exclude outgroup accession if present
        remaining = fasta_df[
            ~fasta_df["Assembly_name"].isin([outgroup_id] if outgroup_id else [])
        ]

        for _, row in remaining.iterrows():
            assembly_id = row["Assembly_name"]
            gca_id = assembly_id.replace("GCF_", "GCA_")

            rem_row = {
                "NCBI Accession": assembly_id,
                "genome": "",   # intentionally blank
                "Strain": "",   # intentionally blank
                "Genome": assembly_id,
                "Label": f"<i></i><n> </n><sup size='0.5'>T</sup><n> ({gca_id})<n>",
                "Pos": -1,
                "Color": "#ff0000",   # Red color for remaining genomes
                "Style": "bold",
                "Size": 1,
                "Rotation": 0
            }

            # Append remaining genome row
            df = pd.concat([df, pd.DataFrame([rem_row])], ignore_index=True)

    # Write annotation CSV file
    df.to_csv(output_tsv, index=False)

# ---------------------------
# Main execution loop
# ---------------------------
def main():
    """
    Main driver function.

    Iterates through all genus-level result directories,
    identifies valid GToTree outputs, and generates
    corresponding annotation CSV files.
    """

    # Load genome summary to retrieve outgroup information
    genome_summary = pd.read_csv(genome_summary_file)

    # Create annotation output directory
    annotation_dir = os.path.join(results_dir, "Annotation")
    os.makedirs(annotation_dir, exist_ok=True)

    # Iterate over all folders in results directory
    for folder in os.listdir(results_dir):
        genus_path = os.path.join(results_dir, folder)

        # Expected GToTree summary file path
        run_files_path = os.path.join(
            genus_path,
            "run_files",
            "NCBI_genomes_summary_info.tsv"
        )

        # Process only valid genus directories
        if os.path.isdir(genus_path) and os.path.exists(run_files_path):
            genus = folder.replace("_GToTree", "")

            # Retrieve outgroup for this genus
            outgroup_id = genome_summary.loc[
                genome_summary["Genus"] == genus,
                "Outgroup"
            ].values[0]

            # Handle missing outgroup values
            if pd.isna(outgroup_id):
                outgroup_id = None

            # Define output annotation file path
            output_tsv = os.path.join(
                annotation_dir,
                f"{genus}_annotation.csv"
            )

            print(f"Processing genus: {genus}, Outgroup: {outgroup_id}")

            # Generate annotation file
            generate_annotation(run_files_path, output_tsv, outgroup_id)

            print(f"[OK] Annotation file written to: {output_tsv}")

# ---------------------------
# Script entry point
# ---------------------------
if __name__ == "__main__":
    main()

