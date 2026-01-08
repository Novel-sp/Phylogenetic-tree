#!/bin/bash
###############################################################################
# Master workflow launcher for phylogenetic tree construction
#
# This script runs the complete phylogenomic pipeline in sequence:
# 1. Input preparation and genome grouping
# 2. GToTree execution (marker gene extraction + alignment)
# 3. IQ-TREE execution (phylogenetic inference)
# 4. Tree annotation file generation
#
# USAGE:
#   bash run_phylogenetic_workflow.sh
#
# REQUIREMENTS:
# - Snakemake installed and available in PATH
# - Conda/Mamba properly configured
# - config.yaml present in the working directory
###############################################################################

# Exit immediately if:
# - Any command fails (-e)
# - Any unset variable is used (-u)
# - Any command in a pipeline fails (-o pipefail)
set -euo pipefail

# ---------------------------------------------------------------------------
# Workflow start message
# ---------------------------------------------------------------------------
echo "=== Starting overall Phylogenetic Tree workflow ==="

# ---------------------------------------------------------------------------
# Step 1: Prepare inputs for GToTree
# ---------------------------------------------------------------------------
# This step:
# - Reads genome_summary
# - Groups genomes by genus
# - Copies genome FASTA files
# - Downloads outgroups
# - Creates FASTA and GCF lists per genus
echo "=== Step 1: Running pre_gtotree_snakefile ==="
snakemake --snakefile pre_gtotree_snakefile --use-conda --cores 20

# ---------------------------------------------------------------------------
# Step 2: Run GToTree
# ---------------------------------------------------------------------------
# This step:
# - Extracts marker genes
# - Aligns single-copy genes
# - Prepares inputs required for IQ-TREE
echo "=== Step 2: Running gtotree_snakefile ==="
snakemake --snakefile gtotree_snakefile --use-conda --cores 20

# ---------------------------------------------------------------------------
# Step 3: Run IQ-TREE
# ---------------------------------------------------------------------------
# This step:
# - Uses GToTree alignments and partitions
# - Selects best-fit models
# - Performs bootstrap-supported phylogenetic inference
echo "=== Step 3: Running iqtree_snakefile ==="
snakemake --snakefile iqtree_snakefile --use-conda --cores 20

# ---------------------------------------------------------------------------
# Step 4: Generate tree annotation files
# ---------------------------------------------------------------------------
# This step:
# - Collects ingroup, outgroup, and remaining genomes
# - Fetches metadata from NCBI when needed
# - Produces annotation CSV files for tree visualization tools
echo "=== Step 4: Running tree_annotation.py ==="
python3 tree_annotation.py

# ---------------------------------------------------------------------------
# Workflow completion message
# ---------------------------------------------------------------------------
echo "=== All steps completed successfully! ==="

