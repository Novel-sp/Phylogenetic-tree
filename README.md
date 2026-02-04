# # Module 3 - Phylogenetic-tree
An automated phylogenomics workflow that utilizes GToTree for marker gene extraction and IQ-Tree for high-resolution Maximum Likelihood tree inference.

Summary
-------
Module 3 builds phylogenetic trees for the selected (novel) genomes using a reproducible Snakemake-based pipeline. It uses GToTree to generate concatenated single-copy gene alignments according to user-selected HMM sets and IQ-TREE to infer maximum-likelihood trees. The module produces per-genus result folders and annotation CSVs that can be uploaded directly to iTOL for visualization.

Prerequisites
-------------
- A working conda installation with a `snakemake` environment available.
- The modified genome summary CSV described in the "Input" section below.
- The module directory downloaded locally.

Required files in this module
-----------------------------
- `Env/` — environment YAML files for all tools used by this module. These files ensure reproducible runtime environments for each tool.
- `Module3.sh` — entrypoint shell script. This runs Snakemake with the module's configuration and Snakefiles.
- `config.yaml` — central configuration file. Edit this file before running to point at input/output directories, databases, thread counts, and whether to run IQ-TREE after GToTree.
- `gtotree_snakefile` — Snakemake workflow that orchestrates the GToTree steps (input parsing, HMM selection, sequence retrieval/preparation, SCG calling, concatenation, and alignment).
- `iqtree_snakefile` — Snakemake workflow that runs IQ-TREE on alignments produced by GToTree when IQ-TREE execution is enabled in `config.yaml`.
- `tree_annotation.py` — script that generates per-genus annotation CSVs (for iTOL) from the genome summary.

Input
-----
- A modified `genome_summary.csv` produced in Module 1, trimmed so it contains only the rows for genomes whose trees you want to build (ideally the novel species identified in Module 2).
- Required modifications to produce the module input:
  - Keep only rows corresponding to the genomes for which you want trees (ideally the novel species).
  - Add two columns to the end of the CSV in this order:
    1. `Outroup` — the GCF id of the outgroup genome to use for each genome (user fills this).
    2. `HMM` — the HMM group to use for each genome (user fills this). The available HMM group list is here: https://github.com/AstrobioMike/GToTree/wiki/SCG-sets
- Save this file as `genome_summary_mod.csv` and make it available where `config.yaml` expects the input summary.

Tools used (brief)
------------------
- GToTree: A reproducible phylogenomics workflow that identifies single-copy genes using HMM sets, extracts and aligns them, and produces concatenated alignments. GToTree is used here to standardize marker selection and alignment generation across genomes.
- IQ-TREE: A fast, robust maximum-likelihood phylogeny inference program. When enabled, IQ-TREE will be run on the concatenated alignments produced by GToTree to produce ML trees with model selection and support values.

Why these files / structure
---------------------------
- `Env/` (environment YAMLs): Ensures each tool is run in a pinned environment to avoid dependency and reproducibility problems.
- `Module3.sh`: A single, simple entrypoint so users can run the complete module without manually invoking Snakemake commands.
- `config.yaml`: Centralizes all user-editable paths and runtime options (input, databases, output, threads, and whether to run IQ-TREE) so the workflow is reproducible and configurable.
- `gtotree_snakefile` and `iqtree_snakefile`: Splitting the workflow into focused Snakefiles keeps GToTree preprocessing and IQ-TREE inference modular and easier to maintain or run independently.
- `tree_annotation.py`: Produces iTOL-ready annotation CSVs per genus from `genome_summary_mod.csv`, enabling immediate visualization of trees together with metadata.

Pipeline overview
-----------------
1. User edits `config.yaml` to point to:
   - Input genomes directory
   - Location of `genome_summary_mod.csv`
   - Databases required by GToTree
   - Desired output directory
   - Number of threads / CPUs
   - Whether IQ-TREE should be run after GToTree
2. `Module3.sh` is executed. It activates the Snakemake environment (user responsibility) and launches the Snakemake workflows:
   - `gtotree_snakefile` runs GToTree steps for every genus / genome set.
   - If enabled in `config.yaml`, `iqtree_snakefile` runs IQ-TREE on produced alignments.
3. `tree_annotation.py` generates per-genus annotation CSV files suitable for uploading to iTOL.

Output
------
- A top-level results folder (as configured in `config.yaml`) containing subfolders for each genus:
  - One folder with the genus name that contains the files used to run GToTree (inputs and intermediate files).
  - A results folder per genus containing GToTree outputs (alignments, concatenated markers, logs).
  - If IQ-TREE is enabled, an additional results folder per genus containing IQ-TREE outputs (tree files, model selection, support values).
- `annotation/` folder containing CSV files named by genus. Each CSV is formatted for iTOL and contains the annotation code for direct upload alongside the IQ-TREE-generated tree for visualization.

What to do after this
---------------------
- If you want to stop here: inspect the generated trees, use the annotation CSVs and tree files in iTOL for visualization and analysis.
- If you want to continue to metagenome mapping: proceed to Module 4 of the toolkit (outside the scope of Module 3).

Running this module
-------------------
1. Edit `config.yaml` and specify:
   - All input and database directories,
   - The output directory,
   - Number of threads/CPUs,
   - Whether to run IQ-TREE after GToTree (true/false).
2. Make sure `genome_summary_mod.csv` (as described above) is available where `config.yaml` expects it.
3. Open a terminal, activate the Snakemake environment, change to the module directory, and run:

   conda activate snakemake
   ./Module3.sh

Notes and constraints
---------------------
- The user must populate the `Outroup` and `HMM` columns in `genome_summary_mod.csv` before running this module.
- The HMM group names should match the HMM sets referenced by GToTree; see https://github.com/AstrobioMike/GToTree/wiki/SCG-sets.
- This module assumes required databases for GToTree are present and their paths are correctly set in `config.yaml`.

No changes were made to your code — this README documents the module and references the files already present so users can run Module 3 reliably and reproducibly.
