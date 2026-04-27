# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains the code for the manuscript: *"Benchmarking 13 tools for mutational signature attribution, including a new and improved algorithm"* (Nanhai Jiang, Yang Wu, Steve Rozen). It benchmarks 13 signature attribution tools — including the new PASA algorithm — on synthetic and real cancer genomic data across three mutation types: SBS (Single Base Substitutions), DBS (Dinucleotide Base Substitutions), and ID (Insertions/Deletions).

The PASA algorithm itself lives in the `mSigAct` R package (v3.0.1-branch), not in this repo.

## How to Run

All scripts assume the working directory is the project root. Scripts are run with `Rscript` or sourced from within R.

**Generate synthetic data:**
```r
Rscript synthetic_data/data_gen_code/generate_SBS_data.R
Rscript synthetic_data/data_gen_code/generate_DBS_data.R
Rscript synthetic_data/data_gen_code/generate_ID_data.R
```

**Run all attribution tools for a mutation type:**
```r
Rscript analysis/code/SBS/run_all_SBS.R   # also sets mut_type env var
Rscript analysis/code/DBS/run_all_DBS.R
Rscript analysis/code/ID/run_all_ID.R
```

**Run a single tool (example):**
```r
Rscript analysis/code/SBS/run_pasa_SBS.R
```

**MSA tool runs via Nextflow/shell (not R):**
```bash
bash analysis/code/SBS/run_msa_syn_SBS.sh
```

**Gather statistics and generate paper outputs:**
```r
Rscript analysis/code/gather_all_stats_and_cpu_time.R
Rscript analysis/code/rank_tools.R
Rscript analysis/code/plot_main_text_and_sup_figs_all_cancer_types_combined.R
```

## Architecture

### Core execution pattern

`generic_analysis.R::run_generic_syn()` is the central function. Every tool (except MSA) plugs into it:

1. `get_all_input.R` loads synthetic spectra, ground-truth signatures, and per-cancer-type signature universes from `synthetic_data/<TYPE>/`.
2. For each cancer type, `run_non_msi_and_msi()` splits samples into MSI-H and non-MSI cohorts, then calls the tool's `attribute_function`.
3. Results are written to `analysis/raw_output/<type>/<tool>/inferred_exposures.csv`.

Each tool has two files:
- `<tool>_analysis.R` — wraps the tool's API in `call_<tool>()`, then calls `run_generic_syn()` via `run_<tool>()`
- `analysis/code/<TYPE>/run_<tool>_<TYPE>.R` — entry point that sources the analysis file and invokes `run_<tool>()`

### Key globals (`analysis/code/common_utils.R`)

- `global_random_seed = 145879`
- `global_raw_tools_to_plot` — canonical internal names for all 13 tools
- `pretty_tool_names()` — maps internal → display names
- `global_measures` — performance metrics: `Combined`, `one_minus_smd`, `prec`, `sens`, `spec`, `scaled_L2`, `KL`
- `plot_output_directory` / `global_output_for_paper` = `"output_for_paper"`

### MSI handling

Stomach and colorectal cancer datasets contain MSI-H samples (column names contain `"MSI"`). These are automatically split from non-MSI samples in `run_non_msi_and_msi()`. MSI signatures are `SBS6, SBS14, SBS15, SBS20, SBS21, SBS26, SBS44`.

### FitMS variants

FitMS is tested at 8 rare-signature thresholds (`fitms_0.001` through `fitms_0.200`). The threshold `fitms_0.010` is treated as the canonical FitMS result in comparisons.

### Python-based tools

MuSiCal and SigPro have Python wrappers (`run_musical.py`, `run_sigpro.py`) and require separate conda environments (`conda_env_for_musical.yml`).

### Output locations

- `analysis/raw_output/<TYPE>/<tool>/` — raw per-tool results (`inferred_exposures.csv`, `time_used.Rds`)
- `output_for_paper/` — aggregated CSVs, PDF figures, XLSX supplementary tables

## Key R Packages

- `mSigAct` (v3.0.1-branch) — PASA algorithm (`PresenceAttributeSigActivity`)
- `mSigTools` — exposure I/O and utilities
- `ICAMS` — mutation catalog handling
- `SynSigGen` — synthetic data generation
- `cosmicsig` / `PCAWG7` — reference signatures and PCAWG data
- `signature.tools.lib` (v2.4.5) — FitMS
