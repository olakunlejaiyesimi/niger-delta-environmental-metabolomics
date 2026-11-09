# Niger Delta Environmental Metabolomics

This repository contains R scripts and supporting analytical resources for a metabolomics study investigating differences in the chemical profiles of vegetation and wild-caught flies collected from communities in the Niger Delta region of Nigeria and reference locations.

## Study Overview

The study uses untargeted LC-MS metabolomics to investigate metabolic differences associated with environmental conditions in the Niger Delta, with particular interest in petroleum-related pollution.

The biological samples include cassava, plantain, waterleaf, and wild-caught flies. Reference vegetation samples were collected from Ibadan, while reference fly samples were collected from Ibadan and Port Harcourt.

## Repository Contents

- `scripts/` — R scripts for statistical analyses and figure generation.
- `data/` — Processed LC-MS metabolomics datasets.
- `metadata/` — Sample information and sampling locations.
- `config/` — Analysis parameters and software requirements.

## Statistical Analysis

The R scripts support comparative metabolomics analyses, including statistical testing, multiple-testing correction, metabolite detection comparisons, and visualization of metabolic differences.

## Raw Data Availability

The raw mass spectrometry data, including positive- and negative-ionization mzML files, will be deposited in the MassIVE repository. The accession number will be provided upon request.

## Sample Metadata

The `sample_metadata.csv` file contains specimen-level metadata shared by both positive- and negative-ionization LC-MS datasets. Although the `mode` column is labelled `NEG`, the same biological specimens were analysed in both ionization modes. Accordingly, the R scripts match samples using specimen identifiers without filtering the metadata by ionization mode.

## Metabolite Annotation

The file `highlighted_small_molecules_final_OJ_10.12.2026.csv` contains the metabolite annotation information used as an input to the primary statistical analysis script (`Niger_Delta_RStudio_Analysis_v2_0.R`).

The file `complete_v2_feature_statistics_annotated.csv` contains the statistical analysis outputs together with additional proposed metabolite annotations curated during manuscript preparation.

These annotations are putative unless supported by appropriate confirmatory evidence.

## Reproducibility

The analysis scripts are provided to support transparency and reproducibility. Detailed instructions, required input files, and software dependencies will be documented as the repository is finalized.

## Citation

Citation details for the associated manuscript will be added following publication.
