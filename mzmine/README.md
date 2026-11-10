# MZmine Data Processing

The LC-MS data associated with the Niger Delta environmental metabolomics study were processed using **MZmine version 4.50**.

The batch files provided in this directory document the processing workflows used for the positive- and negative-ionization datasets.

These files contain the processing parameters and module configurations used during feature detection, chromatographic processing, feature alignment, filtering, and feature-table generation, as applicable to each workflow.

## Batch Files

- **Positive-ionization mode:** See the corresponding positive-mode MZmine batch file in this directory.
- **Negative-ionization mode:** See the corresponding negative-mode MZmine batch file in this directory.

## Reproducibility

To reproduce the processing workflow:

1. Download and install MZmine version 4.50.
2. Obtain the corresponding raw LC-MS mzML files from MassIVE, accession **MSV000103910**.
3. Open MZmine and load the appropriate batch file.
4. Update the raw-data file paths and output directories to match the local computing environment.
5. Execute the batch workflow for the corresponding ionization mode.

The resulting feature tables were used for downstream statistical analyses in R.

## Related Resources

**Raw LC-MS data:** MassIVE accession MSV000103910.

**Processed feature matrices and statistical analyses:** See the `data/` and `scripts/` directories in the parent GitHub repository.

**Note:** Successful execution may require adjustments to local input/output paths and appropriate computational resources.
