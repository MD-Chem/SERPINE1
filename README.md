# SERPINE1 / RESIST-M

R analysis scripts accompanying **[Modeling oxaliplatin resistance in colorectal cancer reveals a SERPINE1-based gene signature (RESIST-M) and therapeutic strategies for pro-metastatic CMS4 subtype](https://doi.org/10.1038/s41419-025-07855-y)**, published in *Cell Death & Disease* **16**, 529 (2025).

[Read the paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12264272/) · [Single-cell dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299427) · [Report a code or documentation issue](https://github.com/MD-Chem/SERPINE1/issues)

## Research context

The study examines acquired oxaliplatin resistance and metastatic behavior in colorectal cancer. Transcriptomic and experimental analyses connect cholesterol biosynthesis, TGF-β signaling, and SERPINE1, and identify the nine-gene RESIST-M signature, evaluated across patient cohorts.

These scripts cover downstream single-cell analyses and cohort-level plots, starting from processed objects or expression and clinical tables. They do not include raw-read alignment or the complete study workflow.

## Scripts

| File | Analyses included | Outputs named in the script |
| --- | --- | --- |
| [`PAIscRNAseq_MD_07052024.R`](PAIscRNAseq_MD_07052024.R) | HCT116 single-cell quality-metric plots; cluster markers; cell-state annotation; UMAP, feature, violin, and dot plots; cell-state proportions by dose; UCell scores for the RCC and RPS signatures defined in the script | PNG/PDF plots, including `FP_RESISTM.png`, `cellprop_bydose.png`, and `RCC_RPS_UCell.pdf`; marker tables `hct116markers.csv` and `bb.csv`; an updated RDS object |
| [`Figure5.R`](Figure5.R) | TCGA-COADREAD signature comparisons across CMS groups; Kaplan–Meier plots and Cox models; an annotated gene-expression heatmap for the cohort labeled “Nant” in the script | `RESIST_M1_boxplot_25032025.pdf`, `rcc_All CMS.pdf`, and `Overall.pdf` |

See **[published Figure 5 in the paper](https://doi.org/10.1038/s41419-025-07855-y)** for representative cohort analyses. Generated plots are not bundled here. Script comments map box plots to 5A/B/E/F, the heatmap to 5G, and survival plots to 5H–M; each requires selecting the relevant signature/settings. PETACC-3 analyses for 5C/D were performed by a collaborator and are not implemented here.

## Data and required inputs

Datasets are stored separately. The paper lists [single-cell RNA-seq: GSE299427](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299427), [bulk RNA-seq: GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428), and [ATAC-seq: GSE298524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298524). The latter two accessions are not directly loaded by these scripts.

| Script / section | Expected local input | Contents used |
| --- | --- | --- |
| Single-cell | `HCT116.rds` | Processed Seurat object; the script uses cluster identities, quality metrics, gene expression, and stored reductions |
| Single-cell | `hct116_md.rds` | Annotated Seurat object; the script uses `Dose`, `seurat_clusters`, and a stored UMAP reduction |
| Figure5 / TCGA | `fromxena_survival_COADREAD_survival.txt` | Survival metadata, including patient/sample identifiers, redaction status, and survival endpoints |
| Figure5 / TCGA | `fromsynapse_clinical_molecular_public_all.txt` | Clinical/molecular annotations joined by `sample`, including CMS labels |
| Figure5 / TCGA | `HiSeqV2.txt` | Gene-expression table with genes in rows; the first column is named `sample` and supplies the row names, and the remaining columns identify patient samples |
| Figure5 / heatmap | `SG-BULK_patient_clinical_information.csv` | `patient_id` and annotations `iCMS`, `CMS`, `TGFBR2`, `group3`, and `group5` |
| Figure5 / heatmap | `SG-BULK_salmonTPM.csv` | Expression matrix with an Ensembl-ID first column and patient columns |
| Figure5 / heatmap | `features.tsv` | Headerless, three-column feature table; the first two columns provide Ensembl IDs and gene symbols |

**Processed single-cell objects:** the [GSM9039076 record within GSE299427](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9039076) provides `GSM9039076_HCT116.rds` and `GSM9039076_hct116_md.rds`, each listed at approximately 2.5 GB. Update the script's input paths to those filenames, or rename local copies to the names expected above. Retain an untouched copy of each download.

For patient-cohort inputs, consult the paper's methods and dataset references. Exact download links and preparation steps for these filenames are not recorded in the scripts. Check identifiers, annotations, and expression format before running.

## R packages

Explicit imports are listed below; not every package is used in every section.

| Script | Imported packages |
| --- | --- |
| `PAIscRNAseq_MD_07052024.R` | `Seurat`, `stringr`, `ggplot2`, `tidyverse`, `dplyr`, `patchwork`, `sctransform`, `scDblFinder`, `DoubletFinder`, `UCell`, `CytoTRACE2` |
| `Figure5.R` | `tidyverse`, `circlize`, `ComplexHeatmap`, `ggplot2`, `ggpubr`, `survival`, `survminer` |

Install packages according to their maintainers' instructions. A tested R/package-version manifest or lockfile is not included. Record `sessionInfo()` when reproducing an analysis.

## Preparing and running an analysis

The files are interactive research scripts and require local preparation before execution.

1. Download the repository and the inputs for the analysis you need.
2. Open a separate R session, install/load the imported packages, and work with copies of the data.
3. Replace the Windows-specific `setwd()`, `readRDS()`, and `saveRDS()` paths with your local locations. Choose output filenames before running.
4. Read the notes below, then run sections interactively, checking identities, metadata, signature, and cohort.
5. Compare results with the paper; retain your settings and `sessionInfo()`.

### Execution notes

- **Session initialization:** both scripts call `rm(list = ls())`, which clears objects in the current R workspace, and call `dev.off()` before opening a plot device. In a fresh session with no active device, omit that initial call or guard it with `if (dev.cur() > 1) dev.off()`.
- **Single-cell object order:** `bb` is used in plotting and marker calls before the later `bb <- RenameIdents(...)` block creates it. Review the intended cluster annotation and run that preparation before sections using `bb`. The five annotation labels must match the input object's cluster levels.
- **Preserve input data:** the single-cell `saveRDS()` call writes to the same `hct116_md.rds` path used as an input. Change the output path to retain the original object.
- **Selected survival analysis:** the supplied survival block uses the `rps` gene set, the full `merged_data` cohort, and overall survival. Its output filename is `rcc_All CMS.pdf`. Check the signature, cohort filter, endpoint, labels, and filename together when selecting another analysis.
- **Interactive plotting:** several plots are displayed through bare object expressions. When adapting the scripts for `source()` or batch execution, use explicit rendering calls where needed, such as `print()` for ggplot/ggsurvplot objects and `ComplexHeatmap::draw()` for heatmaps.
- **Exploratory block:** the single-cell script includes a standalone `ggplot(hct116_md, ...)` call referencing `p.adjust` and `Count` without constructing the corresponding plotting table. That block needs a suitable data frame before use.

Use these scripts section by section with the inputs and checks above; they are not configured as an unattended, end-to-end pipeline.

## Citation

Wong, S. Q. R., Das, M., Tenzin, K., et al. **Modeling oxaliplatin resistance in colorectal cancer reveals a SERPINE1-based gene signature (RESIST-M) and therapeutic strategies for pro-metastatic CMS4 subtype.** *Cell Death & Disease* **16**, 529 (2025). https://doi.org/10.1038/s41419-025-07855-y

