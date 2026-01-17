# psoriasis-colon-scRNA-MR
psoriasis-colon-scRNA-MR
# Integrative scRNA-seq + Mendelian Randomization framework for psoriasis and colon cancer

This repository contains the analysis code and derived results for an integrative single-cell transcriptomics and two-sample Mendelian randomization (scRNA-seq + MR) study that prioritizes shared immune genes between psoriasis and colon cancer and maps them to disease-relevant immune cell states.

## Study overview

We first identified shared candidate genes from scRNA-seq overlap between psoriasis and colon cancer within shared immune compartments (T cells, dendritic cells, and monocytes). We then used two-sample Mendelian randomization (MR) with cis-eQTL instruments (eQTLGen via IEU OpenGWAS) to evaluate whether genetically proxied gene expression is associated with:
- Colon adenocarcinoma (FinnGen endpoint: `finn-b-C3_COLON_ADENO_EXALLC`)
- Psoriasis (GWAS Catalog / IEU OpenGWAS study ID: `ebi-a-GCST90019017`)

A Steiger directionality test was applied when the outcome dataset provided the allele frequency fields required for variance-explained calculations.

## Data sources (public)

- scRNA-seq (GEO):
  - Colon cancer: `GSE231559`
  - Psoriasis: `GSE151177`
- eQTL instruments:
  - eQTLGen (whole-blood cis-eQTLs), accessed via IEU OpenGWAS (gene-specific datasets: `eqtl-a-ENSG*`)
- Outcome GWAS:
  - Colon adenocarcinoma: FinnGen endpoint `finn-b-C3_COLON_ADENO_EXALLC`
  - Psoriasis: GWAS Catalog / IEU OpenGWAS `ebi-a-GCST90019017`

> Note: This repository does not redistribute large third-party datasets (e.g., raw sequencing matrices or GWAS VCF files). All third-party resources are accessed via the identifiers above.

## Repository structure

- `code/`  
  Analysis scripts for scRNA-seq processing, candidate selection, MR, Steiger directionality testing, sensitivity analyses, and visualization.

- `data/derived/`  
  Small derived files needed to reproduce key results and figures (e.g., candidate gene lists and MR summary tables).  
  **Large raw files (e.g., `.vcf.gz`) are not stored here.**

- `outputs/` (optional)  
  Generated plots and intermediate outputs produced by the scripts below.

## Analysis workflow (high level)

### Step 1. Cross-disease candidate gene definition (415 shared genes)
We first identified cross-disease candidate genes by integrating scRNA-seq results from colon cancer and psoriasis. This step yields **415 shared candidate genes** for downstream MR analyses.

**Key output**
- `data/derived/IntersectionGenes_Combined.csv`  
  A list of cross-disease candidate genes used as exposures in MR.

### Step 2. Two-sample MR in two outcomes (colon cancer and psoriasis) with Steiger testing when available
We performed two-sample MR using eQTLGen cis-eQTL instruments as exposures and two independent outcomes:
- Colon adenocarcinoma (FinnGen)
- Psoriasis (OpenGWAS / GWAS Catalog)

Steiger directionality testing was performed when outcome allele frequency information was available; when unavailable, Steiger results are reported as not available.

**Script(s)**
- `code/03_MR_pipeline_colon.R`  
  MR analysis for colon outcome (`finn-b-C3_COLON_ADENO_EXALLC`) including Steiger directionality testing (when available) and sensitivity analyses.
- `code/04_MR_pipeline_psoriasis.R`  
  MR analysis for psoriasis outcome (`ebi-a-GCST90019017`) including sensitivity analyses; Steiger directionality may be limited depending on outcome data fields.

**Key outputs**
- `data/derived/all_odds_ratios.csv`  
  Main MR results (OR and confidence intervals) across genes and MR methods.
- `data/derived/all_heterogeneity_results.csv`  
  Cochran’s Q heterogeneity results.
- `data/derived/all_pleiotropy_results.csv`  
  MR-Egger intercept test results.
- `data/derived/all_steiger_results.csv`  
  Steiger directionality summary (computed when possible; otherwise flagged as not available with reasons).

### Step 3. Filtering MR-prioritized genes (47 for colon cancer; 39 for psoriasis)
From the MR results, we applied predefined criteria to prioritize genes with putative causal effects for each outcome.

**Filtering criteria (summary)**
- Primary IVW association: nominal significance (e.g., IVW P < 0.05)
- Consistent effect direction across methods (risk OR > 1 or protective OR < 1)
- No strong evidence of heterogeneity (Q_pval > 0.05)
- No evidence of directional pleiotropy (MR-Egger intercept P > 0.05)
- Steiger directionality retained when Steiger was computable

**Key outputs**
- `data/derived/COLON_CANCER_MR_genes.csv` (or similar)  
  Final prioritized gene list for colon cancer (**n = 47**).
- `data/derived/PSORIASIS_MR_genes.csv` (or similar)  
  Final prioritized gene list for psoriasis (**n = 39**).

> If your project uses different file names for the final gene lists, update this section accordingly.

### Step 4. Per-gene visualization and diagnostic plots
For each MR-prioritized gene (or a subset), we generated multiple diagnostic plots to assess robustness, including:
- Scatter plot
- Funnel plot
- Leave-one-out plot
- Single-SNP plot
- Radial MR plot (outlier screening)

**Script(s)**
- `code/05_mr_diagnostics_plots.R`

**Key outputs**
- `outputs/diagnostics/`  
  Per-gene plots for colon cancer and psoriasis MR results.

### Step 5. Intersection analysis between colon and psoriasis MR-prioritized genes
We computed the overlap between the colon and psoriasis MR-prioritized gene sets and visualized the intersection (e.g., Venn diagram or UpSet plot).

**Script(s)**
- `code/06_intersection_plots.R`

**Key outputs**
- `outputs/intersection/venn_or_upset_plot.*`
- `data/derived/overlap_genes.csv` (if saved)

## How to reproduce (typical order)

1. Generate cross-disease candidate genes (415 genes) from scRNA-seq overlap  
   - Run: `code/01_scRNA_processing_and_overlap.R` (or your equivalent script)
2. Run MR for colon cancer outcome  
   - Run: `code/03_MR_pipeline_colon.R`
3. Run MR for psoriasis outcome  
   - Run: `code/04_MR_pipeline_psoriasis.R`
4. Apply filtering rules and export prioritized gene lists (47 and 39)  
   - Run: `code/07_filter_mr_results.R` (or the filtering section within your MR scripts)
5. Generate diagnostic plots (scatter, funnel, leave-one-out, radial MR, etc.)  
   - Run: `code/05_mr_diagnostics_plots.R`
6. Visualize overlap between two prioritized gene sets  
   - Run: `code/06_intersection_plots.R`

## Software and environment

Analyses were conducted in R. Key packages include:
- `TwoSampleMR`, `RadialMR`, `VariantAnnotation`, `gwasglue`, `dplyr`, `ggplot2`

For reproducibility, please refer to:
- `requirements/sessionInfo.txt`

## Contact

For questions about code or reproducibility, please open an issue in this repository.
