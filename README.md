# RISE-Meta
Repository for the development and implementation of the **Rank-Based Identification of High-Dimensional Surrogate Markers Meta-Analysis (RISE-Meta)** method.
The methodology is described in the paper *Meta-Analysis of High-Dimensional Surrogate Markers*, available as a preprint from [https://arxiv.org/abs/2605.03819](https://arxiv.org/abs/2605.03819).

---

## Quick start

1. Clone this repository and open the project `RISE-meta.Rproj`. Ensure your working directory is the project root. 
2. Create the required project structure:
```r
source("analysis/create_project_structure.R")
```
3. Download the raw data (described below) and place it in `data-raw/`.
4. Install dependencies and recreate software environment
```r
if (!requireNamespace("renv", quietly = TRUE)){
  install.packages("renv")
}
if (!requireNamespace("BiocManager", quietly = TRUE)){  
  install.packages("BiocManager") # Required to restore bioconductor packages
}
renv::restore()
```
5. Run the full analysis pipeline:
```r
source("analysis/analysis_master.R")
```
Note that due to extensive simulations, this may take a long time (see [Runtime information](#runtime-information)). 

---

## Reproducing the software environment

This repository uses [`renv`](https://rstudio.github.io/renv/) for package management. To install this package andrestore the exact package versions used in this analysis:

```r
if (!requireNamespace("renv", quietly = TRUE)){
  install.packages("renv")
}
if (!requireNamespace("BiocManager", quietly = TRUE)){  
  install.packages("BiocManager") # Required to restore bioconductor packages
}
renv::restore()
```

The complete package environment is recorded in `renv.lock`.

---

## Data access

No raw or processed data is stored in this repository.

To run the project:

1. Download required raw data (see below)
2. Place all files in `data-raw/`
3. Run scripts in `analysis/`

---

### Low-dimensional application

Data is accessed via the **Surrogate** R package which is installed when using `renv::restore()`.

---

### High-dimensional application

#### ImmuneSpace data

1. Go to https://immunespace.org
2. Click **Resources → Immune Signature Projects**
3. Open **Project 2** (Immune Signatures 2) 
4. Click **Proceed to study data** (account required)
5. Open **Clinical and Assay Data**

Download:

- `all_noNorm_eset.rds`

From the **Assay** section download:
- Neutralizing antibody titer
- Hemagglutination inhibition (HAI)

Export each as:

1. Click dataset
2. Click **Export**
3. Select **Excel Workbook (.xlsx)**
4. Download file

Expected filenames:

- `neut_ab_titer_tag.xlsx`
- `hai_tag.xlsx`

where `tag` is replaced by a unique tag generated at the time of download. 
This is dealt with automatically by the scripts when loading in the data. 

---

### Gene set data

- **Blood Transcriptional Modules (BTMs)**
  - https://github.com/shuzhao-li/BTM/tree/master/BTM/datasets
  - File: `BTM_for_GSEA_20131008.gmt`

- **BloodGen3 Modules (BG3M)**
  - DOI: 10.1093/bioinformatics/btab121
  - File: `Suppl_File_1_BIOINF.xls`

---

## Project structure

- `data-raw/`
  Raw data files (not stored in this repository)

- `data/`
  Processed data (not stored in this repository)

- `R/`
  Helper functions

- `analysis/`
  Analysis scripts:
  - preprocessing
  - descriptive
  - application
  - simulation

- `output/`
  Outputs (not stored in this repository):
  - `results/`
  - `figures/`

- `manuscript/`
  Manuscript-related files:
  - `diagrams/` manually created figures

---

## Analysis

### Full pipeline

```r
source("analysis/analysis_master.R")
```

Runs:

- package installation and software environment reproduction
- project folder structure
- preprocessing
- descriptive analysis
- real-data application
- simulation study 

**NOTE:** due to the extensive number of simulations performed, running the simulation files may take a long time.
See [Runtime information](#runtime-information) for details and suggestions. 

---

### Preprocessing

```r
source("analysis/preprocessing/master_preprocessing.R")
```

Steps:

- clinical data preprocessing
- gene expression preprocessing
- immune response preprocessing
- gene set preprocessing
- dataset merging

Outputs of preprocessing go in `data/`

---

### Descriptive analysis

Runs descriptive analyses of the data and produces two supplementary figures.

```r
source("analysis/descriptive/descriptive_master.R")
```

Outputs:

- visual description of dataset (Web Figure 1)
- data for flowchart of samples after each preprocessing step (Web Figure 2)

Saved in:

- `output/results/descriptive/`
- `output/figures/descriptive/`

---

### Application

Runs main and supplementary analyses on real data. 

```r
source("analysis/application/application_master.R")
```
Outputs: 

- Low-dimensional data application (Figure 3)
- High-dimensional data application (Figures 3, 4)
- Ebola vaccine studies application (ebovac2 and prevac; hamburg excluded)
- Supplementary analyses (Supplementary Web Figures 9-16)

Saved in:

- `output/figures/application/`

---

### Simulation

Runs scripts to produce simulation results. 

```r
source("analysis/simulation/simulation_master.R")
```

Outputs: 

- Main simulation results (Figures 1, 2) 
- Supplementary simulation results (Web Figures 4-8)

Saved in:

- `output/results/simulation/` (raw results data)
- `output/figures/simulation/` (figures)

--- 

## Runtime information

The full analysis pipeline took roughly 1 day on a standard laptop with the following specs.

- Dell Inc. Precision 5470 
- 12th Gen Intel® Core™ i7-12700H processor × 20 cores
- 16GB RAM

However, the grand majority of this computational time is taken by the simulation study. 
When excluding the simulation study from the `analysis_master.R` script, the runtime reduces to around 5 minutes. 

To reduce this computational time, one may reduce the number of repetitions in the simulation studies. 
In the simulation files, these are controlled by the variable `J`, which is set to 100,000. To reduce this, you
may replace this with a smaller number. 
In RStudio, on default key bindings this can be achieved by 

- pressing `Ctrl+Shift+F` (find in files)
- typing exactly `J = 100000` and clicking "find"
- clicking "Replace" in the "Find in Files" window
- typing e.g. `J = 10000` and clicking "Replace All"
- Running the master script

Otherwise, one may exclude the simulation files producing results for the Supplementary Information to speed up the computational time. This involves

- opening the file `analysis/simulation/simulation_master.R` 
- commenting out with `#` the lines following `# Supporting information`
- running `analysis/analysis_master.R`

This pipeline can be run in a few hours on a standard PC. 

Note that computational time may be significantly longer for Windows users due to incompatibility of this OS with our parallel computing setup (see [Parallelisation](#parallelisation)).

---

## Parallelisation

Parallel computing has been used to speed up parts of the analysis.
By default, the number of cores is detected automatically using
`parallel::detectCores()` and half of the available cores are used.

### Adjusting the number of cores

To increase or decrease the number of cores used (e.g. if RStudio crashes
or becomes unresponsive), replace the line:

```r
n.cores = parallel::detectCores(all.tests = FALSE, logical = TRUE) / 2
```

with a fixed number, e.g.:

```r
n.cores = 1
```

using the Find in Files approach described in the Runtime section above.

### Operating system caveats

The parallelisation in this repository uses `pbmcapply`, which is based on
forking via `parallel::mclapply`. **This is only supported on Unix-based
systems (Linux and macOS).** Windows is not supported, and the code will fall back to sequential computation (i.e. `n.cores = 1`), 
so the runtime may be significantly slower for windows users. 

---

## Software and version information

Below is a summary of the output of 

```r
sessionInfo()
```

after loading all the libraries necessary for this repository. 

- **OS:** Ubuntu 24.04.4 LTS (x86_64-pc-linux-gnu)
- **R version:** 4.5.3 (2026-03-11)
- **BLAS/LAPACK:** OpenBLAS (LAPACK version 3.12.0)

### Attached packages

```
knitr_1.51          kableExtra_1.4.0    digest_0.6.39       viridis_0.6.5       viridisLite_0.4.3  
scales_1.4.0        janitor_2.2.1       Biobase_2.70.0      BiocGenerics_0.56.0 generics_0.1.4     
readxl_1.4.5        GSA_1.03.3          Surrogate_3.4.1     patchwork_1.3.2     boot_1.3-32        
lubridate_1.9.5     forcats_1.0.1       stringr_1.6.0       purrr_1.2.2         readr_2.2.0        
tidyr_1.3.2         tibble_3.3.1        ggplot2_4.0.3       tidyverse_2.0.0     SurrogateRank_3.0  
dplyr_1.2.1         renv_1.2.2          fs_2.1.0
```

### Packages loaded via namespace

```
tidyselect_1.2.1    farver_2.1.2        S7_0.2.2            fastmap_1.2.0       timechange_0.4.0   
lifecycle_1.0.5     survival_3.8-6      magrittr_2.0.5      compiler_4.5.3      rlang_1.2.0        
tools_4.5.3         xml2_1.5.2          RColorBrewer_1.1-3  withr_3.0.2         grid_4.5.3         
colorspace_2.1-2    MASS_7.3-65         cli_3.6.6           rmarkdown_2.31      miscTools_0.6-30   
reformulas_0.4.4    rstudioapi_0.18.0   tzdb_0.5.0          minqa_1.2.8         splines_4.5.3      
BiocManager_1.30.27 cellranger_1.1.0    vctrs_0.7.3         Matrix_1.7-5        sandwich_3.1-1     
hms_1.1.4           pbmcapply_1.5.1     systemfonts_1.3.2   glue_1.8.1          ggVennDiagram_1.5.7
nloptr_2.2.1        cowplot_1.2.0       stringi_1.8.7       gtable_0.3.6        lme4_2.0-1         
ComplexUpset_1.3.3  pillar_1.11.1       htmltools_0.5.9     R6_2.6.1            maxLik_1.5-2.2     
textshaping_1.0.5   Rdpack_2.6.6        evaluate_1.0.5      lattice_0.22-9      rbibutils_2.4.1    
snakecase_0.11.1    Rcpp_1.1.1-1.1      svglite_2.2.2       gridExtra_2.3       nlme_3.1-169       
xfun_0.57           zoo_1.8-15          pkgconfig_2.0.3 
```

Full version details are recorded in `renv.lock`. Run `renv::restore()` to reproduce the exact package environment.

---

## Citation

The methodology is described in: *Meta-Analysis of High-Dimensional Surrogate Markers* and implemented in the R package `SurrogateRank`, available on the CRAN.