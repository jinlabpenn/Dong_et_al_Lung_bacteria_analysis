# Lung Microbiome Analysis

## Project Overview

This repository contains the R code used to analyze lung microbiome sequencing data and generate the figures presented in our manuscript. The pipeline includes taxonomy assignment using multiple reference databases, decontamination analysis, and visualization of microbiome composition in different sample groups.

## Author
QIANG DONG (Original code)

## Last Updated
March 31, 2025

## Software Requirements

This analysis was performed using R under development (unstable) (2023-12-29 r85751 ucrt) on Windows 11 x64 (build 22631).

### R Packages

The following R packages (with versions) were used:

#### Base Packages
- grid, stats4, stats, graphics, grDevices, utils, datasets, methods, base

#### Analysis-specific Packages
- DECIPHER (3.2.0)
- ggnewscale (0.5.0)
- ggtreeExtra (1.16.0)
- ggtree (3.14.0)
- DESeq2 (1.46.0)
- SummarizedExperiment (1.36.0)
- Biobase (2.66.0)
- MatrixGenerics (1.18.0)
- matrixStats (1.4.1)
- GenomicRanges (1.58.0)
- tidyr (1.3.1)
- pheatmap (1.0.12)
- dada2 (1.34.0)
- Rcpp (1.0.13-1)
- ComplexHeatmap (2.22.0)
- decontam (1.26.0)
- ggrepel (0.9.6)
- UpSetR (1.4.0)
- here (1.0.1)
- plotly (4.10.4)
- circlize (0.4.16)
- openxlsx (4.2.7.1)
- cowplot (1.1.3)
- knitr (1.49)
- dplyr (1.1.4)
- RColorBrewer (1.1-3)
- ape (5.8-1)
- gridExtra (2.3)
- stringr (1.5.1)
- vegan (2.6-8)
- lattice (0.22-6)
- permute (0.9-7)
- ggplot2 (3.5.1)
- reshape2 (1.4.4)
- phyloseq (1.50.0)
- Biostrings (2.74.1)
- GenomeInfoDb (1.42.1)
- XVector (0.46.0)
- IRanges (2.40.1)
- S4Vectors (0.44.0)
- BiocGenerics (0.52.0)

## Directory Structure

```
.
├── data/
│   ├── raw/                  # Raw input data
│   └── processed/            # Processed phyloseq objects
├── results/
│   ├── Figures/              # Output figures (PDF, PNG, etc.)
│   └── Tables/               # Output data tables (CSV, Excel, etc.)
│       ├── contamination_analysis/   # Contamination analysis results
│       ├── extracted_asvs/           # Extracted ASV data
│       └── phylotree/                # Phylogenetic tree data
└── ref_databases/            # Reference databases for taxonomy
```

## Reference Databases
The code automatically downloads the following reference databases:

1. **Silva v138.2**
   - `silva_nr99_v138.2_toGenus_trainset.fa.gz` [Download Link](https://www.dropbox.com/scl/fi/lnwwdbqm5eyqgj8xktw3f/silva_nr99_v138.2_toGenus_trainset.fa.gz?rlkey=tjwi2bb4lkejemvie93dsxeu9&st=slcdc9h3&dl=0)
   - `silva_nr99_v138.2_toSpecies_trainset.fa.gz` [Download Link](https://www.dropbox.com/scl/fi/erhd1gcwuncxkkd5nxeku/silva_nr99_v138.2_toSpecies_trainset.fa.gz?rlkey=ilo06ftnjl9we6pfcgaypy5a1&st=bcywl890&dl=0)

2. **RDP v19**
   - `rdp_19_toGenus_trainset.fa.gz` [Download Link](https://www.dropbox.com/scl/fi/uvr2u3nwzxunt72cfnjyl/rdp_19_toGenus_trainset.fa.gz?rlkey=kdlmuuxcu3hdyhmcoe0l2rbfz&st=u8m032tg&dl=0)
   - `rdp_19_toSpecies_trainset.fa.gz` [Download Link](https://www.dropbox.com/scl/fi/7cfeowmxfw19ftouomzth/rdp_19_toSpecies_trainset.fa.gz?rlkey=sjex7itpw766df6gn82t23f3q&st=9n0sp2ye&dl=0)

3. **GreenGenes 2024_09**
   - `gg2_2024_09_toGenus_trainset.fa.gz` [Download Link](https://www.dropbox.com/scl/fi/qc82i4q8ocfjapjoytw6f/gg2_2024_09_toGenus_trainset.fa.gz?rlkey=0j0lhzf5wnarv8xq2lb1wkmyf&st=lghz38r0&dl=0)
   - `gg2_2024_09_toSpecies_trainset.fa.gz` [Download Link](https://www.dropbox.com/scl/fi/zcd24frb4idegpe3r0fuv/gg2_2024_09_toSpecies_trainset.fa.gz?rlkey=7fxzvqmoego7v24051zgouhbv&st=2svmowvv&dl=0)

## Manuscript Figures Reproduction

This code is specifically designed to reproduce all figures presented in our manuscript. The main figures are:

### Main Figures

- **Figure 2B**: Beta diversity analysis with PCoA plot
  - Shows sample clustering based on Bray-Curtis dissimilarity
  - Generated in section 'Figure 2B: Beta Diversity Analysis'

- **Figure 2C**: Microbiome composition visualization
  - Includes genus-level stacked bar plots and presence/absence visualization
  - Generated in section 'Figure 2C: Stacked Bar Plot and Presence/Absence Plot'

- **Figure 2D**: Phylogenetic tree visualization
  - Displays taxa with prevalence and abundance data
  - Generated in section 'Figure 2D: Phylogenetic plot'

### Supplementary Figures

- **Figure S2B**: Phylum stacked bar plot
  - Shows phylum-level composition across samples
  - Generated in section 'Figure S2B: Phylum stacked bar plot'

- **Figure S2C**: Genus-level differential abundance
  - Volcano plot showing differentially abundant genera
  - Generated in section 'Figure S2C: DESeq2 genus level differential abundance'

- **Figure S2D**: ASV-level differential abundance
  - Volcano plot showing differentially abundant ASVs
  - Generated in section 'Figure S2D: DESeq2 ASV level differential abundance'

## Usage

### Data Requirements

The pipeline requires:
1. A DADA2-processed phyloseq object (`dada2_phyloseq.rds`)
2. Corresponding metadata file (`metadata.csv`)
3. Reference databases for taxonomy assignment (automatically downloaded)

### Execution

To reproduce the figures in the manuscript:

1. **Environment Setup**:
   ```R
   source("main_script.R")
   ```
   This will:
   - Install and load all required packages
   - Create the project directory structure
   - Download reference databases from Dropbox

2. **Raw Data Requirements**:
   - `data/raw/dada2_phyloseq.rds`: DADA2 processed phyloseq object
   - `data/raw/metadata.csv`: Sample metadata with samples as rows

## Key Parameters

The pipeline uses several important parameters that can be customized:

- `freq_ratio_cutoff` = 0.297: Frequency ratio threshold for contamination
- `sample_freq_cutoff` = 0.002: Sample frequency threshold for detection
- `alpha` (DESeq2) = 0.05: Significance threshold for differential abundance

## Outputs

### Processed Data Files
- `DADA2_with_new_metadata.rds`: Initial phyloseq with metadata
- `reference_db_configs.rds`: Configuration details for reference databases
- `silva_taxonomy_raw.rds`, `rdp_taxonomy_raw.rds`, `gg_taxonomy_raw.rds`: Raw taxonomy assignments
- `combined_taxonomy.rds`: Combined taxonomy from all databases
- `DADA2_with_combined_taxonomy.rds`: Final phyloseq with combined taxonomy

### Tables

- `Table1A/B/C_*_taxonomy.csv/rds`: Taxonomy assignments from different databases
- `Table1D_TAXannotate_full.csv`: Consolidated taxonomy with database recommendations
- `Table2A_per_sample_contamination_analysis.xlsx`: Detailed contamination metrics
- `Table2B_extracted_asvs_*.xlsx/csv`: Extracted ASVs passing contamination filters
- `Table3A-F_*.csv`: Summaries, beta diversity stats, and compositional data

### Figures

- `Figure2B_beta_diversity.pdf/jpg`: PCoA plot showing sample clustering
- `Figure2C_microbiome_composition.pdf`: Genus-level composition visualization
- `Figure2D_phylogenetic_tree.pdf/png`: Phylogenetic tree with abundance data
- `FigureS2B_phylum_*.pdf`: Phylum-level composition visualizations
- `FigureS2C_genus_volcano.pdf`: Differential abundance at genus level
- `FigureS2D_asv_volcano.pdf`: Differential abundance at ASV level

## Reference Database Structure

The taxonomy assignment integrates results from three different reference databases:
- **Silva**: Standard taxonomy format
- **RDP**: Standard taxonomy format
- **GreenGenes**: Uses prefixes (d__, p__, c__, etc.) which are removed during processing

## Citation

If you use this code in your research, please cite our manuscript:

```
[Citation to be added after publication]
```

## Contact

For questions or issues related to this code, please contact: jinlabpenn@gmail.com
