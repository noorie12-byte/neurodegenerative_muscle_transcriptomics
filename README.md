# Integrative Transcriptomics of Skeletal Muscle in Neurodegenerative Disease

**Does Skeletal Muscle Reflect Neurodegenerative Pathology Independently of the CNS?**

MSc Bioinformatics Thesis — University of Birmingham, 2025/26  
**Author:** Noor-Ul-Haya Ahmed  
**Supervisor:** Dr. Ekene Anakor  

---

## Overview

This repository contains the full analysis pipeline for an integrative transcriptomic study of skeletal muscle across three neurodegenerative diseases: amyotrophic lateral sclerosis (ALS), Parkinson's disease (PD) and VCP/IBMPFD myopathy.

The study tests whether skeletal muscle harbours a transcriptional signature in neurodegeneration that is independent of CNS pathology — and whether that signature could support peripheral biomarker discovery.

---

## Key Findings

- **43 genes** shared between ALS and PD skeletal muscle at FDR < 0.05, enriched for muscle-contractile biology
- **Direction discordance** — sarcomeric genes downregulated in ALS but upregulated in PD, reflecting denervation atrophy vs compensatory remodelling
- **MGST3 and RNASEL** identified as the strongest neurodegeneration-enriched candidates (Tier 3 — pan-ND, partially specific)
- **8 ALS-specific genes** — FOSL1, TPRG1, PILRB, ACSS1, RFLNA, SLC40A1, B3GLCT, PDE7A — not significant in any non-neurodegenerative muscle disease comparator
- **59 of 60 gene-tissue comparisons** classified as muscle-specific — top candidates not significantly dysregulated in ALS spinal cord, ALS motor cortex or PD substantia nigra across two independent CNS datasets

---

## Datasets

| GSE ID | Disease | Tissue | Role |
|---|---|---|---|
| GSE3307 | ALS | Skeletal muscle | Primary ALS |
| GSE26276 | ALS | Skeletal muscle | Directional only |
| GSE122261 | ALS | Myotubes (cell line) | Secondary |
| GSE100188 | ALS | Skeletal muscle | Small RNA track |
| GSE140089 | PD | Skeletal muscle | Primary PD |
| GSE128177 | PD | Skeletal muscle | Shared cohort |
| GSE30806 | VCP-myopathy | Skeletal muscle | Primary VCP |
| GSE12685 | AD | Brain (PFC) | Brain comparator |
| GSE28146 | AD | Brain (hippocampus) | Brain comparator |
| GSE26927 | ALS/PD/AD | CNS (multi-region) | CNS comparison |
| GSE122649 | ALS | Motor cortex | CNS comparison |

All datasets are publicly available from the [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/).

---

## Repository Structure

```
nd_project/
│
├── 1_download_datasets.Rmd          # Download raw data from GEO
├── 2_QC_Preprocessing_Microarray.Rmd # Microarray QC and normalisation
├── 3_QC_Preprocessing_RNAseq.Rmd    # RNA-seq QC and preprocessing
├── 4_Differential_Expression.Rmd    # DE analysis (limma + DESeq2)
├── 5_Cross_Study_Analysis.Rmd       # Consistency, overlap, pathways, biomarkers
├── 6_visualisations.Rmd             # Volcano plots
├── 7_GSEA.Rmd                       # Gene set enrichment analysis
├── 8_hallmark_pathways.Rmd          # Hallmark MSigDB enrichment
├── 9_brain_vs_muscle_comparison.Rmd # CNS vs muscle comparison
├── gse3307_specificity_test.Rmd     # Internal specificity analysis
│
├── metadata/                        # Sample metadata files
├── raw_data/                        # Raw expression matrices (not tracked)
├── Normalised_Data/                 # Normalised expression matrices
├── DE_Results/                      # Differential expression results
├── Consistency/                     # Cross-study consistency results
├── Pathway/                         # GO, KEGG enrichment results
├── GSEA/                            # GSEA results and plots
├── Hallmark/                        # Hallmark enrichment results
├── Biomarkers/                      # Biomarker scoring and tier tables
├── Specificity/                     # Internal specificity results
├── CNS_Data/                        # CNS dataset results
├── Volcano_Plots/                   # Volcano plot figures
└── figures/                         # Thesis figures
```

---

## Analysis Pipeline

```
Raw GEO data
     │
     ▼
QC + Normalisation
(log2 transform, quantile normalisation, PCA, boxplots)
     │
     ▼
Differential Expression
(limma for microarray, DESeq2 for RNA-seq)
     │
     ├──► Cross-Study Consistency (363 robust ALS genes, 2 PD)
     │
     ├──► Cross-Disease Comparison (43 shared ALS/PD genes)
     │         │
     │         ▼
     │    Pathway Enrichment (GO, KEGG, Hallmark, GSEA)
     │
     ├──► Internal Specificity Testing
     │    (ALS vs JDM, SPG, DMD, FSHD, AQM within GSE3307)
     │
     ├──► CNS vs Muscle Comparison
     │    (GSE26927 + GSE122649)
     │
     └──► Biomarker Prioritisation
          (Weighted scoring + Tiered framework)
```

---

## Software and Dependencies

All analyses performed in **R version 4.5.2**

| Package | Version | Use |
|---|---|---|
| GEOquery | 2.78.0 | Dataset download |
| limma | 3.66.0 | Microarray DE |
| DESeq2 | 1.50.2 | RNA-seq DE |
| clusterProfiler | 4.18.4 | Pathway enrichment |
| org.Hs.eg.db | 3.22.0 | Gene ID mapping |
| AnnotationDbi | 1.72.0 | Annotation |
| msigdbr | 26.1.0 | Hallmark gene sets |
| UpSetR | 1.4.1 | Set visualisation |
| pheatmap | 1.0.13 | Heatmaps |
| ggplot2 | 4.0.3 | Visualisation |
| dplyr | 1.2.1 | Data manipulation |

Install all packages with:

```r
# CRAN packages
install.packages(c("ggplot2", "dplyr", "UpSetR", "pheatmap", "ggrepel"))

# Bioconductor packages
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery", "limma", "DESeq2", "clusterProfiler",
  "org.Hs.eg.db", "AnnotationDbi", "msigdbr",
  "hgu133a.db", "hgu133plus2.db",
  "hugene10sttranscriptcluster.db",
  "huex10sttranscriptcluster.db"
))
```

---

## How to Run

1. Clone this repository
2. Open `nd_project.Rproj` in RStudio
3. Run scripts in order (1 → 10)
4. Raw data will be downloaded automatically from GEO in script 1

> **Note:** Raw expression matrices and large intermediate files are not tracked in this repository due to file size. Running script 1 will download all required data from GEO.

---

## Biomarker Tier Framework

| Tier | Label | Description | n genes |
|---|---|---|---|
| 1 | Pan-ND and ALS-specific | Significant across 3 diseases AND not generic | 0 |
| 2 | ALS-specific only | ALS-specific, not in non-ND muscle diseases | 8 |
| 3 | Pan-ND, partially specific | Across 3 diseases, partially specific | 2 |
| 4 | Partially specific | Some specificity, not pan-ND | 10 |
| 5 | Pan-ND but generic | Across 3 diseases but generic | 5 |
| 6 | Generic muscle disease | Responds to any muscle pathology | 18 |

---

## Citation

If you use this code or data, please cite:

> Ahmed N-U-H (2026). Integrative Transcriptomics Analysis of Skeletal Muscle in Neurodegenerative Disease. MSc Bioinformatics Thesis, University of Birmingham.

---

## Contact

**Noor-Ul-Haya Ahmed**  
MSc Bioinformatics, University of Birmingham  
Supervised by Dr. Ekene Anakor

Analysis scripts available upon request.
