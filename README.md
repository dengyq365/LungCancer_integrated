# LungCancer_integrated

# NSCLC_multiomics_reference_atlas

A harmonized multi-omics transcriptomic resource for non-small-cell lung cancer (NSCLC), integrating bulk RNA sequencing, single-cell RNA sequencing (scRNA-seq), and spatial transcriptomics datasets from publicly available cohorts.

This repository provides standardized computational workflows, processed datasets, and annotation resources for studying tumor microenvironment (TME) heterogeneity, immune–stromal interactions, and immunotherapy-associated transcriptional programs in NSCLC.

---

# Overview

This project integrates:

* Bulk transcriptomic datasets from TCGA and GEO cohorts
* Single-cell RNA sequencing datasets from multiple NSCLC studies
* Spatial transcriptomic datasets for spatial validation
* Clinical annotations including:

  * tumor versus adjacent tissue
  * histological subtype
  * survival information
  * immunotherapy response

The integrated dataset contains:

* 547,360 QC-passed single cells
* 217 tumor and adjacent tissue samples
* harmonized cell-type annotations
* refined myeloid and stromal subpopulation labels

The processed dataset may serve as:

* a reference atlas for NSCLC TME analysis
* a reference for spatial deconvolution (e.g., cell2location)
* a resource for biomarker discovery
* a benchmark for computational integration methods

---

# Repository Structure

```bash
.
├── data/
│   ├── bulk/
│   ├── scRNA/
│   ├── spatial/
│   └── metadata/
│
├── scripts/
│   ├── preprocessing/
│   ├── integration/
│   ├── annotation/
│   ├── cell2location/
│   └── visualization/
│
├── results/
│   ├── figures/
│   ├── DEG/
│   ├── annotation/
│   └── spatial/
│
├── README.md
└── requirements.txt
```

---

# Included Public Datasets

## Bulk transcriptomic datasets

* TCGA-LUAD
* TCGA-LUSC
* GSE19804
* GSE151101
* GSE32863
* GSE31210
* GSE44077
* GSE207422

## Single-cell transcriptomic datasets

* GSE131907
* GSE148071
* GSE127465
* GSE171145
* GSE149655
* GSE117570
* GSE201333
* GSE207422
* GSE179994
* GSE123904
* E-MTAB-6653
* E-MTAB-6149
* Lung Cancer Explorer

## Spatial transcriptomic datasets

* E-MTAB-13530

---

# Computational Workflow

The single-cell integration workflow includes:

1. Data loading and assembly
2. Doublet removal using Scrublet
3. Quality-control filtering
4. Normalization and highly variable gene selection
5. Batch correction using Harmony
6. UMAP dimensionality reduction
7. Leiden clustering
8. Differential expression analysis
9. Cell-type annotation

The main integration pipeline is implemented in:

```bash
scripts/integration/runpipeline.py
```

Key functions include:

* `assemble_data()`
* `run_scrublet()`
* `before_integrated()`
* `run_harmony()`
* `run_deg()`

---

# Quality Control

Cells were filtered according to the following criteria:

* > 200 detected genes
* <25% mitochondrial gene expression
* <8,000 detected genes

Potential doublets were identified using Scrublet and removed before downstream analyses.

Batch effects across cohorts and sequencing platforms were corrected using Harmony integration.

---


# Spatial Transcriptomics

Spatial transcriptomic datasets were used as an orthogonal validation modality to provide spatial context for cell-type localization within NSCLC tissues.

The integrated scRNA-seq atlas may additionally support downstream spatial deconvolution approaches including:

* cell2location
* Tangram
* SpaOTsc

---

# Clinical Annotations

Where available, the following clinical metadata were harmonized:

* age
* sex
* AJCC stage
* overall survival
* progression-free survival
* immunotherapy response

Immunotherapy datasets include response annotations:

* CR
* PR
* SD
* PD

---



# Citation

If you use this resource, please cite:

Deng YQ et al.
*A clinically annotated multi-omics transcriptomic resource of the NSCLC tumor microenvironment.*
(Under review)

---

# Data Availability

Processed datasets and metadata are publicly available through:

* Zenodo: [DOI_LINK]
* GEO: [ACCESSION]
* Original public repositories as referenced in the manuscript

---


