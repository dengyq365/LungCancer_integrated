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


```bash
scripts/integration/runpipeline.py --min_sample_size=1000 --path=h5ad 
```


---

# Quality Control

Cells were filtered according to the following criteria:

* >200 detected genes
* <25% mitochondrial gene expression
* <8,000 detected genes
* >1000 cells per sample

Potential doublets were identified using Scrublet and removed before downstream analyses.

Batch effects across cohorts and sequencing platforms were corrected using Harmony integration.

---


# Spatial Transcriptomics

Spatial transcriptomic datasets were used as an orthogonal validation modality to provide spatial context for cell-type localization within NSCLC tissues.

The integrated scRNA-seq atlas may additionally support downstream spatial deconvolution approaches including:

* cell2location


---

# Clinical Annotations

Where available, the following clinical metadata were harmonized:

* age
* sex
* AJCC stage
* overall survival
* progression-free survival
* immunotherapy response


---




