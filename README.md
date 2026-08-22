# NSCLC Integrated Transcriptomic Resource

This repository provides the major computational workflows used to construct an integrated transcriptomic resource of non-small-cell lung cancer (NSCLC), including scRNA-seq preprocessing, multi-cohort integration, integration benchmarking, and spatial deconvolution.

## Repository structure

```text
LungCancer_integrated/
├── scripts/
│   ├── runpipeline.py
│   ├── bbknn_integration.py
│   ├── combat_integration.py
│   ├── liger_integration.py
│   ├── scanorama_integration.py
│   └── scvi_scanvi_integration.py
├── tutorial/
│   ├── scRNA_plot.ipynb
│   ├── cell2location.ipynb
│   └── ST_plot.ipynb
└── LICENSE
```

## 1. scRNA-seq quality control and preprocessing

The main preprocessing workflow is implemented in `scripts/runpipeline.py` using Scanpy v1.9.6 and Scrublet v0.2.3.

Parameters used in this study:

| Parameter | Value |
|---|---:|
| Minimum sample size | 1,000 cells |
| Maximum sample size | 40,000 cells |
| Minimum genes per cell | 200 |
| Maximum genes per cell | 8,000 |
| Maximum mitochondrial percentage | 25% |
| Minimum cells per gene | 3 |
| Normalization target | 10,000 counts/cell |
| PCA dimensions | 30 |
| Maximum Harmony iterations | 20 |
| Leiden resolution | 1.5 |

Potential doublets were identified using Scrublet with automatically determined thresholds. After quality control, 547,360 cells were retained.

Example:

```bash
python scripts/runpipeline.py \
    --path /path/to/input_h5ad/ \
    --min_sample_size 1000 \
    --max_sample_size 40000 \
    --min_genes 200 \
    --max_genes 8000 \
    --min_cells 3 \
    --mt_pct 25 \
    --resolution 1.5 \
    --max_harmony_iter 20
```

## 2. scRNA-seq integration and benchmarking

Harmony was used as the primary integration method. Harmony-corrected PCA embeddings were used for neighborhood construction, UMAP visualization, and Leiden clustering.

Harmony was benchmarked against BBKNN, ComBat, LIGER, Scanorama, scVI, and scANVI. The corresponding scripts are provided in `scripts/`.

Integration performance was quantitatively assessed using **scIB** metrics for both batch correction and biological conservation.

Example benchmark script:

```text
scripts/bbknn_integration.py
```

The BBKNN workflow uses `orig.ident` as the batch key, 3,000 highly variable genes, 30 principal components.

## 3. Spatial deconvolution using cell2location

Spatial deconvolution was performed using cell2location v0.1.4 with the integrated NSCLC scRNA-seq atlas as reference.

The workflow is provided in:

```text
tutorial/cell2location.ipynb
```

### Reference model

The `RegressionModel` was trained using:

| Parameter | Value |
|---|---:|
| Reference dataset | `NSCLC_ALL.h5ad` |
| Batch key | `orig.ident` |
| Cell-type label | `celltype_LV2` |
| Highly variable genes | 3,000 |
| Training epochs | 250 |
| Training batch size | 4,096 |
| Posterior samples | 1,000 |
| Posterior batch size | 2,500 |
| Accelerator | GPU |

Cell-type-specific signatures were derived from `means_per_cluster_mu_fg`.

### Spatial model

The Cell2location model was trained using:

| Parameter | Value |
|---|---:|
| `N_cells_per_location` | 5 |
| `detection_alpha` | 200 |
| Training epochs | 1,000 |
| `train_size` | 1 |
| Posterior samples | 1,000 |
| Accelerator | GPU |

The `q05_cell_abundance_w_sf` estimates were used for downstream spatial visualization and spatial-region analyses.

The analysis was additionally repeated using the complete integrated scRNA-seq atlas as the reference.
