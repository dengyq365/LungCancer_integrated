#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import sys
import time
import warnings
import logging

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
import scanpy as sc
import anndata as ad

# ==================== Config ====================
DATA_PATH = "./NSCLC_DataSet/NSCLC_ALL.h5ad"
OUT_DIR   = "./integration_results"

BATCH_KEY = "orig.ident"
LABEL_KEY = "celltype_LV1"        

N_HVG       = 3000                # number of HVGs
N_LATENT    = 30                  # latent dimension
N_LAYERS    = 2                   # encoder layers (1 is faster)
MAX_EPOCHS  = 250                 # ~15 min per model on GPU; reduce to 100 for speed
BATCH_SIZE  = 2048
N_NEIGHBORS = 30
LEIDEN_RES  = 0.5

# ==================== Setup ====================
os.makedirs(OUT_DIR, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(OUT_DIR, "scvi_scanvi_run.log")),
    ],
)
log = logging.getLogger("scvi_scanvi")



def get_device():
    import torch
    return "gpu" if torch.cuda.is_available() else "cpu"


def train_scvi(adata, device):
    import scvi
    scvi.settings.num_threads = 40
    scvi.model.SCVI.setup_anndata(adata, batch_key=BATCH_KEY, layer="counts")
    model = scvi.model.SCVI(adata, n_latent=N_LATENT, n_layers=N_LAYERS, gene_likelihood="nb")
    model.train(
        max_epochs=MAX_EPOCHS, batch_size=BATCH_SIZE,
        early_stopping=True, early_stopping_patience=20,
        accelerator=device, devices=1,
    )
    adata.obsm["X_scvi"] = model.get_latent_representation()
    return model


def train_scanvi(adata, device):
    import scvi
    # NaN labels -> "Unknown" (if fully labeled, "Unknown" is just an unused placeholder)
    if adata.obs[LABEL_KEY].isna().any():
        adata.obs[LABEL_KEY] = adata.obs[LABEL_KEY].astype("object").fillna("Unknown").astype("category")
    unlabeled = "Unknown"

    scvi.model.SCANVI.setup_anndata(
        adata, batch_key=BATCH_KEY, labels_key=LABEL_KEY,
        unlabeled_category=unlabeled, layer="counts",
    )
    model = scvi.model.SCANVI(adata, n_latent=N_LATENT, n_layers=N_LAYERS, gene_likelihood="nb")
    model.train(
        max_epochs=MAX_EPOCHS, batch_size=BATCH_SIZE,
        early_stopping=True, early_stopping_patience=20,
        accelerator=device, devices=1,
    )
    adata.obsm["X_scanvi"] = model.get_latent_representation()
    return model


def downstream(adata, name, use_rep):
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=N_NEIGHBORS)
    sc.tl.umap(adata)
    adata.obsm[f"X_umap_{name}"] = adata.obsm["X_umap"].copy()
    sc.tl.leiden(adata, resolution=LEIDEN_RES, key_added=f"leiden_{name}")


def main():
    t0 = time.time()
    device = get_device()
    log.info(f"Device: {device}")

    # ---- 1. Load ----
    log.info(f"Reading {DATA_PATH} ...")
    adata = sc.read_h5ad(DATA_PATH)
    log.info(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes, {time.time()-t0:.0f}s")

    # ---- 2. Raw counts -> layer ----
    if adata.raw is not None and list(adata.raw.var_names) == list(adata.var_names):
        adata.layers["counts"] = adata.raw.X.copy()
        log.info("counts extracted from adata.raw")
    else:
        adata.layers["counts"] = adata.X.copy()
        log.warning("raw missing / gene order mismatch, using X as counts; make sure it is integer!")
    if not np.issubdtype(adata.layers["counts"].dtype, np.integer):
        log.warning("counts are not integer -> scVI/scANVI may fail or give unreliable results!")
    adata.raw = None

    # ---- 3. Highly variable genes ----
    adata.X = adata.layers["counts"]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG)
    adata = adata[:, adata.var.highly_variable].copy()
    log.info(f"After preprocessing: {adata.shape[0]} cells x {adata.shape[1]} HVGs, {time.time()-t0:.0f}s")

    # ---- 4. scVI ----
    log.info("===== Training scVI =====")
    scvi_model = train_scvi(adata, device)
    scvi_model.save(os.path.join(OUT_DIR, "scvi_model"), overwrite=True)
    log.info(f"scVI done, {time.time()-t0:.0f}s")

    # ---- 5. scANVI ----
    log.info("===== Training scANVI =====")
    scanvi_model = train_scanvi(adata, device)
    scanvi_model.save(os.path.join(OUT_DIR, "scanvi_model"), overwrite=True)
    log.info(f"scANVI done, {time.time()-t0:.0f}s")

    # ---- 6. Downstream: neighbors -> UMAP -> leiden ----
    downstream(adata, "scvi", "X_scvi")
    downstream(adata, "scanvi", "X_scanvi")

    # ---- 7. Save ----
    out = os.path.join(OUT_DIR, "scvi_scanvi_integrated.h5ad")
    adata.write_h5ad(out)
    log.info(f"Done! total {time.time()-t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
