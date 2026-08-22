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

# ==================== 配置 ====================
DATA_PATH = "NSCLC_DataSet/NSCLC_ALL.h5ad"
OUT_DIR   = "integration_results"

BATCH_KEY = "orig.ident"

N_HVG       = 3000             
N_PCS       = 30               
LEIDEN_RES  = 0.5
SEED        = 42

os.makedirs(OUT_DIR, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(OUT_DIR, "bbknn_run.log")),
    ],
)
log = logging.getLogger("bbknn")
sc.settings.n_jobs = 40
sc.settings.seed = SEED
np.random.seed(SEED)


def main():
    t0 = time.time()
    log.info(f"load {DATA_PATH} ...")
    adata = sc.read_h5ad(DATA_PATH)
    log.info(f"load completed: {adata.shape[0]} cells x {adata.shape[1]} genes, {time.time()-t0:.0f}s")


    if adata.raw is not None and list(adata.raw.var_names) == list(adata.var_names):
        adata.X = adata.raw.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG)
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=N_PCS, svd_solver="arpack")
    log.info(f"preprocess: {adata.shape[0]} cells x {adata.shape[1]} HVG, {time.time()-t0:.0f}s")

    try:
        sc.external.pp.bbknn(adata, batch_key=BATCH_KEY)
    except AttributeError:
        import bbknn
        bbknn.bbknn(adata, batch_key=BATCH_KEY, use_rep="X_pca", n_pcs=N_PCS, trim=None)
    log.info(f"BBKNN completed, {time.time()-t0:.0f}s")


    sc.tl.umap(adata)
    adata.obsm["X_umap_bbknn"] = adata.obsm["X_umap"].copy()

    out = os.path.join(OUT_DIR, "bbknn_integrated.h5ad")
    adata.write_h5ad(out)
    log.info(f"completed! {time.time()-t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
