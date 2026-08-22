#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import sys
import time
import warnings
import logging
import anndata
import numpy as np

warnings.filterwarnings("ignore")
import scanpy as sc

# ==================== Config ====================
DATA_PATH = "/NSCLC_DataSet/NSCLC_ALL.h5ad"
OUT_DIR   = "/integration_results"

BATCH_KEY = "orig.ident"

N_HVG       = 3000              # number of HVGs selected by select_genes
LIGER_K     = 20                # number of iNMF factors
MINIBATCH   = 5000              # minibatch size
N_NEIGHBORS = 30
MIN_DIST    = 0.3

# ==================== Setup ====================
os.makedirs(OUT_DIR, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(OUT_DIR, "liger_run.log")),
    ],
)
log = logging.getLogger("liger")
sc.settings.n_jobs = 40



def main():
    import pyliger as liger

    t0 = time.time()

    # ---- 1. Load ----
    log.info(f"Reading {DATA_PATH} ...")
    adata = sc.read_h5ad(DATA_PATH)
    log.info(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes, {time.time()-t0:.0f}s")


    # ---- 2. Raw counts -> layer (LIGER requires raw counts) ----
    if adata.raw is not None and list(adata.raw.var_names) == list(adata.var_names):
        adata.layers["counts"] = adata.raw.X.copy()
        adata.raw = None
        log.info("raw counts stored in layers['counts']")
    else:
        adata.layers["counts"] = adata.X.copy()
        log.warning("raw missing / gene order mismatch, using X as counts!")

    # remove empty cells (avoids pyliger remove_missing shape bug)
    sc.pp.filter_cells(adata, min_counts=1)

    # ---- 3. Split into per-batch AnnData list ----
    cats = list(adata.obs[BATCH_KEY].cat.categories)
    adata_list = []
    for b in cats:
        sub = adata[adata.obs[BATCH_KEY] == b].copy()
        sub.X = sub.layers["counts"].copy()   # LIGER requires raw counts
        sub.uns["sample_name"] = b            # pyliger requires a dataset name
        adata_list.append(sub)
    log.info(f"split into {len(cats)} batches, {time.time()-t0:.0f}s")

    # ---- 4. LIGER pipeline (functions modify in place, do NOT reassign) ----
    liger_obj = liger.create_liger(adata_list, remove_missing=False)
    liger.normalize(liger_obj, remove_missing=False)
    liger.select_genes(liger_obj, n_genes=N_HVG)
    liger.scale_not_center(liger_obj)
    liger.optimize_ALS(liger_obj, k = LIGER_K)
    liger.quantile_norm(liger_obj)
    liger.run_umap(liger_obj, distance="cosine", n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST)
    adata = anndata.concat(liger_obj.adata_list)
    log.info(f"LIGER done, {time.time()-t0:.0f}s")

    # ---- 6. Save ----
    out = os.path.join(OUT_DIR, "liger_integrated.h5ad")
    adata.write_h5ad(out)
    log.info(f"Done! total {time.time()-t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
