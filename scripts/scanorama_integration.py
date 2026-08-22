#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys, time, warnings, logging
import numpy as np
warnings.filterwarnings("ignore")
import scanpy as sc

# ==================== Config ====================
DATA_PATH = "./NSCLC_DataSet/NSCLC_ALL.h5ad"
OUT_DIR   = "./integration_results"
BATCH_KEY = "orig.ident"
N_HVG = 3000
N_PCS = 50
N_NEIGHBORS = 30



os.makedirs(OUT_DIR, exist_ok=True)
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)s | %(message)s",
                    handlers=[logging.StreamHandler(sys.stdout),
                              logging.FileHandler(os.path.join(OUT_DIR, "scanorama_run.log"))])
log = logging.getLogger("scanorama")


def main():
    t0 = time.time()
    adata = sc.read_h5ad(DATA_PATH)
    log.info(f"Loaded {adata.shape[0]} x {adata.shape[1]}, {time.time()-t0:.0f}s")
    

    # Preprocess: normalize + log1p + HVG + scale + PCA (scanorama_integrate uses basis='X_pca')
    if adata.raw is not None and list(adata.raw.var_names) == list(adata.var_names):
        adata.X = adata.raw.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG)
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=N_PCS, svd_solver="arpack")

    # Scanorama (scanpy wrapper, auto split/integrate/write to obsm['X_scanorama'])
    log.warning("Scanorama runs all-pairs MNN over 204 batches, may be slow / memory heavy")
    sc.external.pp.scanorama_integrate(adata, key=BATCH_KEY, basis="X_pca")
    log.info(f"Scanorama done, {time.time()-t0:.0f}s")

    # Downstream
    sc.pp.neighbors(adata, use_rep="X_scanorama", n_neighbors=N_NEIGHBORS)
    sc.tl.umap(adata)
    adata.obsm["X_umap_scanorama"] = adata.obsm["X_umap"].copy()
    adata.write_h5ad(os.path.join(OUT_DIR, "scanorama_integrated.h5ad"))
    log.info(f"Done -> scanorama_integrated.h5ad, total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
