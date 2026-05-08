#!/usr/bin/env python3

# Author: Dengyq
# Histroy: 20241121

#
# Describe: This is a simple pipeline for single cell RNA-seq data analysis
# Be careful when employed this script/inner-func to runpipeline
#

import pandas as pd
import numpy as np
import scanpy as sc
import scrublet as scr
import os, sys
from   optparse import OptionParser
import multiprocessing as mp
import anndata
from typing import Union


def load_data_from_path(path: str) -> object:
    hd5s = os.listdir(path)
    anna_list = [sc.read_h5ad(f'{path}/{x}') for x in hd5s]
    ## data tidiness
    anna_list = [sc.AnnData(X = x.X, obs = x.obs.iloc[:,:0], var = x.var.iloc[:,:0]) for x in anna_list]
    return(anna_list)



def run_scrublet(adata: sc.AnnData) -> sc.AnnData:
    adata.var_names_make_unique(join = "_")
    adata.obs_names_make_unique(join = "_")
    adata_new = sc.external.pp.scrublet(adata, verbose= False, copy=True)
    return(adata_new)



def assemble_data(anna_list: list, 
                  sample_name: list,
                  n_min: int = 1000, 
                  n_max: int = 30000,
                  join_type: str = "outer",
                  doDoublet: bool = False) -> sc.AnnData:
    '''
    merge and add batch categories
    [anna1, anna2, ...] -> anna
    '''
    
    sample_s = np.array([ np.logical_and(x.shape[0] > n_min,
                                         x.shape[0] < n_max) for x in anna_list])
    sample_name = [x for x,y in zip(sample_name, sample_s) if y]
    anna_list   = [x for x,y in zip(anna_list, sample_s) if y]
    # add batch var in obs and re_index
    for x,y in zip(anna_list, sample_name):
        x.obs['batch'] = y


    if doDoublet:
        if True: ## use multiple
            num_processes = mp.cpu_count() // 2 + 1
            pool          = mp.Pool(processes=num_processes)
            processed_data_list = pool.map(run_scrublet, anna_list)
            pool.close()
            pool.join()
            anna_list     = processed_data_list
        else: # 
            # # doublet detect
            for x in anna_list:
                sc.external.pp.scrublet(x, verbose= False)
    
    # anna_data = sc.AnnData.concatenate(*anna_list,
    #                                    batch_categories = sample_name,
    #                                    join = join_type)
    anna_data =  anndata.concat(anna_list, join = join_type,fill_value=0)
    anna_data.var_names_make_unique(join = "_")
    anna_data.obs_names_make_unique(join = "_")
    return(anna_data)



def addvar(anna_data: sc.AnnData,
        #    vars: list = ["mt", "ribo"], ## later
           specie: str = "human") -> sc.AnnData:
    '''
    addVar
    mit/RPS/RPL
    '''
    # mitochondrial genes
    if specie == "human":
        anna_data.var['mt'] = anna_data.var_names.str.startswith('MT-')
    elif specie == "mouse":
        anna_data.var['mt'] = anna_data.var_names.str.startswith('Mt-')
    # ribosomal genes
    anna_data.var['ribo'] = anna_data.var_names.str.startswith(("RPS","RPL"))
    sc.pp.calculate_qc_metrics(anna_data, qc_vars=['mt', 'ribo'],
                               percent_top=None, log1p=False, inplace=True)
    return(anna_data)


def stein_gate(anna_data: sc.AnnData, 
               min_genes: int = 200,
               max_genes: int = 8000,
               min_cells: int = 3,
               mt_pct:int = 25) -> sc.AnnData:
    anna_data = anna_data[anna_data.obs['pct_counts_mt'] < mt_pct]
    sc.pp.filter_cells(anna_data, min_genes=min_genes)
    sc.pp.filter_cells(anna_data, max_genes=max_genes)
    sc.pp.filter_genes(anna_data, min_cells=min_cells)


def before_integrated(anna_data: sc.AnnData, batch: Union[str, None] = None, hvg :bool = True):
    # save normalized counts in raw slot.
    anna_data.X = anna_data.X.astype(np.float32)
    anna_data.raw = anna_data
    # normalize to depth 10 000
    sc.pp.normalize_per_cell(anna_data, counts_per_cell_after=1e4)# logaritmize
    sc.pp.log1p(anna_data)
    sc.pp.highly_variable_genes(anna_data, batch_key = batch)
    sc.pp.scale(anna_data)
    sc.tl.pca(anna_data, svd_solver='arpack', use_highly_variable = hvg) ## could set other svd solvers
    sc.pp.neighbors(anna_data, n_pcs=30)
    sc.tl.umap(anna_data)
    # sc.tl.tsne(anna_data, n_jobs = 30)
    return(anna_data)

def  run_harmony(anna_data: sc.AnnData,
                 batch: str = "batch",
                 max_iter_harmony: int = 20,
                 resolution: float = 1.5) -> sc.AnnData:
    sc.external.pp.harmony_integrate(anna_data,
                                     batch,
                                     max_iter_harmony = max_iter_harmony)
    anna_data.obsm['X_pca'] = anna_data.obsm['X_pca_harmony']
    sc.pp.neighbors(anna_data, n_pcs=30)
    sc.tl.umap(anna_data)
    sc.tl.leiden(anna_data, resolution = resolution)
    return(anna_data)


def run_deg(adata: sc.AnnData,
        log2fc_min: float = 1,
        pval_cutoff: float = .05,
        method: str = 't-test_overestim_var',
        col_ind: str = "leiden", 
        use_raw: bool = False) -> pd.DataFrame:
    sc.tl.rank_genes_groups(adata,col_ind , pts = True, rankby_abs=True, method = method, use_raw = use_raw)
    groups =  adata.obs[col_ind].unique().sort_values()
    res_def =  [sc.get.rank_genes_groups_df(adata, group=group, log2fc_min = log2fc_min, pval_cutoff = pval_cutoff) for group in groups]
    for group in range(len(groups)):
        res_def[group]['cluster'] = groups[group]
    res_def = pd.concat(res_def)
    #res_def.to_csv('harmony_cellmarker_leiden.csv')
    return(res_def)


if __name__ == '__main__':
    # Pipline start from here
    # Define Macros
    import logging
    logging.basicConfig(filename='runpipeline.log',
                        encoding='utf-8', level=logging.INFO,
                        format='%(asctime)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    
    parser = OptionParser()
    
    parser.add_option("--path", dest="filename",
                      help="Input Directory", metavar="FILE")
    parser.add_option("--min_sample_size", default = 1000, 
                      help = "min cell size set 1,000 per sample", type = "int") 
    parser.add_option("--max_sample_size", default = 40000, 
                      help = "max cell size set 40,000 per sample", type = "int") 
    parser.add_option("--min_genes", default = 200, 
                      help = "min genes per sample set 200", type = "int")
    parser.add_option("--max_genes", default = 8000,
                      help = "max genes per sample set 8,000", type = "int")
    parser.add_option("--min_cells", default = 3, 
                      help = "min genes per sample set 3", type = "int")
    parser.add_option('--mt_pct', default = 15, help = 'MT-gene set 15%', type = 'int')
    parser.add_option("--doDoublet", action = "store_true",
                      help = "default set Doublet True", default = True)
    parser.add_option("--specie", default = "Human", type = 'string')
    parser.add_option("--join_type", default = "inner", type = "string", help = "default set inner join(inner/outer)")
    parser.add_option("--hvg", action='store_true',  help = "default set hvg as true", default = True)
    parser.add_option("--resolution", default = 1, 
                      help = "resolution set 1",type = "int")
    parser.add_option("--max_harmony_iter", default = 20, 
                      help = "max harmony iter set 20", type = "int") 

    # args = parser.parse_args()
    (options, args) = parser.parse_args()
    assert options.filename is not None, "Input should not set Null"

    # sequential main
    # load h5ad file
    logging.info("Start analysis.")
    logging.info("Load data.")
    adata_list   = load_data_from_path(path = options.filename)
    sample_name  = os.listdir(path=options.filename)
    logging.info(f'Directory: {options.filename}')
    sample_name  = np.char.strip(np.char.strip(sample_name, 'h5ad'), '.')
    logging.info(f'samples: {sample_name}')
    logging.info(f'do Doublet and Assmble.')
    adata        = assemble_data(anna_list   = adata_list,
                          sample_name = sample_name,
                          n_min       = options.min_sample_size,
                          n_max       = options.max_sample_size,
                          join_type   = options.join_type,
                          doDoublet   = options.doDoublet)
    logging.info(f'Delete doublet cells.')
    adata = adata[adata.obs['predicted_doublet'] == False ]
    # set species == "Human" 
    adata = addvar(anna_data = adata,
                #    vars      = None,
                   specie    = "human") 

    # process func
    logging.info(f'Cells and Genes filter start...')
    logging.info(f'#Cell: {adata.shape[0]}, #Gene: {adata.shape[1]}')
    stein_gate(anna_data = adata, 
               min_genes = options.min_genes,
               max_genes = options.max_genes,
               min_cells = options.min_cells,
               mt_pct    = options.mt_pct)
    logging.info(f'Filting is over.')
    logging.info(f'#Cell: {adata.shape[0]}, #Gene: {adata.shape[1]}')
    # adata.write_h5ad('preharmony.h5ad', compression = 'gzip')
    logging.info(f'Harmony integrated start...')
    adata = before_integrated(anna_data = adata, hvg = options.hvg)
    
    adata = run_harmony(anna_data        = adata, 
                        resolution       = options.resolution,
                        max_iter_harmony = options.max_harmony_iter)
    logging.info(f'Harmony integrated Over.')
    adata.write_h5ad("harmony.h5ad", compression  = 'gzip')    
