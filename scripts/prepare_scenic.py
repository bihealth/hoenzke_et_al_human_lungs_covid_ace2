import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import loompy as lp

cells=pd.read_csv('data/seurat/lung_combined_macrophages_barcodes.txt',header=None).squeeze()
samples=cells.str.split('_').str[0].unique()
anno=pd.read_csv('data/seurat/lung_combined_macrophages_meta.csv',header=0,index_col=0)

ads=[]
for lib in samples:
    lib2=lib.replace('-','_') if lib.startswith('lung') else lib
    h5_file=os.path.join('data','explant' if 'lung' in lib2 else 'autopsy','cellbender',lib2+'_filtered.h5')
    if not os.path.isfile(h5_file):
        print(h5_file)
        continue
    ad=sc.read_10x_h5(h5_file)
    ad.X=scipy.sparse.csr_matrix(ad.X)
    ad.obs.index=lib+'_'+ad.obs.index.astype(str)
    ad.obs['sample']=lib
    ad.var_names_make_unique()
    ads.append(ad)

ad=ads[0].concatenate(ads[1:],batch_key=None,index_unique=None)

for orig in ['explant','autopsy']:
    ad_here=ad[anno.index[anno['origin']==orig]]
    row_attrs = {
        "Gene": np.array(ad_here.var_names) ,
    }
    col_attrs = {
        "CellID": np.array(ad_here.obs_names) ,
        "sample": np.array(ad_here.obs['sample']),
        "nGene": np.array( np.sum(ad_here.X.transpose()>0 , axis=0)).flatten() ,
        "nUMI": np.array( np.sum(ad_here.X.transpose() , axis=0)).flatten() ,
    }
    lp.create('data/scenic/momac_'+orig+'.loom', ad_here.X.transpose(), row_attrs, col_attrs)

