import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import loompy as lp
from scipy.spatial.distance import jensenshannon
import json
import zlib
import base64
from math import sqrt

anno=pd.read_csv('data/seurat/lung_combined_macrophages_meta.csv',header=0,index_col=0)

auc_mtx={}
regulons={}
for orig in ['explant','autopsy']:
    lf = lp.connect('data/scenic/momac_'+orig+'_pyscenic.loom', mode='r+', validate=False )
    auc_mtx[orig] = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
    regulons[orig] = lf.ra.Regulons
    lf.close()

for orig in ['explant','autopsy']:
    auc_mtx[orig].columns = auc_mtx[orig].columns.str.replace('\(','_(')
    regulons[orig].dtype.names = tuple( [ x.replace("(","_(") for x in regulons[orig].dtype.names ] )
    auc_mtx[orig].to_csv('data/scenic/momac_scenic_'+orig+'_AUC.csv')
    pd.DataFrame(regulons[orig]).to_csv('data/scenic/momac_scenic_'+orig+'_regulons.csv')

def regulon_specificity_scores(auc_mtx, cell_type_series):
    """
    Calculates the Regulon Specificty Scores (RSS). [doi: 10.1016/j.celrep.2018.10.045]
    :param auc_mtx: The dataframe with the AUC values for all cells and regulons (n_cells x n_regulons).
    :param cell_type_series: A pandas Series object with cell identifiers as index and cell type labels as values.
    :return: A pandas dataframe with the RSS values (cell type x regulon).
    """

    cell_types = list(cell_type_series.unique())
    n_types = len(cell_types)
    regulons = list(auc_mtx.columns)
    n_regulons = len(regulons)
    rss_values = np.empty(shape=(n_types, n_regulons), dtype=np.float)

    def rss(aucs, labels):
        # jensenshannon function provides distance which is the sqrt of the JS divergence.
        return 1.0 - jensenshannon(aucs/aucs.sum(), labels/labels.sum())

    for cidx, regulon_name in enumerate(regulons):
        for ridx, cell_type in enumerate(cell_types):
            rss_values[ridx, cidx] = rss(auc_mtx[regulon_name], (cell_type_series == cell_type).astype(int))

    return pd.DataFrame(data=rss_values, index=cell_types, columns=regulons)


for orig in ['explant','autopsy']:
    rss_subcluster=regulon_specificity_scores(auc_mtx[orig], anno.loc[auc_mtx[orig].index,'subcluster']).T
    rss_contrast=regulon_specificity_scores(auc_mtx[orig], anno.loc[auc_mtx[orig].index,'contrast']).T
    pd.concat([rss_subcluster,rss_contrast],axis=1).to_csv('data/scenic/momac_'+orig+'_regulon_specificity_scores.csv')




