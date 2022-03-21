import scrublet as scr
import pandas as pd
import scanpy as sc
import numpy as np
from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument('-i',dest='input',help="""input file (CellBender h5)""")
parser.add_argument('-o',dest='output',help="""output file (<sample>_scrublet.csv)""")

args=parser.parse_args()

ad=sc.read_10x_h5(args.input)
scrub=scr.Scrublet(ad.X, expected_doublet_rate=0.05)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
pd.DataFrame(dict(scrublet_score=doublet_scores,
                  scrublet_prediction=np.array(['singlet','doublet'])[predicted_doublets.astype(int)]),
             index=ad.obs_names).to_csv(args.output)

