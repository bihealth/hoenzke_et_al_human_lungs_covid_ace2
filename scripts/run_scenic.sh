#!/bin/bash -e

arboreto_with_multiprocessing.py momac_${orig}.loom hs_hgnc_curated_tfs.txt --method grnboost2 --output momac_${orig}_adj.tsv --num_workers 4 --seed 1
pyscenic ctx momac_${orig}_adj.tsv *.feather --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl --output momac_${orig}_reg.csv --mask_dropouts --num_workers 4 --expression_mtx_fname momac_${orig}.loom 
pyscenic aucell momac_${orig}.loom momac_${orig}_reg.csv --output momac_${orig}_pyscenic.loom --num_workers 4 
