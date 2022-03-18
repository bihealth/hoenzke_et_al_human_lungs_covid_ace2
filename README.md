# Code repository for Hönzke et al., "Human lungs show limited permissiveness for SARS-CoV-2 due to scarce ACE2 levels but virus-induced expansion of inflammatory macrophages"

## data access

processed data is available from NCBI GEO under accession [GSE198864](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198864)  

- scRNAseq: cellranger output (raw h5 files) for
  - lung tissue explants (GSM5958261-88)
  - COVID autopsy material (GSM5958253-60)
  - ex vivo infected alveolar macrophages (GSM5958229-52)
  - lung organoids (GSM5958290-99)
- bulk RNAseq: counts and metadata for lung tissue explants (GSM5958300-19)

and should be downloaded and stored in the `cellranger` and `bulk_RNAseq` subfolders like so

```
data
├── AM
│   └── cellranger
├── autopsy
│   └── cellranger
├── explant
│   ├── bulk_RNAseq
│   └── cellranger
└── organoids
    └── cellranger
```

final R objects are also available and should go into `data/seurat`

## CellBender

[CellBender](https://github.com/broadinstitute/CellBender) `remove-background` is used to remove viral RNA background especially from autopsy single-nuc samples

it could be run like so:

```
for material in AM autopsy explant organoids; do
    while read -r sample ncells; do
        bash run_cellbender.sh ${sample} ${ncells} ${material}
    done < expected_cells_${material}.txt
done
```

## scrublet

[scrublet](https://github.com/swolock/scrublet) is used as a first doublet-removal step and could be run like so

```
for material in AM autopsy explant organoids; do
    mkdir -p data/${material}/scrublet
    for h5 in data/cellbender/*.h5; do
        sample=$(basename ${h5} _filtered.h5)
        python run_scrublet.py -i data/${material}/cellbender/${h5} -o data/${material}/scrublet/${sample}_scrublet.csv
    done
done
```

## Seurat processing

R code is split into several markdowns that can be simply knitted

pre-processing:

- `process_explant_controls.Rmd`
- `process_autopsy_controls.Rmd`
- `process_explant_all.Rmd`
- `process_autopsy_all.Rmd`
- `process_alveolar_organoids.Rmd`
- `process_bronchial_organoids.Rmd`
- `process_ex_vivo_AM.Rmd`

integration:

- `combine_lung_tissue_controls.Rmd`
- `combine_lung_tissue_all.Rmd`

differential expression:

- `run_differential_expression.Rmd`

macrophage subclustering:

- `recluster_macrophages.Rmd`

## SCENIC

[SCENIC](https://github.com/aertslab/SCENIC) is run on explant and autopsy macrophages separately

1. prepare loom input files `python prepare_scenic.py`
2. run scenic: `cd data/scenic; bash ../../run_scenic.sh`
3. process results: `python process_scenic.py`

## scDiffCom

[scDiffCom](https://github.com/CyrilLagger/scDiffCom) is used to investigate intercellular receptor-ligand signaling, see `scDiffCom.R`

## external datasets

lung cell atlas reference

- `droplet_normal_lung_seurat_ntiss10x.P2.anno.20191002.RC4.Robj` from [here](https://www.synapse.org/#!Synapse:syn21560412) converted to Seurat V3 and stored as RDS

for the comparison with macrophages from other datasets, we downloaded and processed data as described in the methods section

- Delorey et al.: from the Broad Single Cell Portal under accession [SCP1052](https://singlecell.broadinstitute.org/single_cell/study/SCP1052)
- Grant et al.: from NCBI GEO under accession [GSE155249](http://www.ncbi.nlm.nih.gov/geo/query.acc.cgi?acc=GSE155249) (`GSE155249_supplement.h5ad`)
- Liao et al.: from NCBI GEO under accession [GSE155249](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926)
- Melms et al.: from the Broad Single Cell Portal under accession [SCP1219](https://singlecell.broadinstitute.org/single_cell/study/SCP1219)
- Wendisch et al.: see `BAL_macrophages.rds` [here](https://nubes.helmholtz-berlin.de/s/XrM8igTzFTFSoio)

## read mapping statistics

note that bam files from which read statistics are derived cannot be provided due to data privacy restrictions; statistics has been extracted from bam files using `get_subgenomic_reads.py` and `get_coverage.sh`

## paper figures

all code for paper figures is in `paper_figures.Rmd`
