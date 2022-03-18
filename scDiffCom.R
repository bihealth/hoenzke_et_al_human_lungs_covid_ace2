library(Seurat)
library(scDiffCom)
library(dplyr)
library(plyr)
library(future)
library(dplyr)

options(future.globals.maxSize = 10 * 1024 * 1024^2)

lung_all <- readRDS('data/seurat/lung_tissue_combined_all.rds')
condition_map <- c('control'='control','akut'='acute','chronic'='prolonged','prolonged'='prolonged')
lung_all[['condition']] <- revalue(lung_all$condition, condition_map)
lung_all[['cluster']] <- Idents(lung_all)

DefaultAssay(lung_all) <- 'RNA'

explant <- subset(lung_all, subset=origin=='explant')
autopsy <- subset(lung_all, subset=origin=='autopsy')

plan(multicore, workers = 4)

CCI <- list()
for (cond in c('H3N2','MERS','SCoV1','SCoV2')) {
  scd <- run_interaction_analysis(
    seurat_object = subset(explant, subset=infect %in% c('control',cond)),
    LRI_species = "human",
    seurat_celltype_id = "cluster",
    iterations=10000,
    seurat_condition_id = list(
      column_name = "infect",
      cond1_name = cond,
      cond2_name = "control"
    )
  )
  CCI[[cond]] <- GetTableCCI(scd, "detected", simplified=TRUE) %>%
    dplyr::mutate(condition=cond,
                  origin='explant')
}
for (cond in c('acute','prolonged')) {
  scd <- run_interaction_analysis(
    seurat_object = subset(autopsy, subset=condition %in% c('control',cond)),
    LRI_species = "human",
    seurat_celltype_id = "cluster",
    iterations=10000,
    seurat_condition_id = list(
      column_name = "condition",
      cond1_name = cond,
      cond2_name = "control"
    )
  CCI[[cond]] <- GetTableCCI(scd, "detected", simplified=TRUE) %>%
    dplyr::mutate(condition=cond,
                  origin='autopsy')
}

write.csv(do.call(rbind,CCI),file.path('data','scDiffCom','CCI_tables.csv'))

