library(DoubletFinder)

my_doubletFinder <- function (seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE) {
  require(Seurat)
  require(fields)
  require(KernSmooth)
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, "DF.classifications"] <- classifications
    return(seu)
  }
  if (reuse.pANN == FALSE) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...", 
                sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets) %>%
      NormalizeData(normalization.method = 'LogNormalize',
                    scale.factor = 10000, verbose=FALSE) %>%
      FindVariableFeatures(selection.method = 'vst', nfeatures = 2000, verbose=FALSE) %>%
      ScaleData(verbose=FALSE) %>%
      RunPCA(features = VariableFeatures(.),verbose = FALSE)
    pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, PCs]
    cell.names <- rownames(seu_wdoublets@meta.data)
    nCells <- length(cell.names)
    rm(seu_wdoublets)
    gc()
    dist.mat <- fields::rdist(pca.coord)
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                 ncol = 1))
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      neighbor.names <- rownames(dist.mat)[neighbors]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
    }
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, "pANN"] <- pANN[rownames(seu@meta.data), 1]
    seu@meta.data[, "DF.classifications"] <- classifications
    return(seu)
  }
}