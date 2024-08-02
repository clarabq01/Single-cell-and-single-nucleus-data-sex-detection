library(tidyverse)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(SeuratData)
library(scDblFinder)
library(presto)
library(SingleCellExperiment)
library(dplyr)
library(scran)
library(patchwork)
library(viridis)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(ggforce)
library(gghalves)
library(ggridges)
library(RCurl)
library(glmGamPoi)
# library(BPCells) # BPCells is an R package that allows for computationally efficient single-cell analysis.BPCells allows us to easily analyze these large datasets in memory (ERROR)
library(AnnotationHub)
library(ensembldb)
library(reticulate)
library(sctransform)
library(SingleR)
library(scales)
library(cowplot)

Matson309_ctr_counts <- Read10X_h5(
  filename = "Seurat_Matson309_inj_FPR_0.01.h5",
  use.names = TRUE,
  unique.features = TRUE) 
Matson309_ctr<- CreateSeuratObject(counts = Matson309_ctr_counts, min.cells = 3, min.features = 100)
head(Matson309_ctr@meta.data)

calculate_sex <- function(seurat_object) {
  hembra <- c("Tsix", "Xist")
  macho <- c("Uty", "Ddx3y", "Eif2s3y", "Kdm5d")
  
  gene_hembra <- FetchData(seurat_object, vars = hembra)
  total_gene_hembra <- colSums(gene_hembra)
  total_cell_hembra <- mean(seurat_object@meta.data$nCount_RNA) * length(colnames(seurat_object))
  sexado_hembra <- total_gene_hembra / total_cell_hembra * 1000000
  
  gene_macho <- FetchData(seurat_object, vars = macho)
  total_gene_macho <- colSums(gene_macho)
  total_cell_macho <- mean(seurat_object@meta.data$nCount_RNA) * length(colnames(seurat_object))
  sexado_macho <- total_gene_macho / total_cell_macho * 1000000
  
  result <- -sum(sexado_hembra) - sum(sexado_macho)
  if (result > 0) {
    print("MACHO")
  } else {
    print("HEMBRA")
  }
  
  print("Sexado Hembra:")
  print(sexado_hembra)
  print("Sexado Macho:")
  print(sexado_macho)
}

sex <- calculate_sex(Matson309_ctr)

# [1] "HEMBRA"
# [1] "Sexado Hembra:"
# Tsix      Xist 
# 235.6534 4183.7115 
# [1] "Sexado Macho:"
# Uty 
# 1.382131 




