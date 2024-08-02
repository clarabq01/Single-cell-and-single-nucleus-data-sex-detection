Matson309_ctr_counts <- Read10X_h5(
  filename = "Seurat_Matson309_inj_FPR_0.01.h5",
  use.names = TRUE,
  unique.features = TRUE) 
Matson309_ctr<- CreateSeuratObject(counts = Matson309_ctr_counts, min.cells = 3, min.features = 100)
head(Matson309_ctr@meta.data)


hembra = c("Tsix", "Xist")
macho = c("Uty", "Ddx3y", "Eif2s3y", "Kdm5d")

genehembra <- FetchData(Matson309_ctr, vars= hembra)
totalgenehembra <- colSums(genehembra)
totalgenehembra
totalcellhembra <- mean(Matson309_ctr@meta.data$nCount_RNA) * length(colnames(Matson309_ctr))
totalcellhembra
sexadohembra = totalgenehembra/totalcellhembra * 1000000
sexadohembra
sum(sexadohembra)


genemacho <- FetchData(Matson309_ctr, vars= macho)
totalgenemacho <- colSums(genemacho)
totalgenemacho
totalcellmacho <- mean(Matson309_ctr@meta.data$nCount_RNA) * length(colnames(Matson309_ctr))
totalcellmacho
sexadomacho = totalgenemacho/totalcellmacho * 1000000
sexadomacho
sum(sexadomacho)


result = -sum(sexadohembra)-sum(sexadomacho)
if(result>0){'MACHO'} else {'HEMBRA'}

