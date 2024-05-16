# conda activate /oak/stanford/groups/smontgom/amarder/micromamba/envs/seurat5v2

dir="/oak/stanford/groups/smontgom/amarder/Trisomy"
start="1"
end="4"
DATASET="Trisomy18"
SUFFIX=".PRE_CB"

# load
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(Seurat)
# library(Signac)
# library(GenomicRanges)
# library(future)
library(data.table)
library(harmony)
# library("JASPAR2020")
# library("TFBSTools")
# library(parallel)
library(fgsea)

# library("glmGamPoi") # for SCTransform

source(paste0(dir,"/scripts/single_cell_processing/helper_func/","integration_function_sheet.R"))


f <- paste0(dir,"/output/data/",DATASET,"/RNA_FindClusters",SUFFIX,".rds")
print(paste0("Loading data: ",f,"..."))
dfcombined <- readRDS(file = f)
colName = "clusters_v2"
dfcombined@meta.data[,colName] = NA

RunSubcluster = function(cellName,clusters,kclust,colName) {
  i = dfcombined$seurat_clusters %in% clusters
  dftmp = dfcombined[,i]
  kmeans_results <- kmeans(Embeddings(dftmp,reduction="harmony"),kclust)
  kmeans_results1 <- data.frame(kmeans_results$cluster)
  dfcombined@meta.data[,colName][i] <- paste0(cellName,"_",paste(clusters,collapse = ','),"_",kmeans_results1[,1])
  return(dfcombined)
}

cellName = "HSCs"; clusters=c(5,6,23); kclust = 5
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "unknown"; clusters=c(9); kclust = 4
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "MEMPs"; clusters=c(2,8,11); kclust = 5
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "erythroid"; clusters=c(1); kclust = 3
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "megakaryocytes"; clusters=c(10,18,29); kclust = 4
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "unknown 2"; clusters=c(19,26); kclust = 5
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "mast cells"; clusters=c(25,27); kclust = 3
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "mast cells"; clusters=c(22); kclust = 2
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "pDCs"; clusters=c(24,28); kclust = 5
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "granulocyte progenitors"; clusters=c(12); kclust = 5
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "mono_granulocyte progenitors"; clusters=c(15); kclust = 5
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "pre-pro B cells"; clusters=c(21); kclust = 3
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "pro B cells"; clusters=c(14); kclust = 3
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "pro B cells 2"; clusters=c(16); kclust = 3
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "pro B cells 2"; clusters=c(13); kclust = 3
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "B cells"; clusters=c(13); kclust = 3
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "some sort of B cells 1"; clusters=c(17); kclust = 3
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "some sort of B cells 2"; clusters=c(7); kclust = 3
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "monocytes"; clusters=c(3); kclust = 5
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "cDC2"; clusters=c(4); kclust = 5
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "NK cells"; clusters=c(20); kclust = 3
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

cellName = "T cells"; clusters=c(30); kclust = 3
dfcombined = RunSubcluster(cellName,clusters,kclust,colName)

dfcombined@meta.data[is.na(dfcombined@meta.data[,colName]),colName] = paste0("Previous.",dfcombined@meta.data[is.na(dfcombined@meta.data[,colName]),"seurat_clusters"])

cellTypes = unique(dfcombined@meta.data[,colName])
dir.create(paste0(dir,"/output/data/",DATASET,"/",colName))
print(paste0(dir,"/output/data/",DATASET,"/",colName))
for (cellType in cellTypes) {
  print(cellType)
  cellNames1 <- rownames(dfcombined@meta.data)[dfcombined@meta.data[,colName] == cellType]
  f.out <- paste0(dir,"/output/data/",DATASET,"/",colName,"/RNA_umap.",cellType,".png")
  png(filename=f.out,width = 5000,height=4000,res=500)
  print(DimPlot(dfcombined, label=T,  cells.highlight= list(cellNames1), cols.highlight = c("red"), cols= "grey") & theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
}

genes_to_use = c("CD34","SPINK2","PROM1","MLLT3",
                 "GATA1","KLF1","TESPA1","AHSP",
                 "ALAS2","HBA1","GYPA","GATA2","HDC",
                 "CPA3","ITGA2B","GP9",
                 "MPO","AZU1","SPI1","LYZ","CD14","CD68",
                 "S100A9","MNDA","FCN1","CD163","MS4A7","VCAN","VCAM1","MARCO",
                 "C1QA","C1QB","C1QC",
                 "CTSB","NKG7","PRF1",
                 "GZMA",
                 "TIGIT","TRAC","CD3D","RORC", # T CELLS
                 "IL2RB",
                 "IL7R","DHFR","PAX5","MME","IGLL1","IGHM",
                 "CD79A","CD19","JCHAIN","IRF8","CLEC4C","IL3RA",
                 "CD1C","CLEC4A",
                 "CLEC10A", #cdc1
                 "CLEC9A" , "THBD","XCR1", "BATF3", 
                 "MKI67")

dfcombined@meta.data[,colName] = as.factor(dfcombined@meta.data[,colName])
levels(dfcombined@meta.data[,colName]) = sort(unique(dfcombined@meta.data[,colName]))
Idents(dfcombined) = colName

f.out <- paste0(dir,"/output/data/",DATASET,"/",colName,"/RNA_dotplot.png")
png(filename=f.out,width = 9000,height=9000,res=500)
print(DotPlot(object = dfcombined, features = genes_to_use)) + RotatedAxis()
dev.off()

dir.create(paste0(dir,"/output/data/",DATASET,"/",colName,"/sample_specific_dotplots"))
for (sample in unique(dfcombined$dataset)) {
  print(paste0("Dotplot for sample: ",sample," ..."))
  f.out <- paste0(dir,"/output/data/",DATASET,"/",colName,"/sample_specific_dotplots/RNA_dotplot.",sample,".png")
  png(filename=f.out,width = 9000,height=9000,res=500)
  print(DotPlot(object = subset(dfcombined,dataset==sample), features = genes_to_use) + RotatedAxis())
  dev.off()
}

f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_",colName,SUFFIX,".rds")
print(paste0("Saving data: ",f.out," ..."))
saveRDS(dfcombined,file = f.out)


######################

dfcombined = subset(dfcombined,seurat_clusters!=9)

print("NormalizeData + FindVariableFeatures + ScaleData + RunPCA...")
dfcombined <- NormalizeData(dfcombined,normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  # FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData() %>%
  RunPCA(npcs=50)

print("RunHarmony...")
lambda_val=1; dfcombined <- RunHarmony(dfcombined,"dataset",lambda=lambda_val,dims=1:30)

print("FindNeighbors...")
dfcombined <- FindNeighbors(object = dfcombined, reduction = 'harmony', dims = 1:30)

print("RunUMAP...")
dfcombined <- RunUMAP(dfcombined, reduction = "harmony", dims = 1:30)

Idents(dfcombined) = "seurat_clusters"

# Save data:
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_FindClusters_NoCluster9",SUFFIX,".rds")
print(paste0("Saving data: ",f.out,"..."))
saveRDS(dfcombined,file = f.out)

print("Creating UMAP with seurat clusters...")
# p1 <- DimPlot(dfcombined, group.by = "seurat_clusters")
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_umap.seurat_clusters.cluster_9_removed",SUFFIX,".png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(dfcombined, group.by = "seurat_clusters") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "seurat_clusters"))
dev.off()




