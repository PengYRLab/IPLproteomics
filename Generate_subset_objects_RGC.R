library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(harmony)

setwd("/Users/martinacavallini/Library/CloudStorage/Box-Box/DevelopmentalProteomics/Submission/NatureNeuroscience/UploadtoGithub/Analysis code")

#-------------------------- RGC Dataset ----------------------------

RGC <- readRDS("./Data/p5RGC_final_clusters_20180329.rds") #Load original dataset
count.data <- RGC@count.data #Extract counts from scR object
RGC_obj <- CreateSeuratObject(counts = count.data, project = "p5RGC") #Create Seurat object with counts
RGC_obj$percent.mt <- PercentageFeatureSet(RGC_obj, pattern = "^mt-") #Add mitochondrial gene metadata
RGC_obj$percent.rp <- PercentageFeatureSet(RGC_obj, pattern = "^Rp[sl]") #Add ribosomal gene metadata

RGC_obj <- NormalizeData(RGC_obj)
RGC_obj <- FindVariableFeatures(RGC_obj)
RGC_obj <- ScaleData(RGC_obj)
RGC_obj <- RunPCA(RGC_obj)
ElbowPlot(RGC_obj, ndims = 50)
RGC_obj <- RunHarmony(RGC_obj, group.by.vars = "orig.ident") #Perform Harmony integration for batch effects
RGC_obj <- RunUMAP(RGC_obj, reduction = "harmony", dims = 1:25)
RGC_obj <- FindNeighbors(RGC_obj, reduction = "harmony", dims = 1:25)
RGC_obj <- FindClusters(RGC_obj, resolution = 1)
DimPlot(RGC_obj, reduction = "umap", group.by = "orig.ident", label = T, repel = T) 

#Some batch effects still remain. Remove those samples:
Idents(RGC_obj) <- "orig.ident"
RGC_obj <- subset(RGC_obj, idents = c("p5RGCS1", "p5RGCS10", "p5RGCS11", "p5RGCS12", "p5RGCS13", "p5RGCS2", "p5RGCS3", "p5RGCS4", "p5RGCS5", "p5RGCS6", "p5RGCS7", "p5RGCS8", "p5RGCS9"))
RGC_obj <- NormalizeData(RGC_obj)
RGC_obj <- FindVariableFeatures(RGC_obj)
RGC_obj <- ScaleData(RGC_obj)
RGC_obj <- RunPCA(RGC_obj)
RGC_obj <- RunHarmony(RGC_obj, group.by.vars = "orig.ident")
RGC_obj <- RunUMAP(RGC_obj, reduction = "harmony", dims = 1:25)
RGC_obj <- FindNeighbors(RGC_obj, reduction = "harmony", dims = 1:25)
RGC_obj <- FindClusters(RGC_obj, resolution = 1)
DimPlot(RGC_obj, reduction = "umap", group.by = "orig.ident", label = T)
DimPlot(RGC_obj, reduction = "umap", label = T)

#Perform label transfer with adult RGC dataset as reference
reference <- readRDS("./Data/adultrgc_object.rds")
reference <- FindVariableFeatures(reference) #Required for FindTransferAnchors below
rgc.anchors <- FindTransferAnchors(reference = reference, query = RGC_obj, dims = 1:90, reference.reduction = "pca")
predictions <- TransferData(anchorset = rgc.anchors, refdata = reference$cluster_id, dims = 1:90)
RGC_obj <- AddMetaData(RGC_obj, metadata = predictions)

#Using this mapping strategy, we obtain all RGC types except 39_Novel
DimPlot(RGC_obj, reduction = "umap", group.by = "predicted.id", label = T) +NoLegend()
DotPlot(RGC_obj, group.by = "predicted.id", 
        features = c("Eomes", "Irx3", "Tbr1", "Mafb", "Foxp2", "Neurod2", "Tfap2d", "Bnc2"),
        cluster.idents = T)

saveRDS(RGC_obj, "./Saved_objects/p5rgc_object_allgenes.rds")


#Generate IPL-CSP subset object
IPL_CSP <- read.xlsx("./Data/IPL_CSP_list.xlsx")
genes_vector <- unlist(strsplit(as.character(IPL_CSP$Gene.Names), split = " "))
genes_vector <- unique(trimws(genes_vector))

RGC_subset <- subset(RGC_obj, features = genes_vector) #Subset features to only include IPL-CSP genes
RGC_subset <- NormalizeData(RGC_subset)
RGC_subset <- FindVariableFeatures(RGC_subset)
RGC_subset <- ScaleData(RGC_subset)
RGC_subset <- RunPCA(RGC_subset, npcs = 50)
RGC_subset <- RunHarmony(RGC_subset, group.by.vars = "orig.ident")
RGC_subset <- RunUMAP(RGC_subset, reduction = "harmony", dims = 1:25)
RGC_subset <- FindNeighbors(RGC_subset, reduction = "harmony", dims = 1:25)
RGC_subset <- FindClusters(RGC_subset, resolution = 0.5)
DimPlot(RGC_subset, reduction = "umap", group.by = "orig.ident", label = T)
DimPlot(RGC_subset, reduction = "umap", group.by = "seurat_clusters", label = T)

#Some cells have very low IPL-CSP gene expression. Remove seurat clusters 0 and 17
Idents(RGC_subset) <- "seurat_clusters"
RGC_subset <- subset(RGC_subset, idents = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                            11, 12, 13, 14, 15, 16, 18, 19, 20))

#Final analysis
RGC_subset <- NormalizeData(RGC_subset)
RGC_subset <- FindVariableFeatures(RGC_subset)
RGC_subset <- ScaleData(RGC_subset)
RGC_subset <- RunPCA(RGC_subset, npcs = 50)
RGC_subset <- RunHarmony(RGC_subset, group.by.vars = "orig.ident")
RGC_subset <- RunUMAP(RGC_subset, reduction = "harmony", dims = 1:25)
RGC_subset <- FindNeighbors(RGC_subset, reduction = "harmony", dims = 1:25)
RGC_subset <- FindClusters(RGC_subset, resolution = 0.5)
DimPlot(RGC_subset, reduction = "umap", group.by = "orig.ident", label = T)
DimPlot(RGC_subset, reduction = "umap", group.by = "predicted.id", label = T) + NoLegend()

saveRDS(RGC_subset, "./Saved_objects/p5rgc_object_IPLCSP.rds")


#Generate UMAP plot in Supplementary Figure 5 panel D (all genes but same final cells as p5rgc_object_IPLCSP.rds)
RGC_subset <- readRDS("./Saved_objects/p5rgc_object_IPLCSP.rds")
RGC_obj <- readRDS("./Saved_objects/p5rgc_object_allgenes.rds")
RGC_obj_subset <- RGC_obj[,colnames(RGC_obj) %in% colnames(RGC_subset)]

RGC_obj_subset <- NormalizeData(RGC_obj_subset)
RGC_obj_subset <- FindVariableFeatures(RGC_obj_subset)
RGC_obj_subset <- ScaleData(RGC_obj_subset)
RGC_obj_subset <- RunPCA(RGC_obj_subset)
RGC_obj_subset <- RunHarmony(RGC_obj_subset, group.by.vars = "orig.ident")
RGC_obj_subset <- RunUMAP(RGC_obj_subset, reduction = "harmony", dims = 1:25)
RGC_obj_subset <- FindNeighbors(RGC_obj_subset, reduction = "harmony", dims = 1:25)
RGC_obj_subset <- FindClusters(RGC_obj_subset, resolution = 1)
DimPlot(RGC_obj_subset, reduction = "umap", group.by = "orig.ident", label = T)
DimPlot(RGC_obj_subset, reduction = "umap", group.by = "predicted.id", label = T) + NoLegend()


saveRDS(RGC_obj_subset, "./Saved_objects/p5rgc_object_allgenes_finalcells.rds")

