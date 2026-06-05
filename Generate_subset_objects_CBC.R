library(Seurat)
library(dplyr)
library(tidyr)
library(tidyverse)
library(tibble)
library(readxl)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(harmony)
library(cowplot)
library(patchwork)

setwd("/Users/martinacavallini/Library/CloudStorage/Box-Box/DevelopmentalProteomics/Submission/NatureNeuroscience/UploadtoGithub/Analysis code")

#-------------------------- BC Dataset ----------------------------

load("./Data/seu.BPtype.mouse.lab.harmony.rda")
Idents(seu) <- "cell_type" #Keep only cone bipolar cells (CBCs)
CBC_obj <- subset(seu, idents = c("BC5A", "BC7", "BC6", "BC5C", "BC1A", "BC3B", "BC1B", "BC2", "BC5D", "BC3A", "BC5B", "BC4", "BC8-9"))

CBC_obj <- NormalizeData(CBC_obj)
CBC_obj <- FindVariableFeatures(CBC_obj)
CBC_obj <- ScaleData(CBC_obj)
CBC_obj <- RunPCA(CBC_obj)
ElbowPlot(CBC_obj, ndims = 50)
CBC_obj <- RunHarmony(CBC_obj, group.by.vars = "orig.ident") #Perform Harmony integration for batch effects
CBC_obj <- RunUMAP(CBC_obj, reduction = "harmony", dims = 1:10)
CBC_obj <- FindNeighbors(CBC_obj, reduction = "harmony", dims = 1:10)
CBC_obj <- FindClusters(CBC_obj, resolution = 1)
DimPlot(CBC_obj, reduction = "umap", group.by = "orig.ident", label = T, repel = T)
DimPlot(CBC_obj, reduction = "umap", group.by = "cell_type", label = T, repel = T)

saveRDS(CBC_obj, "./Saved_objects/CBC_object_allgenes.rds")


#Generate IPL-CSP subset object
IPL_CSP <- read.xlsx("./Data/IPL_CSP_list.xlsx")
genes_vector_IPLCSP <- unlist(strsplit(as.character(IPL_CSP$Gene.Names), split = " "))
genes_vector_IPLCSP <- unique(trimws(genes_vector_IPLCSP))

CBC_subset <- subset(CBC_obj, features = genes_vector_IPLCSP) #Subset features to only include IPL-CSP genes
CBC_subset <- NormalizeData(CBC_subset)
CBC_subset <- FindVariableFeatures(CBC_subset)
CBC_subset <- ScaleData(CBC_subset)
CBC_subset <- RunPCA(CBC_subset, npcs = 50)
CBC_subset <- RunHarmony(CBC_subset, group.by.vars = "orig.ident")
CBC_subset <- RunUMAP(CBC_subset, reduction = "harmony", dims = 1:10)
CBC_subset <- FindNeighbors(CBC_subset, reduction = "harmony", dims = 1:10)
CBC_subset <- FindClusters(CBC_subset, resolution = 0.5)
DimPlot(CBC_subset, reduction = "umap", group.by = "orig.ident", label = T)
DimPlot(CBC_subset, reduction = "umap", group.by = "cell_type", label = T)

saveRDS(CBC_subset, "./Saved_objects/CBC_object_IPLCSP.rds")


#Generate Module subset object
CBC_mods <- read.xlsx("./Data/CBC_modules.xlsx") #Load CBC module genes
genes_vector_Mods <- CBC_mods$gene_name

CBC_subset <- subset(CBC_obj, features = genes_vector_Mods) #Subset features to only include CBC module genes
CBC_subset <- NormalizeData(CBC_subset)
CBC_subset <- FindVariableFeatures(CBC_subset)
CBC_subset <- ScaleData(CBC_subset)
CBC_subset <- RunPCA(CBC_subset, npcs = 50)
CBC_subset <- RunHarmony(CBC_subset, group.by.vars = "orig.ident")
CBC_subset <- RunUMAP(CBC_subset, reduction = "harmony", dims = 1:10)
CBC_subset <- FindNeighbors(CBC_subset, reduction = "harmony", dims = 1:10)
CBC_subset <- FindClusters(CBC_subset, resolution = 0.5)
DimPlot(CBC_subset, reduction = "umap", group.by = "orig.ident", label = T)
DimPlot(CBC_subset, reduction = "umap", group.by = "cell_type", label = T)

saveRDS(CBC_subset, "./Saved_objects/CBC_object_Mods.rds")
