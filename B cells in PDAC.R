################################################################
############ Reading in all samples in a shorter script ########
################################################################
library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(hdf5r)
library(Matrix)
rm(list = ls())

####### to read with the files not in zipped gz format

#########################################################
###### For Normal Samples from GSA:PRJCA001063 ##########
#########################################################

dir <- c("/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/N1", 
         "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/N2", 
         "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/N3", 
         "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/N4", 
         "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/N5", 
         "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/N6", 
         "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/N7", 
         "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/N8", 
         "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/N9", 
         "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/N10", 
         "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/N11")
samples_name = c('N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'N10', 'N11')


#########################################################
###### For PDAC Samples from GSA:PRJCA001063   ##########
#########################################################

dir <- c(
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T1",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T2",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T3",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T4",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T5",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T6",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T7",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T8",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T9",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T10",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T11",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T12",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T13",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T14",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T15",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T16",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T17",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T18",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T19",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T20",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T21",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T22",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T23",
  "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/PDAC/T24"
)

# Sample names corresponding to the directories
samples_name <- c('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T20', 'T21', 'T22', 'T23', 'T24')

#########################################################
######   For NAT Samples from GEO:GSE155698    ##########
#########################################################

dir <- c("/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/NAT_pan/NAT_1/filtered_feature_bc_matrix", 
         "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/NAT_pan/NAT_2/filtered_feature_bc_matrix", 
         "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/PDAC_scRNA/NAT_pan/NAT_3/filtered_feature_bc_matrix")
samples_name = c('NAT_1', 'NAT_2', 'NAT_3')

#########################################################
######## Initialize list to store Seurat objects ########
#########################################################
scRNAlist <- list()

for (i in 1:length(dir)) {
  data.dir <- dir[i]
  
  # Check if files are gzipped or not and read them accordingly
  if (file.exists(file.path(data.dir, "barcodes.tsv.gz"))) {
    counts <- Read10X(data.dir = data.dir)
  } else {
    barcode.path <- file.path(data.dir, "barcodes.tsv")
    #features.path <- file.path(data.dir, "genes.tsv") # for pdac
    features.path <- file.path(data.dir, "features.tsv") # for NAT
    matrix.path <- file.path(data.dir, "matrix.mtx")
    
    # Read the matrix
    matrix <- Matrix::readMM(matrix.path)
    
    # Read the features (genes)
    features <- read.delim(features.path, header = FALSE)
    gene_names <- features$V2  # Assuming the second column contains the gene names
    gene_names <- gsub("_", "-", gene_names)  # Replace underscores with dashes
    
    # Ensure gene names are unique
    unique_gene_names <- make.unique(gene_names)
    
    # Read the barcodes
    barcodes <- read.delim(barcode.path, header = FALSE)
    barcodes <- barcodes$V1
    
    # Assign row and column names to the matrix
    rownames(matrix) <- unique_gene_names
    colnames(matrix) <- barcodes
    
    counts <- matrix
  }
  
  # Create Seurat object
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = samples_name[i], min.cells = 3, min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])
  
  # Calculate percentage of mitochondrial genes
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  
  # Calculate percentage of ribosomal genes
  scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")
  
  # Calculate percentage of hemoglobin genes
  HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
  scRNAlist[[i]][["percent.HB"]] <- PercentageFeatureSet(scRNAlist[[i]], features = HB.genes)
}

# Print the list of Seurat objects
print(scRNAlist)

#########################################################
########   Merge and Normalize Seurat objects    ########
#########################################################

##############################################################
###### For 11 Normal Samples from GSA:PRJCA001063   ##########
##############################################################
normal_panc <- Reduce(function(x, y) merge(x, y), scRNAlist)

#saveRDS(normal_panc, file = "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/Nataly-Objects/normal_panc.rds")

######### Normalize data and perform harmony batch effect correction
minGene=200
maxUMI=25000
pctMT=10

VlnPlot(normal_panc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
normal_panc <- subset(normal_panc, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & percent.mt < pctMT)
VlnPlot(normal_panc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

normal_panc <- NormalizeData(normal_panc) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
normal_panc <- RunPCA(normal_panc, verbose = F)
ElbowPlot(normal_panc, ndims = 50)
pc.num=1:25
normal_panc <- normal_panc %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
normal_panc <- FindNeighbors(normal_panc, dims=pc.num)
normal_panc <- FindClusters(normal_panc)
normal_panc <- RunUMAP(normal_panc, dims = 1:25)
DimPlot(normal_panc, label = T)
DimPlot(normal_panc, label = T, group.by = "orig.ident")

normal_panc <- RunPCA(normal_panc, npcs=50, verbose=FALSE)
normal_panc <- RunHarmony(normal_panc, group.by.vars="orig.ident", max.iter.harmony = 20) 
normal_panc <- RunUMAP(normal_panc, reduction = "harmony", dims = 1:25)
normal_panc <- FindNeighbors(normal_panc, reduction = "harmony", dims = 1:25) 
normal_panc <- FindClusters(normal_panc, verbose = FALSE, resolution = 0.4)
DimPlot(normal_panc, label = T)
DimPlot(normal_panc, label = F, group.by = "orig.ident")

saveRDS(normal_panc, file = "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/Nataly-Objects/normal_panc.rds")

##############################################################
######   For 24 PDAC Samples from GSA:PRJCA001063   ##########
##############################################################

# Print the list of Seurat objects
print(scRNAlist)

All_PDAC <- Reduce(function(x, y) merge(x, y), scRNAlist)

#saveRDS(All_PDAC, file = "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/Nataly-Objects/All_PDAC.rds")

######### Normalize data and perform harmony batch effect correction
minGene=200
maxUMI=25000
pctMT=10

VlnPlot(All_PDAC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
All_PDAC <- subset(All_PDAC, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & percent.mt < pctMT)
VlnPlot(All_PDAC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

All_PDAC <- NormalizeData(All_PDAC) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
All_PDAC <- RunPCA(All_PDAC, verbose = F)
ElbowPlot(All_PDAC, ndims = 50)
pc.num=1:30
All_PDAC <- All_PDAC %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
All_PDAC <- FindNeighbors(All_PDAC, dims=pc.num)
All_PDAC <- FindClusters(All_PDAC)
All_PDAC <- RunUMAP(All_PDAC, dims = 1:30)
DimPlot(All_PDAC, label = T)
DimPlot(All_PDAC, label = T, group.by = "orig.ident")

All_PDAC <- RunPCA(All_PDAC, npcs=50, verbose=FALSE)
All_PDAC <- RunHarmony(All_PDAC, group.by.vars="orig.ident", max.iter.harmony = 20) 
All_PDAC <- RunUMAP(All_PDAC, reduction = "harmony", dims = 1:30)
All_PDAC <- FindNeighbors(All_PDAC, reduction = "harmony", dims = 1:30) 
All_PDAC <- FindClusters(All_PDAC, verbose = FALSE, resolution = 0.4)
DimPlot(All_PDAC, label = T)
DimPlot(All_PDAC, label = F, group.by = "orig.ident")

saveRDS(All_PDAC, file = "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/Nataly-Objects/All_PDAC.rds")

#########################################################
######  For 3 NAT Samples from GEO:GSE155698   ##########
#########################################################

NAT_pan <- Reduce(function(x, y) merge(x, y), scRNAlist)

#saveRDS(NAT_pan, file = "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/Nataly-Objects/NAT_pan.rds")

######### Normalize data and perform harmony batch effect correction

minGene=200
maxUMI=25000
pctMT=10

VlnPlot(NAT_pan, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
NAT_pan <- subset(NAT_pan, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & percent.mt < pctMT)
VlnPlot(NAT_pan, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

NAT_pan <- NormalizeData(NAT_pan) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
NAT_pan <- RunPCA(NAT_pan, verbose = F)
ElbowPlot(NAT_pan, ndims = 50)
pc.num=1:30
NAT_pan <- NAT_pan %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
NAT_pan <- FindNeighbors(NAT_pan, dims=pc.num)
NAT_pan <- FindClusters(NAT_pan)
NAT_pan <- RunUMAP(NAT_pan, dims = 1:30)
DimPlot(NAT_pan, label = T)
DimPlot(NAT_pan, label = T, group.by = "orig.ident")

NAT_pan <- RunPCA(NAT_pan, npcs=50, verbose=FALSE)
NAT_pan <- RunHarmony(NAT_pan, group.by.vars="orig.ident", max.iter.harmony = 20) 
NAT_pan <- RunUMAP(NAT_pan, reduction = "harmony", dims = 1:30)
NAT_pan <- FindNeighbors(NAT_pan, reduction = "harmony", dims = 1:30) 
NAT_pan <- FindClusters(NAT_pan, verbose = FALSE, resolution = 0.1)
DimPlot(NAT_pan, label = T)
DimPlot(NAT_pan, label = F, group.by = "orig.ident")

saveRDS(NAT_pan, file = "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/Nataly-Objects/NAT_pan.rds")

#############

saveRDS(normal_panc, file = "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/Nataly-Objects/normal_panc.rds")
saveRDS(All_PDAC, file = "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/Nataly-Objects/All_PDAC.rds")
saveRDS(NAT_pan, file = "/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/Nataly-Objects/NAT_pan.rds")

normal_panc <-readRDS("/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/Nataly-Objects/normal_panc.rds")
All_PDAC<-readRDS("/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/Nataly-Objects/All_PDAC.rds")
NAT_pan <-readRDS("/Volumes/NatalyDrive/Hend_and_Ali_analysis/Hend/Nataly-Objects/NAT_pan.rds")

####################################################
#################### Cell Typing ################### 
####################################################

######################################
############ 24 PDAC samples #########
######################################

All_PDAC <- FindClusters(All_PDAC, verbose = FALSE, resolution = 0.6)
DimPlot(All_PDAC, label = T)

### ductal cell 1 (cluster 6), with res 0.6 (cluster 5)
cluster10Marker=c("AMBP", "CFTR", "MMP7")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### ductal cell 2 (cluster 0,8?), with res 0.6 (cluster 0,10,16)
cluster10Marker=c("KRT19", "KRT7", "TSPAN8", "SLPI")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Acinar (cluster 11 ), with res 0.6 (cluster 17)
cluster10Marker=c("PRSS1", "CTRB1", "CTRB2", "REG1B")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### endocrine cell (cluster 12 ), with res 0.6 (cluster 18)
cluster10Marker=c("CHGB", "CHGA", "INS", "IAPP")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### stellate cell (cluster 3,14 ), with res 0.6 (cluster 4,7,20,23)
cluster10Marker=c("RGS5", "ACTA2", "PDGFRB", "ADIRF")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### fibroblast (cluster 9, 1 ), with res 0.6 (cluster 3,6,14)
cluster10Marker=c("LUM", "DCN", "COL1A1")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

###Endothelia (cluster 2), with res 0.6 (cluster 2,12)
cluster10Marker=c("CDH5", "PLVAP", "VWF", "CLDN5")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### macrophage (cluster 4,14) , with res 0.6 (cluster 1, 20, 11?)
cluster10Marker=c("AIF1", "FCGR1A", "CD14", "CD68", "CSF1R")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),
# 
# ### macrophage (cluster 4,14)
# cluster10Marker=c("ITGAM", "CD68", "HLA-DRA", "ITGAX") ## should be positive all and -ve ITGAX
# DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
#         col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### T cell (cluster 5) , with res 0.6 (cluster 9,13)
cluster10Marker=c("CD3D", "CD3E", "CD4", "CD8A", "CD8B")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### B cell (cluster 7,16,13, 10? ), with res 0.6 (cluster 8,11?,15?,19,22)
cluster10Marker=c("MS4A1", "CD79A", "CD79B", "CD52")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Epithelial cell (cluster 0,6,15 ), with res 0.6 (cluster 21)
cluster10Marker=c("EPCAM", "CDH1", "CLDN3","CLDN4")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### B cell and macrophage marker  (cluster 8? ), with res 0.6 (cluster 11?)
cluster10Marker=c("MS4A1", "CD79A", "CD79B", "CD52", "ITGAM", "CD68", "CSF1R")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Monocytes  (cluster  )
cluster10Marker=c("ITGAM", "CD14", "HLA-DRA", "MRC1", "CD86")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Plasma cell  (cluster  )
cluster10Marker=c("TNFRSF17", "SDC1")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Dendretic cell  (cluster  )
cluster10Marker=c("ITGAX", "HLA-DRA", "ITGAM")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### NK Tcell  (cluster  )
cluster10Marker=c("NCAM1", "CD3D", "CD3G", "CD3E")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Neutrophiles  (cluster  )
cluster10Marker=c("FCGR3A", "CEACAM8", "FUT4")
DotPlot(object = All_PDAC, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### regulation of cell plasticity  (cluster  )
cluster10Marker=c("PAX5", "PRDM1", "SPI1", "CXCR4", "CXCR2", "FOXO1")
DotPlot(object = All_PDAC, features = cluster10Marker, group.by="cell_type",dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

# viridis_plasma_dark_high
# viridis_plasma_light_high
# viridis_magma_dark_high
# viridis_magma_light_high
# viridis_inferno_dark_high
# viridis_inferno_light_high
# viridis_dark_high
# viridis_light_high

####### Add cell types to metadata
# Identify clusters and create a named vector mapping cluster IDs to cell types
cluster_to_cell_type_list <- list(
  "Ductal cell 1" = c("5"),
  "Ductal cell 2" = c("0", "10", "16"),
  "Acinar cell"= c("17"),
  "Endocrine cell"= c("18"),
  "Stellate cell"= c("4","7", "20","23"),
  "Fibroblast"= c("3", "6", "13"),
  "Endothelial cell cell"= c("2", "12"),
  "Macrophage"= c("1", "20"),
  "T cell"= c("9", "13"),
  "B cell"= c("8", "15", "19","22"),
  "Epithelial cell"= c("21"),
  "B cell with macrophage marker"= c("11"))

# Initialize a new metadata column for cell types
All_PDAC$cell_type <- "Unknown"

# Map the clusters to cell types
for (cell_type in names(cluster_to_cell_type_list)) {
  seurat_clusters <- cluster_to_cell_type_list[[cell_type]]
  All_PDAC$cell_type[Idents(All_PDAC) %in% seurat_clusters] <- cell_type
}

DimPlot(All_PDAC, reduction = "umap", label = T, group.by = "cell_type", repel = T) 


FeaturePlot_scCustom(seurat_object = All_PDAC, features = c("PAX5", "PRDM1", "SPI1", "CXCR4", "CXCR2", "FOXO1", "CEBPA", "CEBPB"), reduction = 'umap',
                     colors_use = viridis_plasma_light_high, na_color = "lightgray",na_cutoff=0, raster=FALSE, num_columns =2)
############################################################################################
###### To check the different Bcell states from C. Domínguez Conde et al. ##################
############################################################################################
cluster10Marker=c("CD79A", "SOX4", "IGLL1", "IGHD", "IGHM", "TCL1A", "CD27", "CR2", "MKI67", "POU2AF1", "SUGCT", "JCHAIN", "XBP1")
DotPlot(object = All_PDAC, features = cluster10Marker, group.by="cell_type",dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

####### To subset
All_PDAC_Bcell <- subset(All_PDAC, subset = cell_type %in% c("B cell with macrophage marker", "B cell"))

cluster10Marker=c("CD79A", "SOX4", "IGLL1", "IGHD", "IGHM", "TCL1A", "CD27", "CR2", "MKI67", "POU2AF1", "SUGCT", "JCHAIN", "XBP1")
DotPlot(object = All_PDAC_Bcell, features = cluster10Marker, dot.scale = 10, col.min = -2, col.max = 2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") + coord_flip()

### Immature Bcell  (cluster  )
cluster10Marker=c("CD79A", "SOX4", "IGLL1", "IGHD", "IGHM", "TCL1A")
DotPlot(object = All_PDAC_Bcell, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Naïve Bcell  (cluster 8 )
cluster10Marker=c( "CD79A", "IGHD", "IGHM", "TCL1A")
DotPlot(object = All_PDAC_Bcell, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Memory Bcell  (cluster 15  )
cluster10Marker=c("CD79A", "CD27", "CR2")
DotPlot(object = All_PDAC_Bcell, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Cycling Bcell  (cluster 11)
cluster10Marker=c("TCL1A", "MKI67", "POU2AF1")
DotPlot(object = All_PDAC_Bcell, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Germinal center Bcell  (cluster 15  )
cluster10Marker=c("CD79A", "TCL1A", "CD27", "CR2", "POU2AF1", "SUGCT")
DotPlot(object = All_PDAC_Bcell, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

All_PDAC_Bcell <- FindClusters(All_PDAC_Bcell, verbose = FALSE, resolution = 0.4)
DimPlot(All_PDAC_Bcell, reduction = "umap", label = T, repel = T) 

######################################
######### Normal 11 samples ##########
######################################
#normal_panc <- FindClusters(normal_panc, verbose = FALSE, resolution = 0.6)
DimPlot(normal_panc, label = T)

### with res = 0.4
### ductal cell 1 (cluster 0, 2,9 )
cluster10Marker=c("AMBP", "CFTR", "MMP7")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### ductal cell 2 (cluster 3, 12 )
cluster10Marker=c("KRT19", "KRT7", "TSPAN8", "SLPI")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Acinar (cluster  5,8)
cluster10Marker=c("PRSS1", "CTRB1", "CTRB2", "REG1B")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### endocrine cell (cluster 13 )
cluster10Marker=c("CHGB", "CHGA", "INS", "IAPP")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### stellate cell (cluster 6 )
cluster10Marker=c("RGS5", "ACTA2", "PDGFRB", "ADIRF")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### fibroblast (cluster 4)
cluster10Marker=c("LUM", "DCN", "COL1A1")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

###Endothelial (cluster 11, 1 )
cluster10Marker=c("CDH5", "PLVAP", "VWF", "CLDN5")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### macrophage (cluster 7 ) 
cluster10Marker=c("AIF1", "FCGR1A", "CD14", "CD68", "CSF1R")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### T cell (cluster ) 
cluster10Marker=c("CD3D", "CD3E", "CD4", "CD8A", "CD8B")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### B cell (cluster  )
cluster10Marker=c("MS4A1", "CD79A", "CD79B", "CD52")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Epithelial cell (cluster 2,3,12 )
cluster10Marker=c("EPCAM", "CDH1", "CLDN3","CLDN4")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### B cell and macrophage marker  (cluster  )
cluster10Marker=c("MS4A1", "CD79A", "CD79B", "CD52", "ITGAM", "CD68", "CSF1R")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Monocytes  (cluster  )
cluster10Marker=c("ITGAM", "CD14", "HLA-DRA", "MRC1", "CD86")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Plasma cell  (cluster  )
cluster10Marker=c("TNFRSF17", "SDC1")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Dendretic cell  (cluster  )
cluster10Marker=c("ITGAX", "HLA-DRA", "ITGAM")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### NK Tcell  (cluster  )
cluster10Marker=c("NCAM1", "CD3D", "CD3G", "CD3E")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Neutrophiles  (cluster  )
cluster10Marker=c("FCGR3A", "CEACAM8", "FUT4")
DotPlot(object = normal_panc, features = cluster10Marker,dot.scale = 10,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

####### Add cell types to metadata
# Identify clusters and create a named vector mapping cluster IDs to cell types
cluster_to_cell_type_list <- list(
  "Ductal cell 1" = c("0", "2", "9"),
  "Ductal cell 2" = c("3", "12"),
  "Acinar cell"= c("5", "8"),
  "Endocrine cell"= c("13"),
  "Stellate cell"= c("6"),
  "Fibroblast"= c("4"),
  "Endothelial cell cell"= c("1", "11"),
  "Macrophage"= c("7"))

# Initialize a new metadata column for cell types
normal_panc$cell_type <- "Unknown"

# Map the clusters to cell types
for (cell_type in names(cluster_to_cell_type_list)) {
  seurat_clusters <- cluster_to_cell_type_list[[cell_type]]
  normal_panc$cell_type[Idents(normal_panc) %in% seurat_clusters] <- cell_type
}

DimPlot(normal_panc, reduction = "umap", label = T, group.by = "cell_type", repel = T) 

######################################
########## NAT 3 samples #############
######################################
#NAT_pan <- FindClusters(NAT_pan, verbose = FALSE, resolution = 0.6)
DimPlot(NAT_pan, label = T)

### with res = 0.4
### ductal cell 1 (cluster 6)
cluster10Marker=c("AMBP", "CFTR", "MMP7")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### ductal cell 2 (cluster 6 )
cluster10Marker=c("KRT19", "KRT7", "TSPAN8", "SLPI")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Acinar (cluster  1,12)
cluster10Marker=c("PRSS1", "CTRB1", "CTRB2", "REG1B")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### endocrine cell (cluster  )
cluster10Marker=c("CHGB", "CHGA", "INS", "IAPP")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### stellate cell (cluster 7 )
cluster10Marker=c("RGS5", "ACTA2", "PDGFRB", "ADIRF")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### fibroblast (cluster 14)
cluster10Marker=c("LUM", "DCN", "COL1A1")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

###Endothelial (cluster 10 )
cluster10Marker=c("CDH5", "PLVAP", "VWF", "CLDN5")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### macrophage (cluster 3) 
cluster10Marker=c("AIF1", "FCGR1A", "CD14", "CD68", "ITGAX", "CSF1R")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### T cell (cluster 0 CD8 Tcell (part of 9) 
cluster10Marker=c("CD3D", "CD3E", "CD4", "CD8A", "CD8B")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### B cell (cluster 13 )
cluster10Marker=c("MS4A1", "CD79A", "CD79B", "CD52")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Epithelial cell (cluster 15? )
cluster10Marker=c("EPCAM", "CDH1", "CLDN3","CLDN4")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### B cell and macrophage marker  (cluster  )
cluster10Marker=c("MS4A1", "CD79A", "CD79B", "CD52", "ITGAM", "CD68", "CSF1R")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Monocytes  (cluster 2,5 )
cluster10Marker=c("ITGAM", "CD14", "HLA-DRA", "MRC1", "CD86")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Plasma cell  (cluster )
cluster10Marker=c("TNFRSF17", "SDC1")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Dendretic cell  (cluster  3)
cluster10Marker=c("ITGAX", "HLA-DRA", "ITGAM")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### NK Tcell  (cluster part of 9 )
cluster10Marker=c("NCAM1", "CD3D", "CD3G", "CD3E")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

### Neutrophiles  (cluster  )
cluster10Marker=c("FCGR3A", "CEACAM8", "FUT4")
DotPlot(object = NAT_pan, features = cluster10Marker,dot.scale = 10,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

####### Add cell types to metadata
# Identify clusters and create a named vector mapping cluster IDs to cell types
cluster_to_cell_type_list <- list(
  "Ductal cell" = c("6"),
  "Acinar cell"= c("1","12"),
  "Stellate cell"= c("7"),
  "Fibroblast"= c("14"),
  "Endothelial cell cell"= c("10"),
  "Monocyte"= c("2", "5"),
  "T cell"= c("0", "9"),
  "B cell"= c("13"),
  "Epithelial cell"= c("15"),
  "Macrophage"= c("3"),
  "NK Tcell"=c("9"))

# Initialize a new metadata column for cell types
NAT_pan$cell_type <- "Unknown"

# Map the clusters to cell types
for (cell_type in names(cluster_to_cell_type_list)) {
  seurat_clusters <- cluster_to_cell_type_list[[cell_type]]
  NAT_pan$cell_type[Idents(NAT_pan) %in% seurat_clusters] <- cell_type
}

DimPlot(NAT_pan, reduction = "umap", label = T, group.by = "cell_type", repel = T) 

