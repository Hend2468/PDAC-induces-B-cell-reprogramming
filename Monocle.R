############################# 
######### Monocle ###########
############################# 

##### https://cole-trapnell-lab.github.io/monocle-release/docs/#introduction

####################################################################################### 
########################### scRNA-seq Creating a Seurat Object ########################
####################################################################################### 

library(devtools)
library(Seurat)
library(SeuratData)
library(DESeq2)
library(dplyr)
library(dittoSeq)
library(tidyverse)
library(patchwork)
library(grid)
library(viridis)
library(ggplot2)
library(remotes)
library(future)
library(RColorBrewer)
library(ggridges)
library(data.table)
library(enrichR)
library(scCustomize)
theme_set(theme_bw(16))
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))


NAT1 <- Read10X(data.dir = "/Volumes/NatalyDrive/Hend\ and\ Ali\ analysis/Hend/normal\ pan\ scRNA/AdjNorm_TISSUE_1/filtered_feature_bc_matrix/")
NAT2 <- Read10X(data.dir = "/Volumes/NatalyDrive/Hend\ and\ Ali\ analysis/Hend/normal\ pan\ scRNA/AdjNorm_TISSUE_2/filtered_feature_bc_matrix/")
NAT3 <- Read10X(data.dir = "/Volumes/NatalyDrive/Hend\ and\ Ali\ analysis/Hend/normal\ pan\ scRNA/AdjNorm_TISSUE_3/filtered_feature_bc_matrix/")

NAT1 <- CreateSeuratObject(counts = NAT1, project = "NAT1", min.cells = 3, min.features = 200)
NAT2 <- CreateSeuratObject(counts = NAT2, project = "NAT2", min.cells = 3, min.features = 200)
NAT3 <- CreateSeuratObject(counts = NAT3, project = "NAT3", min.cells = 3, min.features = 200)

PDAC1 <- Read10X(data.dir = "/Volumes/NatalyDrive/Hend\ and\ Ali\ analysis/Hend/PDAC\ scRNA/PDAC/T1/")
PDAC1 <- CreateSeuratObject(counts = PDAC1, project = "PDAC1", min.cells = 3, min.features = 200)


###### NAT1 ######
NAT1[["percent.mt"]] <- PercentageFeatureSet(NAT1, pattern = "^MT-")
NAT1[["percent.rb"]] <- PercentageFeatureSet(NAT1, pattern = "^RP[SL]")

HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(NAT1))
NAT1[["percent.HB"]]<-PercentageFeatureSet(NAT1, features=HB.genes) 

VlnPlot(NAT1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(NAT1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NAT1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

minGene=1200
maxGene=5000
maxUMI=15000
pctMT=10
pctHB=1

NAT1 <- subset(NAT1, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB < pctHB)

VlnPlot(NAT1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
NAT1

######### NAT2 #######
NAT2[["percent.mt"]] <- PercentageFeatureSet(NAT2, pattern = "^MT-")
NAT2[["percent.rb"]] <- PercentageFeatureSet(NAT2, pattern = "^RP[SL]")

HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(NAT2))
NAT2[["percent.HB"]]<-PercentageFeatureSet(NAT2, features=HB.genes) 

VlnPlot(NAT2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(NAT2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NAT2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

minGene=1200
maxGene=5000
maxUMI=15000
pctMT=10
pctHB=1

NAT2 <- subset(NAT2, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & 
                 nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB < pctHB)

VlnPlot(NAT2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
NAT2

######### NAT3 #######
NAT3[["percent.mt"]] <- PercentageFeatureSet(NAT3, pattern = "^MT-")
NAT3[["percent.rb"]] <- PercentageFeatureSet(NAT3, pattern = "^RP[SL]")

HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(NAT3))
NAT3[["percent.HB"]]<-PercentageFeatureSet(NAT3, features=HB.genes) 

VlnPlot(NAT3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(NAT3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NAT3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

minGene=1200
maxGene=5000
maxUMI=15000
pctMT=10
pctHB=1

NAT3 <- subset(NAT3, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & 
                 nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB < pctHB)

VlnPlot(NAT3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
NAT3
############################################################
############ Merge Object and Perform harmony ##############
############################################################

NAT1 <- RenameCells(object = NAT1, add.cell.id = "NAT1")
head(x = colnames(x = NAT1))
NAT2 <- RenameCells(object = NAT2, add.cell.id = "NAT2")
head(x = colnames(x = NAT2))
NAT3 <- RenameCells(object = NAT3, add.cell.id = "NAT3")
head(x = colnames(x = NAT3))

Normal_NAT <- merge(x=NAT1, y =c(NAT2, NAT3),  project = "Normal_NAT")

######### Normalize data and perform harmony batch effect correction
Normal_NAT <- NormalizeData(Normal_NAT) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
Normal_NAT <- RunPCA(Normal_NAT, verbose = F)
ElbowPlot(Normal_NAT, ndims = 50)
pc.num=1:25
Normal_NAT <- Normal_NAT %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
Normal_NAT <- FindNeighbors(Normal_NAT, dims=pc.num)
Normal_NAT <- FindClusters(Normal_NAT)
Normal_NAT <- RunUMAP(Normal_NAT, dims = 1:25)
DimPlot(Normal_NAT, label = T)

Normal_NAT <- RunPCA(Normal_NAT, npcs=50, verbose=FALSE)
Normal_NAT <- RunHarmony(Normal_NAT, group.by.vars="orig.ident", max.iter.harmony = 20) 
Normal_NAT <- RunUMAP(Normal_NAT, reduction = "harmony", dims = 1:25)
Normal_NAT <- FindNeighbors(Normal_NAT, reduction = "harmony", dims = 1:25) 
Normal_NAT <- FindClusters(Normal_NAT, verbose = FALSE, resolution = 0.2)
DimPlot(Normal_NAT, label = T)
DimPlot(Normal_NAT, label = T, split.by = "orig.ident")

top10 <- head(VariableFeatures(Normal_NAT), 10)

saveRDS(Normal_NAT, file = "/Volumes/NatalyDrive/Hend\ and\ Ali\ analysis/Hend/Nataly-Objects/Normal_NAT.rds")


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Normal_NAT)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(NAT1)
NAT1 <- ScaleData(NAT1, features = all.genes)
NAT1 <- RunPCA(NAT1, features = VariableFeatures(object = NAT1))
# Examine and visualize PCA results a few different ways
VizDimLoadings(NAT1, dims = 1:2, reduction = "pca")
ElbowPlot(NAT1)


# Fibroblast: COL1A1, COL3A1, FAP, THY1, PDGFRB, VIM, ACTA2, CDH11, ENG; 
# Epithelial cell: EPCAM, CDH1, CLDN3,CLDN4; 
# B cell: PTPRC, CD19, CD79A, CD79B; 
# T cell: PTPRC, CD4, CTLA4, FOXP3; 
# Vascular endothelial cell: VIM, EMCN,CLEC14A, CDH5, VWF, CD34, PECAM1; 
# Macrophage/monocyte: CD14, CD68, CSF1R, FGFR3, CD33, LYZ; Pericyte: MCAM,CSPG4; 
# Lymphatic endothelial cell: LYVE1, PROX1



########## Clustering #######
NAT1 <- FindNeighbors(NAT1, dims = 1:10)
NAT1 <- FindClusters(NAT1, resolution = 0.5)
NAT1 <- RunUMAP(NAT1, dims = 1:10)
DimPlot(NAT1, reduction = "umap")

############################################
########### Remove Doublets ################
############################################

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("scDblFinder")
# suppressMessages(require(DoubletFinder))
# remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

library(DoubletFinder)
library(scDblFinder)

# define the expected number of doublet cellscells.
nExp <- round(ncol(NAT1) * 0.04)  # expect 4% doublets
NAT1 <- doubletFinder(NAT1, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(NAT1@meta.data)[grepl("DF.classification", colnames(NAT1@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(NAT1, group.by = "orig.ident") + NoAxes(),
                   DimPlot(NAT1, group.by = DF.name) + NoAxes())

VlnPlot(NAT1, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
##Remove all predicted doublets from our data.
NAT1 = NAT1[, NAT1@meta.data[, DF.name] == "Singlet"]

############################################
######### Subset to neg cells only #########
############################################
FeaturePlot(NAT1,"NAT1")

NAT1_expression = GetAssayData(object = NAT1, 
                               assay = "RNA", slot = "data")["NAT1",]
neg_ids = names(which(NAT1_expression==0))
NAT1 = subset(NAT1,cells=neg_ids)
FeaturePlot(NAT1,"NAT1")

############################################
############## Cell Typing #################
############################################

NAT1 <- AddModuleScore(object = NAT1, features = c(Melanoma), name = c("Melanoma"))
NAT1 <- AddModuleScore(object = NAT1, features = c(ColdSig, HotSig), name = c("ColdSig", "HotSig"))

NAT1 <- AddModuleScore(object = NAT1, features = c(IFNG_antigen_pres, IFNG_feedback_sign, IFNG_Chemoattractants, IFNG_Cytotoxic_eff), name = c("IFNG_antigen_pres", "IFNG_feedback_sign", "IFNG_Chemoattractants", "IFNG_Cytotoxic_eff"))

NAT1 <- AddModuleScore(object = NAT1, features = c(CD8_T_Effector_Memory1, CD8_T_Effector_Memory_2, CD8_T_Activated, CD8_T_Progenitor_Exhausted,CD8_T_Proliferating,
                                                   CD8_T_Apoptotic, CD8_NK_like_T, CD8_T_Terminally_Exhausted1, CD8_Treg_like_T, CD8_γδ_like_T, CD8_T_Terminally_Exhausted2, CD8_T_Naive),
                       name = c("CD8_T_Effector_Memory1", "CD8_T_Effector_Memory_2", "CD8_T_Activated", "CD8_T_Progenitor_Exhausted","CD8_T_Proliferating",
                                "CD8_T_Apoptotic", "CD8_NK_like_T", "CD8_T_Terminally_Exhausted1", "CD8_Treg_like_T", "CD8_γδ_like_T", "CD8_T_Terminally_Exhausted2", "CD8_T_Naive"))


NAT1 <- AddModuleScore(object = NAT1, features = c(Naive_CD8_c3_Tn, Naive_CD8_c13_Tn_TCF7, Transit_Effector_CD8_c0_t_Teff, Effector_CD8_c2_Teff, Effector_CD8_c8_Teff_KLRG1,Effector_CD8_c10_Teff_CD244,
                                                   Effector_CD8_c11_Teff_SEMA4A, Central_Memory_CD8_c6_Tcm, Resident_Memory_CD8_c12_Trm, Precursor_Exhausted_CD8_c7_p_Tex, Exhausted_CD8_c1_Tex,
                                                   Stress_Response_CD8_c4_Tstr, IFN_Response_CD8_c5_Tisg, Senescent_like_CD8_c9_Tsen), 
                       name = c("Naive_CD8_c3_Tn", "Naive_CD8_c13_Tn_TCF7", "Transit_Effector_CD8_c0_t_Teff", "Effector_CD8_c2_Teff", "Effector_CD8_c8_Teff_KLRG1","Effector_CD8_c10_Teff_CD244",
                                "Effector_CD8_c11_Teff_SEMA4A", "Central_Memory_CD8_c6_Tcm", "Resident_Memory_CD8_c12_Trm", "Precursor_Exhausted_CD8_c7_p_Tex", "Exhausted_CD8_c1_Tex",
                                "Stress_Response_CD8_c4_Tstr", "IFN_Response_CD8_c5_Tisg", "Senescent_like_CD8_c9_Tsen"))


colnames(NAT1@meta.data)[9:41]
colnames(NAT1@meta.data)[9:41] <- c("Melanoma", "ColdSig", "HotSig", "IFNG_antigen_pres", "IFNG_feedback_sign", "IFNG_Chemoattractants", "IFNG_Cytotoxic_eff",
                                    "CD8_T_Effector_Memory1", "CD8_T_Effector_Memory_2", "CD8_T_Activated", "CD8_T_Progenitor_Exhausted","CD8_T_Proliferating",
                                    "CD8_T_Apoptotic", "CD8_NK_like_T", "CD8_T_Terminally_Exhausted1", "CD8_Treg_like_T", "CD8_γδ_like_T", "CD8_T_Terminally_Exhausted2", "CD8_T_Naive",
                                    "Naive_CD8_c3_Tn", "Naive_CD8_c13_Tn_TCF7", "Transit_Effector_CD8_c0_t_Teff", "Effector_CD8_c2_Teff", "Effector_CD8_c8_Teff_KLRG1","Effector_CD8_c10_Teff_CD244",
                                    "Effector_CD8_c11_Teff_SEMA4A", "Central_Memory_CD8_c6_Tcm", "Resident_Memory_CD8_c12_Trm", "Precursor_Exhausted_CD8_c7_p_Tex", "Exhausted_CD8_c1_Tex",
                                    "Stress_Response_CD8_c4_Tstr", "IFN_Response_CD8_c5_Tisg", "Senescent_like_CD8_c9_Tsen")

############# "ColdSig" vs "HotSig"

data.sig <- NAT1@meta.data[, c("ColdSig", "HotSig")]
data.sig <- data.sig %>% mutate(colnames(data.sig)[apply(data.sig,1,which.max)])
names(data.sig)[3]<-paste("sig.state")
NAT1 <- AddMetaData(NAT1, data.sig$sig.state, col.name = "sig.state")

DimPlot(NAT1, group.by	= c("sig.state"))

############# Melanoma vs non-melanoma

data.mel <- NAT1@meta.data[, c("Melanoma", "sig.state")]

data.mel[data.mel >0] <- "melanoma"
data.mel[data.mel ==0] <- "non_melanoma"
data.mel[data.mel <0] <- "non_melanoma"
NAT1 <- AddMetaData(NAT1, data.mel$Melanoma, col.name = "melanoma_expression")

DimPlot(NAT1, group.by	= c("melanoma_expression"))

############ IFNG signature 

data.sig <- NAT1@meta.data[, c("IFNG_antigen_pres", "IFNG_feedback_sign", "IFNG_Chemoattractants", "IFNG_Cytotoxic_eff")]
data.sig <- data.sig %>% mutate(colnames(data.sig)[apply(data.sig,1,which.max)])
names(data.sig)[5]<-paste("sig.state")
NAT1 <- AddMetaData(NAT1, data.sig$sig.state, col.name = "IFNG.state")

DimPlot(NAT1, group.by	= c("IFNG.state"))

########### CD8clusters1 Oliviera et al., Nature 2021

data.sig <- NAT1@meta.data[, c("CD8_T_Effector_Memory1", "CD8_T_Effector_Memory_2", "CD8_T_Activated", "CD8_T_Progenitor_Exhausted","CD8_T_Proliferating",
                               "CD8_T_Apoptotic", "CD8_NK_like_T", "CD8_T_Terminally_Exhausted1", "CD8_Treg_like_T", "CD8_γδ_like_T", "CD8_T_Terminally_Exhausted2", "CD8_T_Naive")]
data.sig <- data.sig %>% mutate(colnames(data.sig)[apply(data.sig,1,which.max)])
names(data.sig)[13]<-paste("sig.state")
NAT1 <- AddMetaData(NAT1, data.sig$sig.state, col.name = "CD8clusters1")

DimPlot(NAT1, group.by	= c("CD8clusters1"))

########## CD8clusters2 Chu et al., Nat Med 2023

data.sig <- NAT1@meta.data[, c("Naive_CD8_c3_Tn", "Naive_CD8_c13_Tn_TCF7", "Transit_Effector_CD8_c0_t_Teff", "Effector_CD8_c2_Teff", "Effector_CD8_c8_Teff_KLRG1","Effector_CD8_c10_Teff_CD244",
                               "Effector_CD8_c11_Teff_SEMA4A", "Central_Memory_CD8_c6_Tcm", "Resident_Memory_CD8_c12_Trm", "Precursor_Exhausted_CD8_c7_p_Tex", "Exhausted_CD8_c1_Tex",
                               "Stress_Response_CD8_c4_Tstr", "IFN_Response_CD8_c5_Tisg", "Senescent_like_CD8_c9_Tsen")]
data.sig <- data.sig %>% mutate(colnames(data.sig)[apply(data.sig,1,which.max)])
names(data.sig)[15]<-paste("sig.state")
NAT1 <- AddMetaData(NAT1, data.sig$sig.state, col.name = "CD8clusters2")

DimPlot(NAT1, group.by	= c("CD8clusters2"))

FeaturePlot(NAT1, features = c("CD8A", "CD8B"), blend=TRUE)


####################################

saveRDS(NAT1, file = "/Volumes/NatalyDrive/Justin_scRNA/objects/NAT1.rds")
vhl <- readRDS("/Volumes/NatalyDrive/Justin_scRNA/objects/vhl.rds")

