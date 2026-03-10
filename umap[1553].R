library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(hdf5r)
rm(list = ls())
## ??��??ȡ????
### ????????·????????????
setwd("D:/Analysis/scRNA/PDAC3")#ԭʼ???ݴ?????????һĿ¼
assays <- dir("GSE196678/")
dir <- paste0("GSE196678/", assays)
# ???ļ?˳?????????????????Ʋ?Ҫ?????ֿ?ͷ???м䲻???пո?
samples_name = c('L15')
samples_name = c('OS (1)', 'OS (2)', 'OS (3)', 'OS (4)', 'OS (5)', 'OS (6)')
samples_name = c('LiM', 'LuM', 'VM')
samples_name = c('BC2', 'BC3', 'BC5', 'BC6', 'BC10', 'BC11', 'BC16', 'BC17', 'BC20', 'BC21', 'BC22')
samples_name = c('PDAC (1)', 'PDAC (2)', 'PDAC (3)',
                 'PDAC (4)', 'PDAC (5)', 'PDAC (6)',
                 'PDAC (7)', 'PDAC (8)', 'PDAC (9)',
                 'PDAC (10)','PDAC (11)', 'PDAC (12)', 'PDAC (13)',
                 'PDAC (14)', 'PDAC (15)', 'PDAC (16)',
                 'PDAC (17)', 'PDAC (18)', 'PDAC (19)',
                 'PDAC (20)','PDAC (21)', 'PDAC (22)', 'PDAC (23)',
                 'PDAC (24)')

#??ȡ???ݵڶ??ַ???
scRNA <- Read10X(data.dir = "D:/Analysis/scRNA/LRRC15/L15")
####################################??��????seurat????###
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  #??????min.cells???˻????ᵼ??CellCycleScoring???���
  #Insufficient data values to produce 24 bins.  
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i],
                                       min.cells=3, min.features = 200)
  #??ϸ??barcode?Ӹ?ǰ׺????ֹ?ϲ???barcode????
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])   
  #??????��??????????
  if(T){    
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
  }
  #??????????????????
  if(T){
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")
  }
  #??????ϸ??????????
  if(T){
    HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
    scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]], features=HB.genes) 
  }
}
### ???б???????????????
setwd("E:/osteosarcoma/GSE152048")
names(scRNAlist) <- samples_name
#system.time(save(scRNAlist, file = "Integrate/scRNAlist0.Rdata")) 
system.time(saveRDS(scRNAlist, file = "scRNAlist0.rds"))
####1.3 ʹ??merge??????scRNAlist?ϳ?һ??Seurat????#######
scRNA <- readRDS("scRNAlist0.rds")
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
dim(scRNA)
scRNA
#save(scRNA,file='scRNA.Rdata')
# load("scRNA.Rdata")
# An object of class Seurat 
# 18818 features across 19738 samples within 1 assay 
# Active assay: RNA (18818 features, 0 variable features)
table(scRNA$orig.ident)
#save(scRNA,file = 'scRNA_orig.Rdata')
#scRNAlist <- SplitObject(scRNA, split.by = "orig.ident") #?ָ?Seurat????
# # # # # # # # # # # # # ?????ʿ?# # # # # # # # # # # # # 
### ?????ʿ?С????ͼ
# ???ÿ????õ???????
theme.set2 = theme(axis.title.x=element_blank())
# ???û?ͼԪ??
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.HB")
group = "orig.ident"
# ?ʿ?ǰС????ͼ
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)    
dir.create("QC")
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 9, height = 8)
### ?????ʿر?׼
minGene=1200#ԭ??Ϊ500
maxGene=5000
maxUMI=15000
pctMT=10
pctHB=1

### ?????ʿز?????С????ͼ
scRNA <- subset(scRNA, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB < pctHB)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)     
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 10, height = 8) 
#########3. ?鿴????ЧӦ????merge????Seurat???????б?׼???ͽ?ά??###########
# ??ά????
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
scRNA <- RunPCA(scRNA, verbose = F)
ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
scRNA <- scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
scRNA <- FindNeighbors(scRNA, dims=pc.num) %>% FindClusters()
DimPlot(scRNA, label = T)
#?鿴??ͬ?????仹????????ЧӦ?Ĵ???
p <- DimPlot(scRNA, group.by = "orig.ident")
ggsave("UMAP_Samples.pdf", p, width = 8, height = 6)
#???Կ???????ͬ?????仹????????ЧӦ?Ĵ??ڡ?????????ͼ??????????ͼ��??
#??ͼΪÿ???????ľ???ͼ
p <- DimPlot(scRNA, group.by = "orig.ident", split.by = "orig.ident", ncol = 4)
ggsave("UMAP_Samples_Split.pdf", p, width = 18, height = 12)
#saveRDS(scRNA, "scRNA.rds")
#######################harmony????#############################
#rm(list = ls())
scRNA <- readRDS("osteoarthritis.rds")
cellinfo <- subset(scRNA@meta.data, select = c("orig.ident", "percent.mt", "percent.rb", "percent.HB"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
## ??SCT??׼?????ݣ??˷???????NormalizeData?????ǶԵ?????��Ҫ????
memory.limit()
memory.limit(size=156000)
#scRNA <- SCTransform(scRNA)
#?? NormalizeData?????ݽ??б?׼??
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()

###################ʹ??harmony????????#####################
### PCA
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)
### ???Ϸ???1
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony = 20) #max.iter.harmony = 20
### ???Ϸ???2???????????????????ϣ??Ƽ???Ч?????ã?
#scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 
# group.by.vars?????????ð??ĸ?????��????
# max.iter.harmony???õ?????????Ĭ????10??????RunHarmony????????ʾ?ڵ??????ٴκ?????????��??
#??????RunHarmony???????и?lambda??????Ĭ??ֵ??1????????Harmony???ϵ?��?ȡ?lambdaֵ??С??????��?ȱ??󣬷?֮????ֻ??????????Ӱ??????��?ȣ???????Χһ????0.5-2֮?䣩
########################################2.4 ??ά????################################
ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
scRNA <- RunTSNE(scRNA, reduction="harmony", dims=pc.num) %>% RunUMAP(reduction="harmony", dims=pc.num)
#scRNA2 <- RunTSNE(scRNA2, reduction="harmony", dims=pc.num) %>% RunUMAP(reduction="harmony", dims=pc.num)

p <- DimPlot(scRNA, group.by = "orig.ident")
ggsave("monocle  UMAP_Samples_harmony.pdf", p, width = 8, height = 6)
p <- DimPlot(scRNA, group.by = "orig.ident", split.by = "orig.ident", ncol = 4)
ggsave("monocle UMAP_Samples_Split_harmony.pdf", p, width = 18, height = 12)

##save seurat object
saveRDS(scRNA, "monocle SCT_harmony.rds") 
##################################################???Ͻ????��?,??ά????############
setwd("D:/Analysis/scRNA/PDAC3/Epithelial cells")#ԭʼ???ݴ?????????һĿ¼
scRNA <- readRDS("final cancer cell.rds")
##UMAP??????
scRNA <- FindNeighbors(scRNA, dims = 2:40) %>% FindClusters(resolution = 0.004)#resolution????????????
options(repr.plot.height = 4, repr.plot.width = 6)#ԭ??Ϊ4,6
p <-DimPlot(scRNA, reduction = "umap", label = TRUE, pt.size = 0.7)#ԭ??Ϊ0.7
ggsave(" monocle 40-0.05.pdf", p, width = 6.5, height = 5)#ԭ??Ϊwidth = 18, height = 12
mallignant <- read.delim("each cell CNV scores.txt", row.names = 1)#?ļ?Ϊ?ֶ???????СȺ?ϲ?
scRNA <- AddMetaData(scRNA, metadata = mallignant)
p <-DimPlot(scRNA, reduction = "umap", label = TRUE, pt.size = 0.7,group.by = "seurat_clusters",
            split.by = "group")#ԭ??Ϊ0.7

##Ѱ?Ҳ?????????????
scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)#ԭ??Ϊ0.25 

#?ٶ??ϵ???FindAllMarkers,???Ƕ???ϸ??��?Ƚ??ٵ?Ⱥ׼ȷ???Բ?
#remotes::install_github("genecell/COSGR")
library(COSG)
scRNA.markers<- cosg(pbmc,groups='all',assay='RNA',slot='data',mu=1,n_genes_user=100)#ǰ100??????????
#
write.table(scRNA.markers,file="scRNA.markers.xlsx",sep="\t",row.names=F,quote=F)
library(dplyr)
sig.markers=scRNA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
write.table(sig.markers,file="06.markers.xls",sep="\t",row.names=F,quote=F)
top10 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(file="05.DoHeatmap-PDAC.pdf",width=10,height=8)
DoHeatmap(scRNA, features = top10$gene) + NoLegend()#չʾǰ10?????ǻ???????ͼ
dev.off()
pdf(file="05.top20 markers.pdf",width=10,height=8)
VlnPlot(scRNA, features = top10$gene[1:20],pt.size=0)
dev.off()
pdf(file="05.CAF.features-new.pdf",width=8,height=4)

VlnPlot(scRNA, features = c("C0l1a1","Col1a2","Acta2", "Tagln","Igfbp7", "Igfbp", "Igfbp5","Igf1", "Pdgfra","Pdgfb",
                           "Cxcl12", "Il6", "Ccl11","Ccl7","Ccl2","Il1a",
                           "Csf1","Fap", "Ctgf", "Ccl7","Lif", "Cxcl1","Ccr2","Col1a1",
                           "Col6a1","Clec3b","Has1","Col4a1","Cd74","Saa3","H2-Ab1","Sjpi","Fabp4","Car3"),pt.size=0)
VlnPlot(scRNA, features = c("Col1a1","Col1a2","Clec3b","Col14a1","Has1", "Il6","H2-Ab1","Cd74","Saa3","Slpi",
                           "Tagln","Thy1","Col12a1","Thbs2","Col1a1","Col1a2","Pdpn","Dcn","Fap",
                           "Mcam","Cspg4","Pdgfrb","Nes","Rgs5","Igf1"),pt.size=0)
VlnPlot(scRNA, features = c("ACTA2", "TAGLN","IGFBP7", "IGFBP4", "IGFBP5","IGF1", "PDGFRA","PDGFB",
                           "CXCL12", "IL6","IGF1", "CCL11","CCL7","CCL2","CXCR4",
                           "CSF1","FAP", "CTGF", "CCL7","LIF", "CXCL1","CCR2","COL1A1","COL1A2","COL2A2",
                           "COL6A1","IL1A","IL17R","FABP4","CAR3"),pt.size=0)

VlnPlot(scRNA, features = c("COL1A1","COL1A2","FAP","PDPN","DCN","VIM","ACTA2", 
                            "TAGLN","MMP11","MYL9","HOPX","POSTN","TPM1","TPM2",
                           "IL6", "PDGFRA","CXCL12","CFD","DPT","LMNA","AGTR1","HAS1","CXCL1",
                           "CXCL2","CCL2","IL8","IGF1"),pt.size=0)

VlnPlot(scRNA, features = c("COL1A1","COL1A2","FAP","PDPN","DCN","VIM","ACTA2", 
                            "TAGLN","MMP11","MYL9","HOPX","POSTN","TPM1","TPM2",
                            "IL6", "PDGFRA","CXCL12","CFD","DPT","LMNA","AGTR1","HAS1","CXCL1",
                            "CXCL2","CCL2","IL8","SLIT2","LUM"),pt.size=0)
VlnPlot(scRNA, features = c("COL1A1","COL1A2","COL4A1","COL4A2","COL5A2","COL6A3",
                            "SLIT2","LUM"),pt.size=0)
VlnPlot(scRNA, features = c("APOA1","FABP1","DNASE1","CFTR","SOD3","COL18A1"),pt.size=0)
dev.off()
pdf(file="05.iCAF.myCAF-features-new-2.pdf",width=15,height=12)
VlnPlot(scRNA, features = c( "DCN","PDPN","FAP","LUM","PDGFRA","C7","MGP","COL12A1","COL11A1","LOXL2","COL14A1","SDC1","NFIA","IGF1",
                             "SPARCL1","OGN","LGALS1","C5orf46","ANTXR1","C3","ABI3BP","MMP11","SLC16A3","IGFBP4","COL5A2","CTHRC1"),pt.size=0)
dev.off()




pdf(file="05.PDAC-cluster-NEW-NEW.pdf",width=15,height=30)
VlnPlot(scRNA, features = c("COL1A1","COL1A2","ACTA2","MCAM","CSPG4","PDGFRB","NES","RGS5","AMY2A", "AMY2B","KRT19", "KRT18","CLDN18", "MUC5AC","EPCAM", "CDH1", "CLDN3",
                           "PTPRC", "ITGAM", "CCR3","ENPP3","KIT","CD14","APOA1","FABP1","DNASE1","CFTR","SOD3","COL18A1",
                           "FUT4","CD68","HLA-DRA","HLA-DRB1","ITGAX","CD83", "CD86", "CD80","CD163", "MRC1","CD14","CD8A","CD4","EPCAM", "CDH1", "CLDN3",
                           "NCAM1","CD19","SDC1","TNFRSF17","MCAM","CSPG4",
                           "CD3D","IL17R","CD14","LYZ","NKG7","FCGR3A","GNLY","CD79A","MS4A1","TSPAN13","GPR183","PECAM1","CDH5"),pt.size=0)
VlnPlot(scRNA, features = c("Col1a1","Col1a2","Acta2","Amy2a", "Amy2b","Krt19", "Cldn18", "Muc5ac",
                           "Ptprc", "Itgam", "Ccr3","Enpp3","Kit","Cd14",
                           "Fut4","Cd83", "Cd86", "Cd80","Cd163", "Mrc1","Cd14","Cd8a","Cd4",
                           "Ncam1","Cd19","Sdc1","Tnfrsf17","Mcam","Cspg4","Pdgfrb",
                           "Nes","RgS5"),pt.size=0)
dev.off()


#Fibroblast: COL1A1, COL3A1, FAP, THY1, PDGFRB, VIM, ACTA2, CDH11, ENG; Epithelial cell: EPCAM, CDH1, CLDN3,
#CLDN4; B cell: PTPRC, CD19, CD79A, CD79B; T cell: PTPRC, CD4, CTLA4, FOXP3; Vascular endothelial cell: VIM, EMCN,
#CLEC14A, CDH5, VWF, CD34, PECAM1; Macrophage/monocyte: CD14, CD68, CSF1R, FGFR3, CD33, LYZ; Pericyte: MCAM,
#CSPG4; Lymphatic endothelial cell: LYVE1, PROX1

#CD45=PTPRC,CD11B=ITGAM,CD193=CCR3,CD203=ENPP3,CD117=KIT,CD15=FUT4,CD206=MRC1
#CD56=NCAM1 ,CD138=SDC1,BCMA=TNFRSF17
pdf(file="05.ARP.features2.pdf",width=10,height=8)
VlnPlot(scRNA, features = c("Lama5","Thbs2","Lama3","Cd44","Tnc","Thbs1","Sdc1",
                           "Npnt","Lamc2","Itgb6","Itga2","Itga3","Hmmr","Lamb3"),pt.size=0)
VlnPlot(scRNA, features = c("ACTR2", "ACTR3","ARPC1A", "ARPC1B","ARPC2", "ARPC3","ARPC4", "ARPC5"),pt.size=0)
dev.off()
VlnPlot(scRNA, features = c("Sdc4","Lamb3","Itgb8","Thbs4","Thbs2","Thbs1",
                           "Spp1","Npnt","Lamc2","Tnc","Fn1","Col1a2","Col1a1",
                           "Sdc1","Itgb6","Itgb5","Itgb1","Itgav","Itga5",
                           "Itga4","Itga3","Itga2","Cd44","Lama3","Hmmr",
                           "Col5a2",
),pt.size=0)


top10 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#????marker?ڸ???cluster????ͼ
pdf(file="06.tsneHeatmap.pdf",width=12,height=9)
DoHeatmap(object = scRNA, features = top10$gene) + NoLegend()
dev.off()
#??"ACTR2", "ACTR3","ARPC1A", "ARPC1B","ARPC2", "ARPC3","ARPC4", "ARPC5"
#??"Actr2", "Actr3","Arpc1a", "Arpc1b","Arpc2", "Arpc3","Arpc4", "Arpc5"
#????marker??С????ͼ
pdf(file="06.markerViolin-Bcell-2.pdf",width=10,height=12)
VlnPlot(object = scRNA, features = c("CD19","PTPRC","LGM","LGD","CD5","IL2RA","SPN","CXCR4","CXCR5","IL7R","IRF4","IRF8","FOXO1","PAX5","PRDM1","BCL2","BCL6","BCL2L1","LGKAPPA","LG??"))
dev.off()
pdf(file="06.markerViolin-Bcell-4.pdf",width=10,height=12)
VlnPlot(scRNA, features = c("CD19","PTPRC","LGM","LGD","CD5","IL2RA","SPN",
                            "CXCR4","CXCR5","IL7R","IRF4","IRF8","FOXO1",
                            "PAX5","PRDM1","BCL2","BCL6","BCL2L1",
                            "LGKAPPA","LG??"),pt.size=0)
dev.off()
#????marker?ڸ???cluster??ɢ??ͼ
pdf(file="06.markerScatter-Arp complex gava.pdf",width=13,height=12)
FeaturePlot(object = scRNA, features = c("Arp"),
            cols = c("gray","lightyellow", "red"))
dev.off()

pdf(file="06.CAF.markerScatter2.pdf",width=13,height=12)
FeaturePlot(object = scRNA, features = c("ACTA2", "TAGLN","IGFBP7", "IGFBP4", "IGFBP5","IGF1", "PDGFRA","PDGFB",
                                        "CXCL12", "IL6","IGF1", "CCL11","CCL7","CCL2","CXCR4",
                                        "CSF1","FAP", "CTGF", "CCL7","LIF", "CXCL1","COL1A1",
                                        "COL6A1","IL1A"), cols = c("red", "blue"))
dev.off()

#????marker?ڸ???cluster??????ͼ
library(ggplot2)
pdf(file="06.markerBubble-SHUANG-2.pdf",width=8,height=8)
cluster10Marker=c("DOT1L","ROR1","RUNX2")
DotPlot(object = scRNA, features = cluster10Marker,dot.scale = 20,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),

#DotPlot(object = scRNA, features = cluster10Marker,dot.scale = 10,col.min = -2,col.max = 2,cols = c("blue","red"))#col.min = -2,col.max = 2,cols = c("grey","red"),
dev.off()
library(ggplot2)
pdf(file="06.markerBubble-nature cancer-dims2.pdf",width=25,height=10)
cluster10Marker=c("LAMP3","CCL22","TFF1","KRT18","KRT19","KRT8","CLU","MMP7","SPP1","REG1A","CTRB2",
                  "PRSS1","DCN","LUM","CPA3","TPSAB1","CDH5","VWF","PLVAP","IRF7","RGS5","PDGFRB",
                  "CD3E","NCAM1","NKG7","CD3D","CD14","HLA-DRA","GZMB","ITGAX","ITGAM","APOE","LYZ","IGJ","CD79A","MS4A1")
DotPlot(object = scRNA, features = cluster10Marker,dot.scale = 10,col.min = -2,
        col.max = 2,)+ scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),
dev.off()
###################################07.ע??ϸ??????###################################
#????1--ͨ?????ǻ????????ף??????˹?ȷ????????Ⱥ??ϸ?????ͣ????????????ֶ?????ϸ??Ⱥ????
bfreaname.scRNA <- scRNA
#cellmarker??վ????http://biocc.hrbmu.edu.cn/CellMarker/???ֶ????????????ϱ?
#????VlnPlot(scRNA, features = top10$gene[1:20],pt.size=0)?????Ĳ???????????
new.cluster.ids <- c("CD4+ T cells","Granulocytes","CD8+ T cells", "Epithelial cells", "Dendritic cells", 
                     "Macrophages","Epithelial cells","Epithelial cells","Fibroblast","Mast cells",
                     "Macrophages","Plasma Cells","Granulocytes","Epithelial cells","Acinar cells",
                     "B cells","Epithelial cells","Epithelial cells","Pericytes","NK cells",
                     "Epithelial cells","Endothelial cells","Endocrine cells","Dendritic cells","B cells")
#peng data 
new.cluster.ids <- c("Endothelial cells","Pericytes","Macrophages", "Epithelial cells", "Fibroblast", 
                     "T cells","Epithelial cells","Fibroblast","B cells","Epithelial cells",
                     "Epithelial cells","Epithelial cells","Epithelial cells","Epithelial cells","Acinar cells","Epithelial cells",
                     "B cells","Plasma Cells","Epithelial cells","Fibroblast","Macrophages",
                     "Epithelial cells","Epithelial cells","Macrophages","Epithelial cells")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
p<-DimPlot(scRNA, reduction = "umap", label = TRUE, label.size = 6,pt.size = 0.5) + NoLegend()
ggsave("UAMP-cluster-0.25-new.pdf", p, width = 18, height = 12)#ԭ??Ϊwidth = 18, height = 12

#????2-SingleR???????????????Ǵ˷?ע?Ͳ?׼ȷ????Ҫ?ο????ݼ?
library(SingleR)
counts<-scRNA@assays$RNA@counts
clusters<-scRNA@meta.data$seurat_clusters
ann=scRNA@meta.data$orig.ident
singler = CreateSinglerObject(counts, annot = ann, "scRNA", min.genes = 0,
                              species = "Human", citation = "",#??ΪHuman??С??ΪMouse
                              ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = T, clusters = clusters, do.main.types = T,
                              reduce.file.size = T, numCores = 1)
singler$seurat = scRNA
singler$meta.data$xy = scRNA@reductions$tsne@cell.embeddings
clusterAnn=singler$singler[[2]]$SingleR.clusters.main$labels
write.table(clusterAnn,file="07.clusterAnn-PDAC.txt",quote=F,sep="\t",col.names=F)
write.table(singler$other,file="07.cellAnn-PDAC.txt",quote=F,sep="\t",col.names=F)
#?????к??ʵĲο????ݼ?ʱ?????÷??????????Ľ???ע??
#????3???Զ???singleR??ע??

#????4??????seurat???õ?ԭ??????ϸ?????ϵĹ??ܣ????ο?????????ע?????ݽ???ӳ?䴦??
#????5??ʹ??scCATCH???о????????Զ???ע??
devtools::install_github("ZJUFanLab/scCATCH")
devtools::install_github("xuzhougeng/scCATCH")
library(Seurat)
library(scCATCH)
pbmc=scRNA
clu_markers <- findmarkergenes(pbmc,species = "Human",cluster = 'All', 
                               match_CellMatch = FALSE,cancer = NULL,tissue = NULL,
                               cell_min_pct = 0.25,logfc = 0.25,pvalue = 0.05)




#׼??monocle??????Ҫ???ļ?
monocle.matrix=as.matrix(scRNA@assays$RNA@data)
monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
write.table(monocle.matrix,file="07.monocleMatrix-PDAC.txt",quote=F,sep="\t",row.names=F)
monocle.sample=as.matrix(scRNA@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="07.monocleSample-normal and tumorC.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="07.monocleGene-PDAC.txt",quote=F,sep="\t",row.names=F)
write.table(singler$other,file="07.monocleClusterAnn-PDAC.txt",quote=F,sep="\t",col.names=F)
write.table(sig.markers,file="07.monocleMarkers-PDAC.txt",sep="\t",row.names=F,quote=F)
########################monocle????
library(monocle)
setwd("C:/Users/Administrator/Desktop/Analysis/scRNA/CAF/GSE166571")               #???ù???Ŀ¼
monocle.matrix=read.table("07.monocleMatrix-PDAC.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.sample=read.table("07.monocleSample-PDAC.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.geneAnn=read.table("07.monocleGene-PDAC.txt",sep="\t",header=T,row.names=1,check.names=F)
marker=read.table("07.monocleMarkers-PDAC.txt",sep="\t",header=T,check.names=F)

#??Seurat????ת??Ϊmonocle??Ҫ??ϸ????????ϸ??ע?ͱ??ͻ???ע?ͱ???
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)

#??????һ????????????
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

#????ϸ??????????
clusterRt=read.table("07.clusterAnn-PDAC.txt",header=F,sep="\t",check.names=F)
clusterAnn=as.character(clusterRt[,2])
names(clusterAnn)=paste0("cluster",clusterRt[,1])
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

#αʱ??????????
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- setOrderingFilter(cds, marker$gene)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree')
cds <- orderCells(cds)
pdf(file="cluster.trajectory.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
dev.off()
pdf(file="cellType.trajectory.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cell_type2")
dev.off()
#####??????????αʱ??????##############vv
###########RNa???ʽ???ϸ???켣????##############
###########new monocle##############
#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                       'limma', 'S4Vectors', 'SingleCellExperiment',
#                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
#install.packages("devtools")
#devtools::install_github('cole-trapnell-lab/leidenbase')
#install.packages("E:/tumor/Matrix.utils_0.9.7.tar.gz",repos=NULL)
#devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)
library(ggplot2)
##????CDS??????Ԥ????????
scRNA=pbmc
data <- GetAssayData(scRNA, assay = 'RNA', slot = 'counts')
cell_metadata <- scRNA@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds?????൱??seurat??NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)#ԭ??Ϊnum_dim = 50
cds<-align_cds(cds, alighment_group="batch" )
plot_pc_variance_explained(cds)

#umap??ά
cds <- reduce_dimension(cds, preprocess_method = "PCA")
#cds = reduce_dimension(cds, reduction_method="tSNE")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters") + ggtitle('cds.umap')
##??seurat???????Ϲ???umap????
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')

## Monocle3????????
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)

## ʶ???켣
cds <- learn_graph(cds)
p = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = TRUE, 
               label_branch_points = TRUE)
p =plot_cells(cds, genes=c("CREB3L1", "SERPINH1", "P4HB", "IFITM5"))
ggsave("new-celltrace.pdf", p, width = 8, height = 6.6)#ԭ??Ϊwidth = 18, height = 12

#Ѱ??marker ????Find marker genes expressed by each cluster
marker_test_res = top_markers(cds, group_cells_by="cluster", 
                              reference_cells=1000, cores=8)#seurat_clusters
library(dplyr)
top_specific_markers = marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by("group") %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))
plot_genes_by_group(cds,
                    top_specific_markers,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)#seurat_clusters
top_specific_markers = marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="group",
                    ordering_type="cluster_row_col",
                    max.size=3)
#??֪ϸ??????Annotate your cells according to type
colData(cds)$assigned_cell_type=as.character(clusters(cds))
colData(cds)$assigned_cell_type = dplyr::recode(colData(cds)$assigned_cell_type,
                                                "1"="Body wall muscle",
                                                "2"="Germline","3"="Unclassified neurons","4"="Seam cells",
                                                "5"="Coelomocytes",
                                                "6"="Pharyngeal epithelia",
                                                "7"="Vulval precursors",
                                                "8"="Non-seam hypodermis",
                                                "9"="Intestinal/rectal muscle",
                                                "10"="Touch receptor neurons",
                                                "11"="Unclassified neurons",
                                                "12"="flp-1(+) interneurons",
                                                "13"="Canal associated neurons")
plot_cells(cds, group_cells_by="cluster", color_cells_by="assigned_cell_type",group_label_size=4,cell_size=1.5)
#?ֶ?ѡ??ϸ??
cds_subset = choose_cells(cds)
#Warning: package ??shiny?? was built under R version 3.5.3


##ϸ??????ʱ????
 cds <- order_cells(cds)# ????bug??ʹ?ø?????ѡ??rootϸ??
p + geom_vline(xintercept = seq(-7,-6,0.25)) + geom_hline(yintercept = seq(0,1,0.25))
embed <- data.frame(Embeddings(scRNA, reduction = "umap"))
embed <- subset(embed, UMAP_1 > -6.75 & UMAP_1 < -6.5 & UMAP_2 > 0.24 & UMAP_2 < 0.25)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
p1=plot_cells(cds, color_cells_by = "group", label_cell_groups = FALSE, 
           label_leaves = TRUE,  label_branch_points = TRUE)
ggsave("new-celltrace-pseudotime.pdf", p1, width = 6.5, height = 5)#ԭ??Ϊwidth = 18, height = 12
p1=plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = TRUE)
plot_cells(cds, color_cells_by = "cell_type", label_groups_by_cluster=FALSE,
           label_leaves=FALSE, label_branch_points=FALSE)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster=FALSE,
           label_leaves=FALSE, label_branch_points=FALSE)
p2=plot_cells(cds, color_cells_by = "seurat_clusters", label_groups_by_cluster=FALSE,
           label_leaves=TRUE, label_branch_points=TRUE)
#?ֶ?ѡ??rootOrder the cells in pseudotime
cds = order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)
# a helper function to identify the root principal points:
get_earliest_principal_node  <- function(cds, time_bin="Body wall muscle"){
  cell_ids <- which(colData(cds)[, "assigned_cell_type"] == time_bin)
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))


plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)
# 3D trajectories??4???ߣ?
cds_3d = reduce_dimension(cds, max_components = 3)
cds_3d = cluster_cells(cds_3d)
cds_3d = learn_graph(cds_3d)
cds_3d = order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d_plot_obj = plot_cells_3d(cds_3d, color_cells_by="partition")
cds_3d_plot_obj
#????????????ʱ?켣
AFD_genes = c("CREB3L1", "SERPINH1", "P4HB", "IFITM5")
AFD_lineage_cds = cds[AFD_genes,
                      clusters(cds) %in% c(11, 12, 5)]

plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="group",
                         min_expr=0.5)
#Monocle?еĲ???????????
ciliated_genes = top_specific_marker_ids[0:10]
cds_subset = cds[rowData(cds)$gene_short_name %in% ciliated_genes,]
gene_fits = fit_models(cds_subset, model_formula_str = "~seurat_clusters")
fit_coefs = coefficient_table(gene_fits)
fit_coefs
emb_time_terms = fit_coefs %>% filter(term == "seurat_clusters1")
emb_time_terms
emb_time_terms = fit_coefs %>% filter(term == "seurat_clusters1")
emb_time_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)
library(ggplot2)
p=plot_genes_violin(cds_subset[,], group_cells_by="group", ncol=3) +
  theme(axis.text.x=element_text(angle=45, hjust=1))#ncol=3ͼ??????Ŀ
ggsave("????????2.pdf", p, width = 6.5, height = 6.6)#ԭ??Ϊwidth = 18, height = 12



##Ѱ????ʱ?켣????????
#graph_test????????Ҫ?Ľ?????Ī��ָ????morans_I??????ֵ??-1??1֮?䣬0?????˻???û??
#?ռ乲????ЧӦ??1?????˻????ڿռ???????????ϸ???б???ֵ?߶????ơ?
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)
#??ѡtop10??ͼչʾ
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#????????????ͼ
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="predicted.id", 
                         min_expr=0.5, ncol = 2)
#FeaturePlotͼ
plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
##Ѱ?ҹ?????ģ??
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=colData(cds)$predicted.id)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")



##############??λ???Լ?????Ȥ????Ⱥ???????·???################
Idents(scRNA)
levels(scRNA)
head(scRNA@meta.data)
scRNA=scRNA
scRNA = scRNA[,scRNA@meta.data$seurat_clusters %in% c(4,5)]
scRNA = scRNA[,scRNA@meta.data$group %in% c("tumor")]
scRNA<- merge(x = scRNA1, y = scRNA2,
                   add.cell.ids=c("scRNA1","scRNA2"),merge.data=TRUE, project="10X_HSC")
scRNA=Fibroblast
#ʣ?µ??߱?׼???뼴??,?????ݱ?׼????ʼ??Ȼ?󷵻????????з?????ע???޸?·??
scRNA <- NormalizeData(object = Fibroblast1, Normalization.method = "LogPDACize", scale.factor = 10000)
scRNA <- NormalizeData(object = Fibroblast1, scale.factor = 10000)
dim(scRNA)
save(scRNA,file='caner and fibroblast.Rdata')
load("Fibroblast-11.Rdata")
saveRDS(scRNA, "macrohage.rds") 

##################################??????????######################################################
new.cluster.ids <- c("0", "1", "1","0")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
P=DimPlot(scRNA, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("UAMP-cluster.pdf", p, width = 9, height = 6)

table(scRNA$seurat_clusters)
scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, 
                                min.pct = 0.25, logfc.threshold = 0.25)#ԭ??Ϊ0.25 
write.table(scRNA.markers,file="06.markers.xls",sep="\t",row.names=F,quote=F)

new.cluster.ids <- c("0", "1", "1","0")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pdf(file="06.TSNE-PDAC-0.04.pdf",width=6.5,height=5)
TSNEPlot(object = pbmc, pt.size = 0.3, label = TRUE)#TSNE???ӻ?,pt.size???ĵ?????С
dev.off()

#########?Զ???ĳһ?????Ի?????????????####
library(ggplot2)
pbmc<- readRDS('Fibroblast.RDS')
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:20)#dims = 1:20,ԭ??Ϊ30
pbmc <- FindClusters(pbmc, verbose = FALSE,resolution = 0.1)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:30)
DimPlot(pbmc, reduction = "umap", label = TRUE,pt.size=0.8)
scRNA <- pbmc
mallignant <- read.delim("each cell CNV scores.txt", row.names = 1)#?ļ?Ϊ?ֶ???????СȺ?ϲ?
scRNA <- AddMetaData(scRNA, metadata = mallignant)
p1 <-DimPlot(scRNA, reduction = "umap", label = TRUE,pt.size=1.0)
p2 <- DimPlot(scRNA, group.by = "group") + scale_color_manual(values = c("gray","red","green"))
p2 <-FeaturePlot(object = scRNA, features = c("Arp"), #GSVA��????????Ҳ??????
                 cols = c("gray","lightyellow", "red"))
pc <- p1 + p2
ggsave("cnv.pdf", pc, width = 12, height = 5)
pbmc.markers <- FindMarkers(scRNA, only.pos = FALSE,group.by = "group",
                               min.pct = 0.05, logfc.threshold = 0.25,
                            ident.1="tumor",ident.2="normal")

write.table(pbmc.markers,file="06.diff genes-tumor vs normal.xls",sep='\t',quote=F,row.names=T)

pdf(file="06.markerViolin.pdf",width=12,height=12)
VlnPlot(object = scRNA, features = c("APOA1","FABP1","DNASE1","CFTR","SOD3","COL18A1"),group.by = "group",pt.size=0)
dev.off()
markers <- c("ITGA1","ITGA2","ITGA2B","ITGA3","ITGA4","ITGA5","ITGA6","ITGA7","ITGA8","ITGA9",
             "ITGA10","ITGA11","ITGAD","ITGAE","ITGAL","ITGAM","ITGAV","ITGAX","ITGB1","ITGB2",
             "ITGB3","ITGB4","ITGB5","ITGB6","ITGB7","ITGB8")
markers <-data.table::fread("markers.txt", data.table = F)
markers <- as.data.frame(markers)
scRNA=pbmc
markerdata <- ScaleData(scRNA, features = as.character(unique(markers$markers)), assay = "RNA")
DoHeatmap(markerdata,group.by = "seurat_clusters",
          features = as.character(unique(markers$markers)))+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
library(ggplot2)
cluster10Marker=c("ITGA1","ITGA2","ITGA2B","ITGA3","ITGA4","ITGA5","ITGA6","ITGA7","ITGA8",
                  "ITGA9","ITGA10","ITGA11","ITGAD","ITGAE","ITGAL","ITGAM","ITGAV","ITGAX",
                  "ITGB1","ITGB2","ITGB3","ITGB4","ITGB5","ITGB6","ITGB7","ITGB8")
cluster10Marker=c("ITGA2","ITGA3","ITGA4","ITGA5","ITGAV",
                  "ITGB1","ITGB5","ITGB6","ITGB8")
cluster10Marker=c("COL5A2","HMMR","LAMA3","CD44","ITGA2","ITGA3","ITGA4","ITGA5","ITGAV",
                  "ITGB1","ITGB5","ITGB6","SDC1","COL1A1","COL1A2","FN1","TNC","LAMC2","NPNT",
                  "SPP1","THBS1","THBS2","THBS4","ITGB8","LAMB3","SDC4")
DotPlot(object = pbmc, features = cluster10Marker,dot.scale = 15,col.min = -2,
        col.max = 2,)+ RotatedAxis()+scale_color_gradient2(low = "blue", mid = "white", high = "red")#cols = c("blue","red"),
DotPlot(object = scRNA, group.by = "group0.25",features = cluster10Marker,dot.scale = 7,col.min = -2,
        col.max = 2,)+ RotatedAxis()+scale_color_gradient2(low = "blue", mid = "white", high = "red")+ coord_flip()+ theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

##################################????ͼ######################################
library(ggplot2)
markers <- c("ALDH1A3","EMP1","FAM3C","IRS2","MAML2","MCC","PMEPA1","SP100")
markers <- c("Aldha3")
markers <- 
markers <- as.data.frame(markers)
markerdata <- ScaleData(pbmc, features = as.character(unique(markers$markers)), assay = "RNA")
DoHeatmap(markerdata,annotation_col=group,
          features = as.character(unique(markers$markers)),
          group.by = "seurat_clusters",
          assay = 'RNA',
          group.colors = c("#00BFC4","#AB82FF","#00CD00","#C77CFF"))+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

###?Զ???˳??,????Ҫ??ѡ????Ҫ????Ⱥ
markers <- c("Sparcl1","Col15a1","Cxcl14","Serpine2","Col4a2","Col4a1","Col4a2","Eng","Dpp4","Il33",
             "Ackr3","Cd248","Sfrp2","Ly6c1","Scarma3")
markers <- c("Tnc","Cd44","Col1a1","Col1a2","Col5a2","Fn1","Hmmr","Itga2","Itga3","Itga4",
             "Itga5",	"Itgav","Itgb1","Itgb5","Itgb6","Itgb8","Lama3","Lama5","Lamb3","Lamc2",
             "Npnt","Sdc1",	"Sdc4","Spp1","Thbs1","Thbs2","Thbs4")
markers <- c("Aldh1a3","Tgfb1","Actr3","Lrrc15","Itga5","Itgb5")
markers <- as.data.frame(markers)
markerdata <- ScaleData(pbmc, features = as.character(unique(markers$markers)), assay = "RNA")
markerdata$seurat_clusters <- factor(x=markerdata$seurat_clusters,
                              levels = c(2,1,0))
DoHeatmap(markerdata,
          features = as.character(unique(markers$markers)),
          group.by = "seurat_clusters",
          assay = 'RNA',)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
########v????ͼ?л??漰??��????ɫ??һ????��??ʾ????��??��???Ա仯??һ??????չʾ???ࡣ
########??һ????????R?????????ڴ?????ɫ??????Github??ַΪ
#devtools::install_github("caleblareau/BuenColors")
library("BuenColors")
col <- jdb_color_maps[1:9]
names(col) <- levels(cluster_info)
Heatmap(mat,col=c("lightblue","red"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info)

top_anno <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = col), # ????????ɫ
                       labels = levels(cluster_info), 
                       labels_gp = gpar(cex = 0.5, col = "white"))) # ????????

Heatmap(mat,col=c("lightblue","red"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        column_split = cluster_info,
        top_annotation = top_anno, # ????ͼ?ϱ?????ע??
        column_title = NULL ) # ????Ҫ?б???
##############???ڻ????ܶ?ֱ??չʾ??��???????????壬???ǿ???ǿ?????????ǻ???######
######??pbmc.markers??ʼ
mat <- GetAssayData(pbmc, slot = "counts")
mat <- log2(mat + 1)
#??ȡ??????ϸ????????Ϣ
gene_features <- top10
cluster_info <- sort(pbmc$seurat_clusters)
#?Ա???��??????????????ɸѡ,

mat <- as.matrix(mat[pbmc.markers$gene, names(cluster_info)])
mark_gene <- c("Ly6c1","Pi16","Serpine2","Gpx3","Tagln","Cthrc1","Tmem100","Dpp4",
               "Smoc2","Col15a1","Crabp1","Crabp2","Upk3b","Lrrn4","Krt8","Krt18",
               "Cxcl1","Has1","Top2a","Mki67","C1qb","Csf1r","Pnlip")
gene_pos <- which(rownames(mat) %in% mark_gene)

row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, 
                                                 labels = mark_gene))
top_anno <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = col), # ????????ɫ
                       labels = levels(cluster_info), 
                       labels_gp = gpar(cex = 0.5, col = "white"))) # ????????
Heatmap(mat,col=c("lightblue","red"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        top_annotation = top_anno,
        right_annotation = row_anno,
        column_title = NULL)


Heatmap(mat,col=c("lightblue","red"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        top_annotation = top_anno,
        right_annotation = row_anno,
        column_title = NULL,
        heatmap_legend_param = list(
          title = "log2(count+1)",
          title_position = "leftcenter-rot"
        ))
############??ϸ??ת¼??-CellChatϸ???໥???÷???#####################
#??????Ҫ???Ѿ???????ϸ????Ⱥ???õ???seurat object????ʱ????Ҫ????һ??CellChat??????
# Part I: Data input & processing and initialization of CellChat object
#devtools::install_github("sqjin/NMF")
#BiocManager::install("ComplexHeatmap", update = F)
#library(NMF)
#library(ComplexHeatmap)
#library(dplyr)
#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
setwd("C:/Users/Administrator/Desktop/Analysis/scRNA/PDAC2")
scRNA <- readRDS("tumor-macrophage-harmony.rds")
scRNA <- FindNeighbors(scRNA, dims = 1:30) %>% FindClusters(resolution = 0.08)#resolution????????????
options(repr.plot.height = 4, repr.plot.width = 6)#ԭ??Ϊ4,6
DimPlot(scRNA, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#??????????Ϊ???ų?label=0
new.cluster.ids <- c( "Cancer cell 1","Cancer cell 2", "Macrophagy","Cancer cell 3","Cancer cell 4")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
DimPlot(scRNA, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#????CellChat
data.input = scRNA[["RNA"]]@data  # normalized data matrix
labels <- Idents(scRNA)
identity <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
#identity = data.frame(group = scRNA$seurat_clusters, row.names = names(scRNA$seurat_clusters)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels
cellchat <- createCellChat(data.input)
cellchat
#metadata??Ϣ?ӵ?CellChat??????
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")

cellchat <- setIdent(cellchat,ident.use = "labels") # set "labels" as default cell identity--
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
#?????????????ݿ?
CellChatDB <- CellChatDB.human#CellChatDB.mouseΪ?????ݿ?
#???ݿ?????Ϣ?Ǻ?ȫ????
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
showDatabaseCategory(CellChatDB)
#??CellChat?У????ǻ??????????ض?????Ϣ????ϸ???????໥????
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)
#?Ա??????ݽ???Ԥ???���????ϸ??????ͨ?ŷ???????????һ??ϸ??????ʶ??????????????
#?????壬Ȼ?󽫻???????????Ͷ?䵽????-?????໥????(PPI)?????ϡ?????????????????
#?????ʶ??????????????????֮?????໥???á?
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel  ?????ƺ???һЩbug????Linux?Ͼ?Ȼ???С?de??????
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  
#Ȼ????????ͨ??Ϊÿ???໥???÷???һ??????ֵ???????û?????��?ƶ??????????ϵ?ϸ??
#-ϸ??ͨ??
cellchat <- computeCommunProb(cellchat,raw.use = FALSE,population.size = TRUE)
#cellchat <- computeCommunProb(cellchat)#ע?????????????????????þ??ã??????????ߵġ?
#mycomputeCommunProb <-edit(computeCommunProb)  # computeCommunProb?ڲ??ƺ???һЩbug??ͬһ????????window10??û?£?????Linux???б??���??????computeExpr_antagonist?????????????⣬(matrix(1, nrow = 1, ncol = length((group))))????ӦΪ(matrix(1, nrow = 1, ncol = length(unique(group))))?? ??Ȼ???󷵻صĲ??ԡ?de??????
#environment(mycomputeCommunProb) <- environment(computeCommunProb)
#cellchat <- mycomputeCommunProb(cellchat)  # ????????de???ġ?
cellchat <-filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr.csv")
#?ƶ??ź?ͨ·ˮƽ??ϸ??ͨѶ???磨??????????@netP???棬??һ??????ֵ?Ͷ?Ӧ??pval??
#ͨ??????��·????��??????ͨ?Ÿ???��????ϸ?????ľۺ?ͨ??????
cellchat <- computeCommunProbPathway(cellchat)
cellchat@netP$pathways#?鿴????
head(cellchat@LR$LRsig)
df.netp <-subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway.csv")
##########?��ˣ?????ͳ?Ʋ????Ѿ????ɣ???һ???????ݿ??ӻ?###############
#1??????ϸ??Ⱥ?????ۣ?ϸ????????��??ǿ??ͳ?Ʒ???
#ͳ??ϸ????ϸ??֮??ͨ?ŵ???��???ж??ٸ?????-?????ԣ???ǿ?ȣ????ʣ?
cellchat <- aggregateNet(cellchat)
#????ÿ??ϸ?????ж??ٸ?
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf(file="interactions.pdf",width=6.5,height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name ="Interaction weights/strength")
dev.off()
#????ÿ??ϸ?????????źţ?ÿ??ϸ?????θ?????ϸ????????number of interactionͼ??
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) { mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[i, ] <- mat[i, ]
netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                 arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
#ÿ??ϸ?????θ?????ϸ??????????????ǿ??/????ͼ??
mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) { mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[i, ] <- mat[i, ]
netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                 arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
#2???????ź?ͨ·??????-?????鵼??ϸ?????????ӻ???????ͼ??????ͼ??????ͼ????ͼ=չʾ?????ݺʹ???????˼һģһ????
cellchat@netP$pathways#?鿴??????Щ?ź?ͨ·
pathways.show <- c("IGF")
levels(cellchat@idents)# show all celltype
vertex.receiver = c(1,2,4,6)# a numeric vector??vertex.receiver = seq(1,4)
#????????????Ȧͼ
netVisual_aggregate(cellchat, signaling = c("IGF"),
                    layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
#??ͼ
par(mfrow=c(1,1))
netVisual_heatmap(cellchat,signaling = pathways.show,color.heatmap = "Reds")
#3.???׶?TU?????Ϳ??ӻ?ÿ??????-???????????ź?ͨ·?Ĺ??׶?
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.TGFb <- extractEnrichedLR(cellchat,signaling = pathways.show, geneLR.return = FALSE)
#??ȡ??????ͨ·????????????????????��չʾ??Ҳ????ѡ?????????????????ԣ?
LR.show <- pairLR.TGFb[1,] 
vertex.receiver =  c(1,2,4,6) # a numeric vector
netVisual_individual(cellchat,signaling = pathways.show, pairLR.use = LR.show,vertex.receiver = vertex.receiver)
#????ͼ??Circle plot??
netVisual_individual(cellchat,signaling = pathways.show, pairLR.use = LR.show,layout = "circle")
#4???Զ?????��??????ÿ???ź?ͨ·?Ļ???????
pathways.show.all <-cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver = c(1,2,4,6)
dir.create("all_pathways_com_circle")
setwd("all_pathways_com_circle")
for (i in 1:length(pathways.show.all)){
  netVisual(cellchat, signaling= pathways.show.all[i],out.format = c("pdf"),
            vertex.receiver =vertex.receiver,layout = "circle")
  gg <- netAnalysis_contribution(cellchat, signaling =pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"),
         plot=gg, width = 5,height = 2.5, units = 'in',dpi = 300)}
setwd("../")
#5??????????-?????鵼??ϸ????????ϵ???ӻ?
#????ͼ??ȫ?????????壩
levels(cellchat@idents)
p = netVisual_bubble(cellchat, sources.use =c(3,5,7,8,9), 
                     targets.use = c(1,2,4,6),remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble.pdf",p, width = 8, height = 12)#??ϵ???ܰ͵ĵ???
#????ͼ??ָ???ź?ͨ·??????-???壩
#?????ƶ?CCL??CXCL??��???ź?ͨ·
netVisual_bubble(cellchat, sources.use = c(3,5,7,8,9),targets.use = c(1,2,4,6),
                 signaling = c("CCL","CXCL"),remove.isolate = FALSE)
#????ͼ??ָ???ź?ͨ·??????-???岢ָ??ϸ????
pairLR.use <- extractEnrichedLR(cellchat, signaling =c("CCL","CXCL","TGFb"))
netVisual_bubble(cellchat,sources.use = c(3,6,8),targets.use = c(1,4,5),
                 pairLR.use = pairLR.use,remove.isolate = TRUE)
#????ĳ???ź?ͨ·????TGFb???????л?????ϸ??Ⱥ?еı???????չʾ??С????ͼ??????ͼ??
p = plotGeneExpression(cellchat, signaling = "TGFb")
ggsave("TGFb_GeneExpression_vln.pdf", p, width = 8,height = 8)
p = plotGeneExpression(cellchat, signaling = "TGFb", type = "dot")
ggsave("TGFb_GeneExpression_dot.pdf", p, width= 8, height = 6)









#???ݿ??ӻ?
#1.??????Ȧͼ??ʵ??Բ?Ϳ???Բ?ֱ???ʾԴ??Ŀ?ꡣԲ?Ĵ?С??ÿ??ϸ??????ϸ?????ɱ?????
#??Ե??ɫ????Դһ?¡???Խ?֣??ź?Խǿ??????????չʾ??һ??MIF?ź??????????ӣ?
#??????ʾ??Ҫͨ?ŵ?????·????????ͨ??cellchat@netP$pathways???ʡ?
levels(cellchat@idents) 
vertex.receiver = seq(1,4) # a numeric vector
pathways.show <-"IL6" #MIF?ź?ͨ·
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize) 
#????????????Ȧͼ
netVisual_aggregate(cellchat, signaling = c("IL6"),
                    layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
#2.???׶?TU?????Ϳ??ӻ?ÿ??????-???????????ź?ͨ·?Ĺ??׶?
netAnalysis_contribution(cellchat, signaling = pathways.show)
#ʶ??ϸ??Ⱥ???ź?ת?????ã?ͨ??????ÿ??ϸ??Ⱥ????????????ָ?꣬CellChat??????ʱ
#ʶ??ϸ????ͨ???????е???Ҫ?????ߡ??????ߡ??????ߺ?Ӱ???ߡ?
cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
