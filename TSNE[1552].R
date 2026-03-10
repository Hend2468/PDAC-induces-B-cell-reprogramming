#install.packages("devtools")
#library(devtools)
#devtools::install_github('dviraran/SingleR')#ע???ð?
###################################04.????ǰ?ڴ????ͽ???###################################
#??ȡ????
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(data.table)
setwd("C:/Users/Administrator/Desktop/Analysis/scRNA/CAF/GSE166571-2")  #?趨?Լ??Ĺ???Ŀ¼            #???ù???Ŀ¼

#??ȡ?ļ????????ظ?????ȡ??ֵ

rt=read.table("final count_matrix.txt",sep="\t",header=T,check.names=F)
#?????Ǿ???????
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#?????ļ?̫???ֿ???ȡ?ٺϲ?
setwd("C:/Users/Administrator/Desktop/Analysis/scRNA/GSE129455Bai/KPC/fibroblast")
rt2=read.table("gene name-2.txt",sep="\t",header=T,check.names=F)
#?????Ǿ???????
rt2=as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp2=rt2[,2:ncol(rt2)]
dimnames2=list(rownames(exp2),colnames(exp2))
data2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
data2=avereps(data2)
#??һ?????е?pbmc???ٺϲ?

#data=as.data.frame(t(data))#?ϴ??????ô?ת??????????avereps??????external_gene_name.y


#????Ϊ̽????????????bioMart?Ի???????ע?ͣ??????ο?RNAseq????
library('biomaRt')
library("curl")
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(rt)
#listAttributes(mart)
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                    filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(mms_symbols)
#?ϲ?????:res????+mms_symbols?ϲ???һ???ļ?
head(diff_gene_deseq2)
#?ɼ???��???ļ?û?й?ͬ????????????Ҫ?ȸ?'diff_gene_deseq2'????һ??
#??ensembl_gene_id??????????????????:(Ӧ???и??????ķ???)
ensembl_gene_id<-rownames(rt)
exp<-cbind(my_ensembl_gene_id,rt)
colnames(exp)[1]<-c("ensembl_gene_id")
diff_name<-merge(exp,mms_symbols,by="ensembl_gene_id")
head(diff_name)
write.table(diff_name,"gene name.txt",sep ="\t",col.names=NA,quote = FALSE)


#??һ?ֶ?ȡ???ݵķ?ʽ???ļ????е?3???ļ??ֱ??ǣ?barcodes.tsv.gz??
#genes.tsv.gz(or features.tsv.gz)??matrix.mtx.gz???ļ??????ֲ??ܴ���??????ȡ????
library(dplyr)
library(Seurat)
library(patchwork)
data <- Read10X(data.dir = "E:/tumor/GSE196678")
#?????ַ????Ƕ??????ϲ?
library(Seurat)
samples=list.files("GSE162454/")
samples
dir <- file.path('./GSE162454',samples)
names(dir) <- samples
#?ϲ?????
data <- Read10X(data.dir = dir)
#data = CreateSeuratObject(counts, min.cells=1)
dim(data)   #?鿴????????ϸ??????
#table(data@meta.data$orig.ident)  #?鿴ÿ????????ϸ????
#?????ַ???????????loom?ļ?
#devtools::install_github(repo = "hhoeflin/hdf5r")
#devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
setwd("E:/tumor/bone/longbone")
library(loomR)
library(dplyr)
lfile <- connect(filename = "GSM4274191_CS20_longbone_rawdata.loom", mode = "r+")
lfile
lfile$matrix
full.matrix <- lfile$matrix[, ]
dim(x = full.matrix)
A <-lfile$col.attrs$CellID[]
B <-lfile$col.attrs$ClusterName[] 
C <-lfile$col.attrs$ClusterID[]
D <-lfile$row.attrs$Gene[]
#############ת??Seurat?????ķ?ʽ##############
full.matrix <- lfile$matrix[,]
dim(full.matrix)
full.matrix <- t(full.matrix)
gene_name <- lfile$row.attrs$Gene[] #??????
barcode <- lfile$col.attrs$CellID[] #ϸ??ID

colnames(full.matrix) <- barcode
rownames(full.matrix) <-gene_name

S =CreateSeuratObject(counts = full.matrix,project = 'SP',min.cells = 0, min.features = 0)
S$ClusterID  <- lfile$col.attrs$ClusterID[]
S$ClusterName <- lfile$col.attrs$ClusterName[]
S$Regin <- lfile$col.attrs$Region[]
S$Total_molecules  <- lfile$col.attrs$Total_molecules[]
S$Valid <- lfile$col.attrs$Valid[]
S$X <- lfile$col.attrs$X[]
S$Y <-lfile$col.attrs$Y[]
S$size_pix <- lfile$col.attrs$size_pix[]
S$size_um2 <- lfile$col.attrs$size_um2[]

S@assays$RNA@meta.features$Fluorophore <- lfile$row.attrs$Fluorophore[]
S@assays$RNA@meta.features$Hybridization <- lfile$row.attrs$Hybridization[]

#remotes::install_github("aertslab/SCopeloomR")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("LoomExperiment")

library(hdf5r)
library(loomR)
library(LoomExperiment)
library(SCopeLoomR)
path="E:/tumor/bone/longbone"
samples=list.files(path) 
samples
sceList=lapply(samples, function(pro){
  # pro=samples[1]
  print(pro)
  folder=file.path(path,pro)
  print(pro)
  print(folder)
  #????loom?ļ?
  lfile <- connect(filename =folder, mode = "r+")
  ct <- as.data.frame(lfile[["matrix"]][,])
  colnames(ct) <- lfile[["row_attrs/Gene"]][] #?ż?gene_ name
  rownames(ct) <- lfile[["col_attrs/CellID"]][] #?ż?cell_ name
  ct[1:5,1:5]
  ct <- t(ct) #????ת??
  sce=CreateSeuratObject(counts = ct,
                         project = pro )
  return(sce)
})
names(sceList)
samples
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = samples)
pbmc=sce.all
#as.data.frame(sce.all@assays$RNA@counts[1:10, 1:2])

#head(sce.all@meta.data,10)
#table(sce.all@meta.data$orig.ident)
#library(stringr)
#phe=str_split(rownames(sce.all@meta.data),'-',simplify = T)
#head(phe)
#table(phe[,4])
#table(phe[,3])
#sce.all@meta.data$rep=phe[,4]
#sce.all@meta.data$group=phe[,3]
#sce.all@meta.data$orig.ident = paste0(
#sce.all@meta.data$group,
#sce.all@meta.data$rep)
#table(sce.all@meta.data$orig.ident)



#load("pbmc-δ??׼??.Rdata")
dim(data)
dim(pbmc)
setwd("E:/tumor/GSE196678")  #?趨?Լ??Ĺ???Ŀ¼
#??????ת??ΪSeurat???󣬲??????ݽ??й???delim = "_",
pbmc <- CreateSeuratObject(counts = data,project = "seurat", min.cells = 3, min.features = 200,names.delim = "_", )#ԭ��Ϊ3,50
#ʹ??PercentageFeatureSet??????????��???????İٷֱ?#ʹ??PercentageFeatureSet??????????��???????İٷֱ?#ʹ??PercentageFeatureSet??????????��???????İٷֱ?
#save(pbmc,file='pbmc2.Rdata')
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pdf(file="04.featureViolin-PDAc-2.pdf",width=10,height=6)           #????????????С????ͼ
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0.1)
dev.off()
#??????��???????????ʿ?
#load("pbmc-δ??׼??.Rdata")
memory.limit(size=60000)
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 350 & nFeature_RNA < 5000&percent.mt < 10) #?????ݽ??й???,ԭ??Ϊ50??5
dim(pbmc)
#??һ?ַ?????ȡ???????е?PBMC???кϲ???????
#pbmc<- merge(x = pbmc1, y = pbmc2,
#                  add.cell.ids=c("pbmc1","pbmc2"),merge.data=TRUE, project="10X_HSC")
#dim(pbmc)
#save(pbmc,file='pbmc-?ʿ?δ??׼??.Rdata')
#???????ȵ??????Ի?ͼ
pdf(file="04.featureCor-PDAC.pdf",width=10,height=6)              #??????????????????ͼ
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#?????ݽ??б?׼??
pbmc <- NormalizeData(object = pbmc, Normalization.method = "LogPDACize", scale.factor = 10000)
#pbmc <- NormalizeData(object = pbmc, scale.factor = 10000)


#??ȡ??Щ??ϸ????????ϵ???ϴ??Ļ???
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)#ԭ??Ϊ1500
#????????????ͼ
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="04.featureVar-PDAC.pdf",width=10,height=6)              #????????????????ͼ
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

###################################05.PCA???ɷַ???###################################
##PCA????
#PCA??ά֮ǰ?ı?׼Ԥ???���?裬????ȥ??????ЧӦ
#????һ
pbmc=ScaleData(pbmc) #PCA??ά֮ǰ?ı?׼Ԥ???���?裬????ȥ??????ЧӦ
pbmc=RunPCA(object= pbmc,npcs = 30,pc.genes=VariableFeatures(object = pbmc)) 
#??????
#pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nCount_RNA", "percent.mt"))
#summary(pbmc@scale.data[,1])
#pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
#????3??????ϸ?????ڻ???
if(T){
  g2m_genes <- cc.genes$g2m.genes
  g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(pbmc))
  s_genes <- cc.genes$s.genes    
  s_genes <- CaseMatch(search=s_genes, match=rownames(pbmc))
  pbmc <- CellCycleScoring(pbmc, g2m.features=g2m_genes, s.features=s_genes)
  tmp <- RunPCA(pbmc, features = c(g2m_genes, s_genes), verbose = F)
  p <- DimPlot(tmp, reduction = "pca", group.by = "orig.ident")
  ggsave("QC/CellCycle_pca.png", p, width = 8, height = 6)
  rm(tmp)
}
pbmc <- ScaleData(pbmc, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
pbmc <- RunPCA(pbmc, features = c(s_genes, g2m_genes))
DimPlot(pbmc,group.by = 'Phase')
#PCA????
#????ÿ??PCA?ɷֵ????ػ???
pdf(file="05.pcaGene-PDAC.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()

#???ɷַ???ͼ??
pdf(file="05.PCA-PDAC.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()

#???ɷַ?????ͼ
pdf(file="05.pcaHeatmap-PDAC.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#ÿ??PC??pֵ?ֲ??;??ȷֲ?
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="05.pcaJackStraw-PDAC.pdf",width=8,height=6)
JackStrawPlot(object = pbmc, dims = 1:20)
dev.off()
saveRDS(pbmc, "pbmc.rds")
###################################06.TSNE??????????marker????###################################
##TSNE??????
#pbmc=Fibroblast
ElbowPlot(object = pbmc)#???????ʯͼ��ȷ?????ɷ֣???JackStraw???ܲ??
pcSelect=20#ԭ??Ϊ20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)  #?????ڽӾ???,ԭʼΪ1
pbmc <- FindClusters(object = pbmc, resolution = 0.1) # ԭ��resolution = 0.5, ?ɵ??????????٣???ϸ??????,?Ż???׼ģ?黯
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect,check_duplicates=FALSE)#TSNE????.check_duplicates=FALSE
pdf(file="06.TSNE-PDAC-0.5.pdf",width=6.5,height=5)
TSNEPlot(object = pbmc, pt.size = 0.3, label = TRUE)#TSNE???ӻ?,pt.size???ĵ?????С
dev.off()

##UMAP??????
setwd("D:/Analysis/scRNA/PDAC3/Macrophage")
pbmc <- readRDS("Oncocytes.rds")
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:40)#dims = 1:20,ԭ??Ϊ30
pbmc <- FindClusters(pbmc, verbose = FALSE,resolution = 0.04)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:40)
pdf(file="06.umap-40-0.04.pdf",width=6.5,height=5)
DimPlot(pbmc, reduction = "umap", label = TRUE,pt.size=1.2)
dev.off()



write.table(pbmc$seurat_clusters,file="06.tsneCluster-PDAC.txt",quote=F,sep="\t",col.names=F)

##Ѱ?Ҳ?????????????
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, # only.pos = TRUE
                               min.pct = 0.05, logfc.threshold = 0.25)#ԭ??Ϊ0.25 
write.table(pbmc.markers,file="OS.xls",sep="\t",row.names=F,quote=F)
library(dplyr)
sig.markers=pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(file="05.DoHeatmap-PDAC.pdf",width=10,height=8)
DoHeatmap(pbmc, features = top10$gene) #չʾǰ10?????ǻ???????ͼ
dev.off()
pdf(file="05.DoHeatmap-cancer cell.pdf",width=10,height=8)
pbmc@assays$RNA@scale.data <- scale(pbmc@assays$RNA@data, scale = TRUE)
DoHeatmap(pbmc, features = c("Ly6c1","Pi16","Serpine2","Gpx3","Tagln","Cthrc1","Tmem100","Dpp4",
                             "Smoc2","Col15a1","Crabp1","Crabp2","Upk3b","Lrrn4","Krt8","Krt18",
                             "Cxcl1","Has1","Top2a","Mki67","C1qb","Csf1r","Pnlip"),slot = "data")+
dev.off()
pdf(file="05.top20 markers.pdf",width=10,height=8)
VlnPlot(pbmc, features = top10$gene[1:20],pt.size=0)
dev.off()
pdf(file="05.features-cancer cell.pdf",width=15,height=20)
VlnPlot(pbmc, features = c("ALDH1A3","EMP1","FAM3C","IRS2","MAML2","MCC","PMEPA1","SP100"),pt.size=0)
dev.off()
pdf(file="05.cancer cell-marker.pdf",width=15,height=20)
VlnPlot(pbmc, features = c("TGFB2","IFIT3","S100A2","MT2A","CAV1","CD44","COL6A1","TGFBI",
                           "WNT7B","DKK1","KRT6A","TP63","FOXA3","GATA6","PDX1","MUC5B",
                           "MUC5AC","CEACAM6","LYZ","TFF2","EPCAM"),pt.size=0)
dev.off()
pdf(file="05.features-OS-2.pdf",width=15,height=25)

VlnPlot(pbmc, features = c("Acta2", "Tagln","Igfbp7", "Igfbp", "Igfbp5","Igf1", "Pdgfra","Pdgfb",
                           "Cxcl12", "Il6", "Ccl11","Ccl7","Ccl2","Il1a",
                           "Csf1","Fap", "Ctgf", "Ccl7","Lif", "Cxcl1","Ccr2","Col1a1",
                           "Col6a1","Clec3b","Has1","Col4a1","Cd74","Saa3","H2-Ab1","Sjpi","Fabp4","Car3"),pt.size=0)
VlnPlot(pbmc, features = c("Clec3b","Col14a1","Has1", "Il6","H2-Ab1","Cd74","Saa3","Slpi",
                           "Tagln","Thy1","Col12a1","Thbs2","Col1a1","Col1a2","Pdpn","Dcn"),pt.size=0)
VlnPlot(pbmc, features = c("ACTA2", "TAGLN","IGFBP7", "IGFBP4", "IGFBP5","IGF1", "PDGFRA","PDGFB",
                           "CXCL12", "IL6","IGF1", "CCL11","CCL7","CCL2","CXCR4",
                           "CSF1","FAP", "CTGF", "CCL7","LIF", "CXCL1","CCR2","COL1A1",
                           "COL6A1","IL1A","IL17R","FABP4","CAR3","CSPG4","PDGFRB","NES","RGS5",),pt.size=0)

VlnPlot(pbmc, features = c("COL1A1","COL1A2", "FAP","PDPN","DCN","LUM","VIM","ACTA2", "TAGLN","MMP11","MYL9","HOPX","POSTN","TPM1","TPM2",
                           "IL6", "PDGFRA","CXCL12","CFD","DPT","LMNA","AGTR1","HAS1","CXCL1","LRRC15",
                           "CXCL2","CCL2","IL8","COL14A1","CLEC3B","CLEC3B","IGF1","GSN","H2-AB1","CD74","SAA3","SLPI","HLA-DRB1","HLA-DRA","HLA-DPA1","HLA-DOA1",
                           "LY6C1","CSPG4","PDGFRB","NES","RGS5","FABP4","CAR3","PTPRC","EPCAM"),pt.size=0)

VlnPlot(pbmc, features = c("COL1A1","SPP1" ,"LUM","CTSK","IL7R",
                           "CD3","DNKG7","CD74","CD14","FCGR3A","RGS5","ACTA2","CXCL12","MME",
                           "SFRP2","vWF", "MYLPF" ,"RUNX2","ALPP","ALPL","COL2A1","SOX9","ACAN",
                           "CALCR","CTSK","MMP9", "TRAP","TNFRSF11A","SOFAT","SOST","CX43","DMP1","FAP","DPPIV"),pt.size=0)



dev.off()
pdf(file="05.iCAF.myCAF-features-new-3.pdf",width=15,height=12)
VlnPlot(pbmc, features = c( "DCN","PDPN","FAP","LUM","PDGFRA","C7","MGP","COL12A1","COL11A1","LOXL2","COL14A1","SDC1","NFIA",
                             "SPARCL1","OGN","LGALS1","C5orf46","ANTXR1","C3","ABI3BP","MMP11","SLC16A3","IGFBP4","COL5A2","CTHRC1"),pt.size=0)
dev.off()
pdf(file="05.CAF.MOUSE.pdf",width=15,height=12)
VlnPlot(pbmc, features = c("Clec3b","Col14a1","Has1", "CCL2","H2-Ab1","Cd74","Saa3","Slpi","Acta2",
                           "Tagln","Thy1","Col12a1","Thbs2","Col1a1","Col1a2", "Fap","Pdpn","Dcn","Vim,"),pt.size=0)
VlnPlot(pbmc, features = c("Clec3b","C7","Col14a1","Cxcl12", "Ly6a","Dpt","H2-Ab1","Cd74","Saa3","Slpi","Acta2",
                           "Tagln","Thy1","Col12a1","Thbs2","Col1a1","Col1a2", "Fap","Pdpn","Dcn","Vim",
                           "Mki67", "Top2a"),pt.size=0)


dev.off()


pdf(file="05.PDAC-cluster2.pdf",width=20,height=40)
VlnPlot(pbmc, features = c("COL1A1","COL1A2","MCAM","CSPG4","PDGFRB","NES","RGS5","CBF1","TLR1", "TLR10", "TLR2", "TLR6","TLR7",
                           "ACTA2","AMY2A", "AMY2B","KRT19","EPCAM", "CDH1", "CLDN3","CLDN18", "MUC5AC","MUC4","TM4SF1",
                           "PTPRC", "ITGAM", "CCR3","ENPP3","KIT","CD14","FCRL4", "CCR1","NOTCH2","JAM3",
                           "FUT4","CD83", "CD86", "CD80","CD163", "MRC1","CD14","CD8A","CD4",
                           "NCAM1","CD19","CD79A","CD79B","CD27","EIF4EBP1","MS4A1","LRRC15","ALP","STRO-1","Alkaline phosphatase",
                           "CD80","CD84","CD86","SDC1","TNFRSF17","IRF4","PRDM1","XBP1","FKBP11","TNFRSF13C",
                           "RUNX2","SOX9","ALPL","CTSK","VWF","PECAM1","VIM","MYOD1","MYOG",
                           "OGN",),pt.size=0)
VlnPlot(pbmc, features = c("Col1a1","Col1a2","Acta2","Amy2a", "Amy2b","Krt19", "Cldn18", "Muc5ac",
                           "Ptprc", "Itgam", "Ccr3","Enpp3","Kit","Cd14",
                           "Fut4","Cd83", "Cd86", "Cd80","Cd163", "Mrc1","Cd14","Cd8a","Cd4",
                           "Ncam1","Cd19","Sdc1","Tnfrsf17","Mcam","Cspg4","Pdgfrb",
                           "Nes","RgS5"),pt.size=0)
VlnPlot(pbmc, features = c("Ly6c1","Pi16","Serpine2","Gpx3","Tagln","Cthrc1","Tmem100","Dpp4",
                           "Smoc2","Col15a1","Crabp1","Crabp2","Upk3b","Lrrn4","Krt8","Krt18",
                           "Cxcl1","Has1","Top2a","Mki67","C1qb","Csf1r","Pnlip"),pt.size=0)
dev.off()
pdf(file="05.PDAC-cluster-2.pdf",width=15,height=12)
VlnPlot(pbmc, features = c("COL1A1","ACTA2","MUC4","TM4SF1", "MUC5B","CRISP3"),pt.size=0)
dev.off()

#CD45=PTPRC,CD11B=ITGAM,CD193=CCR3,CD203=ENPP3,CD117=KIT,CD15=FUT4,CD206=MRC1
#CD56=NCAM1 ,CD138=SDC1,BCMA=TNFRSF17
pdf(file="05.ARP.features.pdf",width=10,height=8)
VlnPlot(pbmc, features = c("Lama5","Thbs2","Lama3","Cd44","Tnc","Thbs1","Sdc1",
                           "Npnt","Lamc2","Itgb6","Itga2","Itga3","Hmmr","Lamb3"),pt.size=0)
VlnPlot(pbmc, features = c("ACTR2", "ACTR3","ARPC1A", "ARPC1B","ARPC2", "ARPC3","ARPC4", "ARPC5"),pt.size=0)
dev.off()
VlnPlot(pbmc, features = c("Sdc4","Lamb3","Itgb8","Thbs4","Thbs2","Thbs1",
                           "Spp1","Npnt","Lamc2","Tnc","Fn1","Col1a2","Col1a1",
                           "Sdc1","Itgb6","Itgb5","Itgb1","Itgav","Itga5",
                           "Itga4","Itga3","Itga2","Cd44","Lama3","Hmmr",
                           "Col5a2",
                           ),pt.size=0)



write.table(sig.markers,file="06.markers-PDAC.xls",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#????marker?ڸ???cluster????ͼ
pdf(file="06.tsneHeatmap.pdf",width=12,height=9)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()
#??"ACTR2", "ACTR3","ARPC1A", "ARPC1B","ARPC2", "ARPC3","ARPC4", "ARPC5"
#??"Actr2", "Actr3","Arpc1a", "Arpc1b","Arpc2", "Arpc3","Arpc4", "Arpc5"
#????marker??С????ͼ
pdf(file="06.markerViolin-DKK3.pdf",width=10,height=12)
VlnPlot(object = pbmc, features = c("DKK3","ACTR2", "ACTR3","ARPC1A", "ARPC1B","ARPC2", "ARPC3","ARPC4", "ARPC5"))
dev.off()
pdf(file="06.markerViolin-ACRP30.pdf",width=10,height=12)
VlnPlot(object = pbmc, features = c("APN","ACTR2", "ACTR3","ARPC1A", "ARPC1B","ARPC2", "ARPC3","ARPC4", "ARPC5"))
dev.off()

#????marker?ڸ???cluster??ɢ??ͼ
pdf(file="06.markerScatter-IGF1-2.pdf",width=12.5,height=20)
FeaturePlot(object = pbmc, features = c("DCN","PDPN","FAP","LUM","PDGFRA","C7","MGP","COL12A1","COL11A1","LOXL2","COL14A1","SDC1","NFIA",
                                        "SPARCL1","OGN","LGALS1","C5orf46","ANTXR1","C3","ABI3BP","MMP11","SLC16A3","IGFBP4","COL5A2","CTHRC1",
                                        "IGF1","IL6","ACTA2","TPM1"),
            cols = c("darkgray", "red","red"))
dev.off()

pdf(file="06.CAF.markerScatter-ITGA11.pdf",width=6.5,height=5)
FeaturePlot(object = pbmc, features = c("ITGA11"), cols = c( "gray91","red"))
dev.off()

#????marker?ڸ???cluster??????ͼ
DotPlot(pbmc, group.by = 'seurat_clusters',
        features = unique(pbmc$GENE_SYMBOL)) + RotatedAxis()+ scale_color_gradient2(low = "blue", mid = "white", high = "red")+ coord_flip()+ theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
library(ggplot2)
pdf(file="06.markerBubble-ECM-2.pdf",width=16,height=4)
cluster10Marker=c("CD19",	"IL2RA", "TNFRSF8","IgG", "CD27", "CD38", "CD78", "SDC1", "SLAMF7","IL6",
                  "IgA", "IgE", "MS4A1",  "CD40", "CD80", "PDCD1LG2","CXCR3", "CXCR4", "CXCR5", 
                  "CXCR6","CD1A","CD1B", "NOTCH2","IgD", "CR2", "CD22", "FCER2",
                  "CD5", "CD24", "TLR4","IL10", "TGFB1","ITGAM","SDC1","LY6C1")
cluster10Marker=c("Col5a2","Hmmr","Lama3","Cd44","Itga2","Itga3","Itga4","Itga5","Itgav",
                  "Itgb1","Itgb5","Itgb6","Sdc1","Col1a1","Col1a2","Fn1","Tnc","Lamc2","Npnt",
                  "Spp1","Thbs1","Thbs2","Thbs4","Itgb8","Lamb3","Sdc4")
cluster10Marker=c("COL5A2","HMMR","LAMA3","CD44","ITGA2","ITGA3","ITGA4","ITGA5","ITGAV",
                  "ITGB1","ITGB5","ITGB6","SDC1","COL1A1","COL1A2","FN1","TNC","LAMC2","NPNT",
                  "SPP1","THBS1","THBS2","THBS4","ITGB8","LAMB3","SDC4")
cluster10Marker=c("KIT","IL7R","ATXN1","FLT3","PTPRC","SPN","RAG1","RAG2","EBF1","PAX5","SOX4",
                  "LEF1","CD19","RAG7","CD24","EIF4EBP1","TNFRSF13C","TNFRSF17",
                  "MS4A1","CD79A","CD79B","CXCR4")
cluster10Marker=c("ACTA2","TAGLN","MMP11","MYL9","HOPX","POSTN","TPM1","TPM2",
                             "PDGFRA","CXCL12","CFD","DPT","C3","CCL2","CXCL1","COL14A1","AGTR1","HAS1","C7","IGF1","IL6")
cluster10Marker=c("DCN","TAGLN","PDGFRB","COL3A1","ACTA2",
                  "OS9","FBN2","GLI1","MBD6","PTCH1",
                  "SOX9","COL2A1","SOX6","COL9A1","MIA",
                  "SFRP2","NPW","HSPA2","HSPA6","SOHLH1",
                  "C1QA","HLA-DRA","CD14","CD74","HLA-DRB1","CD163",
                  "EPCAM")
cluster10Marker=c("PTPRC","ITGAM","CD14","CD68","HLA-DRA","CD80","CD86","NOS2","CD163","MRC1",
                  "IL6","IL10","IL12","TNFA")
DotPlot(object = pbmc, features = cluster10Marker,dot.scale = 15,col.min = 0,
        col.max = 2,)+ RotatedAxis()+scale_color_gradient2(low = "white", mid = "pink", high = "red")#cols = c("blue","red"),
DotPlot(object = pbmc, features = cluster10Marker,dot.scale = 7,col.min = -2,
        col.max = 2,)+ RotatedAxis()+scale_color_gradient2(low = "blue", mid = "white", high = "red")+ coord_flip()+ theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))



#DotPlot(object = pbmc, features = cluster10Marker,dot.scale = 10)#col.min = -2,col.max = 2,cols = c("grey","red"),
dev.off()
###################################07.ע??ϸ??????###################################
#????1--ͨ?????ǻ????????ף??????˹?ȷ????????Ⱥ??ϸ?????ͣ????????????ֶ?????ϸ??Ⱥ????
bfreaname.pbmc <- pbmc
#cellmarker??վ????http://biocc.hrbmu.edu.cn/CellMarker/???ֶ????????????ϱ?
#????VlnPlot(pbmc, features = top10$gene[1:20],pt.size=0)?????Ĳ???????????
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#????2-SingleR???????????????Ǵ˷?ע?Ͳ?׼ȷ????Ҫ?ο????ݼ?
library(SingleR)
counts<-pbmc@assays$RNA@counts
#write.table(counts,file="07.monoclecounts-PDAC.txt",sep='\t',quote=F,row.names=T)

clusters<-pbmc@meta.data$seurat_clusters
ann=pbmc@meta.data$orig.ident
singler = CreateSinglerObject(counts, annot = ann, "pbmc", min.genes = 0,
  species = "Mouse", citation = "",#??ΪHuman??С??ΪMouse
  ref.list = list(), normalize.gene.length = F, variable.genes = "de",
  fine.tune = F, do.signatures = T, clusters = clusters, do.main.types = T,
  reduce.file.size = T, numCores = 1)
singler$seurat = pbmc
singler$meta.data$xy = pbmc@reductions$tsne@cell.embeddings
clusterAnn=singler$singler[[2]]$SingleR.clusters.main$labels
write.table(clusterAnn,file="07.clusterAnn-PDAC.txt",quote=F,sep="\t",col.names=F)
write.table(singler$other,file="07.cellAnn-PDAC.txt",quote=F,sep="\t",col.names=F)
#?????к??ʵĲο????ݼ?ʱ?????÷??????????Ľ???ע??
#????3???Զ???singleR??ע??

#????4??????seurat???õ?ԭ??????ϸ?????ϵĹ??ܣ????ο?????????ע?????ݽ???ӳ?䴦??

#׼??monocle??????Ҫ???ļ?
monocle.matrix=as.matrix(pbmc@assays$RNA@data)
monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
write.table(monocle.matrix,file="07.monocleMatrix-PDAC.txt",quote=F,sep="\t",row.names=F)
monocle.sample=as.matrix(pbmc@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="07.monocleSample-PDAC.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="07.monocleGene-PDAC.txt",quote=F,sep="\t",row.names=F)
write.table(singler$other,file="07.monocleClusterAnn-PDAC.txt",quote=F,sep="\t",col.names=F)
write.table(sig.markers,file="07.monocleMarkers-PDAC-4.txt",sep="\t",row.names=F,quote=F)
########################monocle????
library(monocle)
monocle.matrix=read.table("07.monocleMatrix-PDAC.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.sample=read.table("07.monocleSample-PDAC.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.geneAnn=read.table("07.monocleGene-PDAC.txt",sep="\t",header=T,row.names=1,check.names=F)
marker=read.table("07.monocleMarkers-PDAC-4.txt",sep="\t",header=T,check.names=F)

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
plot_cell_trajectory(cds, color_by = "Cluster",show_branch_points = FALSE)+  
  scale_colour_manual(
    values =c("indianred1","yellow4","mediumseagreen","dodgerblue","darkorchid1")
    # aesthetics = c("colour", "fill")
  )
dev.off()
pdf(file="cellType.trajectory.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cell_type2")
dev.off()
#??λ???Լ?????Ȥ????Ⱥ???????·???
Idents(pbmc)
levels(pbmc)
head(pbmc@meta.data)
pbmc= pbmc[,pbmc@meta.data$seurat_clusters %in% c(3,4,5,7,10,12)]
Bcells = pbmc[,pbmc@meta.data$seurat_clusters %in% c(17)]
Fibroblast<- merge(x = pbmc, y = pbmc1,
                   add.cell.ids=c("pbmc","pbmc1"),merge.data=TRUE, project="10X_HSC")
pbmc=Fibroblast
#ʣ?µ??߱?׼???뼴??,?????ݱ?׼????ʼ??Ȼ?󷵻????????з?????ע???޸?·??
pbmc <- NormalizeData(object = Fibroblast1, Normalization.method = "LogPDACize", scale.factor = 10000)
pbmc <- NormalizeData(object = Fibroblast1, scale.factor = 10000)
dim(pbmc)
save(pbmc,file='fibroblast.Rdata')
saveRDS(Fibroblast, "all.rds")

##########################????С????ͼ########################
#install.packages("vioplot") 
#???Ƚ???????????ΪArp complex.txt??ʽ
library(vioplot)
data=read.table("Arp complex Epithelial cells.txt", header = T, sep = "\t", row.names = 1, check.names = F)
data =data.frame(data)
table(data$name_clusters)
range(data$StromalScore)
##??????Ϣ????
group<-data$name_clusters
length(group)==dim(data)[[1]]##ȷ????Ϣƥ??

data$group<-group
data=data[order(data$group),]#?????ݽ???????
table(data$group)
library(ggpubr)
my_comparisons <- list(c("G1","G2"), c("G2", "G3"), c("G3", "G4"),c("G4", "GX"))##?????趨
my_comparisons <- list(c( "Dead","Alive"))##?????趨
e<-data %>% 
  #dplyr::filter(group %in% c("i", "ii","iii", "iv")) %>% 
  ggviolin(x = "group", y = c(colnames(data)[3:10]), fill = "group",#ѡ????1??3ԭ??Ϊ????????ֵ??С
           combine = T,
           #palette = c("#00AFBB", "#E7B800", "#FC4E07","#E7B800"),##
           ylab="PDACized Expression",
           add = "boxplot", add.params = list(fill = "white"))
e+stat_compare_means(method = "t.test",
                     #label = "p.signif",##?Ǻ?????
                     comparisons = my_comparisons)
ggsave(file = "arp violin.pdf", width = 10, height = 8.5)




##########??????????
setwd("C:/Users/Administrator/Desktop/Analysis/scRNA/PDAC and PDAC")  #?趨?Լ??Ĺ???Ŀ¼            #???ù???Ŀ¼
library(Seurat)
#??ȡ?????걾??Ϣ
samples=list.files("PDAC/")
samples
nm_e<- file.path('./PDAC',samples)
names(nm_e) <- samples
#?ϲ?????
nm_e<- Read10X(data.dir = nm_e)
row.names(nm_e) <- nm_e[,1]
nm_e <- nm_e[,-1]
nm_e[1:4,1:4]
#data = CreateSeuratObject(counts, min.cells=1)
dim(nm_e)   #?鿴????????ϸ??????
#??ȡ???ٰ??걾??Ϣ
samples=list.files("PDAC/")
samples
tm_e<- file.path('./PDAC',samples)
names(tm_e) <- samples
#?ϲ?????
tm_e<- Read10X(data.dir = tm_e)
#data = CreateSeuratObject(counts, min.cells=1)
row.names(tm_e) <- tm_e[,1]
tm_e <- tm_e[,-1]
tm_e[1:4,1:4]
dim(tm_e)   #?鿴????????ϸ??????
#?????ݺϲ?
test <- cbind(nm_e,tm_e)
test <- as.matrix(test)
#ԭ��?Ļ?????̫??????ȡ???е?symbol??ʽ
rownames(test) <- sapply(strsplit(rownames(test),"_"),"[",2)
test[1:4,1:4]
#??????Ϣ
group_dat <- data.frame(group=c(rep('PDAC',ncol(nm_e)),
                                    rep('tumor',ncol(tm_e))))
rownames(group_dat) <- colnames(test)
#????Seurat????
library(Seurat)
scRNA <- CreateSeuratObject(counts=test, 
                            meta.data=group_dat)
dim(scRNA) #????432????????57241??????
table(scRNA@meta.data$group)
#?ʿ?
minGene=500
maxGene=4000
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
pctMT=30
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
dim(scRNA) #??ʣ57??????
table(scRNA@meta.data$group)




