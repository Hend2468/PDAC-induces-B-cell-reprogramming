

library(Seurat)
library(readxl)
# Read your marker table
marker_table <- read_excel("/Volumes/NatalyDrive/Hend_and_Ali_analysis/Bcell\ states.xlsx")

# Create a list of markers for each cell type
marker_list <- split(marker_table$Gene, marker_table$Cell.type)

All_PDAC_Bcell <- JoinLayers(object = All_PDAC_Bcell) #### to be able to addmodulescore

# Add module scores to the Seurat object based on marker genes
for(cell_typeB in names(marker_list)) {
  All_PDAC_Bcell <- AddModuleScore(
    object = All_PDAC_Bcell,
    features = list(marker_list[[cell_typeB]]),
    name = paste0(cell_typeB, "_score")
  )
}

# Identify the cell type for each cluster based on the highest module score
All_PDAC_Bcell$cell_typeB <- apply(
  All_PDAC_Bcell@meta.data[, grep("_score1$", colnames(All_PDAC_Bcell@meta.data))],
  1,
  function(scores) names(scores)[which.max(scores)]
)

# Rename clusters based on identified cell types
Idents(All_PDAC_Bcell) <- All_PDAC_Bcell$cell_typeB
Idents(All_PDAC_Bcell) <- All_PDAC_Bcell$seurat_clusters

DimPlot(All_PDAC_Bcell, reduction = "umap", label = F, pt.size = 0.5)

dittoBarPlot(All_PDAC_Bcell, var = "cell_typeB", group.by = "seurat_clusters")
dittoBarPlot(All_PDAC_Bcell, var = "cell_typeB", group.by = "cell_type")

