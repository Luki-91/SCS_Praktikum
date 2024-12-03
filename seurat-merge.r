library(Seurat)
library(tidyverse)

# For memory T cells - replace these with your actual T cell subtypes
t_cell_types <- c(
    "Memory_CD4_Treg",
    "Memory_CD4_Tcm",
    "Memory_CD4_Th17",
    "Memory_CD4_Trm",
    "Memory_CD8_Tcm",
    "Memory_CD8_Trm",
    "Memory_ID2_effector",
    "Memory_IFN_I",
    "Memory_Hobit_Trm",
    "Apoptotic",
    "Proliferating"
)

# For stroma cells - replace these with your actual stromal subtypes
stroma_types <- c(
    "ST2",
    "classical_stroma",
    "special_stroma",
    "contamination"
)

# Function to assign cell types based on your clustering or marker expression
assign_t_cell_types <- function(seurat_obj) {
    # Replace this with your actual classification logic
    # This could be based on marker genes, clustering results, or existing metadata

    # Example: if using clustering results
    cell_types <- case_when(
        SO.BM$RNA_snn_res.0.6 == "0" ~ "Memory_CD8_Tcm",
        SO.BM$RNA_snn_res.0.6 == "1" ~ "Memory_CD8_Trm",
        SO.BM$RNA_snn_res.0.6 == "2" ~ "Memory_CD8_Tcm",
        SO.BM$RNA_snn_res.0.6 == "3" ~ "Memory_CD4_Treg",
        SO.BM$RNA_snn_res.0.6 == "4" ~ "Memory_CD4_Tcm",
        SO.BM$RNA_snn_res.0.6 == "5" ~ "Memory_ID2_effector",
        SO.BM$RNA_snn_res.0.6 == "6" ~ "Memory_CD4_Th17",
        SO.BM$RNA_snn_res.0.6 == "7" ~ "Apoptotic",
        SO.BM$RNA_snn_res.0.6 == "8" ~ "Memory_Hobit_Trm",
        SO.BM$RNA_snn_res.0.6 == "9" ~ "Apoptotic",
        SO.BM$RNA_snn_res.0.6 == "10" ~ "Memory_IFN_I",
        SO.BM$RNA_snn_res.0.6 == "11" ~ "Memory_IFN_I",
        SO.BM$RNA_snn_res.0.6 == "12" ~ "Proliferating",)

    return(cell_types)
}

assign_stroma_types <- function(seurat_obj) {
    # Replace this with your actual classification logic for stroma cells
  cell_types <- case_when(
    SO$RNA_snn_res.0.3 == "0" ~ "classical_stroma",
    SO$RNA_snn_res.0.6 == "1" ~ "contamination",
    SO$RNA_snn_res.0.6 == "2" ~ "special_stroma",
    SO$RNA_snn_res.0.6 == "3" ~ "ST2",
    SO$RNA_snn_res.0.6 == "4" ~ "special_stroma",
    SO$RNA_snn_res.0.6 == "5" ~ "contamination",
    SO$RNA_snn_res.0.6 == "6" ~ "special_stroma",
    SO$RNA_snn_res.0.6 == "7" ~ "special_stroma",
    SO$RNA_snn_res.0.6 == "8" ~ "ST2",
    SO$RNA_snn_res.0.6 == "9" ~ "contamination",
    SO$RNA_snn_res.0.6 == "10" ~ "contamination",
    SO$RNA_snn_res.0.6 == "11" ~ "special_stroma",
    SO$RNA_snn_res.0.6 == "12" ~ "ST2",)

    return(cell_types)
}

# Assign cell types to each object
SO.BM$cell_type <- assign_t_cell_types(SO.BM)
SO$cell_type <- assign_stroma_types(SO)

# Merge the objects
combined_seurat <- merge(SO.BM,
                        y = SO,
                        add.cell.ids = c("T", "Stroma"),
                        project = "combined_analysis")

# Verify the merge and cell type assignments
print("Cell type distribution:")
print(table(combined_seurat$cell_type))

# Process the combined dataset
combined_seurat <- NormalizeData(combined_seurat)
combined_seurat <- FindVariableFeatures(combined_seurat)
combined_seurat <- ScaleData(combined_seurat)
combined_seurat <- RunPCA(combined_seurat)
ElbowPlot(SO, ndims = 30, reduction = "pca")

combined_seurat <- RunUMAP(combined_seurat, dims = 1:30)
combined_seurat <- FindNeighbors(SO, reduction = "harmony", dims =1:20)
combined_seurat <- FindClusters(SO, resolution =  c(0.1,0.2,0.3,0.4,0.5,0.6))

# Visualize the combined data with all subtypes
p1 <- DimPlot(combined_seurat,
              group.by = "cell_type",
              label = TRUE,
              repel = TRUE) +
      theme(legend.position = "right") +
      ggtitle("Combined Cell Types")

# Create a simplified grouping for broad categories
combined_seurat$broad_type <- ifelse(
    grepl("^Memory", combined_seurat$cell_type),
    "T cells",
    "Stromal cells"
)

# Visualize with broad grouping
p2 <- DimPlot(combined_seurat,
              group.by = "broad_type",
              label = TRUE,
              repel = TRUE) +
      ggtitle("Broad Cell Types")

p1 + p2

# Save the combined object
saveRDS(combined_seurat, "combined_seurat_object.rds")
