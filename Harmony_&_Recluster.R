library(Seurat)
library(tidyverse)
library(harmony)

load("Lucas_noXCR.integrated_SNN_rpca_50.RData")

# getting rid of contaminations and dying cells
cleaning <- read_csv("w_o_contamination+dying_cells.csv")
cleaned <- cleaning$Barcode
SO.integrated <- subset(SO.integrated, cells = cleaned)

# getting rid of the old experiments
idx <- which((SO.integrated@meta.data$experiment != 'Tcell_vehicle') &
               (SO.integrated@meta.data$experiment != 'Tcells_treated'))
SO.integrated <- SO.integrated[,idx]

DefaultAssay(SO.integrated) <- "RNA"
SO.integrated <- NormalizeData(SO.integrated)
SO.integrated <- FindVariableFeatures(SO.integrated)
SO.integrated <- ScaleData(SO.integrated)
SO.integrated <-RunPCA(SO.integrated)
ElbowPlot(SO.integrated, ndims = 30, reduction = "pca")
SO.integrated <- RunUMAP(SO.integrated, dims = 1:30)
SO.integrated <- FindNeighbors(SO.integrated, dims = 1:30)

# Run Harmony
SO.har <- NormalizeData(SO.integrated)
SO.har <- FindVariableFeatures(SO.har)
SO.har <- ScaleData(SO.har)
SO.har <- RunPCA(SO.har)
ElbowPlot(SO.har, ndims = 30, reduction = "pca")
SO.har <- RunHarmony(SO.har, group.by.vars = "experiment", plot_convergence = TRUE)
SO.har <- RunUMAP(SO.har, reduction = "harmony", dims = 1:30)
SO.har <- FindNeighbors(SO.har, reduction = "harmony", dims =1:30)
SO.har <- FindClusters(SO.har, resolution =  0.4)
DimPlot(SO.har, group.by = c("experiment","seurat_clusters"))
FeaturePlot(SO.har, features = c("Cd4","Cd8a","Isg15","Zfp683","Cd69"),
            cols = c("blue","red"), split.by = "experiment")

SO_cite = Read10X_h5("filtered_feature_bc_matrix.h5")
SO_cite = CreateSeuratObject(counts = SO_cite[["Antibody Capture"]])
SO_cite = NormalizeData(SO_cite, normalization.method = "CLR")
SO.har@assays$Cite = SO_cite@assays$RNA
SO.har@assays[["Cite"]]@key
SO.har@assays[["Cite"]]@key <- "cite_"

## filterring for cells, which do not have any of these barcodes
unchanged <- read_csv("subset_filter.csv")
filter <- unchanged$Barcode

subsetted <- subset(SO.har, cells = filter)
subsetted <- ScaleData(subsetted)
subsetted <- NormalizeData(subsetted)
subsetted <- FindVariableFeatures(subsetted)
subsetted <- RunHarmony(subsetted, "experiment", plot_convergence = TRUE)
subsetted <- FindNeighbors(subsetted, reduction="harmony", dims = 1:30)
subsetted <- FindClusters(subsetted, resolution = 0.4)
subsetted <- RunUMAP(subsetted, reduction="harmony", dims = 1:30)


DimPlot(subsetted, reduction = "umap",
        group.by = "seurat_clusters", split.by = "experiment", ncol = 2)

subsetted@meta.data <- subsetted@meta.data %>%
  mutate(Condition = case_when(
    grepl("vehicle", experiment) ~ "vehicle",
    grepl("treated", experiment) ~ "treated",
    TRUE ~ "Unknown"  # Default if none of the patterns match
  ))

# setting res of 0.4 to the active Ident
Idents(subsetted) <- subsetted@meta.data$RNA_snn_res.0.4
cluster0_con <- FindConservedMarkers(subsetted, ident.1 = 0, only.pos=T,
                                     grouping.var = "Condition")
write.csv2(cluster0_con, file ="cluster0_con.csv")
cluster1_con <- FindConservedMarkers(subsetted, ident.1 = 1, only.pos=T,
                                     grouping.var = "Condition")
write.csv2(cluster1_con, file ="cluster1_con.csv")
cluster2_con <- FindConservedMarkers(subsetted, ident.1 = 2, only.pos=T,
                                     grouping.var = "Condition")
write.csv2(cluster2_con, file ="cluster2_con.csv")
cluster3_con <- FindConservedMarkers(subsetted, ident.1 = 3, only.pos=T,
                                     grouping.var = "Condition")
write.csv2(cluster3_con, file ="cluster3_con.csv")
cluster4_con <- FindConservedMarkers(subsetted, ident.1 = 4, only.pos=T,
                                     grouping.var = "Condition")
write.csv2(cluster4_con, file ="cluster4_con.csv")
cluster5_con <- FindConservedMarkers(subsetted, ident.1 = 5, only.pos=T,
                                     grouping.var = "Condition")
write.csv2(cluster5_con, file ="cluster5_con.csv")


# Adding a new metadata column with condition+cluster
subsetted$comparison <- paste(subsetted$Condition, sep = "_", Idents(subsetted))
# making this the new Ident
Idents(subsetted) <- "comparison"

cluster0 <- FindMarkers(subsetted, ident.1 = 'vehicle_0', ident.2 = 'treated_0')
write.csv2(cluster0, file ="cluster0.csv")
cluster1 <- FindMarkers(subsetted, ident.1 = 'vehicle_1', ident.2 = 'treated_1')
write.csv2(cluster1, file ="cluster1.csv")
cluster2 <- FindMarkers(subsetted, ident.1 = 'vehicle_2', ident.2 = 'treated_2')
write.csv2(cluster2, file ="cluster2.csv")
cluster3 <- FindMarkers(subsetted, ident.1 = 'vehicle_3', ident.2 = 'treated_3')
write.csv2(cluster3, file ="cluster3.csv")
cluster4 <- FindMarkers(subsetted, ident.1 = 'vehicle_4', ident.2 = 'treated_4')
write.csv2(cluster4, file ="cluster4.csv")
cluster5 <- FindMarkers(subsetted, ident.1 = 'vehicle_5', ident.2 = 'treated_5')
write.csv2(cluster5, file ="cluster5.csv")


subsetted <- RenameIdents(subsetted, '0'='Gzmm CD8+', '1'='Dapl1 CD8+', '2'='Izumo1r CD4',
                          '3'='Th17-like', '4'='Gzmk Effector', '5'='Treg-like')

DimPlot(subsetted, reduction = "umap", split.by = "Condition", ncol = 2, pt.size = 0.8)
FeaturePlot(subsetted, features = c("Cd4","Cd8a","Foxp3","Rorc","Cd69"), split.by = "Condition",pt.size = .8)
