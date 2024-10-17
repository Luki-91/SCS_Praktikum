library(Seurat)
library(tidyverse)
library(harmony)

setwd("C:/Users/lukas/Documents/CQ_Praktikum")
load("Lucas_noXCR.integrated_SNN_rpca_50.RData")

# adding the individual replicates to meta data
replicates <- read.csv("Lukas_Analysis.csv")
rownames(replicates) <- replicates$Barcode
repl_subset <- replicates$Barcode
replicates$Barcode <- NULL

SO.integrated <- AddMetaData(SO.integrated, replicates)

# adding the CITE Seq as a seperate assay
SO_cite = Read10X_h5("filtered_feature_bc_matrix.h5")
SO_cite = CreateSeuratObject(counts = SO_cite[["Antibody Capture"]])
SO_cite = NormalizeData(SO_cite, normalization.method = "CLR")
SO.integrated@assays$Cite = SO_cite@assays$RNA
SO.integrated@assays[["Cite"]]@key
SO.integrated@assays[["Cite"]]@key <- "cite_"

# getting rid of contaminations and dying cells
cleaning <- read_csv("w_o_contamination.csv")
cleaned <- cleaning$Barcode
SO.integrated <- subset(SO.integrated, cells = cleaned)

# getting rid of the old experiments
idx <- which((SO.integrated@meta.data$experiment != 'Tcell_vehicle') &
               (SO.integrated@meta.data$experiment != 'Tcells_treated'))
SO.integrated <- SO.integrated[,idx]

# subsetting the SO to the barcodes where I have metadata for Replicates

replicates <- read.csv("Lukas_Analysis.csv")
repl_subset <- replicates$Barcode
SO.integrated <- subset(SO.integrated, cells = repl_subset)

# subsetting the SO to only CD4+ T cells
unchanged <- read_csv("CD8_reclustering.csv")
filter <- unchanged$Barcode

SO.integrated <- subset(SO.integrated, cells = filter)

DefaultAssay(SO.integrated) <- "RNA"

DimPlot(SO.integrated,
        #group.by = "RNA_snn_res.0.6",
        #split.by = "Replicates",
        label.size = 10, pt.size = 2,label = T, ncol = 2)

# Run Harmony
SO.har <- NormalizeData(SO.integrated)
SO.har <- FindVariableFeatures(SO.har)
SO.har <- ScaleData(SO.har)
SO.har <- RunPCA(SO.har)
ElbowPlot(SO.har, ndims = 30, reduction = "pca")
SO.har <- RunHarmony(SO.har, group.by.vars = "Replicates", plot_convergence = TRUE)
SO.har <- RunUMAP(SO.har, reduction = "harmony", dims = 1:20)
SO.har <- FindNeighbors(SO.har, reduction = "harmony", dims =1:20)
SO.har <- FindClusters(SO.har, resolution =  c(0.4,0.5,0.6))
rm("SO.integrated")

SO.har@meta.data <- SO.har@meta.data %>%
  mutate(Condition = case_when(
    grepl("vehicle", experiment) ~ "vehicle",
    grepl("treated", experiment) ~ "treated",
    TRUE ~ "Unknown"  # Default if none of the patterns match
  ))

DimPlot(SO.har, group.by = "RNA_snn_res.0.6",
        split.by = "Condition",
        label.size = 10, pt.size = .5,label = T, ncol = 2)

features  <-  c("Cd4","Cd8a","Isg20","Ms.CD69","Foxp3","Ms.CD4","Ms.CD8a")

for (marker in features) {
  p <- FeaturePlot(SO.har, features = marker,
            label.size = 10, pt.size = .5,label = F
            #, split.by = "Replicates",ncol=2
            )
  print(p)
}
## filterring for cells, which do not have any of these barcodes
# unchanged <- read_csv("subset_filter.csv")
# filter <- unchanged$Barcode
#
# subsetted <- subset(SO.har, cells = filter)
# subsetted <- ScaleData(subsetted)
# subsetted <- NormalizeData(subsetted)
# subsetted <- FindVariableFeatures(subsetted)
# subsetted <- RunHarmony(subsetted, "experiment", plot_convergence = TRUE)
# subsetted <- FindNeighbors(subsetted, reduction="harmony", dims = 1:20)
# subsetted <- FindClusters(subsetted, resolution = 0.4)
# subsetted <- RunUMAP(subsetted, reduction="harmony", dims = 1:20)
#
#
# DimPlot(subsetted, reduction = "umap",
#         group.by = "seurat_clusters", split.by = "experiment", ncol = 2,
#         label.size = 10, pt.size = 2,label = T)

## annotating the dataset ##

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah,
              pattern = c("Mus musculus", "EnsDb"),
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb,
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

write.csv2(annotations, file ="annotations.csv")

# setting res of 0.4 to the active Ident
Idents(SO.har) <- SO.har@meta.data$RNA_snn_res.0.6
cluster0_con <- FindConservedMarkers(SO.har, ident.1 = 0, only.pos=T,
                                     grouping.var = "Condition",
                                     logfc.threshold = 0.25)
write.csv2(cluster0_con, file ="cluster0_con.csv")
cluster1_con <- FindConservedMarkers(SO.har, ident.1 = 1, only.pos=T,
                                     grouping.var = "Condition",
                                     logfc.threshold = 0.25)
write.csv2(cluster1_con, file ="cluster1_con.csv")
cluster2_con <- FindConservedMarkers(SO.har, ident.1 = 2, only.pos=T,
                                     grouping.var = "Condition",
                                     logfc.threshold = 0.25)
write.csv2(cluster2_con, file ="cluster2_con.csv")
cluster3_con <- FindConservedMarkers(SO.har, ident.1 = 3, only.pos=T,
                                     grouping.var = "Condition",
                                     logfc.threshold = 0.25)
write.csv2(cluster3_con, file ="cluster3_con.csv")
cluster4_con <- FindConservedMarkers(SO.har, ident.1 = 4, only.pos=T,
                                     grouping.var = "Condition",
                                     logfc.threshold = 0.25)
write.csv2(cluster4_con, file ="cluster4_con.csv")
cluster5_con <- FindConservedMarkers(SO.har, ident.1 = 5, only.pos=T,
                                     grouping.var = "Condition",
                                     logfc.threshold = 0.25)
write.csv2(cluster5_con, file ="cluster5_con.csv")
cluster6_con <- FindConservedMarkers(SO.har, ident.1 = 6, only.pos=T,
                                     grouping.var = "Condition",
                                     logfc.threshold = 0.25)
write.csv2(cluster6_con, file ="cluster6_con.csv")
cluster7_con <- FindConservedMarkers(SO.har, ident.1 = 7, only.pos=T,
                                     grouping.var = "Condition",
                                     logfc.threshold = 0.25)
write.csv2(cluster7_con, file ="cluster7_con.csv")
cluster8_con <- FindConservedMarkers(SO.har, ident.1 = 8, only.pos=T,
                                     grouping.var = "Condition",
                                     logfc.threshold = 0.25)
write.csv2(cluster8_con, file ="cluster8_con.csv")
cluster9_con <- FindConservedMarkers(SO.har, ident.1 = 9, only.pos=T,
                                     grouping.var = "Condition",
                                     logfc.threshold = 0.25)
write.csv2(cluster9_con, file ="cluster9_con.csv")



# cluster0_ann_markers <- cluster0_con %>%
#   rownames_to_column(var="gene") %>%
#   left_join(y = unique(annotations[, c("gene_name", "gene_biotype")]),
#             by = c("gene" = "gene_name"))
#
# View(cluster0_ann_markers)


# Adding a new metadata column with condition+cluster
SO.har$comparison <- paste(SO.har$Condition, sep = "_", Idents(SO.har))
# making this the new Ident
Idents(SO.har) <- "comparison"

cluster0 <- FindMarkers(SO.har, ident.1 = 'vehicle_0', ident.2 = 'treated_0')
write.csv2(cluster0, file ="cluster0.csv")
cluster1 <- FindMarkers(SO.har, ident.1 = 'vehicle_1', ident.2 = 'treated_1')
write.csv2(cluster1, file ="cluster1.csv")
cluster2 <- FindMarkers(SO.har, ident.1 = 'vehicle_2', ident.2 = 'treated_2')
write.csv2(cluster2, file ="cluster2.csv")
cluster3 <- FindMarkers(SO.har, ident.1 = 'vehicle_3', ident.2 = 'treated_3')
write.csv2(cluster3, file ="cluster3.csv")
cluster4 <- FindMarkers(SO.har, ident.1 = 'vehicle_4', ident.2 = 'treated_4')
write.csv2(cluster4, file ="cluster4.csv")
cluster5 <- FindMarkers(SO.har, ident.1 = 'vehicle_5', ident.2 = 'treated_5')
write.csv2(cluster5, file ="cluster5.csv")
cluster6 <- FindMarkers(SO.har, ident.1 = 'vehicle_6', ident.2 = 'treated_6')
write.csv2(cluster6, file ="cluster6.csv")
cluster7 <- FindMarkers(SO.har, ident.1 = 'vehicle_7', ident.2 = 'treated_7')
write.csv2(cluster7, file ="cluster7.csv")
cluster8 <- FindMarkers(SO.har, ident.1 = 'vehicle_8', ident.2 = 'treated_8')
write.csv2(cluster8, file ="cluster8.csv")
cluster9 <- FindMarkers(SO.har, ident.1 = 'vehicle_9', ident.2 = 'treated_9')
write.csv2(cluster9, file ="cluster9.csv")


Idents(SO.har) <- SO.har@meta.data$RNA_snn_res.0.6
SO.har <- RenameIdents(SO.har, '0'='CD8+ Tcm', '1'='CD4+CD154+', '2'='Treg',
                       '3'='CD8+ Tcm', '4'='CD8+ Trm', '5'='Th17-like',
                       '6'='Nkg7+ CD8+','7'='CXCR6+ CD4+','8'='Apoptotic cells+',
                       '9'='Tcf7+ CD4+','10'='Tfr','11'='IFN I stiumlated',
                       '12'='Mki67+')

DimPlot(SO.har, reduction = "umap", split.by = "Condition",
        ncol = 2, label.size = 10, pt.size = 2,label = F)

FeaturePlot(SO.har, features = c("Cd4","Cd8a","Foxp3","Isg15","Ms.CD69"),
            label.size = 10, pt.size = 2,label = F,
            split.by = "Condition")

FeaturePlot(SO.har, features = c("Gzmm","Gzmk","Izumo1r","Sostdc1","Tmem176a"),
            label.size = 10, pt.size = 2,label = F,
            split.by = "Condition")

SaveSeuratRds(SO.har, file = "./SO.CD8+.RDS")

FeaturePlot(SO.har,
            features = c("mmHashtag5","mmHashtag6","mmHashtag7","mmHashtag8"),
            label.size = 10, pt.size = 2,label = F,
            split.by = "experiment")

