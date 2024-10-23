library(Seurat)
library(tidyverse)
library(harmony)

setwd("C:/Users/lukas/Documents/CQ_Praktikum")
load("Lucas_noXCR.integrated_SNN_rpca_50.RData")


Disc_colors = c("#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72","#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",

                "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#E6AB02","#7570B3","#BEAED4","#666666","#999999",

                "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000","#aeae5c","#1e90ff","#00bfff","#56ff0d","#ffff00",



                "#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72","#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",

                "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#E6AB02","#7570B3","#BEAED4","#666666","#999999",

                "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000","#aeae5c","#1e90ff","#00bfff","#56ff0d","#ffff00"

)
use_these_colors = c(Disc_colors[1:length(Idents(SO.har))])


# # adding the individual replicates to meta data
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

# # getting rid of contaminations and dying cells
cleaning <- read_csv("w_o_contamination.csv")
cleaned <- cleaning$Barcode
SO.integrated <- subset(SO.integrated, cells = cleaned)

# # getting rid of the old experiments
idx <- which((SO.integrated@meta.data$experiment != 'Tcell_vehicle') &
               (SO.integrated@meta.data$experiment != 'Tcells_treated'))
SO.integrated <- SO.integrated[,idx]

# # subsetting the SO to the barcodes where I have metadata for Replicates

replicates <- read.csv("Lukas_Analysis.csv")
repl_subset <- replicates$Barcode
SO.integrated <- subset(SO.integrated, cells = repl_subset)

# adding new meta data with organ
SO.integrated$organ <- ifelse(grepl("BM", SO.integrated$experiment), "BM",
                              ifelse(grepl("Spleen", SO.integrated$experiment), "Spleen", NA))

# Subset for BM samples
SO.BM <- subset(SO.integrated, subset = organ == "BM")

# Subset for Spleen samples
SO.Spleen <- subset(SO.integrated, subset = organ == "Spleen")


SaveSeuratRds(SO.Spleen, file = "./SO.Spleen.RDS")
SaveSeuratRds(SO.BM, file = "./SO.BM.RDS")

