library(Seurat)
library(tidyverse)
library(harmony)
library(scCustomize)

# plot <- FeaturePlot_scCustom(SO.har, features = "Ms.CD69",
#                      label.size = 10, pt.size = .8,label = F) & theme(legend.position = "right")
#
# trm_cells <- CellSelector(plot)
#
# plot <- FeaturePlot_scCustom(SO.har, features = "Ms.CD69",
#                              label.size = 10, pt.size = .8,label = F
#
# ) & theme(legend.position = "right")
#
# trm_cells <- CellSelector(plot)
#
# # Create a new metadata column with all three categories
# SO.har$expression_groups <- NA  # default category
# SO.har$expression_groups[colnames(SO.har) %in% treg_cells] <- "treg"
# SO.har$expression_groups[colnames(SO.har) %in% trm_cells] <- "trm"

FeatureScatter(SO.har,
               feature1 = "Klf2",
               feature2 = "Ms.CD69",
               group.by = "RNA_snn_res.0.4",
               split.by = "Replicates",
               #cols = c("#DC050C","#1965B0","#882E72"),
               slot = "data",
               jitter=T)


plot <- FeatureScatter(SO.har,
                       feature1 = "Klf2",
                       feature2 = "Ms.CD69",
                       group.by = "RNA_snn_res.0.4")

Klf2_neg <- CellSelector(plot)

plot <- FeatureScatter(SO.har,
                       feature1 = "Klf2",
                       feature2 = "Ms.CD69",
                       group.by = "RNA_snn_res.0.4")

Klf2_pos <- CellSelector(plot)
Klf2_double_pos <- CellSelector(plot)
Klf2_double_neg <- CellSelector(plot)

plot <- FeatureScatter(SO.har,
                       feature1 = "Klf2",
                       feature2 = "Ms.CD69",
                       group.by = "RNA_snn_res.0.4")

Klf2_CD69_pos <- CellSelector(plot)

plot <- FeatureScatter(SO.har,
                       feature1 = "Cd69",
                       feature2 = "Ms.CD69",
                       group.by = "RNA_snn_res.0.4")

CD69_double_pos <- CellSelector(plot)

plot <- FeatureScatter(SO.har,
                       feature1 = "Cd69",
                       feature2 = "Ms.CD69",
                       group.by = "RNA_snn_res.0.4")

CD69_cite_pos <- CellSelector(plot)

plot1 <- FeatureScatter(SO.har,
                       feature1 = "Cd69",
                       feature2 = "Ms.CD69",
                       group.by = "RNA_snn_res.0.4")

CD69_neg_cells <- CellSelector(plot1)

# Create a new metadata column with all three categories
SO.har$CD69_expression <- NA  # default category
SO.har$CD69_expression[colnames(SO.har) %in% CD69_double_pos] <- "CD69+_cells"
SO.har$CD69_expression[colnames(SO.har) %in% CD69_cite_pos] <- "CD69+_cite_cells"
SO.har$CD69_expression[colnames(SO.har) %in% CD69_neg_cells] <- "CD69-_cells"

SO.har$Klf2_expression <- NA  # default category
SO.har$Klf2_expression[colnames(SO.har) %in% Klf2_neg] <- "Klf2-_cells"
SO.har$Klf2_expression[colnames(SO.har) %in% Klf2_pos] <- "Klf2+_cells"
SO.har$Klf2_expression[colnames(SO.har) %in% Klf2_double_pos] <- "Klf2+CD69+_cells"
SO.har$Klf2_expression[colnames(SO.har) %in% Klf2_double_neg] <- "Klf2-CD69-_cells"

# creating another metadata column in order to compare the individual groups
Idents(SO.har) <- SO.har@meta.data$Klf2_expression
SO.har$comparison_Klf2 <- paste(SO.har$Condition, sep = "_", Idents(SO.har))
# making this the new Ident
Idents(SO.har) <- "comparison_Klf2"

# finding markers between the different groups
Klf2_positive_vehicle <-
  FindMarkers(SO.har, ident.1 = 'vehicle_Klf2+_cells', ident.2 = 'treated_Klf2+_cells',
              test.use="MAST")
write.csv2(Klf2_positive_vehicle, file ="Klf2_positive_vehicle_vs_wort.csv")

Klf2_negative_vehicle <-
  FindMarkers(SO.har, ident.1 = 'vehicle_Klf2-_cells', ident.2 = 'treated_Klf2-_cells',
              test.use="MAST")
write.csv2(Klf2_negative_vehicle, file ="Klf2_negative_vehicle_vs_wort.csv")

Klf2_double_positive<-
  FindMarkers(SO.har, ident.1 = 'vehicle_Klf2+CD69+_cells', ident.2 = 'treated_Klf2+CD69+_cells',
              test.use="MAST")
write.csv2(Klf2_double_positive, file ="Klf2_double_positive_vehicle_vs_wort.csv")

Klf2_double_negative <-
  FindMarkers(SO.har, ident.1 = 'vehicle_Klf2-CD69-_cells', ident.2 = 'treated_Klf2-CD69-_cells',
              test.use="MAST")
write.csv2(Klf2_double_negative, file ="Klf2_double_negative_vehicle_vs_wort.csv")

Klf2_double_negative_vehicle <-
  FindMarkers(SO.har, ident.1 = 'vehicle_Klf2-_cells', ident.2 = 'vehicle_Klf2-CD69-_cells',
              test.use="MAST")
write.csv2(Klf2_double_negative_vehicle, file ="Klf2-_vs_Klf2-CD69-_veh.csv")

Klf2_double_positive_vehicle <-
  FindMarkers(SO.har, ident.1 = 'vehicle_Klf2+_cells', ident.2 = 'vehicle_Klf2-CD69-_cells',
              test.use="MAST")
write.csv2(Klf2_double_positive_vehicle, file ="Klf2+_vs_Klf2-CD69-_veh.csv")

Klf2_double_negative_wort <-
  FindMarkers(SO.har, ident.1 = 'treated_Klf2-_cells', ident.2 = 'treated_Klf2-CD69-_cells',
              test.use="MAST")
write.csv2(Klf2_double_negative_wort, file ="Klf2-_vs_Klf2-CD69-_wort.csv")

Klf2_double_positive_wort <-
  FindMarkers(SO.har, ident.1 = 'treated_Klf2+_cells', ident.2 = 'treated_Klf2-CD69-_cells',
              test.use="MAST")
write.csv2(Klf2_double_positive_wort, file ="Klf2+_vs_Klf2-CD69-_wort.csv")

Klf2_negative_vs_positive_vehicle <-
  FindMarkers(SO.har, ident.1 = 'vehicle_Klf2-_cells', ident.2 = 'vehicle_Klf2+_cells',
              test.use="MAST")
write.csv2(Klf2_negative_vs_positive_vehicle, file ="Klf2_negative_vs_positive_vehicle.csv")

Klf2_negative_vs_positive_wort <-
  FindMarkers(SO.har, ident.1 = 'treated_Klf2-_cells', ident.2 = 'treated_Klf2+_cells',
              test.use="MAST")
write.csv2(Klf2_negative_vs_positive_wort, file ="Klf2_negative_vs_positive_wort.csv")

# Visualize with three colors
DimPlot(SO.har, group.by = "Klf2_expression",split.by = "Replicates",pt.size = 3) + scale_color_manual(values = c("#1965B0","#7BAFDE","#882E72","#FF7F00")) # for high, low, others

# creating another metadata column in order to compare the individual groups
Idents(SO.har) <- SO.har@meta.data$CD69_expression
SO.har$comparison_CD69 <- paste(SO.har$Condition, sep = "_", Idents(SO.har))
# making this the new Ident
Idents(SO.har) <- "comparison_CD69"

# finding markers between the different groups
CD69_double_positive_vehicle <-
  FindMarkers(SO.har, ident.1 = 'vehicle_CD69+_cells', ident.2 = 'treated_CD69+_cells',
              test.use="MAST")
write.csv2(CD69_double_positive_vehicle, file ="CD69_double_positive_vehicle.csv")

CD69_cite_positive_vehicle <-
  FindMarkers(SO.har, ident.1 = 'vehicle_CD69+_cite_cells', ident.2 = 'treated_CD69+_cite_cells',
              test.use="MAST")
write.csv2(CD69_cite_positive_vehicle, file ="CD69_cite_positive_vehicle.csv")

CD69_negative_vehicle <-
  FindMarkers(SO.har, ident.1 = 'vehicle_CD69-_cells', ident.2 = 'treated_CD69-_cells',
              test.use="wilcox")
write.csv2(CD69_negative_vehicle, file ="CD69_negative_vehicle.csv")



CD69_double_positive_vs_cite_pos_veh <-
  FindMarkers(SO.har, ident.1 = 'vehicle_CD69+_cells', ident.2 = 'vehicle_CD69+_cite_cells',
              test.use="wilcox")
write.csv2(CD69_double_positive_vs_cite_pos_veh, file ="CD69_double_positive_vs_cite_pos_veh.csv")

CD69_double_positive_vs_neg_veh <-
  FindMarkers(SO.har, ident.1 = 'vehicle_CD69+_cells', ident.2 = 'vehicle_CD69-_cells',
              test.use="wilcox")
write.csv2(CD69_double_positive_vs_neg_veh, file ="CD69_double_positive_vs_neg_veh.csv")

CD69_cite_positive_vs_neg_veh <-
  FindMarkers(SO.har, ident.1 = 'vehicle_CD69+_cite_cells', ident.2 = 'vehicle_CD69-_cells',
              test.use="wilcox")
write.csv2(CD69_cite_positive_vs_neg_veh, file ="CD69_cite_positive_vs_neg_veh.csv")



CD69_double_positive_vs_cite_pos_wort <-
  FindMarkers(SO.har, ident.1 = 'treated_CD69+_cells', ident.2 = 'treated_CD69+_cite_cells',
              test.use="wilcox")
write.csv2(CD69_double_positive_vs_cite_pos_wort, file ="CD69_double_positive_vs_cite_pos_wort.csv")

CD69_double_positive_vs_neg_wort <-
  FindMarkers(SO.har, ident.1 = 'treated_CD69+_cells', ident.2 = 'treated_CD69-_cells',
              test.use="wilcox")
write.csv2(CD69_double_positive_vs_neg_wort, file ="CD69_double_positive_vs_neg_wort.csv")

CD69_cite_positive_vs_neg_wort <-
  FindMarkers(SO.har, ident.1 = 'treated_CD69+_cite_cells', ident.2 = 'treated_CD69-_cells',
              test.use="wilcox")
write.csv2(CD69_cite_positive_vs_neg_wort, file ="CD69_cite_positive_vs_neg_wort.csv")





SO.har_subset <- subset(SO.har, cells = CD69_cells)

DimPlot(SO.har_subset,split.by = "Replicates")

# Run Harmony
SO.har <- NormalizeData(SO.har_subset)
SO.har <- FindVariableFeatures(SO.har)
SO.har <- ScaleData(SO.har)
SO.har <- RunPCA(SO.har)
ElbowPlot(SO.har, ndims = 30, reduction = "pca")
SO.har <- RunHarmony(SO.har, group.by.vars = "Replicates", plot_convergence = TRUE)
SO.har <- RunUMAP(SO.har, reduction = "harmony", dims = 1:20)
SO.har <- FindNeighbors(SO.har, reduction = "harmony", dims =1:20)
SO.har <- FindClusters(SO.har, resolution =  c(0.1,0.2,0.3))

Disc_colors = c("#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72","#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",

                "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#E6AB02","#7570B3","#BEAED4","#666666","#999999",

                "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000","#aeae5c","#1e90ff","#00bfff","#56ff0d","#ffff00",



                "#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72","#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",

                "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#E6AB02","#7570B3","#BEAED4","#666666","#999999",

                "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000","#aeae5c","#1e90ff","#00bfff","#56ff0d","#ffff00"

)
use_these_colors = c(Disc_colors[1:length(Idents(SO.har))])

DimPlot(SO.har, group.by = "RNA_snn_res.0.2",
        split.by = "Condition",
        label.size = 10, pt.size = .5,label = T, ncol = 2,
        cols = use_these_colors)

FeaturePlot_scCustom(SO.har, features = "Klf2",
                     label.size = 10, pt.size = .8,label = F
                     , split.by = "Condition",ncol=2
) & theme(legend.position = "right")

SaveSeuratRds(SO.har, file = "./SO.BM.CD4+gating.RDS")
