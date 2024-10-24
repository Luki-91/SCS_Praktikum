library(Seurat)
library(tidyverse)
library(harmony)

setwd("C:/Users/lukas/Documents/CQ_Praktikum/BM_reclustered/CD8+")
SO.har <- readRDS("SO.BM.CD8+.RDS")

setwd("C:/Users/lukas/Documents/CQ_Praktikum/BM_reclustered/CD8+/mast")
Idents(SO.har) <- "comparison_CD8"
DefaultAssay(SO.har)

#### Finding markers using MAST ####
cluster0 <- FindMarkers(SO.har, ident.1 = 'vehicle_0', ident.2 = 'treated_0',
                        test.use="MAST")
write.csv2(cluster0, file ="cluster0.csv")
cluster1 <- FindMarkers(SO.har, ident.1 = 'vehicle_1', ident.2 = 'treated_1',
                        test.use="MAST")
write.csv2(cluster1, file ="cluster1.csv")
cluster2 <- FindMarkers(SO.har, ident.1 = 'vehicle_2', ident.2 = 'treated_2',
                        test.use="MAST")
write.csv2(cluster2, file ="cluster2.csv")
cluster3 <- FindMarkers(SO.har, ident.1 = 'vehicle_3', ident.2 = 'treated_3',
                        test.use="MAST")
write.csv2(cluster3, file ="cluster3.csv")
cluster4 <- FindMarkers(SO.har, ident.1 = 'vehicle_4', ident.2 = 'treated_4',
                        test.use="MAST")
write.csv2(cluster4, file ="cluster4.csv")
cluster5 <- FindMarkers(SO.har, ident.1 = 'vehicle_5', ident.2 = 'treated_5',
                        test.use="MAST")
write.csv2(cluster5, file ="cluster5.csv")
cluster6 <- FindMarkers(SO.har, ident.1 = 'vehicle_6', ident.2 = 'treated_6',
                        test.use="MAST")
write.csv2(cluster6, file ="cluster6.csv")
cluster7 <- FindMarkers(SO.har, ident.1 = 'vehicle_7', ident.2 = 'treated_7',
                        test.use="MAST")
write.csv2(cluster7, file ="cluster7.csv")
cluster8 <- FindMarkers(SO.har, ident.1 = 'vehicle_8', ident.2 = 'treated_8',
                        test.use="MAST")
write.csv2(cluster8, file ="cluster8.csv")


#### Finding markers using roc ####
setwd("C:/Users/lukas/Documents/CQ_Praktikum/BM_reclustered/CD8+/roc")
cluster0 <- FindMarkers(SO.har, ident.1 = 'vehicle_0', ident.2 = 'treated_0',
                        test.use="roc")
write.csv2(cluster0, file ="cluster0.csv")
cluster1 <- FindMarkers(SO.har, ident.1 = 'vehicle_1', ident.2 = 'treated_1',
                        test.use="roc")
write.csv2(cluster1, file ="cluster1.csv")
cluster2 <- FindMarkers(SO.har, ident.1 = 'vehicle_2', ident.2 = 'treated_2',
                        test.use="roc")
write.csv2(cluster2, file ="cluster2.csv")
cluster3 <- FindMarkers(SO.har, ident.1 = 'vehicle_3', ident.2 = 'treated_3',
                        test.use="roc")
write.csv2(cluster3, file ="cluster3.csv")
cluster4 <- FindMarkers(SO.har, ident.1 = 'vehicle_4', ident.2 = 'treated_4',
                        test.use="roc")
write.csv2(cluster4, file ="cluster4.csv")
cluster5 <- FindMarkers(SO.har, ident.1 = 'vehicle_5', ident.2 = 'treated_5',
                        test.use="roc")
write.csv2(cluster5, file ="cluster5.csv")
cluster6 <- FindMarkers(SO.har, ident.1 = 'vehicle_6', ident.2 = 'treated_6',
                        test.use="roc")
write.csv2(cluster6, file ="cluster6.csv")
cluster7 <- FindMarkers(SO.har, ident.1 = 'vehicle_7', ident.2 = 'treated_7',
                        test.use="roc")
write.csv2(cluster7, file ="cluster7.csv")
cluster8 <- FindMarkers(SO.har, ident.1 = 'vehicle_8', ident.2 = 'treated_8',
                        test.use="roc")
write.csv2(cluster8, file ="cluster8.csv")


#### Finding markers using bimod ####
setwd("C:/Users/lukas/Documents/CQ_Praktikum/BM_reclustered/CD8+/bimod")
cluster0 <- FindMarkers(SO.har, ident.1 = 'vehicle_0', ident.2 = 'treated_0',
                        test.use="bimod")
write.csv2(cluster0, file ="cluster0.csv")
cluster1 <- FindMarkers(SO.har, ident.1 = 'vehicle_1', ident.2 = 'treated_1',
                        test.use="bimod")
write.csv2(cluster1, file ="cluster1.csv")
cluster2 <- FindMarkers(SO.har, ident.1 = 'vehicle_2', ident.2 = 'treated_2',
                        test.use="bimod")
write.csv2(cluster2, file ="cluster2.csv")
cluster3 <- FindMarkers(SO.har, ident.1 = 'vehicle_3', ident.2 = 'treated_3',
                        test.use="bimod")
write.csv2(cluster3, file ="cluster3.csv")
cluster4 <- FindMarkers(SO.har, ident.1 = 'vehicle_4', ident.2 = 'treated_4',
                        test.use="bimod")
write.csv2(cluster4, file ="cluster4.csv")
cluster5 <- FindMarkers(SO.har, ident.1 = 'vehicle_5', ident.2 = 'treated_5',
                        test.use="bimod")
write.csv2(cluster5, file ="cluster5.csv")
cluster6 <- FindMarkers(SO.har, ident.1 = 'vehicle_6', ident.2 = 'treated_6',
                        test.use="bimod")
write.csv2(cluster6, file ="cluster6.csv")
cluster7 <- FindMarkers(SO.har, ident.1 = 'vehicle_7', ident.2 = 'treated_7',
                        test.use="bimod")
write.csv2(cluster7, file ="cluster7.csv")
cluster8 <- FindMarkers(SO.har, ident.1 = 'vehicle_8', ident.2 = 'treated_8',
                        test.use="bimod")
write.csv2(cluster8, file ="cluster8.csv")