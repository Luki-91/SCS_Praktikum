library(Seurat)
library(tidyverse)
library(scRepertoire)

setwd("C:/Users/lukas/Documents/CQ_Praktikum")
SO_TCR <- readRDS("SO_TCR.RDS")
setwd("C:/Users/lukas/Documents/CQ_Praktikum/TCR/outs/vdj_t")
tcr_v <- read.csv("filtered_contig_annotations.csv")

# renaming barcodes, since in the SO the first two experiments,
# were the ones, which were removed
tcr_v$barcode <- gsub("-4", "-6", tcr_v$barcode)
tcr_v$contig_id <- gsub("-4", "-6", tcr_v$contig_id)
tcr_v$barcode <- gsub("-3", "-5", tcr_v$barcode)
tcr_v$contig_id <- gsub("-3", "-5", tcr_v$contig_id)
tcr_v$barcode <- gsub("-2", "-4", tcr_v$barcode)
tcr_v$contig_id <- gsub("-2", "-4", tcr_v$contig_id)
tcr_v$barcode <- gsub("-1", "-3", tcr_v$barcode)
tcr_v$contig_id <- gsub("-1", "-3", tcr_v$contig_id)

contig.list <- loadContigs(tcr_v, format = "10X")

contig.list <- createHTOContigList(contig = tcr_v,
                                   sc.data = SO_TCR,
                                   group.by = "Replicates")

head(contig_list[[1]])

combined.TCR <- combineTCR(contig.list,
                           samples = c("BM_Veh_1", "BM_Veh_2",
                                       "BM_Wort_1", "BM_Wort_2",
                                       "Spleen_veh_1", "Spleen_veh_2",
                                       "Spleen_Wort_1", "Spleen_Wort_2"),
                           removeNA = T,
                           removeMulti = FALSE,
                           filterMulti = F)

all.TCR <- combineTCR(contig.list,
                           removeNA = T,
                           removeMulti = FALSE,
                           filterMulti = F)

head(combined.TCR[[1]])

clonalQuant(combined.TCR,
            cloneCall="strict",
            chain = "both",
            scale = TRUE)

clonalAbundance(combined.TCR,
                cloneCall = "gene",
                scale = FALSE)

clonalLength(combined.TCR,
             cloneCall="aa",
             chain = "both")

clonalCompare(combined.TCR,
              top.clones = 10,
              samples = c("BM_Veh_1","BM_Veh_2","BM_Wort_1", "BM_Wort_2"),
              cloneCall="gene",
              graph = "alluvial")

clonalScatter(combined.TCR,
              cloneCall ="strict",
              x.axis = "BM_Veh_1",
              y.axis = "Spleen_veh_1",
              dot.size = "total",
              graph = "proportion")

clonalHomeostasis(combined.TCR,
                  cloneCall = "gene")

clonalProportion(combined.TCR,
                 cloneCall = "gene")

percentAA(combined.TCR,
          chain = "TRA",
          aa.length = 20)

positionalEntropy(combined.TCR,
                  chain = "TRA",
                  aa.length = 20)

clonalDiversity(combined.TCR,
                cloneCall = "gene")

clonalRarefaction(combined.TCR,
                  plot.type = 2,
                  hill.numbers = 1,
                  n.boots = 2)

sub_combined <- clonalCluster(combined.TCR[[2]],
                              chain = "TRA",
                              sequence = "aa",
                              threshold = 0.85,
                              group.by = NULL)

head(sub_combined[,c(1,2,13)])

combined.TCR$BM_Veh_1$barcode <-
  str_remove_all(combined.TCR$BM_Veh_1$barcode, "BM_Veh_1_")
combined.TCR$BM_Veh_2$barcode <-
  str_remove_all(combined.TCR$BM_Veh_2$barcode, "BM_Veh_2_")
combined.TCR$BM_Wort_1$barcode <-
  str_remove_all(combined.TCR$BM_Wort_1$barcode, "BM_Wort_1_")
combined.TCR$BM_Wort_2$barcode <-
  str_remove_all(combined.TCR$BM_Wort_2$barcode, "BM_Wort_2_")
combined.TCR$Spleen_veh_1$barcode <-
  str_remove_all(combined.TCR$Spleen_veh_1$barcode, "Spleen_veh_1")
combined.TCR$Spleen_veh_2$barcode <-
  str_remove_all(combined.TCR$Spleen_veh_2$barcode, "Spleen_veh_2")
combined.TCR$Spleen_Wort_1$barcode <-
  str_remove_all(combined.TCR$Spleen_Wort_1$barcode, "Spleen_Wort_1")
combined.TCR$Spleen_Wort_2$barcode <-
  str_remove_all(combined.TCR$Spleen_Wort_2$barcode, "Spleen_Wort_2")


BM_Veh_1 <- subsetClones(combined.TCR, "sample", "BM_Veh_1")
BM_Veh_2 <- subsetClones(combined.TCR, "sample", "BM_Veh_2")
BM_Wort_1 <- subsetClones(combined.TCR, "sample", "BM_Wort_1")
BM_Wort_2 <- subsetClones(combined.TCR, "sample", "BM_Veh_1")
Spleen_veh_1 <- subsetClones(combined.TCR, "sample", "Spleen_veh_1")
Spleen_veh_2 <- subsetClones(combined.TCR, "sample", "Spleen_veh_2")
Spleen_Wort_1 <- subsetClones(combined.TCR, "sample", "Spleen_Wort_1")
Spleen_Wort_2 <- subsetClones(combined.TCR, "sample", "BM_Veh_1")

#subset_list <- c(BM_Veh_1,BM_Veh_2,BM_Wort_1,BM_Wort_2,
#                 Spleen_veh_1,Spleen_veh_2,Spleen_Wort_1,Spleen_Wort_2)

#### Doesn't work in a loop: the first argument has to be replaced everytime
SO_TCR <- combineExpression(BM_Veh_1,
                            SO_TCR,
                            cloneCall="gene",
                            proportion = TRUE)


Disc_colors = c("#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72","#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",

                "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#E6AB02","#7570B3","#BEAED4","#666666","#999999",

                "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000","#aeae5c","#1e90ff","#00bfff","#56ff0d","#ffff00",



                "#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72","#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",

                "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#E6AB02","#7570B3","#BEAED4","#666666","#999999",

                "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000","#aeae5c","#1e90ff","#00bfff","#56ff0d","#ffff00"

)
use_these_colors = c(Disc_colors[1:length(Idents(SO_TCR))])

DimPlot(SO_TCR, group.by = "cloneSize",
        split.by = "Condition",
        label.size = 10, pt.size = .8,label = F, ncol = 2,
        cols = use_these_colors)
