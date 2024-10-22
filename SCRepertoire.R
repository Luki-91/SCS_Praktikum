library(Seurat)
library(tidyverse)
library(scRepertoire)

setwd("C:/Users/lukas/Documents/CQ_Praktikum")
SO.har <- readRDS("SO.har.RDS")
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
                                   sc.data = SO.har,
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
