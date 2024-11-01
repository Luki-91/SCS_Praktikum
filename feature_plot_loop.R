library(Seurat)
library(tidyverse)
library(scCustomize)

setwd("C:/Users/lukas/Documents/CQ_Praktikum/BM_reclustered/CD4+/res_0.4")
SO.har <- readRDS("SO.BM.CD4+.RDS")
setwd("C:/Users/lukas/Documents/CQ_Praktikum/plots/CD4+_res_0.4")

features  <-  c("Rnaset2a","Cysltr2","Slamf7","Ddx19a","Isoc2b","Ankrd46","Gm10076","Gm12596","Scamp3")

for (marker in features) {
  p <- FeaturePlot_scCustom(SO.har, features = marker,
                   label.size = 10, pt.size = .8,label = F
                   , split.by = "Condition",ncol=2
  ) & theme(legend.position = "right")
  print(p)
  # Step 4: Save the plot with specified width of 1000 pixels
  ggsave(filename = paste0(marker, "_feature_plot.png"),
         plot = p, width = 4500, height = 1800, units = "px",bg="white")

  # Optional: Print a message to track progress
  print(paste("Plot for feature", marker, "saved."))
}

cluster_1_marker_gene_list <- list(c("S100a4","Tnfrsf4","Ikzf2","Stat1","Itm2c","mt-Atp8","9130401M01Rik","Klrg1","Cd74","Cd81","Cish","Ifi27l2a"))
object <- AddModuleScore(object = SO.har, features = cluster_1_marker_gene_list, name = "cluster_1_marker_gene_list")
FeaturePlot_scCustom(seurat_object = object, features = "cluster_1_marker_gene_list1",
            label.size = 10, pt.size = .8,label = F
            , split.by = "Condition",ncol=2
) & theme(legend.position = "right")

cluster_3_marker_gene_list <- list(c("Eomes","Izumo1r","Isg15","Stat1","Junb","mt-Atp8","Atp5k","Sp110","Parp14","Ddx5","Cd6","Ifi27l2a"))
object <- AddModuleScore(object = SO.har, features = cluster_3_marker_gene_list, name = "cluster_3_marker_gene_list")
FeaturePlot_scCustom(seurat_object = object, features = "cluster_3_marker_gene_list1",
            label.size = 10, pt.size = .8,label = F
            , split.by = "Condition",ncol=2
) & theme(legend.position = "right")

cluster_6_marker_gene_list <- list(c("Cysltr2","Pef1","Asb6","Ifitm10","Dusp2","Klrk1","Atp5k","mt-Atp8","Faap24","Ikzf3","Ctsw","Lemd2"))
object <- AddModuleScore(object = SO.har, features = cluster_6_marker_gene_list, name = "cluster_6_marker_gene_list")
FeaturePlot_scCustom(seurat_object = object, features = "cluster_6_marker_gene_list1",
                     label.size = 10, pt.size = .8,label = F
                     , split.by = "Condition",ncol=2
) & theme(legend.position = "right")

Wort_marker_gene_list <- list(c("Il7r","Rnaseh2a","Gas5","Ifngr1","Bcl2","Resf1","Il2ra","Atp1b3","2410006H16Rik","Npm1","Ccl5","Rps27rt"))
object <- AddModuleScore(object = SO.har, features = Wort_marker_gene_list, name = "Wort_marker_gene_list")
FeaturePlot_scCustom(seurat_object = object, features = "Wort_marker_gene_list1",
                     label.size = 10, pt.size = .8,label = F
                     , split.by = "Condition",ncol=2
) & theme(legend.position = "right")
