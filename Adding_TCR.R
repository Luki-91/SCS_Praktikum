library(Seurat)
library(tidyverse)

Disc_colors = c("#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72","#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",

                "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#E6AB02","#7570B3","#BEAED4","#666666","#999999",

                "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000","#aeae5c","#1e90ff","#00bfff","#56ff0d","#ffff00",



                "#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72","#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",

                "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#E6AB02","#7570B3","#BEAED4","#666666","#999999",

                "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000","#aeae5c","#1e90ff","#00bfff","#56ff0d","#ffff00"

)


#####Add TCR info ####
setwd("C:/Users/lukas/Documents/CQ_Praktikum/TCR/outs/vdj_t")

tcr_v <- read.csv("filtered_contig_annotations.csv")
tcr_v <- tcr_v[,c("barcode", "cdr3","cdr3_nt", "chain")]
tcr_v$TRA <- NA
tcr_v$TRB <- NA

# renaming barcodes, since in the SO the first two experiments,
# were the ones, which were removed
tcr_v$barcode <- gsub("-4", "-6", tcr_v$barcode)
tcr_v$barcode <- gsub("-3", "-5", tcr_v$barcode)
tcr_v$barcode <- gsub("-2", "-4", tcr_v$barcode)
tcr_v$barcode <- gsub("-1", "-3", tcr_v$barcode)


tcr_df <- tcr_v %>%
  mutate(TRA = case_when(chain=="TRA" ~ as.character(cdr3))) |>
  mutate(TRB = case_when(chain=="TRB" ~ as.character(cdr3)))

df_TRA <- tcr_df |>  filter(tcr_df$chain == "TRA") |> distinct(barcode, .keep_all = T)
df_TRB <- tcr_df %>% filter(tcr_df$chain == "TRB") |> distinct(barcode, .keep_all = T)

rownames(df_TRA) <- df_TRA[,1]
df_TRA[,1] <- NULL #I repeated until only names and TRA sequence stayed-> 4X
df_TRA[,1] <- NULL
df_TRA[,1] <- NULL
df_TRA[,1] <- NULL
df_TRA[,2] <- NULL

rownames(df_TRB) <- df_TRB[,1]
df_TRB[,1] <- NULL #I repeated until only names and TRB sequence stayed-> 5X
df_TRB[,1] <- NULL
df_TRB[,1] <- NULL
df_TRB[,1] <- NULL
df_TRB[,1] <- NULL

aa_cl <- merge(df_TRA, df_TRB, by = "row.names")
aa_cl$Row.names <- str_replace_all(aa_cl$Row.names, "-\\d","-1")
rownames(aa_cl) <- aa_cl[,1]
aa_cl[,1] <- NULL


aa_cl$TRA_TRB <- paste(aa_cl$TRA, aa_cl$TRB, sep = ":")
aa_cl[,1] <- NULL
aa_cl[,1] <- NULL


SO <- AddMetaData(object=SO,
                  metadata=aa_cl,
                  col.name = "aa_clones")

SO@meta.data$aa_clones = as.factor(SO@meta.data$aa_clones)
levels(SO@meta.data$aa_clones)


####PLOT CELLS WITHOUT TCRs ####

Plot_tmp = as.data.frame(SO@reductions$umap@cell.embeddings)

Plot_tmp$clonotype = SO@meta.data$aa_clones

Plot_tmp$col2 = SO@meta.data$RNA_snn_res.0.6
levels(SO@meta.data$RNA_snn_res.0.6)

Plot_tmp$col = ifelse(is.na(Plot_tmp$clonotype),"No clonotype",as.character(Plot_tmp$col2))
Plot_tmp$col <- factor(Plot_tmp$col, levels = c(levels(SO@meta.data$RNA_snn_res.0.6)))
use_these_colors = c(Disc_colors[1:length(levels(Plot_tmp$col2))],"#000000")
Plot = Plot_tmp[order(Plot_tmp$col),]

p <- ggplot(Plot, aes(umap_1, umap_2, color= col)) +
  geom_point(size=0.3) +
  labs(title= paste("UMAP"," highlighting ","cells without TCR",sep = ""),color="") +
  xlab("Umap 1") +
  ylab("Umap 2") +
  scale_color_manual(values = use_these_colors,na.value = "#000000") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(size=20),
        legend.text=element_text(size=10,face="bold"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)) +
  guides(color = guide_legend(override.aes = list(size = 5)))
p

####Subset no TCRs####

summary(is.na(SO@meta.data$aa_clones))

# Convert aa_clones to character type
SO@meta.data$aa_clones <- as.character(SO@meta.data$aa_clones)

# Assign 'NoTCR' to missing values
SO@meta.data$aa_clones[is.na(SO@meta.data$aa_clones)] <- 'NoTCR'

# Convert aa_clones to character type
SO@meta.data$aa_clones <- as.factor(SO@meta.data$aa_clones)

# Check if it worked
summary(SO@meta.data$aa_clones[SO@meta.data$aa_clones == "NoTCR"])


SO_TCR <- subset(x = SO, subset = (aa_clones == "NoTCR"), invert = T)

summary(SO_TCR@meta.data$aa_clones[SO_TCR@meta.data$aa_clones == "NoTCR"])

summary(is.na(SO_TCR@meta.data$aa_clones))

SaveSeuratRds(SO_TCR, file = "./SO_TCR.RDS")

####Bubble plots of only top genes####

# CD4h_FAM <- FindAllMarkers(SO.harmony, only.pos = TRUE,
#                            min.pct = 0.5,
#                            logfc.threshold = 0.25)
#
# CD4h_FAM_top10 <- CD4h_FAM %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# CD4h_FAM_top10 <- as.factor(unique(CD4h_FAM_top10$gene))