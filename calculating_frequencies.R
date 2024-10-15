# load libraries
library(Seurat)
library(tidyverse)

setwd("C:/Users/lukas/Documents/CQ_Praktikum")

# Read RDS created in Seurat Video Tutorials--video 7
subsetted <- readRDS("subsetted.RDS")
DimPlot(subsetted, reduction = "umap", label = TRUE)

view(subsetted@meta.data)


# cluster_frequency
# Step 1: Calculate total cell count per replicate
cluster_total <- subsetted@meta.data %>%
  group_by(Replicates) %>%
  summarise(total_count = n())

# Step 2: Calculate the count per cluster and replicate
cluster_frequency <- subsetted@meta.data %>%
  group_by(Replicates, seurat_clusters) %>%
  summarise(count = n())

# Step 3: Join the total cell count to the cluster_frequency table
cluster_frequency <- cluster_frequency %>%
  left_join(cluster_total, by = "Replicates")

# Step 4: Calculate the relative frequency
cluster_frequency <- cluster_frequency %>%
  mutate(relative_frequency = count / total_count) %>%
  mutate(data_set = "subsetted")

ggplot(cluster_frequency, aes(x = as.factor(seurat_clusters), y = relative_frequency, fill = Replicates)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster Number", y = "Relative Frequency", title = "Relative Frequency of Clusters by Replicate") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")

cluster_frequency <- cluster_frequency %>%
  mutate(Organ = case_when(
    grepl("BM", Replicates, ignore.case = TRUE) ~ "BM",
    grepl("Spleen", Replicates, ignore.case = TRUE) ~ "Spleen"))

cluster_frequency <- na.omit(cluster_frequency)

ggplot(cluster_frequency, aes(x = as.factor(seurat_clusters), y = relative_frequency, fill = Replicates)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Organ) +  # Split the plot by the new 'Group' column
  labs(x = "Cluster Number", y = "Relative Frequency", title = "Relative Frequency of Clusters by Organ") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")

write.csv(cluster_frequency, file="../Seurat object/cluster_frequency.CSV")

ggplot(cluster_frequency, aes(x = data_set, y = count,
                              fill=seurat_clusters)) + geom_col()

