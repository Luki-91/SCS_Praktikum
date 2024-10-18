library(ggplot2)
library(ggrepel)

setwd("C:/Users/lukas/Documents/CQ_Praktikum/full_DEG+CON")
# List of cluster numbers
clusters <- seq(1,12)  # Add all your cluster numbers here

for (cluster in clusters) {

  # Step 1: Read the CSV for the current cluster
  df <- read.csv2(paste0("cluster", cluster, ".csv"))

  # Step 2: Process the data
  df$diffexpressed <- "NO"
  df$diffexpressed[df$avg_log2FC > 0.5 & df$p_val_adj < 0.05] <- "Vehicle"
  df$diffexpressed[df$avg_log2FC < -0.5 & df$p_val_adj < 0.05] <- "Wortmannin"

  df$delabel <- NA
  df$delabel[df$diffexpressed != "NO"] <- df$X[df$diffexpressed != "NO"]

  # Step 3: Generate the plot
  p <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
    geom_point() +
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=c("black", "blue",  "red")) +
    geom_vline(xintercept=c(-0.5, 0.5), col="red") +
    ggtitle(paste("Cluster", cluster, "Vehicle vs. Wortmannin"))

  # Step 4: Save the plot with specified width of 1000 pixels
  ggsave(filename = paste0("Cluster_", cluster, "_Vehicle_vs_Wortmannin.png"),
         plot = p, width = 3000, height = 1800, units = "px",bg="white")

  # Optional: Print a message to track progress
  print(paste("Plot for Cluster", cluster, "saved."))
}
