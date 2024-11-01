library(ggrepel)
setwd("C:/Users/lukas/Documents/CQ_Praktikum/full_DEG+CON")


df<-read.csv2('cluster_cluster_2357_ml_signatures.csv')

# ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val_adj))) + geom_point()

# add a column of NAs
df$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$log2FC > 0.2] <- "Vehicle"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$log2FC < -0.2] <- "Wortmannin"

#p <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed)) + geom_point() + theme_minimal()

df$delabel <- NA
df$delabel[df$diffexpressed != "NO"] <- df$gene[df$diffexpressed != "NO"]

#ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
#  geom_point() +
#  theme_minimal() +
#  geom_text()

# Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
# plot adding up all layers we have seen so far
ggplot(data=df, aes(x=log2FC, y=importance_score, col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("black", "blue",  "red"))+
  ggtitle('Cluster 2357 Vehicle vs. Wortmannin')
