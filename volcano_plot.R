library(ggrepel)
setwd("C:/Users/lukas/Documents/CQ_Praktikum/full_DEG+CON")


df<-read.csv2('cluster12.csv')

# ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val_adj))) + geom_point()

# add a column of NAs
df$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$avg_log2FC > 0.5 & df$p_val_adj < 0.05] <- "Vehicle"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$avg_log2FC < -0.5 & df$p_val_adj < 0.05] <- "Wortmannin"

#p <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed)) + geom_point() + theme_minimal()

df$delabel <- NA
df$delabel[df$diffexpressed != "NO"] <- df$X[df$diffexpressed != "NO"]

#ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
#  geom_point() +
#  theme_minimal() +
#  geom_text()

# Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
# plot adding up all layers we have seen so far
ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("black", "blue",  "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  ggtitle('Cluster 12 Vehicel vs. Wortmannin')
