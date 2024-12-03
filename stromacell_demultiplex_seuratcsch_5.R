#source("Z:/TempIzar/Emila/seurat_scripts_recent.3.R")
library("Seurat")

setwd("C:/Users/lukas/Documents/CQ_Praktikum/stroma_cells")

loc_name = "Lucas_stroma"


SO = Read10X_h5("LUCAS_StromaAggr.h5")
SO = CreateSeuratObject(counts = SO[["Antibody Capture"]])
library_id <- read.csv("aggregation.csv")

print(unique(experiment))
print(rownames(SO@assays$RNA@counts))
msHashtags <- grep("^mmHashtag", row.names(SO@assays$RNA@features), value = T)
if(length(msHashtags) == 0){
  msHashtags <- grep("^hsHashtag", row.names(SO@assays$RNA@features), value = T)
}
if(length(msHashtags) == 0){
  stop("no hashtags found")
}
SO_hs <- subset(SO, features = msHashtags)
SO_hs = NormalizeData(SO_hs)

experiment <- colnames(SO_hs@assays$RNA)
for (suffix_pos in c(length(library_id[,1]):1)) {
  suffix = paste0("-",suffix_pos, "$")
  h1 <- grep(suffix, colnames(SO_hs@assays$RNA), value = T)
  h1.ID <- match(h1,colnames(SO_hs@assays$RNA))
  experiment[h1.ID] <- library_id[suffix_pos,1]
}

SO_hs@meta.data$experiment <- experiment
#SO_hs = NormalizeData(SO_hs)
#SO_hs = FindVariableFeatures(SO_hs, do.plot = FALSE)
#SO_hs <- ScaleData(SO_hs, verbose = FALSE)
#loc_name = paste0(loc_name, ".hs")
#SO_hs <- RunPCA(SO_hs, npcs = 50, verbose = FALSE, features = msHashtags)
#SO_hs <- RunTSNE(SO_hs, npcs = 50, verbose = FALSE, features = hsHashtags, check_duplicates = FALSE)
#SO_hs <- RunUMAP(SO_hs, npcs = 50, verbose = FALSE, features = hsHashtags, check_duplicates = FALSE)


#SO_hs <- FindNeighbors(SO_hs, dims = 1:dim(SO_hs@reductions$pca)[2], reduction = "pca")
#SO_hs <- FindClusters(SO_hs, resolution = 1, reduction = "pca")
#SO_hs <- FindClusters(SO_hs, resolution = 0.9, reduction = "pca")
#SO_hs <- FindClusters(SO_hs, resolution = 0.8, reduction = "pca")
# SO_hs <- FindClusters(SO_hs, resolution = 0.7, reduction = "pca")
# SO_hs <- FindClusters(SO_hs, resolution = 0.6, reduction = "pca")
# SO_hs <- FindClusters(SO_hs, resolution = 0.5, reduction = "pca")
# SO_hs <- FindClusters(SO_hs, resolution = 0.4, reduction = "pca")
# SO_hs <- FindClusters(SO_hs, resolution = 0.3, reduction = "pca")
# SO_hs <- FindClusters(SO_hs, resolution = 0.2, reduction = "pca")
# SO_hs <- FindClusters(SO_hs, resolution = 0.1, reduction = "pca")
#
# writeJustPCA(SO=SO_hs, raw_name = paste0(loc_name,".SNN.rpca"), no_pca = 50, do_pca = TRUE, do_tsne = TRUE, do_umap = TRUE)
# write_mwtadata(SO_hs, paste0(loc_name,".SNN.rpca", "int_meta"))
#
# write_dgCMatrix_tab(mat=SO_hs@assays$RNA@data, filename = paste0(loc_name, "_expr_all.txt"))
# write_dgCMatrix_tab(mat=SO_hs@assays$RNA@counts, filename = paste0(loc_name, "_counts_all.txt"))
# write_dgCMatrix_tab(mat=SO_hs@assays$RNA@scale.data, filename = paste0(loc_name, "_scale_all.txt"))



check = FALSE
if(check){ #individual lucas
  unq_exp <- unique(experiment)
  exp_list = list()
  for (i in 1:length(unq_exp)) {
    #i <- 0
    exp_list[[i]] = which(SO_hs@meta.data$experiment == unq_exp[i])
  }

  argh <- function(x, scalearg = 1){
    x = x/scalearg
    return(log(x + sqrt(x * x + 1)))
  }


  val_x = SO_hs@assays$RNA@layers$counts
  dim(val_x)
  rnames <- rownames(SO_hs@assays$RNA@features)
  cnames <- rownames(SO_hs@assays$RNA@cells)
  colnames(val_x)
  anno <- rep("unknown", length(cnames))

  #length(anno)

  table(anno)

  ################
  # Experiment 1
  ################

  cur_exp = 1
  idx_list = list()
  unique(SO_hs@meta.data$experiment[exp_list[[cur_exp]]])
  cur_idx = 1
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 200)
  res <- hist(cur_val, breaks = 120, xlim = c(5,10))
  axis(side=1,at=seq(5,10,0.1), labels=seq(5,10,0.1))
  cbind(res$breaks[c(-1)], res$counts) #4.5
  idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 7.3)]]

  cur_idx = 2
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 100)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.5
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.0)]]

  cur_idx = 3
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  cbind(res$breaks[c(-1)], res$counts) #5.0
  idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 6.3)]]

  cur_idx = 4
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  cbind(res$breaks[c(-1)], res$counts) #5.0
  idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 6.5)]]

  cur_idx = 5
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  cur_idx = 6
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  cbind(res$breaks[c(-1)], res$counts) #5.0
  idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.5)]]

  cur_idx = 7
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.5)]]

  cur_idx = 8
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.9)]]

  cur_idx = 9
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  cur_idx = 10
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  table(anno)
  overlaps <- c()
  for(idx_1 in 1:(length(idx_list)-1)){
    for(idx_2 in (idx_1 + 1):(length(idx_list))){
      overlaps <- c(overlaps,intersect(idx_list[[idx_1]],idx_list[[idx_2]]))
      #print(length(overlaps))
    }
  }
  print(length(overlaps))
  overlaps <- unique(overlaps)
  tail(overlaps)
  for(idx_1 in 1:(length(idx_list))){
    print(paste0(rnames[idx_1], ":", length(idx_list[[idx_1]])))
  }

  table(anno)
  for(idx_1 in 1:(length(idx_list))){
    anno[which(cnames %in% idx_list[[idx_1]])] <- rnames[idx_1]
  }
  anno[which(cnames %in% overlaps)] <- "ambigous"
  table(anno)

  ################
  # Experiment 2
  ################

  cur_exp = 2
  idx_list = list()
  unique(SO_hs@meta.data$experiment[exp_list[[cur_exp]]])
  cur_idx = 1
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 200)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #4.5
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.6)]]

  cur_idx = 2
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 100)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.5
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.0)]]

  cur_idx = 3
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  cur_idx = 4
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  cbind(res$breaks[c(-1)], res$counts) #5.0
  idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.6)]]

  cur_idx = 5
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  cur_idx = 6
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  cbind(res$breaks[c(-1)], res$counts) #5.0
  idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.4)]]

  cur_idx = 7
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.5)]]

  cur_idx = 8
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.9)]]

  cur_idx = 9
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  cur_idx = 10
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  table(anno)
  overlaps <- c()
  for(idx_1 in 1:(length(idx_list)-1)){
    for(idx_2 in (idx_1 + 1):(length(idx_list))){
      overlaps <- c(overlaps,intersect(idx_list[[idx_1]],idx_list[[idx_2]]))
      #print(length(overlaps))
    }
  }
  #print(length(overlaps))
  overlaps <- unique(overlaps)
  #tail(overlaps)
  for(idx_1 in 1:(length(idx_list))){
    print(paste0(rnames[idx_1], ":", length(idx_list[[idx_1]])))
  }

  table(anno)
  for(idx_1 in 1:(length(idx_list))){
    anno[which(cnames %in% idx_list[[idx_1]])] <- rnames[idx_1]
  }
  anno[which(cnames %in% overlaps)] <- "ambigous"
  table(anno)

  ################
  # Experiment 3
  ################

  cur_exp = 3
  idx_list = list()
  #unique(SO_hs@meta.data$experiment[exp_list[[cur_exp]]])
  cur_idx = 1
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 200)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #4.5
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.6)]]

  cur_idx = 2
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 100)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.5
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.0)]]

  cur_idx = 3
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  cur_idx = 4
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.8)]]

  cur_idx = 5
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  cbind(res$breaks[c(-1)], res$counts) #5.0
  idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 3.8)]]

  cur_idx = 6
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  cbind(res$breaks[c(-1)], res$counts) #5.0
  idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.2)]]

  cur_idx = 7
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.5)]]

  cur_idx = 8
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.9)]]

  cur_idx = 9
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  cur_idx = 10
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  table(anno)
  overlaps <- c()
  for(idx_1 in 1:(length(idx_list)-1)){
    for(idx_2 in (idx_1 + 1):(length(idx_list))){
      overlaps <- c(overlaps,intersect(idx_list[[idx_1]],idx_list[[idx_2]]))
      #print(length(overlaps))
    }
  }
  #print(length(overlaps))
  overlaps <- unique(overlaps)
  #tail(overlaps)
  for(idx_1 in 1:(length(idx_list))){
    print(paste0(rnames[idx_1], ":", length(idx_list[[idx_1]])))
  }

  table(anno)
  for(idx_1 in 1:(length(idx_list))){
    anno[which(cnames %in% idx_list[[idx_1]])] <- rnames[idx_1]
  }
  anno[which(cnames %in% overlaps)] <- "ambigous"
  table(anno)

  ################
  # Experiment 4
  ################

  cur_exp = 4
  idx_list = list()
  #unique(SO_hs@meta.data$experiment[exp_list[[cur_exp]]])
  cur_idx = 1
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 200)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #4.5
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.6)]]

  cur_idx = 2
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 100)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.5
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.0)]]

  cur_idx = 3
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  cur_idx = 4
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.8)]]

  cur_idx = 5
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 3.8)]]

  cur_idx = 6
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.2)]]

  cur_idx = 7
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  cbind(res$breaks[c(-1)], res$counts) #5.0
  idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.6)]]

  cur_idx = 8
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  cbind(res$breaks[c(-1)], res$counts) #5.0
  idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.3)]]

  cur_idx = 9
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  cur_idx = 10
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  table(anno)
  overlaps <- c()
  for(idx_1 in 1:(length(idx_list)-1)){
    for(idx_2 in (idx_1 + 1):(length(idx_list))){
      overlaps <- c(overlaps,intersect(idx_list[[idx_1]],idx_list[[idx_2]]))
      #print(length(overlaps))
    }
  }
  print(length(overlaps))
  overlaps <- unique(overlaps)
  #tail(overlaps)
  for(idx_1 in 1:(length(idx_list))){
    print(paste0(rnames[idx_1], ":", length(idx_list[[idx_1]])))
  }

  table(anno)
  for(idx_1 in 1:(length(idx_list))){
    anno[which(cnames %in% idx_list[[idx_1]])] <- rnames[idx_1]
  }
  anno[which(cnames %in% overlaps)] <- "ambigous"
  table(anno)

  ################
  # Experiment 5
  ################

  cur_exp = 5
  idx_list = list()
  #unique(SO_hs@meta.data$experiment[exp_list[[cur_exp]]])
  cur_idx = 1
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 200)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #4.5
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.6)]]

  cur_idx = 2
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 100)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.5
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.0)]]

  cur_idx = 3
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  cur_idx = 4
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.8)]]

  cur_idx = 5
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  cbind(res$breaks[c(-1)], res$counts) #5.0
  idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.7)]]

  cur_idx = 6
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  cbind(res$breaks[c(-1)], res$counts) #5.0
  idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.8)]]

  cur_idx = 7
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.6)]]

  cur_idx = 8
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.3)]]

  cur_idx = 9
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  cur_idx = 10
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  table(anno)
  overlaps <- c()
  for(idx_1 in 1:(length(idx_list)-1)){
    for(idx_2 in (idx_1 + 1):(length(idx_list))){
      overlaps <- c(overlaps,intersect(idx_list[[idx_1]],idx_list[[idx_2]]))
      #print(length(overlaps))
    }
  }
  print(length(overlaps))
  overlaps <- unique(overlaps)
  #tail(overlaps)
  for(idx_1 in 1:(length(idx_list))){
    print(paste0(rnames[idx_1], ":", length(idx_list[[idx_1]])))
  }

  table(anno)
  for(idx_1 in 1:(length(idx_list))){
    anno[which(cnames %in% idx_list[[idx_1]])] <- rnames[idx_1]
  }
  anno[which(cnames %in% overlaps)] <- "ambigous"
  table(anno)

  ################
  # Experiment 6
  ################

  cur_exp = 6
  idx_list = list()
  #unique(SO_hs@meta.data$experiment[exp_list[[cur_exp]]])
  cur_idx = 1
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 200)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #4.5
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.6)]]

  cur_idx = 2
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 100)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.5
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.0)]]

  cur_idx = 3
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  cur_idx = 4
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.8)]]

  cur_idx = 5
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.7)]]

  cur_idx = 6
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.8)]]

  cur_idx = 7
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  cbind(res$breaks[c(-1)], res$counts) #5.0
  idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 6.1)]]

  cur_idx = 8
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  cbind(res$breaks[c(-1)], res$counts) #5.0
  idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 6.1)]]

  cur_idx = 9
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  cur_idx = 10
  cur_val = argh(val_x[cur_idx,exp_list[[cur_exp]]])
  res <- hist(cur_val, breaks = 120)
  #res <- hist(cur_val, breaks = 120, xlim = c(3,8))
  #axis(side=1,at=seq(3,8,0.1), labels=seq(3,8,0.1))
  #cbind(res$breaks[c(-1)], res$counts) #5.0
  #idx_list[[cur_idx]] <- cnames[exp_list[[cur_exp]][which(cur_val > 5.65)]]

  table(anno)
  overlaps <- c()
  for(idx_1 in 1:(length(idx_list)-1)){
    for(idx_2 in (idx_1 + 1):(length(idx_list))){
      overlaps <- c(overlaps,intersect(idx_list[[idx_1]],idx_list[[idx_2]]))
      #print(length(overlaps))
    }
  }
  print(length(overlaps))
  overlaps <- unique(overlaps)
  #tail(overlaps)
  for(idx_1 in 1:(length(idx_list))){
    print(paste0(rnames[idx_1], ":", length(idx_list[[idx_1]])))
  }

  table(anno)
  for(idx_1 in 1:(length(idx_list))){
    anno[which(cnames %in% idx_list[[idx_1]])] <- rnames[idx_1]
  }
  anno[which(cnames %in% overlaps)] <- "ambigous"
  table(anno)

  getwd()
  write.table(file = paste0("Lukasanalyzed",loc_name, ".hash.", "csv"), cbind(barcodes = cnames, hashtags = anno), sep = ",", quote = FALSE, row.names = F)
  table(anno) / sum(table(anno))
  SO_hs@meta.data$hashtags <- anno
  colnames(SO_hs@meta.data)
  write_metadata(SO_hs, paste0("Lukasanalyzed",loc_name,".SNN.rpca", "int_meta"))
  writeCounts(SO_hs,ident.2 = "hashtags", sort.2 = FALSE)

}

writeCounts <- function(SO, ident.1 = "experiment", ident.2, sort.2 = TRUE) {
  uq.1 <- unique(SO@meta.data[[ident.1]])
  uq.2 <- unique(SO@meta.data[[ident.2]])
  overlap <- matrix(nrow = length(uq.1),ncol = length(uq.2))
  colnames(overlap) = uq.2
  rownames(overlap) = uq.1
  for(i in 1:length(uq.1)){
    #i = 1
    idx.1 <- row.names(SO@meta.data)[which(SO@meta.data[[ident.1]] == uq.1[i])]
    for(ii in 1:length(uq.2)){
      #   ii = 1
      idx.2 <- row.names(SO@meta.data)[which(SO@meta.data[[ident.2]] == uq.2[ii])]
      overlap[i,ii] <- length(intersect(idx.1, idx.2))
    }
  }
  if(sort.2){
    overlap <- overlap[,sort(as.numeric(colnames(overlap)), index.return = TRUE)$ix]
  } else {
    overlap <- overlap[,sort(colnames(overlap), index.return = TRUE)$ix]
  }
  #overlap <- overlap[sort(as.numeric(rownames(overlap)), index.return = TRUE)$ix,]
  write.table(cbind("Name" = row.names(overlap),overlap), file = paste0(ident.2,".counts.tsv"), sep = "\t", quote = F, row.names = F)
}
