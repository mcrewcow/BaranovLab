umap = cbind("Barcode" = rownames(Embeddings(object = FD82, reduction = "umap")), Embeddings(object = FD82, reduction = "umap"))
write.table(umap, file="C://Users/Emil/10X/umap_FD82.csv", sep = ",", quote = F, row.names = F, col.names = T)
m <- read.csv("C://Users/Emil/10X/umap_FD82.csv",header=TRUE) #change the coords with FA
rownames(m) <- m$Barcode
m <- m[,-1]
m <- as.matrix(m)
FD82[['FA']] <- CreateDimReducObject(embeddings = m, key = 'X_draw_graph_fa_', assay = 'RNA')
FD82
DimPlot(FD82, reduction = 'FA', label = TRUE, repel = TRUE, label.box = TRUE)
colnames(x = FD82[["FA"]]@cell.embeddings) <- paste0("X_draw_graph_fa_", 1:2)
FD82@reductions[["FA"]]@global <- TRUE
colnames(x = FD82[["FA"]]@cell.embeddings) <- paste0("Xdrawgraphfa_", 1:2)
