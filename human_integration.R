ProcessInt <- function(data.integrated){
+     data.integrated <- ScaleData(data.integrated, verbose = FALSE)
+     data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = FALSE)
+     data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
+     data.integrated <- FindClusters(data.integrated, resolution = 0.5)
+     data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:20)
+     data.integrated <- RunTSNE(data.integrated,  dims.use = 1:10 )
+ }
adult$annotation <- adult@active.ident
FD59$annotation <- FD59@active.ident
FD82$annotation <- FD82@active.ident
FD125$annotation <- FD125@active.ident
integration_list <- list(adult, FD59, FD82, FD125)
features <- SelectIntegrationFeatures(object.list = integration_list)

data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)
data.combined <- IntegrateData(anchorset = data.anchors)
human.combined <- ProcessInt(data.combined)
