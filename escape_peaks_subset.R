p <- data.frame(rgcs@meta.data$GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE[rgcs$stage == 'Adult'])
pvis <- ggplot(p, aes(x=rgcs.meta.data.GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE.rgcs.stage....)) + geom_density()
ggplotly(pvis)
cellidspos <- names(which(rgcs$stage == 'Adult' & rgcs$GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE >= 0.45))
cellidsneg <- names(which(rgcs$stage == 'Adult' & rgcs$GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE < 0.45))
rgcs <- SetIdent(rgcs, cells = cellidspos, value = 'Maintenance of syn structure +')
rgcs <- SetIdent(rgcs, cells = cellidsneg, value = 'Maintenance of syn structure -')

#Same for FD125,FD82,FD59

rgcs$maint_syn_str <- rgcs@active.ident

cellidspos <- names(which(rgcs$stage == 'Adult' & rgcs$GOBP_MAINTENANCE_OF_PRESYNAPTIC_ACTIVE_ZONE_STRUCTURE >= 0.714))
cellidsneg <- names(which(rgcs$stage == 'Adult' & rgcs$GOBP_MAINTENANCE_OF_PRESYNAPTIC_ACTIVE_ZONE_STRUCTURE < 0.471))
cellidsmid <- names(which(rgcs$stage == 'Adult' & rgcs$GOBP_MAINTENANCE_OF_PRESYNAPTIC_ACTIVE_ZONE_STRUCTURE >= 0.471 & 
rgcs$GOBP_MAINTENANCE_OF_PRESYNAPTIC_ACTIVE_ZONE_STRUCTURE < 0.714))
rgcs <- SetIdent(rgcs, cells = cellidspos, value = 'Maintenance of presyn act zone structure +')
rgcs <- SetIdent(rgcs, cells = cellidsneg, value = 'Maintenance of presyn act zone structure -')
rgcs <- SetIdent(rgcs, cells = cellidsmid, value = 'Maintenance of presyn act zone structure +-')

rgcsfd59 <- subset(rgcs, subset = stage == 'FD59')

ProcessSeu <- function(Seurat){
Seurat <- NormalizeData(Seurat)
Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
Seurat <- ScaleData(Seurat)
Seurat <- RunPCA(Seurat)
Seurat <- FindNeighbors(Seurat, dims = 1:25)
Seurat <- FindClusters(Seurat, resolution = 2)
Seurat <- RunUMAP(Seurat, dims = 1:25)
Seurat <- RunTSNE(Seurat,  dims.use = 1:10 )
DimPlot(object = Seurat, reduction = "umap")
return (Seurat)
}

rgcsfd59 <- ProcessSeu(rgcsfd59)
 
 FeaturePlot(rgcsfd59, reduction = 'umap', features = c('POU4F2','RBPMS'))
 DimPlot(rgcsfd59, reduction = 'umap', label = TRUE, label.box = TRUE)
 DimPlot(rgcsfd59, reduction = 'umap', label = TRUE, label.box = TRUE, split.by = 'presyn_act_zone_org') + ggtitle('FD59')
 
 #maint_syn_strmaint_presyn_act_zone_strmaint_postsyn_spec_strpos_reg_syn_plastpostsyn_act_cytosk_orgpresyn_act_zone_org
 
 rgcsfd59@active.ident <- rgcsfd59$maint_syn_str
 
 markers <- FindAllMarkers(rgcsfd59, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
    
    markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top20
DoHeatmap(rgcsfd59, features = top20$gene) + NoLegend()

#After merging labels:
rgcsfd59 <- RenameIdents(rgcsfd59, 'Maintenance of syn structure - Maintenance of postsyn spec structure - Maintenance of presyn act zone structure -
 Pos reg syn plast +- Postsyn act cytosk org +- Presyn act zone org -' = '1', 'Maintenance of syn structure + Maintenance of postsyn spec structure -
 Maintenance of presyn act zone structure + Pos reg syn plast - Postsyn act cytosk org - Presyn act zone org +' = '3', 'Maintenance of syn structure
 + Maintenance of postsyn spec structure - Maintenance of presyn act zone structure +- Pos reg syn plast + Postsyn act cytosk org +- Presyn act zone org +-' = '3.5')
