#download the datasets from Hammond et al, 2019 https://github.com/samuel-marsh/Hammond-et-al_2019_Microglia_scRNAseq
wget -O Hammond_et-al-2019_Seurat_Converted_v4.qs https://figshare.com/ndownloader/files/37624052 #use terminal
install.packages('qs')
library(qs)

hammond_all_samples@active.ident <- hammond_all_samples$Age
hammond <- subset(hammond_all_samples, subset = Paper_Cluster == 'Mono/Mac', invert = TRUE)
table(hammond$Paper_Cluster)

hammond.markers <- FindAllMarkers(hammond, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

hammond.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

hammond.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top20
DoHeatmap(hammond, features = top20$gene) + NoLegend()

hammond_normal <- subset(hammond, idents = c('E14','P4/P5','P30','P100','P540'))
#rerun DEG for normal only
