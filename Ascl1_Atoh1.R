#dogs Levi late dataset
levi.data <- Read10X_h5("G://Levi_Jon/late/Ascl1_Atoh1_late_timepoint/outs/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
levi <- CreateSeuratObject(counts = levi.data, project = "late", min.cells = 3, min.features = 200)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
  Seurat <- ScaleData(Seurat)
  Seurat <- RunPCA(Seurat)
  Seurat <- FindNeighbors(Seurat, dims = 1:20)
  Seurat <- FindClusters(Seurat, resolution = 0.5)
  Seurat <- RunUMAP(Seurat, dims = 1:20)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:20 )
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}

levi [["percent.rb"]] <- PercentageFeatureSet(levi , pattern = "^Rps|^Rpl|^Mrps|^Mrpl", assay = 'RNA')

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "www", host = "dec2021.archive.ensembl.org")
  mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "www", host = "dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

levi <- CellCycleScoring(levi, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = FALSE)
levi[["percent.mt"]] <- PercentageFeatureSet(levi, pattern = "^mt-")

levi <- subset(levi, subset = nCount_RNA > 300 & nCount_RNA < 10000 & nFeature_RNA > 600 & nFeature_RNA < 4000 & percent.mt < 15 & percent.rb < 15)

levi <- ScaleData(levi, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))

levi <- ProcessSeu(levi)

RDoublet <- function(tmp){
  sweep.res.list <- paramSweep_v3(tmp, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pKopt <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  pKopt <- pKopt[order(pKopt, decreasing = TRUE) ]
  pKopt <- pKopt[1]
  homotypic.prop <- modelHomotypic(tmp$seurat_clusters) 
  nExp_poi <- round(0.1*length(colnames(tmp)))  ## Assuming 10% doublet formation rate 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi, reuse.pANN = FALSE)
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",pKopt,nExp_poi, sep="_"))
  return (tmp) 
}

remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

levi <- RDoublet(levi)
levi <- subset(levi, subset = DF.classifications_0.25_0.04_241 == 'Singlet') #remove the doublets
levi <- subset(levi, subset = DF.classifications_0.25_0.04_191 == 'Singlet') #remove the doublets

levi <- ProcessSeu(levi)

levi_late <- levi

#repeat for early timepoint

levi_early <- levi

ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = FALSE)
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = FALSE)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
  data.integrated <- FindClusters(data.integrated, resolution = 0.5)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:20)
  data.integrated <- RunTSNE(data.integrated,  dims.use = 1:20 )
}

#Choose the objects for integration

integration_list <- list(levi_late, levi_early)

features <- SelectIntegrationFeatures(object.list = integration_list)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)

data.combined <- IntegrateData(anchorset = data.anchors)

levi.combined <- ProcessInt(data.combined)

levi.combined.markers <- FindAllMarkers(levi.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

levi.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

levi.combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(levi.combined, features = top10$gene) + NoLegend()

SaveH5Seurat(levi.combined, 'G://Levi_Jon/combined.h5Seurat')
Convert("G://Levi_Jon/combined.h5Seurat", dest = "h5ad") #for further in-python processing

library(escape)
gene.sets1 <- getGeneSets(library = "C5", gene.sets = c("GOBP_GLIAL_CELL_MIGRATION",'GOBP_NEURON_MIGRATION'
                                                        ,'GOBP_NEURON_MATURATION','GOBP_GLIAL_CELL_DEVELOPMENT',
                                                        'GOBP_GLIAL_CELL_PROLIFERATION'),species = 'Mus musculus')
levi.combined <- SetIdent(levi.combined, value = 'orig.ident')

ES <- enrichIt(obj = levi.combined,
gene.sets = gene.sets1,
groups = 1000, cores = 8)
levi.combined <- AddMetaData(levi.combined, ES)
ES2 <- data.frame(levi.combined[[]], Idents(levi.combined))
colnames(ES2)[ncol(ES2)] <- "cluster"
ridgeEnrichment(ES2, gene.set = "GOBP_GLIAL_CELL_MIGRATION", group = 'orig.ident', add.rug = TRUE) 

levi.combined <- SetIdent(levi.combined, value = 'seurat_clusters')
levi.combined <- RenameIdents(levi.combined, '0'='Transitory',
                              '1'='Transitory',
                              '2'='Transitory',
                              '3'='Transitory',
                              '4'='Partially muller glia',
                              '5'='Microglia',
                              '6'='RGC',
                              '7'='Cones',
                              '8'='Mature muller glia',
                              '9'='AC',
                              '10'='Early neurons')

levi.combined$EK_anno <- levi.combined@active.ident


levi.subset.glia <- subset(levi.combined, idents = c('Partially muller glia','Mature muller glia'))

levi.subset.glia$EK_anno <- levi.subset.glia@active.ident

ES <- enrichIt(obj = levi.subset.glia,
               gene.sets = gene.sets1,
               groups = 1000, cores = 8)
levi.subset.glia <- AddMetaData(levi.subset.glia, ES)
ES2 <- data.frame(levi.subset.glia[[]], Idents(levi.subset.glia))
colnames(ES2)[ncol(ES2)] <- "cluster"
ridgeEnrichment(ES2, gene.set = "GOBP_GLIAL_CELL_MIGRATION", group = 'orig.ident', add.rug = TRUE) + facet_wrap(~EK_anno) 

levi.subset.RGC <- subset(levi.combined, idents = c('RGC'))

levi.subset.RGC$EK_anno <- levi.subset.RGC@active.ident

ES <- enrichIt(obj = levi.subset.RGC,
               gene.sets = gene.sets1,
               groups = 1000, cores = 8)
levi.subset.RGC <- AddMetaData(levi.subset.RGC, ES)
ES2 <- data.frame(levi.subset.RGC[[]], Idents(levi.subset.RGC))
colnames(ES2)[ncol(ES2)] <- "cluster"
ridgeEnrichment(ES2, gene.set = "GOBP_GLIAL_CELL_MIGRATION", group = 'orig.ident', add.rug = TRUE)

DefaultAssay(levi.subset.RGC) <- 'RNA'
plot <- DotPlot(levi.subset.RGC, features = c('Adgrg1','Adgrl3','Ccr4','Celsr1','Celsr2','Celsr3',
                                              'Cxcr4','Dcc','Drd1','Drd2','Erbb4','Esr2','Fzd3','Gfra3','Gpr173','Il1r1','Nr2f1','Nr2f2',
                                              'Nr4a2','Nsmf','Ntrk2','Ptprz1','Robo1','Robo2','Robo3','Unc5c','Unc5d'),  cols = 'RdBu', dot.scale = 10, split.by = 'orig.ident')
plot + theme(axis.text.x = element_text(angle = 315, family = 'Arial'), axis.text.y = element_text(family = 'Arial'))



library(CellChat)
levi.combined <- RenameIdents(levi.combined, 'AC'='Transitory')
levi.combined <- subset(levi.combined, idents = c('Microglia'), invert = TRUE)
table(levi.combined@active.ident)

levi.combined$EK_anno2 <- levi.combined@active.ident
levi.combined.late <- subset(levi.combined, subset = orig.ident == 'late')
levi.combined.early <- subset(levi.combined, subset = orig.ident == 'early')

DefaultAssay(levi.combined.early) <- 'RNA'
cellchat <- createCellChat(object = levi.combined.early , group.by = "EK_anno2")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat_levi.early <- cellchat

DefaultAssay(levi.combined.late) <- 'RNA'
cellchat <- createCellChat(object = levi.combined.late , group.by = "EK_anno2")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentraliy(cellchat, slot.name = "netP")
cellchat_levi.late <- cellchat

netAnalysis_contribution(cellchat_levi.early, signaling = 'NT')
pairLR.NT <- extractEnrichedLR(cellchat_levi.early, signaling = 'NT', geneLR.return = FALSE)
LR.show <- pairLR.NT[1:3,]
netVisual_individual(cellchat_levi.early, signaling = 'NT', pairLR.use = LR.show, layout = "chord")

mp_genes <- read.table(file = "clipboard", sep = "\t", header=TRUE)
mp_genes <- mp_genes$MP
st_genes <- read.table(file = "clipboard", sep = "\t", header=TRUE)
st_genes <- st_genes$ST
DefaultAssay(levi.subset.RGC) <- 'RNA'
gene.sets <- list(MP = mp_genes, ST = st_genes)

ESr <- enrichIt(obj = levi.subset.RGC,
               gene.sets = gene.sets,
               groups = 1000, cores = 8)
levi.subset.RGC <- AddMetaData(levi.subset.RGC, ESr)
ES2r <- data.frame(levi.subset.RGC[[]], Idents(levi.subset.RGC))
colnames(ES2r)[ncol(ES2r)] <- "cluster"
ES2r$orig.ident <- factor(ES2r$orig.ident , levels = c('late','early'))
ridgeEnrichment(ES2r, gene.set = "MP", group = 'orig.ident', add.rug = TRUE)
ridgeEnrichment(ES2r, gene.set = "ST", group = 'orig.ident', add.rug = TRUE)


levi.combined.markers <- FindAllMarkers(levi.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
levi.combined.markers %>%
group_by(cluster) %>%
slice_max(n = 2, order_by = avg_log2FC)
levi.combined.markers %>%
group_by(cluster) %>%
top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(levi.combined, features = top10$gene) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))


plot <- DotPlot(levi.subset.RGC, features = c('Adgrg1','Adgrl3','Ccr4','Celsr1','Celsr2','Celsr3',
'Cxcr4','Dcc','Drd1','Drd2','Erbb4','Esr2','Fzd3','Gfra3','Gpr173','Il1r1','Nr2f1','Nr2f2',
'Nr4a2','Nsmf','Ntrk2','Ptprz1','Robo1','Robo2','Robo3','Unc5c','Unc5d'),  cols = 'RdBu', dot.scale = 10, split.by = 'orig.ident')
plot + theme(axis.text.x = element_text(angle = 315, family = 'Arial'), axis.text.y = element_text(family = 'Arial'))

FD59  <- LoadH5Seurat('C://Users/Emil/10X/scretina/FD59.h5Seurat')
FD59$stage <- 'FD59'
FD82  <- LoadH5Seurat('C://Users/Emil/10X/scretina/FD82.h5Seurat')
FD82$stage <- 'FD82'
FD125  <- LoadH5Seurat('C://Users/Emil/10X/scretina/FD125.h5Seurat')
FD125$stage <- 'FD125'
adult  <- LoadH5Seurat('C://Users/Emil/10X/scretina/adult.h5Seurat')
adult$stage <- 'Adult'
FD59RGC <- subset(FD59, idents = c('RGC'))
FD82RGC <- subset(FD82, idents = c('RGC'))
FD125RGC <- subset(FD125, idents = c('RGC'))
adultRGC <- subset(adult, idents = c('RGC'))
humanRGC <- merge(FD59RGC, y = c(FD82RGC, FD125RGC, adultRGC))
humanRGC$stage <- factor(humanRGC$stage, levels = c('Adult','FD125','FD82','FD59'))
mp_genes <- lapply(mp_genes, toupper) 
mp_genes <- unlist(mp_genes)
st_genes <- lapply(st_genes, toupper) 
st_genes <- unlist(st_genes)


gene.sets <- list(MP = mp_genes, ST = st_genes)

ESr2 <- enrichIt(obj = humanRGC,
                gene.sets = gene.sets,
                groups = 1000, cores = 8)
humanRGC <- AddMetaData(humanRGC, ESr2)
ES2r2 <- data.frame(humanRGC[[]], Idents(humanRGC))
colnames(ES2r2)[ncol(ES2r2)] <- "cluster"

ridgeEnrichment(ES2r2, gene.set = "MP", group = 'orig.ident', add.rug = TRUE)

ESr3 <- enrichIt(obj = human.combined.RGC,
                 gene.sets = gene.sets,
                 groups = 1000, cores = 8)
human.combined.RGC <- AddMetaData(human.combined.RGC, ESr3)
ESr3 <- data.frame(human.combined.RGC[[]], Idents(human.combined.RGC))
colnames(ESr3)[ncol(ESr3)] <- "cluster"

ridgeEnrichment(ESr3, gene.set = "MP", group = 'orig.ident', add.rug = TRUE)


ht1 <- netAnalysis_signalingRole_heatmap(cellchat_levi.late, pattern = "outgoing", height = 2, signaling = c('SEMA5','MK'))
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_levi.late, pattern = "incoming", height = 2, signaling = c('SEMA5','MK'))

ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat_levi.early, pattern = "outgoing", height = 6, signaling = c('NEGR','LAMININ','NGL','NT','SEMA3','ESAM'))
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_levi.early, pattern = "incoming", height = 6, signaling = c('NEGR','LAMININ','NGL','NT','SEMA3','ESAM'))

ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat_levi.early, pattern = "outgoing", height = 3, width = 3, signaling = c('NEGR','LAMININ','NGL','NT','SEMA3','ESAM'))
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_levi.early, pattern = "incoming", height = 3, width = 3, signaling = c('NEGR','LAMININ','NGL','NT','SEMA3','ESAM'))

ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat_levi.late, pattern = "outgoing", height = 1, width = 3, signaling = c('SEMA5','MK'))
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_levi.late, pattern = "incoming", height = 1, width = 3, signaling = c('SEMA5','MK'))

ht1 + ht2

strwidth <- function(x) {0.5}
netVisual_aggregate(cellchat_levi.early, signaling = c('LAMININ'), layout = "chord", show.legend = FALSE)
netVisual_aggregate(cellchat_levi.early, signaling = c('NT'), layout = "chord", show.legend = FALSE)
