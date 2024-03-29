#Install and load the packages for the analysis

library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reticulate)

library(EnhancedVolcano)
library(DoubletFinder)
library(CytoTRACE)
library(biomaRt)
library(SeuratDisk)
library(SeuratWrappers)
library(Matrix)

library(dittoSeq)

#Download and import the data
data_dir <- 'C:/Users/Emil/BaranovLab/Retina/D59_fetal_filtered_gene_bc_matrices/GRCh38'
list.files(data_dir)
D59fetal <- Read10X(data.dir = data_dir)
D59fetalS = CreateSeuratObject(counts = D59fetal)

#This is important in future during ScaleData
memory.limit(size=56000)

#Import S/G2M cell cycle genes
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

#Create the function for creating and processing the Seurat objects

ProcessSeu <- function(Seurat){
Seurat <- NormalizeData(Seurat)
Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
Seurat <- ScaleData(Seurat) #could be replaced with SCTransform
Seurat <- RunPCA(Seurat)
Seurat <- FindNeighbors(Seurat, dims = 1:10)
Seurat <- FindClusters(Seurat, resolution = 0.5)
Seurat <- RunUMAP(Seurat, dims = 1:10)
Seurat <- RunTSNE(Seurat,  dims.use = 1:10 )
DimPlot(object = Seurat, reduction = "umap")
return (Seurat)
}

#Quantify the RB/MT, cell cycle expression levels and process the data
D59fetalS[["percent.rb"]] <- PercentageFeatureSet(D59fetalS, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
D59fetalS[["percent.mt"]] <- PercentageFeatureSet(D59fetalS, pattern = "^MT-")
D59fetalS <- CellCycleScoring(D59fetalS, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
D59fetalS <- subset(D59fetalS, subset = nCount_RNA > 1000 & nCount_RNA < 20000 & nFeature_RNA > 600 & nFeature_RNA < 6000 & percent.mt < 15 & percent.rb < 30)
D59fetalS <- ScaleData(D59fetalS, verbose = T, vars.to.regress = c('nCount_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score")) 
#If receive the error for nbins during the ScaleData, run NormalizeData before it
D59fetalS <- ProcessSeu(D59fetalS)

#Basic visualisation
FeaturePlot(D59fetalS, features = c('RBPMS'))
DimPlot(D59fetalS, reduction = 'umap', label = TRUE, repel = TRUE)


###Metadata input
#D59metadata <- read.table(file = "clipboard", sep = "\t", header=TRUE)
#rownames(D59metadata) <- D59metadata$barcode
#D59fetalS1 <- AddMetaData(object = D59fetalS1, metadata = D59metadata)
#D59fetalS2 <- subset(D59fetalS2, subset = type != 'NA')
#D59fetalS2 <- RenameIdents(object = D59fetalS2, '0' = 'RGC', '1' = 'T1', '2' = 'RGC','3' = 'Progenitors', '4' = 'T1', '5' = 'Amacrine', '6' = 'RGC', '7' = 'Progenitors', '8' = 'Progenitors')

#To find the DEGs
D59fetalS.markers <- FindAllMarkers(D59fetalS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

D59fetalS.markers %>%
     group_by(cluster) %>%
     slice_max(n = 2, order_by = avg_log2FC)

D59fetalS.markers %>%
     group_by(cluster) %>%
     top_n(n = 10, wt = avg_log2FC) -> top10
 DoHeatmap(D59fetalS, features = top10$gene) + NoLegend()

#For merging P and C datasets
D82PCfetalS <- merge(D82cfetalS , y = D82pfetalS)

#FD59 clusters annotation
#D59fetalS1 <- RenameIdents(object = D59fetalS1, '0' = 'RGC', '1' = 'T1A', '2' = 'RGC','3' = 'Progenitors', '4' = 'T1B', '5' = 'AC/HC', '6' = 'RGC', '7' = 'Progenitors', '8' = 'Progenitors')

ProcessInt(my_data_frame) #This line is to choose the dataset for Integrated Processing in Seurat

ProcessInt <- function(data.integrated){
data.integrated <- ScaleData(data.integrated, verbose = FALSE)
data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = FALSE)
data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
data.integrated <- FindClusters(data.integrated, resolution = 0.5)
data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:20)
data.integrated <- RunTSNE(data.integrated,  dims.use = 1:10 )
}

#Choose the objects for integration

integration_list <- list(poc, pocGH, soc, socGH, pocV, socV)

features <- SelectIntegrationFeatures(object.list = integration_list)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)

data.combined <- IntegrateData(anchorset = data.anchors)

finaldata <- ProcessInt(data.combined)

#For saving/loading the Seurat objects
SaveH5Seurat(D125PCfetalS1, 'C://Users/Emil/10X/scretina/FD125.h5Seurat', overwrite = TRUE)
seuratobject <- LoadH5Seurat('C://Users/Emil/10X/scretina/FD125.h5Seurat')

#For subsetting the clusters / identification
D82PCfetalS11 <- subset(D82PCfetalS, subset = seurat_clusters == 11)
names <- colnames(D82PCfetalS11)
D82PCfetalS <- SetIdent(D82PCfetalS, cells = names, value = 'Amacrine')
DimPlot(D82PCfetalS, reduction = 'umap', label = TRUE, label.box = TRUE)

#For subcluster/expression analysis
HLAexpr <- GetAssayData(adultretinaSsub, assay = 'RNA', slot = 'data')['HLA-DRA',]
poshla <- names(which(HLAexpr>0))
neghla <- names(which(HLAexpr=0))
poscells <- subset(adultretinaSsub, cells = poshla)
negcells <- subset(adultretinaSsub, cells = neghla)
adultretinaS1 <- SetIdent(adultretinaS1, cells = poshla, value = 'Microglia')

library(htmlwidgets)
library(plotly)
library(ggplotly)

#For Violin Plots
#Subset RGCs and merge the data adding the $stage metadata
RGCCXCR <- merge(D45orgSrgc, y = c(D60orgSrgc, D59fetalSRGC, D82PCfetalSRGC, D125PCfetalSRGC, adultretinaSRGC))
RGCCXCR$stage <- factor(RGCCXCR$stage, levels = c('FD59','FD82','FD125','Adult','OD45','OD60'))
library(ggplot2)
vln_df = data.frame(CXCR4 = RGCCXCR[["RNA"]]@data["CXCR4",], cluster = RGCCXCR$stage)
ggplot(vln_df, aes(x = cluster, y = CXCR4)) + geom_violin(aes(fill = cluster), trim=TRUE, scale = "width") + geom_jitter(width = 0.2) + theme_minimal()

#For bars
ggplot(bars, aes(x = Stage, y = CXCR4, fill = Stage)) + geom_bar(stat = 'identity') + theme_minimal()
ggplot(bars, aes(x = Stage, y = CXCR4, color = Stage)) + geom_bar(stat = 'identity', fill = 'white', size = 1.2) +
 geom_text(aes(label = round(CXCR4, digits = 1)), vjust = -1, color = 'black', size = 3.5) + theme_minimal() 
ggplot(bars, aes(x = Stage, y = CXCR4, color = Stage)) + geom_bar(stat = 'identity', fill = 'white', size = 1.2) +
 geom_text(aes(label = Amount), vjust = -1, color = 'black', size = 3.5) + theme_minimal()

#For the dotplots merge the data with the $analysis metadata
#The genes include: CXCR4, DCC, PTCH1, ROBO1 etc
analysis <- RenameIdents(analysis, 'Rod' = 'PR','Cone' = 'PR', 'ON-bipolar' = 'Bipolar', 'OFF-cone bipolar' = 'Bipolar', 'GABA-Amacrine' = 'Amacrine')
CXCR.FD59 <- analysis[['RNA']]@data["CXCR4",] * (analysis$an == "FD59")
CXCR.FD82 <- analysis[['RNA']]@data["CXCR4",] * (analysis$an == "FD82")
CXCR.FD125 <- analysis[['RNA']]@data["CXCR4",] * (analysis$an == "FD125")
CXCR.Adult <- analysis[['RNA']]@data["CXCR4",] * (analysis$an == "Adult")
analysis[['NEW']] <- CreateAssayObject(data = rbind(CXCR.FD59, CXCR.FD82, CXCR.FD125, CXCR.Adult))
DefaultAssay(analysis) <- "NEW"
DotPlot(analysis, features = c("CXCR.FD59", "CXCR.FD82",'CXCR.FD125','CXCR.Adult'), dot.scale = 10) 

analysisrgcs <- subset(analysis, idents = c('RGC'))
#analysisrgcs$an <- factor(RGCCXCR$stage, levels = c('FD59','FD82','FD125','Adult'))
DotPlot(analysisrgcs, features = c('ACKR3','ADGRG1','ADGRL3','CCR4','CELSR1','CELSR2','CELSR3',
'CXCR4','DCC','DRD1','DRD2','ERBB4','ESR2','FGFR1','FZD3','GFRA3','GPR173','IL1R1','NR2F1','NR2F2',
'NR4A2','NSMF','NTRK2','PTCH1','PTPRZ1','ROBO1','ROBO2','ROBO3','UNC5C','UNC5D'), split.by = 'an', cols = c('blue','blue','blue','blue')) 
#Or cols = 'RdBu' for red/blue scale colour
plot <- DotPlot(D59fetalSre06, features = c('ADGRG1','ADGRL3','CCR4','CELSR1','CELSR2','CELSR3',
                                    'CXCR4','DCC','DRD1','DRD2','ERBB4','ESR2','FZD3','GFRA3','GPR173','IL1R1','NR2F1','NR2F2',
                                    'NR4A2','NSMF','NTRK2','PTPRZ1','ROBO1','ROBO2','ROBO3','UNC5C','UNC5D'),  cols = 'RdBu', dot.scale = 20)
plot + theme(axis.text.x = element_text(angle = 315, family = 'Arial'), axis.text.y = element_text(family = 'Arial'))
#For FeatureScatter
RGCCXCR <- merge(D59fetalSRGC, y = c(D82PCfetalSRGC, D125PCfetalSRGC, adultretinaSRGC))
RGCCXCR$stage <- factor(RGCCXCR$stage, levels = c('FD59','FD82','FD125','Adult'))
FeatureScatter(RGCCXCR, feature1 = 'CXCR4', feature2 = 'RBPMS', group.by = 'stage', pt.size = 5)

#For the multiple gene analysis
mp_genes <- read.table(file = "clipboard", sep = "\t", header=TRUE) #Read the file with list of the genes
mp_genes <- as.list(mp_genes)
mp_genes <- lapply(mp_genes, toupper) 
mp_genes <- unlist(mp_genes)
mp <- unique(mp_genes)
DotPlot(analysisrgcs, features = mp, split.by = 'an', cols = c('blue','blue','blue','blue'))

#For RBPMS/POU4F2 populations division
RBPMSexpr <- GetAssayData(D82PCfetalSRGC, assay = 'RNA', slot = 'data')['RBPMS',]
posrbpms <- names(which(RBPMSexpr>0))
poscellsRBPMS <- subset(D82PCfetalSRGC, cells = posrbpms)

POU4F2expr <- GetAssayData(D82PCfetalSRGC, assay = 'RNA', slot = 'data')['POU4F2',]
pospou4f2 <- names(which(POU4F2expr>0))
poscellspou4f2 <- subset(D82PCfetalSRGC, cells = pospou4f2)

#Escape analysis
library(escape)
library(GSEABase)
library(dittoSeq)
gene.sets <- list(Multipolar = c("DCX",     "SLC17A6", "NEUROD1", "ACTR2",   "L1CAM",   "STK25",   "CDKN1B",  "CDK5",   
                                  "NRP1",    "NUAK1",   "KIF21B",  "NHLH1",   "STK24",   "NDEL1",   "PLXNB2",  "RND2",   
                                 "NEUROD2", "UNC5D",   "PLXNA2",  "AXIN1",   "PLXNA4",  "NEUROD6", "NEUROG2", "PLXND1", 
                                  "GJA1",    "HTR6"   ),Somal = c("RTN4",     "PAFAH1B1", "RNF7",     "CDH2",     "ZMIZ1",    "CUL5",     "RELN",    
                                                                  "CRK",      "DIXDC1",   "NCK2",     "ARF6",    "ADGRG1",   "WDR47",    "MBOAT7",  
                                                                  "CTNNB1",   "PIK3CA",   "LRP8",     "NDEL1",    "SOCS7",    "RAPGEF1",  "LAMB1",   
                                                                   "DAB2IP",   "FBXO45",   "SRGAP2",   "LIMK1",    "POU3F2",   "DAB1",     "SOCS3",   
                                                                 "GLI3",     "NR2E1",    "MDGA1",    "FOXG1",    "FUT10",    "CCDC141",  "DISC1",   
                                                                  "POU3F3",   "LHX6",     "COL3A1",   "P2RY12"), Migration = mp)

#mp is http://www.informatics.jax.org/go/term/GO:0001764

#Download the genesets from MSigDb
gene.sets1 <- getGeneSets(library = "C5", gene.sets = "GOBP_NEURON_MIGRATION")
#Perform the matrix enrichment
ES <- enrichIt(obj = allmergedpourbp, 
               gene.sets = gene.sets1, 
               groups = 1000, cores = 4)
ES2 <- data.frame(allmergedpourbp1[[]], Idents(allmergedpourbp1))
colnames(ES2)[ncol(ES2)] <- "cluster"
ES2$stage <- factor(ES2$stage, levels = c('FD59','FD82','FD125','Adult'))
ES2$labeling <- factor(ES2$labeling, levels = c('RBPMS+','POU4F2+'))
ridgeEnrichment(ES2, gene.set = "GOBP_NEURON_MIGRATION", group = "labeling", add.rug = TRUE) + facet_wrap(~stage, ncol = 1)

#allmergedpourbp is the Seurat object with merged RBPMS/POU4F2 subsets of RGCs for FD59,FD82,FD125 and Adult
allmergedpourbp1 <- allmergedpourbp
allmergedpourbp1 <- AddMetaData(allmergedpourbp1, ES)
D82poscellspou4f2 <- RenameIdents(D82poscellspou4f2, 'RGC' = 'POU4F2+')
D82poscellsRBPMS <- RenameIdents(D82poscellsRBPMS, 'RGC' = 'RBPMS+')

#To generate the heatmap
ggplot(migration_allgenes1, aes(x = Gene, y = Timepoint, fill= Amount)) +
geom_tile(color = "white",
lwd = 1.5,linetype = 1) + geom_text(aes(label = round(Amount,1), colours = 'black')) + 
scale_fill_gradientn(colours=c('white','lightblue','blue','darkblue'),guide="colorbar", trans = 'log') + 
facet_grid(~Migrtype, scale = 'free', space = 'free')

