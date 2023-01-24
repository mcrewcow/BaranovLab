library(SingleCellExperiment)


exprMat <- as.matrix(GetAssayData(object = tmp, slot = "data"))
cellInfo <- colnames(x = tmp)
cellInfo <- data.frame(seuratCluster=Idents(tmp))

library(SCENIC)
db='C:/Users/nelso/Downloads/Ristarget/'
list.files(db)
scenicOptions <- initializeScenic(org="mgi", dbDir= db , nCores=10)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
scenicOptions@settings[["nCores"]] <- 1 
runGenie3(exprMat_filtered_log, scenicOptions)

exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

# Optional: Binarize activity
# aucellscenicOptions, aucType="AUC") # choose settings
App <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
# savedSelections <- shiny::runApp(aucellApp)
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

# Export:
# saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
export2loom(scenicOptions, exprMat)

# To save the current status, or any changes in settings, save the object again:
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Exploring output 
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
viewMotifs(tableSubset) 

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 

# Cell-type specific regulators (RSS): 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

library(AUCell)
cellInfo <- data.frame(seuratCluster=Idents(tmp))
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
tfs <- c("E2f1","E2f2","E2f3", "E2f8", "Sox9", "Sox2")
tfs <- c("rna_Hif1a","rna_Arnt","Epas1", "rna_Arnt2", "rna_Hif3a")
tfs <- c("Stat1","Stat3","Stat6", "Nfkb1")
tfs <- c("Sox9","Sox6","Prdm1", "Sox4")
tfs <- c("Lef1","Tcf7","Prdm1", "Sox4")

regulonAUC <- selectRegulons(regulonAUC, tfs, onlyNonDuplicatedExtended = FALSE)

#succcessful
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

# regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
#    function(cells) rowMeans(data.frame(regulon_data[,cells])))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                        cluster_columns =FALSE,column_order = sort(colnames(regulonActivity_byCellType_Scaled)))
library(Seurat)
dr_coords <- Embeddings(CPOC, reduction="tsne")

tfs <- c("Sox10","Irf1","Sox9", "Dlx5")
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "AUC")
