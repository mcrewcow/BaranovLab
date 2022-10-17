#load the packages for the analysis
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
library(ggalt) #this one we need for geom_xspline() - smooth version of geom_line()
library(escape) #GSEA analysis for BDNF/GDNF-apoptosis correlation
#load h5seurat files 
FD59  <- LoadH5Seurat('C://Users/Emil/10X/scretina/FD59.h5Seurat')
FD59$stage <- 'FD59' #add stage label

#generate the dual Y axis plot for GDNF-BDNF expression
gg1 <- ggplot(factorsall, aes(x = Stage, y = GDNF, fill = Gene, color = Gene, group = Gene)) +
 geom_point() + geom_xspline() + scale_y_continuous("BDNF", sec.axis = sec_axis(~./75, name = "GDNF")) + theme_bw()
 
gg2 <- gg2 <- ggplot(factorsall, aes(x = Stage, y = GDNF, fill = Gene, color = Gene, group = Gene)) + geom_point() + geom_xspline() + scale_y_continuous("BDNF", sec.axis = sec_axis(~./75, name = "GDNF")) +
 theme_bw() + theme(axis.line=element_blank(),axis.text.x=element_blank(),
+     axis.text.y=element_blank(),axis.ticks=element_blank(),
+     axis.title.x=element_blank(),
+     axis.title.y=element_blank(),legend.position="none",
+     panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
+     panel.grid.minor=element_blank(),plot.background=element_blank())
gg1 + gg2

#cellchat analysis
FD59$labels <- FD59@active.ident #save current cluster annotation
cellchat <- createCellChat(object = FD59, group.by = "labels")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0, raw.use = FALSE, population.size = TRUE)

cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

computeAveExpr(cellchat, features = c("BDNF","CNTF",'GDNF'), type =  "truncatedMean", trim = 0)
pairLR.NT <- extractEnrichedLR(cellchat, signaling = 'NT', geneLR.return = FALSE) 
netAnalysis_contribution(cellchat, signaling = pathways.show) #check the pathway for L-R pairs and their contribution
LR.show <- pairLR.NT[1,] #choose the L-R part needed
LR.show
netVisual_individual(cellchat, signaling = 'NT', pairLR.use = LR.show, layout = "chord", show.legend = TRUE) #generate the chord plot


ESCAPE
library(escape)
RGCfetal <- merge(FD59RGC, y = c(FD82RGC, FD125RGC))
ES <- enrichIt(obj = RGCfetal,
gene.sets = gene.sets1,
groups = 1000, cores = 4)
RGCfetalescape <- AddMetaData(RGCfetal, ES)
ES2 <- data.frame(RGCfetalescape [[]], Idents(RGCfetalescape))
colnames(ES2)[ncol(ES2)] <- "cluster"
ES2$stage <- factor(ES2$stage, levels = c('FD59','FD82','FD125'))
fetalescapefd59 <- subset(fetalescape, subset = stage == 'FD59')
fetalescapefd82 <- subset(fetalescape, subset = stage == 'FD82')
fetalescapefd125 <- subset(fetalescape, subset = stage == 'FD125')

genessurvival <- c('GDNF','BDNF','NTRK2','GFRA1')

pathways <- c('HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY','HALLMARK_APOPTOSIS','HALLMARK_DNA_REPAIR')

survivalfun <- function(escape_object) {
    for(i in 1:length(pathways)) {
        for(j in 1:length(genessurvival)) {
            p <- FeatureScatter(escape_object,
                                feature1 = genessurvival[j], feature2 = pathways[i], group.by = 'stage')
            ifelse(i == 1 & j == 1, ggsum <- p, ggsum <- ggsum + p) } }
    print(ggsum) }

survivalfun(fetalescapefd59)
survivalfun(fetalescapefd82)
survivalfun(fetalescapefd125)
          
fetalescapefd59GFRA1 <- GetAssayData(fetalescapefd59, assay = 'RNA', slot = 'data')['GFRA1',]
poshla <- names(which(fetalescapefd59GFRA1>0))
poscells <- subset(fetalescapefd59, cells = poshla)
FeatureScatter(poscells,
               feature1 = 'GFRA1', feature2 = 'HALLMARK_APOPTOSIS', group.by = 'stage', cols = 'Orange', span = TRUE, pt.size = 2)

gg1 <- FeatureScatter(poscells,
               feature1 = 'GFRA1', feature2 = 'HALLMARK_APOPTOSIS', group.by = 'stage', cols = 'Orange', pt.size = 2) + ylim(0,3500) + xlim(0.5,3)

gg1 + stat_smooth(method = "lm",
                formula = y ~ x,
                geom = "smooth", color = 'Orange', size = 2)
#save the files after analysis: cellchat and escape
SaveH5Seurat(FD59cellchatrawusefalsepopsizetrue, 'C://Users/Emil/10X/scretina/FD59cellchat.h5Seurat', overwrite = TRUE)
saveRDS(FD59cellchatrawusefalsepopsizetrue, file = "C://Users/Emil/10X/scretina/FD59cellchat.rds")
