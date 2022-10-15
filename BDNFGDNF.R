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

#load h5seurat files 
FD59  <- LoadH5Seurat('C://Users/Emil/10X/scretina/FD59.h5Seurat')
FD59$stage <- 'FD59'

gg1 <- ggplot(factorsall, aes(x = Stage, y = GDNF, fill = Gene, color = Gene, group = Gene)) +
 geom_point() + geom_xspline() + scale_y_continuous("BDNF", sec.axis = sec_axis(~./75, name = "GDNF")) + theme_bw()
 
gg2 <- gg2 <- ggplot(factorsall, aes(x = Stage, y = GDNF, fill = Gene, color = Gene, group = Gene)) + geom_point() + geom_xspline() + scale_y_continuous("BDNF", sec.axis = sec_axis(~./75, name = "GDNF")) +
 theme_bw() + theme(axis.line=element_blank(),axis.text.x=element_blank(),
+     axis.text.y=element_blank(),axis.ticks=element_blank(),
+     axis.title.x=element_blank(),
+     axis.title.y=element_blank(),legend.position="none",
+     panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
+     panel.grid.minor=element_blank(),plot.background=element_blank())


FD59$labels <- FD59@active.ident
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


pairLR.NT <- extractEnrichedLR(cellchat, signaling = 'NT', geneLR.return = FALSE)
netAnalysis_contribution(cellchat, signaling = pathways.show)
LR.show <- pairLR.NT[1,]
LR.show
netVisual_individual(cellchat, signaling = 'NT', pairLR.use = LR.show, layout = "chord", show.legend = TRUE)
