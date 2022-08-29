RBPMSexpr <- GetAssayData(D125PCfetalSRGC, assay = 'RNA', slot = 'data')['RBPMS',]
posrbpms <- names(which(RBPMSexpr>0))
negrbpms <- names(which(RBPMSexpr == 0))
D125PCfetalS1cellchat <- SetIdent(D125PCfetalS1cellchat, cells = posrbpms, value = 'Mature RGCs')
D125PCfetalS1cellchat  <- SetIdent(D125PCfetalS1cellchat , cells = negrbpms, value = 'Immature RGCs')
DimPlot(D125PCfetalS1cellchat, label = TRUE)
D125PCfetalS1cellchat$labels <- D125PCfetalS1cellchat@active.ident
cellchat <- createCellChat(object = D125PCfetalS1cellchat, group.by = "labels")
CellChatDB <- CellChatDB.human
> CellChatDB.use <- CellChatDB
> cellchat@DB <- CellChatDB.use
> cellchat <- subsetData(cellchat)
> future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
> cellchat <- identifyOverExpressedInteractions(cellchat)
> 
> cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
> cellchat <- computeCommunProbPathway(cellchat)
> cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
> ht1
> ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
> ht2
> ht1+ht2
