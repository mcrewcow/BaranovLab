p <- data.frame(rgcs@meta.data$GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE[rgcs$stage == 'Adult'])
pvis <- ggplot(p, aes(x=rgcs.meta.data.GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE.rgcs.stage....)) + geom_density()
ggplotly(pvis)
cellidspos <- names(which(rgcs$stage == 'Adult' & rgcs$GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE >= 0.45))
cellidsneg <- names(which(rgcs$stage == 'Adult' & rgcs$GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE < 0.45))
rgcs <- SetIdent(rgcs, cells = cellidspos, value = '1')
rgcs <- SetIdent(rgcs, cells = cellidsneg, value = '0')

#Same for FD125,FD82,FD59

rgcs$maint_syn_str <- rgcs@active.ident

cellidspos <- names(which(rgcs$stage == 'Adult' & rgcs$GOBP_MAINTENANCE_OF_PRESYNAPTIC_ACTIVE_ZONE_STRUCTURE >= 0.714))
cellidsneg <- names(which(rgcs$stage == 'Adult' & rgcs$GOBP_MAINTENANCE_OF_PRESYNAPTIC_ACTIVE_ZONE_STRUCTURE < 0.471))
cellidsmid <- names(which(rgcs$stage == 'Adult' & rgcs$GOBP_MAINTENANCE_OF_PRESYNAPTIC_ACTIVE_ZONE_STRUCTURE >= 0.471 & 
rgcs$GOBP_MAINTENANCE_OF_PRESYNAPTIC_ACTIVE_ZONE_STRUCTURE < 0.714))
rgcs <- SetIdent(rgcs, cells = cellidspos, value = '1')
rgcs <- SetIdent(rgcs, cells = cellidsneg, value = '0')
rgcs <- SetIdent(rgcs, cells = cellidsmid, value = '0.5')

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


#Repetative part for DimPlots + markers
rgcsfd125$stat <- paste(rgcsfd125$maint_syn_str, rgcsfd125$maint_postsyn_spec_str, rgcsfd125$maint_presyn_act_zone_str, rgcsfd125$pos_reg_syn_plast, rgcsfd125$postsyn_act_cytosk_org, rgcsfd125$presyn_act_zone_org)

DimPlot(rgcsfd125, reduction = 'umap', label = TRUE, label.box = TRUE, split.by = 'maint_syn_str') + ggtitle('FD125')
DimPlot(rgcsfd125, reduction = 'umap', label = TRUE, label.box = TRUE, split.by = 'maint_postsyn_spec_str') + ggtitle('FD125')
DimPlot(rgcsfd125, reduction = 'umap', label = TRUE, label.box = TRUE, split.by = 'pos_reg_syn_plast') + ggtitle('FD125')
DimPlot(rgcsfd125, reduction = 'umap', label = TRUE, label.box = TRUE, split.by = 'postsyn_act_cytosk_org') + ggtitle('FD125')
DimPlot(rgcsfd125, reduction = 'umap', label = TRUE, label.box = TRUE, split.by = 'presyn_act_zone_org') + ggtitle('FD125')
DimPlot(rgcsfd125, reduction = 'umap', label = TRUE, label.box = TRUE, split.by = 'maint_presyn_act_zone_str') + ggtitle('FD125')

#active.ident line + markers
rgcsfd125@active.ident <- rgcsfd125$maint_syn_str
rgcsfd125@active.ident <- rgcsfd125$maint_presyn_act_zone_str
rgcsfd125@active.ident <- rgcsfd125$maint_postsyn_spec_str
rgcsfd125@active.ident <- rgcsfd125$pos_reg_syn_plast
rgcsfd125@active.ident <- rgcsfd125$postsyn_act_cytosk_org
rgcsfd125@active.ident <- rgcsfd125$presyn_act_zone_org
 markers <- FindAllMarkers(rgcsfd125, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
    
    markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top20
write.csv(markers, 'C://Users/Emil/10X/scretina/markers.csv')



gene.sets1 <- getGeneSets(library = "C5", gene.sets = c('GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE','GOBP_PRESYNAPTIC_ACTIVE_ZONE_ORGANIZATION',
                                                        'GOBP_POSITIVE_REGULATION_OF_SYNAPTIC_PLASTICITY','GOBP_MAINTENANCE_OF_PRESYNAPTIC_ACTIVE_ZONE_STRUCTURE',
                                                        'GOBP_MAINTENANCE_OF_POSTSYNAPTIC_SPECIALIZATION_STRUCTURE',
                                                        'GOBP_POSTSYNAPTIC_ACTIN_CYTOSKELETON_ORGANIZATION'), species = "Mus musculus")

ES <- enrichIt(obj = mousefetalRGC, 
               gene.sets = gene.sets1, 
               groups = 1000, cores = 4, 
               min.size = NULL)
mousefetalRGC <- AddMetaData(mousefetalRGC, ES)
ES2 <- data.frame(mousefetalRGC [[]], Idents(mousefetalRGC ))
colnames(ES2)[ncol(ES2)] <- "cluster"
#after it visualisation from the 1 line

mousefetalRGC$maint_syn_str  <- as.numeric(as.character(mousefetalRGC$maint_syn_str))
mousefetalRGC$presyn_act_zone_org  <- as.numeric(as.character(mousefetalRGC$presyn_act_zone_org))
mousefetalRGC$pos_reg_syn_plast  <- as.numeric(as.character(mousefetalRGC$pos_reg_syn_plast))
mousefetalRGC$maint_presyn_act_zone_str  <- as.numeric(as.character(mousefetalRGC$maint_presyn_act_zone_str))
mousefetalRGC$maint_postsyn_spec_str   <- as.numeric(as.character(mousefetalRGC$maint_postsyn_spec_str))
mousefetalRGC$postsyn_act_cytosk_org   <- as.numeric(as.character(mousefetalRGC$postsyn_act_cytosk_or))
mousefetalRGC$postsyn_act_cytosk_org   <- as.numeric(as.character(mousefetalRGC$postsyn_act_cytosk_org))

mousefetalRGC$stat <- mousefetalRGC$presyn_act_zone_org  + mousefetalRGC$maint_syn_str +
    mousefetalRGC$pos_reg_syn_plast + mousefetalRGC$maint_presyn_act_zone_str +
    mousefetalRGC$maint_postsyn_spec_str + mousefetalRGC$postsyn_act_cytosk_org

mousefetalRGC <- SetIdent(mousefetalRGC, value = mousefetalRGC$stat)
mousefetalRGC <- RenameIdents(mousefetalRGC, '0' = 'low','0.5'='low','1'='low','1.5'='mid','2'='mid','2.5'='mid','3'='mid','3.5'='mid','4'='mid','4.5'='mid','5'='high','5.5'='high','6'='high')
table(mousefetalRGC$stat[mousefetalRGC$stage == 'E14'])

ggplot(my_data1, aes(x = Day, y = Percentage, fill = Amount)) + geom_bar(position="stack", stat="identity") 

#for mouse
gene_list <- c('ADORA1',	'ADRA2A',	'ADRB2',	'AGTR2',	'AHSG',	'AKT1',	'AKT2',	'ALK',	'AR',	'AREG',	'ARNT',	'AXL',	'BDKRB2',
               'CEACAM1',	'BLK',	'BRAF',	'DDR1',	'RUNX2',	'CD4',	'CD7',	'CD8A',	'CD63',	'CDH3',	'CDH13',	'CHN1',	'CHRNA3',
               'AP3S1',	'CSF1R',	'CSPG4',	'CTNNB1',	'DOK1',	'EFNA1',	'EFNA2',	'EFNA3',	'EFNA4',	'EFNA5',	'EFNB1',	'EFNB2',
               'EFNB3',	'EGFR',	'EPHA2',	'EPHA1',	'EPHA3',	'EPHA4',	'EPHA5',	'EPHA7',	'EPHA8',	'EPHB1',	'EPHB2',	'EPHB3',	'EPHB4',
               'EPHB6',	'ERBB2',	'ERBB3',	'ERBB4',	'EFEMP1', 'FER',	'FES',	'FGFR1',	'FGFR3',	'FGFR2',	'FGFR4',	'FGR',	'FLT1',	'FLT3',
               'FLT4',	'FRK',	'FUT7',	'GATA3',	'GFRA1',	'GFRA2',	'GFRA3',	'GHR',	'GHRHR',	'GHSR',	'GPR21',	'GPER1',	'FFAR3',	'GRB2',
               'GRB7',	'GRB10',	'GRB14',	'HCK',	'NRG1',	'NDST1',	'IGF1R',	'IGF2R',	'IGFBP1',	'IGFBP2',	'IGFBP3',	'IGFBP4',	'IGFBP5',
               'IGFBP6',	'RBPJ',	'INSR',	'INSRR',	'IRS1',	'ITGA1',	'ITGA5',	'ITGB3',	'JAK2',	'JAK3',	'KDR',	'KIT',	'STMN1', 'LCK',	'LEP',
               'LRP1',	'LTK',	'LYN',	'MET',	'FOXO4',	'MST1R',	'MUSK',	'NDN',	'NEDD9',	'NGFR',	'NKX3-1',	'NTRK1',	'NTRK2',	'NTRK3',	'ROR1',
               'ROR2',	'DDR2',	'PAK1',	'PAK2',	'PAK3',	'PDGFRA',	'PDGFRL',	'PDGFRB',	'PIGR',	'PIK3C2A',	'PIK3CA',	'PIK3CB',	'PIK3CD',	'PIK3R1',
               'PIK3R2',	'PLAUR',	'PRLR',	'PSEN1',	'PTGIR',	'PTPN1',	'PTPN2',	'PTPN3',	'PTPN11',	'PTPN12',	'PTPRA',	'PTPRE',	'PTPRG',	'PTPRJ',
               'PTPRR',	'RAC1',	'RAF1',	'RARRES2',	'RET',	'ROBO1',	'ROS1',	'RYK', 'SORT1',	'SORL1',	'SOS1',	'SOX9',	'SRC',	'SREBF1',	'SRMS',
               'STAT3',	'STAT5A',	'STAT5B',	'STAT6',	'TEK',	'TIAM1',	'TIE1',	'TYRO3',	'VTN',	'WNT1',	'WNT5A',	'YES1',	'FZD4',	'PIK3R3',
               'SOCS1',	'IRS2',	'TNK1',	'NRP2',	'NRP1',	'SOCS2',	'HAP1',	'MPZL1',	'SOCS3',	'HIP1R',	'GPRC5A',	'FIBP',	'REPS2',	'NOG',	'KL', 'MVP',	
               'NR1H4',	'FGFBP1',	'NAMPT',	'TNK2',	'SPRY3',	'SPRY1',	'SPRY2',	'CNKSR1',	'STUB1',	'EFS',	'MERTK',	'SH2B2',	'TXNIP',	'RGS14',	'NRG3',
               'FRS3',	'FRS2',	'NEU3',	'EMILIN1',	'PTPRT',	'PTP4A3',	'LMTK2',	'SIRT2',	'SETX',	'SIK2',	'SIRT1',	'LEPROTL1',	'FLRT3',	'FLRT2',
               'FLRT1',	'SHC2',	'DSTYK',	'NGEF',	'PRKD2',	'ADGRA2',	'GREM1',	'CYFIP2',	'NPTN',	'SNX5',	'DLL1',		'ADIPOR1',	'CRIM1',
               'GHRL',	'FGFRL1',	'ERRFI1',	'WNT4',	'LEPROT',	'STYK1',	'SMPD3',	'SULF2',	'RTN4',	'SEMA6A',		'NCOA5',	'SNX6',	'SMOC2',
               'GFRA4',	'HHIP',	'GIGYF1',	'GKAP1',	'NDEL1',	'MVB12B',	'OSBPL8',	'RASGRP4',	'NUS1',	'CLNK',	'ZFYVE27',	'SOCS4',	'OTOL1',	'IL31RA',	
               'SAMD10',	'SOGA1',	'FGFBP3',	'SESN3',	'DAB2IP',	'CLEC14A',	'CADM4',	'MUC20',	'DOK6',	'FAM83B',	'STXBP4',	'EPHA10',	'EPHA6',	'LRIT3',
               'GFRAL', 'ADRB2','ANGPT1','ANGPT4','CHRNA3','DGKQ','EFNA5','EGF','GREM1','NRG1','NRG3','PDGFC','PRLR','TAL1', 'AGT','AGTR2',
               'AKT1S1','BCAR1','BDNF','CASP3','CORO1A','CYFIP1','CYFIP2','DDIT4','DOK5','GFRA1','HAP1','KIDINS220','MAGI2','NDN','NGF','NTF3',
               'NTRK1','NTRK2','NTRK3','PPP2R5B','PTPN11','RAF1','RAP1A','RAPGEF1','RAPGEF2','SORT1','SOS1','SPRY1','SPRY2','SRC','TMEM108','WASF1',
               'ZDHHC17','ZFYVE27', 'ADAM17','AGR2','APP','AREG','ARF4','ARTN','ATXN2','BTC','CADM4','CBLC','CCDC88A','CD2AP','CD300LF','CDH5',
               'CNTF','CSF2','CSF3','DAB2IP','ECM1','EFEMP1','EGF','EPGN','ERAP1','ERBB4','EREG','ERN1','ESM1','FAM83B','FER','FGF1','FGF10','FGF16',
               'FGF17','FGF18','FGF2','FGF20','FGF21','FGF22','FGF23','FGF3','FGF4','FGF5','FGF6','FGF7','FGF8','FGF9','FLRT2','FLRT3','FRS2',
               'FRS3','FYN','GATA3','GDNF','GLMN','GRB2','GREM1','HBEGF','HIP1','IL10','IL11','IL12A','IL12B','IL12RB1','IL1A','IL1B','IL1F10','IL1R1',
               'IL1RAP','IL1RN','IL2','IL21','IL23R','IL27RA','IL3','IL4','IL5','IL6','IL6ST','IL7','IL9',
               'IRAK4','ITGA5','ITGB3','JAK2','KL','KLB','LINGO1','LYN','MS4A1','MYD88','NCSTN','NPTN','NRTN','PDCL3','PDGFA','PDGFB','PDGFC','PDGFD',
               'PDGFRA','PDGFRB','PGF','PIBF1','PLSCR1','PSEN1','PSPN','PTPRJ','PYCARD','RNF126','RNF41','SDCBP','SHC1','SLC9A3R1','SNX1','SNX2','SNX4',
               'SOCS5','SRC','TGFA','TIMM50','TLR5','TLR9','TNK2','TOLLIP','TRIP6','TSLP','VAV2','VAV3','VEGFA','VEGFB','VEGFC', 'AGT','AGTR2',
               'AKT1S1','BCAR1','BDNF','CASP3','CORO1A','CYFIP1','CYFIP2','DDIT4','DOK5','GFRA1','HAP1','KIDINS220','MAGI2','NDN','NGF','NTF3',
               'NTRK1','NTRK2','NTRK3','PPP2R5B','PTPN11','RAF1','RAP1A','RAPGEF1','RAPGEF2','SORT1','SOS1','SPRY1','SPRY2','SRC','TMEM108','WASF1',
               'ZDHHC17','ZFYVE27', 'ARF6', 'RAB3A', 'CTBP2', 'ERC1', 'PCLO', 'DBN1', 'DBNL', 'INA', 'RAB3A','CBLN2','CTBP2','DBN1','ARF6','TAGLN3',
               'ERC1','CTBP2','PCLO','MALAT1', 'DNAJC6','STAT3','RALBP1',
             #from here after cellchat 
               'Ncam1', 'L1cam', 'Nrxn1', 'Nlgn1', 'Mpzl1', 'Efna5', 'Nrxn3', 'Nrxn2', 'Ephb2', 'Cntn2', 'Cntnap2', 'Mstn', 'Tgfbr1', 'Acvr2b',
               #from here after individual pathway DEGs 
               'Rab3a','Sdf4','Sncg','Tubb2a','Mapt','Meg3','Ly6h','Uchl1','Stmn3','Pcsk1n','Ina','Stmn4','Syt4','Nefl','Snrpn','Aplp1',
               'S100a10','Rab6b','Tppp3','Apbb2','Rbpms','Celf4','Rac3','Map1lc3a','Islr2','Nsg1','Calb2','Snap25','Runx1t1','Mgst3','Pou4f1','Gng3','Gadd45a', 
               'Erc1','Ly6h','Meg3','Sncg','Tubb2a','Snhg11','Mapt','Pcsk1n','Snrpn','Nell2','Rab6b','Syt4','Celf4','Rbpms','Islr2','Apbb2','Tppp3','Map1lc3a','Calb2',
               'Dbn1','Cplx2','Islr2','Ap2m1','Oaz2','Snrpn','Mtch1','Prelid1','Fscn1','Map1lc3a','Uchl1','Arhgdia','Rab3a','Cotl1','Ctxn1',
               'Erc1','Erc2','Ly6h','Mapt','Sncg','Tubb2a','Snhg11','Meg3','Pcsk1n','Snrpn','Rab6b','Islr2','Syt4','Rbpms','Tppp3','Apbb2','Map1lc3a','Celf4','Calb2',
               'Arf6','Snhg11','Nell2','Ly6h','Snrpn','Rab3a','Meg3','Pcsk1n','Calb2',
               'Dbnl','Actg1','Dbn1','Ina','Rab3a','Islr2','Snrpn','Uchl1','Tubb2a','Thra','Sncg','Mapt','Pcsk1n','Meg3','Map1lc3a','Tppp3','Ly6h','Apbb2','Calb2',
               'Snhg11','Gadd45a','Mdk','Prdx1','Gm10260','Hes6','H2afv','Dlx2','Hes5','Btg2','Mfng','Btbd17',
               'Rassf4','Ccnd1','Itm2a','Hmgn3','Tead2','Alyref','Sfrp2','Fos')

#for human

gene_list <- c('ADORA1',	'ADRA2A',	'ADRB2',	'AGTR2',	'AHSG',	'AKT1',	'AKT2',	'ALK',	'FASLG',	'AR',	'AREG',	'ARNT',	'AXL',	'BDKRB2',
               'CEACAM1',	'BLK',	'BRAF',	'DDR1',	'RUNX2',	'CD4',	'CD7',	'CD8A',	'CD8B',	'CD63',	'CDH3',	'CDH13',	'CHN1',	'CHRNA3',
               'AP3S1',	'CSF1R',	'CSPG4',	'CTNNB1',	'DOK1',	'EFNA1',	'EFNA2',	'EFNA3',	'EFNA4',	'EFNA5',	'EFNB1',	'EFNB2',
               'EFNB3',	'EGFR',	'EPHA2',	'EPHA1',	'EPHA3',	'EPHA4',	'EPHA5',	'EPHA7',	'EPHA8',	'EPHB1',	'EPHB2',	'EPHB3',	'EPHB4',
               'EPHB6',	'ERBB2',	'ERBB3',	'ERBB4',	'EFEMP1', 'FER',	'FES',	'FGFR1',	'FGFR3',	'FGFR2',	'FGFR4',	'FGR',	'FLT1',	'FLT3',
               'FLT4',	'FRK',	'FUT7',	'GATA3',	'GFRA1',	'GFRA2',	'GFRA3',	'GHR',	'GHRHR',	'GHSR',	'GPR21',	'GPER1',	'FFAR3',	'GRB2',
               'GRB7',	'GRB10',	'GRB14',	'HCK',	'NRG1',	'NDST1',	'IGF1R',	'IGF2R',	'IGFBP1',	'IGFBP2',	'IGFBP3',	'IGFBP4',	'IGFBP5',
               'IGFBP6',	'RBPJ',	'INSR',	'INSRR',	'IRS1',	'ITGA1',	'ITGA5',	'ITGB3',	'JAK2',	'JAK3',	'KDR',	'KIT',	'STMN1', 'LCK',	'LEP',
               'LRP1',	'LTK',	'LYN',	'MET',	'FOXO4',	'MST1R',	'MUSK',	'NDN',	'NEDD9',	'NGFR',	'NKX3-1',	'NTRK1',	'NTRK2',	'NTRK3',	'ROR1',
               'ROR2',	'DDR2',	'PAK1',	'PAK2',	'PAK3',	'PDGFRA',	'PDGFRL',	'PDGFRB',	'PIGR',	'PIK3C2A',	'PIK3CA',	'PIK3CB',	'PIK3CD',	'PIK3R1',
               'PIK3R2',	'PLAUR',	'PRLR',	'PSEN1',	'PTGIR',	'PTPN1',	'PTPN2',	'PTPN3',	'PTPN11',	'PTPN12',	'PTPRA',	'PTPRE',	'PTPRG',	'PTPRJ',
               'PTPRR',	'RAC1',	'RAF1',	'RARRES2',	'RET',	'ROBO1',	'ROS1',	'RYK', 'SORT1',	'SORL1',	'SOS1',	'SOX9',	'SRC',	'SREBF1',	'SRMS',
               'STAT3',	'STAT5A',	'STAT5B',	'STAT6',	'TEK',	'TIAM1',	'TIE1',	'TYRO3',	'VTN',	'WNT1',	'WNT5A',	'YES1',	'RAB7A',	'FZD4',	'PIK3R3',
               'SOCS1',	'IRS2',	'TNK1',	'NRP2',	'NRP1',	'SOCS2',	'HAP1',	'MPZL1',	'SOCS3',	'HIP1R',	'GPRC5A',	'FIBP',	'REPS2',	'NOG',	'KL', 'MVP',	
               'NR1H4',	'FGFBP1',	'NAMPT',	'TNK2',	'SPRY3',	'SPRY1',	'SPRY2',	'CNKSR1',	'STUB1',	'EFS',	'MERTK',	'SH2B2',	'TXNIP',	'RGS14',	'NRG3',
               'FRS3',	'FRS2',	'NEU3',	'EMILIN1',	'PTPRT',	'PTP4A3',	'LMTK2',	'SIRT2',	'SETX',	'SIK2',	'ANKS1A',	'SIRT1',	'LEPROTL1',	'FLRT3',	'FLRT2',
               'FLRT1',	'SHC2',	'DSTYK',	'NGEF',	'PRKD2',	'ADGRA2',	'GREM1',	'CYFIP2',	'DNAI1',	'NPTN',	'SNX5',	'DLL1',	'PILRB',	'ADIPOR1',	'CRIM1',
               'GHRL',	'FGFRL1',	'ERRFI1',	'WNT4',	'LEPROT',	'STYK1',	'SMPD3',	'SULF2',	'RTN4',	'SEMA6A',	'NCOA5',	'SNX6',	'SMOC2',
               'GFRA4',	'HHIP',	'GIGYF1',	'GKAP1',	'NDEL1',	'MVB12B',	'OSBPL8',	'RASGRP4',	'NUS1',	'CLNK',	'ZFYVE27',	'SOCS4',	'OTOL1',	'IL31RA',	
               'SAMD10',	'SOGA1',	'FGFBP3',	'SESN3',	'DAB2IP',	'CLEC14A',	'CADM4',	'MUC20',	'DOK6',	'FAM83B',	'STXBP4',	'EPHA10',	'EPHA6',	'LRIT3',
               'GFRAL', 'ADRB2','ANGPT1','ANGPT4','CHRNA3','DGKQ','EFNA5','EGF','GREM1','NRG1','NRG3','PDGFC','PILRB','PRLR','TAL1', 'AGT','AGTR2',
               'AKT1S1','BCAR1','BDNF','CASP3','CORO1A','CYFIP1','CYFIP2','DDIT4','DOK5','GFRA1','HAP1','KIDINS220','MAGI2','NDN','NGF','NTF3','NTF4',
               'NTRK1','NTRK2','NTRK3','PPP2R5B','PTPN11','RAF1','RAP1A','RAPGEF1','RAPGEF2','SORT1','SOS1','SPRY1','SPRY2','SRC','TMEM108','WASF1',
               'ZDHHC17','ZFYVE27', 'ADAM17','AGR2','APP','AREG','ARF4','ARTN','ATXN2','BTC','CADM4','CBLC','CCDC88A','CD2AP','CD300LF','CDH5',
               'CNTF','CSF2','CSF3','DAB2IP','ECM1','EFEMP1','EGF','EPGN','ERAP1','ERBB4','EREG','ERN1','ESM1','FAM83B','FER','FGF1','FGF10','FGF16',
               'FGF17','FGF18','FGF19','FGF2','FGF20','FGF21','FGF22','FGF23','FGF3','FGF4','FGF5','FGF6','FGF7','FGF8','FGF9','FLRT2','FLRT3','FRS2',
               'FRS3','FYN','GATA3','GDNF','GLMN','GRB2','GREM1','HBEGF','HIP1','IL10','IL11','IL12A','IL12B','IL12RB1','IL1A','IL1B','IL1F10','IL1R1',
               'IL1RAP','IL1RN','IL2','IL21','IL23R','IL27RA','IL3','IL36A','IL36B','IL36G','IL36RN','IL37','IL4','IL5','IL6','IL6R','IL6ST','IL7','IL9',
               'IRAK4','ITGA5','ITGB3','JAK2','KL','KLB','LINGO1','LYN','MS4A1','MYD88','NCSTN','NPTN','NRTN','PDCL3','PDGFA','PDGFB','PDGFC','PDGFD',
               'PDGFRA','PDGFRB','PGF','PIBF1','PLSCR1','PSEN1','PSPN','PTPRJ','PYCARD','RNF126','RNF41','SDCBP','SHC1','SLC9A3R1','SNX1','SNX2','SNX4',
               'SOCS5','SRC','TGFA','TIMM50','TLR5','TLR9','TNK2','TOLLIP','TRIP6','TSLP','VAV2','VAV3','VEGFA','VEGFB','VEGFC', 'AGT','AGTR2',
               'AKT1S1','BCAR1','BDNF','CASP3','CORO1A','CYFIP1','CYFIP2','DDIT4','DOK5','GFRA1','HAP1','KIDINS220','MAGI2','NDN','NGF','NTF3','NTF4',
               'NTRK1','NTRK2','NTRK3','PPP2R5B','PTPN11','RAF1','RAP1A','RAPGEF1','RAPGEF2','SORT1','SOS1','SPRY1','SPRY2','SRC','TMEM108','WASF1',
               'ZDHHC17','ZFYVE27', 'ARF6', 'RAB3A', 'CTBP2', 'ERC1', 'PCLO', 'DBN1', 'DBNL', 'INA', 'RAB3A','CBLN2','CTBP2','DBN1','ARF6','TAGLN3',
               'ERC1','CTBP2','PCLO','MALAT1','MT-CYB', 'DNAJC6','SAMD4A','STAT3','RALBP1',
             #from here after cellchat 
               'MDK', 'NCL', 'CADM1', 'NCAM1', 'PTN', 'PTPRZ1', 'ITGA6', 'ITGB1', 'ITGA4', 'L1CAM', 'SPP1', 'ITGA8', 'ITGA5', 'ITGA4', 'CNTN1', 'NRCAM',
               'SEMA6A', 'PLXNA2',
               'ITGAV', 'FGFR1', 'DLL3', 'ITGB5', 'IGF1', 'IGF1R', 'EPHB2', 'EFNB3', 

               #from here after individual pathway DEGs 
               'RAB3A','CBLN2','RPS4Y1',
               'CPLX2','RPS4Y1','DBN1','CLEC2D','EIF5A','RPL17','AC009501.4','RPS10','TUBB3','NBEAL1','MIF','FKBP1A','ERC1','CTBP2','PCLO','RPS4Y1','ATP6V1G2','FKBP1B',
               'ARF6', 'RPS4Y1','FARP1','MEG3','NEFL','APLP2','APP','INA','SNAP25','LY6H','ATP1B1','ARL6IP5','PRNP','ITM2C','SERPINE2','GABBR2','PCP4',
               'REEP5','SPTBN1','PIK3R1','NEUROD1','FABP7','ERC1','CTBP2','PCLO','RPS4Y1','ATP6V1G2','FKBP1B')


)

gene_list <- as.list(gene_list)
gene_list  <- lapply(gene_list , tolower) #or toupper and skip next two lines if human data
library(stringr)
gene_list <- str_to_title(gene_list) 
gene_list <- unlist(gene_list)
gene_list  <- unique(gene_list)
mousefetalRGC$score <- mousefetalRGC@active.ident
avexpr <- AverageExpression(mousefetalRGC, features = gene_list, assays = 'RNA', group.by = c('stage','score'))


avexprfun <- function(average_expression_table) {
    average_expression_table <- as.data.frame(average_expression_table)
    average_expression_table$barcode <- ''
    for(i in 1:length(rownames(avexpr$RNA))) {
        if(average_expression_table[i,1] > average_expression_table[i,2] &
 average_expression_table[i,2] > average_expression_table[i,3] | average_expression_table[i,1] 
< average_expression_table[i,2] & average_expression_table[i,2] < average_expression_table[i,3]) {average_expression_table$barcode[i] <- 'E14'}
        if(average_expression_table[i,4] > average_expression_table[i,5] &
 average_expression_table[i,5] > average_expression_table[i,6] | average_expression_table[i,4]
 < average_expression_table[i,5] & average_expression_table[i,5] < average_expression_table[i,6]) {average_expression_table$barcode[i] <- paste(average_expression_table$barcode[i], 'E16', sep = ' ')}
        if(average_expression_table[i,7] > average_expression_table[i,8] &
 average_expression_table[i,8] > average_expression_table[i,9] | average_expression_table[i,7] <
 average_expression_table[i,8] & average_expression_table[i,8] < average_expression_table[i,9]) {average_expression_table$barcode[i] <- paste(average_expression_table$barcode[i], 'E18', sep = ' ')}
    }
    return(average_expression_table)
}


test <- avexprfun(avtest)
table(test$barcode)
test$genes <- rownames(test)
barcodes <- test[, c("genes", "barcode")]
barcodes <- barcodes[!(is.na(barcodes$barcode) | barcodes$barcode==""), ]

install.packages('VennDiagram')               

library('VennDiagram')  

install.packages('venneuler')

library('venneuler')

ggplot(my_data, aes(x = barcode, y = SHIFT, label = GENE)) + 
geom_text_repel(data=my_data, aes(label=GENE),position = position_jitter(seed = 1),
                color = c('blue','blue','gray','gray','red','blue','blue','blue','gray','gray',
                          'blue','gray','blue','gray','gray','gray','blue','blue','gray','gray','gray',
                          'blue','blue','blue','blue','blue','blue','red','gray','red','red','red','gray','gray')) +
geom_jitter(position = position_jitter(seed = 1)) + geom_hline(yintercept=0) + geom_hline(yintercept=50, linetype="dashed", color = "red") + 
geom_hline(yintercept=-50, linetype="dashed", color = "red") + theme_bw()
