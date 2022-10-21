low <- function(x) {
    0 + runif(1, 0, 0.7)
}

high <- function(x) {
    2.3 + runif(1, 0, 0.7)
}

up <- function(x) {
    exp(x/2)-1 + runif(1, -0.4, 0.4)
}

down <- function(x) {
    1/x + runif(1, 0, 0.4)
}

library(ggplot2)
library(patchwork)


AA <- c('AA','Aa','aA','aa')
dim(AA) <- c(1,4)
BB <- c('BB','Bb','bB','bb')
dim(BB) <- c(4,1)
AABB <- cbind(rep(AA,each=nrow(BB)),rep(BB,nrow(AA)))
write.csv(AABB, 'C://Users/Emil/10X/matric.csv')
Excel =B2&""&B3
CC <- c('CC','Cc','cC','cc')
dim(CC) <- c(1,4)
DD <- c('DD','Dd','dD','dd')
dim(DD) <- c(4,1)
CCDD <- rbind(rep(CC,each=nrow(DD)),rep(DD,nrow(CC)))
write.csv(CCDD, 'C://Users/Emil/10X/matric1.csv')
Excel =B2&""&B3
AABB <- read.table(file = "clipboard", sep = "\t", header=FALSE)
dim(AABB) <- c(1,16)
AABB <- as.matrix(AABB)
CCDD <- read.table(file = "clipboard", sep = "\t", header=FALSE)
CCDD <- as.matrix(CCDD)
dim(CCDD) <- c(16,1)
abcd <- rbind(rep(AABB,each=nrow(CCDD)),rep(CCDD,nrow(AABB)))
write.csv(abcd, 'C://Users/Emil/10X/matricall.csv')
Excel =B2&""&B3
matabcd <- read.table(file = "clipboard", sep = "\t", header=FALSE)

tabfun <- function(matgenerated) {
    for(i in 1:length(matgenerated)) {
        tmp <- matgenerated[i]
        tmp  <- unlist(strsplit(as.character(tmp), ""))
        ifelse(tmp[1] == 'A',ifelse(tmp[2] == 'A',function1 <- high, function1 <- down),
               ifelse(tmp[2] == 'A',function1 <- up, function1 <- low)) 
        ifelse(tmp[3] == 'B',ifelse(tmp[4] == 'B',function2 <- high, function2 <- down),
               ifelse(tmp[4] == 'B',function2 <- up, function2 <- low)) 
        ifelse(tmp[5] == 'C',ifelse(tmp[6] == 'C',function3 <- high, function3 <- down),
               ifelse(tmp[6] == 'C',function3 <- up, function3 <- low)) 
        ifelse(tmp[7] == 'D',ifelse(tmp[8] == 'D',function4 <- high, function4 <- down),
               ifelse(tmp[8] == 'D',function4 <- up, function4 <- low)) 
        
        p <- ggplot(data = data.frame(x = c(0,pi)), mapping = aes(x = x))
        
        p <- p +
            stat_function(fun=function1,aes(colour="RGC Receptor"), size=1.5) +
            stat_function(fun=function2,aes(colour="Retina Ligand"), size=1.5) +
            stat_function(fun=function3,aes(colour="Retina Receptor"), size=1.5) +
            stat_function(fun=function4, aes(colour='RGC Ligand'), size=1.5) +
            ylim(0,3) + theme_bw()
         ifelse(i == 1, ggsum <- p +xlab('Age') + ylab('Expression') +
                   theme(axis.text.y=element_blank(),
                         axis.text.x=element_blank()), ggsum <- ggsum + p +xlab('Age') + ylab('Expression') +
                   theme(axis.text.y=element_blank(),
                         axis.text.x=element_blank()) + plot_layout(guides = "collect") & theme(legend.position = "bottom")) } 
    print(ggsum)}




CellChatDB <- CellChatDB.human
interaction_input <- CellChatDB$interaction
library(dplyr)
library(tidyr)
int1 <- interaction_input %>% separate(interaction_name, c('Ligand','R1','R2'))

int2 <- int1[c(1:3)]
int2 <- head(int2, - 4) 
lrpairs <- function(tablet) {
    for(i in 1:length(tablet$Ligand)) {
        if(is.na(tablet$R2[i]) == FALSE) {
            tablet <- rbind(tablet, data.frame(Ligand = tablet$Ligand[i], R1 = tablet$R2[i], R2 = tablet$R2[i])) } }
    return(tablet) }

int2 <- lrpairs(int2)
df_added <- data.frame(Ligand = c('ITGA4','ITGA9','ITGB1','ITGB7','VSIR'), R1 = c('VCAM1','VCAM1','VCAM1','VCAM1','IGSF11'))
int3 <- int2[c(1:2)]
int3 <- rbind(int3, df_added)	 
human <- merge(FD59, y = c(FD82, FD125, adult))
humanRGC <- subset(human, idents = c('RGC'))
human$labels <- human@active.ident
DefaultAssay(humanRGC) <- 'RNA'
p <- AverageExpression(humanRGC, features = int3$Ligand, assays = 'RNA', group.by = 'stage')
p <- as.data.frame(p)
int3 <- int3 %>% distinct(Ligand, R1, .keep_all = TRUE)

exprassign <- function(tabletexpr, tabletall) {
    
    for(i in 1:length(tabletexpr$RNA.Adult)) {
        for(j in 1:length(tabletall$Ligand)) {
            if(tabletall$Ligand[j] == rownames(tabletexpr)[i]) {
                tabletall[j,4] <- tabletexpr[i,1]
                tabletall[j,5] <- tabletexpr[i,2]
                tabletall[j,6] <- tabletexpr[i,3]
                tabletall[j,7] <- tabletexpr[i,4]
            } } }
    return(tabletall) } 
int4$LigandRGCAdult <- 'NA'
int4$LigandRGCFD125 <- 'NA'
int4$LigandRGCFD59 <- 'NA'
int4$LigandRGCFD82 <- 'NA'
int4 <- exprassign(p, int3)
	
p <- AverageExpression(humanRGC, features = int3$R1, assays = 'RNA', group.by = 'stage')
p <- as.data.frame(p)
int4$RecRGCAdult <- 'NA'
int4$RecRGCFD125 <- 'NA'
int4$RecRGCFD59 <- 'NA'
int4$RecRGCFD82 <- 'NA'
exprassign1 <- function(tabletexpr, tabletall, tabletfinal) {
    
    for(i in 1:length(tabletexpr$RNA.Adult)) {
        for(j in 1:length(tabletall$R1)) {
            if(tabletall$R1[j] == rownames(tabletexpr)[i]) {
                tabletfinal[j,8] <- tabletexpr[i,1]
                tabletfinal[j,9] <- tabletexpr[i,2]
                tabletfinal[j,10] <- tabletexpr[i,3]
                tabletfinal[j,11] <- tabletexpr[i,4]
            } } }
    return(tabletfinal) } 
int4 <- exprassign1(p, int3, int4)

p <- AverageExpression(human, features = int3$Ligand, assays = 'RNA', group.by = 'stage')
p <- as.data.frame(p)

int4$LigandAdult <- 'NA'
int4$LigandFD125 <- 'NA'
int4$LigandFD59 <- 'NA'
int4$LigandFD82 <- 'NA'

exprassign2 <- function(tabletexpr, tabletall, tabletfinal) {
    
    for(i in 1:length(tabletexpr$RNA.Adult)) {
        for(j in 1:length(tabletall$Ligand)) {
            if(tabletall$Ligand[j] == rownames(tabletexpr)[i]) {
                tabletfinal[j,12] <- tabletexpr[i,1]
                tabletfinal[j,13] <- tabletexpr[i,2]
                tabletfinal[j,14] <- tabletexpr[i,3]
                tabletfinal[j,15] <- tabletexpr[i,4]
            } } }
    return(tabletfinal) } 
int4 <- exprassign2(p, int3, int4)

p <- AverageExpression(human, features = int3$R1, assays = 'RNA', group.by = 'stage')
p <- as.data.frame(p)

int4$RecAdult <- 'NA'
int4$RecFD125 <- 'NA'
int4$RecFD59 <- 'NA'
int4$RecFD82 <- 'NA'

exprassign3 <- function(tabletexpr, tabletall, tabletfinal) {
    
    for(i in 1:length(tabletexpr$RNA.Adult)) {
        for(j in 1:length(tabletall$R1)) {
            if(tabletall$R1[j] == rownames(tabletexpr)[i]) {
                tabletfinal[j,16] <- tabletexpr[i,1]
                tabletfinal[j,17] <- tabletexpr[i,2]
                tabletfinal[j,18] <- tabletexpr[i,3]
                tabletfinal[j,19] <- tabletexpr[i,4]
            } } }
    return(tabletfinal) } 
int4 <- exprassign3(p, int3, int4)
