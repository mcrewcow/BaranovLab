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
        ifelse(i == 1, ggsum <- p, ggsum <- ggsum + p + plot_layout(guides = "collect") & theme(legend.position = "bottom")) } 
    print(ggsum)}
