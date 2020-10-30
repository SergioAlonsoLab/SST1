library(magrittr)



try(setwd("/imppc/labs/mplab/salonso/SST1/"))
source("/imppc/labs/mplab/salonso/AnalysisFCS/src/analyzeFCS.R")


fcsSet <- read.flowSet(path = "data/Ciclo Celular DAPI/",pattern = "*.fcs")

# 149 experiments
length(fcsSet)

fcsSet <- as.list(fcsSet@frames)

data.frame(filename=names(fcsSet),
           tubeName=sapply(fcsSet,function(x) x@description$`TUBE NAME`),
           date=sapply(fcsSet,function(x) x@description$`$DATE`),
           events=sapply(fcsSet,function(x) nrow(exprs(x)))) -> foo


foo$betterDate <- as.Date(foo$date,format=c("%d-%B-%y"))

write.table(foo,file = "sandbox/fcs.tsv",row.names=F,quote=F,sep="\t")

# pdf("sandbox/fcs_graphs.pdf",20,30)
# par(mfrow=c(6,4))
# for(i in order(foo$betterDate)) {
#   
#   x <- exprs(fcsSet[[i]]) %>% as.data.frame
#   name <- paste(fcsSet[[i]]@description$`TUBE NAME`,fcsSet[[i]]@description$`$DATE`,sep="\n")
#   mySmoothScatter(x$`FSC-A`,x$`SSC-A`)
#   title(name,xlab="FSC-A",ylab="SSC-A",cex.main=2)
#   polygon(c(4e4,26e4,26e4,4e4),c(1e3,1e3,26e4,26e4),col="#FFFF0033")
#   mtext(sprintf("n=%i",nrow(x)),3,-2)
#   x <- subset(x,`FSC-A` > 4e4 & `FSC-A` < 26e4 & `SSC-A` > 1e3 & `SSC-A` < 26e4)
#   mySmoothScatter(x$`Pacific Blue-A`,x$`Pacific Blue-W`,ylim=range(x$`Pacific Blue-W`)*c(0,1))
#   title(name,xlab="Pacific Blue-A",ylab="Pacific Blue-W",cex.main=2)
#   mtext(sprintf("n=%i",nrow(x)),3,-2)
# }
# dev.off()

# the experiments for cell cycle of LS174T clones were run on 27-MAY-2019

fcsSet <- fcsSet[which(foo$date=="27-MAY-2019")]


for(i in 1:length(fcsSet)) {
  
  x <- exprs(fcsSet[[i]]) %>% as.data.frame
  
  #mySmoothScatter(x[,"FSC-A"],x[,"SSC-A"])
  #title(names(fcsSet)[i])
  
  x <- subset(x,`FSC-A` > 4e4 & `FSC-A` < 26e4 & `SSC-A` < 26e4)
  exprs(fcsSet[[i]]) <- as.matrix(x)
  
}


e0 <- lapply(fcsSet,function(x) as.data.frame(exprs(x)))

names(e0) <- names(fcsSet)

e0 <- e0[2]


lapply(e0,function(x) {
  gate0 <- gating(x$`Pacific Blue-A`,x$`Pacific Blue-W`)
  G0 <- locator(n=1)
  
  x1 <- x[gate0,]
  mix2 <- mixtools::normalmixEM(x1$`Pacific Blue-A`,mu=G0$x*c(1,1.8))
  
  gate1 <- gate2 <- rep(F,length(gate0))
  gate1[gate0] <- mix2$posterior[,1] > .5 | x1$`Pacific Blue-A` < mix2$mu[1]
  gate2[gate0] <- ! gate1[gate0]
  
  col <- ifelse(gate1,"lightgreen","lightgrey")
  col[gate2] <- "pink"

  plot(x$`Pacific Blue-A`,x$`Pacific Blue-W`,pch=20,col=col,cex=.3)
  locator()
  return(list(g1=gate1,g2=gate2))
}) -> gates

for(i in 1:length(e0)) {
  x <- e0[[i]]
  
  group <- ifelse(gates[[i]]$g1,"G1","NA")
  group[gates[[i]]$g2] <- "G2"
  x$group <- group
  ggplot(x) + geom_point(aes(`Pacific Blue-A`,`Pacific Blue-W`,col=group))
  ggplot(subset(x,group != "NA")) + geom_density(aes(`Pacific Blue-A`,col=group))
  ggplot(subset(x,group != "NA")) + geom_density(aes(`FSC-A`,col=group))
  
}

