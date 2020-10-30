# read and process nuclei size

library(gdata)
library(magrittr)

source("src/auxiliary.R")

nucSizeLS <- read.xls("data/NucleoSizeImageJ.xlsx",3)
nucSizeOV <- read.xls("data/OV90nucleusSize_BG.xlsx")[,c("Sample","areas")]
colnames(nucSizeOV) <- c("sample","area")

ggplot(nucSizeLS,aes(sample,area)) + geom_boxplot() + geom_hline(yintercept = 1000) + vxt
ggplot(nucSizeOV,aes(sample,area)) + geom_boxplot() + geom_hline(yintercept = 1000) + vxt

# scale areas

nucSizeLS$area <- nucSizeLS$area * (50/155) ^ 2 # 155 pixels correspond to 50 microns
nucSizeOV$area <- nucSizeOV$area * (50/155) ^ 2 # 155 pixels correspond to 50 microns

# edit clone names

nucSizeLS$sample <- factor(nucSizeLS$sample,c("LS174_23",paste("CLONE",1:14,sep="_")))
levels(nucSizeLS$sample)[1] <- "LS174T"
nucSizeLS <- na.omit(nucSizeLS)

nucSizeOV$sample %>% gsub("OV90-C","CLONE_",.) %>% gsub("_CL","",.) %>% factor(.,c("OV90",paste("CLONE",1:17,sep="_"))) -> nucSizeOV$sample

ggplot(nucSizeLS,aes(sample,area)) + geom_boxplot() + vxt
ggplot(nucSizeOV,aes(sample,area)) + geom_boxplot() + vxt


write.table(nucSizeLS,file = "data/NucleiSizeLS.tsv",sep="\t",quote=F,row.names = F)
write.table(nucSizeOV,file = "data/NucleiSizeOV.tsv",sep="\t",quote=F,row.names = F)
