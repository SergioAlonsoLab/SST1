# Prepare data from xlsx to tsv tables 
# This script is not intended to be shared in the final version
# Only tsv tables will be shared

library(gdata)
library(magrittr)

source("src/auxiliary.R")

subValue <- function(find,sub,where) {
  where[where==find] <- sub
  return(where)
}

try(setwd("~/Documents/WORK/SST1/"))
try(setwd("~/Documents/SST1/"))

# Prepare SST1 MSQPCR data ----

sst1 <- read.xls(("data/SST1_MARIA_data.xlsx"),1,stringsAsFactors=F)
sst1$TYPE <- 
  gsub(" ","",sst1$TYPE) %>%
  gsub("LS-174T","CL",.) %>% factor()

levels(sst1$TYPE) <- c("CL","Normal","Tumor")

sst1$Case.number <- 
  gsub("Cell Line","LS174T",sst1$Case.number) %>% 
  gsub("CLONE","LS174T Clone ",.)

sst1$Normalized.SYBR <- as.numeric(sst1$Normalized.SYBR) %>% log2

sst1 <- sst1[,c("Case.number","TYPE","Normalized.SYBR")]
colnames(sst1) <- c("Case.number","Type","SST1.RDL")
write.table(sst1,file="data/SST1.MSQCR.tsv",sep="\t",quote=F,row.names = F)

rm(sst1)

# Prepare LINE1 MSQPCR data ----

line1 <- read.xls("data/LINE1_summary_Bea.xlsx",2,stringsAsFactors=F)
line1 <- line1[,c("Case.number","TYPE","average.RDR")]
colnames(line1) <- c("Case.number","Type","LINE1.RDL")

line1 <- subset(line1,! is.na(LINE1.RDL))
line1$LINE1.RDL <- log2(line1$LINE1.RDL)

line1$Type <- factor(line1$Type)
levels(line1$Type) <- c("Met.Liver","Normal","Normal.Liver","Tumor")

write.table(line1,file="data/LINE1.MSQCR.tsv",sep="\t",quote=F,row.names = F)

rm(line1)

# Prepare SST1 methylation data of LS174T clones ----

methylation <- read.xls("data/BisulfiteClonesData.xlsx",skip=1,stringsAsFactors=F)
names(methylation)[1:2] <- c("Sample","Sequence")
coors <- methylation[1,paste("CpG",1:29,sep=".")]
methylation <- methylation[-1,]
methylation <- subset(methylation, !Sequence %in% c(".",".."))
methylation$Sample <- factor(methylation$Sample,c("LS174T cell line",paste("LS174T Clone",1:14))) # reorder the clones
levels(methylation$Sample) <- gsub("LS174T ","",levels(methylation$Sample)) %>% gsub("cell line","LS174T",.)

write.table(methylation,file="data/BisulfiteClonesData.tsv",sep="\t",row.names=F,quote=F)
write.table(coors,file="data/coors.tsv",sep="\t",row.names=F,quote=F)
# Prepare Nuclei size data ----

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


# Prepare data from old cases. File used in NAR 2018 ----

allCases <- read.xls("data/SST1-all data.xls",3)
allCases$Classification <- factor(allCases$Hypo.divided.Demethylation.and.Severe.demethylation)
levels(allCases$Classification)[3] <- "NC"
allCases$Classification <- factor(allCases$Classification,c("NC","D","SD"))

allCases$Classification2 <- allCases$Classification
allCases$Classification2[allCases$Difference >= 0.1] <- "SD"
allCases$Classification2[allCases$Difference < 0.1] <- "D"
allCases$Classification2[allCases$Difference < 0.05] <- "NC"

table(allCases$Classification,allCases$Classification2)

names(allCases)[1] <- "Case.number"

write.table(allCases,file="data/allCases.tsv",sep="\t",quote=F,row.names = F)

rm(allCases)

# Prepare IMPPC database


imppc <- read.xls("data/IMPPC_DB.xls")
write.table(imppc,file="data/IMPPC.tsv",sep="\t",quote=F,row.names = F)
rm(imppc)
