# Analyze beta values

library(ggplot2)
library(magrittr)
library(ggrepel)
library(ggsci)

betas0 <- readRDS("data/beta_values_BMIQ.rds")
rownames(betas0) <- betas0$ID

samples <- colnames(betas0)[6:ncol(betas0)]
sampleType <- gsub("[0-9bR_]","",samples) %>% factor(.,c("N","T","NL","ML"))


# The data has been previosuly normalized using BMIQ in BMIQ_normalization.R
# No need to check it again

# XY and chromosomes have been removed
table(betas0$Chromosome)

# PCA and UMAP analyses

library(umap)

pca0 <- prcomp(t(betas0[,samples])) # unscaled because beta values MUST range between 0 and 1
umap0 <- umap(t(betas0[,samples]))

# PCA shows that FFPE samples are different from the rest, and 458T is an outlayer among tumors

gg0 <- ggplot(data.frame(sampleType=sampleType,UMAP=umap0$layout,pca0$x))

gg1 <- gg0 +  geom_point(aes(fill=sampleType),alpha=.8,pch=21,size=4) +
  geom_label_repel(label=samples,alpha=.8) + scale_fill_aaas() + theme(axis.text = element_text(size=12)) 

gg1 + aes(x=PC1,y=PC2) 
gg1 + aes(x=PC2,y=PC3) 
gg1 + aes(x=PC3,y=PC4)

gg1 + aes(x=UMAP.1,y=UMAP.2)

gg1 + aes(x=PC1,y=UMAP.1)
gg1 + aes(x=PC2,y=UMAP.1)
gg1 + aes(x=PC1,y=UMAP.2)
gg1 + aes(x=PC2,y=UMAP.2)

