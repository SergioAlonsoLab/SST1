# Normalize beta values using BMIQ method from wateRmelon
# This script takes some time to complete

library(wateRmelon)
library(IlluminaHumanMethylation450kmanifest)

original_file <- "/imppc/labs/mplab/share/Illumina450/RnBeads_report2/export_data/csv/betas_1.csv"

betas1 <- as.data.frame(data.table::fread(original_file))
probeType <- ifelse(betas1$ID %in% IlluminaHumanMethylation450kmanifest@data$TypeI$Name,1,2)

for(i in 6:ncol(betas1)) {
  betas1[,i] <- BMIQ(betas1[,i],probeType)$nbeta
}

saveRDS(betas1,"data/beta_values_BMIQ.rds")


