

library(tidyverse)
library(ggplot2)
library(corncob)


#gitdir <- "C:/1Files/2Git/projects2/3microbiome/danone/1final/supplement"
#setwd(gitdir)



############################################################################################################


infant_data <- readRDS(file = "simulation code and data/data/infant_data_filtered.RDS")
infantWide <- infant_data %>% spread(taxon, count)
infantWideDF <- data.frame(infantWide)

dim(infantWideDF)
#[1] 272  53

cntsMat <- as.matrix(infantWideDF[, c(-1:-2)])
#dim(cntsMat)
#[1] 272  51


############################################################################################################


metaFile <- read.table(file = "data example/Sample ids data.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)


#--check
#sum(infantWideDF[, 1] != metaFile[, 1])
#0


metData <- data.frame(infantWideDF[, c(1:2)])
metData$study <- factor(metData$study)
metData$ArmSimplified <- metaFile$ArmSimplified


rm(infantWideDF, infantWide, infant_data, metaFile)
#ls()
#[1] "cntsMat" "gitdir"  "metData"
