
rm(list = ls())
library(vegan)
library(ggplot2)
library(corncob)
library(gridExtra)
library(grid)


gitdir <- "C:/1Files/2Git/projects2/3microbiome/danone/1final/supplement"
setwd(gitdir)


####################################################################################################


source("data example/functions.R")
ls()


####################################################################################################


#---load data and get subset for study B
load("data example/resStudyB.RData")
ls()
#[1] "numDenom" "ratioRes" "relRes" 
#out2


#sum(ratioRes$fdr.1vs4 < 0.05, na.rm = TRUE)
#6
#sum(ratioRes$fdr.2vs4 < 0.05, na.rm = TRUE)
#107
#sum(ratioRes$fdr.3vs4 < 0.05, na.rm = TRUE)
#100



####################################################################################################


outSorted1 <- resSorted(numDenomIn = numDenom, adj = ratioRes[, "fdr.1vs4"], alpha = 0.05)
outSorted2 <- resSorted(numDenomIn = numDenom, adj = ratioRes[, "fdr.2vs4"], alpha = 0.05)
outSorted3 <- resSorted(numDenomIn = numDenom, adj = ratioRes[, "fdr.3vs4"], alpha = 0.05)



summary(outSorted2)
#--ratioResult
#--ratioResultFilt
#--taxaResult
#--taxaResultFilt
#--taxaResultFiltRanked


#head(outSorted2$ratioResult)
#--> Result per ratio (all)
#--> This object rankNums and rankDenoms are the taxa ranked by the number of hypothesis with a FDR below alpha
#--> nums: taxon used as numerator (a in a / (a + b))
#--> denoms: taxon used as denominator
#--> adjusted p-value as used as input. These are here equal:
#--------data.frame(outSorted2$ratioResult$adj, ratioRes[, "fdr.2vs4"])
#--> rankNums: this is the rank of the numerator taxon when ranked the number of rejected hypothesis (sumReject, see taxaResult)
#--> rankDenoms: this is the rank of the demerator taxon when ranked the number of rejected hypothesis (sumReject, see taxaResult)


#head(outSorted2$ratioResultFilt)
#--> Result per ratio filtered for taxa without rejections


#head(outSorted2$taxaResult)
#--> results sorted to taxa.
#--------sumReject: total number of rejections per taxa. Each ratios is counted twice here.
#--------sumRejRankAve: Average rank of taxa, ranked on sumReject
#--------sumRejRank: Rank of taxa, ranked on sumReject, ties = "last"
#--------sumRejSeq: number of rejections per taxa, "unwrapped". Each ratios is counted once here. 
#Each ratio is assiged to the taxa with highest number of rejections
#--------sumRejSeqRank: rank based on sumRejSeq
#--------seqTied: diagnostic to see if taxa are tied for sumRejSeq
#--note that sumReject is double sumRejSeq, and double the total number of rejections
#sum(outSorted2$taxaResult$sumReject)
#sum(outSorted2$taxaResult$sumRejSeq)
#sum(outSorted2$ratioResult$adj < 0.05)


#head(outSorted2$taxaResultFilt)
#--> Result per taxa filtered for taxa without rejections


#head(outSorted2$taxaResultFiltRanked)
#--> As taxaResultFilt but ranked


####################################################################################################


#----select all taxa that were involved in a significant ratio in any of the comparisons
taxaPlot <- sort(unique(c(outSorted1$taxaResultFilt$numberTaxon, outSorted2$taxaResultFilt$numberTaxon, outSorted3$taxaResultFilt$numberTaxon)))
length(taxaPlot)


#---figPairs has three options for type
#------type 1 includes all taxa
#------type 2 includes only taxa with significant ratio
#------type 3 includes the taxa specified with taxaPlot
#------In this figure the aim is to include all taxa signifcant in any of the 3 treatment comparisons
#------The numbers in the diagonal come from numberTaxon, which correspond to the colnumber of the count matrix used as input


pout1 <- figPairs(sortedResIn = outSorted1, type = 3, taxaPlot = taxaPlot)
p1 <- pout1$p
pout2 <- figPairs(sortedResIn = outSorted2, type = 3, taxaPlot = taxaPlot)
p2 <- pout2$p
pout3 <- figPairs(sortedResIn = outSorted3, type = 3, taxaPlot = taxaPlot)
p3 <- pout3$p


#--add titles
p1 <- p1 + ggtitle("Treatment 1 versus Treatment 4")
p1 <- p1 + theme(plot.title = element_text(size = 16, face = "plain", hjust = 0.5, vjust = -8))
p2 <- p2 + ggtitle("Treatment 2 versus Treatment 4")
p2 <- p2 + theme(plot.title = element_text(size = 16, face = "plain", hjust = 0.5, vjust = -8))
p3 <- p3 + ggtitle("Treatment 3 versus Treatment 4")
p3 <- p3 + theme(plot.title = element_text(size = 16, face = "plain", hjust = 0.5, vjust = -8))


#---add letters
fz <- 14
yloc <- 0.98
grobA <- grobTree(textGrob("(A)", x=0.08,  y=yloc, gp=gpar(fontsize = fz, face = "bold")))
grobB <- grobTree(textGrob("(B)", x=0.08,  y=yloc, gp=gpar(fontsize = fz, face = "bold")))
grobC <- grobTree(textGrob("(C)", x=0.08,  y=yloc, gp=gpar(fontsize = fz, face = "bold")))
p1 <- p1 + annotation_custom(grobA)
p2 <- p2 + annotation_custom(grobB)
p3 <- p3 + annotation_custom(grobC)


#--add legend
p1 <- p1 + theme(legend.position = c(0.15, 0.85), 
	legend.text = element_text(colour = "black", size = 11),
	legend.key = element_blank(),
	legend.title = element_text(colour = "black", size = 13))
p1 <- p1 + theme(legend.key.height = unit(0.52, "cm"), legend.key.width = unit(0.93, "cm"))




########################################################################
#---add dots for taxa significant from relative betabin
########################################################################
#---ggplot code
#https://stackoverflow.com/questions/12409960/ggplot2-annotate-outside-of-plot


sigRel1 <- which(relRes$fdr.1vs4 < 0.05)
ww1 <- which(pout1$dfPoint$numberTaxon %in% sigRel1)

sigRel2 <- which(relRes$fdr.2vs4 < 0.05)
ww2 <- which(pout2$dfPoint$numberTaxon %in% sigRel2)

sigRel3 <- which(relRes$fdr.3vs4 < 0.05)
ww3 <- which(pout3$dfPoint$numberTaxon %in% sigRel3)


p1 <- p1 + annotation_custom(
	grob = pointsGrob(pch = 20, gp = gpar(cex = 1.2)),
      ymin = -1, ymax = -1, xmin = ww1, xmax = ww1)
gt1 <- ggplot_gtable(ggplot_build(p1))
gt1$layout$clip[gt1$layout$name == "panel"] <- "off"



for(i in 1:length(ww2))
{
	#i <- 1

	p2 <- p2 + annotation_custom(
		grob = pointsGrob(pch = 20, gp = gpar(cex = 1.2)),
	      ymin = -1, ymax = -1, xmin = ww2[i], xmax = ww2[i])
	gt2 <- ggplot_gtable(ggplot_build(p2))
}
gt2$layout$clip[gt2$layout$name == "panel"] <- "off"



for(i in 1:length(ww3))
{
	#i <- 1

	p3 <- p3 + annotation_custom(
		grob = pointsGrob(pch = 20, gp = gpar(cex = 1.2)),
	      ymin = -1, ymax = -1, xmin = ww3[i], xmax = ww3[i])
	gt3 <- ggplot_gtable(ggplot_build(p3))
}
gt3$layout$clip[gt3$layout$name == "panel"] <- "off"



########################################################################
#---add dots for taxa from ratio analysis
########################################################################


sigRatio1 <- outSorted1$taxaResultFiltRanked$numberTaxon[outSorted1$taxaResultFiltRanked$sumRejSeq > 1]
sigRatio2 <- outSorted2$taxaResultFiltRanked$numberTaxon[outSorted2$taxaResultFiltRanked$sumRejSeq > 1]
sigRatio3 <- outSorted3$taxaResultFiltRanked$numberTaxon[outSorted3$taxaResultFiltRanked$sumRejSeq > 1]


ww1 <- which(pout1$dfPoint$numberTaxon %in% sigRatio1)
ww2 <- which(pout2$dfPoint$numberTaxon %in% sigRatio2)
ww3 <- which(pout3$dfPoint$numberTaxon %in% sigRatio3)


p1 <- p1 + annotation_custom(
	grob = pointsGrob(pch = 17, gp = gpar(cex = 0.8)),
      ymin = -0.2, ymax = -0.2, xmin = ww1, xmax = ww1)
gt1 <- ggplot_gtable(ggplot_build(p1))
gt1$layout$clip[gt1$layout$name == "panel"] <- "off"
#grid.draw(gt1)


for(i in 1:length(ww2))
{
	#i <- 1

	p2 <- p2 + annotation_custom(
		grob = pointsGrob(pch = 17, gp = gpar(cex = 0.8)),
	      ymin = -0.2, ymax = -0.2, xmin = ww2[i], xmax = ww2[i])
	gt2 <- ggplot_gtable(ggplot_build(p2))
}
gt2$layout$clip[gt2$layout$name == "panel"] <- "off"
#grid.draw(gt2)


for(i in 1:length(ww3))
{
	#i <- 1

	p3 <- p3 + annotation_custom(
		grob = pointsGrob(pch = 17, gp = gpar(cex = 0.8)),
	      ymin = -0.2, ymax = -0.2, xmin = ww3[i], xmax = ww3[i])
	gt3 <- ggplot_gtable(ggplot_build(p3))
}
gt3$layout$clip[gt3$layout$name == "panel"] <- "off"
#grid.draw(gt3)


########################################################################



dev.new(width = 15, height = 7)
grid.arrange(
	gt1, gt2, gt3,
nrow = 1, ncol = 3)



#--prevent accidental overwriting
doWrite <- FALSE
if(doWrite)
{
	nameFile <- "data example/dataExample.pdf"
	pdf(nameFile,width=15, height=7)
		grid.arrange(
			gt1, gt2, gt3,
		nrow=1, ncol=3)
	dev.off() 
}



