

rm(list = ls())

library(ggplot2)
library(tidyverse)
library(corncob)


gitdir <- "C:/1Files/2Git/projects2/3microbiome/danone/1final/supplement"
setwd(gitdir)

source("data example/functions.R")
source("data example/loadData.R")
ls()


####################################################################################################


indexSample <- metData$study == "B"
cnts <- cntsMat[indexSample, ]
treat <- factor(metData$ArmSimplified[indexSample])


ww <- which(colSums(cnts == 0) == 0)
#ww[2:3]


####################################################################################################


nn <- dim(cnts)[1]
pp <- dim(cnts)[2]
nums <- rep(1:(pp-1), times = (pp-1):1)
denoms <- c()
for(j in 2:pp) denoms <- c(denoms, j:pp)
numDenom <- data.frame(nums, denoms)
pRatio <- dim(numDenom)[1]


####################################################################################################


nSim <- 50
saveRel = saveRatio = saveRatio2 <- list()

set.seed(11)
generateSeeds <- sample(1:1e6, nSim, replace = FALSE)

nSigRatio = nSigRatioMod <- numeric()
estRelSave = estRatioSave <- list()
pvalsRelSave = pvalsRatioSave <- list()

for(mm in 1:nSim)
{

	#mm <- 2
	print(mm)

	set.seed(generateSeeds[mm])
	treatSim <- rep(c(0, 1), times = c(74, 75))[sample(nn)]

	cntsSim <- cnts
	cntsSim[treatSim == 1, ww[2]] <- cntsSim[treatSim == 1, ww[2]] * 8
	cntsSim[treatSim == 1, ww[3]] <- cntsSim[treatSim == 1, ww[3]] * 8


	####################################################################################################


	pvalsRel = estRel <- numeric()
	for(i in 1:pp)
	{
		#i <- 1
		#print(i)

		dfIn <- data.frame(succes = cntsSim[, i], fail = rowSums(cntsSim) - cntsSim[, i], x = treatSim)

		corncob <- bbdml(formula = cbind(succes, fail) ~ x, phi.formula = ~1, data = dfIn)
		corncobNull <- bbdml(formula = cbind(succes, fail) ~ 1, phi.formula = ~1, data = dfIn)

		summ <- summary(corncob)
		estRel[i] <- summ$coeff[2,1]
		pvalsRel[i] <- lrtest(corncob, corncobNull)

	}


	adjRel <- p.adjust(pvalsRel, method = "BH")
	saveRel[[mm]] <- which(adjRel < 0.05)	


	####################################################################################################


	pvalsRatio = estRatio <- numeric()
	for(ij in 1:pRatio)
	{

		#ij <- 1
		#ij <- 342

		#--select indices for numerator / denominator
		i <- numDenom[ij, 1]
		j <- numDenom[ij, 2]


		#--fit main model for scenario
		dfIn <- data.frame(succes = cntsSim[, i], fail = cntsSim[, j], x = treatSim)
		fit <- fitCorncobLRT(dfIn)

		if(!is.na(fit$pval))
		{

			pvalsRatio[ij] <- fit$pval
			estRatio[ij] <- fit$corncob$b.mu[2]

		}
	}

	estRelSave[[mm]] <- estRel
	estRatioSave[[mm]] <- estRatio

	pvalsRelSave[[mm]] <- pvalsRel
	pvalsRatioSave[[mm]] <- pvalsRatio

	pvalsRatio[is.na(pvalsRatio)] <- 1
	adjRatio <- p.adjust(pvalsRatio, method = "BH")
	nSigRatio[mm] <- sum(adjRatio < 0.05)

	indexTmp <- (numDenom$nums == ww[2] | numDenom$denoms == ww[2] | numDenom$nums == ww[3] | numDenom$denoms == ww[3]) & !(numDenom$nums == ww[3] & numDenom$denoms == ww[2] | numDenom$nums == ww[2] & numDenom$denoms == ww[3])
	nSigRatioMod[mm] <- sum(adjRatio[indexTmp] < 0.05)

	if(sum(adjRatio < 0.05) > 1)
	{
		out1 <- resSorted(numDenom, adj = adjRatio, alpha = 0.05)
		saveRatio[[mm]] <- out1$taxaResultFilt$numberTaxon[out1$taxaResultFilt$sumRejSeq > 1]
		saveRatio2[[mm]] <- out1$taxaResultFilt$numberTaxon[out1$taxaResultFilt$sumRejSeq > 0]

	} else {

		saveRatio[[mm]] <- NULL
	}

}


shifted <- ww[2:3]
tabRel <- table(unlist(saveRel))
tabRatio <- table(unlist(saveRatio))


shiftedRatio <- numDenom$nums %in% shifted | numDenom$denoms %in% shifted
shiftedRatio[numDenom$nums %in% shifted & numDenom$denoms %in% shifted] <- FALSE
#table(shiftedRatio)
#table(shiftedRatio)



#--sort estimates
nSim <- length(estRelSave)
estRelShift = estRelNull <- c()
estRatioShift = estRatioNull <- c()
for(m in 1:nSim)
{
	estRelShift <- c(estRelShift, estRelSave[[m]][shifted])
	estRelNull <- c(estRelNull, estRelSave[[m]][-shifted])

	estRatioShift <- c(estRatioShift, estRatioSave[[m]][shiftedRatio])
	estRatioNull <- c(estRatioNull, estRatioSave[[m]][!shiftedRatio])
}




#--total identified
sum(tabRel)
#--incorrect identified
sum(tabRel[!names(tabRel) %in% shifted])
#--correct identified
sum(tabRel[names(tabRel) %in% shifted])


#---Number rejected ratios, and how many are in modified taxa
#--total
sum(nSigRatio)
#---correct
sum(nSigRatioMod)
#---correct rate
1 - sum(nSigRatioMod) / sum(nSigRatio)



save <- FALSE
if(save) save(saveRel, saveRatio, saveRatio2,
		estRelSave, estRatioSave,
		pvalsRelSave, pvalsRatioSave,
		nSigRatio, nSigRatioMod,
		shifted, numDenom,
		tabRel, tabRatio,
		estRelShift, estRelNull, estRatioShift, estRatioNull, shiftedRatio,
		file = "supplement example/supplCompareRelRatio.RData")




################################################################################
################################################################################





#---selected taxa per sim in a list
saveRel
saveRatio
saveRatio2

#---estimates per simulation in a list
estRelSave
estRatioSave

#---p-value per simulation in a list
pvalsRelSave
pvalsRatioSave

#---Number rejected ratios, and how many are in modified taxa
sum(nSigRatioMod)
sum(nSigRatio)
sum(nSigRatioMod) / sum(nSigRatio)

(sum(nSigRatio) - sum(nSigRatioMod)) / sum(nSigRatio)
