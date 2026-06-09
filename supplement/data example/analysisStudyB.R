

rm(list = ls())


#gitdir <- "C:/1Files/2Git/projects2/3microbiome/danone/1final/supplement"
gitdir <- "M:/2projects/1completed/2017-2022/1microbiome/danone/1final/supplement/data example"
setwd(gitdir)


####################################################################################################


source("data example/functions.R")

ls()
#[1] "figPairs"         "figPairs2"        "fitCorncobLRT"
#[5] "gitdir"           "numDenomFig"      "rejectedPerTaxon" "seqSumRejected"  


####################################################################################################


#---load data and get subset for study B
source("data example/loadData.R")
ls()


indexSample <- metData$study == "B"
cnts <- cntsMat[indexSample, ]
zeroes <- colSums(cnts == 0) / dim(cnts)[1]
cnts <- cntsMat[indexSample, ]
treat <- factor(metData$ArmSimplified[indexSample])
nSamp <- length(treat)


table(treat)
#Test_1 Test_2 Test_3 Test_4 
#    34     37     42     36 



####################################################################################################


#---Create data.frame with ratio combinations
pp <- dim(cnts)[2]
nums <- rep(1:(pp-1), times = (pp-1):1)
denoms <- c()
for(j in 2:pp) denoms <- c(denoms, j:pp)
numDenom <- data.frame(nums, denoms)
pRatio <- dim(numDenom)[1]


#--model matrix for treatment comparisons
#--used for selecting subsets of data
comps <- cbind(1:3, 4)
nComp <- dim(comps)[1]
mm <- model.matrix(~treat - 1)

#---these are the treatment comparisons
#     [,1] [,2]
#[1,]    1    4
#[2,]    2    4
#[3,]    3    4


####################################################################################################


pvalsBB = estBB <- matrix(nrow = pRatio, ncol = nComp)
for(ij in 1:pRatio)
{
	if(ij %in% 10 == 0) print(ij)
	#ij <- 1140

	#--select indices for numerator / denominator
	i <- numDenom[ij, 1]
	j <- numDenom[ij, 2]

	#---loop through the pairwise comparissons
	for(k in 1:nComp)
	{
		#k <- 1

		#--select relevant columns of the design matrix
		xPW <- mm[, comps[k, ]]
		index <- rowSums(xPW) > 0

		#--put stuff in DF and fit
		dfIn <- data.frame(succes = cnts[index, i], fail = cnts[index, j], x = factor(treat[index]))
		fitPC <- fitCorncobLRT(dfIn)

		if(!is.na(fitPC$pval))
		{

			pvalsBB[ij, k] <- fitPC$pval
			estBB[ij, k] <- summary(fitPC$corncob)$coef[2, 1]

		}
	}
}


#---pvals and FDR
pvalsBB[is.na(pvalsBB)] <- 1
fdrBB <- apply(pvalsBB, 2,  function(x) p.adjust(x, method = "BH"))


#fracZero <- colSums(cnts == 0) / dim(cnts)[1]
#numDenom$numsZero <- fracZero[numDenom$nums]
#numDenom$denomsZero <- fracZero[numDenom$denoms]



####################################################################################################
#---Analyse met library size as denominator
####################################################################################################



pvalsRel = tvalsRel = estRel <- matrix(nrow = pp, ncol = nComp)
for(i in 1:pp)
{
	print(i)
	#i <- 1

	for(k in 1:nComp)
	{
		#k <- 1

		xPW <- mm[, comps[k, ]]
		index <- rowSums(xPW) > 0

		#--prepare DF 
		dfIn <- data.frame(succes = cnts[index, i], fail = rowSums(cnts[index, ]) - cnts[index, i], x = factor(treat[index]))

		#--fit corncub
		fitPC <- fitCorncobLRT(dfIn)

		pvalsRel[i, k] <- fitPC$pval
		estRel[i, k] <- summary(fitPC$corncob)$coef[2, 1]
		tvalsRel[i, k] <- summary(fitPC$corncob)$coef[2, 3]

	}

}


fdrRel <- apply(pvalsRel, 2,  function(x) p.adjust(x, method = "BH"))



####################################################################################################
#-----Sort result of ratios
####################################################################################################


txt1 <- paste("pval", paste(comps[1, ], collapse = "vs"), sep = ".")
txt2 <- paste("pval", paste(comps[2, ], collapse = "vs"), sep = ".")
txt3 <- paste("pval", paste(comps[3, ], collapse = "vs"), sep = ".")
colnames(pvalsBB) <- c(txt1, txt2, txt3)
colnames(pvalsRel) <- c(txt1, txt2, txt3)


txt1 <- paste("fdr", paste(comps[1, ], collapse = "vs"), sep = ".")
txt2 <- paste("fdr", paste(comps[2, ], collapse = "vs"), sep = ".")
txt3 <- paste("fdr", paste(comps[3, ], collapse = "vs"), sep = ".")
colnames(fdrBB) <- c(txt1, txt2, txt3)
colnames(fdrRel) <- c(txt1, txt2, txt3)


txt1 <- paste("est", paste(comps[1, ], collapse = "vs"), sep = ".")
txt2 <- paste("est", paste(comps[2, ], collapse = "vs"), sep = ".")
txt3 <- paste("est", paste(comps[3, ], collapse = "vs"), sep = ".")
colnames(estBB) <- c(txt1, txt2, txt3)
colnames(estRel) <- c(txt1, txt2, txt3)

#head(estBB)


ratioRes <- data.frame(pvalsBB, fdrBB, estBB)
relRes <- data.frame(pvalsRel, fdrRel, estRel)

#head(ratioRes)
#head(relRes)


if(FALSE) save(numDenom, ratioRes, relRes, file = "data example/resStudyB.RData")



