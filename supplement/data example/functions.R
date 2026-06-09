


############################################################################################################
#---LRT with corncub
############################################################################################################


fitCorncobLRT <- function(dfIn)
{

	suppressWarnings(corncobNull <- tryCatch(
	{
		corncob <- bbdml(formula = cbind(succes, fail) ~ 1, phi.formula = ~1, data = dfIn)
	},
		error=function(cond) {
		return(NULL)
	},
		finally={}
	))


	suppressWarnings(corncob <- tryCatch(
	{
		corncob <- bbdml(formula = cbind(succes, fail) ~ x, phi.formula = ~1, data = dfIn)
	},
		error=function(cond) {
		return(NULL)
	},
		finally={}
	))

	if(!is.null(corncob) & !is.null(corncobNull))
	{
		pval <- lrtest(corncob, corncobNull)

	} else {

		pval <- NA
	}

	return(list(pval = pval, corncobNull = corncobNull, corncob = corncob))

}




########################################################################################################################
#----This function counts the number of ratios below alpha per taxa
#----To exclude certain taxa, a a subset of resLong can also be input
########################################################################################################################


rejectedPerTaxon <- function(resLong, pp = NULL, alpha = NULL)
{
	#resLong <- numDenom
	#resLong <- numDenomIn

	sumHypothesis <- numeric(pp) + NA
	for(i in 1:pp)
	{
		subset <- resLong$nums == i | resLong$denoms == i
		if(sum(subset) > 0) sumHypothesis[i] <- sum(resLong$adj[subset] < alpha, na.rm = TRUE)
	}

	return(list(sumHypothesis = sumHypothesis))

}


########################################################################################################################
#----This function counts the number of ratios below alpha per taxa
########################################################################################################################


seqSumRejected <- function(resLong, alpha = NULL)
{
	#resLong = numDenomIn

	p <- max(resLong$denoms)
	tied = sumRejSeq <- numeric(p) + NA
	nTies <- 0
	include <- !logical(dim(resLong)[1])

	for(i in 1:(p - 1))
	{
		#i <- 1
		#i <- 2
		#i <- 3

		if(sum(resLong[include, ]$adj < alpha) > 0)
		{

			#---counts the number of rejections per taxon
			outReject <- rejectedPerTaxon(resLong[include, ], pp = p, alpha = alpha)
			if(i == 1) sumRej <- outReject$sumHypothesis

			#---Next we want to exclude the taxon with the most rejections
			#---if there are multiple maxima, pick both
 			whichMax <- which(outReject$sumHypothesis == max(outReject$sumHypothesis, na.rm = TRUE))

			#---Save the number of rejections
			sumRejSeq[whichMax] <- outReject$sumHypothesis[whichMax]

			#--the taxon that has the maximum is excluded in the next iteration
			include[resLong$nums %in% whichMax | resLong$denoms %in% whichMax] <- FALSE

			#--diagnostic, to see how ofthen ties in maximum occur
			if(length(whichMax) > 1) 
			{
				nTies <- nTies + 1
				tied[whichMax] <- nTies
			}

		}

	}

	sumRejSeq[is.na(sumRejSeq)] <- 0

	return(list(sumRej = sumRej, sumRejSeq = sumRejSeq, tied = tied))

}





########################################################################################################################
#--this function sorts the ratio betabinomial result
########################################################################################################################


resSorted <- function(numDenomIn, adj = NULL, alpha = NULL)
{
	#numDenomIn <- numDenom
	#adj <- ratioRes[, "fdr.1vs4"]
	#adj <- adjRatio
	#alpha <- 0.05


	#--------------------------------------------------------------------------------
	#----get number of rejections, given adj
	#--------------------------------------------------------------------------------


	ratioResult <- numDenomIn
	ratioResult$adj <- adj
	outSumRej <- seqSumRejected(ratioResult, alpha = alpha)


	#--------------------------------------------------------------------------------
	#----rank taxa
	#--------------------------------------------------------------------------------


	rr1 <- length(outSumRej$sumRej) + 1 - rank(outSumRej$sumRej, ties = "average")
	rr2 <- length(outSumRej$sumRej) + 1 - rank(outSumRej$sumRej, ties = "last")

	ratioResult$rankNums <- rr2[ratioResult$nums]
	ratioResult$rankDenoms <- rr2[ratioResult$denoms]

	rrSeq <- length(outSumRej$sumRejSeq) + 1 - rank(outSumRej$sumRejSeq, ties = "average")
	rrSeq[outSumRej$sumRejSeq == 0] <- NA


	#--------------------------------------------------------------------------------
	#----res per taxon
	#--------------------------------------------------------------------------------


	taxaResult <- data.frame(
			numberTaxon = 1:length(outSumRej$sumRej), 
			sumReject = outSumRej$sumRej, 
			sumRejRankAve = rr1, 
			sumRejRank = rr2, 
			sumRejSeq = outSumRej$sumRejSeq, 
			sumRejSeqRank = rrSeq, 
			seqTied = outSumRej$tied)


	#----------------------------------------------------------------
	#----Filter result for taxa without rejections
	#----------------------------------------------------------------

	zeroReject <- taxaResult$sumReject == 0
	filt <- which(zeroReject)
	index <- !(ratioResult$nums %in% filt | ratioResult$denoms %in% filt)

	ratioResultFilt <- ratioResult[index, ]
	taxaResultFilt <- taxaResult[!zeroReject, ]


	#----------------------------------------------------------------


	#--add also a ranked version of the taxa result
	taxaResultFiltRanked <- taxaResultFilt[order(taxaResultFilt$sumRejRank), ]


	attr(ratioResult, "alpha") <- alpha
	attr(ratioResult, "alpha") <- alpha
	attr(ratioResultFilt, "alpha") <- alpha

	return(list(ratioResult = ratioResult, ratioResultFilt = ratioResultFilt, taxaResult = taxaResult, taxaResultFilt = taxaResultFilt, taxaResultFiltRanked = taxaResultFiltRanked))

}


########################################################################################################################
#---This is the function that makes the taxa versus taxa plot
########################################################################################################################



figPairs <- function(sortedResIn, type = 1, taxaPlot = NULL)
{

	alpha <- attr(sortedResIn$ratioResult, "alpha")
	taxTot <- length(sortedResIn$taxaResult$numberTaxon)

	if(type == 1)
	{

		taxaPlot <- sortedResIn$taxaResult$numberTaxon
		ratioResPlot <- sortedResIn$ratioResult

		taxaResultPlot <- sortedResIn$taxaResult
		taxaResultPlot <- taxaResultPlot[order(taxaResultPlot$sumRejRank), ]

	}

	if(type == 2)
	{

		taxaPlot <- sortedResIn$taxaResultFilt$numberTaxon
		ratioResPlot <- sortedResIn$ratioResultFilt

		taxaResultPlot <- sortedResIn$taxaResultFilt
		taxaResultPlot <- taxaResultPlot[order(taxaResultPlot$sumRejRank), ]

	}

	if(type == 3) 
	{ 

		#--take subset ratios
		index <- (sortedResIn$ratioResult$nums %in% taxaPlot & sortedResIn$ratioResult$denoms %in% taxaPlot)
		ratioResPlot <- sortedResIn$ratioResult[index, ]

		#--take subset taxa
		index2 <- sortedResIn$taxaResult$numberTaxon %in% taxaPlot
		taxaResultPlot <- sortedResIn$taxaResult[index2 , ]

		#--new ranking for ordering in the plot
		taxaResultPlot$sumRejRankNew <- as.numeric(factor(taxaResultPlot$sumRejRank))

		#--add new ranking to ratios
		tmp <- numeric(taxTot) + NA
		tmp[taxaResultPlot$numberTaxon] <- taxaResultPlot$sumRejRankNew
		ratioResPlot$rankNums <- tmp[ratioResPlot$nums]
		ratioResPlot$rankDenoms <- tmp[ratioResPlot$denoms]

		taxaResultPlot <- taxaResultPlot[order(taxaResultPlot$sumRejRank), ]

	}


	cc <- cut(ratioResPlot$adj, c(0, alpha, 1))


      dfPlot <- data.frame(
      	numeratorRank = factor(c(ratioResPlot[, "rankNums"], ratioResPlot[, "rankDenoms"])), 
      	nums = factor(c(ratioResPlot[, "nums"], ratioResPlot[, "nums"])), 
      	denominatorRank = factor(c(ratioResPlot[, "rankDenoms"], ratioResPlot[, "rankNums"])), 
      	denoms = factor(c(ratioResPlot[, "denoms"], ratioResPlot[, "denoms"])), 
      	FDR = factor(c(as.character(cc), as.character(cc))))


	nTaxPlot <- length(taxaPlot)


      p <- ggplot(data = dfPlot, aes(x = numeratorRank, y = denominatorRank))
      p <- p + geom_tile(col = "grey60", aes(fill = FDR))
      p <- p + scale_x_discrete(expand = c(0, 0))
      p <- p + scale_y_discrete(expand = c(0, 0))
      p <- p + ggtitle("")
      p <- p + theme_classic()


	tabcc <- table(cc)
	if(tabcc[1] == 0) p <- p + scale_fill_manual(values=c("grey60"))
	if(tabcc[2] == 0) p <- p + scale_fill_manual(values=c("red"))
	if(all(tabcc > 0)) p <- p + scale_fill_manual(values=c("red", "grey60"))


	dfPoint <- data.frame(numeratorRank = 1:nTaxPlot, denominatorRank =  1:nTaxPlot, numberTaxon = taxaResultPlot$numberTaxon)
	p <- p + geom_text(data = dfPoint, label = taxaResultPlot$numberTaxon, cex = 2.5, fontface = "bold")
	p <- p + theme(legend.position = "none")
	p <- p + xlab("numerator taxon")
	p <- p + xlab("numerator taxon")
	p <- p + ylab("denominator taxon")
	p <- p + theme(axis.title=element_text(size = 15, face = "plain"))
	p <- p + theme(axis.text = element_blank())
	p <- p + theme(axis.ticks = element_blank())


	#---margins
	p <- p + theme(plot.margin = unit(c(-15, 3, 2, 2), "pt"))
	#--default all are 5.5
	#top, right*, bot, left


	return(list(p = p, dfPlot = dfPlot, dfPoint = dfPoint))

}

