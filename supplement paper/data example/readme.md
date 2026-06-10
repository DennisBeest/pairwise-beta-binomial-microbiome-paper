
# code data example

The file analysisStudyB.R does the analysis. 

This files sources loadData.R and functions.R for the needed data and functions.  
It fits both pairwise ratio beta-binomial model, and the MWW beta-binomial model (relative abundances).

The data is loaded from the "simulation code and data" directory. The treatments are loaded from Sample ids data.csv. 
The result of the analysis is stored in the resStudyB.RData object.

The file figRatio.R provides an overview of the result and makes the visualisation of the paper. 
The figRatio.R file also sources the functions.R.
The functions sourced by figRatio.R create an overview of the result, and make a basic figure.
These are then modified and combined figRatio.R file.







