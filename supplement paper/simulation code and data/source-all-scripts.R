# Open R in command line and source this script to run simulation study. 
# Note that parallel loops are used, the functions in those do not work well
# with RStudio. 6 cores are used, if not available, differences in random-stream
# are to be expected.
source("01-data-filtering-and-dispersion-estimation.R")
source("02-parametric-simulations-bifido.R")
source("03-parametric-simulations-bifido-extended.R")
source("04-parametric-simulations-anaerostipes.R")
source("05-sampling-simulations-bifido.R")
source("06-sampling-simulations-anaero.R")
source("07-process-simulation-output.R")
source("08-create-result-visualizations.R")
source("09-create-figure-1-raw-vs-fitted.R")