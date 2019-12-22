# set working directory. if you don't use that directory for your 
# data files, you will need to use full path for them 
# ex : "C:/windows/sucks/data.txt" instead of just "data.txt"
setwd("E:/Dissertation/library_dcm2")

# read the script file that calibrates the DC4
source("utils/probitUtils.R")
source("utils/create_X.R")
source("utils/models.R")
source("utils/dc_utils.R")
source("models/logit.R")
source("models/dc4.R")

# read your data. if your data is separated by spaces just use read.table
# it it's csv use read.csv
# if your data is not in your working directory use full path to your file
D = read.table("data_2009.txt", header = TRUE)

# common is a list whose elements are vector of the same length than the
# number of alternatives. The elements of each vector corresponds to
# variables that must have the same coefficients for each alternative
common = list(c("zero","logsum1","logsum2","logsum3","logsum4"))

# specific is a list that contains as many vectors as there are alternatives
# each vector contains variables that have their own coefficient for the
# corresponding alternative
specific = list(c(), c("const", "HHFAMINC", "DRVRCNT", "URSIZE",  
"HHR_SEX", "HTRESDN_1000"), c("const", "HHFAMINC", "DRVRCNT", 
"URSIZE",  "HHR_SEX", "HTRESDN_1000"), c("const", "HHFAMINC", 
"DRVRCNT", "URSIZE",  "HHR_SEX", "HTRESDN_1000"), c("const", 
"HHFAMINC", "DRVRCNT", "URSIZE", "HHR_SEX", "HTRESDN_1000"))

# reg is the vector that contains the variables of your regression
reg = c("const", "HHFAMINC", "HOMEOWN", "HHR_SEX", 
"HTRESDN_1000", "MEAN_COST")

# YPrString and YRegString contain the name of the Discrete and the Continuous variable
YPrString = "HHVEHCNT"
YRegString = "MILES_10k"


specPmv = list(
	probit = list(common = common, specific = specific),
	reg = reg,
	YPr = YPrString,
	YReg = YRegString,
	method = "pmvnorm",
	delta = 1e-5,
	reltol = 1e-5,
	SD = "bootstrap",
	nboot = 100,
	verbose = TRUE,
	start = "start_dc4_2009_noPTe-8.txt"
	#output = "~/Dropbox/dc/results_dc4/PmvPlots_dc4_2009_noPTe-5.pdf",
	#plotrange = 0.8,
	# number of points we use for the plots
	# needs to be relatively big
	#numpoints = 200
	)

modelFns = dc4

Sys.time()
result_dc4_2009_nb100 = model(modelFns, specPmv, D)
result_dc4_2009_repara100nb = dc4$reparam(result_dc4_2009_nb100)
Sys.time()

#write.table(result_dc4_2009, "~/Dropbox/dc/results_dc4/result_dc_2009_noPTe-5_bootstrap.txt",
#quote = F, row.names = F)
	
