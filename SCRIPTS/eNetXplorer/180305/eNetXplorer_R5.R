library(missForest)

rm(list=ls())

source("/Users/candiajm/ACTIVE/H5N1/March_2018/SCRIPTS/eNetXplorer/eNetXplorer_devel3.R")

path_infile = "/Users/candiajm/ACTIVE/H5N1/March_2018/DATA/eNet/"
path_outfile = "/Users/candiajm/ACTIVE/H5N1/March_2018/RESULTS/eNetXplorer/"

#################################################################
# choose parameter set to reproduce results reported in the paper
run_index = 5 # integer value in the 1:5 range
if (run_index==1) {
    model = "FP_GPgset"
    donor_group = "ALL"
    modifier = ""
} else if (run_index==2) {
    model = "FP_GPgset"
    donor_group = "NEG"
    modifier = ""
} else if (run_index==3) {
    model = "FP_GPgset"
    donor_group = "POS"
    modifier = "wo44"
} else if (run_index==4) {
    model = "FB_GBpax_GBpbmc"
    donor_group = "NEG"
    modifier = ""
} else if (run_index==5) {
    model = "SOMAscan"
    donor_group = "ADJ"
    modifier = "wo44"
}
#################################################################

endpoint = "titer_I_MN_d28"
# Subject metadata, endpoints, predictors.
data = read.table(paste0(path_infile,"eNet_InputData_NEW.txt"),header=T,stringsAsFactors=F,sep="\t")

data[data[,"sex"]!="Male Gende","sex"] = 0 # F=0
data[data[,"sex"]=="Male Gende","sex"] = 1 # M=1

donor = data[,"donor_short"]
n_donor = length(donor)

# We select the predictor matrix
if (model=="FP_GPgset") {
    predictor = c("sex","age",paste0("Fp0",1:8),paste0("Gp0",1:9),paste0("Gp",10:14))
}
if (model=="FB_GBpax_GBpbmc") {
    predictor = c("sex","age",paste0("Fb0",1:9),paste0("Fb",10:83),paste0("Gb0",1:9),paste0("Gb",10:23),paste0("GbWB0",1:9),paste0("GbWB",10:19))
}
if (model=="SOMAscan") {
    predictor = c("sex","age",colnames(data)[171:1475])
}
x = data[,predictor]

# We select the endpoint
if (endpoint=="titer_I_MN_d28-d0") {
    y = data[,"titer_I_MN_d28"]-data[,"titer_I_MN_d0"]
}
if (endpoint=="titer_I_MN_d28") {
    y = data[,"titer_I_MN_d28"]
}

# We select the cohort
if (donor_group=="ALL") {
    donor_select = rep(T,n_donor)
} else if ((donor_group=="POS")|(donor_group=="NEG")) {
    donor_select = data[,"blinded_flag"]==donor_group
} else if ((donor_group=="ADJ")|(donor_group=="NONADJ")) {
    donor_select = data[,"unblinded_flag"]==donor_group
} else {
    stop("Donor group not properly specified!\n")
}
if (modifier=="wo44") {
    donor_select[donor=="s44"]=F
}

donor_remove = apply(is.na(x),1,sum)>10 # we exclude subjects with missing values on many (>10) predictor variables
donor_select = donor_select&(!donor_remove)

donor = donor[donor_select]
x = x[donor_select,]
descriptor = colnames(x)
x = matrix(as.numeric(unlist(x)),ncol=ncol(x))
rownames(x) = donor
colnames(x) = descriptor

set.seed(12345) # to fix seed of stochastic (random-forest-based) imputation
x = missForest(x)$ximp # to impute missing values
y = y[donor_select]

#################################################################

version = "v1"
label = paste0("R",run_index)
result = eNetXplorer(x=x,y=y,family="gaussian",alpha=seq(0,1,by=0.1),
seed=123,nlambda.ext=1000,n_run=500, n_perm_null=100)
save(result,file=paste0(path_outfile,label,"_",version,".Robj"))
pdf.summary(result,path=path_outfile,filename=paste0(label,"_",version,".pdf"))


