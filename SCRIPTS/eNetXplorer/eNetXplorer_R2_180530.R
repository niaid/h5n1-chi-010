# rm(list=ls())
source(file.path("SCRIPTS/0_initialize.r"))
library(missForest)
library(eNetXplorer)

version = "v4"

path_infile = file.path(PROJECT_DIR, "DATA_PROCESSED/eNet")
path_outfile = file.path(PROJECT_DIR, "RESULTS/eNet", version)
dir.create(path_outfile, recursive = T,  showWarnings = F)

data = read.table(file.path(path_infile,"eNet_InputData_r2.txt"), header=T, stringsAsFactors=F, sep="\t")
run_index = 2 # integer value in the 1:9 range
endpoint = "titer_I_MN_d28"

data[data[,"sex"]=="F","sex"] = 0
data[data[,"sex"]=="M","sex"] = 1

donor = data[,"donor_short"]
n_donor = length(donor)

predictor_rm = c("Gp09", "Gp10")

#################################################################
# choose parameter set to reproduce results reported in the paper
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
    model = "FP_GPgset"
    donor_group = "NONADJ"
    modifier = ""
} else if (run_index==5) {
    model = "FP_GPgset"
    donor_group = "ADJ"
    modifier = "wo44"
} else if (run_index==6) {
    model = "FB_GBpax_GBpbmc"
    donor_group = "NEG"
    modifier = ""
} else if (run_index==7) {
    model = "FB_GBpax_GBpbmc"
    donor_group = "NONADJ"
    modifier = ""
} else if (run_index==8) {
    model = "SOMAscan"
    donor_group = "ADJ"
    modifier = "wo44"
} else if (run_index==9) {
    model = "SOMAscan_wCounts"
    donor_group = "ADJ"
    modifier = "wo44"
}
#################################################################


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
if (model=="SOMAscan_wCounts") {
    predictor = c("sex","age",colnames(data)[171:1482]) # adds 7 columns with cell subset counts at the end
}

pred_rm_index = predictor %in% predictor_rm
predictor = predictor[!pred_rm_index]

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

label = paste0("R",run_index)
result = eNetXplorer(x=x,y=y,family="gaussian",alpha=seq(0,1,by=0.1), seed=123) # ,nlambda.ext=1000,n_run=1000, n_perm_null=250)
                                                                                # Comment
# out from here to use default parameters.
save(result,file=file.path(path_outfile, paste0(label,"_",version,".Robj")))
summaryPDF(result,path=path_outfile,filename=paste0(label,"_",version,".pdf"))


