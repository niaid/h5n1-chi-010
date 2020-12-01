# PURPOSE: Correlation between individual markers from our signature and d42 titer in Bali cohort.

# rm(list=ls())
source(file.path("SCRIPTS/0_initialize.r"))
library(missForest)
library(eNetXplorer)

version = "v42"

path_infile = file.path(PROJECT_DIR, "DATA_PROCESSED/eNet")
path_outfile = file.path(PROJECT_DIR, "RESULTS/eNet", version)
dir.create(path_outfile, recursive = T, showWarnings = F)

data = read.table(file.path(path_infile,"eNet_InputData_r9.txt"), header=T, stringsAsFactors=F, sep="\t")
endpoint = "Day.42"

data[data[,"sex"]=="F","sex"] = 0
data[data[,"sex"]=="M","sex"] = 1

donor = data[,"SubjectID"]
n_donor = length(donor)

# Generate feature table.
predictor = c("sex","age",colnames(data)[-c(1:6)])
x = data[,predictor]

# Target    
y = data[,endpoint]

donor_select = rep(T,n_donor)

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

# Check if target is NA and remove the subject. 
na_idx <- is.na(y)
y <- y[!na_idx]
x <- x[!na_idx,]

result = eNetXplorer(x=x,y=y,family="gaussian",alpha=seq(0,1,by=0.1), seed=123,nlambda.ext=1000,n_run=1000, n_perm_null=250)
                                                                                # comment
# out to use default parameters.
save(result,file=file.path(path_outfile, paste0("emory_adj","_",version,".Robj")))
summaryPDF(result,path=path_outfile,filename=paste0("emory_adj","_",version,".pdf"))


