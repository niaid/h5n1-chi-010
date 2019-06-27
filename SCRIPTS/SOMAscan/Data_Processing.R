library(calibrate)

# rm(list=ls())

# PROJECT_DIR = "/Users/candiajm/ACTIVE/H5N1/WORKFLOW"

setwd(file.path(PROJECT_DIR,"DATA_PROCESSED","SOMAscan"))

norm = "Hyb.Cal.MedNorm"

sample = as.matrix(read.table("Samples.txt",header=T,stringsAsFactors=F,sep="\t"))
somamer =  as.matrix(read.table("Somamers.txt",header=T,stringsAsFactors=F,sep="\t"))
RFU = as.matrix(read.table(paste0(norm,"_RFU.txt"),header=F,stringsAsFactors=F,sep="\t"))

# We focus on samples only.
sample_select = sample[,3]=="Sample"
sample = sample[sample_select,]
sample[,2] = paste0("s",substr(sample[,2],7,8)) # to simplify sample labels
RFU = RFU[sample_select,]

# We average duplicates.
tmp = table(sample[,2])
duplicate = names(tmp[tmp==2])
n_duplicate = length(duplicate)

select_sample = !sample[,2]%in%duplicate
tmp_sample = sample[select_sample,]
tmp_RFU = RFU[select_sample,]
for (i_duplicate in 1:n_duplicate) {
    tmp_sample = rbind(tmp_sample,sample[which(sample[,2]==duplicate[i_duplicate])[1],])
    tmp_RFU = rbind(tmp_RFU,apply(RFU[sample[,2]==duplicate[i_duplicate],],2,mean))
}
sample = tmp_sample
RFU = tmp_RFU

RFU = log10(RFU) # we log the RFUs

# We save the data for downstream analysis
write(t(RFU),ncol=ncol(RFU),file="Final_RFU.txt",sep="\t")
write(t(rbind(colnames(sample),sample)),ncol=ncol(sample),file="Final_Samples.txt",sep="\t")
write(t(rbind(colnames(somamer),somamer)),ncol=ncol(somamer),file="Final_Somamers.txt",sep="\t")

