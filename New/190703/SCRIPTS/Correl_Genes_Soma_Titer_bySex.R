library(ppcor)

rm(list=ls())

abs_path = "/Users/candiajm/ACTIVE/H5N1/190703/"
path_infile = paste0(abs_path,"DATA/")
path_outfile = paste0(abs_path,"RESULTS/")

pax_data = readRDS(paste0(path_infile,"PAX_dat0.in_pbmc.isv.0.7_8144g_170622.rds")) # 42 x 8144
gene = colnames(pax_data)
n_gene = length(gene)
sample_metadata_pax = readRDS(paste0(path_infile,"PAX_info0.in_170622.rds"))
donor_pax = paste0("s",substr(sample_metadata_pax[,1],7,8))

pbmc_data = readRDS(paste0(path_infile,"PBMC_dat0.in_isv.0.7_8144g_170216.rds")) # 41 x 8144 (s10 missing)
cat(paste0("Integrity check passed = ",sum(gene==colnames(pbmc_data))==n_gene),"\n")
sample_metadata_pbmc = readRDS(paste0(path_infile,"PBMC_info0.in_170216.rds"))
donor_pbmc = paste0("s",substr(sample_metadata_pbmc[,1],7,8))

soma_data = read.table(paste0(path_infile,"SOMA.txt"),header=T,stringsAsFactors=F,sep="\t")

sex = rep(0,nrow(soma_data))
sex[soma_data[,"sex"]=="Male Gende"] = 1 # F=0, M=1

pax_in_soma_ALL = donor_pax%in%soma_data[,"donor_short"]
soma_in_pax_ALL = soma_data[,"donor_short"]%in%donor_pax
pax_in_soma_ADJ = donor_pax%in%soma_data[soma_data[,"unblinded_flag"]=="ADJ","donor_short"]
soma_in_pax_ADJ = soma_data[,"donor_short"]%in%donor_pax[pax_in_soma_ADJ]
pbmc_in_soma_ALL = donor_pbmc%in%soma_data[,"donor_short"]
soma_in_pbmc_ALL = soma_data[,"donor_short"]%in%donor_pbmc
pbmc_in_soma_ADJ = donor_pbmc%in%soma_data[soma_data[,"unblinded_flag"]=="ADJ","donor_short"]
soma_in_pbmc_ADJ = soma_data[,"donor_short"]%in%donor_pbmc[pbmc_in_soma_ADJ]

target_soma = "RSPO3"
target_RFU = soma_data[,target_soma]
target_titer = "d28"
titer = soma_data[,"titer_I_MN_d28"]

# pax genes
RES_ALL = c("r","pVal","p_r","p_pVal")
for (i_gene in 1:n_gene) {
    RES = NULL
    correl = cor.test(pax_data[pax_in_soma_ADJ,i_gene],target_RFU[soma_in_pax_ADJ])
    RES = c(RES,correl$estimate,correl$p.value)
    partial_correl = pcor.test(pax_data[pax_in_soma_ADJ,i_gene],target_RFU[soma_in_pax_ADJ],sex[soma_in_pax_ADJ])
    RES = c(RES,partial_correl$estimate,partial_correl$p.value)
    RES_ALL = rbind(RES_ALL,RES)
}
output = cbind(c("gene",gene),RES_ALL)
write(t(output),ncol=ncol(output),file=paste0(path_outfile,"pax_ADJ_",target_soma,"_bySex.txt"),sep="\t")

# pbmc genes
RES_ALL = c("r","pVal","p_r","p_pVal")
for (i_gene in 1:n_gene) {
    RES = NULL
    correl = cor.test(pbmc_data[pbmc_in_soma_ADJ,i_gene],target_RFU[soma_in_pbmc_ADJ])
    RES = c(RES,correl$estimate,correl$p.value)
    partial_correl = pcor.test(pbmc_data[pbmc_in_soma_ADJ,i_gene],target_RFU[soma_in_pbmc_ADJ],sex[soma_in_pbmc_ADJ])
    RES = c(RES,partial_correl$estimate,partial_correl$p.value)
    RES_ALL = rbind(RES_ALL,RES)
}
output = cbind(c("gene",gene),RES_ALL)
write(t(output),ncol=ncol(output),file=paste0(path_outfile,"pbmc_ADJ_",target_soma,"_bySex.txt"),sep="\t")
