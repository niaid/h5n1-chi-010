# rm(list=ls())
source("SCRIPTS/0_initialize.r")
# PROJECT_DIR = "/Users/candiajm/ACTIVE/H5N1/WORKFLOW"
path_infile = file.path(PROJECT_DIR,"DATA_ORIGINAL","SOMAscan")
path_outfile = file.path(PROJECT_DIR,"DATA_PROCESSED","SOMAscan")

norm = c("Raw","Hyb","Hyb.MedNorm","Hyb.MedNorm.Cal","Hyb.Cal","Hyb.Cal.MedNorm")
n_norm = length(norm)

dilution = c("0.005","1","40")
n_dilution = length(dilution)

seq_id_remove = c("2795-23","3590-8","5071-3","5118-74","5073-30") # taken down by Somalogic as of 12/2016.

# somamers
tmp = t(as.matrix(read.table(paste0(path_infile,"/Somamers.txt"),header=F,sep="\t")))
colnames(tmp) = tmp[1,]
tmp = tmp[-1,]
seq_id = strsplit(tmp[,1],"_")
seq_id = matrix(unlist(seq_id),ncol=2,byrow=T)[,1] # to remove the sub-sequence info (after the underscore)
tmp[,1] = seq_id
somamer_remove = tmp[,1]%in%seq_id_remove
tmp = tmp[!somamer_remove,]
metadata_analyte=tmp

# RFU (already Hyb normalized)
tmp = as.matrix(read.table(paste0(path_infile,"/RFU.txt"),header=F,sep="\t"))
tmp = tmp[,!somamer_remove]
RFU_measured = tmp

# samples
tmp = as.matrix(read.table(paste0(path_infile,"/Samples.txt"),header=T,sep="\t"))
metadata_sample = tmp[,c(1,6,7)]

# we fix naming inconsistencies
metadata_sample[(metadata_sample[,3]=="CHI QC")|(metadata_sample[,3]=="QC_CHI")|(metadata_sample[,3]=="CHI_QC")|(metadata_sample[,3]=="QC_CHI1"),3] = "QC_CHI"
metadata_sample[metadata_sample[,3]=="QC_CHI",2] = "QC_CHI"
metadata_sample[metadata_sample[,3]=="sample",3] = "Sample"
plate_label = unique(metadata_sample[,1])
n_plate = length(plate_label)

# we remove HCE (data is already Hyb normalized, so we don't need them)
select_somamer = metadata_analyte[,9]!="Hybridization Control Elution"
metadata_analyte = metadata_analyte[select_somamer,]
n_analyte=nrow(metadata_analyte)
RFU_measured = RFU_measured[,select_somamer]

# we write metadata to output
write(t(rbind(colnames(metadata_analyte),metadata_analyte)),ncol=ncol(metadata_analyte),file=paste0(path_outfile,"/Somamers.txt"),sep="\t")
write(t(rbind(colnames(metadata_sample),metadata_sample)),ncol=ncol(metadata_sample),file=paste0(path_outfile,"/Samples.txt"),sep="\t")

RFU = vector("list",n_norm) # fluorescence intensity matrix

# 1. raw data
# RFU[[1]] = Not Available!

# 2. hybridization control normalization
RFU[[2]] = RFU_measured

# 3. we add median normalization
sample_type_id = metadata_sample[,3]
sample_type_id[sample_type_id=="QC"] = metadata_sample[sample_type_id=="QC",2]
RFU[[3]] = RFU[[2]]
for (i_plate in 1:n_plate) {
    sample_type = unique(sample_type_id[metadata_sample[,1]==plate_label[i_plate]])
    n_sample_type = length(sample_type)
    for (i_sample_type in 1:n_sample_type) {
        select_sample = (metadata_sample[,1]==plate_label[i_plate]) & (sample_type_id==sample_type[i_sample_type])
        for (i_dilution in 1:n_dilution) {
            select_somamer = metadata_analyte[,10] == dilution[i_dilution]
            data = RFU[[3]][select_sample,select_somamer,drop=F]
            data2 = data
            for (i in 1:ncol(data2)) {
                data2[,i] = data2[,i]/median(data2[,i])
            }
            for (i in 1:nrow(data)) {
                data[i,] = data[i,]/median(data2[i,])
            }
            RFU[[3]][select_sample,select_somamer] = data
        }
    }
}

# 4. we add inter-plate calibration
RFU[[4]] = RFU[[3]]
calibrator ="SL16689"
reference_all_plates = apply(RFU[[4]][metadata_sample[,2]==calibrator,],2,median)
for (i_plate in 1:n_plate) {
    sample_select = metadata_sample[,1]==plate_label[i_plate]
    reference = apply(RFU[[4]][sample_select&(metadata_sample[,2]==calibrator),],2,median)/reference_all_plates
    for (i_analyte in 1:n_analyte) {
        RFU[[4]][sample_select,i_analyte] = RFU[[4]][sample_select,i_analyte]/reference[i_analyte]
    }
}

# 5. We perform inter-plate calibration without median normalization (except for non-samples).
RFU[[5]] = RFU[[2]]
sample_select = metadata_sample[,3]!="Sample"
RFU[[5]][sample_select,] = RFU[[3]][sample_select,]
reference_all_plates = apply(RFU[[5]][metadata_sample[,2]==calibrator,],2,median)
for (i_plate in 1:n_plate) {
    sample_select = metadata_sample[,1]==plate_label[i_plate]
    reference = apply(RFU[[5]][sample_select&(metadata_sample[,2]==calibrator),],2,median)/reference_all_plates
    for (i_analyte in 1:n_analyte) {
        RFU[[5]][sample_select,i_analyte] = RFU[[5]][sample_select,i_analyte]/reference[i_analyte]
    }
}

# 6. we add multiplate median normalization on Samples.
RFU[[6]] = RFU[[5]]
select_sample = metadata_sample[,3]=="Sample"
for (i_dilution in 1:n_dilution) {
    select_somamer = metadata_analyte[,10] == dilution[i_dilution]
    data = RFU[[6]][select_sample,select_somamer,drop=F]
    data2 = data
    for (i in 1:ncol(data2)) {
        data2[,i] = data2[,i]/median(data2[,i])
    }
    for (i in 1:nrow(data)) {
        data[i,] = data[i,]/median(data2[i,])
    }
    RFU[[6]][select_sample,select_somamer] = data
}

# we write RFU data to output
for (i_norm in 2:n_norm) {
    write(t(RFU[[i_norm]]),ncol=ncol(RFU[[i_norm]]),
    file=paste0(path_outfile,"/",norm[i_norm],"_RFU.txt"),sep="\t")
}
