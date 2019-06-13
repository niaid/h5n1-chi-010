require(gplots)
require(RColorBrewer)

rm(list=ls())

source(file.path("SCRIPTS/0_initialize.r"))
path_infile = file.path(PROJECT_DIR,"DATA_PROCESSED","Flow_10c")
path_outfile = file.path(PROJECT_DIR,"RESULTS","Flow_10c")
dir.create(path_outfile, showWarnings = F)

################################
measure_type = "PercentOfParent"
#measure_type = "PercentOfTotal"
thres = 0.50
donor_thres = 0.30
################################

cell_type = c("Bcell","DCMonoNK","Tcell","Thelper","Treg")

infile_root = c("H5N1_Bcell-H5N1_Bcell","H5N1_DCMonNK-DCMonNK",
"H5N1_Tcell-H5N1-Tcell","H5N1_Thelper-Thelper3","H5N1_Treg-H5N1-Treg")
n_infile = length(infile_root)

if (measure_type == "PercentOfParent") {
    infile = paste0("/",infile_root,"_freq.V1.1")
}
if (measure_type == "PercentOfTotal") {
    infile = paste0("/",infile_root,"_counts.V1.1")
}

donor = paste0("H5N1-0",c(paste0("0",1:9),paste0("1",0:9),paste0("2",0:6),
paste0("2",8:9),paste0("3",0:9),"40","42","43","44"))
donor_short = paste0("s",c(paste0("0",1:9),paste0("1",0:9),paste0("2",0:6),
        paste0("2",8:9),paste0("3",0:9),"40","42","43","44"))
n_donor = length(donor)
TP = c("d0_0h","d0_4h","d0_12h","d1","d7","d21_0h","d21_4h","d21_12h","d22","d28","d42")
n_TP = length(TP)

##### For use in "PercentOfTotal" mode #####
cellpop_ref = vector("list",n_infile)
cellpop_ref[[1]] = "LYM"
cellpop_ref[[2]] = "MNC"
cellpop_ref[[3]] = "LYM"
cellpop_ref[[4]] = "LYM"
cellpop_ref[[5]] = "LYM"
cellpop_remove = vector("list",n_infile)
cellpop_remove[[2]] = c("CD45+")
############################################

traj_matrix_ALL = NULL
traj_labels_ALL = NULL
traj_matrix_filtered = NULL
traj_labels_filtered = NULL
for (i_infile in 1:n_infile) {
    
    dataset = read.table(paste0(path_infile,infile[i_infile]),header=T,stringsAsFactors=F,check.names=F,sep="\t")
    dataset = dataset[dataset[,1]!="CHI-002",] # we remove the bridge sample
    dataset = dataset[,(colnames(dataset)!="58")&(colnames(dataset)!="58 (MFI)")] # we remove cellpop 58 (DCMono)
    cellpop = colnames(dataset)[-(1:2)]
    n_cellpop = length(cellpop)
    
    # we "fix" TP labels
    TP_original = c("d8","d29","d31","d43")
    TP_fixed = c("d7","d28","d28","d42")
    for (i in 1:length(TP_original)) {
        select = dataset[,2]==TP_original[i]
        if (sum(select)>0) {
            dataset[select,2]=rep(TP_fixed[i],sum(select))
        }
    }
    
    if (measure_type == "PercentOfTotal") {
        # we transform counts in percent of total
        col_ref = which(colnames(dataset)==cellpop_ref[[i_infile]])
        col_remove = which(colnames(dataset)%in%cellpop_remove[[i_infile]])
        dataset2 = dataset[,1:2]
        cellpop2 = NULL
        for (i_col in 3:ncol(dataset)) {
            if (!(i_col%in%c(col_ref,col_remove))) {
                dataset2 = cbind(dataset2,100*dataset[,i_col]/dataset[,col_ref])
                cellpop2 = c(cellpop2,colnames(dataset)[i_col])
            }
        }
        dataset = dataset2
        cellpop = cellpop2
        n_cellpop = length(cellpop)
    }
    
    # we filter out CPs
    for (i_cellpop in 1:n_cellpop) {
        tmp = dataset[,c(1,2,2+i_cellpop)]
        traj_matrix = matrix(rep(NA,n_donor*(n_TP-1)),ncol=(n_TP-1))
        traj_labels = matrix(rep("",n_donor*5),ncol=5)
        traj_flag = rep(F,n_donor)
        for (i_donor in 1:n_donor) {
            tmp2 = tmp[tmp[,1]==donor[i_donor],]
            baseline = tmp2[tmp2[,2]==TP[1],3]
            median_ref = median(tmp2[,3],na.rm=T)
            if (is.na(baseline)||is.na(median_ref)||median_ref<0.10) {  # to avoid illegal divisions
                traj_flag[i_donor] = T
            }
            for (i_TP in 2:n_TP) {
                select = tmp2[,2]==TP[i_TP]
                if (sum(select)>0) {
                    traj_matrix[i_donor,i_TP-1] = (mean(tmp2[select,3])-baseline)/median_ref # NOTE: we normalize relative to the median
                }
            }
            traj_labels[i_donor,] = c(cell_type[i_infile],cellpop[i_cellpop],donor_short[i_donor],paste(cellpop[i_cellpop],donor_short[i_donor],sep="_"),paste(cell_type[i_infile],cellpop[i_cellpop],donor_short[i_donor],sep="_"))
        }
        select = !traj_flag
        if (sum(select)>1) {
            traj_matrix = traj_matrix[select,]
            traj_labels = traj_labels[select,]
            traj_matrix_ALL = rbind(traj_matrix_ALL,traj_matrix)
            traj_labels_ALL = rbind(traj_labels_ALL,traj_labels)
            cellpop_select = F
            for (i_col in 1:ncol(traj_matrix)) {
                if (((sum(traj_matrix[,i_col]<(-thres),na.rm=T)/n_donor)>donor_thres)|((sum(traj_matrix[,i_col]>thres,na.rm=T)/n_donor)>donor_thres)) {
                    cellpop_select = T
                    break
                }
            }
            if (cellpop_select) {
                traj_matrix_filtered = rbind(traj_matrix_filtered,traj_matrix)
                traj_labels_filtered = rbind(traj_labels_filtered,traj_labels)
            }
        }
    }
}

# we save the trajectory matrix
ver = "_v3"
outfile = paste(path_outfile,"/TrajMatrix_AllCPs_",measure_type,ver,".txt",sep="")
write(t(traj_matrix_ALL),ncol=ncol(traj_matrix_ALL),file=outfile,sep="\t")
outfile = paste(path_outfile,"/TrajMatrixLabels_AllCPs_",measure_type,ver,".txt",sep="")
write(t(traj_labels_ALL),ncol=ncol(traj_labels_ALL),file=outfile,sep="\t")

outfile = paste(path_outfile,"/TrajMatrix_DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,ver,".txt",sep="")
write(t(traj_matrix_filtered),ncol=ncol(traj_matrix_filtered),file=outfile,sep="\t")
outfile = paste(path_outfile,"/TrajMatrixLabels_DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,ver,".txt",sep="")
write(t(traj_labels_filtered),ncol=ncol(traj_labels_filtered),file=outfile,sep="\t")
