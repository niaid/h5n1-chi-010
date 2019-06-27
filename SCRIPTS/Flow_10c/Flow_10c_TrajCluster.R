require(gplots)
require(RColorBrewer)
library(grid)
library(gridExtra)
library(fpc)
library(dynamicTreeCut)

rm(list=ls())

source(file.path("SCRIPTS/0_initialize.r"))
path_infile1 = file.path(PROJECT_DIR,"DATA_PROCESSED","Flow_10c")
path_infile2 = file.path(PROJECT_DIR,"RESULTS","Flow_10c")
path_outfile = file.path(PROJECT_DIR,"RESULTS","Flow_10c")

################################
measure_type = "PercentOfParent"
#measure_type = "PercentOfTotal"
thres = 0.50
donor_thres = 0.30
################################

cell_type = c("Bcell","DCMonoNK","Tcell","Thelper","Treg")
n_cell_type = length(cell_type)

TP = c("d0_4h","d0_12h","d1","d7","d21_0h","d21_4h","d21_12h","d22","d28","d42") # we remove baseline TP
n_TP = length(TP)

# we read the trajectory matrix
traj_matrix = data.matrix(read.table(paste(path_infile2,"/TrajMatrix_DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,".txt",sep=""),header=F,sep="\t"))
cellpop_donor = read.table(paste(path_infile2,"/TrajMatrixLabels_DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,".txt",sep=""),header=F,colClass="character",stringsAsFactors=F,sep="\t")
cellpop = unique(cellpop_donor[,2])
n_cellpop = length(cellpop)
donor = unique(cellpop_donor[,3])
n_donor = length(donor)

# we add cellpop labels
cellpop_labels = cellpop
for (i_cell_type in 1:n_cell_type) {
    dataset_labels = read.table(paste(path_infile1,"/GateList_",cell_type[i_cell_type],".txt",sep=""),header=F,colClass="character",stringsAsFactors=F,check.names=F,sep="\t")
    for (i_cellpop in 1:n_cellpop) {
        index = which(dataset_labels[,1]==cellpop[i_cellpop])
        if (length(index)>0) {
            cellpop_labels[i_cellpop] = dataset_labels[index,2]
        }
    }
}

# we calculate distances between trajectory pairs
method = "CORpearson" # other methods: DTW, INT.PER, etc
clustering_method = "cutreeHybrid"
n_instance = nrow(traj_matrix)
similarity_matrix = matrix(rep(0,n_instance**2),ncol=n_instance)
for (i_instance in 1:n_instance) {
    for (j_instance in 1:i_instance) {
        x = traj_matrix[i_instance,]
        y = traj_matrix[j_instance,]
        select = (is.finite(x))&(is.finite(y)) # to check for NA's and Inf's
        x = x[select]
        y = y[select]
        similarity_matrix[i_instance,j_instance] = sqrt(2*(1.-cor(x,y,method="pearson")))
        similarity_matrix[j_instance,i_instance] = similarity_matrix[i_instance,j_instance]
    }
}

hc <- hclust(as.dist(similarity_matrix), method="average")

# we run several instances with different numbers and choose among them using silhouette
deepSplit_param = c(0,1,2,3,4,0,1,2,3,4)
minClusterSize_param = round(c(0.05,0.05,0.05,0.05,0.05,0.1,0.1,0.1,0.1,0.1)*n_instance)
n_param = length(deepSplit_param)
clust_id = vector("list",n_param)
output = c("deepSplit","minClusterSize","n_clust","intra","inter","ratio","sil","pearsongamma","dunn2","ch")
for (i_param in 1:n_param) {
    treeClust = cutreeHybrid(dendro=hc,distM=similarity_matrix,deepSplit=deepSplit_param[i_param],minClusterSize = minClusterSize_param[i_param])
    clust_id[[i_param]] = treeClust$labels # the "noise cluster" is cluster ID=0
    n_clust = max(clust_id[[i_param]])
    if (n_clust>1) {
        select = clust_id[[i_param]]!=0
        cluster_stats = cluster.stats(d=similarity_matrix[select,select],clustering=clust_id[[i_param]][select],silhouette=T,G2=F,G3=F)
        intra = cluster_stats$average.within
        inter = cluster_stats$average.between
        ratio = cluster_stats$wb.ratio # intra/inter
        sil = cluster_stats$avg.silwidth # avg silhouette
        pearsongamma = cluster_stats$pearsongamma
        dunn2 = cluster_stats$dunn2
        ch = cluster_stats$ch
        output = rbind(output,c(deepSplit_param[i_param],minClusterSize_param[i_param],n_clust,intra,inter,ratio,sil,pearsongamma,dunn2,ch))
    }
}
outfile = paste(path_outfile,"/DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,"_",method,"_",clustering_method,"_ClusterMeasures.txt",sep="")
write(t(output),ncol=ncol(output),file=outfile)

outfile = paste(path_outfile,"/DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,"_",method,"_",clustering_method,"_ClusterMeasures.pdf",sep="")
pdf(outfile,width=7,height=5)
x_label = paste(output[-1,1],output[-1,2],sep="_")
y_label = output[1,-c(1,2)]
output = output[-1,-c(1,2)]
x_plot = 1:nrow(output)
for (i_meas in 1:ncol(output)) {
    y_plot = as.numeric(output[,i_meas])
    plot(range(x_plot),range(y_plot),xaxt="n",type="n",xlab="parameters [deepSplit_minClusterSize]",ylab=y_label[i_meas],cex.lab=0.85,main="")
    points(x_plot,y_plot,cex=1,col=1,pch=16,lty=1,type="b")
    axis(side = 1, at = x_plot, labels = x_label, las = 2, cex.axis=1)
}
dev.off()

run_select = 1 # we select the first set of clustering parameters
clust_id = clust_id[[run_select]]
clust_out = data.frame(node = cellpop_donor$V4, cluster=clust_id)
write.table(clust_out, file=paste(path_outfile,"/DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,"_",method,"_",clustering_method,"_clusters.txt",sep=""),
            sep="\t", col.names=T, row.names=F, quote = F)
n_clust = max(clust_id)
delta_median = matrix(rep(0,n_clust*n_TP),ncol=n_TP)
delta_mad = matrix(rep(0,n_clust*n_TP),ncol=n_TP)
clust_size = rep(0,n_clust)
for (i_clust in 1:n_clust) {
    select = clust_id==i_clust
    delta_median[i_clust,1:n_TP] = apply(traj_matrix[select,],2,median,na.rm=T)
    delta_mad[i_clust,1:n_TP] = apply(traj_matrix[select,],2,mad,na.rm=T)
    clust_size[i_clust] = sum(select)
}

x = 1:n_TP
eps = 0.2
range_x = range(x)
if (measure_type=="PercentOfParent") {ylab = "(TP-TP0)/median(TP) (% parent)"}
if (measure_type=="PercentOfTotal") {ylab = "(TP-TP0)/median(TP) (% total)"}

outfile = paste(path_outfile,"/DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,"_",method,"_",clustering_method,"_deepSplit",deepSplit_param[run_select],"_minClusterSize",minClusterSize_param[run_select],".pdf",sep="")
pdf(outfile)

myCol = rev(redblue(16)) # blue-white-red
for (i_clust in 1:n_clust) {
    select = clust_id==i_clust
    heatmap_panel = traj_matrix[select,]
    rownames(heatmap_panel) = cellpop_donor[select,5]
    colnames(heatmap_panel) = TP
    myBreaks = quantile(abs(heatmap_panel),probs=c(seq(0.40,0.90,by=0.10),0.95,1),na.rm=T)
    myBreaks = c(-rev(myBreaks),0,myBreaks)
    heatmap.2(heatmap_panel, scale="none", Rowv=T, Colv=F, na.rm=T, na.color="black",
    col = myCol, breaks = myBreaks, dendrogram = "none", cexRow=0.1, cexCol=0.6, trace="none", key=T)
}

# we plot the trajectory clusters by median/mad
range_y = range(c(delta_median-delta_mad,delta_median+delta_mad))
title = paste0("Cluster size: ",clust_size[1])
for (i_clust in 2:n_clust) {
    title = paste0(title,", ",clust_size[i_clust])
}
plot(range_x,range_y,xaxt="n",type="n",xlab="",ylab=ylab,cex.lab=1.1,main=title)
for (i_clust in 1:n_clust) {
    offset = 2*eps*(i_clust-1)/(n_clust-1)-eps
    points(x+offset,delta_median[i_clust,],cex=1,col=i_clust,pch=16,lty=2,type="b")
    arrows(x+offset,delta_median[i_clust,]-delta_mad[i_clust,],x+offset,delta_median[i_clust,]+delta_mad[i_clust,],col=i_clust,length=0.025,angle=90,code=3)
}
axis(side = 1, at = x, labels = TP, las = 2, cex.axis=0.8)
abline(h=0,col="black",lty=3)

# we clean up cluster IDs
remove = clust_id==0
# if (sum(remove)>0) {
    clust_id = clust_id[!remove]
    cellpop_donor2 = cellpop_donor[!remove,]
# }
cellpop2 = unique(cellpop_donor2[,2])
n_cellpop2 = length(cellpop2)
cellpop_labels2 = rep("",n_cellpop2)
cellpop_labels2n = rep("",n_cellpop2)
for (i in 1:n_cellpop2) {
    cellpop_labels2[i] = cellpop_labels[which(cellpop==cellpop2[i])]
    cellpop_labels2n[i] = paste0("ID",cellpop2[i],": ",cellpop_labels[which(cellpop==cellpop2[i])])
}
donor2 = unique(cellpop_donor2[,3])
n_donor2 = length(donor2)

# we generate heatmap
cluster_heatmap = matrix(rep(NA,n_donor2*n_cellpop2),ncol=n_cellpop2)
rownames(cluster_heatmap) = donor2
colnames(cluster_heatmap) = cellpop_labels2
for (i in 1:nrow(cellpop_donor2)) {
    donor_index = which(donor2==cellpop_donor2[i,3])
    cellpop_index = which(cellpop2==cellpop_donor2[i,2])
    cluster_heatmap[donor_index,cellpop_index] = clust_id[i]
}
hm <- try(
heatmap.2(cluster_heatmap, scale="none", Rowv=T, Colv=T, na.color="white", col=1:n_clust,
dendrogram = "both", margins=c(5,5), cexRow=0.4, cexCol=0.35, key=F, trace="none",main="")
)
if(inherits(hm,"try-error")) { # to catch errors due to clustering on a very sparse matrix
    try(
    heatmap.2(cluster_heatmap, scale="none", Rowv=F, Colv=F, na.color="white", col=1:n_clust,
    dendrogram = "none", margins=c(5,5), cexRow=0.4, cexCol=0.35, key=F, trace="none",main="")
    )
}

donor_clust_score = matrix(rep(0,n_donor*n_clust),ncol=n_clust)
rownames(donor_clust_score) = donor
colnames(donor_clust_score) = paste0("Fp0",1:n_clust)
module = matrix(rep(0,n_cellpop2*n_clust),ncol=n_clust)
rownames(module) = cellpop_labels2
colnames(module) = paste0("Fp0",1:n_clust)
for (i_clust in 1:n_clust) {
    cluster_heatmap2 = matrix(rep(0,n_donor2*n_cellpop2),ncol=n_cellpop2)
    index = which(cluster_heatmap==i_clust,arr.ind=T)
    for (i in 1:nrow(index)) {
        cluster_heatmap2[index[i,1],index[i,2]]=1
    }
    n_instance = ncol(cluster_heatmap2)
    similarity_matrix = matrix(rep(0,n_instance**2),ncol=n_instance)
    for (i_instance in 1:n_instance) {
        for (j_instance in 1:i_instance) {
            x1 = cluster_heatmap2[,i_instance]
            x2 = cluster_heatmap2[,j_instance]
            similarity_matrix[i_instance,j_instance] = sum(abs(x1-x2))/n_donor2
            similarity_matrix[j_instance,i_instance] = similarity_matrix[i_instance,j_instance]
        }
    }
    KM_id = kmeans(similarity_matrix,2,nstart=50)$cluster
    
    fraction_KM_id1 = sum(cluster_heatmap2[,KM_id==1])/sum(KM_id==1)
    fraction_KM_id2 = sum(cluster_heatmap2[,KM_id==2])/sum(KM_id==2)
    if (fraction_KM_id1>fraction_KM_id2) {
        KM_id_target = 1
    } else {
        KM_id_target = 2
    }
    CP_clust = KM_id==KM_id_target
    # module components
    module[,i_clust] = CP_clust*i_clust
    cellpop_labels2[CP_clust] # labels for CPs in module
    cellpop_ref = cellpop2[CP_clust]
    x_ref = delta_median[i_clust,]
    for (i_donor in 1:n_donor) {
        select = (cellpop_donor[,3]==donor[i_donor])&(cellpop_donor[,2]%in%cellpop_ref)
        if (sum(select)>=1) {
            x_traj = apply(traj_matrix[select,,drop=F],2,mean)
        } else {
            x_traj = rep(0,ncol(traj_matrix))
        }
        select_cor = (is.finite(x_traj))&(is.finite(x_ref))
        donor_clust_score[i_donor,i_clust] = cor(x_traj[select_cor],x_ref[select_cor],method="pearson")
    }
}

module = module[apply(module,1,sum)>0,] # to remove CPs that don't participate from any module
module2 = module
rownames(module2) = cellpop_labels2n
heatmap.2(module2, scale="none", Rowv=T, Colv=F, col=0:n_clust,
dendrogram = "none", margins=c(10,20), cexRow=0.5, labCol=F,trace="none", key=F,
add.expr=mtext(side = 1, text = colnames(module2), at = 1:ncol(module2), las = 2, line = 0.5,col = 1:ncol(module2), cex = 0.8))

# heatmap with donor vs cluster scores
myCol = rev(redblue(16)) # blue-white-red
myBreaks = c(seq(0,1,len=length(myCol)+1))
heatmap.2(donor_clust_score, scale="none", Rowv=T, Colv=F, na.rm=T, na.color="black",
col = myCol, breaks = myBreaks, dendrogram = "row", cexRow=0.5, labCol=F,trace="none", key=T,
add.expr=mtext(side = 1, text = colnames(donor_clust_score), at = 1:ncol(donor_clust_score), las = 2, line = 0.5,col = 1:ncol(donor_clust_score), cex = 0.8))

outfile = paste(path_outfile,"/DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,"_",method,"_",clustering_method,"_deepSplit",deepSplit_param[run_select],"_minClusterSize",minClusterSize_param[run_select],"_DonorScores.txt",sep="")
output = rbind(c("Donor",colnames(donor_clust_score)),cbind(rownames(donor_clust_score),donor_clust_score))
write(t(output),ncol=ncol(output),file=outfile,sep=",")

#########################################################
# trajectory clusters (filtered by CP and donor)
cor_thres = 0.50
delta_median_filt = matrix(rep(0,n_clust*n_TP),ncol=n_TP)
delta_mad_filt = matrix(rep(0,n_clust*n_TP),ncol=n_TP)
clust_size_filt = rep(0,n_clust)
donor_clust_score_filt = donor_clust_score*0
for (i_clust in 1:n_clust) {
    CP_this_clust = rownames(module[module[,i_clust]>0,,drop=F])
    select = donor_clust_score[,i_clust]>cor_thres
    select = (select==T)&(!is.na(select))
    donor_clust_score_filt[select,i_clust]=i_clust
    donor_this_clust = rownames(donor_clust_score[select,])
    traj_matrix_filt = NULL
    for (i_CP in 1:length(CP_this_clust)) {
        for (i_donor in 1:length(donor_this_clust)) {
            select = (cellpop_donor[,2]==cellpop[which(cellpop_labels==CP_this_clust[i_CP])])&(cellpop_donor[,3]==donor_this_clust[i_donor])
            if (sum(select)>0) {
                traj_matrix_filt = rbind(traj_matrix_filt,traj_matrix[select,])
            }
        }
    }
    delta_median_filt[i_clust,1:n_TP] = apply(traj_matrix_filt,2,median,na.rm=T)
    delta_mad_filt[i_clust,1:n_TP] = apply(traj_matrix_filt,2,mad,na.rm=T)
    clust_size_filt[i_clust] = nrow(traj_matrix_filt)
}

#####################################################

# we plot the trajectory clusters by median/mad
range_y = range(c(delta_median-delta_mad,delta_median+delta_mad))
ncol = min((n_clust+1)%/%2,3)
nrow = n_clust%/%ncol+1
par(mfrow=c(nrow,ncol),oma=c(0,0,5,0))
for (i_clust in 1:n_clust) {
    leftmost_col = i_clust%%ncol==1
    bottom_row = i_clust>((nrow-1)*ncol)
    if (leftmost_col) {ylabel=ylab} else {ylabel=""}
    if (bottom_row) {xlabel="TP"} else {xlabel=""}
    plot(range_x,range_y,xaxt="n",yaxt="n",type="n",xlab=xlabel,ylab=ylabel,cex.lab=0.95)
    if (leftmost_col) {
        axis(side=2,cex.axis=0.6)
    } else {
        #axis(side=2,labels=F)
        axis(side=2,cex.axis=0.6) # we include tick labels regardless of quadrant
    }
    if (bottom_row) {
        axis(side=1,at=x,labels=TP,las=2,cex.axis=0.6)
    } else {
        #axis(side=1,at=x,labels=F)
        axis(side=1,at=x,labels=TP,las=2,cex.axis=0.6) # we include tick labels regardless of quadrant
    }
    offset = 0
    points(x+offset,delta_median[i_clust,],cex=1,col=i_clust,pch=16,lty=2,type="b")
    arrows(x+offset,delta_median[i_clust,]-delta_mad[i_clust,],x+offset,delta_median[i_clust,]+delta_mad[i_clust,],col=i_clust,length=0.025,angle=90,code=3)
    abline(h=0,col="black",lty=3)
}
title = paste0("Cluster size: ",clust_size[1])
for (i_clust in 2:n_clust) {
    title = paste0(title,", ",clust_size[i_clust])
}
title(title,outer=T)

# we plot the trajectory clusters by median/mad (filtered by CP/donor)
ncol = min((n_clust+1)%/%2,3)
nrow = ceiling(n_clust/ncol)
par(mfrow=c(nrow,ncol))
for (i_clust in 1:n_clust) {
    # each plot is independently rescaled
    range_y = range(c(delta_median_filt[i_clust,]-delta_mad_filt[i_clust,],delta_median_filt[i_clust,]+delta_mad_filt[i_clust,]))
    leftmost_col = i_clust%%ncol==1
    bottom_row = i_clust>((nrow-1)*ncol)
    if (leftmost_col) {ylabel=ylab} else {ylabel=""}
    if (bottom_row) {xlabel="TP"} else {xlabel=""}
    plot(range_x,range_y,xaxt="n",yaxt="n",type="n",xlab=xlabel,ylab=ylabel,cex.lab=0.95)
    if (leftmost_col) {
        axis(side=2,cex.axis=0.6)
    } else {
        #axis(side=2,labels=F)
        axis(side=2,cex.axis=0.6) # we include tick labels regardless of quadrant
    }
    if (bottom_row) {
        axis(side=1,at=x,labels=TP,las=2,cex.axis=0.6)
    } else {
        #axis(side=1,at=x,labels=F)
        axis(side=1,at=x,labels=TP,las=2,cex.axis=0.6) # we include tick labels regardless of quadrant
    }
    offset = 0
    points(x+offset,delta_median_filt[i_clust,],cex=1,col=i_clust,pch=1,lty=3,type="b")
    arrows(x+offset,delta_median_filt[i_clust,]-delta_mad_filt[i_clust,],x+offset,delta_median_filt[i_clust,]+delta_mad_filt[i_clust,],col=i_clust,length=0.025,angle=90,code=3)
    abline(h=0,col="black",lty=3)
}
title = paste0("Cluster size: ",clust_size_filt[1])
for (i_clust in 2:n_clust) {
    title = paste0(title,", ",clust_size_filt[i_clust])
}
title = paste0(title," [filtered by CP/donor; each plot independently rescaled]")
title(title,outer=T)

dev.off()

# We save processed data for downstream analysis

# 1. Cell Populations vs Module
outfile = paste(path_outfile,"/DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,"_",method,"_",clustering_method,"_deepSplit",deepSplit_param[run_select],"_minClusterSize",minClusterSize_param[run_select],"_CPs.txt",sep="")
CP_desc = rownames(module)
CP_code = rep(NA,length(CP_desc))
for (i_CP in 1:length(CP_desc)) {
    CP_code[i_CP] = cellpop[which(cellpop_labels==CP_desc[i_CP])]
}
output = rbind(c("CP_desc","CP_code",colnames(module)),cbind(CP_desc,CP_code,module))
write(t(output),ncol=ncol(output),file=outfile,sep=",")

# 2. Donors vs Module
outfile = paste(path_outfile,"/DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,"_",method,"_",clustering_method,"_deepSplit",deepSplit_param[run_select],"_minClusterSize",minClusterSize_param[run_select],"_Donors.txt",sep="")
donor_desc = rownames(donor_clust_score_filt)
output = rbind(c("Donor_desc",colnames(donor_clust_score_filt)),cbind(donor_desc,donor_clust_score_filt))
write(t(output),ncol=ncol(output),file=outfile,sep=",")

# 3. Module profiles
outfile = paste(path_outfile,"/DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,"_",method,"_",clustering_method,"_deepSplit",deepSplit_param[run_select],"_minClusterSize",minClusterSize_param[run_select],"_Modules.txt",sep="")
output = rbind(c("timepoint",colnames(module)),cbind(TP,t(delta_median_filt)))
write(t(output),ncol=ncol(output),file=outfile,sep=",")

