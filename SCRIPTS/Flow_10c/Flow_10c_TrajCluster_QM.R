require(gplots)
require(RColorBrewer)
library(grid)
library(gridExtra)
library(fpc)
library(dynamicTreeCut)
library(ComplexHeatmap)
# rm(list=ls())

source(file.path("SCRIPTS/0_initialize.r"))
path_infile1 = file.path(PROJECT_DIR,"DATA_PROCESSED","Flow_10c")
path_infile2 = file.path(PROJECT_DIR,"RESULTS","Flow_10c")
path_outfile = file.path(PROJECT_DIR,"RESULTS","Flow_10c")

dn.fig = file.path(PROJECT_DIR, "FIGURES/Flow_QC")
dir.create(dn.fig, showWarnings = F)

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
# outfile = paste(path_outfile,"/DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,"_",method,"_",clustering_method,"_ClusterMeasures.txt",sep="")
# write(t(output),ncol=ncol(output),file=outfile)

# silhouette analysis
select = clust_id[[1]]!=0
sil_comb = cluster::silhouette(clust_id[[1]][select], dmatrix=similarity_matrix[select,select])
sil_comb.asw = as.data.frame(sil_comb[,]) %>% 
  arrange(cluster, desc(sil_width)) %>% 
  mutate(x = seq_along(cluster)) %>% 
  mutate(cluster = factor(cluster)) %>% 
  group_by(cluster) %>% 
  mutate(ASW = mean(sil_width)) %>% 
  ungroup()
sil_comb.df = data.frame(pattern = 1:8, ASW = summary(sil_comb)$clus.avg.width) %>% 
  mutate(pattern = paste0("Fp0",pattern))
cm.pat = palette()
# factoextra::fviz_silhouette(sil_comb, label = FALSE, print.summary = F) +
ggplot(sil_comb.asw, aes(x, sil_width, fill=cluster, col=cluster)) +
  geom_col() +
  scale_color_manual(values = cm.pat) +
  scale_fill_manual(values = cm.pat) +
  guides(col=F, fill=F) + 
  xlab("") +
  ylab("Silhouette width") +
  geom_errorbar(data=sil_comb.asw, aes(x, ymax = ASW, ymin = ASW, group = cluster),
                size=0.5, linetype = "longdash", inherit.aes = F, width = 1) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())
fn.fig = file.path(dn.fig, "Flow_patterns_silhouette")
ggsave(paste0(fn.fig, ".png"), w=5, h=2.5)
ggsave(paste0(fn.fig, ".pdf"), w=5, h=2.5)

# outfile = paste(path_outfile,"/DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,"_",method,"_",clustering_method,"_ClusterMeasures.pdf",sep="")
# pdf(outfile,width=7,height=5)
# x_label = paste(output[-1,1],output[-1,2],sep="_")
# y_label = output[1,-c(1,2)]
# output = output[-1,-c(1,2)]
# x_plot = 1:nrow(output)
# for (i_meas in 1:ncol(output)) {
#     y_plot = as.numeric(output[,i_meas])
#     plot(range(x_plot),range(y_plot),xaxt="n",type="n",xlab="parameters [deepSplit_minClusterSize]",ylab=y_label[i_meas],cex.lab=0.85,main="")
#     points(x_plot,y_plot,cex=1,col=1,pch=16,lty=1,type="b")
#     axis(side = 1, at = x_plot, labels = x_label, las = 2, cex.axis=1)
# }
# dev.off()

run_select = 1 # we select the first set of clustering parameters
clust_id = clust_id[[run_select]]
clust_out = data.frame(node = cellpop_donor$V4, cluster=clust_id)
# write.table(clust_out, file=paste(path_outfile,"/DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,"_",method,"_",clustering_method,"_clusters.txt",sep=""),
#             sep="\t", col.names=T, row.names=F, quote = F)
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

# x = 1:n_TP
# eps = 0.2
# range_x = range(x)
# if (measure_type=="PercentOfParent") {ylab = "(TP-TP0)/median(TP) (% parent)"}
# if (measure_type=="PercentOfTotal") {ylab = "(TP-TP0)/median(TP) (% total)"}
# 
# outfile = paste(path_outfile,"/DeltaThr",thres,"_DonorThr",donor_thres,"_",measure_type,"_",method,"_",clustering_method,"_deepSplit",deepSplit_param[run_select],"_minClusterSize",minClusterSize_param[run_select],".pdf",sep="")
# pdf(outfile)

# myCol = rev(redblue(16)) # blue-white-red
# for (i_clust in 1:n_clust) {
#     select = clust_id==i_clust
#     heatmap_panel = traj_matrix[select,]
#     rownames(heatmap_panel) = cellpop_donor[select,5]
#     colnames(heatmap_panel) = TP
#     myBreaks = quantile(abs(heatmap_panel),probs=c(seq(0.40,0.90,by=0.10),0.95,1),na.rm=T)
#     myBreaks = c(-rev(myBreaks),0,myBreaks)
#     heatmap.2(heatmap_panel, scale="none", Rowv=T, Colv=F, na.rm=T, na.color="black",
#     col = myCol, breaks = myBreaks, dendrogram = "none", cexRow=0.1, cexCol=0.6, trace="none", key=T)
# }

# we plot the trajectory clusters by median/mad
# range_y = range(c(delta_median-delta_mad,delta_median+delta_mad))
# title = paste0("Cluster size: ",clust_size[1])
# for (i_clust in 2:n_clust) {
#     title = paste0(title,", ",clust_size[i_clust])
# }
# plot(range_x,range_y,xaxt="n",type="n",xlab="",ylab=ylab,cex.lab=1.1,main=title)
# for (i_clust in 1:n_clust) {
#     offset = 2*eps*(i_clust-1)/(n_clust-1)-eps
#     points(x+offset,delta_median[i_clust,],cex=1,col=i_clust,pch=16,lty=2,type="b")
#     arrows(x+offset,delta_median[i_clust,]-delta_mad[i_clust,],x+offset,delta_median[i_clust,]+delta_mad[i_clust,],col=i_clust,length=0.025,angle=90,code=3)
# }
# axis(side = 1, at = x, labels = TP, las = 2, cex.axis=0.8)
# abline(h=0,col="black",lty=3)

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
# hm <- try(
# heatmap.2(cluster_heatmap, scale="none", Rowv=T, Colv=T, na.color="white", col=1:n_clust,
# dendrogram = "both", margins=c(5,5), cexRow=0.4, cexCol=0.35, key=F, trace="none",main="")
# )
# if(inherits(hm,"try-error")) { # to catch errors due to clustering on a very sparse matrix
#     try(
#     heatmap.2(cluster_heatmap, scale="none", Rowv=F, Colv=F, na.color="white", col=1:n_clust,
#     dendrogram = "none", margins=c(5,5), cexRow=0.4, cexCol=0.35, key=F, trace="none",main="")
#     )
# }

donor_clust_score = matrix(rep(0,n_donor*n_clust),ncol=n_clust)
rownames(donor_clust_score) = donor
colnames(donor_clust_score) = paste0("Fp0",1:n_clust)
module = matrix(rep(0,n_cellpop2*n_clust),ncol=n_clust)
rownames(module) = cellpop_labels2
colnames(module) = paste0("Fp0",1:n_clust)

df.qm = data.frame()
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
    
    # hm = Heatmap(similarity_matrix, name = "Hamming")
    # fn.fig = file.path(dn.fig, sprintf("Flow_patterns_similarity_matrix_Fp%02d", i_clust))
    # png(paste0(fn.fig, ".png"), w=500, h=400)
    # draw(hm)
    # dev.off()
    # pdf(paste0(fn.fig, ".pdf"), w=5, h=4)
    # draw(hm)
    # dev.off()
    
    KM_id = kmeans(similarity_matrix,2,nstart=50)$cluster
    
    fraction_KM_id1 = sum(cluster_heatmap2[,KM_id==1])/sum(KM_id==1)
    fraction_KM_id2 = sum(cluster_heatmap2[,KM_id==2])/sum(KM_id==2)
    if (fraction_KM_id1>fraction_KM_id2) {
        KM_id_target = 1
        KM_id_out = KM_id
        KM_fraction1 = fraction_KM_id1
        KM_fraction2 = fraction_KM_id2
    } else {
        KM_id_target = 2
        KM_id_out = 3 - KM_id
        KM_fraction1 = fraction_KM_id2
        KM_fraction2 = fraction_KM_id1
    }
    
    hm = Heatmap(similarity_matrix, name = "Hamming", split = KM_id_out)
    hm2 = Heatmap(similarity_matrix, name = "Hamming", split = factor(KM_id_out,level=c(2,1)), cluster_columns = F,
                  column_order = unlist(lapply(row_order(hm),rev))) +
          Heatmap(KM_id_out, col=c("black", "grey"), show_heatmap_legend = F, name="")
    fn.fig = file.path(dn.fig, sprintf("Flow_patterns_similarity_matrix_km2_Fp%02d", i_clust))
    png(paste0(fn.fig, ".png"), w=500, h=400)
    draw(hm2)
    dev.off()
    pdf(paste0(fn.fig, ".pdf"), w=5, h=4)
    draw(hm2)
    dev.off()
    
    kdist1 = similarity_matrix[KM_id_out==1, KM_id_out==1] %>% as.dist()
    kdist2  = similarity_matrix[KM_id_out==2, KM_id_out==2] %>% as.dist()
    
    hm3 = Heatmap(t(cluster_heatmap2), split=factor(KM_id_out,level=c(2,1)), 
                  cluster_rows = F, row_order = unlist(row_order(hm2)),
                  col=c("grey90","darkblue"),
                  show_heatmap_legend = F, row_title_side = "right",
                  name=paste0("Fp0",i_clust), column_title = paste0("Fp0",i_clust)) +
         Heatmap(KM_id_out, col=c("black", "grey"), show_heatmap_legend = F, name="")
    fn.fig = file.path(dn.fig, sprintf("Flow_patterns_pattern_heatmap_km2_Fp%02d", i_clust))
    png(paste0(fn.fig, ".png"), w=500, h=400)
    draw(hm3)
    dev.off()
    pdf(paste0(fn.fig, ".pdf"), w=5, h=4)
    draw(hm)
    dev.off()
    
    library(factoextra)
    sil = cluster::silhouette(KM_id_out, dmatrix=similarity_matrix)
    df.sil = as.data.frame(sil[,]) %>% 
      arrange(desc(cluster), desc(sil_width)) %>% 
      mutate(x=seq_along(cluster), cluster=factor(cluster), neighbor=factor(neighbor))
    df.sil = df.sil %>% 
      group_by(cluster) %>% 
      mutate(ASW = mean(sil_width, na.rm=T))
    ggplot(df.sil, aes(x, sil_width, fill=cluster)) +
      geom_col() +
      scale_color_manual(values = c("black", "grey")) +
      scale_fill_manual(values = c("black", "grey")) +
      guides(col=F, fill=F) + 
      # coord_flip() +
      # ggtitle(sprintf("Fp%02d",i_clust)) +
      xlab("") +
      ylab("Silhouette width") +
      geom_errorbar(aes(x, ymax = ASW, ymin = ASW, group = cluster),
                    size=0.5, col="red", linetype = "dashed", inherit.aes = F, width = 1) +
      # geom_hline(data=df.sil.asw, aes(yintercept=sil_width), col="red", lty=2) +
      theme_minimal() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank())
    fn.fig = file.path(dn.fig, sprintf("Flow_patterns_silhouette_horz_Fp%02d", i_clust))
    ggsave(paste0(fn.fig, ".png"), w=3,h=2)
    ggsave(paste0(fn.fig, ".pdf"), w=2,h=3)
    
    factoextra::fviz_silhouette(sil, label = FALSE, print.summary = F) +
      scale_color_manual(values = c("black", "grey")) +
      scale_fill_manual(values = c("black", "grey")) +
      guides(col=F, fill=F) + 
      coord_flip() +
      ggtitle(sprintf("Fp%02d",i_clust)) +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank())
    fn.fig = file.path(dn.fig, sprintf("Flow_patterns_silhouette_Fp%02d", i_clust))
    ggsave(paste0(fn.fig, ".png"), w=4,h=3)
    ggsave(paste0(fn.fig, ".pdf"), w=4,h=3)
    
    sil.asw = summary(sil)$clus.avg.width
    df.sil = data.frame(asw.1 = sil.asw[1], asw.2 = sil.asw[2])
    df.qm = rbind(df.qm,
                  cbind(
                    data.frame(pattern=i_clust,
                               nCP1=sum(KM_id_out==1), nCP2=sum(KM_id_out==2),
                               KM_fraction1, KM_fraction2,
                    dist.mean.1 = mean(kdist1), dist.mean.2 = mean(kdist2)),
                    df.sil
                  ))
    
    
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


df.qm = df.qm %>% mutate(pattern=fct_rev(sprintf("Fp%02d", pattern)))

p0 = df.qm %>% 
  dplyr::select(pattern, nCP1, nCP2) %>% 
  gather("stat", "value", -pattern) %>% 
  mutate(stat = sub("nCP", "cluster", stat) %>% fct_rev()) %>% 
  ggplot(aes(pattern, value, fill=stat)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("grey", "black")) +
  xlab("") +
  ylab("No. CP") +
  coord_flip() +
  guides(fill = guide_legend(reverse=T)) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = "bottom")


p1 = df.qm %>% 
  dplyr::select(pattern, matches("KM_fraction")) %>% 
  gather("stat", "value", -pattern) %>% 
  mutate(stat = sub("KM_fraction", "cluster", stat) %>% fct_rev()) %>% 
  ggplot(aes(pattern, value, fill=stat)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("grey", "black")) +
  xlab("") +
  ylab("Membership Fraction") +
  coord_flip() +
  guides(fill = guide_legend(reverse=T)) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = "bottom")

p2 = df.qm %>% 
  dplyr::select(pattern, matches("dist.mean")) %>% 
  gather("stat", "value", -pattern) %>%
  mutate(stat = sub("dist.mean.", "cluster", stat) %>% fct_rev()) %>% 
  ggplot(aes(pattern, value, fill=stat)) +
  geom_col(col=NA, position = "dodge") +
  scale_fill_manual(values = c("grey", "black")) +
  xlab("") +
  ylab("Average Humming distance") +
  coord_flip() +
  guides(fill = guide_legend(reverse=T)) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = "bottom")

p3 = df.qm %>% 
  dplyr::select(pattern, matches("asw")) %>% 
  gather("stat", "ASW", -pattern) %>% 
  mutate(stat = sub("asw.", "cluster", stat) %>% fct_rev()) %>% 
  ggplot(aes(pattern, ASW, group=stat, fill=stat)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("grey", "black")) +
  xlab("") +
  ylab("Average silhouette width") +
  coord_flip() +
  guides(fill = guide_legend(reverse=T)) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = "bottom")

p4 = sil_comb.df %>% 
  mutate(pattern=fct_rev(pattern)) %>% 
ggplot(aes(pattern, ASW, fill=pattern)) +
  geom_col() +
  scale_fill_manual(values = rev(cm.pat)) +
  xlab("") +
  ylab("ASW") +
  coord_flip() +
  guides(fill = F) +
  theme_bw()
fn.fig = file.path(dn.fig, "Flow_patterns_ASW")
ggsave(paste0(fn.fig, ".png"), plot=p4, w=3, h=4.4)  
ggsave(paste0(fn.fig, ".pdf"), plot=p4, w=3, h=4.4)  

g0 <- ggplot_gtable(ggplot_build(p0))
g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))
g3 <- ggplot_gtable(ggplot_build(p3))

l = list(g0, g1, g2, g3)

library(gridExtra)

mg = grid.arrange(arrangeGrob(g0), arrangeGrob(g1), arrangeGrob(g2), arrangeGrob(g3), nrow=1)#, width=c(0.5, 0.5))
fn.fig = file.path(dn.fig, "Flow_patterns_QMs")
ggsave(file=paste0(fn.fig, ".png"), plot=mg, w=9, h=5)
ggsave(file=paste0(fn.fig, ".pdf"), plot=mg, w=9, h=5)

fwrite(df.qm, file.path(path_outfile, "Flow_pattern_QM.txt"), sep="\t")



