
# was d0_wgcna_170319.r
# Initialize.
source("SCRIPTS/0_initialize.r")

dn.in = file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC/baseline")
dat0.in = readRDS(file.path(dn.in, "dat0.in_isv.0.7_8144g.rds"))
info0.in = readRDS(file.path(dn.in, "info0.in.rds")) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("subject.id")

dn.fig = file.path(PROJECT_DIR, "FIGURES/Baseline_PBMC")
dir.create(dn.fig, showWarnings = F)

library(WGCNA)
enableWGCNAThreads()
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(dat0.in, powerVector = powers, verbose = 5)

# Plot the results:

png(file.path(dn.fig, "wgcna_select_power.png"), width=600, height=400)
# sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red") # this line corresponds to using an R^2 cut-off of h
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

power.sel = 6
net = blockwiseModules(dat0.in, power = power.sel,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = sprintf("d0TOM.%d",power.sel), 
                       verbose = 3) 

png(file.path(dn.fig, sprintf("dat0.in-networkConstruction-auto_HCLUST_power.%d.png", power.sel)),
         width=1000, height=600)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05) 
dev.off()

table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

dn.out = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/baseline")
dir.create(dn.out, showWarnings = F)
save(net, MEs, moduleLabels, moduleColors, geneTree, 
     file = file.path(dn.out, sprintf("dat0.in-networkConstruction-auto_power.%d.RData", power.sel)))
