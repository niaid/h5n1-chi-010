library(eNetXplorer)
library(tidyverse)


# MAIN

ver <- "pear"
results_dir <- file.path("RESULTS", "rspo3", ver)

# Read eNetXplorer training results.
load(file.path("RESULTS", "rspo3", "new_enet_pearson.Robj"))
fit_summary <- summary(eNet)
wd <- getwd()

# Export data.
export(eNet, dest_dir=results_dir, dest_dir_create=TRUE, delim="tsv", input.data=TRUE, summary.data=TRUE, output.data=TRUE)
setwd(wd)

# Select genes at a particular alpha value. 
a = 0.2
ai = which.min(abs(eNet$alpha - a))
albl = paste0("a",a)

feature = eNet$feature

freq.mean = eNet$feature_freq_mean[,ai]
freq.sd = eNet$feature_freq_sd[,ai]
freq.pval = eNet$feature_freq_model_vs_null_pval[,ai]
coef.mean = eNet$feature_coef_wmean[,ai]
coef.sd = eNet$feature_coef_wsd[,ai]
coef.pval = eNet$feature_coef_model_vs_null_pval[,ai]
DF.model = data.frame(feature, freq.mean, freq.sd, freq.pval, coef.mean, coef.sd, coef.pval, model="model")

freq.mean = eNet$null_feature_freq_mean[,ai]
freq.sd = eNet$null_feature_freq_sd[,ai]
coef.mean = eNet$null_feature_coef_wmean[,ai]
coef.sd = eNet$null_feature_coef_wsd[,ai]
DF.null = data.frame(feature, freq.mean, freq.sd, freq.pval, coef.mean, coef.sd, coef.pval, model="null")

DF.stat = rbind(DF.model, DF.null) %>% 
          mutate(model = factor(model, levels=c("null","model"))) %>% 
            mutate(feature = factor(feature, levels=sort(unique(feature), decreasing=T)))

features.in = NULL
if(is.null(features.in)){
  fi = DF.model$freq.mean > max(DF.null$freq.mean, na.rm=T) & 
    abs(DF.model$coef.mean) > max(abs(DF.null$coef.mean[!is.na(DF.model$coef.mean)]), na.rm=T) & DF.model$coef.pval < 0.2
} else {
  fi = DF.stat$feature %in% features.in
}
features.in = feature[fi]

saveRDS(features.in, file.path(results_dir, paste(albl, "_", ver, ".rds", sep="")))

stop()

# Plots
# Summary PDF.
summaryPDF(eNet, dest_dir=file.path("FIGURES", "rspo3"), dest_file=paste("summary_", ver ,".pdf", sep =""))

# Model performance.
png(filename = file.path("FIGURES","rspo3", paste("model_performance_", ver, ".png", sep ="")), width = 6, height = 4, units = "in", pointsize = 12, res = 300)
plot(eNet, plot.type="summary")
dev.off()

# Top features.
png(filename = file.path("FIGURES","rspo3", paste("top_features_", ver, ".png", sep="")), width = 6, height = 4, units = "in", pointsize = 12, res = 300)
plot(eNet, alpha.index = which.max(eNet$model_QF_est), plot.type = "featureCaterpillar",stat = c("coef"))
dev.off()

# Selectd features coefficients.
png(filename = file.path("FIGURES","rspo3", paste("feature_coeff_", ver, ".png", sep="")), width = 6, height = 4, units = "in", res = 300, pointsize = 8)
plot(eNet, alpha.index = which.max(eNet$model_QF_est), plot.type = "featureHeatmap", stat = c("coef"), notecex = 1.5)
dev.off()


