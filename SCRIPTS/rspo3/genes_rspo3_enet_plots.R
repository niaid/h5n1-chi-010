library(eNetXplorer)
library(tidyverse)


# MAIN

# Read eNetXplorer training results.
ver <- "v4"
fit <- readRDS(file.path("RESULTS", "rspo3", paste("genes_rspo3_enet_fit_", ver , ".rds", sep ="")))
fit_summary <- summary(fit)
wd <- getwd()

# Plots
# export(fit, file.path("RESULTS", "rspo3"), delim="tab", input.data=F, summary.data=T, output.data=T)
# setwd(wd)

# Summary PDF.
summaryPDF(x=fit, path=file.path("FIGURES", "rspo3"), file=paste("summary_", ver ,".pdf", sep =""))

# Model performance.
png(filename = file.path("FIGURES","rspo3", paste("model_performance_", ver, ".png", sep ="")), width = 6, height = 4, units = "in", pointsize = 12, res = 300)
plot(fit, plot.type="summary")
dev.off()

# Top features.
png(filename = file.path("FIGURES","rspo3", paste("top_features_", ver, ".png", sep="")), width = 6, height = 4, units = "in", pointsize = 12, res = 300)
plot(fit, alpha.index = which.max(fit$model_QF_est), plot.type = "featureCaterpillar",stat = c("coef"))
dev.off()

# Selectd features coefficients.
png(filename = file.path("FIGURES","rspo3", paste("feature_coeff_", ver, ".png", sep="")), width = 6, height = 4, units = "in", res = 300, pointsize = 8)
plot(fit, alpha.index = which.max(fit$model_QF_est), plot.type = "featureHeatmap", stat = c("coef"), notecex = 1.5)
dev.off()


