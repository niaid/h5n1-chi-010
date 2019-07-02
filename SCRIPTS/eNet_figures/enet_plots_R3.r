source("SCRIPTS/0_initialize.r")
library(Matrix)

dn.enet = file.path(PROJECT_DIR, "RESULTS/eNet")
run.id = "R3"
run.ver = "v4"
fn.run = file.path(dn.enet, run.ver, glue::glue("{run.id}_{run.ver}.Robj"))

load(fn.run, verbose = T)

a = 0.3
ai = which(round(result$alpha, digits = 2) == a)
albl = paste0("a",a)

feature = result$feature

freq.mean = result$feature_freq_mean[,ai]
freq.sd = result$feature_freq_sd[,ai]
freq.pval = result$feature_freq_model_vs_null_pval[,ai]
coef.mean = result$feature_coef_wmean[,ai]
coef.sd = result$feature_coef_wsd[,ai]
coef.pval = result$feature_coef_model_vs_null_pval[,ai]
DF.model = data.frame(feature, freq.mean, freq.sd, freq.pval, coef.mean, coef.sd, coef.pval, model="model")

freq.mean = result$null_feature_freq_mean[,ai]
freq.sd = result$null_feature_freq_sd[,ai]
coef.mean = result$null_feature_coef_wmean[,ai]
coef.sd = result$null_feature_coef_wsd[,ai]
DF.null = data.frame(feature, freq.mean, freq.sd, freq.pval, coef.mean, coef.sd, coef.pval, model="null")

DF.stat = rbind(DF.model, DF.null) %>% 
  mutate(model = factor(model, levels=c("null","model"))) %>% 
  mutate(feature = factor(feature, levels=sort(unique(feature), decreasing=T)))

dn.fig = file.path(PROJECT_DIR, "FIGURES/eNet")
dir.create(dn.fig, showWarnings = F)

ggplot(DF.stat, aes(coef.mean, freq.mean, label=feature)) +
  geom_text(vjust=0) + facet_wrap(~model, nrow=1) +
  ggtitle(paste0(run.id, ", alpha = ", a))
fn.fig = file.path(dn.fig, sprintf("enet_%s_%s_%s_coef.mean_vs_freq_with_null.png", run.id, run.ver, albl))
ggsave(fn.fig, w=8, h=4)

ggplot(DF.model, aes(coef.mean, coef.pval, label=feature)) +
  geom_text(vjust=0) +
  geom_hline(yintercept = 0.2) +
  scale_y_log10() +
  ggtitle(paste0(run.id, ", alpha = ", a))
fn.fig = file.path(dn.fig, sprintf("enet_%s_%s_%s_coef.mean_vs_coef.pval.png", run.id, run.ver, albl))
ggsave(fn.fig, w=4, h=4)

ggplot(DF.model, aes(freq.pval, coef.pval, label=feature)) +
  geom_text(vjust=0) +
  geom_hline(yintercept = 0.2) +
  scale_y_log10() +
  ggtitle(paste0(run.id, ", alpha = ", a))
fn.fig = file.path(dn.fig, sprintf("enet_%s_%s_%s_coef.pval_vs_freq.pval.png", run.id, run.ver, albl))
ggsave(fn.fig, w=4, h=4)


# scatter plot of predictin performance
DF.plot = data.frame(
  subject = rownames(result$predictor),
  measured = result$response,
  pred_mean = result$predicted_values[[ai]][,1],
  pred_sd = result$predicted_values[[ai]][,2]
)
p0 = ggplot(DF.plot, aes(measured, pred_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin=pred_mean-pred_sd, ymax=pred_mean+pred_sd), width=0.1) +
  geom_smooth(method="lm", alpha=0.2, col="red", lty=2) +
  geom_text(aes(label=sub("s","",subject)), hjust=-0.3, size=3) +
  xlab("MN titer, measured at day 28") +
  ylab("MN titer, out-of-bag predicted") +
  theme_bw()

# plots of selected features
# features.in = c("Gp01","Gp02","Gp03","Fp01","Fp05")
features.in = NULL

if(is.null(features.in)){
  fi = DF.model$freq.mean > max(DF.null$freq.mean) & 
    abs(DF.model$coef.mean) > max(abs(DF.null$coef.mean)) & DF.model$coef.pval < 0.2
} else {
  fi = DF.stat$feature %in% features.in
}
features.in = feature[fi]

p1 = ggplot(DF.stat %>% filter(feature %in% features.in), aes(feature, freq.mean*100, fill=model)) +
  geom_col(position="dodge") +
  geom_errorbar(aes(ymin=pmax(0, (freq.mean-freq.sd)*100), ymax=pmin(100,(freq.mean+freq.sd)*100)),
    position=position_dodge(width = 0.9), width=0.25) +
  scale_fill_manual(values = c("grey60","red")) +
  coord_flip() +
  xlab("Predictor") +
  ylab("% selected") +
  guides(fill=F) +
  theme_classic() + theme(axis.line.y = element_blank())
p2 = ggplot(DF.stat %>% filter(feature %in% features.in), aes(feature, coef.mean, fill=model)) +
  geom_col(position="dodge") +
  geom_errorbar(aes(ymin=coef.mean-coef.sd, ymax=coef.mean+coef.sd),
                position=position_dodge(width = 0.9), width=0.25) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("grey60","red")) +
  coord_flip(ylim=c(-1.1,1.1)) +
  xlab("") +
  ylab("Coefficient mean") +
  guides(fill=F) +
  theme_classic() + theme(axis.line.y = element_blank())
g0 <- ggplot_gtable(ggplot_build(p0))
g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))
l = list(g0, g1, g2)

library(gridExtra)
mg = grid.arrange(arrangeGrob(g0), arrangeGrob(g1,g2, nrow=1), nrow=1)#, width=c(0.5, 0.5))
fn.fig = file.path(dn.fig, sprintf("eNet_%s_%s_%s_selected", run.id, run.ver, albl))
ggsave(file=paste0(fn.fig, ".png"), plot=mg, w=9, h=3)
ggsave(file=paste0(fn.fig, ".pdf"), plot=mg, w=9, h=3)
