source("SCRIPTS/0_initialize.r")
source("SCRIPTS/functions/today.r")
source("SCRIPTS/functions/get_score.r")
library(ComplexHeatmap)
library(circlize)

eset.gene = readRDS("DATA_PROCESSED/Emory/eset.gene.rds")

dat = exprs(eset.gene)
info = pData(eset.gene)
iord = info$day %>% unique %>% sub("Day","",.) %>% as.numeric() %>% order()
info = info %>% 
  mutate(day=factor(day, levels=info$day[iord]))

z = matrix(nrow=nrow(info),ncol=3)
colnames(z) = paste0("G0",1:3)

for(k in 1:3) {
  fn.sig = sprintf("GE_pattern_%02d_gene.sig.contr_ann.txt",k)
  pat.gene = read.table(fn.sig, sep="\t", header=T, row.names=NULL, 
                        stringsAsFactors = F, quote = "")$gene
  gi = rownames(dat) %in% pat.gene
  sum(gi)
  
  tmp = dat[gi,]
  z[,k] = get_score(tmp)
}

infz = info %>% 
  cbind(data.frame(z)) %>% 
  gather("pattern","score",matches("^G0")) %>% 
  mutate(pattern = sub("G","Gp",pattern)) %>% 
  mutate(day = sub("Day","D",day)) %>% 
  mutate(day = factor(day, levels=levels(info$day) %>% sub("Day","D",.)))

ss = c("HCAF017", "HCAF026")
ggplot(infz %>% filter(subject %in% ss), aes(day, score, col=pattern, group=pattern)) + geom_line() +
  facet_wrap(~subject, nrow=1) +
  geom_vline(xintercept=c(2,5), lty=2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        strip.background = element_blank())
ggsave(sprintf("emory/emory_profiles_wrong.subjects_3patterns_%s.png",k,today()), w=5,h=2)
  
