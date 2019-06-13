# was import_luminex_ip10_170426.r

source(file.path(PROJECT_DIR, "SCRIPTS/functions/hours2days.r"))

dn = file.path(PROJECT_DIR, "DATA_ORIGINAL/Luminex/")
fn = dir(dn,"[^_]_CIR.csv")
df.lumi = data.frame()
for(f in seq_along(fn)) {
  # tmp = read.csv(file.path(dn,fn[f]), sep="\t", header=T, row.names=NULL, stringsAsFactors = F) %>% 
  tmp = fread(file.path(dn,fn[f]), sep="\t") %>% 
    mutate(batch = f)
  df.lumi = rbind(df.lumi, tmp)
}

df.lumi = df.lumi %>% 
  dplyr::rename(subject = `Sample ID`) %>% 
  arrange(Time) %>% 
  mutate(time = hours2days(Time)) %>% 
  mutate(time = factor(time, levels=unique(time)))

cidx = (!names(df.lumi) %in% c("subject", "Time", "time", "batch")) %>% which()
df.lumi.long = df.lumi %>% 
  gather("analyte", "value", cidx) %>% 
  mutate(value = as.numeric(value)) %>% 
  mutate(analyte = sub("^Hu ", "", analyte))

dn.out = file.path(PROJECT_DIR, "DATA_PROCESSED/Luminex/")
dir.create(dn.out, showWarnings = F)
saveRDS(df.lumi.long, file.path(dn.out, "luminex_data.rds"))
