# was import_luminex_ip10_170426.r

dn.in = file.path(PROJECT_DIR, "DATA_PROCESSED/Luminex/")
df.lumi.long = readRDS(file.path(dn.in, "luminex_data.rds"))

df.ip10 = df.lumi.long %>% 
  filter(analyte=="IP-10")
  # mutate(label=as.numeric(sub("H5N1-","",Sample.ID))) %>%

df.ip10.score = df.ip10 %>% 
  group_by(subject, analyte) %>% 
  mutate(fc.d1 = value[Time==24] / value[Time==0 & batch==1]) %>%  
  mutate(fc.d22 = value[Time==528] / value[Time==0 & batch==2]) %>%  
  mutate(fc = max(value[Time %in% c(24,516,528)],na.rm=T)/
           median(value[!Time %in% c(24,516,528)],na.rm=T)) %>% 
  mutate(fc.mean = mean(c(value[Time==24] / value[Time==0 & batch==1], 
                   value[Time==528] / value[Time==0 & batch==2]), na.rm=T) ) %>% 
  mutate(fc.max = max(c(value[Time==24] / value[Time==0 & batch==1], 
                          value[Time==528] / value[Time==0 & batch==2]), na.rm=T) ) %>% 
  ungroup() %>% 
  dplyr::select(-Time, -time, -batch, -value) %>% 
  distinct() #%>% 
  # mutate(subject = factor(subject, levels=unique(subject)[order(fc)]))

dn.out = file.path(PROJECT_DIR, "RESULTS/Luminex/")
dir.create(dn.out, showWarnings = F)
saveRDS(df.ip10.score, file.path(dn.out, "ip10_2peak_score.rds"))
fwrite(df.ip10.score, file.path(dn.out, "ip10_2peak_score.txt"), sep="\t")
