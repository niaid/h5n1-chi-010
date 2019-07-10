library(rlang)
library(lazyeval)
source("SCRIPTS/0_initialize.r")
source(file.path(PROJECT_DIR, "SCRIPTS/functions/color_functions.r"))

pop10 = "CXCR3+"
pop15.list = c("PD-1+ICOS+ CXCR3+ Tfh", "PD-1+ICOS+ CXCR3- Tfh")
tm.sel = list(c("d0", "d1", "d7"), c("d21", "d22", "d28"))

dn.fig = file.path(PROJECT_DIR, "FIGURES/Tfh")
dir.create(dn.fig, showWarnings = F)

### read and process 10-color flow
fn.tfh = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_10c/CXCR3_freq.txt")
tfh = fread(fn.tfh, data.table=F, header = T) %>% 
  rename(subject=`PATIENT ID`, time=`SAMPLE ID`) %>% 
  rename(T_CD3CD4 = `4.2`, T_CD3 = `1.2`, Tfh = `24`) %>% 
  filter(subject!="")

tfh = tfh %>% 
  mutate(time.point = case_when(
    time=="d8" ~ "d7",
    time=="d21_2h" ~ "d21_4h",
    time=="d29" ~ "d28",
    time=="d31" ~ "d28",
    time=="d43" ~ "d42",
    time=="d0_0h" ~ "d0",
    time=="d21_0h" ~ "d21",
    TRUE ~ time
  ))

### read and process 15-color flow
fn.tfh15 = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_15c/H5N1-all-T4 021518.txt")
tfh15 = fread(fn.tfh15, data.table=F) %>% 
  filter(!Sample %in% c("Mean","StdDev"))
tfh15$Sample = tfh15$`$FIL`
tfh15 = tfh15 %>% select(Sample, matches("^Lym"))

names(tfh15)[-1] = c("CXCR3+ Tfh","PD-1+ICOS+ CXCR3+ Tfh", 
                     "Tfh2","PD-1+ICOS+ Tfh2", 
                     "CXCR3- Tfh","PD-1+ICOS+ CXCR3- Tfh")
tfh15.levels = names(tfh15)[-1]

tfh15 = tfh15 %>% 
  filter(!grepl("H5N1-CHI", Sample)) %>% 
  filter(!Sample %in% c("Mean","StdDev"))
tfh15 = tfh15 %>% 
  mutate(subject = str_extract(.$Sample, "H5N1-\\d{3}")) %>% 
  mutate(time.point = str_extract(.$Sample, "day \\d{1,2}") %>% sub("ay ","", .))

### demographics and adjuvant status
fn.clin = file.path(PROJECT_DIR, "DATA_ORIGINAL/Clinical/clinical_info_adj.txt")
df.clin = fread(fn.clin) %>% 
  rename(subject=`Subject ID`)


### read titer d28
fn.mn = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Titres", "titre.txt")
df.mn = read.table(fn.mn, sep="\t", header=T, row.names=NULL, stringsAsFactors=T)
df.mn = df.mn %>% 
  mutate_at(-(1:2), function(x) sub("<20","10",x) %>% as.numeric()) %>% 
  mutate(day = as.numeric(sub("[dD]","",TimePt))) %>% 
  mutate(day = ifelse(day %in% c(29,31), 28, day) %>% factor()) %>% 
  mutate(subject = sprintf("H5N1-%03d",Sample.ID)) %>% 
  filter(day=="28") %>% 
  select(subject, titer.d28=A.Indonesia) %>% 
  mutate(titer.d28.log2 = log2(titer.d28))


for (tm in tm.sel) {
    
  df.tfh = tfh %>% 
    filter(time.point %in% tm[c(1,2)]) %>% 
    mutate(time.point = factor(time.point, levels = unique(time.point))) %>% 
    select(subject, time.point, one_of(pop10)) %>% 
    gather("population","pop", pop10)
  
  df.tfh = df.tfh %>% 
    group_by(subject, population) %>% 
    mutate(non0base = any(time.point==tm[1] & pop>0)) %>% 
    ungroup() %>% 
    filter(subject != "H5N1-044") %>%
    filter(non0base) %>% 
    select(-non0base)
  
  df.tfh.change = df.tfh %>% 
    group_by(subject,time.point,population) %>% mutate(n=n()) %>% ungroup() %>% filter(n==1) %>% select(-n) %>% #remove duplicates
    group_by(subject, population) %>% 
    mutate(FC = pop / pop[time.point==tm[1]]) %>% 
    select(-pop) %>% 
    filter(time.point != tm[1]) %>% 
    gather(type, value, "FC") %>% 
    unite(population, population, time.point, type, sep=".") %>% 
    spread(population, value)
  
  
  for(pop15 in pop15.list) {

    df.tfh15 = tfh15 %>% 
      filter(time.point %in% tm[c(1,3)]) %>% 
      mutate(time.point = factor(time.point, levels = unique(time.point))) %>% 
      select(subject, time.point, one_of(pop15)) %>%
      gather("population","pop", one_of(pop15))
    
    df.tfh15 = df.tfh15 %>% 
      group_by(subject, population) %>% 
      mutate(non0base = any(time.point==tm[1] & pop>0)) %>% 
      ungroup() %>% 
      filter(subject != "H5N1-044") %>%
      filter(non0base) %>% 
      select(-non0base)
    
    df.tfh15.change = df.tfh15 %>% 
      group_by(subject, population) %>% 
      mutate(FC = pop / pop[time.point==tm[1]]) %>% 
      select(-pop) %>% 
      filter(time.point != tm[1]) %>% 
      gather(type, value, "FC") %>% 
      unite(population, population, time.point, type, sep=".") %>% 
      spread(population, value)
    
    
    df.comb = inner_join(df.tfh.change, df.tfh15.change, by=c("subject")) %>% 
      inner_join(df.mn, by="subject") %>% 
      inner_join(df.clin %>% select(subject, Adjuvant), by="subject")
    
    DF = rbind(mutate(df.comb, Adjuvant="All"), df.comb) %>% 
      mutate(Adjuvant = factor(Adjuvant, levels = c("All","Adj","NonAdj")))
    
    cols = names(DF)

    v1 = 2
    for(v2 in c(3,5)) {
    col1 = cols[v1]
    col2 = cols[v2]
    
    df.cc = DF %>% 
      group_by(Adjuvant) %>% 
      summarise_(cor = interp(~cor.test(var1, var2, method="spearman", exact=F)$estimate, 
                              var1=as.name(col1), var2=as.name(col2)),
                 p = interp(~cor.test(var1, var2, method="spearman", exact=F)$p.value, 
                            var1=as.name(col1), var2=as.name(col2))) %>% 
      mutate(sign = case_when(
        p<=0.001 ~ "***",
        p<=0.01 ~ "**",
        p<=0.05 ~ "*",
        TRUE ~ ""
      )) %>% 
      tibble::add_column(var1=col1, var2=col2, .before = 1)
    
    xrange = range(pull(DF, col1), na.rm = T, finite = T)
    yrange = range(pull(DF, col2), na.rm = T, finite = T)
    xlbl = xrange[1] + diff(xrange)*0.7
    ylbl = yrange[1] + diff(yrange)*seq(0.8,1,length.out=3) %>% rev()
    df.cc = df.cc %>% 
      mutate(y=Inf, x=xlbl)
    
    
    col1 = glue::glue("`{cols[v1]}`")
    col2 = glue::glue("`{cols[v2]}`")
    p = ggplot(DF %>% filter(Adjuvant!="All"), aes_string(col1, col2, col="Adjuvant")) +
      geom_point(size=2) +
      geom_smooth(col = "black", method="lm", fill=NA, show.legend = F) +
      geom_text(data=df.cc %>% filter(Adjuvant=="All"), 
                aes(x=x,y=y,label=glue::glue("p = {format(p, digit=2)}{sign}")), hjust=0, vjust=1.2, 
                inherit.aes = F,
                show.legend = F) +
      theme_bw()
    fn.fig = glue::glue("{dn.fig}/{cols[v1]}_vs_{cols[v2]}_correlation.png")
    ggsave(fn.fig, w=5, h=3)
    
  }
  }
}
