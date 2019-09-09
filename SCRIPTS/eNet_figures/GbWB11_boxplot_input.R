# Initialize
source("SCRIPTS/0_initialize.r")

# DEMOGRAPHICS

fn.demo = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Clinical", "clinical_info.txt")
df.demo = fread(fn.demo) %>% 
  mutate(donor_short = sub("H5N1-0", "s", `Subject ID`)) %>% 
  dplyr::select(donor = `Subject ID`, donor_short, sex= Gender, age = Age)

# ADJUAVNT PREDICTION

fn.pred = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Clinical", "clinical_info_adj.txt")
df.pred = fread(fn.pred)
df.pred$subject = sub("H5N1-0", "s", df.pred$`Subject ID`) 
si = match(df.demo$donor_short, df.pred$subject)
df.demo$blinded_flag = df.pred$Adjuvant[si]
df.demo = df.demo %>% mutate (blinded_flag = sub("NonAdj", "NEG", blinded_flag)) %>%
        mutate (blinded_flag = sub("Adj", "POS", blinded_flag))

# MN TITERS

fn.mn = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Titres", "titre.txt")
df.mn = fread(fn.mn) %>% 
  mutate_at(-(1:2), function(x) sub("<20","10",x) %>% as.numeric()) %>% 
  # mutate(day = as.numeric(sub("[dD]","",TimePt))) %>% 
  mutate(day = sub("D","d",TimePt)) %>% 
  mutate(day = ifelse(day %in% c("d29","d31"), "d28", day)) %>% 
  mutate(donor = sprintf("H5N1-%03d",`Sample ID`)) %>% 
  mutate(titer = log2(`A/Indonesia`/10)) %>% 
  # mutate(day = paste0("titer_I_MN_", day)) %>% 
  dplyr::select(donor, day, titer) %>% 
  spread(day, titer)
names(df.mn)[-1] = sub("^d", "titer_I_MN_d", names(df.mn)[-1])

# fwrite(df.mn, "mn.txt", sep="\t")
# H5N1-019 at d42 has different MN titer in mn (1280) and ab (320) data


# FLOW BASELINE FREQUENCIES

fn.gates = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_10c", "GateList_10c.txt")
df.gates = fread(fn.gates, header = T)

fn.b = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_10c", "H5N1_Bcell-H5N1_Bcell_freq.V1.1")
df.b = fread(fn.b, header = T) %>% 
  dplyr::rename(donor = V1, time = Timepoint) %>% 
  dplyr::select(-LYM) %>% 
  dplyr::filter(time=="d0_0h") %>% 
  gather("Id", "value", -c(donor, time))

fn.t = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_10c", "H5N1_Tcell-H5N1-Tcell_freq.V1.1")
df.t = fread(fn.t, header = T) %>% 
  dplyr::rename(donor = V1, time = Timepoint) %>% 
  dplyr::select(-LYM) %>% 
  dplyr::filter(time=="d0_0h") %>% 
  gather("Id", "value", -c(donor, time))

fn.th = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_10c", "H5N1_Thelper-Thelper_freq.V1.1")
df.th = fread(fn.th, header = T) %>% 
  dplyr::rename(donor = V1, time = Timepoint) %>% 
  dplyr::select(-LYM) %>% 
  dplyr::filter(time=="d0_0h") %>% 
  gather("Id", "value", -c(donor, time))

fn.tr = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_10c", "H5N1_Treg-H5N1-Treg_freq.V1.1")
df.tr = fread(fn.tr, header = T) %>% 
  dplyr::rename(donor = V1, time = Timepoint) %>% 
  dplyr::select(-LYM) %>% 
  dplyr::filter(time=="d0_0h") %>% 
  gather("Id", "value", -c(donor, time))

fn.dc = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_10c", "H5N1_DCMonNK-DCMonNK_freq.V1.1")
df.dc = fread(fn.dc, header = T) %>% 
  dplyr::rename(donor = V1, time = Timepoint) %>% 
  dplyr::select(-`CD45+`, -MNC, -`58`) %>% 
  dplyr::filter(time=="d0_0h") %>% 
  gather("Id", "value", -c(donor, time))

df.flow = bind_rows(df.b, df.t, df.th, df.tr, df.dc) %>% 
  mutate(Id = as.numeric(Id))
df.flow = df.flow %>% 
  inner_join(df.gates %>% dplyr::select(Code, Id), by="Id")
df.fb = df.flow %>% 
  dplyr::select(-time, -Id) %>% 
  spread(Code, value)

# GE BASELINE WGCNA modules scores 

fn.gb = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/baseline/GE.pbmc_d0_WGCNA_ME_scores.txt")
df.gb = fread(fn.gb) %>% 
  dplyr::rename(donor=subject)

fn.gbw = file.path(PROJECT_DIR, "RESULTS/Microarrays/PAXgene/baseline/GE.pax_d0_WGCNA_ME_scores.txt")
df.gbw = fread(fn.gbw) %>% 
  dplyr::rename(donor=subject)

# COMBINE ALL TOGETHER

DF = inner_join(df.demo, df.mn, by="donor") %>% 
  left_join(df.fb, by="donor") %>% 
  left_join(df.gb, by="donor") %>% 
  left_join(df.gbw, by="donor") %>% 
  dplyr::filter(blinded_flag == "NEG")
fn.out = file.path(PROJECT_DIR, "DATA_PROCESSED/eNet/eNet_InputData_r6.txt")
fwrite(DF, fn.out, sep="\t")

# fn.old = file.path(PROJECT_DIR, "DATA_PROCESSED/eNet", "eNet_InputData.txt")
# DF.old = fread(fn.old)
# DF2 = DF.old[,1:ncol(DF)]
# 
# all.equal(DF, DF2)
# which(DF$blinded_flag != DF2$blinded_flag)

