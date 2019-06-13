# DEMOGRAPHICS

fn.demo = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Clinical", "clinical_info.txt")
df.demo = fread(fn.demo) %>% 
  mutate(donor_short = sub("H5N1-0", "s", `Subject ID`)) %>% 
  dplyr::select(donor = `Subject ID`, donor_short, sex= Gender, age = Age)

# ADJUAVNT PREDICTION

fn.pred = file.path(PROJECT_DIR, "RESULTS/Adjuvant_prediction/adjuvant_predicted_subjects.txt")
df.pred = fread(fn.pred)
si = match(df.demo$donor_short, df.pred$subject)
df.demo$blinded_flag = df.pred$predict[si]

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

# FLOW PATTERN SCORES

fn.fp = file.path(PROJECT_DIR, "RESULTS/Flow_10c", "DeltaThr0.5_DonorThr0.3_PercentOfParent_CORpearson_cutreeHybrid_deepSplit0_minClusterSize31_DonorScores.txt")
df.fp = fread(fn.fp) %>% 
  dplyr::rename(donor_short = Donor) %>% 
  arrange(donor_short)

# GE PATTERN SCORES

fn.gp = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery/", "GE_pattern_scores_42subj.txt")
df.gp = fread(fn.gp) %>% 
  dplyr::rename(donor=subject)

# COMBINE ALL TOGETHER

DF = inner_join(df.demo, df.mn, by="donor") %>% 
  # left_join(df.ab, by="donor") %>% 
  left_join(df.fp, by="donor_short") %>% 
  left_join(df.gp, by="donor") %>% 
  filter(blinded_flag == "NEG")
fn.out = file.path(PROJECT_DIR, "DATA_PROCESSED/eNet/eNet_InputData_r2.txt")
fwrite(DF, fn.out, sep="\t")

# fn.old = file.path(PROJECT_DIR, "DATA_PROCESSED/eNet", "eNet_InputData.txt")
# DF.old = fread(fn.old)
# DF2 = DF.old[,1:ncol(DF)]
# 
# all.equal(DF, DF2)
# which(DF$blinded_flag != DF2$blinded_flag)

