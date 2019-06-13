# DEMOGRAPHICS

fn.demo = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Clinical", "clinical_info_adj.txt")
df.demo = fread(fn.demo) %>% 
  mutate(donor_short = sub("H5N1-0", "s", `Subject ID`)) %>% 
  mutate(Adjuvant = toupper(Adjuvant)) %>% 
  dplyr::select(donor = `Subject ID`, donor_short, sex= Gender, age = Age, unblinded_flag = Adjuvant)

# ADJUAVNT PREDICTION

# fn.pred = file.path(PROJECT_DIR, "RESULTS/Adjuvant_prediction/adjuvant_predicted_subjects.txt")
# df.pred = fread(fn.pred)
# si = match(df.demo$donor_short, df.pred$subject)
# df.demo$blinded_flag = df.pred$predict[si]

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


# SOMASCAN BASELINE DATA

fn.rfu = file.path(PROJECT_DIR, "DATA_PROCESSED/SOMAscan", "Final_RFU.txt")
fn.samp = file.path(PROJECT_DIR, "DATA_PROCESSED/SOMAscan", "Final_Samples.txt")
fn.soma = file.path(PROJECT_DIR, "DATA_PROCESSED/SOMAscan", "Final_Somamers.txt")

df.soma = fread(fn.rfu, header = F, data.table = F)
df.samp = fread(fn.samp, data.table = F)
df.ann = fread(fn.soma, data.table = F)

names(df.soma) = df.ann$Target
df.soma$donor_short = df.samp$SampleId


# COMBINE ALL TOGETHER

DF = inner_join(df.demo, df.mn, by="donor") %>% 
  left_join(df.soma, by="donor_short") %>% 
  dplyr::filter(unblinded_flag == "ADJ")
fn.out = file.path(PROJECT_DIR, "DATA_PROCESSED/eNet/eNet_InputData_r8.txt")
fwrite(DF, fn.out, sep="\t")

# CHECK
# fn.old = file.path(PROJECT_DIR, "DATA_PROCESSED/eNet", "eNet_InputData_NEW.txt")
# DF.old = fread(fn.old, data.table = F)
# names(DF) = make.names(names(DF)) %>% gsub("\\.\\.","\\.",.)
# i1 = names(DF) %in% names(DF.old)
# sum(i1)
# names(DF)[!i1]
# i2 = names(DF.old) %in% names(DF)
# i2 = match(names(DF), names(DF.old))
# 
# DF2 = DF.old[,i2]
# 
# all.equal(DF, DF2)
# which(DF$blinded_flag != DF2$blinded_flag)
# k = 11
# all.equal(DF[,k], DF2[,k])
# names(DF)[k]
# names(DF2)[k]
