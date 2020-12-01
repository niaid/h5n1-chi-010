source("SCRIPTS/0_initialize.r")

# Load Demographics Data
df_demo <- read.table(file.path(PROJECT_DIR, "DATA_ORIGINAL", "Emory", "emory_demo.txt"), sep = "\t", header = T)
df_demo$Sex <- gsub( "female", "F", df_demo$Sex)
df_demo$Sex <- gsub( "male", "M", df_demo$Sex)
df_demo$SubjectID <- as.character(df_demo$SubjectID)
colnames(df_demo)[4:5] <- c("sex", "age")

# Load Titre Data
df_mn_emory <- read.table(file.path(PROJECT_DIR, "DATA_PROCESSED", "Emory", "Emory_Adj_MN_D0_21_42.txt"), sep = "\t", header = T)
df_mn_emory[2:4] <- log2(df_mn_emory[2:4]/10)
colnames(df_mn_emory)[1] <- "SubjectID"
df_mn_emory$SubjectID <- as.character(df_mn_emory$SubjectID)

# SOMASCAN BASELINE DATA

fn.rfu = file.path(PROJECT_DIR, "DATA_PROCESSED", "Emory", "somascan","Final_RFU_MedNorm_Cal.txt")
fn.samp = file.path(PROJECT_DIR, "DATA_PROCESSED","Emory", "somascan", "Final_Samples_MedNorm_Cal.txt")
fn.soma = file.path(PROJECT_DIR, "DATA_PROCESSED", "Emory", "somascan", "Final_Somamers_MedNorm_Cal.txt")

df.soma = fread(fn.rfu, header = F, data.table = F)
df.samp = fread(fn.samp, data.table = F)
df.ann = fread(fn.soma, data.table = F)

names(df.soma) = df.ann$Target
df.soma$SubjectID = df.samp$SampleId

# COMBINE ALL TOGETHER
DF <- inner_join(df_demo[, c("SubjectID", "sex", "age")], df_mn_emory, by ="SubjectID") %>%
        left_join(df.soma, by = "SubjectID")

fn.out = file.path(PROJECT_DIR, "DATA_PROCESSED/eNet/eNet_InputData_r9.txt")
fwrite(DF, fn.out, sep="\t")



