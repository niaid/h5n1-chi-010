# Purpose: To find out if the difference between NEG and POS titers
# on day 42 is greater than on day 28.

library(tidyverse)
library(stringr)

adj_status <- read.table(file.path("DATA_ORIGINAL", "Clinical","clinical_info_adj.txt"), sep ="\t", header = TRUE)
titre <- read.table(file.path("DATA_ORIGINAL", "Titres","titre.txt"), sep ="\t", header = TRUE)

d28 <- titre[titre$TimePt == "D28", c("Sample.ID", "A.Indonesia")]
d28$A.Indonesia <- as.numeric(gsub("<", "", d28$A.Indonesia))

d42 <- titre[titre$TimePt == "D42", c("Sample.ID", "A.Indonesia")]
d42$A.Indonesia <- as.numeric(gsub("<", "", d42$A.Indonesia))

d28$Sample.ID <- paste0("H5N1-",str_pad(d28$Sample.ID, 3, pad = "0"))
d42$Sample.ID <- paste0("H5N1-",str_pad(d42$Sample.ID, 3, pad = "0"))

colnames(d28) <- c("Subject.ID", "A.Indonesia")
colnames(d42) <- c("Subject.ID", "A.Indonesia")

d28 <- merge(d28, adj_status)
d42 <- merge(d42, adj_status)

diff_28 <- median(d28[d28$Adjuvant == "Adj", "A.Indonesia"]) - median(d28[d28$Adjuvant == "NonAdj", "A.Indonesia"]) 
diff_42 <- median(d42[d42$Adjuvant == "Adj", "A.Indonesia"]) - median(d42[d42$Adjuvant == "NonAdj", "A.Indonesia"]) 
