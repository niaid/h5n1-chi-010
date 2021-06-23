# 
suppressMessages(library(here)) 
# load gene_pattern_genes
gpath = list.files(path = here("signature_curation/GE_pattern_genes/"), full.names = TRUE)
gpnames = list.files(path = here("signature_curation/GE_pattern_genes/"))
gpnames = str_sub(gpnames, 1,13)
gp = lapply(gpath, function(x){read_delim(x, delim = "\t") %$% gene })
names(gp) = gpnames 

# save pattern gene list 
saveRDS(gp, file = here("signature_curation/H5_GenePatterns.rds"))
