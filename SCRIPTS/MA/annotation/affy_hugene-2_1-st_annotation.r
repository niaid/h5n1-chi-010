file.ann = file.path(PROJECT_DIR, "DATA_ORIGINAL/Microarrays/annotation", "HuGene-2_1-st-v1.na35.hg19.transcript.csv")
ann = read_csv(file.ann, na = "---", progress = TRUE, skip = 21)

ann2 <- ann %>%
  select(one_of("transcript_cluster_id","gene_assignment")) %>%
  filter(gene_assignment != "---")

ann3 <- ann2 %>%
  mutate(gene_assignment = strsplit(gene_assignment, " /// ")) %>%
  unnest(gene_assignment) %>%
  filter(!grepl("WARNING",gene_assignment))
  
ann4 <- ann3 %>%
  separate(gene_assignment, c("accession","symbol","title","cytoband","entrez"), sep=" // ")

ann5 <- ann4 %>%
  group_by(transcript_cluster_id) %>%
  summarise(n=n_distinct(symbol), symbol=first(symbol), title=first(title), cytobind=first(cytoband), entrez=first(entrez)) %>%
  filter(n==1) %>%
  select(-n)

fn.out = file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/annotation/affy_hugene-2_1-st_unique_all.txt")
write.table(ann5, file=fn.out, sep="\t", col.names = TRUE, row.names=FALSE)

ann6 <- ann5 %>%
  select(one_of(c("transcript_cluster_id","symbol", "entrez"))) %>%
  rename(probeset=transcript_cluster_id, hgnc_symbol=symbol, llid=entrez) %>%
  mutate(llid=sub("---","NA", llid))

fn.out = file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/annotation/affy_hugene-2_1-st_unique.txt")
write.table(ann6, file=fn.out, sep="\t", col.names = TRUE, row.names=FALSE, quote = FALSE)
