# create hallmark termdf  
reactome = readRDS(here("signature_curation/reactome.rds"))

# create combined signature termdf
term2gene = list()
for (i in 1:length(reactome)) {
  
  # get genes in module 
  gene = reactome[[i]] %>% unlist(use.names = FALSE) %>% as.character()
  module = rep(names(reactome[i]), length(gene))
  df = data.frame(module, gene) %>% mutate_if(is.factor, as.character)
  
  entrez = AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, 
                                 keys = gene, 
                                 columns = c("ENTREZID", "SYMBOL"),
                                 keytype = "SYMBOL") %>% 
    select(gene = SYMBOL , ENTREZID) %>% 
    filter(ENTREZID %ni% '<NA>')
  
  term2gene[[i]] = full_join(df, entrez, by = "gene") %>% select(module, entrezid = ENTREZID)
}
# merge results list 
term_df = term2gene %>% bind_rows()
saveRDS(term_df, file = here("signature_curation/reactome_cluster_profiler_termdf.rds"))

