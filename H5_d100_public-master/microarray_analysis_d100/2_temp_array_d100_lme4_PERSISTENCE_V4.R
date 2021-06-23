#  bulk fold change delta PERSISTENCE SIGNATURES 
suppressMessages(library(here))
suppressMessages(library(tidyverse))
# suppressMessages(library(variancePartition))
# suppressMessages(library(magrittr))
'%ni%' = Negate('%in%')
# suppressMessages(library(scglmmr))

# save path 
# datapath = here("microarray_analysis_d100/generated_data_V4/"); dir.create(datapath)
#  figpath = here("microarray_analysis_d100/figures_V4/"); dir.create(figpath)

# pattern gene object createed in dir signature_curation 
gp = readRDS(here("signature_curation/H5_GenePatterns.rds"))
pattern_genes = unlist(gp, use.names = FALSE)
# persist_pattern = gp[c(4,6:8, 10:12)] 
# expression dataframe 
edf = readRDS(here("microarray_analysis_d100/generated_data_V4/array_d0_d100_patterngenes_andmetadata_dataframe.rds"))
met = read_delim(file = "data/metadata/clinical_info_adj.txt",delim = "\t")


df_t_out <- data.frame()
for(pp in 1:14){
        persist_pattern = gp[pp] 
        persist_genes = unlist(x = persist_pattern, use.names = FALSE) %>% unique()

# gene indexes for tidy 
        index1 = persist_genes[1]; index2 = persist_genes[length(persist_genes)]

# tidy 
        eddf = edf %>% 
          # remove subject without day 100 measurement 
          filter(subject.id %ni% c("H5N1-005", "H5N1-003", "H5N1-006")) %>% 
          select(subject.id, time.point, all_of(persist_genes)) %>% 
          gather(gene, expression, all_of(index1):all_of(index2)) %>%  
          # calculate average expression of the multiple baseline timepoints 
          group_by(subject.id, time.point, gene) %>% 
          summarize(expression = mean(expression, na.rm = TRUE)) %>% 
          arrange(gene, subject.id) %>% 
          ungroup() %>% 
          spread(gene, expression) %>% 
          mutate(sample = paste(subject.id, time.point, sep = "_")) %>% 
          select(sample,  everything())

# gene expressiondata matrix 
        exprs_mat = eddf[ ,c(1, 4:ncol(eddf))] %>%
          column_to_rownames("sample")

        genes <- colnames(exprs_mat)
        sub_time <- rownames(exprs_mat)
        d0_sub <- sub_time[str_detect(sub_time, "d000")]
        d100_sub <- sub_time[str_detect(sub_time, "d100")]

        for(g in genes){
                d0_dat <- exprs_mat[d0_sub, g]
                d100_dat <- exprs_mat[d100_sub, g]
                t_res <- t.test(d0_dat, d100_dat, paired = TRUE)
                temp_df <- data.frame("gene" = g, "pattern" = pp, "t_p" = t_res$p.value)
                df_t_out <- rbind(df_t_out, temp_df)
        }




}
df_t_adj <- df_t_out %>% mutate("p_val_adj" = p.adjust(t_p, method="fdr"))
df_t_sig <- df_t_adj[df_t_adj$p_val_adj < 0.05,]

all_per_pat <- data.frame(table(df_t_out$pattern))
sig_per_pat <- data.frame(table(df_t_sig$pattern))

df_frac <- data.frame("pattern" = all_per_pat$Var1, "fraction" = sig_per_pat$Freq / all_per_pat$Freq)


h5res = readRDS(file = here('microarray_analysis_d100/generated_data_V4/formattedresults_d100_lme4model_array.rds'))


stop()
df_w_out <- data.frame()
for(g in genes){
        d0_dat <- exprs_mat[d0_sub, g]
        d100_dat <- exprs_mat[d100_sub, g]
        w_res <- wilcox.test(d0_dat, d100_dat)
        temp_df <- data.frame("gene" = g, "w_p" = w_res$p.value)
        df_w_out <- rbind(df_w_out, temp_df)
}

df_w_adj <- df_w_out %>% mutate("p_val_adj" = p.adjust(w_p, method="fdr"))
df_w_sig <- df_w_adj[df_w_adj$p_val_adj < 0.05,]


