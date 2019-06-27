dn.in = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_10c/")
fn.ann = file.path(dn.in, "flow_10c_info.txt")
ann = fread(fn.ann, data.table = F)

files = c("H5N1_Bcell-H5N1_Bcell_freq.V1.1",
          "H5N1_DCMonNK-DCMonNK_freq.V1.1",
          "H5N1_Tcell-H5N1-Tcell_freq.V1.1",
          "H5N1_Thelper-Thelper_freq.V1.1",
          "H5N1_Treg-H5N1-Treg_freq.V1.1") %>% 
  file.path(dn.in, .)

flow.d0 = data.frame()
for(f in files) {
  tmp = fread(f, header=T, data.table=F) %>% 
    dplyr::rename(donor=V1) %>% 
    filter(Timepoint == "d0_0h") %>% 
    dplyr::select(-Timepoint) %>% 
    gather(ID, freq, -donor) %>% 
    mutate(ID = as.numeric(ID)) %>% 
    filter(!is.na(ID))
  flow.d0 = rbind(flow.d0, tmp)
}

flow.d0 = flow.d0 %>% 
  inner_join(ann %>% dplyr::select(Id, Code), by=c("ID"="Id")) %>% 
  dplyr::select(-ID) %>% 
  spread(Code, freq)

fn.out = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_10c/flow_10c_matrix_Code_d0.txt")
fwrite(flow.d0, file=fn.out, sep="\t")

# ff = fread(file.path(PROJECT_DIR, "DATA_PROCESSED/eNet/eNet_InputData_NEW.txt"), header=T, data.table = F) %>% 
#   select(donor, matches("^Fb\\d\\d"))
# 
# all.equal(flow.d0$donor, ff$donor)
# plot(flow.d0$Fb61, ff$Fb61)
