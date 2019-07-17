# PURPOSE: To plot titer response rate against time.

# Initialize the environment with PROJETCT_DIR and some commonly used packages.
source("SCRIPTS/0_initialize.r")

# Read demographics and adjuvant status data.
fn.demo = file.path(PROJECT_DIR, "DATA_ORIGINAL/Clinical/clinical_info_adj.txt")
df.demo = fread(fn.demo, data.table = F) %>% 
  select(subject=`Subject ID`, Adjuvant)

# Read serology data.
fn.hai = file.path(PROJECT_DIR, "DATA_ORIGINAL", "HAI", "H5N1_serology.txt")
df.hai = read.table(fn.hai, sep="\t", 
                      header=T, row.names=NULL, stringsAsFactors=F) %>% 
  mutate_at(-(1), function(x) sub("<10","5",x) %>% as.numeric()) %>% 
  gather("tmp","hai",-subject) %>% 
  separate("tmp",c("day","replicate"))

df.hai3 = df.hai %>% 
  group_by(subject, day) %>% 
  summarise(hai.median=median(hai, na.rm=T), 
            hai.min=min(hai, na.rm=T), hai.max=max(hai, na.rm=T) ) %>% 
  ungroup() %>% 
  mutate(day=factor(day,levels=c("day0","day1","day7","day21","day22","day28","day42","day100")))

df.hai3.time = df.hai3 %>% 
  select(-hai.min, -hai.max) %>% 
  filter(day %in% c("day0", "day21", "day28", "day42")) %>% 
  spread("day","hai.median")

df.hai3.fc = df.hai3.time %>% 
  mutate_at(-1, function(x)ifelse(is.na(x),5,x)) %>% 
  mutate_at(-1, function(x)x/.$day0)

df.resp = df.hai3.fc %>% 
  inner_join(df.demo, by="subject") %>% 
  mutate(Adjuvant=factor(Adjuvant, levels=c("NonAdj","Adj"))) %>% 
  select(-day0) %>% 
  gather("Time","hai",-subject, -Adjuvant) %>% 
  group_by(Time, Adjuvant) %>% 
  summarize(Response.rate=sum(hai>=4)/n()*100) %>% 
  ungroup()

# Invoke custom function for color palette. 
source(file.path(PROJECT_DIR, "SCRIPTS/functions/color_functions.r"))
cm = gg_color_hue(2) %>% rev()

# Plot and save the figure file
plot_rate <- ggplot(df.resp, aes(Time, Response.rate, fill=Adjuvant, col=Adjuvant)) +
  geom_bar(stat="identity", position=position_dodge()) +
  # geom_errorbar(aes(ymin=Response.rate.CI.low, ymax=Response.rate.CI.high), 
  #               col="black", position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = cm) +
  scale_color_manual(values = cm) +
  xlab("Days after first vaccination") +
  ylab("Response rate, %") +
  theme_bw()

## Check if the figure output directory already exists else create it.
dn.fig = file.path(PROJECT_DIR, "FIGURES/titers")
if(dir.exists(dn.fig)){
        print("Figure directory already exists.")
}else {
        print("Figure directory doesn't exists, creating it now.")
        dir.create(dn.fig, recursive = T, showWarnings = T)
}
print("Saving figure file.")
fn.fig = file.path(dn.fig, "titer_response_4FC")
ggsave(paste0(fn.fig, ".png"), plot = plot_rate, w=5, h=3)
ggsave(paste0(fn.fig, ".pdf"), plot = plot_rate, w=5, h=3)
unlink(plot_rate)
