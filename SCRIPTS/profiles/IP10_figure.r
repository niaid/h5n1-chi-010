# PURPOSE: To plot IP-10 data.

# Initialize the environment with PROJECT_DIR and commonly used packages.
source("SCRIPTS/0_initialize.r")

# MAIN
# Load luminex data.
dn.in = file.path(PROJECT_DIR, "DATA_PROCESSED/Luminex")
df = readRDS(file.path(dn.in, "luminex_data.rds")) %>% 
  filter(analyte == "IP-10")

df.mean = df %>% 
  group_by(time) %>% 
  summarise(value=mean(value, na.rm=T)) %>% 
  ungroup()

# Plot and save figure files.
plot_ip  <- ggplot(df, aes(time, value, group=subject)) +
  geom_line(alpha=0.2) +
  geom_line(data=df.mean, group=1, size=1) +
  # geom_vline(xintercept = c(5, 11), lty=2, col="black") +
  # geom_vline(xintercept = c(1, 7), lty=1, col="blue") +
  xlab("Time points") +
  ylab("IP-10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

## Check if the figure output folder exists.
dn.fig = file.path(PROJECT_DIR, "FIGURES/profiles")
if(dir.exists(dn.fig)){
        print("Figure output folder exists.")
}else {
        print("Fgiure output folder doesn't exits. Creating it now.")
        dir.create(dn.fig, recursive = T, showWarnings = T)
}
fn.fig = file.path(dn.fig, "IP-10_profiles")
ggsave(paste0(fn.fig, ".png"), plot = plot_ip, h=3.5, w=4)
ggsave(paste0(fn.fig, ".pdf"), plot = plot_ip, h=3.5, w=4)
