# PURPOSE: To plot titer response against time (Day 0, Day 21, Day 28, Day 42 &
# Day 100)

# Initialize the environment with path to PROJECT_DIR and some commonly used
# packages.
source(file.path("SCRIPTS/0_initialize.r"))

# FUNCTIONS
factor.date <- function(x, frmt="%m/%d/%y") {
        factor(x,
        levels = unique(x[
        order(as.Date(as.character(x),format=frmt))]))
}

# MAIN
# Load data.
fn.mn = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Titres", "titre.txt")
df.mn = read.table(fn.mn, sep="\t", header=T, row.names=NULL, stringsAsFactors=T) %>% 
  mutate_at(-(1:2), function(x) sub("<20","10",x) %>% as.numeric()) %>% 
  mutate(day = sub("D","d",TimePt)) %>% 
  mutate(day = ifelse(day %in% c("d29","d31"), "d28", day) %>% factor.date(), subject=factor(Sample.ID))

yticks = sort(unique(df.mn$A.Indonesia))
yticks = yticks[seq(2,length(yticks),2)]

# Check if figure output directory already exists else create it.
dn.fig = file.path(PROJECT_DIR, "FIGURES/titers")
if(dir.exists(dn.fig)){
        print ("Figure directory already exists.")
} else {
        dir.create(dn.fig, recursive = T, showWarnings = T)
}

# Plot and save the figure file.
print ("Plotting and saving figure.")
plot_out  <- ggplot(df.mn, aes(day, log2(A.Indonesia), group=subject, col=subject)) +
  geom_line(position = position_jitter(width=0, height=0.1)) +
  scale_y_continuous(name="MN titer (A/Indonesia)", breaks=log2(yticks), labels = yticks) +
  xlab("Time points") +
  theme_bw() + theme(legend.position = "none", panel.grid.minor=element_blank(),
                     axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) 

ggsave(file.path(dn.fig, "MN_titer_profiles_all_subjects_indonesia.png"), plot = plot_out, h=3.3, w=4)
unlink(plot_out)
