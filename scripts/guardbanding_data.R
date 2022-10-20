# Guardbanding Study

# input file data

a <- list.files(path='/Volumes/ivd/raw', pattern='AH2YKTBGXK')
sequencing_data <- read.csv('/Users/kkwock/OneDrive - Guardant Health/Data/Guardbanding/H2YKTBGXK_autoqc_sample_qc.hdr.tsv', sep = '\t') 
seq.filepath <- paste0('/Volumes/ivd/raw/', a)

fcid <- c("AH2YKTBGXK")  
sample_metrics <- list()
metric_values <- c(sample_coverage_exceptions, sample_contamination_pct, sample_gc_bias, sample_non_singleton_families, sample_avg_family_size, sample_on_target_rate, sample_gender_status_mismatch,  sample_germline_contamination)

sample_metrics <-  sequencing_data %>% 
  mutate(FCID = fcid)  %>%
  select(FCID, run_sample_id, ) %>% 
  unique()

colnames(sample_metrics) <- c("FCID", "run_sample_id", "on_target_reads", "total_reads", "raw_reads", "MAPD", "Coverage_Exceptions", "Family_contamination", "GC_IQR", "Non_singleton_families","Average_family_size", "On_target_rate", "Gender_status_mismatch", "Germline_contamination", "Batch_type", "FC_type")
sample_metrics


sample_metrics_ce_gciqr_select_long <-  sample_metrics_ce_gciqr_select %>% 
  gather(trait, value, c(on_target_reads:Germline_contamination)) %>% 
  mutate(value = as.numeric(value)) 

sample_metrics_ce_gciqr_select_long

## Plot the raw reads, total reads, and on target reads

traits <- c("on_target_reads", "total_reads", "raw_reads")
sample_metrics_ce_gciqr_select_long_selected_traits <- sample_metrics_ce_gciqr_select_long %>% filter(trait %in% traits) %>% 
  filter(Batch_type == "P360Q-669")
sample_metrics_ce_gciqr_select_long_selected_traits
sample_metrics_ce_gciqr_select_long_selected_traits$trait <- factor(sample_metrics_ce_gciqr_select_long_selected_traits$trait, levels = c("on_target_reads", "total_reads", "raw_reads"))

reads_barplot <- ggplot(sample_metrics_ce_gciqr_select_long_selected_traits, aes(x = run_sample_id, y =value, fill = FC_type )) +
  geom_bar(stat="identity")+
  facet_wrap(FC_type ~ trait, nrow = 2, scales = 'free_x') +
  scale_fill_manual(values=c("steelblue", "coral")) +
  labs(title = "LDT50_P360Q-669", colour = "Sample Type") +
  ylab("Reads value") +
  xlab("Sample ID") +
  theme_Publication() +
  theme(legend.position = 'bottom',
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 29),
        axis.text.x = element_text(angle = 45, size = 10 , hjust = 1))
reads_barplot
