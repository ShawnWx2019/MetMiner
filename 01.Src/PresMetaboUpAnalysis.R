#######################################################################
#         Prj: MetMiner
#         Assignment: UHPLC analysis of pseudotargeted metabolomics
#         Author: Shawn Wang
#         Init Date: Jan 03, 2023
#         Update: Jul 14, 2023
#######################################################################
#> log: 
#> Jul 14. fix bugs for ms2 selection.
#> Jul 14. ready to push
#> getopt
TEST = "FALSE"
options(stringsAsFactors = F)
options(warn = -1)
suppressMessages(if (!require('getopt')) BiocManager::install('getopt'))
suppressMessages(if (!require('crayon')) BiocManager::install('crayon'))



msg_yes = green$bold$italic
msg_no = red$bold$italic
msg_run = blue$bold$italic$underline
msg_warning = yellow$bold$italic


# args --------------------------------------------------------------------

command=matrix(c(
  'help', 'h', 0, 'logic', 'help information',
  'MS1', 'x', 1, 'character', 'File Path of MS1". ',
  'MS2', 'y', 1, 'character','File Path of MS2',
  'sample_info', 'z', 1, 'character','sample information',
  'column', 'c', 1, 'character','column type rp or hilic',
  'ppm', 'm', 1, 'double', 'Related to m/z\n                         fluctuation of m/z value (ppm) from scan to scan - depends on the mass spectrometer accuracy\n                         eg:10',
  'threads', 't', 1, 'integer', 'How many threads do you want to use.\n                          eg:20\n                          default:8',
  'show_output', 'o', 2, 'Logic', 'Need tic, bpc and rt_correction_plot or not\n                          default: TRUE',
  'p_min', 'a', 2, 'integer', 'Minimal of peakwidth for peak detection.\n                          eg:10\n                          default:10',
  'p_max', 'b', 2, 'integer', 'Max of peakwidth for peak detection.\n                          eg:10\n                          default:10',
  'min_fraction', 'f', 2, 'double', 'Related to Samples,\n                         to be valid, a group must be found in at least minFraction*n samples, with n=number of samples for each class of samples. A minFraction=0.5 corresponds to 50%.\n                         n=10, minFraction=0.5 => found in at least 5 samples\n                         default:0.5',
  'QC_tag', 'g', 2, 'character', 'File path of QCn.mzXML file placed.\n                         eg: QC.mzXML files placed in NEG/QC, set -g "QC", QC.mzXML files placed in NEG/qc, set -g "qc"',
  'snthresh', 's', 2, 'double', 'Related to intensity,\n                         signal/noise ratio threshold\n                         default: 5',
  'noise', 'n', 2, 'double', 'Related to intensity,\n                         each centroid must be greater than the “noise” value\n                         default:500'
),byrow = T, ncol = 5)
args = getopt(command)

## help information
if (!is.null(args$help)) {
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status=1)
}

## default value
if (is.null(args$MS1)){
  message(msg_no("-p error \nPlease fill in the file path of MS1 files correctly!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}
if (is.null(args$MS2)){
  message(msg_no("-n error \nPlease fill in the file path of MS2 files correctly!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}
if (is.null(args$sample_info)){
  message(msg_no("-n error \nPlease fill in the file path of sample_info files correctly!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}
if (is.null(args$QC_tag)){
  message(msg_no("-g error \nPlease fill in the file path of QC files correctly!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}

if (is.null(args$ppm)){
  args$ppm = 10
}

if (is.null(args$column)){
  args$column = 'rp'
}

if (is.null(args$threads)){
  args$threads = 8
}

if (is.null(args$snthresh)){
  args$snthresh = 5
}

if (is.null(args$noise)){
  args$noise = 500
}

if (is.null(args$min_fraction)){
  args$min_fraction = 0.5
}


if (is.null(args$p_min)){
  args$p_min = 10
}

if (is.null(args$p_max)){
  args$p_max = 60
}

if (is.null(args$show_output)){
  args$show_output = FALSE
}


# test --------------------------------------------------------------------
##test
if (TEST == "TRUE") {
  T.MS1 = "02.DemoData/MS1";
  T.MS2 = "02.DemoData/MS2";
  T.sample_info = "02.DemoData/sample_info.xlsx"
  T.ppm = 10;
  T.threads = 14;
  T.show_output = FALSE;
  T.p_min = 10;
  T.p_max = 60;
  T.min_fraction = 0.5;
  T.QC_tag = "QC";
  T.column = 'rp';
  T.snthresh = 5;
  T.noise = 500;
  ms1.pos <- paste0(T.MS1,"/POS/")
  ms1.neg <- paste0(T.MS1,"/NEG/")
  ms2.pos <- paste0(T.MS2,"/POS/")
  ms2.neg <- paste0(T.MS2,"/NEG/")
} else {
  T.MS1 = args$MS1;
  T.MS2 = args$MS2;
  T.ppm = args$ppm;
  T.sample_info = args$sample_info;
  T.threads = args$threads;
  T.show_output = args$show_output;
  T.p_min = args$p_min;
  T.p_max = args$p_max;
  T.min_fraction = args$min_fraction;
  T.QC_tag = args$QC_tag;
  T.snthresh = args$snthresh;
  T.noise = args$noise
  T.column = args$column
  ms1.pos <- paste0(T.MS1,"/POS/")
  ms1.neg <- paste0(T.MS1,"/NEG/")
  ms2.pos <- paste0(T.MS2,"/POS/")
  ms2.neg <- paste0(T.MS2,"/NEG/")
}
#> 
library(tidymass)
suppressMessages(if (!require('devtools')) BiocManager::install('devtools'))
suppressMessages(if (!require('conflicted')) BiocManager::install('conflicted'))
suppressMessages(if (!require('tidyverse')) BiocManager::install('tidyverse'))
suppressMessages(if (!require('patchwork')) BiocManager::install('patchwork'))
suppressMessages(if (!require('bruceR')) install.packages("bruceR", dep=TRUE))
suppressMessages(if (!require('MDAtoolkits')) install_github(repo = "ShawnWx2019/MDAtoolkits",ref = 'master'))
suppressMessages(if (!require('IMOtoolkits')) install_github(repo = "ShawnWx2019/IMOtoolkits"))

suppressMessages(conflict_prefer_all("dplyr"))
#tidymass::update_tidymass(packages = "all",from = "gitlab")

# process data ------------------------------------------------------------

message(msg_run("Step1. Peak peaking in positive model. This may take a long time. Please be patient!"))
process_data(
  path = ms1.pos,
  polarity = "positive",
  ppm = T.ppm,
  threads = T.threads,
  snthresh = T.snthresh,
  noise = T.noise,
  peakwidth = c(T.p_min,T.p_max),
  output_tic = T.show_output,
  output_bpc = T.show_output,
  output_rt_correction_plot = T.show_output,
  min_fraction = T.min_fraction,
  fill_peaks = FALSE,
  group_for_figure = T.QC_tag
)

##> NEG model
##>
message(msg_run("Step2. Peak peaking in negative model. This may take a long time. Please be patient!"))
process_data(
  path = ms1.neg,
  polarity = "negative",
  ppm = T.ppm,
  threads = T.threads,
  snthresh = T.snthresh,
  noise = T.noise,
  peakwidth = c(T.p_min,T.p_max),
  output_tic = T.show_output,
  output_bpc = T.show_output,
  output_rt_correction_plot = T.show_output,
  min_fraction = T.min_fraction,
  fill_peaks = FALSE,
  group_for_figure = T.QC_tag
)

message(msg_yes("Peak picking steps fininsed. Please check result in /NEG/Result or /POS/Result"))



# functions ---------------------------------------------------------------
##> theme

theme1 =   theme(
  panel.border = element_rect(size = 1.5),
  axis.title = element_text(size = 14,color = 'black'),
  axis.text = element_text(size = 12,color = 'black')
)

##> QC sample batch effect detection.
batch_detect = function(object) {
  plt_batch =
    object %>%
    extract_expression_data() %>%
    dplyr::select(contains("QC")) %>%
    pivot_longer(contains("QC"),values_to = "value",names_to = "QC_samples") %>%
    drop_na() %>%
    mutate(
      value = log2(value),
      QC_samples = paste0(
        "QC",str_pad(
          string = gsub("QC","",QC_samples),width = 3,side = 'left',pad = "0"
        )
      )
    ) %>%
    arrange(QC_samples) %>%
    ggplot(data = .,mapping = aes(x = QC_samples,y = value,color = QC_samples))+
    geom_boxplot()+
    ylab("log2(Raw peak area)")+
    xlab("")+
    theme_bw()+
    theme1+
    theme(
      axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1),
      legend.position = "none"
    )
  return(plt_batch)
}

## peak distribution
peak_distribution <- function(object) {
  plt_peak_dis =
    object %>%
    `+`(1) %>%
    log(10) %>%
    show_mz_rt_plot() +
    scale_size_continuous(range = c(0.01, 2))+
    theme1
  return(plt_peak_dis)
}

##> working dir
dir.create("workdir",showWarnings = F,recursive = T)
# 1.0 positive model ----------------------------------------------------------

##> import mass_dataset
load(paste0(ms1.pos,"Result/object"))

object.pos <- object
##> load sample_info
sample_info <- import(T.sample_info)
object.pos <-
  object.pos %>% 
  activate_mass_dataset('sample_info') %>% 
  select(sample_id) %>% left_join(sample_info)

dir.create("workdir/01.raw",showWarnings = F,recursive = T)
save(object.pos,file = "workdir/01.raw/object.pos")

##> batch detection
batch_detect_plt.pos <- 
  batch_detect(object = object.pos)

ggsave(filename = "workdir/01.raw/raw_peak_boxplot.pos.png",plot = batch_detect_plt.pos,width = 10,height = 5)
ggsave(filename = "workdir/01.raw/raw_peak_boxplot.pos.pdf",plot = batch_detect_plt.pos,width = 10,height = 5)

##> peak distribution

peak_distribution_plt.pos <-
  peak_distribution(object = object.pos)
ggsave(filename = "workdir/01.raw/peak_distribution_plt.pos.png",plot = peak_distribution_plt.pos,width = 10,height = 5)
ggsave(filename = "workdir/01.raw/peak_distribution_plt.pos.pdf",plot = peak_distribution_plt.pos,width = 10,height = 5)

##> missing value
plt_mv_raw.pos<-
  show_sample_missing_values(object = object.pos, percentage = TRUE,color_by = 'group',order_by = 'injection.order')+theme1+
  geom_hline(yintercept = 80,color = "red",linetype="dashed")+
  ylim(c(0,100))+
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()+
  theme(axis.text.x = element_text(size = 3))

ggsave(filename = "workdir/01.raw/missing_value_distribution.pos.png",plot = plt_mv_raw.pos,width = 10,height = 5)
ggsave(filename = "workdir/01.raw/missing_value_distribution.pos.pdf",plot = plt_mv_raw.pos,width = 10,height = 5)

# 1.2 Outlier detect and remove outliers ----------------------------------

qc_id = object.pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)

object.pos <-
  object.pos %>%
  mutate_variable_na_freq(according_to_samples = qc_id)
##> na frequence is less than 0.2 in qc samples.
object.pos.mv <-
  object.pos %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(na_freq < 0.2)

dir.create("workdir/02.remove_noise/",showWarnings = F,recursive = T)
save(object.pos.mv,file = "workdir/02.remove_noise/object.pos.mv")

plt_mv_remove_noise.pos<-
  show_sample_missing_values(object = object.pos.mv, percentage = TRUE,color_by = 'group',order_by = 'injection.order')+theme1+
  geom_hline(yintercept = 80,color = "red",linetype="dashed")+
  ylim(c(0,100))+
  theme(axis.text.x = element_text(size = 8,angle = 90)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()+
  theme(axis.text.x = element_text(size = 3))

ggsave(filename = "workdir/02.remove_noise/plt_mv_remove_noise.pos.png",plot = plt_mv_remove_noise.pos,width = 10,height = 5)
ggsave(filename = "workdir/02.remove_noise/plt_mv_remove_noise.pos.pdf",plot = plt_mv_remove_noise.pos,width = 10,height = 5)

outlier_samples.pos <-
  object.pos.mv %>% 
  `+`(1) %>% 
  log(2) %>% 
  scale() %>% 
  detect_outlier(na_percentage_cutoff = 0.8)

outlier_table.pos <-
  extract_outlier_table(outlier_samples.pos)

out_name.pos <-
  outlier_table.pos %>%
  filter(according_to_na == TRUE) %>% 
  rownames()

##> remove outlier based on 
if(length(out_name.pos) != 0) {
  object.pos.mv <- 
    object.pos.mv %>% 
    activate_mass_dataset('expression_data') %>% 
    select(-all_of(out_name.pos))
}

plt_mv_remove_outlier.pos<-
  show_sample_missing_values(object = object.pos.mv, percentage = TRUE,color_by = 'group',order_by = 'injection.order')+theme1+
  geom_hline(yintercept = 80,color = "red",linetype="dashed")+
  ylim(c(0,100))+
  theme(axis.text.x = element_text(size = 8,angle = 90)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()+
  theme(axis.text.x = element_text(size = 3))

ggsave(filename = "workdir/02.remove_noise/plt_mv_remove_outlier.pos.png",plot = plt_mv_remove_outlier.pos,width = 10,height = 5)
ggsave(filename = "workdir/02.remove_noise/plt_mv_remove_outlier.pos.pdf",plot = plt_mv_remove_outlier.pos,width = 10,height = 5)


# impute mv ---------------------------------------------------------------
object.pos.impute =
  object.pos.mv %>%
  impute_mv(method = 'knn')
dir.create("workdir/03.impute_mv/",showWarnings = F,recursive = T)
save(object.pos.impute,file = "workdir/03.impute_mv/object.pos.impute")


# normalization -----------------------------------------------------------
dir.create("workdir/04.normalization/",recursive = T,showWarnings = F)

object_inte_pos <-
  integrate_data(object.pos.impute,'qc_mean')
save(object_inte_pos,file = "workdir/04.normalization/object_inte_pos")


# rsd ---------------------------------------------------------------------

rsd <- 
  object_inte_pos %>% 
  activate_mass_dataset("sample_info") %>%
  dplyr::filter(class == T.QC_tag) %>%
  extract_expression_data() %>%
  rownames_to_column("ID") %>%
  pivot_longer(contains(T.QC_tag),names_to = "tag",values_to = "value") %>%
  select(-tag) %>%
  group_by(ID) %>%
  summarise(
    raw.rsd = (sd(value,na.rm = T)/mean(value,na.rm = T))*100
  ) %>% arrange(raw.rsd) %>% 
  mutate(
    check = case_when(
      raw.rsd > 30 ~ "Failed",
      TRUE ~ "PASS"
    )
  )

rsd_check <- 
  rsd %>% 
  group_by(check) %>% 
  summarise(num = n())

writexl::write_xlsx(x = list(
  rsd_val = rsd,
  rsd_summary = rsd_check
),"workdir/04.normalization/feature_rsd.pos.xlsx")


# add MS2 -----------------------------------------------------------------
#> add ms2 info
object_pos.ms2 <-
  object_inte_pos %>%
  mutate_ms2(
    object = .,
    column = T.column,
    polarity = 'positive',
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = ms2.pos
  )
#> extract features with ms2
object_pos_MRM <- 
  oneStepMRMselection(obj_ms2 = object_pos.ms2)

# feature annotation ------------------------------------------------------

##> load database 
##> MS2


load("~/.HPC_tidymass/MS_db/mona_database0.0.3.rda")
load("~/.HPC_tidymass/MS_db/RIKEN_PlaSMA_database0.0.1.rda")
load("~/.HPC_tidymass/MS_db/knapsack_ath_db.rda")
load("~/.HPC_tidymass/MS_db/kegg_ms1_database0.0.3.rda")
load("~/.HPC_tidymass/MS_db/RPLC.database.rda")

object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos_MRM,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = T.column,
    threads = 5,
    database = kegg_ms1_database0.0.3,
    candidate.num = 2
  )



object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos_anno,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = T.column,
    threads = 5,
    database = knapsack_ath_db,
    candidate.num = 2
  )
# for ms2

object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos_anno,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = T.column,
    threads = T.threads,
    database = snyder_database_rplc0.0.3,
    candidate.num = 2
  )

object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos_anno,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = T.column,
    threads = T.threads,
    database = mona_database0.0.3,
    candidate.num = 2
  )


object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos_anno,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = T.column,
    threads = T.threads,
    database = RIKEN_PlaSMA_database0.0.1,
    candidate.num = 2
  )

save(object_pos_anno,file = "workdir/05.Annotation/object_pos_anno.rds")

anno_tmp <- 
object_pos_anno %>% extract_annotation_table()


# neg ---------------------------------------------------------------------

load(paste0(ms1.neg,"Result/object"))

object.neg <- object
##> load sample_info
sample_info <- import(T.sample_info)
object.neg <-
  object.neg %>% 
  activate_mass_dataset('sample_info') %>% 
  select(sample_id) %>% left_join(sample_info)


save(object.neg,file = "workdir/01.raw/object.neg")

##> batch detection
batch_detect_plt.neg <- 
  batch_detect(object = object.neg)

ggsave(filename = "workdir/01.raw/raw_peak_boxplot.neg.png",plot = batch_detect_plt.neg,width = 10,height = 5)
ggsave(filename = "workdir/01.raw/raw_peak_boxplot.neg.pdf",plot = batch_detect_plt.neg,width = 10,height = 5)

##> peak distribution

peak_distribution_plt.neg <-
  peak_distribution(object = object.neg)
ggsave(filename = "workdir/01.raw/peak_distribution_plt.neg.png",plot = peak_distribution_plt.neg,width = 10,height = 5)
ggsave(filename = "workdir/01.raw/peak_distribution_plt.neg.pdf",plot = peak_distribution_plt.neg,width = 10,height = 5)

##> missing value
plt_mv_raw.neg<-
  show_sample_missing_values(object = object.neg, percentage = TRUE,color_by = 'group',order_by = 'injection.order')+theme1+
  geom_hline(yintercept = 80,color = "red",linetype="dashed")+
  ylim(c(0,100))+
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()+
  theme(axis.text.x = element_text(size = 3))

ggsave(filename = "workdir/01.raw/missing_value_distribution.neg.png",plot = plt_mv_raw.neg,width = 10,height = 5)
ggsave(filename = "workdir/01.raw/missing_value_distribution.neg.pdf",plot = plt_mv_raw.neg,width = 10,height = 5)

# 1.2 Outlier detect and remove outliers ----------------------------------

qc_id = object.neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)

object.neg <-
  object.neg %>%
  mutate_variable_na_freq(according_to_samples = qc_id)
##> na frequence is less than 0.2 in qc samples.
object.neg.mv <-
  object.neg %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(na_freq < 0.2)


save(object.neg.mv,file = "workdir/02.remove_noise/object.neg.mv")

plt_mv_remove_noise.neg<-
  show_sample_missing_values(object = object.neg.mv, percentage = TRUE,color_by = 'group',order_by = 'injection.order')+theme1+
  geom_hline(yintercept = 80,color = "red",linetype="dashed")+
  ylim(c(0,100))+
  theme(axis.text.x = element_text(size = 8,angle = 90)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()+
  theme(axis.text.x = element_text(size = 3))

ggsave(filename = "workdir/02.remove_noise/plt_mv_remove_noise.neg.png",plot = plt_mv_remove_noise.neg,width = 10,height = 5)
ggsave(filename = "workdir/02.remove_noise/plt_mv_remove_noise.neg.pdf",plot = plt_mv_remove_noise.neg,width = 10,height = 5)

outlier_samples.neg <-
  object.neg.mv %>% 
  `+`(1) %>% 
  log(2) %>% 
  scale() %>% 
  detect_outlier(na_percentage_cutoff = 0.8)

outlier_table.neg <-
  extract_outlier_table(outlier_samples.neg)

out_name.neg <-
  outlier_table.neg %>%
  filter(according_to_na == TRUE) %>% 
  rownames()

##> remove outlier based on 
if(length(out_name.neg) != 0) {
  object.neg.mv <- 
    object.neg.mv %>% 
    activate_mass_dataset('expression_data') %>% 
    select(-all_of(out_name.neg))
}

plt_mv_remove_outlier.neg<-
  show_sample_missing_values(object = object.neg.mv, percentage = TRUE,color_by = 'group',order_by = 'injection.order')+theme1+
  geom_hline(yintercept = 80,color = "red",linetype="dashed")+
  ylim(c(0,100))+
  theme(axis.text.x = element_text(size = 8,angle = 90)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()+
  theme(axis.text.x = element_text(size = 3))

ggsave(filename = "workdir/02.remove_noise/plt_mv_remove_outlier.neg.png",plot = plt_mv_remove_outlier.neg,width = 10,height = 5)
ggsave(filename = "workdir/02.remove_noise/plt_mv_remove_outlier.neg.pdf",plot = plt_mv_remove_outlier.neg,width = 10,height = 5)


# impute mv ---------------------------------------------------------------
object.neg.impute =
  object.neg.mv %>%
  impute_mv(method = 'knn')
save(object.neg.impute,file = "workdir/03.impute_mv/object.neg.impute")


# normalization -----------------------------------------------------------

object_inte_neg <-
  integrate_data(object.neg.impute,'qc_mean')
save(object_inte_neg,file = "workdir/04.normalization/object_inte_neg")


# rsd ---------------------------------------------------------------------

rsd <- 
  object_inte_neg %>% 
  activate_mass_dataset("sample_info") %>%
  dplyr::filter(class == T.QC_tag) %>%
  extract_expression_data() %>%
  rownames_to_column("ID") %>%
  pivot_longer(contains(T.QC_tag),names_to = "tag",values_to = "value") %>%
  select(-tag) %>%
  group_by(ID) %>%
  summarise(
    raw.rsd = (sd(value,na.rm = T)/mean(value,na.rm = T))*100
  ) %>% arrange(raw.rsd) %>% 
  mutate(
    check = case_when(
      raw.rsd > 30 ~ "Failed",
      TRUE ~ "PASS"
    )
  )

rsd_check <- 
  rsd %>% 
  group_by(check) %>% 
  summarise(num = n())

writexl::write_xlsx(x = list(
  rsd_val = rsd,
  rsd_summary = rsd_check
),"workdir/04.normalization/feature_rsd.neg.xlsx")


# add MS2 -----------------------------------------------------------------
#> add ms2 info
object_neg.ms2 <-
  object_inte_neg %>%
  mutate_ms2(
    object = .,
    column = T.column,
    polarity = 'negative',
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = ms2.neg
  )
#> extract features with ms2
object_neg_MRM <- 
  oneStepMRMselection(obj_ms2 = object_neg.ms2)

# feature annotation ------------------------------------------------------

##> load database 
##> MS2

object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg_MRM,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = T.column,
    threads = 5,
    database = kegg_ms1_database0.0.3,
    candidate.num = 2
  )


object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg_anno,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = T.column,
    threads = 5,
    database = knapsack_ath_db,
    candidate.num = 2
  )
# for ms2

object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg_anno,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = T.column,
    threads = T.threads,
    database = snyder_database_rplc0.0.3,
    candidate.num = 2
  )

object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg_anno,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = T.column,
    threads = T.threads,
    database = mona_database0.0.3,
    candidate.num = 2
  )

object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg_anno,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = T.column,
    threads = T.threads,
    database = RIKEN_PlaSMA_database0.0.1,
    candidate.num = 2
  )

save(object_neg_anno,file = "workdir/05.Annotation/object_neg_anno.rds")



# merge object ------------------------------------------------------------

dir.create("workdir/05.Annotation/",recursive = T,showWarnings = F)
load("~/.HPC_tidymass/MS_db/labID2INCHIKEY.rda")
object_merge_original <- 
  merge_mass_dataset(
    x = object_neg_anno,
    y = object_pos_anno,
    sample_direction = "inner",
    variable_direction = 'full',
    sample_by = 'sample_id',
    variable_by = c("variable_id","mz","rt","ms2_spectrum_id","precursor","product")
  ) %>% 
  activate_mass_dataset("sample_info") 


ori_vari_info <-
  object_merge_original %>% 
  extract_variable_info() %>% 
  select(variable_id,mz,rt,Compound.name,Lab.ID,ms2_spectrum_id,precursor,product,Adduct,mz.error,Total.score,Level) %>% 
  distinct() %>% 
  filter(Adduct == "(M-H)-" |
           Adduct == "(M-H2O-H)-"|
           Adduct == "(M+NH4-2H)-"|
           Adduct == "(2M-H)-"|
           Adduct == "(M+F)-" |
           Adduct == "(M+H)+" |
           Adduct == "(2M+H)+"|
           Adduct == "(M+H-H2O)+"|
           Adduct == "(M+NH4)+"|
           Adduct == "(M+H-2H2O)+" |
           Adduct == "(M+H-2H2O)+"|
           Adduct == "(2M+NH4)+"
  ) %>% 
  mutate(Adduct.judge = case_when(
    Adduct == "(M+H)+" ~ 1,
    Adduct == "(M-H)-" ~ 1,
    Adduct == "(M+F)-" ~ 2,
    TRUE ~ 3
  )) %>% 
  group_by(variable_id) %>% 
  slice_min(Adduct.judge) %>% 
  slice_max(Total.score) %>% 
  arrange(Compound.name) %>% 
  mutate(query = URLencode(Compound.name)) %>% 
  ungroup
MRM_method <- ori_vari_info %>% 
  select(variable_id,mz,rt,precursor,product)
writexl::write_xlsx(list(annotation = ori_vari_info,TQMS_method = MRM_method),"workdir/05.Annotation/pseudotargeted_QC_anno.xlsx")



save.image(paste0("workdir/",Sys.time() %>% str_replace_all("-|:| ","_"),".rda"))

#load("workdir/2023_07_16_09_54_55.035713.rda")
