#######################################################################
#         Prj: MetMiner
#         Assignment: TQMS Data Normalization
#         Author: Shawn Wang
#         Init Date: Jan 03, 2023
#         Update: Jul 14, 2023
#######################################################################
TEST = "FALSE"
options(stringsAsFactors = F)
options(warn = -1)
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
suppressMessages(if (!require('getopt')) BiocManager::install('getopt'))
suppressMessages(if (!require('crayon')) BiocManager::install('crayon'))

msg_yes = green$bold$italic
msg_no = red$bold$italic
msg_run = blue$bold$italic$underline
msg_warning = yellow$bold$italic


# args --------------------------------------------------------------------

command=matrix(c(
  'help', 'h', 0, 'logic', 'help information',
  'expmat', 'e', 1, 'character', 'File Path of MS1". ',
  'sample_info', 's', 1, 'character','sample information',
  'color_by','c',1, 'In PCA biplot, samples will be colored by which column in sample_info',
  'dimension','d',0,'how many dimentions will be used to generate PCA biplot | 2 or 3 ',
  'heterogeneous', 'g', 1, 'character', "yes or no",
  'norm_method','n',1,'character','normalization method: svr,serrf,total,median,mean,pqn,loess'
),byrow = T, ncol = 5)
args = getopt(command)

## help information
if (!is.null(args$help)) {
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status=1)
}

## default value
if (is.null(args$expmat)){
  message(msg_no("-p error \nPlease fill in the file path of expmat files correctly!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}
if (is.null(args$sample_info)){
  message(msg_no("-n error \nPlease fill in the file path of sample_info files correctly!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}
if (is.null(args$heterogeneous)){
  args$heterogeneous = 'no'
}

if (is.null(args$color_by)){
  args$color_by = 'group'
}

if (is.null(args$dimension)){
  args$dimension = 2
}

if (is.null(args$norm_method)){
  args$norm_method = 'svr'
}

##test
if (TEST == "TRUE") {
  T.expmat = "02.DemoData/exp_mat_v2.csv";
  T.sample_info = "02.DemoData/sample_info_v2.csv";
  T.norm_method = "serrf"
  heterogeneous = 'yes'
  color_by = 'group'
  dimension = 2
} else {
  T.expmat = args$expmat;
  T.sample_info = args$sample_info;
  T.norm_method = args$norm_method
  heterogeneous = args$heterogeneous
  color_by = args$color_by
  dimension = args$dimension
}

library(tidymass)
suppressMessages(if (!require('devtools')) install.packages('devtools'))
suppressMessages(if (!require('conflicted')) install.packages('conflicted'))
suppressMessages(if (!require('tidyverse')) install.packages('tidyverse'))
suppressMessages(if (!require('progressr')) install.packages('progressr'))
suppressMessages(if (!require('PCAtools')) install_github('kevinblighe/PCAtools'))
suppressMessages(if (!require('conflict')) install.packages('conflict'))
suppressMessages(if (!require('patchwork')) install.packages('patchwork'))
#suppressMessages(if (!require('ComplexHeatmap')) install_github('jokergoo/ComplexHeatmap'))
suppressMessages(if (!require('MDAtoolkits')) install_github(repo = "ShawnWx2019/MDAtoolkits",ref = 'master'))
suppressMessages(if (!require('IMOtoolkits')) install_github(repo = "ShawnWx2019/IMOtoolkits"))
# tryCatch(
#   {pacman::p_load(affy, parallel,ranger,caret,pcaMethods,ggplot2,tidyr,graphics,grDevices,Hmisc,gtools,cowplot,RColorBrewer,readr,plotly,stringr,GGally,dplyr,e1071,officer,bootstrap,pdist,metabolomics)
# }
# )
suppressMessages(if (!require('parallel')) install.packages("parallel"))
suppressMessages(if (!require('ranger')) install.packages("ranger"))
suppressMessages(conflict_prefer_all('dplyr'))
# importdata --------------------------------------------------------------
handlers(handler_pbcol(
  adjust = 1.0,
  complete = function(s) cli::bg_red(cli::col_black(s)),
  incomplete = function(s) cli::bg_cyan(cli::col_black(s))
))

raw_exp_mat <- read.csv(
  T.expmat
)

expression_data <- 
  raw_exp_mat %>% 
  select(-rt,-mz) %>% 
  column_to_rownames("variable_id")

variable_info <- 
  raw_exp_mat %>% 
  select(variable_id:rt)

sample_info <- read.csv(
  T.sample_info 
)

## 1.2 => construct a massdataset by tidymass.

object <- 
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

## 1.3 => save rawdata object

dir.create("workdir2/01.datacleaning",showWarnings = F,recursive = T)

save(object,file = "workdir2/01.datacleaning/object_raw.rds")

# step2.  data cleaning ---------------------------------------------------

## 2.2 => remove noisy features.

object %>% 
  activate_mass_dataset(
    what = "sample_info"
  ) %>% 
  count(class)


## noise remove method: QC and WT 0.2 and MT 0.5

## 2.2.1 Get sample ID of each group.
qc_id <- 
  object %>% 
  activate_mass_dataset(
    what = "sample_info"
  ) %>% 
  dplyr::filter(class == "QC") %>% 
  pull(sample_id)

subject_id <- 
  object %>% 
  activate_mass_dataset(
    what = "sample_info"
  ) %>% 
  dplyr::filter(class == "Subject") %>% 
  pull(sample_id)

## 2.2.2 Calculate missing frequency.
object <- 
  object %>% 
  mutate_variable_na_freq(according_to_samples = qc_id) %>% 
  mutate_variable_na_freq(according_to_samples = subject_id) 
## 2.2.3 remove noisy feature
if(heterogeneous=='no') {
  object.mv <- 
    object %>% 
    activate_mass_dataset(
      what = "variable_info"
    ) %>% 
    filter(na_freq < 0.2 & na_freq.1 < 0.5 )
} else {
  object.mv <- 
    object %>% 
    activate_mass_dataset(
      what = "variable_info"
    ) %>% 
    filter(na_freq < 0.2)
}

if(heterogeneous=='no') {
  outlier_samples <- 
    object %>% 
    `+`(1) %>% 
    log(2) %>% 
    scale() %>% 
    detect_outlier()
  outlier_table <-
    extract_outlier_table(outlier_samples) 
  outlier_name <- 
  outlier_table %>% 
    apply(1, function(x){
      sum(x)
    }) %>% 
    `>`(0) %>% 
    which() %>% names
  object.mv <- 
    object.mv %>% 
    activate_mass_dataset('expression_data') %>% 
    select(-all_of(outlier_name))
}

##> QC samples were concentrated for 5x, so the outliers were ignored

## 2.4 => missing value imputation 

object.imput <- 
  impute_mv(
    object = object.mv,
    method = "knn"
  )
get_mv_number(object.imput)
dir.create("workdir2/01.datacleaning/data_quality_before_cleaning/",showWarnings = F,recursive = T)
save(object.mv,file = "workdir2/01.datacleaning/data_quality_before_cleaning/object_impute_mv")
save(object.imput,file = "workdir2/01.datacleaning/data_quality_before_cleaning/object_impute_imput")
##> PCA for raw data
##> 
pca_obj.raw <- IMO_plt_pca(obj = object.imput,tag = "group",center = T,scale = T,removeVar = .1,interactive = F,showImage = F)

##> The batch effect alignment of tidymass needs to provide rt and mz of features of different batches. Since there is no record, SERRF is used for QC-based normalization.
##> SERRF function was simplified.
## load functions:


if (T.norm_method == "serrf") {
  object.norm<- IMOtoolkits::run_serrf(obj = object.imput,
                                         QC_tag = "QC",
                                         cluster_num = 7,
                                         seed = 126)
} else {
  object_svr <-
    object.imput %>% 
    normalize_data(method = T.norm_method)
  
  object.norm <-
    integrate_data(object = object_svr,method = 'qc_mean')
  
}

dir.create("workdir2/02.normalized/")
save(object.norm,file = "workdir2/02.normalized/object.norm.rds")
plt.rsd <- IMO_plt_rsd(obj_old = object.imput,obj_new = object.norm,QC_tag = "QC")
ggsave(filename = "workdir2/02.normalized/rsd.pdf",plot = plt.rsd$plot,width = 10,height = 9.5)

# pca ---------------------------------------------------------------------

pca_new <- IMO_plt_pca(obj = object.norm,tag = color_by,center = T,scale = T,removeVar = .1,interactive = F)
pca_old <- IMO_plt_pca(obj = object.imput,tag = color_by,center = T,scale = T,removeVar = .1,interactive = F)

if(dimension == 2) {
  pca_plt <- pca_old$plot + ggtitle("before")+theme(plot.title = element_text(hjust = .5)) +
    pca_new$plot +ggtitle("after")+theme(plot.title = element_text(hjust = .5))+
    plot_annotation(tag_levels = "A")
  
  ggsave(filename = "workdir2/02.normalized/pca.pdf",plot = pca_plt,width = 10,height = 7)
} else if (dimension == 3) {
  
}


norm_expmat <- object.norm %>% 
  extract_expression_data() %>% 
  rownames_to_column("variable_id")
writexl::write_xlsx(norm_expmat,path = "workdir2/02.normalized/Normalized_expmat.xlsx")

message(msg_run("Finish"))