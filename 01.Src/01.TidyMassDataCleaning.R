####################################################################################
#               Prj: MetMiner
#               Assignment: feature re annotation with tiydmass
#               Author: Shawn Wang
#               Date: Apr 11, 2022
####################################################################################
library(tidyverse)
library(tidymass)
library(MDAtoolkits)
library(PCAtools)
library(ComplexHeatmap)
library(progressr)
select = dplyr::select
# importdata --------------------------------------------------------------
handlers(handler_pbcol(
  adjust = 1.0,
  complete = function(s) cli::bg_red(cli::col_black(s)),
  incomplete = function(s) cli::bg_cyan(cli::col_black(s))
))

raw_exp_mat <- read.csv(
  "../01.Data/exp_mat_v2.csv"
)

expression_data <- 
  raw_exp_mat %>% 
  select(-rt,-mz) %>% 
  column_to_rownames("variable_id")

variable_info <- 
  raw_exp_mat %>% 
  select(variable_id:rt)

sample_info <- read.csv(
  "../01.Data/sample_info_v2.csv"
)

## 1.2 => construct a massdataset by tidymass.

object <- 
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

## 1.3 => save rawdata object

dir.create("../02.Progress_new/01.datacleaning",showWarnings = F,recursive = T)

save(object,file = "../02.Progress_new/01.datacleaning/object_raw.rds")

# step2.  data cleaning ---------------------------------------------------

## 2.1 => quality validation
dir.create("../02.Progress_new/01.datacleaning/data_quality_before_cleaning",showWarnings = F,recursive = T)
massqc_report(
  object = object,
  path = "../02.Progress_new/01.datacleaning/data_quality_before_cleaning/"
)


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
  dplyr::filter(group == "QC") %>% 
  pull(sample_id)

wt_id <- 
  object %>% 
  activate_mass_dataset(
    what = "sample_info"
  ) %>% 
  dplyr::filter(group == "WT") %>% 
  pull(sample_id)

mt_id <- 
  object %>% 
  activate_mass_dataset(
    what = "sample_info"
  ) %>% 
  dplyr::filter(group == "MT") %>% 
  pull(sample_id)


## 2.2.2 Calculate missing frequency.
object <- 
  object %>% 
  mutate_variable_na_freq(according_to_samples = qc_id) %>% 
  mutate_variable_na_freq(according_to_samples = wt_id) %>% 
  mutate_variable_na_freq(according_to_samples = mt_id)
## 2.2.3 remove noisy feature
object.mv <- 
  object %>% 
  activate_mass_dataset(
    what = "variable_info"
  ) %>% 
  filter(na_freq < 0.2 & na_freq.1 < 0.2 & na_freq.2 < 0.5)

outlier_samples <- 
  object %>% 
  `+`(1) %>% 
  log(2) %>% 
  scale() %>% 
  detect_outlier()

outlier_table <-
  extract_outlier_table(outlier_samples) 

outlier_table %>% 
  apply(1, function(x){
    sum(x)
  }) %>% 
  `>`(0) %>% 
  which()
##> QC samples were concentrated for 5x, so the outliers were ignored

## 2.4 => missing value imputation 

object.imput <- 
  impute_mv(
    object = object.mv,
    method = "knn"
  )
get_mv_number(object.imput)

save(object.mv,file = "../02.Progress_new/01.datacleaning/data_quality_before_cleaning/object_impute_mv")
save(object.imput,file = "../02.Progress_new/01.datacleaning/data_quality_before_cleaning/object_impute_imput")
##> PCA for raw data
##> 
library(IMOtoolkits)
pca_obj.raw <- IMO_plt_pca(obj = object.imput,tag = "group",center = T,scale = T,removeVar = .1,interactive = F,showImage = F)

##> The batch effect alignment of tidymass needs to provide rt and mz of features of different batches. Since there is no record, SERRF is used for QC-based normalization.
##> SERRF function was simplified.
## load functions:
# source('~/My_Repo/IMOtoolkits/R/run_serrf.R')

object_svr <-
  object.imput %>% 
  normalize_data(method = 'svr')

object_svr_inter <-
  integrate_data(object = object_svr,method = 'qc_mean')

rsd.obj = IMO_plt_rsd(
  obj_old = object.imput,obj_new = object_svr_inter,QC_tag = "QC"
)
IMO_plt_pca(obj = object_svr_inter,tag = "group",center = T,scale = T,removeVar = .1,interactive = T)


object.serrf <- IMOtoolkits::run_serrf(obj = object.imput,
                                       QC_tag = "QC",
                                       cluster_num = 7,
                                       seed = 126)



IMO_plt_rsd(obj_old = object.imput,obj_new = object.serrf,QC_tag = "QC")
IMO_plt_pca(obj = object.serrf,tag = "group",center = T,scale = T,removeVar = .1,interactive = T)
dir.create("../02.Progress_new/02.normalized/")
save(object.serrf,file = "../02.Progress_new/02.normalized/obj.serrf.rds")
