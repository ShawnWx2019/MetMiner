####################################################################################
#               Prj: MetMiner
#               Assignment: feature re annotation with tiydmass
#               Author: Shawn Wang
#               Date: May 11, 2023
####################################################################################

# part1 load packages -----------------------------------------------------

library(tidyverse)
library(tidymass)
library(conflicted)
conflict_prefer('filter',"dplyr")
conflict_prefer('select',"dplyr")
column = "rp"
##> load normalized massdatasets
load("../01.Data/obj.serrf.rds")
object <- object.serrf %>% update_mass_dataset()
vari_info <- object %>% 
  extract_variable_info()

exp_mean <- object %>% 
  extract_expression_data() %>% 
  select(contains("QC")) %>% 
  mutate(mean = rowMeans(.)) %>% 
  select(mean) %>% 
  rownames_to_column('variable_id')
checklist <- 
vari_info %>% left_join(exp_mean) %>% 
  mutate(mw = case_when(
    str_detect(variable_id,"neg") ~ mz+1,
    str_detect(variable_id,"pos") ~ mz-1,
  )) %>% 
  select(variable_id,mw,mz,rt,mean) %>% 
  arrange(mw)
library(furrr)
plan(multisession,workers = 3)
round1_filter.mw <- furrr::future_map_dfr(.x = 1:nrow(checklist),.f = function(.x) {
  furrr::future_map_dfr(1:nrow(checklist),.f = function(.y){
    data.frame(
      feature_1 = checklist[.x,1] %>% as.character(),
      feature_2 = checklist[.y,1] %>% as.character(),
      mw_diff = abs(checklist[.x,2] %>% as.numeric() - checklist[.y,2] %>% as.numeric()),
      rt_diff = abs(checklist[.x,4] %>% as.numeric() - checklist[.y,4] %>% as.numeric()),
      mean_diff = checklist[.x,5] %>% as.numeric() - checklist[.y,5] %>% as.numeric()
    )
  })
}) %>% filter(feature_1 != feature_2 & mw_diff < 1 & rt_diff < 0.3,mean_diff > 0)

anno_all <- readxl::read_xlsx("../02.Progress/annotation.xlsx")

feature_redundant = round1_filter.mw %>% 
  select(feature_2 ) %>% 
  distinct() %>% left_join(
    anno_all,by = c('feature_2' = "variable_id")
  ) %>% group_by(feature_2) %>% 
  slice_head(n = 1)

writexl::write_xlsx(feature_redundant,"../02.Progress_new/reduandent.xlsx")


# remove redundant --------------------------------------------------------

object <- 
object %>% 
  activate_mass_dataset("variable_info") %>% 
  left_join(feature_redundant %>% select(feature_2) %>% setNames("variable_id") %>% mutate(tags = 1)) %>% 
  filter(is.na(tags)) %>% 
  select(-tags)


object <- 
  object %>% 
  activate_mass_dataset('variable_info') %>% 
  mutate(rt = rt*60) 
##> generated neg/pos mode.
object_neg <-
  object %>% activate_mass_dataset('variable_info') %>% 
  filter(str_detect(variable_id,"neg"))
object_pos <-
  object %>% activate_mass_dataset('variable_info') %>% 
  filter(str_detect(variable_id,"pos"))

##> add ms2 infomation

object_neg.ms2 <- 
  object_neg %>% 
  mutate_ms2(
    object = .,
    column = 'rp',
    polarity = 'negative',
    ms1.ms2.match.mz.tol = 10,
    ms1.ms2.match.rt.tol = 35,
    path = "../01.Data/NEG/"
  )

object_pos.ms2 <- 
  object_pos %>% 
  mutate_ms2(
    object = .,
    column = 'rp',
    polarity = 'positive',
    ms1.ms2.match.mz.tol = 10,
    ms1.ms2.match.rt.tol = 35,
    path = "../01.Data/POS/"
  )


# part2 feature annotation -----------------------------------------------------------

#> load database
load("~/.HPC_tidymass/MS_db/ReSpect_database2.0.rda")
load("~/.HPC_tidymass/MS_db/PlaSMA_database2.0.rda")
load("~/.HPC_tidymass/MS_db/knapsack_ath_db.rda")
load("~/.HPC_tidymass/MS_db/mona_database0.0.4.rda")
load("~/.HPC_tidymass/MS_db/massbank_database0.0.4.rda")
load("~/.HPC_tidymass/MS_db/ath_plantcyc.database.rda")
load("~/.HPC_tidymass/MS_db/kegg_ms1_database0.0.3.rda")
#> feature annotation

### for positive model ------------------------------------------------------
object_pos_anno.kegg<- 
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 10,
    column = column,
    threads = 7,
    database = kegg_ms1_database0.0.3,
    candidate.num = 10)

anno_kegg.pos <- 
  object_pos_anno.kegg %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M+H)+") 

anno_kegg.pos %>% dplyr::select(Compound.name) %>% unique()

object_pos_anno.ReSpect <- 
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 10,
    column = column,
    threads = 7,
    database = ReSpect_database,
    candidate.num = 10)

anno_respect.pos <- 
  object_pos_anno.ReSpect %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M+H)+") 

object_pos_anno.PlaSMA <- 
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 10,
    column = column,
    threads = 7,
    database = PlaSMA_database,
    candidate.num = 10)

anno_plasma.pos <- 
  object_pos_anno.PlaSMA  %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M+H)+") 

object_pos_anno.MoNA <- 
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 10,
    column = column,
    threads = 7,
    database = mona_database0.0.3,
    candidate.num = 10)

anno_mona.pos <- 
  object_pos_anno.MoNA  %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M+H)+") 

object_pos_anno.massbank <- 
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 10,
    column = column,
    threads = 7,
    database = massbank_database0.0.3,
    candidate.num = 10)

anno_massbank.pos <- 
  object_pos_anno.massbank  %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M+H)+")

object_pos_anno.knapsack <- 
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 10,
    column = column,
    threads = 2,
    database = knapsack_ath_db,
    candidate.num = 10)

anno_knapsack.pos <- 
  object_pos_anno.knapsack  %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M+H)+") 

object_pos_anno.planycyc <- 
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 10,
    column = column,
    threads = 7,
    database = ath_plantcyc.database,
    candidate.num = 10)

anno_plantcyc.pos <- 
  object_pos_anno.planycyc %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M+H)+") 



anno_pos_merge <- 
  bind_rows(anno_kegg.pos,anno_mona.pos,anno_plasma.pos,anno_respect.pos,anno_knapsack.pos,anno_plantcyc.pos) %>% 
  arrange(variable_id)
anno_pos_merge %>% pull(variable_id) %>% unique() %>% length()
# neg ---------------------------------------------------------------------

object_neg_anno.kegg <- 
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 10,
    column = column,
    threads = 7,
    database = kegg_ms1_database0.0.3,
    candidate.num = 10)

anno_kegg.neg <- 
  object_neg_anno.kegg %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M-H)-") 


object_neg_anno.ReSpect <- 
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 10,
    column = column,
    threads = 7,
    database = ReSpect_database,
    candidate.num = 10)

anno_respect.neg <- 
  object_neg_anno.ReSpect %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M-H)-") 


object_neg_anno.PlaSMA <- 
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 10,
    column = column,
    threads = 7,
    database = PlaSMA_database,
    candidate.num = 10)

anno_plasma.neg <-
  object_neg_anno.PlaSMA  %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M-H)-") 

object_neg_anno.MoNA <- 
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 10,
    column = column,
    threads = 7,
    database = mona_database0.0.3,
    candidate.num = 10)

anno_mona.neg<-
  object_neg_anno.MoNA  %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M-H)-") 

object_neg_anno.massbank <- 
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 10,
    column = column,
    threads = 7,
    database = massbank_database0.0.3,
    candidate.num = 10)

anno_massbank.neg <- 
  object_neg_anno.massbank  %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M-H)-") 

object_neg_anno.knapsack <- 
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 10,
    column = column,
    threads = 2,
    database = knapsack_ath_db,
    candidate.num = 10)

anno_knapsack.neg <- 
  object_neg_anno.knapsack  %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M-H)-") 

object_neg_anno.planycyc <- 
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 10,
    column = column,
    threads = 7,
    database = ath_plantcyc.database,
    candidate.num = 10)

anno_plantcyc.neg <- 
  object_neg_anno.planycyc %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M-H)-") 

ms2_plot_mass_dataset(object = object_neg_anno.ReSpect, 
                      variable_id = "Met_neg_1585", polarity = 'negative',
                      database = ReSpect_database, 
                      interactive_plot = FALSE)
ms2_plot_mass_dataset(object = object_pos_anno.massbank, 
                      variable_id = "Met_pos_19", polarity = 'positive',
                      database = massbank_database0.0.3, 
                      interactive_plot = FALSE)

anno_neg_merge <- 
  bind_rows(anno_knapsack.neg,anno_massbank.neg,anno_mona.neg,anno_plasma.neg,anno_respect.neg,anno_plantcyc.neg) %>% 
  arrange(variable_id)
anno_neg_merge %>% pull(variable_id) %>% unique() %>% length()

# merge annotation-------------------------------------------------------------------

anno_merge <- bind_rows(anno_neg_merge,anno_pos_merge) %>% inner_join(
  object %>% extract_variable_info()
)

writexl::write_xlsx(anno_merge,"../02.Progress_new/annotation.xlsx")


# annotation rm duplicated ------------------------------------------------
anno_merge_spec <- 
  anno_merge %>% 
  mutate(Compound.name = str_split(Compound.name,";",2,T)[,1]) %>% 
  group_by(variable_id,Compound.name) %>% 
  filter(Level == min(Level)) %>% 
  filter(Total.score == max(Total.score)) %>% 
  slice_head(n = 1) %>% 
  arrange(variable_id)

##> export for MetDNA
# dir.create("../02.Progress_new/MetDNA/NEG/",showWarnings = F,recursive = T)
# dir.create("../02.Progress_new/MetDNA/POS/",showWarnings = F,recursive = T)
# 
# export_mass_dataset4metdna(object = object_neg.ms2,path = "../02.Progress_new/MetDNA/NEG/")
# 
# export_mass_dataset4metdna(object = object_pos.ms2,path = "../02.Progress_new/MetDNA/POS/")
writexl::write_xlsx(anno_merge_spec,"../02.Progress_new/annotation_spec.xlsx")

anno_merge_rm_reduandant <- 
  anno_merge_spec %>% 
  group_by(variable_id) %>% 
  filter(Level == min(Level)) %>% 
  filter(Total.score == max(Total.score)) %>% 
  slice_head(n = 1)

writexl::write_xlsx(anno_merge_rm_reduandant,"../02.Progress_new/annotation_clean.xlsx")
## export image
save.image(file = "../02.Progress_new/02.Progress/Annotation.rda")
