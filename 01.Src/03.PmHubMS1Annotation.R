####################################################################################
#               Prj: MetMiner
#               Assignment: PMhub-ms1 annotation
#               Author: Shawn Wang
#               Date: Oct 11, 2023
####################################################################################
####################################################################################
##>                     Update at Nov 12,2023.
##>
##>              WANNING:  PMhub annotation were processed online, 
##>          so this part was removed from final annotation results!
##>
####################################################################################

##> Run PMhub ms1 annotation by tidymass.
library(tidyverse)
library(tidymass)

# load rawannotation ------------------------------------------------------
anno_table <- readxl::read_xlsx("../02.Progress_new/Final_annotation.xlsx")
##> mark unannotated features with MW and RT
anno_table <-
  anno_table %>% 
  mutate(
    Compound_name = case_when(
      is.na(Compound_name) ~ paste0("MW:",round(mw,3),"@RT:",round(RT,3)),
      TRUE  ~ Compound_name
    )
  ) %>% 
  mutate()

## pmhub ms1 annotation database.

load("~/.HPC_tidymass/MS_db/pmhub.rda")

object_pos_anno.pmhub<- 
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 10,mz.ppm.thr = 20,
    column = column,
    threads = 7,
    database = pmhub_db_v0.01,
    candidate.num = 10)

anno_pmhub.pos <- 
  object_pos_anno.pmhub %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M+H)+") 
anno_pmhub.pos %>% arrange(Compound.name) %>% tibble()

object_neg_anno.pmhub<- 
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 10,mz.ppm.thr = 20,
    column = column,
    threads = 7,
    database = pmhub_db_v0.01,
    candidate.num = 10)

anno_pmhub.neg <- 
  object_neg_anno.pmhub %>% 
  extract_annotation_table() %>% 
  filter(Adduct == "(M-H)-") 
pmhub_tbl <- rbind(anno_pmhub.neg,anno_pmhub.pos)
writexl::write_xlsx(pmhub_tbl,"../02.Progress_new/pmhub.xlsx")
pmhub_clean <- 
pmhub_tbl %>% 
  group_by(variable_id) %>% 
  filter(Total.score == max(Total.score)) %>% 
  slice_head(n = 1) %>% 
  arrange(Compound.name)
writexl::write_xlsx(pmhub_clean,"../02.Progress_new/pmhub_clean.xlsx")

anno_pmhub_tbl <- 
  pmhub_clean %>% 
  select(
    variable_id,Compound.name,Adduct,Total.score,Lab.ID,KEGG.ID
  ) %>% 
  setNames(c("Compound_id","Compound_name","Adduct","Total.score","Lab.ID.iso","KEGG.ID.iso")) %>% 
  left_join(anno_table %>% select(Compound_id,mw,RT,CD)) %>% 
  select(colnames(anno_table))

anno_table2 <- rbind(anno_table,anno_pmhub_tbl)

anno_table2 %>% group_by(Compound_id) %>% 
  arrange(Compound_id) %>% 
  filter(!is.na(Compound_name))

writexl::write_xlsx(anno_table2,"../02.Progress_new/annotation_clean2.xlsx")
  
anno_table2 <- readxl::read_xlsx("../02.Progress_new/annotation_clean2.xlsx")

anno_table3 <- 
anno_table2 %>% 
  group_by(Compound_id) %>% 
  slice_min(level) %>% 
  filter(Total.score == max(Total.score)) %>% 
  slice_head()

writexl::write_xlsx(anno_table3,"../02.Progress_new/annotation_clean3.xlsx")sz
