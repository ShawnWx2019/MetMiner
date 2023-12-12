####################################################################################
#               Prj: MetMiner
#               Assignment: PMhub-ms2 annotation
#               Author: Shawn Wang
#               Date: Oct 11, 2023
####################################################################################

## use raw HPLC-MS mgf file to PMhub online annotation. and generate "new_anno.pos.pmhub.xlsx" and "new_anno.pos.pmhub.xlsx"

# load packages and previous annotation object ----------------------------

load("../02.Progress/Annotation.rda")
library(tidyverse)
library(tidymass)
conflicted::conflict_prefer('select','dplyr')
conflicted::conflict_prefer('filter','dplyr')
anno <- readxl::read_xlsx("../02.Progress_new/annotation_clean2.xlsx",sheet = 2)

# extract variable informations -------------------------------------------
variable_all_anno <- object.serrf |> update_mass_dataset() |> 
  extract_variable_info() |> dplyr::select(1:3)
##> add mz and rt information for each feature
anno_add_level1 <- variable_all_anno |> left_join(anno |> dplyr::select(-mz,-rt)) |> 
  select(variable_id,mz,rt,Compound.name,Lab.ID,Adduct,KEGG.ID) |> 
  mutate(Formula = NA)

anno_add_tidymass <- anno_add_level1 |> filter(!is.na(Compound.name))
##> pick out unannotated features by tidymass
anno_rest <- anno_add_level1 |> dplyr::filter(is.na(Compound.name)) |> 
  select(variable_id,mz,rt)

# match PMhub online annotation -------------------------------------------

pmhub_new <- readxl::read_xlsx("../02.Progress_new/pmHub.pos.xlsx")

new_anno.pos.pmhub <- 
  pos.info |> select(variable_id,mz,rt) |> mutate(mz = round(mz,2)) |> left_join(pmhub_new |> mutate(mz = round(mz,2)),by = 'mz') |> 
  mutate(drt = abs(rt/60 - RT)) |> tibble() |> 
  relocate(drt,.after = RT) |> 
  group_by(variable_id) |> 
  filter(drt == min(drt))

writexl::write_xlsx(new_anno.pos.pmhub,"../02.Progress_new/new_anno.pos.pmhub.xlsx")

neg.info <- object_neg.ms2 |> extract_variable_info()
pmhub_new <- readxl::read_xlsx("../02.Progress_new/pmHub.neg.xlsx")

new_anno.neg.pmhub <- 
  neg.info |> select(variable_id,mz,rt) |> mutate(mz = round(mz,2)) |> left_join(pmhub_new |> mutate(mz = round(mz,2)),by = 'mz') |> 
  mutate(drt = abs(rt/60 - RT)) |> tibble() |> 
  relocate(drt,.after = RT) |> 
  group_by(variable_id) |> 
  filter(drt == min(drt))
writexl::write_xlsx(new_anno.neg.pmhub,"../02.Progress_new/new_anno.neg.pmhub.xlsx")
##> remove compounds not belongs to plant manually, and re-import.
pmhub.pos <- readxl::read_xlsx("../02.Progress_new/new_anno.pos.pmhub.xlsx") 
pmhub.neg <- readxl::read_xlsx("../02.Progress_new/new_anno.neg.pmhub.xlsx")

pmhub.anno <- bind_rows(pmhub.neg,pmhub.pos)
colnames(pmhub.anno)
pmhub.anno <- 
  pmhub.anno |> select(variable_id,pubchemID,libraryMetaboliteID,molecular,source,keggID) |> 
  setNames(
    c('variable_id','Compound.name','Lab.ID','Formula','Database','KEGG.ID')
  )
##> add adduct information
pmhub.anno2 <- inner_join(anno_rest,pmhub.anno,by = 'variable_id') |> 
  mutate(Adduct = case_when(
    str_detect(variable_id,'neg') ~ '(M-H)-',
    str_detect(variable_id,'pos') ~ '(M+H)+')
  ) |> select(-Database)
##> merge final annotation file
anno_all_final <- bind_rows(anno_add_tidymass,pmhub.anno2)

writexl::write_xlsx(anno_all_final,"../02.Progress_new/annotation_clean4.xlsx")
