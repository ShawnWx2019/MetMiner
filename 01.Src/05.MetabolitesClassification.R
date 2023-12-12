####################################################################################
#               Prj: MetMiner
#               Assignment: feature re annotation with tiydmass
#               Author: Shawn Wang
#               Date: Oct 25, 2023
####################################################################################

new_anno <- readxl::read_xlsx("../02.Progress_new/annotation_clean3.xlsx")

# load pacakges -----------------------------------------------------------
library(tidyverse)
library(MDAtoolkits)
library(progressr)
handlers(handler_pbcol(
  adjust = 1.0,
  complete = function(s) cli::bg_red(cli::col_black(s)),
  incomplete = function(s) cli::bg_cyan(cli::col_black(s))
))

# get cid and InChiKey via compound name ----------------------------------

tmp.name_list <- new_anno |> pull(Compound_name) |> unique()
##> get inchikey
tmp.InChIKey <- MDAtoolkits::mda_get_cid_fast(query = tmp.name_list,core_num = 12)
tmp.done <- tmp.InChIKey |> pull(query) |> sort()
##> round2 get inchikey
query2 = setdiff(tmp.name_list |> sort(),tmp.done)
tmp.InChIKey2 <- MDAtoolkits::mda_get_cid_fast(query = query2 ,core_num = 12)
##> merge results
tmp.inchikeyall <- bind_rows(tmp.InChIKey,tmp.InChIKey2)

# class -------------------------------------------------------------------
tmp.inchi.query <- tmp.inchikeyall |> pull(InChIKey)
tmp.class <- cbf_crawler(query = tmp.inchi.query,delay_max = 0.5)
tmp.x = left_join(tmp.inchikeyall,tmp.class,by = 'InChIKey')
##> merge classififcation 
anno_final <- left_join(new_anno,tmp.x |> select(query,InChIKey,superclass,class,subclass,parent_levels,description),by = c('Compound_name'= 'query'))
writexl::write_xlsx(anno_final,"../02.Progress_new/anno_class.final.xlsx")

##################################################
##>             Update Oct 30,2023
##>     Pmhub online annotation classification    
##################################################
library(tidyverse)
library(tidymass)
library(MDAtoolkits)
##> import pmhub annotations
anno_tbl <- readxl::read_xlsx("../02.Progress_new/annotation_clean4.xlsx")
pmhub.smile <- readxl::read_xlsx("../02.Progress_new/pmhub_smile.xlsx") |> 
  setNames(c('Lab.ID','smiles','Compound.name'))
## clean compound name
anno_tbl.clean <- 
  anno_tbl |> 
  mutate(
    Compound.name = str_remove_all(Compound.name,"\\?|\\'") 
  )|> 
  mutate(
    Compound.name = str_remove_all(Compound.name,'\\"') 
  )
pmhub.smile.clean <-
  pmhub.smile |> 
  mutate(
    Compound.name = str_remove_all(Compound.name,"\\?|\\'") 
  )|> 
  mutate(
    Compound.name = str_remove_all(Compound.name,'\\"') 
  )

# round 1 get inchikey based on compound name -----------------------------
## compound name as key
query_name = anno_tbl.clean |> pull(Compound.name) |> unique()

name2Inchi = mda_get_cid_fast(query = query_name,core_num = '14')

# round 2 get inchikey based on SMILE -------------------------------------
query_rest <- anti_join(data.frame(query = query_name),name2Inchi |> dplyr::select(query)) |> 
  setNames('Compound.name') |> 
  inner_join(pmhub.smile.clean) 
rest_smile2inchi = mda_get_cid_fast(query = (query_rest |> pull(smiles)),core_num = 14,probe = 'smiles')
## merge results
name2Inchi2 <- inner_join(rest_smile2inchi,query_rest,by = c("query" = "smiles")) |> 
  dplyr::select(Compound.name,cid,formula,mw,InChIKey)
name2Inchi.final = name2Inchi |> setNames(colnames(name2Inchi2)) |> 
  rbind(name2Inchi2)
##> get compound classification by inchikey
compound_class <- MDAtoolkits::cbf_crawler(query = name2Inchi.final |> pull(InChIKey),delay_max = 0.5)
rest_inchi_class <-compound_class |> filter(is.na(kingdom)) |> pull(InChIKey)
compound_class2 <- MDAtoolkits::cbf_crawler(query = rest_inchi_class,delay_max = 0.5 )
writexl::write_xlsx(compound_class2,"../02.Progress_new/failed.xlsx")
##> merge data
final_anno_with_class <- left_join(anno_tbl.clean,name2Inchi.final) |> left_join(compound_class)

writexl::write_xlsx(final_anno_with_class,"../02.Progress_new/final_anno_with_class.xlsx")

# Get KEGG annotation ---------------------------------------------------------------------
##> add InChiKey information manually by other database.
final_anno_with_class_man <- readxl::read_xlsx("../02.Progress_new/final_anno_with_class.xlsx")
##> new inchikeys
new_inch_query <- final_anno_with_class_man |> 
  filter(is.na(kingdom)) |> pull(InChIKey) |> unique() |> 
  str_remove_all(" ")
##> new classification information
class_new <- MDAtoolkits::cbf_crawler(query = new_inch_query,delay_max = 0.5 )

class_old <- final_anno_with_class_man |> dplyr::select(colnames(class_new)) |> 
  filter(!is.na(kingdom))
class_all <- rbind(class_new,class_old)
##> merge all data
final_anno_with_class_new <- 
  final_anno_with_class_man |> dplyr::select(-c(13:18)) |> 
  left_join(class_all) |> 
  distinct()
writexl::write_xlsx(final_anno_with_class_new,"../02.Progress_new/final_anno_with_class_new.xlsx")

