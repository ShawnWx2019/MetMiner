####################################################################################
#               Prj: MetMiner
#               Assignment: Get kegg annotation information
#               Author: Shawn Wang
#               Date: Oct 30, 2023
####################################################################################

# import packages ---------------------------------------------------------

library(clusterProfiler)
library(tidyverse)
library(tidymass)


# get KEGG cid based on compound name and InChiKey --------------------------

class <- readxl::read_xlsx("../02.Progress_new/final_anno_with_class_new.xlsx")
##> convert compound name to kegg cid by KEGG database
cname2kegg <- MDAtoolkits::mda_name2kegg(query = class |> pull(Compound.name) |> unique(),core_num = 14)
##> convert inchikey to kegg cid
inchi2kegg <- mda_CTS_kegg(query = class |> pull(InChIKey) |> unique() |> str_remove_all(" "),key = "InChIKey",core_num = 14)
inchi2kegg <- 
  inchi2kegg |> 
  filter(KEGG != "") |> 
  mutate(KEGG = str_split(KEGG,",",2,T)[,1])
##> convert compound name to kegg cid by CTS database
ctsname2kegg <- mda_CTS_kegg(query = class |> pull(Compound.name) |> unique(),key = "name",core_num = 14)
ctsname2kegg <- 
  ctsname2kegg |> 
  filter(KEGG != "") |> 
  mutate(KEGG = str_split(KEGG,",",2,T)[,1])
##> Merge KEGG information with raw KEGG CID which compound database provided
class_add_kegg <- 
  class |> left_join(cname2kegg |> filter(KEGG != "") |> setNames(c("Compound.name","KEGG"))) |> 
  mutate(KEGG.ID = case_when(
    is.na(KEGG.ID) ~ KEGG,
    TRUE ~ KEGG.ID
  )) |> dplyr::select(-KEGG)|> 
  left_join(inchi2kegg |> filter(KEGG != "")) |> 
  mutate(KEGG.ID = case_when(
    is.na(KEGG.ID) ~ KEGG,
    TRUE ~ KEGG.ID
  )) |> dplyr::select(-KEGG) |> 
  left_join(ctsname2kegg |> filter(KEGG != "") |> setNames(c("Compound.name","KEGG"))) |> 
  mutate(KEGG.ID = case_when(
    is.na(KEGG.ID) ~ KEGG,
    TRUE ~ KEGG.ID
  )) |> dplyr::select(-KEGG)
##> Export compound annotations with classification and complately KEGG cid
writexl::write_xlsx(class_add_kegg,"../02.Progress_new/annotation_clean5.xlsx")
