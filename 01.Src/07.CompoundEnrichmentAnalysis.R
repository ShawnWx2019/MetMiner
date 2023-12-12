####################################################################################
#               Prj: MetMiner
#               Assignment: Classification and KEGG enrichment
#               Author: Shawn Wang
#               Date: Oct 30, 2023
####################################################################################

# import library ----------------------------------------------------------

library(clusterProfiler)
library(RCurl)
library(XML)
library(tidyverse)
library(tidymass)

# load data-------------------------------------------------------------------------
##> module 2 compounds
C2M <- read.delim("../02.Progress_new/WGCNA/round4/01.Gene2Module.xls",header = T,sep = "\t")

##> clean network file for gephi

# path = "../03.progress/WGCNA/round4/"
# 
# edge_file <- dir(path = path,pattern = "04.*.edge.xls")
# 
# edge <- purrr::map_dfr(
#   .x = edge_file,.f = function(.x) {
#     a = read.delim(paste0(path,.x),header = T,sep = "\t") %>% 
#       select(fromNode,	toNode,	weight) %>% 
#       setNames(c("source","target",'weight'))
#   }
# ) 
# edge.new <- 
#   edge %>% inner_join(node %>% select(source)) %>% inner_join(node %>% select(source) %>% setNames("target"))
# 


# convert classification to t2g and t2n -----------------------------------
##> To run enrichment ananlysis, classification data MUST be convert as TERM2GENE (Classifiction term to metabolites) and TERM2NAME (Classification Term to desciription)

conflict_prefer(name = "filter",winner = "dplyr")
conflict_prefer(name = "select",winner = "dplyr")
##> import classification table
class <- readxl::read_xlsx("../02.Progress_new/final_anno_with_class_new.xlsx")

## reformat, remove unclassified metabolites and convert wider table to long.
class2 <- 
  class %>% 
  filter(!is.na(superclass)) %>% 
  select(variable_id,Compound.name,superclass,class,subclass,parent_levels) %>% 
  separate_longer_delim(cols = parent_levels,delim = " | ") %>% ## seprate parent levels 
  pivot_longer(!all_of(c("variable_id","Compound.name")),names_to = "levels",values_to = "NAME") %>% 
  filter(!is.na(NAME))


TERM2NAME = 
  class2 %>% 
  select(NAME) %>% 
  distinct() %>% 
  mutate(TERM = paste0("MC:",str_pad(c(1:nrow(.)),5,"left",'0'))) %>% ## add terms
  filter(NAME != "NA") %>% 
  select(TERM,NAME)

TERM2COM = 
  class2 %>% 
  inner_join(TERM2NAME) %>% 
  select(TERM,variable_id) %>% 
  setNames(c("TERM","COMPOUND")) %>% 
  distinct()

writexl::write_xlsx(list(
  db = class2,
  TERM2NAME = TERM2NAME,
  TERM2COM = TERM2COM
),path = "../01.Data/classyfire_database.xlsx")


# summary of classification results ----------------------------------------------------------
##> count compound numbers of each term
class3 <-
  class %>% 
  select(-Compound.name,-Lab.ID,-description,-InChIKey) %>% 
  filter(!superclass == "NA") %>% 
  mutate(
    class = dplyr::case_when(
      is.na(class) | class == "NA" ~ "not classified",
      TRUE ~ class
    ),
    subclass = dplyr::case_when(
      is.na(subclass) | subclass == "NA" ~ "not classified",
      TRUE ~ subclass
    )
  ) %>% group_by(superclass,class) %>% 
  summarise(n = n())


# Classification enrichment for each module -------------------------------

conflict_prefer("dotplot","clusterProfiler")
res.cate.black<- enricher(
  gene = C2M %>% dplyr::filter(Module == "black") %>% pull(GID),
  TERM2GENE = TERM2COM ,TERM2NAME = TERM2NAME,pvalueCutoff = 0.5
)

dotplot(res.cate.black,showCategory = 10)

res.cate.blue<- enricher(
  gene = C2M %>% dplyr::filter(Module == "blue") %>% pull(GID),
  TERM2GENE = TERM2COM,TERM2NAME = TERM2NAME,pvalueCutoff = 0.5
)

dotplot(res.cate.blue,showCategory = 10)

res.cate.brown <- enricher(
  gene = C2M %>% dplyr::filter(Module == "brown") %>% pull(GID),
  TERM2GENE = TERM2COM,TERM2NAME = TERM2NAME,pvalueCutoff = 0.5
)

dotplot(res.cate.brown,showCategory = 10)

res.cate.green<- enricher(
  gene = C2M %>% dplyr::filter(Module == "green") %>% pull(GID),
  TERM2GENE = TERM2COM,TERM2NAME = TERM2NAME,pvalueCutoff = 0.5
)

dotplot(res.cate.green,showCategory = 10)

res.cate.greenyellow<- enricher(
  gene = C2M %>% dplyr::filter(Module == "greenyellow") %>% pull(GID),
  TERM2GENE = TERM2COM,TERM2NAME = TERM2NAME,pvalueCutoff = 0.5
)

dotplot(res.cate.greenyellow,showCategory = 10)

res.cate.magenta<- enricher(
  gene = C2M %>% dplyr::filter(Module == "magenta") %>% pull(GID),
  TERM2GENE = TERM2COM,TERM2NAME = TERM2NAME,pvalueCutoff = 0.5
)

dotplot(res.cate.magenta,showCategory = 10)

res.cate.pink<- enricher(
  gene = C2M %>% dplyr::filter(Module == "pink") %>% pull(GID),
  TERM2GENE = TERM2COM,TERM2NAME = TERM2NAME,pvalueCutoff = 0.5
)

dotplot(res.cate.pink,showCategory = 10)

res.cate.purple<- enricher(
  gene = C2M %>% dplyr::filter(Module == "purple") %>% pull(GID),
  TERM2GENE = TERM2COM,TERM2NAME = TERM2NAME,pvalueCutoff = 0.5
)

dotplot(res.cate.purple,showCategory = 10)

res.cate.red <- enricher(
  gene = C2M %>% dplyr::filter(Module == "red") %>% pull(GID),
  TERM2GENE = TERM2COM,TERM2NAME = TERM2NAME,pvalueCutoff = 0.5
)

dotplot(res.cate.red,showCategory = 10)

res.cate.turquoise<- enricher(
  gene = C2M %>% dplyr::filter(Module == "turquoise") %>% pull(GID),
  TERM2GENE = TERM2COM,TERM2NAME = TERM2NAME,pvalueCutoff = 0.5
)

dotplot(res.cate.turquoise,showCategory = 10)

res.cate.tan<- enricher(
  gene = C2M %>% dplyr::filter(Module == "tan") %>% pull(GID),
  TERM2GENE = TERM2COM,TERM2NAME = TERM2NAME,pvalueCutoff = 0.5
)

dotplot(res.cate.tan,showCategory = 10)

res.cate.yellow<- enricher(
  gene = C2M %>% dplyr::filter(Module == "yellow") %>% pull(GID),
  TERM2GENE = TERM2COM,TERM2NAME = TERM2NAME,pvalueCutoff = 0.5
)

dotplot(res.cate.yellow,showCategory = 10)


result_all <- clusterProfiler::merge_result(
  enrichResultList = list(
    black = res.cate.black,
    blue = res.cate.blue,
    brown = res.cate.brown,
    green = res.cate.green,
    greenyellow = res.cate.greenyellow,
    magenta = res.cate.magenta,
    pink = res.cate.pink,
    purple = res.cate.purple,
    red = res.cate.red,
    tan = res.cate.tan,
    turquoise = res.cate.turquoise,
    yellow = res.cate.yellow
  )
)

writexl::write_xlsx(as.data.frame(result_all),"../02.Progress_new/classificationEnrich.xlsx")
enrich_class <- 
  dotplot(result_all,showCategory = 15,label_format = 100) + 
  coord_flip() +
  xlab("Module")+
  guides(size = guide_legend(title = "Metabolite\nRatio"))+
  ggtitle("Classification enrichment")+
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(hjust = .5)
  )
ggsave(plot = enrich_class,filename = "../02.Progress_new/ClassificationEnrichment.pdf",
       width = 13,height = 8)

# treemap for class -----------------------------------------------------------------
class_treemap_data <- 
  class |> 
  filter(superclass != 'NA') |> 
  select(superclass,class,subclass,parent_levels) |> 
  group_by(superclass) |> 
  mutate(num = n()) |> 
  separate_longer_delim(cols = parent_levels,delim = " | ") |> 
  mutate(class = case_when(
    class == "NA" ~ NA,
    TRUE ~ class
  ),subclass = case_when(
    subclass == "NA" ~ NA,
    TRUE ~ subclass
  ),parent_levels = case_when(
    parent_levels == "NA" ~ NA,
    TRUE ~ parent_levels
  )) |> 
  arrange(num)
pdf("../02.Progress_new/class_treemap.pdf",width = 12,height = 14)
treemap::treemap(
  class_treemap_data,
  index = c("superclass","class","subclass","parent_levels"),
  vSize = "num",
  vColor = 'num',
  type = "dens",
  palette = "RdYlBu", 
  fontsize.labels=c(16, 10),
  align.labels=list(c("center", "center"), c("left", "top")),
  title = 'ClassyFire',
)
dev.off()

load("../01.Data/obj.serrf.rds")
variable_info_final <- 
  object.serrf |> extract_variable_info() |> 
  left_join(class |> select(variable_id,Compound.name)) |> 
  mutate(Compound.name = case_when(
    is.na(Compound.name) ~ paste0("MW:",round(mz,3),"@RT:",round(rt,3)),
    TRUE ~ Compound.name
  ))

writexl::write_xlsx(variable_info_final,"../02.Progress_new/variable_info_final.xlsx")

# for kegg ---------------------------------------------------------------
##> KEGG plant database provided by TBtools
plantKeggback <- read.delim("../02.data/TBtools.Plant.KEGG.backEnd.20200724",header = F,quote = "")

t2n <- 
  plantKeggback %>% 
  select(V1,V3) %>% 
  mutate(term = str_extract(V3,"\\d++(?= )"),
         name = str_remove(V3,"\\d++(?= )")) %>% 
  rename("gene" = "V1") %>% 
  select(term,name) %>% 
  distinct()

t2n2 <- 
  plantKeggback %>% 
  select(V1,V4) %>% 
  mutate(V4 = str_remove(V4,"B  ")) %>% 
  mutate(term = str_extract(V4,"\\d++(?= )"),
         name = str_remove(V4,"\\d++(?= )")) %>% 
  rename("gene" = "V1") %>% 
  select(term,name) %>% 
  distinct()
##> remove useless terms
plant_kegg <- rbind(t2n,t2n2) %>% 
  mutate(term = paste0("map",term)) %>% 
  filter(!str_detect(name,"Other|Unclassified|Drug|Bacterial|Not")) %>% 
  arrange(name)
## Get cid via kegg mapid
plant_kegg2name <- map_dfr(.x = plant_kegg$term,.f = function(.x){
  url = paste0("https://rest.kegg.jp/link/cpd/",.x)
  tbl2 = data.frame(
    TERM = .x,
    COMPOUND = ""
  )
  Sys.sleep(sample(runif(n = 50, 1, 1.5) %>% round(., digits = 2), 1))
  tryCatch(
    {
      xx2 <- getURL(url) 
      if(xx2 == "\n") {
        tbl2 = tbl2
      }else {
        tbl2 = xx2 %>% 
          read.table(text = .,sep = "\t") %>% 
          setNames(c("TERM","COMPOUND")) %>% 
          mutate(
            TERM = str_remove(TERM,"path\\:"),
            COMPOUND = str_remove(COMPOUND,"cpd\\:"))
      }
    }
  )
  return(tbl2)
})

kegg.match <- class |> select(variable_id,Compound.name,KEGG.ID) |> drop_na() 
plant_kegg2name2 = plant_kegg2name %>% dplyr::filter(!COMPOUND == "")
##》 retain kegg cids match with plant compounds
kegg2gene = inner_join(
  plant_kegg2name,kegg.match,by =c("COMPOUND" = "KEGG.ID")
) %>% distinct()
x = 
  inner_join(plant_kegg,kegg2gene,by = c('term'='TERM')) |> 
  distinct()
writexl::write_xlsx(x,"../02.Progress_new/KEGG_anno_with_map.xlsx")
k2c.plant_nb <- kegg2gene %>% 
  dplyr::select(TERM,variable_id) |> 
  k2n.plant_db <- inner_join(
    k2c.plant_nb,plant_kegg,by = c("TERM" = 'term')
  ) %>% dplyr::select(TERM,name) %>% distinct()

k2c.plant_nb = plant_kegg2name %>% 
  inner_join(kegg.match %>% dplyr::select(variable_id,KEGG.ID) %>% setNames(c("GENE","KEGGID"))) %>% 
  select(PathwayID,GENE) %>% 
  distinct()

# kegg pathway analysis ---------------------------------------------------

kegg_all <- enricher(
  gene = k2c.plant_nb |> pull(variable_id) |> unique(),
  TERM2GENE = k2c.plant_nb,TERM2NAME = k2n.plant_db,pvalueCutoff = 1,qvalueCutoff = 1
)
kegg_all.tbl <- 
  kegg_all |> as.data.frame() |> arrange(desc(Count)) |> 
  mutate(Description = factor(Description,levels = Description))
kegg_summary <-  
  ggplot(data = kegg_all.tbl,mapping = aes(x = Count,y = Description,fill = Count))+
  geom_col(color = 'black') +
  xlab('Metabolite number') + 
  ylab('')+
  scale_fill_viridis_c(alpha = .7)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 12,colour = 'black'),
    axis.ticks = element_line(size = 1,color = 'black'),
    panel.border = element_rect(linewidth = 1,colour = "black"),
    legend.position = ""
  )

ggsave(kegg_summary,filename = "../02.Progress_new/Figure4.KEGGsummay.pdf",width = 10,height = 8)


# kegg enrichment for each wgcna module -----------------------------------

res.KEGG.black<- enricher(
  gene = C2M %>% dplyr::filter(Module == "black") %>% pull(GID),
  TERM2GENE = k2c.plant_nb,TERM2NAME = k2n.plant_db,pvalueCutoff = 1
)

dotplot(res.KEGG.black,showCategory = 10)

res.KEGG.blue<- enricher(
  gene = C2M %>% dplyr::filter(Module == "blue") %>% pull(GID),
  TERM2GENE = k2c.plant_nb,TERM2NAME = k2n.plant_db,pvalueCutoff = 1
)

dotplot(res.KEGG.blue,showCategory = 10)

res.KEGG.brown <- enricher(
  gene = C2M %>% dplyr::filter(Module == "brown") %>% pull(GID),
  TERM2GENE = k2c.plant_nb,TERM2NAME = k2n.plant_db,pvalueCutoff = 1
)

dotplot(res.KEGG.brown,showCategory = 10)

res.KEGG.green<- enricher(
  gene = C2M %>% dplyr::filter(Module == "green") %>% pull(GID),
  TERM2GENE = k2c.plant_nb,TERM2NAME = k2n.plant_db,pvalueCutoff = 1
)

dotplot(res.KEGG.green,showCategory = 10)

res.KEGG.greenyellow<- enricher(
  gene = C2M %>% dplyr::filter(Module == "greenyellow") %>% pull(GID),
  TERM2GENE = k2c.plant_nb,TERM2NAME = k2n.plant_db,pvalueCutoff = 1
)

dotplot(res.KEGG.greenyellow,showCategory = 10)

res.KEGG.magenta<- enricher(
  gene = C2M %>% dplyr::filter(Module == "magenta") %>% pull(GID),
  TERM2GENE = k2c.plant_nb,TERM2NAME = k2n.plant_db,pvalueCutoff = 1
)

dotplot(res.KEGG.magenta,showCategory = 10)

res.KEGG.pink<- enricher(
  gene = C2M %>% dplyr::filter(Module == "pink") %>% pull(GID),
  TERM2GENE = k2c.plant_nb,TERM2NAME = k2n.plant_db,pvalueCutoff = 1
)

dotplot(res.KEGG.pink,showCategory = 10)

res.KEGG.purple<- enricher(
  gene = C2M %>% dplyr::filter(Module == "purple") %>% pull(GID),
  TERM2GENE = k2c.plant_nb,TERM2NAME = k2n.plant_db,pvalueCutoff = 1
)

dotplot(res.KEGG.purple,showCategory = 10)

res.KEGG.red <- enricher(
  gene = C2M %>% dplyr::filter(Module == "red") %>% pull(GID),
  TERM2GENE = k2c.plant_nb,TERM2NAME = k2n.plant_db,pvalueCutoff = 1
)

dotplot(res.KEGG.red,showCategory = 10)

res.KEGG.tan <- enricher(
  gene = C2M %>% dplyr::filter(Module == "tan") %>% pull(GID),
  TERM2GENE = k2c.plant_nb,TERM2NAME = k2n.plant_db,pvalueCutoff = 1
)

dotplot(res.KEGG.tan,showCategory = 10)

res.KEGG.turquoise<- enricher(
  gene = C2M %>% dplyr::filter(Module == "turquoise") %>% pull(GID),
  TERM2GENE = k2c.plant_nb,TERM2NAME = k2n.plant_db,pvalueCutoff = 1
)

dotplot(res.KEGG.turquoise,showCategory = 10)




result_all <- clusterProfiler::merge_result(
  enrichResultList = list(
    black = res.KEGG.black,
    blue = res.KEGG.blue,
    brown = res.KEGG.brown,
    green = res.KEGG.green,
    greenyellow = res.KEGG.greenyellow,
    magenta = res.KEGG.magenta,
    pink = res.KEGG.pink,
    purple = res.KEGG.purple,
    red = res.KEGG.red,
    tan = res.KEGG.tan,
    turquoise = res.KEGG.turquoise
  )
)
writexl::write_xlsx(result_all %>% as.data.frame(),"../03.progress/05.Network/KEGGenrichment.xlsx")
enrich_class <- 
  dotplot(result_all,showCategory = 100) + 
  coord_flip() +
  xlab("Module")+
  guides(size = guide_legend(title = "Metabolite\nRatio"))+
  ggtitle("Classification enrichment")+
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(hjust = .5)
  )
ggsave(plot = enrich_class,filename = "../03.progress/05.Network/KEGGenrichment.pdf",
       width = 13,height = 8)
