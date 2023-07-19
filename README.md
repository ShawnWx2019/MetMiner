# Users Manual of MetMiner

[![R version](https://img.shields.io/badge/R-v4.1.1-salmon)](https://www.r-project.org) ![lifecycle](https://img.shields.io/badge/lifecycle-Experimental-lightcyan) [![license](https://img.shields.io/badge/license-MIT-red)](https://opensource.org/licenses/MIT) [![Myblog](https://img.shields.io/badge/Blog-ShanwLearnBioinfo-purple)](https://shawnwx2019.github.io/)

An Integrated Pipeline for Large-Scale Metabolomics Data Processing and Data Mining.

<center>
<img src="https://github.com/ShawnWx2019/MetMiner/tree/main/04.www/Figure1.png" alt="Image" style="max-width:80%">
</center>

# Getting started

[English](https://github.com/ShawnWx2019/MetMiner/blob/main/README.md) \| [中文](https://github.com/ShawnWx2019/MetMiner/blob/main/README.CN.md)

------------------------------------------------------------------------

R version: `>4.1.1`

OS: `MacOS > 10.10`, `Win 7-11`, `Ubuntu 20.04`

*Only Ubuntu 20.04 has passed the test. Other Linux distributions need to be tested.*

# Dependence

LC-MS data analysis framwork: [**TidyMass**](https://www.tidymass.org/) developed by [Dr. Xiaotao Shen](https://www.shenxt.info/). [@Citation](#refer-anchor-1)

```r
if(!require(remotes)){
install.packages("remotes")
}
remotes::install_gitlab("tidymass/tidymass")
```

Metabolomics Downstream Analysis toolkits: [**MDAtoolkits**](https://github.com/ShawnWx2019/MDAtoolkits/tree/master) 

```r
## install from github
suppressMessages(if (!require('MDAtoolkits')) install_github(repo = "ShawnWx2019/MDAtoolkits",ref = 'master'))
## two functions from another package called IMOtoolkits. will be intergreted with MDAtoolkits.
suppressMessages(if (!require('IMOtoolkits')) install_github(repo = "ShawnWx2019/IMOtoolkits"))
```

Untargeted metabolomics upstream analysis pipeline based on tidyMass: [**HPC-tidymass** ](https://github.com/ShawnWx2019/HPC_tidymass) 

```r
##>  you can download this pipeline from github or just clone the repo to your server
##>  1.0 clone the repo to your server
git clone https://github.com/ShawnWx2019/HPC_tidymass.git

##> Run the initialization script to configure the runtime environment.

cd HPC_tidymass && chmod +x init.sh && bash init.sh

##> test
hpc-runTidymass -h
```

A user-friendly WGCNA Shiny app: [**WGCNA-shinyApp** ](https://github.com/ShawnWx2019/WGCNA-shinyApp) follow the steps of [WGCNA-shinyapp readme file](https://github.com/ShawnWx2019/WGCNA-shinyApp/blob/main/README.md)



# Step 1. Upstream analysis

Here we refer to the process of obtaining the metabolomics data from LC-MS raw data to the generation of a metabolite accumulation matrix as <font color=blue>**「upstream data analysis data」**</font>. It includes processes such as <font color=green> 「raw data format conversion」, 「peak picking」, 「data cleaning」, and 「data normalization」</font>. In this pipeline, this part is typically performed by TidyMass or Compound Discoverer. If you have used other software for peak picking, you can refer to https://www.tidymass.org/start/create_mass_dataset/ to create a mass_dataset, which can then be used for subsequent analysis with TidyMass.

## 1.1 Pseudotargeted Metabolomics

### 1.1.1 Generate TQMS Method.

For pseudotargeted metabolomics, in the initial stage, HPLC experiments can be conducted using 3-10 QC (quality control) samples. Since the data volume is relatively small, and considering the completeness of commercial software databases, we recommend utilizing dedicated analysis software for peak detection and library searching. Alternatively, you can also opt to perform upstream analysis directly with TidyMass, followed by MRM (Multiple Reaction Monitoring) selection to generate a TQMS (Triple Quadrupole Mass Spectrometry) method.

**Generate TQMS method from raw data**


If using commercial software, it is recommended to follow the instructions provided in the software's manual for step-by-step analysis. If using TidyMass, it is advised to refer to the tutorial available at https://www.tidymass.org/start/ and follow the provided guidelines. Alternatively, you can also execute the [PresMetaboUpAnalysis.R](https://github.com/ShawnWx2019/MetMiner/blob/main/01.Src/PresMetaboUpAnalysis.R) provided by the pipeline.

1. Convert `.raw` data to `.mzXML` and `.mgf` by [MSCovert](https://proteowizard.sourceforge.io/download.html)

2. prepare your files according to the file storage locations shown below.

```dir
02.DemoData
├── MS1
│   ├── NEG
│   │   └── QC
│   │       ├── QC_01.mzXML
│   │       ├── QC_02.mzXML
│   │       ├── QC_03.mzXML
│   │       └── QC_04.mzXML
│   └── POS
│       └── QC
│           ├── QC_01.mzXML
│           ├── QC_02.mzXML
│           ├── QC_03.mzXML
│           └── QC_04.mzXML
├── MS2
│   ├── NEG
│   │   ├── QC_01.mgf
│   │   ├── QC_02.mgf
│   │   ├── QC_03.mgf
│   │   └── QC_04.mgf
│   └── POS
│       ├── QC_01.mgf
│       ├── QC_02.mgf
│       ├── QC_03.mgf
│       └── QC_04.mgf
├── sample_info.csv
```

3. run `PreMetaboUpAnalysis.R` under terminal or in Rstudio.

    3.1 Note that you should add or remove database according to your own needs, starting from line 480 of the script. Also, make sure to modify the location where the database file is stored.
    
    ```r
    ##> load database
    load("~/.HPC_tidymass/MS_db/mona_database0.0.3.rda")
    load("~/.HPC_tidymass/MS_db/RIKEN_PlaSMA_database0.0.1.rda")
    load("~/.HPC_tidymass/MS_db/knapsack_ath_db.rda")
    load("~/.HPC_tidymass/MS_db/kegg_ms1_database0.0.3.rda")
    load("~/.HPC_tidymass/MS_db/RPLC.database.rda")
    ##> MS1
        object_pos_anno <-
          annotate_metabolites_mass_dataset(
            object = object_pos_MRM,
            polarity = 'positive',
            ms1.match.ppm = 15,
            column = T.column,
            threads = 5,# for ms1 database, DO NOT recommand use all threads, maybe half of cores act better.
            database = kegg_ms1_database0.0.3,
            candidate.num = 2
          )
    ##> MS2
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
    ```
    
    3.2 The column name for 'sample information' must be consistent with that in the [sample_info.csv](https://github.com/ShawnWx2019/MetMiner/blob/main/02.DemoData/sample_info.csv). You can add more columns, such as 'treat', 'tissue', 'day', and so on.
    
    3.3 Steps of MRM selection. See `MDAtoolkits::extract_fragment`, **a)**. order fragment by intensity of each feature. **b)**. Add tags to fragments which have mz gap with precursor ions larger than 15. **c)**. pick the first tagged fragment as product ion. 
    
  

```bash
  # run Rscript ../01.Src/PreMetaboUpAnalysis.R -h check the explanation of each parameters
  Rscript ../01.Src/PreMetaboUpAnalysis.R \
          -x 02.DemoData/MS1 \ ## MS1 .mzXML file
          -y 02.DemoData/MS2 \ ## MS2 .mgf file
          -z 02.DemoData/sample_info.csv \ ## sample information
          -c 'rp' \ 
          -m 5 \
          -t 6 \
          -o T
          -g "QC"
```

4. check result

`PreMetaboUpAnalysis.R`carries out data cleansing and metabolite annotation on QC samples, selecting features that possess MS2 spectra. Then, following the MRM selection steps mentioned above, it generates a method for TQMS. You can check the results in the file named [pseudotargeted_QC_anno.xlsx](https://github.com/ShawnWx2019/MetMiner/blob/main/workdir/05.Annotation/pseudotargeted_QC_anno.xlsx). Additionally, please refer to the [workdir](https://github.com/ShawnWx2019/MetMiner/blob/main/workdir/) for intermediate files pertaining to data cleaning and metabolite annotation. If you find the default steps unsatisfactory, or if you need to conduct additional analyses, you can load all intermediate variables into R or Rstudio for subsequent analyses using the following code:

```r
# 'EnvName' is generated based on the current time for convenient traceback
load("workdir/EnvName.rda")
```
This will allow you to further manipulate and analyze the data as necessary. 

**Generate TQMS method from clean data**

For LC-MS data processed by Compound Discoverer (CD) or other software, you first need to export the intensity data for both the precursor and fragment ions. The order of the fragments is arranged based on their intensity. This data should be organized into the table format: [CDLC-MS.xlsx](https://github.com/ShawnWx2019/MetMiner/blob/main/02.DemoData/CDLC-MS.xlsx). The column name in bold-red is necessary. Then, you can generate a Triple Quadrupole Mass Spectrometry (TQMS) method using the `MDAtoolkits::mrm_selection_cd` function.


### 1.1.2 Data Cleaning of Metabolite TQMS Quantitative Results

After acquiring the TQMS method, we proceed to perform quantitative detection on all samples. Once the quantitative results are obtained, we can construct a mass_dataset following the steps outlined on the tidymass official website. This allows us to carry out subsequent data cleaning and standardization. Alternatively, you can use the Rscript:  [PresMetaboDataCleaning.R](https://github.com/ShawnWx2019/MetMiner/blob/main/01.Src/PresMetaboDataCleaning.R) to conduct a one-step data processing, directly yielding a standardized expression matrix, as well as the feature RSD and PCA results both before and after standardization. Results as [workdir2](https://github.com/ShawnWx2019/MetMiner/blob/main/workdir2/) 

**NOTICE：** When using the script for data cleaning, one important aspect to pay attention to is the removal of outliers. In metabolomics experiments involving heterogeneous samples (encompassing different tissues or species), substantial differences between samples can lead to a high proportion of missing values. Some of these might be mistakenly identified and removed as outliers. Therefore, when processing data that includes different tissues or species, add the `-g | --heterogeneous 'yes'` parameter. In addition, if your QC samples have been concentrated, they could also be erroneously removed as outliers. In such a case, it would be necessary to include the `-g | --heterogeneous ' yes` parameter as well. Additionally, if discrepancies exist among your samples, modifications to the script will be required. For example, in our study, we have randomly incorporated seven col-0 into the mutants, setting our filtering parameters at QC < 0.2, WT < 0.2, and MT < 0.5. The results within workdir2 were generated using default parameters, specifically applying only QC < 0.2. Consequently, the quantity of metabolites identified significantly exceeded the final count of metabolites selected in the publication.

```bash
 Rscript 01.Src/PresMetaboDataCleaining.R \
             -e 02.DemoData/exp_mat_v2.csv \
             -s 02.DemoData/sample_info_v2.csv \
             -n 'svr' \
             -g 'yes' &
```

## 1.2 untargeted Metabolomics

When analyzing large-scale untargeted metabolomics data, peak picking and library searching for hundreds or even thousands of samples simultaneously pose a challenge for personal computers. At this point, many commercial software solutions or PC-dependent LC-MS data analysis tools may struggle. 

After exploring various options, we chose to use tidyMass to complete the full process from data format conversion to peak picking, data cleaning, and metabolite annotation for large-scale untargeted metabolomics upstream data analysis. To facilitate usage on server-side, we have also set up an analysis workflow known as HPC-tidyMass.

**NOTIC:** To use this HPC-tidymass, several conditions must be met:

1. The user must be added to the Docker group, as the initial format conversion relies on MSconvert within Docker.
2. The naming convention for original files: QC as `QC_001`, sample as `S_0001`, with zeros used to pad the digit length.
3. When using the script, the absolute path must be used for the original file path, as Docker only recognizes absolute paths when running.

```bash
##> run
runTidymass \
          -i /home/data/shawn/project/fbox/01.rawdata \
          -t 1 \ ## DO NOT use 2, it dose not work now！
          -c 'rp' \
```
when runTidymass finish, all files, including original data and format-converted files such as .mzXML and .mgf, can be found in the `working_dir` directory. This also includes all mass_datasets produced from the peak picking to data normalization steps, as well as the mass_datasets which contains feature annotation. Moreover, in order to best retain the annotation results, we've categorized the annotation outcomes into four types. The first one is the unfiltered original annotation file, named `Original_annotation`. The second is the result after redundancy removal from the metabolites, named `rm_redundant`. The third type, `Only_MS2`, retains only the features that match with the database via MS2. The third category often has higher confidence; however, due to the limitations of the metabolites annotated in the database, many metabolites could be filtered out, particularly in plant materials abundant in secondary metabolites. 

**NOTIC:** The script has a mechanism in place for checkpoint resumption, which is useful when encountering errors. When an error arises, identify its cause, then check the result files in the `working_dir` to see if they are complete. If they are not, remove the corresponding result directory and rerun the script.


# Step 2. Downstream Analysis

Downstream data analysis of metabolomics is a critical step in metabolomics data mining. At present, most tools developed for metabolomics downstream analysis mainly serve animals, humans, and diseases. However, the data mining for plant metabolomics is relatively underdeveloped. To address this, we've developed `MDAtoolkits (beta)` to perform some conventional downstream metabolomics analyses. These include typical analyses such as differential metabolite (DAM) analysis, PCA analysis, metabolite classification, metabolite enrichment analysis, and more.

In addition, for large-scale metabolomics data, the swift, efficient, and accurate extraction of useful information from complex datasets has always been a challenge. To tackle this, we use the WGCNA method to quickly associate metabolites with samples. Through iterative removal of compounds that cannot be classified, we can rapidly and efficiently construct a weighted co-accumulation network of metabolites. At the same time, by integrating sample attributes, we can analyze the underlying biological issues at the metabolic level in the population scale.

## 2.1 Conventional Analysis and Visualization of Metabolomics Data

**Compound classification**   

ClassyFire is a comprehensive, flexible, and computable, purely structure-based chemical taxonomy (ChemOnt), developed by chemists, along with a computer program (ClassyFire) that uses only chemical structures and structural features to automatically assign all known chemical compounds to a taxonomy consisting of more than 4800 different categories [Djoumbou et.al](#refer-anchor-2). 

The ClassyFire database facilitates the retrieval of compound classifications by utilizing the compound's unique InChIKey. We generally execute batch conversions of metabolites using the [ClassyFire Batch by Fiehn Lab (cbf)](https://cfb.fiehnlab.ucdavis.edu/#/). Although some local databases offer corresponding InChIKeys for their compounds, a majority do not. For these unaccommodated metabolites, we must first procure the InChIKey from PubChem using the compound's name, then subsequently obtain the classification data via the CBF.

To streamline this process, we've developed a series of web crawler scripts to fetch the information. The functions are as follows:

1. `MDAtoolkits::mda_get_cid & MDAtoolkits::mda_pubchem_crawler`, convert compound name to pubchem cid and InChIKey via [webchem package](https://github.com/ropensci/webchem) [Szöcs E et.al](#refer-anchor-3).

2. `MDAtoolkits::mda_get_cid_fast`, convert compound name to pubchem InChIKey via [PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest).  **[Recommend]**

3. `MDAtoolkits::mda_classfire`, convert InChIKey to classyfire data via [classyfireR](https://github.com/aberHRML/classyfireR).

4. `MDAtoolkits::cbf_crawler`, convert InChIKey to classyfire data via [cbf api](https://cfb.fiehnlab.ucdavis.edu/). **[Recommend]**

**NOTICE** We propose that metabolite classification analysis be performed directly on the database, and that the classification information be incorporated into the database during its construction. We have already carried out classification analyses on the Arabidopsis Thaliana Knapsack and KEGG databases, embedding the classification data within them. Going forward, we will complete the classification annotations for other database.

**Different Accumulation Metabolites Analysis (DAM analysis)**

DAM constitutes a critical step in metabolomics data analysis, enabling direct comparison of metabolic differences between two type of samples. Common methods for identifying DAMs include using a t-test, analysis of variance (ANOVA), or Wilcoxon rank-sum test for p-value combined with log2 fold change for judgment. Alternatively, one could employ methods such as partial least squares-discriminant analysis (PLS-DA) or orthogonal partial least squares-discriminant analysis (OPLS-DA) to derive VIP (Variable Importance in Projection) values and log2 fold change as criteria. The `MDAtoolkits::DM_analysis` function allows for swift completion of differential testing as well as PLS-DA (or OPLS-DA). It also provides essential statistical metrics such as the p-value for differential significance testing, q-value for multiple testing corrections, log2 fold change representing the difference multiple, and the VIP (Variable Importance in Projection) value.

```r
mat <- 
object %>% 
  extract_expression_data()
res_case_vs_control = 
  DM_analysis(x = mat,
              left_index = c("S_001","S_002","S_003"),
              right_index = c("S_004","S_005","S_006"),
              left = "case",
              right = "control",
              method = 't-test',
              method2 = 'opls-da') 
```

In addition, to better integrate with tidymass, we've designed a set of visualization functions that directly operate on 'mass_datamass'. These functions facilitate the visualization of Principal Component Analysis (PCA)`mda_pca`, volcano plots for differential analysis`DAM_volcano `, as well as visual representation of the relative standard deviation (RSD) results for metabolite stability in quality control (QC) samples before and after normalization `mda_rsd`.

**KEGG or ClassyFire enrichment analysis**

Enrichment analysis is a powerful tool for exploring how metabolite sets influence metabolic pathways. For instance, conducting KEGG and ClassyFire enrichment analyses on differential metabolites or metabolites with similar expression patterns allows us to discern which metabolic pathways these biologically significant metabolite sets participate in, as well as their classification characteristics. Theoretically, we should use all metabolites contained in the species as a background, then confirm the significance of enrichment through hypergeometric distribution testing or Fisher's exact test. However, in practice, it's difficult to know all the metabolites a species contains, and the common approach is to use detected metabolites as the background for enrichment analysis. With this in mind, we can generate the `Compound2Term` and `Term2Name` files using `MDAtoolkits::mda_make_Keggdb & MDAtoolkits::mda_make_enrichdb`, and then perform enrichment analysis using `clusterProfiler::enricher`. 

*a) run KEGG enrichment analysis*

Currently, we have constructed the KEGG database for major plants and crops (rice, corn, wheat, soybeans, rapeseed, cotton, Arabidopsis, TBtools-keggbackend), which can be found in the [03.Document](https://github.com/ShawnWx2019/MetMiner/blob/main/03.Document) folder.

The detailed steps are as follows: 

First, you need to search for the corresponding Organisms code for the species through the website https://www.genome.jp/kegg/catalog/org_list.html. Then, you can construct the species metabolite kegg database through the following code.
 
```r
library(MDAtoolkits)
library(progressr)
handlers(handler_pbcol(
      adjust = 1.0,
    complete = function(s) cli::bg_red(cli::col_black(s)),
  incomplete = function(s) cli::bg_cyan(cli::col_black(s))
))
##> For Brassica napus (rape)
bna_kegg_db <- mda_make_keggdb(organism = 'bna')
##> save the long-time used bna kegg database.
save(bna_kegg_db,file = "03.Document/bna_kegg_db.rda")
```

Second, 、 acquire the corresponding KEGG CID for the compounds through databases like KEGG, [CTS](http://cts.fiehnlab.ucdavis.edu/batch), etc. We can accomplish the automated conversion by using `MDAtoolkits::mda_name2kegg & MDAtoolkits::CTS_kegg`. 

```r
##> import annotation result
anno <- readxl::read_xlsx("02.DemoData/compound_anno.xlsx",sheet = 1)
##> filter out the compounds without KEGG.id (some database contains keggid, others are not.)
query <- anno %>% 
  mutate(KEGG.ID = str_split(KEGG.ID.iso,"|",2,T)[,1]) %>% ## get the represent KEGG.id
  filter(is.na(KEGG.ID) & !str_detect(Compound_name,"^MW:")) %>% ## remove features already have KEGG.ID and features without compound name.
  pull(Compound_name) %>% URLencode() %>% unique()## encode the Compound name as URL for web wrawler.
library(MDAtoolkits)
name2kegg1 = mda_name2kegg(query = query[1:10]) %>% filter(KEGG != "")
name2kegg2 = mda_CTS_kegg(query = query[1:10],key = 'name') %>% filter(KEGG != "")
name2kegg3 = rbind(name2kegg1,name2kegg2) 
name2kegg4 <- anno %>% 
  mutate(KEGG.ID = str_split(KEGG.ID.iso,"|",2,T)[,1]) %>% 
  filter(is.na(KEGG.ID)) %>% 
  select(Compound_ID,Compound_name) %>% 
  inner_join(name2kegg3,by = "Compound_name") 
ID2kegg <- anno %>% 
  mutate(KEGG.ID = str_split(KEGG.ID.iso,"|",2,T)[,1]) %>% 
  filter(!is.na(KEGG.ID)) %>% select(Compound_id,Compound_name,KEGG.ID) %>% 
  rename("KEGG" = "KEGG.ID") %>% 
  rbind(name2kegg4) %>% 
  select(Compound.id,KEGG) %>% 
  distinct()
##> load ath databse
load("03.Document/ath_kegg_db.rda")
##> customized Term 2 Compound_id
TERM2COMPOUND = ath_kegg_db$TERM2GENE %>% inner_join(ID2TERM,by = c("CID"="KEGG"),multiple = "all") %>% 
  distinct() %>% select(TERM,Compound_id)
##> KEGG enrichment analysis
test_c = anno$Compound_id[1:400] ## head 400 compounds.
library(clusterProfiler)
res.kegg <- enricher(
  gene = test_c,
  pvalueCutoff = 1, ##> Randomly selected metabolites are likely not to show significant enrichment, hence a p-value of 1 was chosen. When doing this on your own, you can set a stricter threshold for the p-value, such as 0.05.
  qvalueCutoff = 1,
  TERM2GENE = TERM2COMPOUND,TERM2NAME = ath_kegg_db$TERM2NAME
)
##> visualize
dotplot(res.kegg,showCategory = 15) ## dotplot
enrichplot::cnetplot(res.kegg) ## term - compound network
res.sim <- enrichplot::pairwise_termsim(res.kegg)
enrichplot::emapplot(res.sim) ## KEGG term network
```

*b) run CLassyFire enrichment analysis*

ClassyFire has not created a species specific database, so we have to build a database derived from the classification results, and then carrying out enrichment analysis. 

Code:

```r
##> classyfire result via MDAtoolkits::cbf_crawler
class <- readxl::read_xlsx("02.DemoData/compound_anno.xlsx",sheet = 2)
##> transform
class2 <- 
  class %>% 
  filter(!is.na(superclass)) %>% 
  select(-description,-InChIKey) %>% 
  separate_longer_delim(cols = `parent levels`,delim = " | ") %>% 
  select(-`Lab ID`) %>% 
  pivot_longer(!all_of(c("Compound ID","Compound name")),names_to = "levels",values_to = "NAME") %>% 
  filter(!is.na(NAME))
  
##> construct database
TERM2NAME = 
  class2 %>% 
  select(NAME) %>% 
  distinct() %>% 
  mutate(TERM = paste0("MC:",str_pad(c(1:nrow(.)),5,"left",'0'))) %>% 
  filter(NAME != "NA") %>% 
  select(TERM,NAME)

TERM2COM = 
  class2 %>% 
  left_join(TERM2NAME) %>% 
  select(TERM,`Compound ID`) %>% 
  setNames(c("TERM","COMPOUND")) %>% 
  distinct()

writexl::write_xlsx(list(
  db = class2,
  TERM2NAME = TERM2NAME,
  TERM2COM = TERM2COM
),path = "02.DemoData/classyfire_database.xlsx")

##> enrichment

library(clusterProfiler)

res.class <- enricher(
  gene = test_c,pvalueCutoff = 1,
  qvalueCutoff = 1,
  TERM2GENE = TERM2COM,
  TERM2NAME = TERM2NAME
)

dotplot(res.class,showCategory = 15) 
enrichplot::cnetplot(res.class)

res.sim <- enrichplot::pairwise_termsim(res.class)

enrichplot::emapplot(res.sim)
```

**Complete DAM analysis pipeline**

We have developed a streamlined, one-step workflow tailored for paired DAM analyses. This comprehensive process encompasses differential analysis, volcano plot generation, PCA, heatmap construction, and KEGG enrichment analysis. Once the enrichment analysis is completed, we proceed by extracting the metabolites from significantly enriched pathways amidst the differential metabolites, followed by their visualization in a heatmap. This notably enhances the ease and efficiency of data mining within our enrichment analysis results. [run_DAM_analysis.R](https://github.com/ShawnWx2019/MetMiner/blob/main/01.Src/run_DAM_analysis.R)

```bash
	Rscript run_DAM_analysis_v3.R  \
		--peakfile expmat.txt \
		--group group.txt \
		--meta_anno annotation.txt \ # must have compound_id compound_name and KEGG.ID 
		--kegg_db ath_kegg.rds \
		--left "left" \
		--right 'right' \
		--VIP 1 \
		--pvalue 0.05 \
		--qvalue 1 \
		--log2fc 0.26 \
    --test.method "t-test" \
    --pls.method "pls-da"
```

group file:

group|sample_id
---|---
left|s_01
left|s_02
left|s_03
right|s_04
right|s_05
right|s_06

## 2.2 Advanced Analysis of Metabolomics Data

In the face of metabolomics data mining involving large sample sizes, the extraction of useful information from complex datasets remains a considerable challenge. After testing clustering techniques like k-means and hclust, we discovered that the iterative Weighted Gene Co-expression Network Analysis (iWGCNA) represents an efficient method for data mining. It swiftly and accurately aggregates metabolites exhibiting shared accumulation trends. By delving into the types of compounds encapsulated within each module, alongside their affected metabolic pathways, we are able to uncover a wealth of valuable information.

While we can adhere to the WGCNA manual for a step-by-step analysis, we have developed a user-friendly Shiny app for conducting WGCNA analyses to enhance convenience. This software has become highly stable after several updates. Moreover, we've crafted a one-step WGCNA script, incorporating interactive elements through readline and plotly. We are now capable of excluding outliers and selecting suitable power values, dependent on the progress of our analysis. Instructions for operating the Shiny app can be found here: [WGCNA ShinyApp](https://github.com/ShawnWx2019/WGCNA-shinyApp/blob/main/README.md). The one-step WGCNA can be executed with the code provided below, and the final results can be previewed within [workdir3](https://github.com/ShawnWx2019/MetMiner/blob/main/workdir3/).

```r
library(ShinyWGCNA)
library(tidyverse)
library(tidymass)

expmat <- read.delim("02.DemoData/exp_mat_final.txt")
category <- read.delim("02.DemoData/category.txt")
dir.create(path = "workdir3")
ShinyWGCNA::oneStepWGCNA(
  exp.mat = expmat,
  iterative = T,
  pheofile = category,
  lazy_modle = F,
  save_fig = T,
  F.path = "workdir3/",
  datatype = "peak area",r_cutoff = 0.15
)

```

## Reference

<div id=#refer-anchor-1>

1. Shen, X., Yan, H., Wang, C., Gao, P., Johnson, C.H. & Snyder, M.P. (2022) TidyMass an object-oriented reproducible analysis framework for LC--MS data. Nature Communications, 13, 4365. DOI: [https://doi.org/10.1038/s41467-022-32155-w](https://www.nature.com/articles/s41467-022-32155-w)

</div>

<div id=#refer-anchor-2>

2. Djoumbou Feunang, Y., Eisner, R., Knox, C., et al. (2016) ClassyFire: automated chemical classification with a comprehensive, computable taxonomy. Journal of Cheminformatics, 8, 61. DOI: [https://doi.org/10.1186/s13321-016-0174-y](https://doi.org/10.1186/s13321-016-0174-y)

</div>

<div id=#refer-anchor-3>

3. Szöcs E, Stirling T, Scott ER, et al. (2020) webchem: An R Package to retrieve chemical information from the web. Journal of statistical software 93:. DOI: [https://doi.org/10.18637/jss.v093.i13](https://doi.org/10.18637/jss.v093.i13)

</div>
