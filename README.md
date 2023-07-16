# Users Manual of MetMiner

[![R version](https://img.shields.io/badge/R-v4.1.1-salmon)](https://www.r-project.org) ![lifecycle](https://img.shields.io/badge/lifecycle-Experimental-lightcyan) [![license](https://img.shields.io/badge/license-MIT-red)](https://opensource.org/licenses/MIT) [![Myblog](https://img.shields.io/badge/Blog-ShanwLearnBioinfo-purple)](https://shawnwx2019.github.io/)

An Integrated Pipeline for Large-Scale Metabolomics Data Processing and Data Mining.

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

For LC-MS data processed by Compound Discoverer (CD) or other software, you first need to export the intensity data for both the precursor and fragment ions. The order of the fragments is arranged based on their intensity. This data should be organized into the table format: [CDLC-MS.xlsx](https://github.com/ShawnWx2019/MetMiner/blob/main/02.DemoDat6a/CDLC-MS.xlsx). The column name in bold-red is necessary. Then, you can generate a Triple Quadrupole Mass Spectrometry (TQMS) method using the `MDAtoolkits::mrm_selection_cd` function.


### 1.1.2 Data Cleaning of Metabolite TQMS Quantitative Results

After acquiring the TQMS method, we proceed to perform quantitative detection on all samples. Once the quantitative results are obtained, we can construct a mass_dataset following the steps outlined on the tidymass official website. This allows us to carry out subsequent data cleaning and standardization. Alternatively, you can use the Rscript:  [PresMetaboDataCleaning.R](https://github.com/ShawnWx2019/MetMiner/blob/main/01.Src/PresMetaboDataCleaning.R) to conduct a one-step data processing, directly yielding a standardized expression matrix, as well as the feature RSD and PCA results both before and after standardization. Results as [workdir2](https://github.com/ShawnWx2019/MetMiner/blob/main/workdir2/) 

**NOTICE：** When using the script for data cleaning, one important aspect to pay attention to is the removal of outliers. In metabolomics experiments involving heterogeneous samples (encompassing different tissues or species), substantial differences between samples can lead to a high proportion of missing values. Some of these might be mistakenly identified and removed as outliers. Therefore, when processing data that includes different tissues or species, add the `-g | --heterogeneous 'yes'` parameter. In addition, if your QC samples have been concentrated, they could also be erroneously removed as outliers. In such a case, it would be necessary to include the `-g | --heterogeneous ' yes` parameter as well.

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

## 2.1 Conventional  Analysis and Visualization of Metabolomics Data

**Compound classification**   

ClassyFire is a comprehensive, flexible, and computable, purely structure-based chemical taxonomy (ChemOnt), developed by chemists, along with a computer program (ClassyFire) that uses only chemical structures and structural features to automatically assign all known chemical compounds to a taxonomy consisting of more than 4800 different categories [Djoumbou et.al](#refer-anchor-2). 

The ClassyFire database facilitates the retrieval of compound classifications by utilizing the compound's unique InChIKey. We generally execute batch conversions of metabolites using the [ClassyFire Batch by Fiehn Lab (cbf)](https://cfb.fiehnlab.ucdavis.edu/#/). Although some local databases offer corresponding InChIKeys for their compounds, a majority do not. For these unaccommodated metabolites, we must first procure the InChIKey from PubChem using the compound's name, then subsequently obtain the classification data via the CBF.

To streamline this process, we've developed a series of web crawler scripts to fetch the information. The functions are as follows:

1. `MDAtoolkits::mda_get_cid & MDAtoolkits::mda_pubchem_crawler`, convert compound name to pubchem cid and InChIKey via [webchem package](https://github.com/ropensci/webchem)[Szöcs E et.al](#refer-anchor-3)

2. `MDAtoolkits::mda_get_cid_fast`, convert compound name to pubchem InChIKey via [PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)  **[Recommend]**

3. `MDAtoolkits::mda_classfire`, convert InChIKey to classyfire data via [classyfireR](https://github.com/aberHRML/classyfireR)

4. `MDAtoolkits::cbf_crawler`, convert InChIKey to classyfire data via [cbf api](https://cfb.fiehnlab.ucdavis.edu/) **[Recommend]**

**[NOTICE]** We propose that metabolite classification analysis be performed directly on the database, and that the classification information be incorporated into the database during its construction. We have already carried out classification analyses on the Arabidopsis Thaliana Knapsack and KEGG databases, embedding the classification data within them. Going forward, we will complete the classification annotations for other database.

**Different Accumulation Metabolites Analysis (DAM analysis)**








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
