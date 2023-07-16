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

# Step 1. Upstream analysis

Here we refer to the process of obtaining the metabolomics data from LC-MS raw data to the generation of a metabolite accumulation matrix as <font color=blue>**「upstream data analysis data」**</font>. It includes processes such as <font color=green> 「raw data format conversion」, 「peak picking」, 「data cleaning」, and 「data normalization」</font>. In this pipeline, this part is typically performed by TidyMass or Compound Discoverer. If you have used other software for peak picking, you can refer to https://www.tidymass.org/start/create_mass_dataset/ to create a mass_dataset, which can then be used for subsequent analysis with TidyMass.

## 1.1 Pseudotargeted metabolomics

For pseudotargeted metabolomics, in the initial stage, HPLC experiments can be conducted using 5-10 QC (quality control) samples. Since the data volume is relatively small, and considering the completeness of commercial software databases, we recommend utilizing dedicated analysis software for peak detection and library searching. Alternatively, you can also opt to perform upstream analysis directly with TidyMass, followed by MRM (Multiple Reaction Monitoring) selection to generate a TQMS (Triple Quadrupole Mass Spectrometry) method.

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

For LC-MS data processed by Compound Discoverer (CD) or other software, you first need to export the intensity data for both the precursor and fragment ions. The order of the fragments is arranged based on their intensity. This data should be organized into the table format: [CDLC-MS.xlsx](https://github.com/ShawnWx2019/MetMiner/blob/main/CDLC-MS.xlsx). The column name in bold-red is necessary. Then, you can generate a Triple Quadrupole Mass Spectrometry (TQMS) method using the `MDAtoolkits::mrm_selection_cd` function.




## Reference

<div id=#refer-anchor-1></div>

1. Shen, X., Yan, H., Wang, C., Gao, P., Johnson, C.H. & Snyder, M.P. (2022) TidyMass an object-oriented reproducible analysis framework for LC--MS data. Nature Communications, 13, 4365. DOI: [https://doi.org/10.1038/s41467-022-32155-w](https://www.nature.com/articles/s41467-022-32155-w)


