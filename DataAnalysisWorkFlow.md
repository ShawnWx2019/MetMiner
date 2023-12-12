# Data analysis code of MetMiner paper

## Section 1. Data cleaning

**Code: ** [01.TidyMassDataCleaning.R](https://github.com/ShawnWx2019/MetMiner/blob/main/01.Src/01.TidyMassDataCleaning.R)

- 1.1 Construct massdataset by tidymass, introduce TQMS data matrix and sample information to tidymass.  

- 1.2 Remove noise features: Missing values (mv) in QC `0.2`, Col-0 plant `0.2` and F-box mutants `0.5`  

- 1.3 Missing value imputation: `Knn` method  

- 1.4 Normalization: `SERRF` method

## Section 2. Re-annotation with tidymass

**Code: **[02.ReannotationWithTidymass.R](https://github.com/ShawnWx2019/MetMiner/blob/main/01.Src/02.ReannotationWithTidymass.R)  

- 2.1 Remove redundance between Pos model and Neg model: mw `1 dalton` and rt in `0.3 min`  

- 2.2 Import MS2 information to massdataset from HPLC QC samples. 

    a. MSconvert [.raw ==> .mgf]  
    
    b. Tidymass [`mutate_ms2`]  
    
    c. mz.tol [10 ppm] & rt.tol [35 second]  
    
- 2.3 MS1 database
    a. knapsack_ath_db.rda;  
    
    b. ath_plantcyc.database.rda;  
    
    c. kegg_ms1_database0.0.3.rda; 

- 2.4 MS2 database
    a. ReSpect_database2.0.rda;
    
    b. PlaSMA_database2.0.rda;
    
    c. mona_database0.0.4.rda;
    
    d. massbank_database0.0.4.rda
    
## Section 3. PmHub annotation and merge results with tidymass annotations

**Code: **  ~~[03.PmHubMS1Annotation.R](https://github.com/ShawnWx2019/MetMiner/blob/main/01.Src/03.PmHubMS1Annotation.R)~~  
[04.MergePMHubOnlineAnnotation.R](https://github.com/ShawnWx2019/MetMiner/blob/main/01.Src/04.MergePMHubOnlineAnnotation.R)

> [!WARNING]
> The annotation of PMhub MS1 was abandoned. Due to too much redundancy, we subsequently conducted an online annotation of P Mhub based on the original MS2 information, which achieved better results.  

## Section 4. Metabolites classification and KEGG CID annotation

**Code: **  [05.MetabolitesClassification.R](https://github.com/ShawnWx2019/MetMiner/blob/main/01.Src/05.MetabolitesClassification.R)
[06.CompoundKEGGcidAnnotation.R](https://github.com/ShawnWx2019/MetMiner/blob/main/01.Src/06.CompoundKEGGcidAnnotation.R)  

- Get pubchem cid and InChIKey by compound names.  
    a. MDAtoolkits [`mda_get_cid_fast()`];  
    
    b. The rest compounds InChIKey were obtained from searching other online database, such as Knapsack, metacyc etc.. manually  
    
- Get classification information via InChIKey   MDAtoolkits [`cbf_crawler()`]  

- Get KEGG cid via compound name, InChIKey or SMILE MDAtoolkits 
    a. KEGG database [`MDAtoolkits::mda_name2kegg`]  
    
    b. KEGG database [`MDAtoolkits::mda_CTS_kegg`]
    
## Section 5. Enrichment analysis 
 
**Code: **  [07.CompoundEnrichmentAnalysis.R](https://github.com/ShawnWx2019/MetMiner/blob/main/01.Src/07.CompoundEnrichmentAnalysis.R)

- Construct classification enrichment database,
     
     a. **TERM2NAME**: [TERM] Class ID (add by myself) to [Name] Decription: eg: C_0001 \t flavanol；
     
     b. **TERM2GENE**: [TERM] Class ID to [Gene] Metabolite ID: eg: C_0001 \t Com_00001;
     
- Enrichment analysis clusterProfiler [`enricher`]

- Construct KEGG database

     a. Extract plant metabolic pathways: use [TBtools KEGG plant backend].
     
     b. According to the plant metabolic pathway, crawl the KEGG CID of the metabolites in the pathway from the kegg database. 
     
     c. Filter plant specific KEGG cid annotation
     
     d. KEGG pathway summary and enrichment analysis
