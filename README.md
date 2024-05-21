# User' Manual
## Overview
In this work, we integrate large GWAS data for breast, colorectal, lung, ovarian, pancreatic, and prostate cancers and population-sale proteomics data from 96,991 participants combing from Atherosclerosis Risk in Communities study (ARIC), deCODE genetics, and UK Biobank Pharma Proteomics Project (UKB-PPP) to identify risk proteins associated with each cancer. We prioritize therapeutic drugs for potential repurposing based on druggable risk proteins targeted by approved drugs or undergoing clinical trials for cancer treatment or other indications. We further evaluate the effect of cancer risk for those drugs approved for the indications, using over 3.5 million EHR database at Vanderbilt University Medical Center (VUMC). Findings from this study offers novel insights into therapeutic drugs targeting risk proteins for cancer prevention and intervention.

![My Image](Figure1.png)

**Step1:** Identify the most significant genetic variants associated with six cancer types (breast, ovary, prostate, colorectum, lung, and pancreas) using the latest GWAS data. Perform statistical fine-mapping to pinpoint the lead variants at each risk locus for each cancer type. Integrate GWAS and Fine-mapping results to pinpoint independent cancer risk variants. Map these genetic independent cancer risk variants to pQTLs from ARIC, deCODE, and UKB-PPP, to pinpoint putative cancer risk proteins (*P* < 0.05, Bonferroni correction). Perform colocalization analyses to identify cancer risk proteins with high confidence by evaluating the likelihood of shared causal variants between pQTLs and GWAS.

**Step2:** Annotate the identified cancer risk proteins using drug information from DrugBank, ChEMBL, the Therapeutic Target Database, and OpenTargets. Identify druggable proteins that are either therapeutic targets of approved drugs or are undergoing clinical trials.

**Step3** Utilize electronic health records from VUMC to emulate treated-control drug trials for drugs approved for indications other than cancer. Apply the Inverse Probability of Treatment Weighting framework and use the Cox proportional hazard model to evaluate the cancer risk.

## Methods
### Meta
A fixed effect meta-analyses was applied between pQTLs from ARIC and deCODE to increase the power of pQTLs.
- Executive code: meta_pQTLs.R
- Example files: meta_example_file.tsv

### coloc
A genetic colocalisation analysis of protein expression and cancer risk was conducted to check whether they share common genetic causal variants. Protein with posterior probability (PP.H4) > 0.5 are considered as cancer risk proteins.

- Executive code: coloc_pQTLs_eQTLs_GWAS.R
- Example files:
1) GWAS summary statistics: coloc_breast_dbSNPs_impute_summary_statistic_polyfun.txt
2) pQTLs: coloc_ANXA4.all.pqtls.csv

### SMR+HEIDI
SMR(Summary-data-based Mendelian Randomization) is a method that uses summary data from GWAS and protein QTL studies to test for a causal relationship between protein expression and cancer risk. A protein Bonferroni-adjusted SMR *P* < 0.05 is considered as significant. A followed HEIDI test is performed on significant SMR results to determine if the colocalized signals can be explained by one single causal variant or by multiple causal variants in the locus. HEIDI *P* >= 0.05 (no obvious evidence of heterogeneity of estimated effects or linkage).

- Executive code:
```
/method/smr-1.3.1-linux-x86_64/smr-1.3.1
--bfile /hg19/1kg.chr1.phase3.20130502.Eur
--gwas-summary /BRCA/CA14.rs12048493.smr.ma
--beqtl-summary /BRCA/CA14.rs12048493.smr
--out BRCA.CA14.rs12048493
--target-snp rs12048493 --thread-num 20
```
- Example files:
1) GWAS summary statistics: CA14.rs12048493.smr.ma
2) pQTLs: CA14.rs12048493.smr
3) 1000G reference files (1kg.chr1.phase3.20130502.Eur) are downloaded from [1000G project](https://www.internationalgenome.org/category/genotypes/)

### Inverse probability of treatment weighting (IPTW)


### Cox proportional hazard model

### Figures
`./Figures/Figure4.tgz`
`./Figures/Figure5.tgz`
`./Figures/FigureS3.tgz`

## Contact
Qing Li: qing.li@vumc.org
Qingyuan Song:
Zhijun Yin:
Xingyi Guo: xingyi.guo@vumc.org
