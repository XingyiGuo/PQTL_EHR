# User' Manual
## Overview
In this work, we integrate large GWAS data for breast, colorectal, lung, ovarian, pancreatic, and prostate cancers and population-sale proteomics data from 96,991 participants combing from Atherosclerosis Risk in Communities study (ARIC), deCODE genetics, and UK Biobank Pharma Proteomics Project (UKB-PPP) to identify risk proteins associated with each cancer. We prioritize therapeutic drugs for potential repurposing based on druggable risk proteins targeted by approved drugs or undergoing clinical trials for cancer treatment or other indications. We further evaluate the effect of cancer risk for those drugs approved for the indications, using over 3.5 million EHR database at Vanderbilt University Medical Center (VUMC). Findings from this study offers novel insights into therapeutic drugs targeting risk proteins for cancer prevention and intervention.

![My Image](Figure1.png)

**Step1:** Identify the most significant genetic variants associated with six cancer types (breast, ovary, prostate, colorectum, lung, and pancreas) using the latest GWAS data. Perform statistical fine-mapping to pinpoint the lead variants at each risk locus for each cancer type. Integrate GWAS and Fine-mapping results to pinpoint independent cancer risk variants. Map these genetic independent cancer risk variants to pQTLs from ARIC, deCODE, and UKB-PPP, to pinpoint putative cancer risk proteins (*P* < 0.05, Bonferroni correction). Perform colocalization analyses to identify cancer risk proteins with high confidence by evaluating the likelihood of shared causal variants between pQTLs and GWAS.

**Step2:** Annotate the identified cancer risk proteins using drug information from DrugBank, ChEMBL, the Therapeutic Target Database, and OpenTargets. Identify druggable proteins that are either therapeutic targets of approved drugs or are undergoing clinical trials.

**Step3** Utilize electronic health records from VUMC to emulate treated-control drug trials for drugs approved for indications other than cancer. Apply the Inverse Probability of Treatment Weighting framework and use the Cox proportional hazard model to evaluate the cancer risk.

## Methods
### Meta

### coloc


### SMR+HEIDI

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
