#module R/4.2.0

library(coloc)
library(hash)
library(data.table)
library(stats)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

ea.prots <- read.csv("SixCancersIndexVariants_PQTLs_sig_surrogate_20230508.csv")
ea.sig=ea.prots[ea.prots$Below_Pcutoff=="True",]

d <- hash() 
d[["BRCA"]] <- "breast"
d[["Lung"]] <- "lung"
d[["OV"]] <- "ovarian"
d[["PAAD"]] <- "pancreatic"
d[["PRAD"]] <- "prostate"

gwas_sample <- hash()
gwas_sample[["BRCA"]] <- "228951"
gwas_sample[["Lung"]] <- "112610"
gwas_sample[["OV"]] <- "63347"
gwas_sample[["PAAD"]] <- "14998"
gwas_sample[["PRAD"]] <- "140306"

gwas_case_s <- hash()
gwas_case_s[["BRCA"]] <- "0.537"
gwas_case_s[["Lung"]] <- "0.291"
gwas_case_s[["OV"]] <- "0.354"
gwas_case_s[["PAAD"]] <- "0.552"
gwas_case_s[["PRAD"]] <- "0.564"

###specific for gwas
gwas_exsits=FALSE
cur_tsname=ea.sig[1,1]
if (file.exists(paste0(d[[cur_tsname]], '_dbSNPs_impute_summary_statistic_polyfun.txt'))){gwas_exsits=TRUE}
if(gwas_exsits){
	print(paste0("reading gwas for ", cur_tsname))
	gwas = fread(paste0(d[[cur_tsname]], '_dbSNPs_impute_summary_statistic_polyfun.txt'),sep="\t")
        gwas <- gwas %>% mutate(minor_freq = ifelse(A1FREQ<=0.5, A1FREQ, 1-A1FREQ))
	gwas = gwas[complete.cases(gwas), ]
	gwas = gwas[gwas$minor_freq<1, ]
	gwas = gwas[gwas$minor_freq>=0.01, ]
        gwas = gwas[gwas$PVALUE<0.5, ]

}else{
	print(paste0("check gwas for ", cur_tsname))
	stop
}


##tsname= args[1]#for non tissue specific 
system.time(for(i in 1:dim(ea.sig)[1]){
	print(i)
	pnames <- ea.sig[i,19]
	tsname <- ea.sig[i,1] #tissue specific proteins
	print(paste0(tsname, "_", pnames))
	
	eqtls_exsits=FALSE
	pqtls_exsits=FALSE
	
	pph4_pqtlseqtls <- NULL
	pph4_pqtlsgwas <- NULL
	
	##Part 1: pQTLs and eQTLs colocalization
	if (file.exists(paste0("/eqlts_pqtls_tss500k_refine/",tsname,"/",pnames,".all.pqtls.csv"))){pqtls_exsits=TRUE}
	
	if (file.exists(paste0("/eqtls_193proteins/",tsname, "/",pnames,".eqtls.csv"))){eqtls_exsits=TRUE}
	
	if(eqtls_exsits && pqtls_exsits){
		print("reading eqtls")
		eqtls = fread(paste0("/eqtls_193proteins/",tsname, "/",pnames,".eqtls.csv"),sep=",")
		eqtls = eqtls[complete.cases(eqtls), ]
		eqtls = eqtls[eqtls$maf<1, ]
		eqtls = eqtls[eqtls$maf>0.01, ]
		
	    print("reading pqtls")
		pqtls = fread(paste0("/eqlts_pqtls_tss500k_refine/",tsname,"/",pnames,".all.pqtls.csv"),sep=",")
		pqtls = pqtls[complete.cases(pqtls), ]
		pqtls = pqtls[pqtls$ImpMAF<1, ]
		pqtls = pqtls[pqtls$ImpMAF>=0.01, ]
		
		#order by variant_id
		eqtls= eqtls[order(eqtls$rsids),]
		pqtls= pqtls[order(pqtls$rsids),]
	
		common = intersect(as.character(pqtls$rsids),as.character(eqtls$rsids))
		print(length(common))
		if(length(common)>0){
			eqtls.1 = eqtls[match(as.character(common),as.character(eqtls$rsids)),]
			pqtls.1 = pqtls[match(as.character(common),as.character(pqtls$rsids)),]

			# print(eqtls.1[!complete.cases(eqtls.1),])
			# print(pqtls.1[!complete.cases(pqtls.1),])
			pqtls.1$N <- gsub("nocommon.*", pqtls.1$N[1], pqtls.1$N)
			
			ds1 = list(eqtls.1$rsids, as.numeric(eqtls.1$N), as.numeric(eqtls.1$pval_nominal), as.numeric(eqtls.1$maf))
			names(ds1) = c("snp","N","pvalues","MAF")
			ds2 = list(pqtls.1$rsids,as.numeric(pqtls.1$N), as.numeric(pqtls.1$meta_Pvalue), as.numeric(pqtls.1$ImpMAF))
			names(ds2) = c("snp","N","pvalues","MAF")

			ds1$snp = as.character(ds1$snp)
			ds2$snp = as.character(ds2$snp);
			ds1$type = ds2$type = "quant"
			c1 = coloc.abf(ds1, ds2)
			pph4_pqtlseqtls <- rbind(pph4_pqtlseqtls,c1$summary[-1])
			write.table(pph4_pqtlseqtls,paste0("/coloc_res_eqtls_pqtls_tss500k/",tsname, ".",pnames,".pqtlseqtls.coloc.txt"),col.names = T,row.names = F,quote = F)
		}else{
			print("No common rs found among two datasets")
		}
	}else if (pqtls_exsits){
		print("eqtls does not exist")
	}else if(eqtls_exsits){
		print("pqtls does not exist")
	}else{
		print("Neither pqtls nor eqtls exist")
	}

	##Part 2: pQTLs and GWAS colocalization
	if (file.exists(paste0("/gwas_pqtls_rs50k_refine/",tsname,"/",pnames,".all.pqtls.csv"))){pqtls_exsits=TRUE}
	if (file.exists(paste0(d[[cur_tsname]], '_dbSNPs_impute_summary_statistic_polyfun.txt'))){gwas_exsits=TRUE}
	
	if (gwas_exsits && pqtls_exsits){
		if(cur_tsname == tsname){
			print("gwas has been loaded successfully")
		}else if (tsname != "COREAD"){
			print("reading gwas")
			gwas = fread(paste0(d[[tsname]], '_dbSNPs_impute_summary_statistic_polyfun.txt'),sep="\t")
                        gwas <- gwas %>% mutate(minor_freq = ifelse(A1FREQ<=0.5, A1FREQ, 1-A1FREQ))
						gwas = gwas[complete.cases(gwas), ]
						gwas = gwas[gwas$minor_freq<1, ]
						gwas = gwas[gwas$minor_freq>=0.01, ]
						gwas = gwas[gwas$PVALUE<0.5, ]
			cur_tsname=tsname
		}else{
			print("COREAD gwas not exists")
			cur_tsname=tsname
			gwas_exsits=FALSE
			next
		}
		
		print("reading pqtls")
		pqtls = fread(paste0("/gwas_pqtls_rs50k_refine/",tsname,"/",pnames,".all.pqtls.csv"),sep=",")
		pqtls = pqtls[complete.cases(pqtls), ]
		pqtls = pqtls[pqtls$ImpMAF<1, ]
		pqtls = pqtls[pqtls$ImpMAF>=0.01, ]
		
		common.1 = intersect(as.character(pqtls$rsids),as.character(gwas$SNP))
		print(length(common.1))
		if(length(common.1)>0){
			gwas.1 = gwas[match(as.character(common.1),as.character(gwas$SNP)),]
			pqtls.2 = pqtls[match(as.character(common.1),as.character(pqtls$rsids)),]
			pqtls.2$N <- gsub("nocommon.*", pqtls.2$N[1], pqtls.2$N)
			
			ds3 = list(gwas.1$SNP, as.numeric(rep(gwas_sample[[tsname]],dim(gwas.1)[1])), as.numeric(gwas.1$PVALUE), as.numeric(gwas.1$A1FREQ))
			names(ds3) = c("snp","N","pvalues","MAF")
			ds4 = list(pqtls.2$rsids,as.numeric(pqtls.2$N), as.numeric(pqtls.2$meta_Pvalue), as.numeric(pqtls.2$ImpMAF))
			names(ds4) = c("snp","N","pvalues","MAF")

			ds3$snp = as.character(ds3$snp)
			ds4$snp = as.character(ds4$snp)
			ds3$type = "cc"
			ds4$type = "quant"
			ds3$s=as.numeric(gwas_case_s[[tsname]])
			c2 = coloc.abf(ds4, ds3)
			pph4_pqtlsgwas <- rbind(pph4_pqtlsgwas,c2$summary[-1])
			write.table(pph4_pqtlsgwas,paste0("/coloc_res_gwas_pqtls_rs50k/",tsname,".",pnames,".pqtlsgwas.coloc.txt"),col.names = T,row.names = F,quote = F)
		}else{
			print("No common rs found among two datasets")
		}
	}
})