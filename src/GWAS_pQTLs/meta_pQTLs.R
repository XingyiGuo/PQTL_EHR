library(meta)
args <- commandArgs(trailingOnly=T)
prefix <- args[1]
input <- paste0("./both/",prefix,".PMID35501419_pQTL_cis_NG.out")
df <- read.table(input,header=F,stringsAsFactors=F,sep="\t",fill=T)
n <- dim(df)[1]


sink(paste0(input,"_meta_04072023"))
for(i in 1:n)
{
	if(!grepl("nocommon",df[i,11]) && ! grepl("nocommon",df[i,29])){
	
		# the alleles for beta value in both datasets are the same 
		if(df[i,7] == df[i,27])
		{
		beta <- as.numeric(c(df[i,11],df[i,29]))
		se <- as.numeric(c(df[i,12],df[i,32]))
		meta <- metagen(beta,se,common=gs("common"))
		cat(cat(c(as.character(df[i,]), meta$TE.common, meta$seTE.common, meta$pval.common,"meta",df[i,27]),sep="\t"),"\n")
		}
	
		# the alleles for beta value in both datasets are differents
		if(df[i,7] != df[i,27])
		{
		x <- as.numeric(df[i,11])*(-1)
		beta <- as.numeric(c(x,df[i,29]))
		se <- as.numeric(c(df[i,12],df[i,32]))
		meta <- metagen(beta,se,common=gs("common"))
		cat(cat(c(as.character(df[i,]), meta$TE.common, meta$seTE.common, meta$pval.common,"meta",df[i,27]),sep="\t"),"\n")
		}
	}

	if(!grepl("nocommon",df[i,11]) && grepl("nocommon",df[i,29])){
		# update 04072023; convert beta as the beta of  alternative allele
		if(df[i,7] == df[i,6]){
			cat(cat(c(as.character(df[i,]), df[i,11], df[i,12], df[i,14],"PMID35501419",df[i,6]),sep="\t"),"\n")
		}
		if(df[i,7] != df[i,6]){
			x <- as.numeric(df[i,11])*(-1)
			cat(cat(c(as.character(df[i,]), x, df[i,12], df[i,14],"PMID35501419",df[i,6]),sep="\t"),"\n")
		}
	}

	if(grepl("nocommon",df[i,11]) && ! grepl("nocommon",df[i,29])){	
		cat(cat(c(as.character(df[i,]), df[i,29], df[i,32], df[i,30],"pQTL_cis_NG",df[i,27]),sep="\t"),"\n")
	}
}	
sink()