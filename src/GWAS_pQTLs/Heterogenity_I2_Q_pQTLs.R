library(metafor)
library(dplyr)
library(data.table)
args = commandArgs(trailingOnly = TRUE)

meta_files=c("aa","ab","ac","ad","ae","af","ag","ah","ai","aj","ak","al","am","an",
            "ao","ap","aq","ar","as","at","au","av","aw","ax","ay","az",
            "ba","bb","bc","bd","be","bf","bg","bh","bi")

 
# Function: run FE meta for two studies and return Q
rma_two <- function(b1, se1, b2, se2) {
  m <- rma(yi = c(b1, b2), sei = c(se1, se2), method = "FE")
  c(
    beta_FE = as.numeric(coef(m)),
    se_FE   = sqrt(vcov(m)),
    z_FE    = as.numeric(coef(m))/sqrt(vcov(m)),
    p_FE    = 2*pnorm(abs(as.numeric(coef(m))/sqrt(vcov(m))), lower.tail = FALSE),
    Q       = as.numeric(m$QE),      # Cochran's Q
    Q_p     = as.numeric(m$QEp),     # p-value for Q (df = 1 when k=2)
  )
}

meta_file=meta_files[as.numeric(args[1])]
meta_df = fread(paste0("/meta_res_hg38_rmnocommon/x", meta_file))
# drop rows containing "nocommon" in each column
set1 <- unique(meta_df$entrezgenesymbol[!grepl("nocommon", meta_df$entrezgenesymbol)])
set2 <- unique(meta_df$Gene_Name[!grepl("nocommon", meta_df$Gene_Name)])
# union of two character vectors
GeneNames <- union(set1, set2)

# length of the union
print(length(GeneNames))
for(index in 1:length(GeneNames)){
	Gene = GeneNames[index]
	meta_df_per_gene=meta_df[(meta_df$entrezgenesymbol==Gene) & (meta_df$Tag=="meta"), ]  #pQTL from both datasets
	#hamonize the sign of beta based on effect allele
	# Subset only the needed columns
	pqtls_two_std <- meta_df_per_gene[, c("ID", "A1", "BETA", "SE", "effectAllele", "Beta", "SE.1", "Counted_allele_in_regression")]
	pqtls_two_std$BETA <- as.numeric(pqtls_two_std$BETA)
	pqtls_two_std$SE   <- as.numeric(pqtls_two_std$SE)
	pqtls_two_std$Beta <- as.numeric(pqtls_two_std$Beta)
	pqtls_two_std$SE.1 <- as.numeric(pqtls_two_std$SE.1)
	# Study 1 aligned to Counted_allele_in_regression
	pqtls_two_std$beta1_count_allele <- ifelse(
	  pqtls_two_std$A1 == pqtls_two_std$Counted_allele_in_regression,
	  pqtls_two_std$BETA,
	  -pqtls_two_std$BETA
	)
	pqtls_two_std$se1_count_allele <- pqtls_two_std$SE
	# Study 2 aligned to Counted_allele_in_regression
	pqtls_two_std$beta2_count_allele <- ifelse(
	  pqtls_two_std$effectAllele == pqtls_two_std$Counted_allele_in_regression,
	  pqtls_two_std$Beta,
	  -pqtls_two_std$Beta
	)
	pqtls_two_std$se2_count_allele <- pqtls_two_std$SE.1

	# Apply row-wise across all SNPs
	res_mat <- mapply(
	  rma_two,
	  pqtls_two_std$beta1_count_allele, pqtls_two_std$se1_count_allele,
	  pqtls_two_std$beta2_count_allele, pqtls_two_std$se2_count_allele
	)
	# Bind results back (transpose mapply output)
	res_df <- as.data.frame(t(res_mat))
	pqtls_two_std <- cbind(pqtls_two_std, res_df)
	fwrite(pqtls_two_std, paste0("/meta_res_hg38_rmnocommon/Cochran_Q/perProtein/x", meta_file,"_",Gene,".csv"))
}



