library(data.table)
library("remotes")
library("readxl")
library(openxlsx)

library("enrichR")
listEnrichrSites()
setEnrichrSite("Enrichr")

dbs <- listEnrichrDbs()

# 101 risk proteins 
df <- read_excel("/AJHG_revision/SupplementaryTables_20250903.xlsx", sheet = "S6", skip = 1)
genes <- unique(df[["Gene...2"]])
enriched <- enrichr(genes, dbs$libraryName)

header=c("Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value", "Old.Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes")
# Initialize empty list to collect filtered data frames
filtered_list <- list()
# Loop through each enrichment result
for (sheet_name in names(enriched)) {
  df <- enriched[[sheet_name]]
  # Check if required column exists
  if ("Adjusted.P.value" %in% colnames(df)) {
    df_sig <- df[df$Adjusted.P.value < 0.1, ]
    
    # Proceed only if non-empty
    if (nrow(df_sig) > 0) {
      df_sig$Source <- sheet_name
      filtered_list[[sheet_name]] <- df_sig
    }
  }
}

# Combine all filtered data frames into one
if (length(filtered_list) > 0) {
  combined_df <- do.call(rbind, filtered_list)
  # Reorder columns if needed
  combined_df <- combined_df[, c(header, "Source")]
  # Save to Excel
  write.xlsx(combined_df, 
             file = "/AJHG_revision/EnrichR/enriched_genes.xlsx", 
             rowNames = FALSE)
} else {
  message("No significant results (Adjusted.P.value < 0.05) found in any enrichment set.")
}

## Stratify by cancers
df_all <- read_excel("/AJHG_revision/SupplementaryTables_20250903_QY_LQ.xlsx", sheet = "S6", skip = 1)
for(cancer in unique(df_all$Cancer)){
	df = df_all[df_all["Cancer"]==cancer, ]
	genes <- unique(df[["Gene"]])
	enriched <- enrichr(genes, dbs$libraryName)

	header=c("Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value", "Old.Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes")
	# Initialize empty list to collect filtered data frames
	filtered_list <- list()
	# Loop through each enrichment result
	for (sheet_name in names(enriched)) {
	  df <- enriched[[sheet_name]]
	  # Check if required column exists
	  if ("Adjusted.P.value" %in% colnames(df)) {
		df_sig <- df[df$Adjusted.P.value < 0.1, ]
		
		# Proceed only if non-empty
		if (nrow(df_sig) > 0) {
		  df_sig$Source <- sheet_name
		  filtered_list[[sheet_name]] <- df_sig
		}
	  }
	}

	# Combine all filtered data frames into one
	if (length(filtered_list) > 0) {
	  combined_df <- do.call(rbind, filtered_list)
	  # Reorder columns if needed
	  combined_df <- combined_df[, c(header, "Source")]
	  # Save to Excel
	  write.xlsx(combined_df, 
				 file = paste0("/AJHG_revision/EnrichR/enriched_genes_",cancer,".xlsx"), 
				 rowNames = FALSE)
	} else {
	  message("No significant results (Adjusted.P.value < 0.05) found in any enrichment set.")
	}
}
