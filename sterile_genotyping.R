## INITIALIZATION ##
# Load required packages
library(tidyverse)
# Define input / output files
args <- commandArgs(trailingOnly=TRUE)
high_eff_file <- args[1]
gene_list_file <- args[2]
high_eff_annot_file <- args[3]
bed_file <- args[4]
indivs_file <- args[5]

## MERGE DATAFRAMES ##
# Load data
gene_list <- read.table(gene_list_file, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
high_eff <- read.table(high_eff_file, header = TRUE, sep = "\t", 
                    stringsAsFactors = FALSE)
# Remove "chr" column (duplicate)
gene_list <- gene_list[,2:9]
# Remove non-standard genotypes
high_eff <- high_eff[, -c(grep("LIS", names(high_eff)))]
high_eff <- high_eff[, -c(grep("Female", names(high_eff)))]
# Split "effects" column for variants which effect multiple genes on ";"
high_eff <- high_eff %>% separate(effect.s., c("eff1","eff2","eff3","eff4",NA), 
                              sep = ";", fill = "right")
# Replace empty strings resulting from tidyr separate with NA
high_eff[high_eff == ""] <- NA
# Split each effect and gene column on "-"
high_eff <- high_eff %>% separate(eff1, c("eff1","gene1"), 
                                  sep = "-", fill = "right")
high_eff <- high_eff %>% separate(eff2, c("eff2","gene2"), 
                                  sep = "-", fill = "right")
high_eff <- high_eff %>% separate(eff3, c("eff3","gene3"), 
                                  sep = "-", fill = "right")
high_eff <- high_eff %>% separate(eff4, c("eff4","gene4"), 
                                  sep = "-", fill = "right")
# Keep only high_eff genes in gene.list
for (i in 1:(dim(high_eff)[1]))
{
  if (!high_eff$gene1[i] %in% unique(gene_list$gene_ID))
  {
    high_eff$gene1[i] <- NA
    high_eff$eff1[i] <- NA
  }
  if (!high_eff$gene2[i] %in% unique(gene_list$gene_ID))
  {
    high_eff$gene2[i] <- NA
    high_eff$eff2[i] <- NA
  }
  if (!high_eff$gene3[i] %in% unique(gene_list$gene_ID))
  {
    high_eff$gene3[i] <- NA
    high_eff$eff3[i] <- NA
  }
  if (!high_eff$gene4[i] %in% unique(gene_list$gene_ID))
  {
    high_eff$gene4[i] <- NA
    high_eff$eff4[i] <- NA
  }
}
# Collapse effects and genes into one column, filtering out NA values
col_ids <- sort(c(grep(pattern = "eff\\d", names(high_eff)), 
                  grep(pattern = "gene\\d", names(high_eff))))
high_eff$effect <- apply(high_eff[, col_ids], 
                         1, function(x) paste(x[!is.na(x)], collapse = "-"))
high_eff <- high_eff[,-c(col_ids)]
high_eff <- high_eff %>% separate(effect, c("effect","gene_ID"), 
                                  sep = "-", fill = "right")
# Annotate high eff variant table with gene list using merge
high_eff_annot <- merge(gene_list, high_eff, by = "gene_ID", all.y = TRUE)
# Add a flag column to filter for unique variants (needed because VCF is decomposed)
high_eff_annot$flag <- paste(high_eff_annot$POS, high_eff_annot$REF, 
                             high_eff_annot$ALT, sep = "_")
# Add column to quantify amount of support for a variant's gene identity
high_eff_annot <- add_column(high_eff_annot, 
                             Support = rep(1, times = dim(high_eff_annot)[1]),
                             .after = "E_value")
# Collapse variants with multiple sources by combining "Source" and "E_value" columns
for (i in unique(high_eff_annot$flag[duplicated(high_eff_annot$flag)]))
{
  high_eff_annot[which(high_eff_annot$flag == i)[1], "Support"] <- 
    length(which(high_eff_annot$flag == i))
  high_eff_annot[which(high_eff_annot$flag == i)[1], "Source"] <- 
    paste(high_eff_annot[which(high_eff_annot$flag == i), "Source"], collapse = ";")
  high_eff_annot[which(high_eff_annot$flag == i)[1], "E_value"] <- 
    paste(high_eff_annot[which(high_eff_annot$flag == i), "E_value"], collapse = ";")
  high_eff_annot <- high_eff_annot[-c(which(high_eff_annot$flag == i)[-c(1)]),]
}
# Delete flag column
high_eff_annot$flag <- NULL
# Subset only rows where at least one female and one male have a "yes"
females <- high_eff_annot[, grep("FG", names(high_eff_annot))]
males <- high_eff_annot[, grep("MG", names(high_eff_annot))]
rows <- list()
for (row in 1:dim(females)[1])
{
  fems <- which(females[row,]=="yes")
  mens <- which(males[row,]=="yes")
  if (identical(fems, integer(0)) || identical(mens, integer(0)))
  {
    next
  } else {
    rows <- append(rows, row)
  }
}
rows <- unlist(rows)
high_eff_annot <- high_eff_annot[rows,]
# Keep only genotypes that have a "yes" in at least one variant
genos <- c(grep("FG", names(high_eff_annot)), grep("MG", names(high_eff_annot)))
high_eff_annot <- cbind(high_eff_annot[, -c(genos)],
                        high_eff_annot[,names(high_eff_annot 
                                              %>% select_if(~any(. == "yes", na.rm = TRUE)))])


## EXTRACT GENE REGION INFO ##
# BED-style dataframe of regions ~500 bp around remaining genes ##
bed <- paste(high_eff_annot$CHROM, high_eff_annot$Start, 
             high_eff_annot$End, high_eff_annot$gene_ID, sep = "-")
bed <- data.frame(bed = unique(bed))
bed <- bed %>% separate(bed, c("#CHR","START","END","GENE"), sep = "-")
bed$START <- (as.numeric(bed$START) - 500)
bed$END <- (as.numeric(bed$END) + 500)
bed[bed < 0] <- 0
## EXTRACT INDIV. GENOTYPE FOR EACH GENE ##
indivs <- data.frame(GENE = character(), INDIV = character())
for (row in 1:dim(high_eff_annot)[1])
{
  indices <- which(high_eff_annot[row,]=="yes")
  ids <- gsub(x = names(high_eff_annot[,indices]), pattern = "\\.", replacement = "-")
  gene <- rep(high_eff_annot[row, "gene_ID"], length(ids))
  temp <- data.frame(GENE = gene, INDIV = ids)
  indivs <- rbind(indivs, temp)
}
indivs <- unique(indivs)


## EXPORT TO FILES ##
# Write BED file
write.table(bed, bed_file, row.names = FALSE, quote = FALSE, sep = "\t")
# Write dataframe of genes with individual genotypes for secondary validation
write.table(indivs, indivs_file, col.names = FALSE, 
            row.names = FALSE, quote = FALSE, sep = "\t")
# Clean up high effect dataframe with tidyverse
high_eff_annot <- relocate(high_eff_annot, "effect", .after = "protein_product")
high_eff_annot <- arrange(high_eff_annot, desc(Support), desc(E_value))
# Fix individual IDs in header
names(high_eff_annot) <- gsub(x = names(high_eff_annot), 
                              pattern = "\\.", replacement = "-")
write.table(high_eff_annot, high_eff_annot_file, row.names = FALSE,
            quote = FALSE, sep = "\t")
