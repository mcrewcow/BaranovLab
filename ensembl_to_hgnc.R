#Seurat utils based

ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get the Ensembl gene IDs from your Seurat object
ensembl_gene_ids <- rownames(hippo_sten)


gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters = "ensembl_gene_id",
                   values = ensembl_gene_ids,
                   mart = ensembl)

# Create a mapping data frame with NA values for all HGNC symbols
gene_mapping <- data.frame(Ensembl_Gene_ID = ensembl_gene_ids, HGNC_Symbol = NA)

# Get the HGNC symbols for Ensembl gene IDs that have mappings
matching_indices <- match(gene_info$ensembl_gene_id, ensembl_gene_ids)
gene_mapping$HGNC_Symbol[matching_indices] <- gene_info$hgnc_symbol

hippo_sten <- Seurat.utils::RenameGenesSeurat(hippo_sten, newnames = gene_mapping$HGNC_Symbol)
