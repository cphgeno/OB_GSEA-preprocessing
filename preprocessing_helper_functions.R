filter_geneset_file <- function(geneset.input, geneset.output, genes_in_dataset){
  # function to fetch geneset file and filter out the genes not in the pre-processed input counts
  print('--------Reading geneset list------------')

  if (!file.exists(geneset.input)){
    stop(paste('Cannot find GeneSet input file:', geneset.input))
  }

  # Parse gene sets and only keep those present in pre-processed dataset
  genesets <- read.gmt(geneset.input)
  genesets_filt <- list()
  for (geneset in names(genesets)){
    genes <- genesets[[geneset]]
    overlap <- genes[genes %in% genes_in_dataset]
    if (length(overlap) >= 10) {
      genesets_filt[[geneset]] <- overlap
    } else {
      print(paste('GeneSet:', geneset, 'does not contain enough (>=10) genes present in input df ---> Rejected'))
    }
  }

  if (length(genesets_filt) > 0){
    save_gmt(genesets_filt, geneset.output)
  } else {
    stop('All GeneSets rejected! ---> Terminating run')
  }
}


save_gmt <- function(gene_sets, filepath) {
  # gene_sets: a named list where each element is a vector of gene symbols
  # filepath: output .gmt file path
  
  gmt_output <- file(filepath, "wt")
  for (set_name in names(gene_sets)) {
    genes <- gene_sets[[set_name]]
    line <- paste(c(set_name, "NA", genes), collapse = "\t")
    writeLines(line, gmt_output)
  }
  close(gmt_output)
}


filter_genes_by_name <- function(df, min_expression_threshold=1) {
  # function to remove genes that are not protein-coding
  # additionally, it can remove genes based on minimal gene mean expression

  gene_symbols <- rownames(df)
  
  # Define regular expressions for filtering
  mito_genes <- grepl("^MT-", gene_symbols)  # Mitochondrial genes
  ribo_genes <- grepl("^RPL|^RPS", gene_symbols)  # Ribosomal genes
  lncRNA_genes <- grepl("^LINC|^LOC", gene_symbols)  # Long non-coding RNAs
  pseudo_genes <- grepl("-P[0-9]+$", gene_symbols)  # Pseudogenes ending in "-P#"
  small_RNAs <- grepl("^SNORD|^SCARNA|^MIR", gene_symbols)  # Small RNAs like snoRNAs, miRNAs
  
  # Combine all filters
  genes_to_remove <- mito_genes | ribo_genes | lncRNA_genes | pseudo_genes | small_RNAs
  
  # Apply filtering
  filtered_df <- df[!genes_to_remove, ]
  
  # Filter out lowly expressed genes
  # filtered_df <- filtered_df[rowMeans(filtered_df) > min_expression_threshold, ]
  
  return(filtered_df)
}

counts_filtering <- function(counts_path, metadata = '', reference = FALSE){
  # function to convert counts to ranks after filtering out only protein coding and non constant/null values
  counts <- read.csv(counts_path, sep = '\t', check.names = FALSE)
  counts <- aggregate(. ~ Geneid, data = counts, FUN = sum) # sum up duplicated gene names
  counts <- column_to_rownames(counts, 'Geneid')
  # make sure they are all numbers to avoid issues with filtering
  counts[] <- lapply(counts, function(x) if (is.factor(x) | is.character(x)) as.numeric(as.character(x)) else x)
  counts[is.na(counts)] <- 0
  print('df original dims:')
  print(dim(counts))

  # subset counts df if too big to avoid memory limit errors
  if (dim(counts)[2] > 5000) {
    set.seed(219)
    random_cols = sample(colnames(counts), 5000)
    counts <- counts[, random_cols]
    print('df filtered dims:')
    print(dim(counts))
  }

  # intersect samples for which you have metadata for - not for reference
  if (!reference) {
    metadata <- metadata %>% filter(!annotation == '')
    samples_wmetadata <- intersect(metadata$filename, names(counts))
    metadata <- metadata %>% filter(filename %in% samples_wmetadata)
    counts <- counts[names(counts) %in% samples_wmetadata]
    print('df only samples with metadata dims:')
    print(dim(counts))
  }

  counts <- counts[!apply(counts == 0, 1, all), ]
  # remove constant value in all samples
  counts <- counts %>% filter(apply(., 1, n_distinct) > 1)
  print('df filtered - no constant/0 value dims:')
  print(dim(counts))


  # only protein coding genes included
  counts_filt <- filter_genes_by_name(counts)
  print("df filtered - only 'protein' genes dims:")
  print(dim(counts_filt))

  if (!reference) {
    return(list(df = counts_filt, meta = metadata))
  } else {
    return(counts_filt)
  }
}


rank_dataframe <- function(data_to_rank){
    # rank all columns independently in a dataframe
    data_ranked <- data.frame(row.names = row.names(data_to_rank))
    for (samplename in names(data_to_rank)) {
        data_wsample_ranking <- data_to_rank %>%
            mutate(ranking = rank(.[,samplename], ties.method = "min", na.last = 'keep'))
        data_ranked[, samplename] <- as.numeric(data_wsample_ranking$ranking)
    }
    return(data_ranked)
}

