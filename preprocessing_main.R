suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(edgeR)
  library(qusage)
  library(fs)
})

# source helper function script
get_script_dir <- function() {
  commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(
      col = value, into = c("key", "value"), sep = "=", fill = "right"
    ) %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value) %>%
    dirname(.)
}
source(file.path(get_script_dir(), "preprocessing_helper_functions.R"))

# define option list 
option_list <- list(
  make_option("--output_dir", "-o",
    dest="output_dir", type="character",
    help="output directory where files will be saved"),
  make_option("--name",
    type = "character",
    help = "name of the dataset"),
  make_option("--data.counts",
    dest="data_counts", type="character",
    help="input directory where datasets and genesets are found"),
  make_option("--data.metadata",
    dest="metadata", type="character",
    help="input directory where datasets and genesets are found"),
  make_option("--data.genesets",
    dest="genesets", type="character",
    help="input directory where datasets and genesets are found"),
  make_option("--data.referencecounts",
    dest="reference_counts", type="character",
    help="input directory where datasets and genesets are found")
)
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)
analysis_name <- opts$name

# metadata checkpoint
metadata_df <- read.csv(opts$metadata, sep = '\t', check.names = FALSE)
if (!'annotation' %in% names(metadata_df)){
  stop('Missing required "annotation" column in metadata table')
}
if (!'true_label' %in% names(metadata_df)){
  stop('Missing required "true_label" column in metadata table')
}
metadata_df$annotation <- ifelse(metadata_df$annotation == '', 'Unknown', metadata_df$annotation)

# Read data and filter genes as well as samples
input_data <- counts_filtering(opts$data_counts, metadata_df)
counts_filt <- input_data$df
metadata_df <- input_data$meta
reference_filt <- counts_filtering(opts$reference_counts, '', reference = TRUE)
common_genes <- intersect(rownames(counts_filt), rownames(reference_filt))
input_ranked <- rank_dataframe(counts_filt[common_genes, ])
reference_ranked <- rank_dataframe(reference_filt[common_genes, ])

reference_ranked_mean <- data.frame( # use for rankdiff
  centroid = apply(reference_ranked, 1, mean)
)
rownames(reference_ranked_mean) <- rownames(reference_ranked)

# Get rank difference and then re rank it for consistency downstream with non rank-diffed data
rank_diff_df <- input_ranked - reference_ranked_mean$centroid
rank_diff_df_ranked <- as.data.frame(apply(rank_diff_df, 2, function(col) rank(col, ties.method = "min")))

if (!dir.exists(opts$output_dir)) {
  dir.create(opts$output_dir, recursive = TRUE)
}

write.table(rownames_to_column(input_ranked, var='Geneid'), path(opts$output_dir, paste0(analysis_name, '-sampleranks.tsv')), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(rownames_to_column(rank_diff_df_ranked, var='Geneid'), path(opts$output_dir, paste0(analysis_name, '-samplediffranks.tsv')), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(rownames_to_column(reference_ranked, var='Geneid'), path(opts$output_dir, paste0(analysis_name, '-referenceranks.tsv')), sep = '\t', quote = FALSE, row.names = FALSE)

write.table(metadata_df, path(opts$output_dir, paste0(analysis_name, '-metadata_preprocessed.tsv')), sep = '\t', quote = FALSE, row.names = FALSE)


# Copy GeneSet file to the tmp dir for input
print(paste0('Fetching ', analysis_name, ' gmt file for input'))
filter_geneset_file(opts$genesets, path(opts$output_dir, paste0(opts$name, "-GenesetsOI_preprocessed.gmt")), row.names(input_ranked))

print('-----------✅ Input preprocessing complete ✅---------------')
