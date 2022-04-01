# -----------------------------------------------------------------------------
#
# Join aligned domain architecture dataframes
#
# Jason Jiang - Created: 2022/01/25
#               Last edited: 2022/01/25
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Combine dataframes of domain architecture alignments for each
#       microsporidia species.
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))

################################################################################

main <- function() {
  # -----------------------------------------------------------------------------
  # $1 = domain arch folder
  # $2 = list of essential yeast genes
  # $3 = output folder
  # $4 - $n = microsporidia species names
  # -----------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = TRUE)
  DA_dir <- args[1]
  essential_list <- readLines(args[2])
  output_dir <- args[3]
  microsp_names <- args[4 : length(args)]
  
  merged_aligned_DA <- merge_DA_dfs(DA_dir, essential_list, microsp_names)
  
  outfile = '/merged_DA_alignments.csv'
  write_csv(merged_aligned_DA, str_c(output_dir, outfile))
}

################################################################################

# Helper functions

merge_DA_dfs <- function(DA_dir, essential_list, microsp_names) {
  # ----------------------------------------------------------------------------
  # Collects all csv files of aligned domain architectures for each microsporidia
  # species in DA_dir, and combines these csv files into one dataframe.
  #
  # DA_dir: directory with csv files of aligned domain architectures
  # essential_list: character vector of the uniprot protein names for essential
  #                 genes in yeast
  # microsp_names: character vector of abbreviated microsporidia names, used in
  # `              files in DA_dir`
  # ----------------------------------------------------------------------------
  # Get filepaths to all domain arch csv files
  DA_files <- sapply(microsp_names,
                     function(x) {str_c(DA_dir, '/', x, '_S_cere/', x,'_S_cere_DA_aligned.csv')})
  
  # initialize dataframe for accumulating all of these csvs into
  merged_DA_df <- initialize_merged_df(DA_files[1])
  
  # loop over all domain arch csv files and add to merged_DA_df
  for (i in 1 : length(DA_files)) {
    DA_df <- read_csv(DA_files[i], show_col_types = F)
    
    DA_df <- DA_df %>%
      mutate(species = microsp_names[i],
             essential = yeast %in% essential_list,
             num_conserved_DA = sum(is.na(DA_df$DA_diffs)),
             DA_conservation_rate = num_conserved_DA / nrow(DA_df),
             num_lost_doms_DA = sum(!is.na(DA_df$lost_doms)),
             dom_loss_rate = num_lost_doms_DA / nrow(DA_df)) %>%
      # arrange columns so they're in same order as in merged_DA_df
      select(species, essential, everything())
    
    merged_DA_df <- rbind(merged_DA_df, DA_df)
  }
  
  return(merged_DA_df)
}


initialize_merged_df <- function(DA_file) {
  # ----------------------------------------------------------------------------
  # Returns an empty dataframe for accumulating domain architecture csv files
  # into.
  #
  # DA_file: filepath to any domain architecture csv file, to use for initializing
  #          size + column names of the empty dataframe.
  # ----------------------------------------------------------------------------
  example_df <- read_csv(DA_file, show_col_types = F)
  
  col_names <- c('species', 'essential', colnames(example_df),
                 'num_conserved_DA', 'DA_conservation_rate', 'num_lost_doms_DA',
                 'dom_loss_rate')
  
  # add 6 to ncol to account for new columns to add (species, essential, etc)
  merged_DA_df <- data.frame(matrix(nrow = 0, ncol = ncol(example_df) + 6))
  colnames(merged_DA_df) <- col_names
  
  return(merged_DA_df)
}

################################################################################

main()  # run this script