# -----------------------------------------------------------------------------
#
# Get orthogroup protein summary statistics
#
# Jason Jiang - Created: 2022/01/19
#               Last edited: 2022/01/24
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: For each microsporidia-yeast pair passed in, write a new csv file
#       with the median ortholog lengths between each pair, median ortholog
#       length differences and signficance of length difference.
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))

################################################################################

main <- function() {
  # $1 = output directory for summary stats csvs
  # $2 - $n = all csvs of orthogroups to process
  args <- commandArgs(trailingOnly = T)
  output_dir <- args[1]
  species_pairs <- args[2 : length(args)]
  microsp_names = str_extract(species_pairs, '[:upper:]_[:lower:]{2,4}\\.?')
  
  prt_len_df <- get_prt_len_stats(species_pairs, microsp_names)
  write_csv(prt_len_df, file = str_c(output_dir, '/SCO_protein_lengths_summary.csv'))
}

################################################################################

# Helper functions

get_prt_len_stats <- function(species_pairs, microsp_names) {
  # initialize empty dataframe to accumulate rows from each species pair
  combined_df <- data.frame(matrix(ncol = 10, nrow = 0))
  colnames(combined_df) <- c('species', 'num_orthogroups', 'orthogroup',
                             'microsp_ortho', 'yeast_ortho', 'microsp_len',
                             'yeast_len', 'microsp_yeast_ratio',
                             'species_median_ratio', 'p_species')
  
  # initialize vector to hold all individual microsp/yeast ortholog lengths
  for (i in 1 : length(species_pairs)) {
    pair = species_pairs[i]
    name = microsp_names[i]
    combined_df <- rbind(combined_df, get_species_summary(pair, name))
  }
  
  # apply bonferroni p-value correction for multiple testing, multiplying
  # each p-value by the number of species tested
  # calculate overall p-value for median difference between microsporidia and
  # yeast ortholog lengths
  combined_df <- combined_df %>%
    mutate(p_species_corrected = p_species * length(species_pairs),
           p_overall = wilcox.test(microsp_len, yeast_len, paired = TRUE,
                                   alternative = 'two.sided')$p.value)

  return(combined_df)
}


get_species_summary <- function(species_pair, microsp_name) {
  # Load in csv for species_pair
  # rename orthogroup, microsp, S_cere column, microsp_len column, S_cere_len
  # column
  # add in species column, num_orthogroups column, p_species column
  sco_csv <- read_csv(species_pair, show_col_types = FALSE)
  
  # rename columns named with the microsporidia species name to more generic
  # names
  names(sco_csv)[names(sco_csv) == microsp_name] <- 'microsp_ortho'
  names(sco_csv)[names(sco_csv) == str_c(microsp_name, '_len')] <- 'microsp_len'
  
  sco_csv <- sco_csv %>%
    rename(orthogroup = Orthogroup,
           yeast_ortho = S_cere,
           yeast_len = S_cere_len) %>%
    mutate(species = microsp_name,
           num_orthogroups = nrow(sco_csv),
           microsp_yeast_ratio = microsp_len / yeast_len,
           species_median_ratio = median(microsp_yeast_ratio),
           p_species = wilcox.test(microsp_len, yeast_len, paired = TRUE,
                                    alternative = 'two.sided')$p.value) %>%
    select(species, num_orthogroups, orthogroup, microsp_ortho, yeast_ortho,
           microsp_len, yeast_len, microsp_yeast_ratio, species_median_ratio,
           p_species)
  
  return(sco_csv)
}

################################################################################

main()