# -----------------------------------------------------------------------------
#
# Plot microsporidia protein lengths vs. %identity with yeast orthologs
#
# Jason Jiang - Created: 2022/02/02
#               Last edited: 2022/03/14
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Create scatterplots comparing lengths of microsporidia proteins to their
#       percent identity with their yeast orthologs, for each individual
#       microsporidia species and all species together.
#
# Thanks to Brandon Murareanu for this RScript layout
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(seqinr))

################################################################################

main <- function() {
  # ----------------------------------------------------------------------------
  # Arguments taken:
  #   $1 = filepath to csv of single-copy orthogroups + ortholog lengths,
  #        generated while making Fig 1B
  #
  #   $2 = filepath to directory with pairwise ortholog alignments for all
  #        single-copy ortholog pairs, for all microsporidia species w/ yeast
  #
  #   $3 = filepath to directory to save outputs in
  # ----------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  sco_lengths <- read_csv(args[1], show_col_types = F)
  ortho_alignments_dir <- args[2]
  out_dir <- args[3]
  
  # Calculate %identity (pid) for each single-copy ortholog pair, across all
  # microsporidia species
  sco_lengths <- sco_lengths %>%
    mutate(pid = calculate_percent_id(name, yeast, ortho_alignments_dir)) %>%
    # Split sco_lengths into list of dataframes, with each dataframe holding
    # ortholog pair length data + %identity for 1 species only
    split(sco_lengths$name)
  
  # draw plots for each individual microsporidia and save plots to out_dir
  Map(plot_percent_identities, sco_lengths, names(sco_lengths),
      rep(out_dir, length(sco_lengths)))
  
  # draw plot for for all microsporidia species together and save plot to
  # out_dir
  all_species <- bind_rows(sco_lengths)  # merge all dataframes into 1
  plot_percent_identities(all_species, 'All species', out_dir, TRUE)
  
  # write csv of all microsporidia species and their %conservation with yeast
  # orthologs
  write_csv(all_species, str_c(out_dir, '/ortholog_percent_identities.csv'))
}

################################################################################

### Helper functions

calculate_percent_id <- function(species, yeast, ortho_alignments_dir) {
  # ----------------------------------------------------------------------------
  # Returns %identity between a microsporidia - yeast ortholog pair, specified
  # by microsporidia species, yeast ortholog name and directory
  # with pairwise ortholog alignments.
  #
  # Percent identity is calculated as
  #   100 * [(identical aligned residues) / (number of aligned residues)]
  # ----------------------------------------------------------------------------
  ortho_alignment <- read.fasta(file = str_c(ortho_alignments_dir, species,
                                             str_c(yeast, '.fa'), sep = '/'),
                                seqtype = 'AA', as.string = T)
  
  # NOTE: aligned sequences should be of same length due to introduced gaps
  s1 <- str_split(toString(ortho_alignment[[1]]), '')[[1]]
  s2 <- str_split(toString(ortho_alignment[[2]]), '')[[1]]
  
  # Identical residues = alignment positions with exact same amino acid
  # Aligned residues = residues in both sequences not aligned with any gaps in
  # the other sequence
  return(100 * (sum(s1 == s2) / sum(s1 != '-' & s2 != '-')))
}

calculate_percent_id <- Vectorize(calculate_percent_id,
                                  vectorize.args = c('species', 'yeast',
                                                     'ortho_alignments_dir'))


plot_percent_identities <- function(sco_df, name, out, many_points=FALSE) {
  # ----------------------------------------------------------------------------
  # Makes a scatterplot of microsporidia protein lengths to their %identity
  # with their yeast orthologs.
  #
  # Arguments:
  #   sco_df = dataframe of orthogroup data, with %id between ortholog pairs
  #
  #   name = name of microsporidia species
  #
  #   out = output directory to save plot in
  #
  #   many_points = boolean indicating if we want to adjust the scatterplot to
  #                 accomodate many points
  # ----------------------------------------------------------------------------
  if (many_points) {
    # decrease opacity of points to accommodate a larger number of points
    # to plot
    pid_plot <- ggplot(data = sco_df, aes(x = microsp_yeast_ratio, y = pid)) +
      geom_point(alpha = 0.1) +
      geom_rug(alpha = 0.1) +
      labs(x = 'Relative length of microsporidia ortholog to yeast', y = '%Identity with yeast',
           title = str_c(name, '\n', str_c('n = ', nrow(sco_df), ' ortholog pairs'))) +
      theme(axis.title = element_text(color = 'black', size = 18),
            axis.text = element_text(size = 14),
            title = element_text(color = 'black', size = 24)) +
      theme_bw()
      
  } else {
    pid_plot <- ggplot(data = sco_df, aes(x = microsp_yeast_ratio, y = pid)) +
      geom_point(alpha = 0.5) +
      geom_rug(alpha = 0.2) +
      labs(x = 'Relative length of microsporidia ortholog to yeast', y = '%Identity with yeast',
           title = str_c(name, '\n', str_c('n = ', nrow(sco_df), ' ortholog pairs'))) +
      theme(axis.title = element_text(color = 'black', size = 18),
            axis.text = element_text(size = 14),
            title = element_text(color = 'black', size = 24)) +
      theme_bw()
  }
  
  ggsave(filename = str_c(out, '/', name, '.tiff'), plot = pid_plot, units = 'in',
         width = 7, height = 7, dpi = 600)
}

################################################################################

main()