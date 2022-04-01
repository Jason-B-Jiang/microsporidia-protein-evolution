# -----------------------------------------------------------------------------
#
# Parse yeast gene intron coordinates
#
# Jason Jiang - Created: 2022/03/07
#               Last edited: 2022/03/07
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Parse in SGD gene feature text files for yeast genes with introns, and
#       output a csv file listing intron coordinates for these yeast genes.
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))

################################################################################

main <- function() {
  # ----------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = filepath to directory containing SGD gene feature text files for all
  #        yeast genes w/ introns
  #
  #   $2 = filepath to save csv file of yeast gene intron coordinates to
  # ----------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  intron_coordinates <- get_intron_coordinates(args[1])
  write_csv(intron_coordinates, args[2])
}

################################################################################

### Helper functions

get_intron_coordinates <- function(yeast_introns_dir) {
  # ----------------------------------------------------------------------------
  # Parses SGD gene feature files into a dataframe of intron coordinates for each
  # yeast gene with introns.
  #
  # Arguments:
  #   yeast_introns_dir = filepath to directory with SGD gene feature files
  # ----------------------------------------------------------------------------
  # NOTE: this code block produces a lot of warnings because the SGD feature files
  # are formatted a bit weirdly
  intron_coordinates <- Sys.glob(str_c(yeast_introns_dir, '/*.txt')) %>%
    as.list() %>%
    # load in each SGD gene feature file as a dataframe
    lapply(function(x) {read_table(x, comment = '!', show_col_types = F)}) %>%
    # combine all dataframes in the list into a single dataframe
    bind_rows() %>%
    # only select 2nd - 5th columns (SGD gene name, gene feature,
    # gene feature coordinates and location of gene in chromosome)
    select(2, 3, 4, 5) %>%
    # rename columns to more descriptive names
    rename(sgd_name = Feature_1, feature = Systematic, coordinates = Name,
           genome_coordinates = Feature_2) %>%
    rowwise() %>%
    mutate(genome_coordinates = str_split(genome_coordinates, ':')[[1]][2])
  
  return(intron_coordinates)
}

################################################################################

main()