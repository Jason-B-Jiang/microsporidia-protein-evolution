# -----------------------------------------------------------------------------
#
# Make microsporidia clade hashtable
#
# Jason Jiang - Created: 2022/01/26
#               Last edited: 2022/01/26
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Take in a spreadsheet labelling each microsporidia species by their
#       evolutionary clades, and output a hashtable mapping each species
#       to their respective clade.
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))

################################################################################

main <- function() {
  # $1 = filepath to csv file of microsporidia to clade labels
  # $2 = filepath + filename of where we want to save the clade hashtable
  args <- commandArgs(trailingOnly = T)
  clades <- read_csv(args[1], show_col_types = F)
  outfile <- args[2]
  
  clade_hash <- new.env()
  for (i in 1 : nrow(clades)) {
    species <- clades$species[i]
    clade <- clades$clade[i]
    
    # set microsporidia species name as key, and clade as value
    clade_hash[[species]] <- clade
  }
  
  saveRDS(clade_hash, outfile)
}

################################################################################

main()