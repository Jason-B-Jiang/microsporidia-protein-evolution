# -----------------------------------------------------------------------------
#
# Map SGD to Uniprot names
#
# Jason Jiang - Created: 2022/03/07
#               Last edited: 2022/03/07
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Create a hashtable mapping SGD yeast gene names to their corresponding
#       Uniprot yeast protein names
#
# Thanks to Brandon Murareanu for this RScript layout
#
# NOTE: yeast.txt was downloaded from https://www.uniprot.org/docs/yeast.txt
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))

################################################################################

main <- function() {
  # ----------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = filepath to yeast.txt, a file from Uniprot mapping SGD yeast gene
  #        names to their Uniprot protein names
  #
  #   $2 = filepath to save the hashtable of SGD to Uniprot names to, as an RDS
  #        object
  # ----------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  
  sgd_uniprot_hash <- map_sgd_uniprot_names(args[1])
  saveRDS(sgd_uniprot_hash, args[2])
}

################################################################################

### Helper functions

map_sgd_uniprot_names <- function(yeast_names) {
  # ----------------------------------------------------------------------------
  # Parses yeast.txt (yeast_names) and returns hashtable mapping SGD yeast gene
  # names to Uniprot protein names.
  #
  # Arguments:
  #   yeast_names: filepath to yeast.txt file
  # ----------------------------------------------------------------------------
  yeast_names <- readLines(yeast_names) %>%
    as.list() %>%
    # turn each row in yeast.txt into a list element
    # each list element is a size 2 vector, with first vector element being SGD
    # name for yeast gene, and second vector element being Uniprot protein name
    lapply(function(x) {str_split(substr(x[1], 76, str_length(x[1])), "\\s+")[[1]][c(1, 2)]})
  
  # Now turn information from this list into a hashtable, with keys being SGD
  # gene names and values being Uniprot names
  yeast_hash <- new.env()
  for (pair in yeast_names) {
    yeast_hash[[pair[1]]] <- pair[2]
  }
  
  return(yeast_hash)
}

################################################################################

main()