# -----------------------------------------------------------------------------
#
# Parse Pfam2GO domain annotations
#
# Jason Jiang - Created: 2022/03/23
#               Last edited: 2022/03/23
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Format text file of Pfam2GO Pfam domain GO annotations from
#       http://current.geneontology.org/ontology/external2go/pfam2go as a tidy
#       dataframe.
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = filepath to text file of pfam2go data
  #   $2 = filepath to directory to save formatted pfam2go data in
  # ---------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  pfam2go_df <- format_pfam2go(args[1])
  out <- args[2]
  
  write_csv(pfam2go_df, str_c(out, '/pfam2go.csv'))
}

################################################################################

format_pfam2go <- function(pfam2go_txt) {
  # ---------------------------------------------------------------------------
  # Loads in pfam2go text data and formats as a tidy dataframe.
  # ---------------------------------------------------------------------------
  pfam2go_df <- data.frame(pfam_ID = character(), pfam_dom = character(),
                           GO_term = character(), GO_ID = character())
  
  pfam2go_lines <- readLines(pfam2go_txt)
  
  # parsing in each line of the pfam2go file is messy af because each field uses
  # a different delimiter :/
  for (line in pfam2go_lines) {
    pfam_ID <- str_remove(str_split(line, ' ')[[1]][1], 'Pfam:')
    pfam_dom <- str_split(line, ' ')[[1]][2]
    GO_term <- str_remove(str_remove(str_split(line, ' > ')[[1]][2], ' ; .+'), 'GO:')
    GO_ID <- str_split(line, ' ; ')[[1]][2]
    
    pfam2go_df[nrow(pfam2go_df) + 1,] = c(pfam_ID, pfam_dom, GO_term, GO_ID)
  }
  
  return(pfam2go_df)
}