#!/usr/bin/Rscript

# -----------------------------------------------------------------------------
#
# Make Pfam family id:clan id has table
#
# Jason Jiang - Created: 2022/01/14
#               Last edited: 2022/01/14
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Create and save a hashtable of Pfam families to their corresponding
#       Pfam clans
#
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

library(readr)


main <- function() {
  # load in Pfam-A.clans.tsv spreadsheet as dataframe
  args <- commandArgs(trailingOnly=TRUE)
  clan_data <- args[1]  # tsv mapping families to clans
  output <- args[2]  # name to save resulting hashtable as
  
  clan_df <- read_tsv(clan_data)
  
  # initialize hash table for mapping pfam families to pfam clans
  fam_to_clan_table <- new.env()
  
  # set keys as family ids, and values as clan ids in fam_to_clan_table
  for (i in 1 : nrow(clan_df)) {
    family_id <- as.character(clan_df[i, 'Family_ID'])
    clan_id <- as.character(clan_df[i, 'Clan_name'])
    fam_to_clan_table[[family_id]] <- clan_id
  }
  
  saveRDS(fam_to_clan_table, file = output)
}


main()