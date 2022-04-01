# -----------------------------------------------------------------------------
#
# Get single copy orthogroups
#
# Jason Jiang - Created: 2022/01/04
#               Last edited: 2022/01/04
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: From OrthoFinder orthogroups for a microsporidia species and yeast,
#       output a tsv with only the single copy orthogroups, with length
#       annotations from each ortholog
#
#
# Thanks to Brandon Murareanu for this RScript layout
#
# ------------------------------------------------------------------------------

library(dplyr)
library(readr)
library(seqinr)
library(tidyr)
library(stringr)


main <- function() {
  orthofinder_results <- commandArgs(trailingOnly = T)
  
  if (length(args) != 1) {
    stop("Usage: [orthofinder orthogroups results]")
  }
  
  sco_tsv <- get_single_copy_orthogroups(orthofinder_results)
  
  # Name output csv file as "../results/{species}_S_cere_SCOs.csv
  output_csv = str_c("../results/",
                     str_extract(orthofinder_results, "[:upper:]_[:lower:]{2,4}\\.?_S_cere"),
                     "_SCOs.csv")
  
  write_csv(sco_tsv, output_csv)
}


get_ortho_lengths <- function(orthogroup, orthofinder_results) {
  # load in fasta file with orthogroup sequences
  og_fa <- read.fasta(file = paste0(orthofinder_results,
                                    '/Single_Copy_Orthologue_Sequences/',
                                    orthogroup, '.fa'))
  
  # return vector of protein lengths for each ortholog
  return(paste0(getLength(og_fa[[1]]), ',', getLength(og_fa[[2]])))
}

get_ortho_lengths <- Vectorize(get_ortho_lengths,
                               vectorize.args = c("orthogroup",
                                                  "orthofinder_results"))


get_single_copy_orthogroups <- function(orthofinder_results) {
  # Load in tsv of orthogroups from orthofinder results, as well as text file
  # of single copy orthogroups
  orthogroups_tsv <- read_tsv(
    paste0(orthofinder_results, "/Orthogroups/Orthogroups.tsv")
    )
  
  sco_names <- scan(paste0(orthofinder_results,
                           "/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"),
                    what = character())
  
  # initialize column names for ortholog lengths from each species
  sp_1 = paste0(colnames(orthogroups_tsv)[2], "_len")
  sp_2 = paste0(colnames(orthogroups_tsv)[3], "_len")
  
  cat(paste0(length(sco_names), " SCOs for ", sp_1, " and ", sp_2), sep = '\n')
  
  # filter orthogroups tsv to only single copy orthogroups, and add new columns
  # for protein lengths of ortholog from each species
  orthogroups_tsv <- orthogroups_tsv %>%
    filter(Orthogroup %in% sco_names) %>%
    mutate(ortholog_lengths = get_ortho_lengths(Orthogroup,
                                                orthofinder_results)) %>%
    separate(ortholog_lengths, into = c(sp_1, sp_2), sep = ',')
  
  return(orthogroups_tsv)
}


main()