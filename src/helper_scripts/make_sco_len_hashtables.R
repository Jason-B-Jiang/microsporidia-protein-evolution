# -----------------------------------------------------------------------------
#
# Create orthogroup length hashtables
#
# Jason Jiang - Created: 2022/01/19
#               Last edited: 2022/01/19
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: For each microsporidia-yeast pair passed in, create a hash table
#       for that pair mapping orthogroups to their microsporidia and yeast
#       ortholog lengths
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

library(stringr)
library(readr)

################################################################################

main <- function() {
  # $1 = output directory for hashtables
  # $2 - $n = all csvs of orthogroups to process
  args <- commandArgs(trailingOnly = T)
  output_dir <- args[1]
  species_pairs <- args[2 : length(args)]
  microsp_names = str_extract(species_pairs, '[:upper:]_[:lower:]{2,4}\\.?')
  
  write_sco_hashes(species_pairs, microsp_names, output_dir)
}

################################################################################

# Helper functions

write_sco_hashes <- function(species_pairs, microsp_names, output_dir) {
  # species_pairs and microsp_names should be character vectors of the same
  # length, and with each position corresponding to the other
  for (i in 1 : length(species_pairs)) {
    microsp_name = microsp_names[i]
    sco_csv = read_csv(species_pairs[i], show_col_types = FALSE)
    sco_hash = new.env()  # hashtable of orthogroup ortholog lengths
    
    # get vectors for microsp and yeast lengths + orthogroups
    # iterate over vectors in parallel and save in sco_hash
    # save sco_hash, and move onto the next species pair
    orthogroups = sco_csv$Orthogroup
    microsp_lens = as.numeric(sco_csv[[str_c(microsp_name, '_len')]])
    yeast_lens = as.numeric(sco_csv[['S_cere_len']])
    
    for (j in 1 : length(orthogroups)) {
      og = orthogroups[j]
      microsp_len = microsp_lens[j]
      yeast_len = yeast_lens[j]
      
      sco_hash[[og]] <- c(microsp_len, yeast_len)
    }
    
    # save the hashtable in specified output directory
    saveRDS(sco_hash, file = str_c(output_dir, '/', microsp_name, '_SCO_hash.rds'))
    cat(str_c('Created ', str_c(output_dir, '/', microsp_name, '_SCO_hash.rds')),
        sep = '\n')
  }
}

################################################################################

main()