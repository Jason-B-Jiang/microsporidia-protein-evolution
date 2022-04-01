#!/usr/bin/Rscript

# -----------------------------------------------------------------------------
#
# Merge cath-resolve-hits domain architectures
#
# Jason Jiang - Created: 2022/01/14
#               Last edited: 2022/01/14
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Merge all cath-resolve-hits outputs for ortholog pairs in a species
#       into a single csv file
#
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse, quietly = T))


# Define global variables
UNIPROT_REGEX = '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'
CRH_HEADER = c('query_id', 'match_id', 'score', 'boundaries',
               'resolved', 'cond_evalue', 'indp_evalue')


main <- function() {
  # args taken in:
  # $1: directory with orthogroup domain archs for a species pair
  # $2: filepath to hashtable of pfam families to clans
  
  args <- commandArgs(trailingOnly = T)
  og_domain_archs <- args[1]
  clan_hash <- readRDS(args[2])  # load in hashtable for pfam families to clans
  
  # Get name of microsporidia species for this set of orthogroup domain archs
  microsp_name <- str_extract(og_domain_archs, "[:upper:]_[:lower:]{2,4}\\.?")
  
  cat(str_c('Merging cath-resolve-hits results for ', microsp_name, ' and S_cere.'),
      sep = '\n')
  
  merged_archs <- merge_archs(og_domain_archs, microsp_name, clan_hash)
  
  # Combine domain architectures for both species into one dataframe
  matched_archs <- combine_archs(merged_archs[[microsp_name]],
                                 merged_archs[['S_cere']])
  
  # Write csvs for cath-resolve-hits outputs for all orthologs, one csv for
  # the microsporidia species and another for the yeast orthologs
  write_csv(merged_archs[[microsp_name]],
            str_c(og_domain_archs, '/', microsp_name, '_crh.csv'))
  
  write_csv(merged_archs[['S_cere']],
            str_c(og_domain_archs, '/S_cere_crh.csv'))
  
  # Write csvs for domain architectures of each species in pair, as well
  # as the combined domain architecture dataframe
  write_csv(matched_archs,
            str_c(og_domain_archs, '/', microsp_name, '_S_cere_DA.csv'))
}


# Make sure to vectorize this function for multiple of_domain_archs
merge_archs <- function(og_domain_archs, microsp_name, clan_hash) {
  # Returns a list of vectors, with each vector holding the cath-resolve-hits
  # outputs for the orthologs of a single species
  #
  # Inputs:
  # og_domain_archs: filepath to orthogroup domain architectures for a species
  # microsp_name: name of microsporidia species in the orthogroups
  # clan_hash: hashtable mapping pfam domain families to their clans
  #
  # Return:
  # Named list of vectors, each vector being the ortholog cath-resolve-hits
  # domain architectures for yeast or the microporidia species
  
  og_folders = list.files(path = og_domain_archs, pattern = 'OG0*',
                          full.names = T)
  
  crh_files = list.files(path = str_c(og_folders, '/crh'), full.names = T)
  
  # Get file paths to microsporida and yeast domain architecture crh outputs
  yeast_crh <-
    crh_files[str_detect(crh_files, UNIPROT_REGEX) & !str_detect(crh_files, '\\.1')]
  
  microsp_crh <- crh_files[!(crh_files %in% yeast_crh)]
  
  microsp_crh_df <- get_crh_df(microsp_crh, clan_hash)
  yeast_crh_df <- get_crh_df(yeast_crh, clan_hash)
  
  return_list <- list(microsp_crh_df, yeast_crh_df)
  names(return_list) <- c(microsp_name, 'S_cere')
  
  return(return_list)
}


get_crh_df <- function(crh_files, clan_hash) {
  # Returns the dataframe corresponding to the all crh output file specified by
  # crh_files
  # Add columns for orthogroup, clans and domain lengths
  
  # Initialize empty dataframe for holding all domain architecture data
  merged_crh_df <- data.frame(matrix(nrow = 0, ncol = length(CRH_HEADER) + 2))
  colnames(merged_crh_df) <- c('orthogroup', CRH_HEADER, 'dom_length')
  
  # Add rows for each ortholog's domain architecture to crh_df
  for (crh_file in crh_files) {
    # added as.character and [1] because of some weird csv writing issues
    orthogroup = as.character(str_match(crh_file, 'OG[:digit:]{7}'))[1]
    
    crh_df <- read_delim(file = crh_file, comment = '#', col_names = CRH_HEADER,
                         delim = ' ', show_col_types = FALSE) %>%
      mutate(orthogroup = orthogroup, dom_length = get_domain_length(resolved)) %>%
      rowwise() %>%
      select(orthogroup, everything())
  
    merged_crh_df <- rbind(merged_crh_df, crh_df)
  }
  
  # Add column for clans of each domain
  clans <- unlist(mget(as.character(merged_crh_df$match_id), envir = clan_hash,
                       ifnotfound = NA))
  merged_crh_df$clan <- clans
  
  return(merged_crh_df)
}


get_domain_length <- function(resolved_boundaries) {
  domain_segments <- str_split(resolved_boundaries, ',')[[1]]
  total_len = 0
  
  for (segment in domain_segments) {
    boundaries = as.numeric(str_split(segment, '-')[[1]])
    total_len = total_len + (boundaries[2] - boundaries[1] + 1)
  }
  
  return(total_len)
}

get_domain_length <- Vectorize(get_domain_length)


combine_archs <- function(microsp_df, yeast_df) {
  # Combines two dataframes from get_crh_df into one dataframe
  
  microsp_df <- microsp_df %>%
    group_by(query_id) %>%
    # collapse domain architectures into strings for each orthogroup
    mutate(match_id = paste(match_id, collapse = ', '),
           clan = paste(clan, collapse = ', '),
           dom_length = str_c(dom_length, collapse = ', ')) %>%
    # keep only one row for each orthogroup
    distinct(orthogroup, .keep_all = TRUE) %>%
    select(orthogroup, query_id, match_id, clan, dom_length)
  
  yeast_df <- yeast_df %>%
    group_by(query_id) %>%
    # collapse domain architectures into strings for each orthogroup
    mutate(match_id = paste(match_id, collapse = ', '),
           clan = paste(clan, collapse = ', '),
           dom_length = str_c(dom_length, collapse = ', ')) %>%
    # keep only one row for each orthogroup
    distinct(orthogroup, .keep_all = TRUE) %>%
    select(orthogroup, query_id, match_id, clan, dom_length)
  
  # Filter dataframes to only orthogroups in common
  shared_og <- intersect(microsp_df$orthogroup, yeast_df$orthogroup)
  
  microsp_df <- filter(microsp_df, orthogroup %in% shared_og)
  yeast_df <- filter(yeast_df, orthogroup %in% shared_og)
  
  combined_df <- data.frame(orthogroup = microsp_df$orthogroup,
                            microsp = microsp_df$query_id,
                            yeast = yeast_df$query_id,
                            microsp_DA = microsp_df$match_id,
                            microsp_DA_clan = microsp_df$clan,
                            microsp_dom_lengths = microsp_df$dom_length,
                            yeast_DA = yeast_df$match_id,
                            yeast_DA_clan = yeast_df$clan,
                            yeast_dom_lengths = yeast_df$dom_length)
  
  return(combined_df)
}


main()