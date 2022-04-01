# -----------------------------------------------------------------------------
#
# Get orthogroup linker summary statistics
#
# Jason Jiang - Created: 2022/01/21
#               Last edited: 2022/01/24
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Compare linker lengths between microsporidia and yeast orthologs
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))

################################################################################

main <- function() {
  # $1: output directory
  # $2: path to directory with hashtables of ortholog lengths for each
  # microsporidia - yeast species pair
  # $3: path to csvs of domain architectures in orthologs amongst each
  # microsporidia - yeast species pair
  # $4 - $n: names of microsporidia species to analyze
  args <- commandArgs(trailingOnly = T)
  out_dir <- args[1]
  hash_dir <- args[2]
  DA_dir <- args[3]
  microsp_names = args[4 : length(args)]
  
  result_dfs <- compare_domains_and_linkers(DA_dir, hash_dir, microsp_names)
  
  write_csv(result_dfs[['linkers']], str_c(out_dir, '/SCO_linker_lengths_summary.csv'))
  write_csv(result_dfs[['domains']], str_c(out_dir, '/SCO_domain_lengths_summary.csv'))
}

################################################################################

# Helper functions for main

compare_domains_and_linkers <- function(DA_dir, hash_dir, microsp_names) {
  # Initialize empty domain length dataframe, and empty linker length df
  linker_df <- data.frame(matrix(ncol = 10, nrow = 0))
  colnames(linker_df) <-
    c('species', 'microsp_ortho', 'yeast_ortho', 'microsp_len', 'yeast_len',
      'microsp_linker_len', 'yeast_linker_len', 'p_species',
      'p_species_corrected', 'p_overall')
  
  # Get linker df for each species, and accumulate to overall linker df
  # Get domain df for each species, and accumulate to overall domain df
  domain_df <- data.frame(matrix(ncol = 13, nrow = 0))
  colnames(domain_df) <-
    c('species', 'orthogroup', 'microsp_ortho', 'yeast_ortho', 'microsp_dom',
      'yeast_dom', 'microsp_clan', 'yeast_clan', 'microsp_dom_len',
      'yeast_domain_len', 'p_species', 'p_species_corrected', 'p_overall')
  
  for (i in 1 : length(microsp_names)) {
    DA_df <- read_csv(str_c(DA_dir, '/', microsp_names[i], '_S_cere/',
                            microsp_names[i], '_S_cere_DA.csv'),
                      show_col_types = FALSE)
    
    ortho_len_hash <-
      readRDS(str_c(hash_dir, '/', microsp_names[i], '_SCO_hash.rds'))
    
    sp_linkers <- analyze_sp_linkers(DA_df, ortho_len_hash, microsp_names[i])
    sp_domains <- analyze_sp_domains(DA_df, microsp_names[i])
    
    linker_df <- rbind(linker_df, sp_linkers)
    domain_df <- rbind(domain_df, sp_domains)
  }
  
  # add columns for corrected per-species p-values (using bonferroni) and
  # overall p-values for all species
  linker_df$p_species_corrected <- linker_df$p_species * length(microsp_names)
  linker_df$p_overall <- wilcox.test(
    linker_df$microsp_linker_len,
    linker_df$yeast_linker_len,
    paired = TRUE,
    alternative = 'two.sided'
  )$p.value
  
  domain_df$p_species_corrected <- domain_df$p_species * length(microsp_names)
  domain_df$p_overall <- wilcox.test(
    domain_df$microsp_dom_len,
    domain_df$yeast_dom_len,
    paired = TRUE,
    alternative = 'two.sided'
  )$p.value
  
  
  # return hashmap of linker and domain dataframes
  return_hash <- new.env()
  return_hash[['linkers']] <- linker_df
  return_hash[['domains']] <- domain_df
  
  return(return_hash)
}


analyze_sp_linkers <- function(DA_df, ortho_len_hash, sp) {
  sp_linker_df <- DA_df %>%
    rowwise() %>%
    mutate(microsp_len = ortho_len_hash[[orthogroup]][1],
           yeast_len = ortho_len_hash[[orthogroup]][2]) %>%
    ungroup() %>%
    mutate(microsp_linker_len = get_linker_len(microsp_len, microsp_dom_lengths),
           yeast_linker_len = get_linker_len(yeast_len, yeast_dom_lengths),
           species = sp)
  
  sp_linker_df$p_species <- wilcox.test(sp_linker_df$microsp_linker_len,
                                        sp_linker_df$yeast_linker_len,
                                        paired = TRUE,
                                        alternative = 'two.sided')$p.value
  sp_linker_df$p_species_corrected <- NA
  sp_linker_df$p_overall <- NA
  
  return(format_sp_linker_df(sp_linker_df))
}

get_linker_len <- function(prt_len, dom_len) {
  dom_len <- sum(as.numeric(str_split(dom_len, ', ')[[1]]))
  return(as.numeric(prt_len) - dom_len)
}

get_linker_len <- Vectorize(get_linker_len,
                                vectorize.args = c('prt_len', 'dom_len'))

format_sp_linker_df <- function(sp_linker_df) {
  sp_linker_df <- sp_linker_df %>%
    select(species, microsp, yeast, microsp_len, yeast_len,
           microsp_linker_len, yeast_linker_len, p_species,
           p_species_corrected, p_overall) %>%
    rename(microsp_ortho = microsp, yeast_ortho = yeast)
  
  return(sp_linker_df)
}


analyze_sp_domains <- function(DA_df, sp) {
  # Create two separate dataframes for yeast and microsporidia domains
  microsp_df <- DA_df %>%
    select(-yeast, -yeast_DA, -yeast_dom_lengths, -yeast_DA_clan) %>%
    separate_rows(microsp_DA, microsp_dom_lengths, microsp_DA_clan, sep = ', ') %>%
    group_by(orthogroup, microsp_DA) %>%
    # merge together entries for repeat domains in the same protein
    mutate(microsp_merged_dom_length = sum(as.numeric(microsp_dom_lengths))) %>%
    distinct(orthogroup, microsp_DA, .keep_all = TRUE) %>%
    select(-microsp_dom_lengths) %>%
    ungroup() %>%
    mutate(id = str_c(orthogroup, microsp_DA, sep = ','))
    
  
  yeast_df <- DA_df %>%
    select(-microsp, -microsp_DA, -microsp_dom_lengths, -microsp_DA_clan) %>%
    separate_rows(yeast_DA, yeast_dom_lengths, yeast_DA_clan, sep = ', ') %>%
    group_by(orthogroup, yeast_DA) %>%
    # merge lengths of tandem repeat domains together
    mutate(yeast_merged_dom_length = sum(as.numeric(yeast_dom_lengths))) %>%
    distinct(orthogroup, yeast_DA, .keep_all = TRUE) %>%
    select(-yeast_dom_lengths) %>%
    ungroup() %>%
    mutate(id = str_c(orthogroup, yeast_DA, sep = ','))
  
  # use orthogroup + domain as identifier
  # filter both dfs to only identifiers in common
  # verify that all orthogroups and domains are in the same order
  # calculate p value for domain lengths
  # get median difference for this particular species
  # merge into one dataframe and return
  shared_doms <- intersect(microsp_df$id, yeast_df$id)
  
  microsp_df <- microsp_df %>%
    filter(id %in% shared_doms) %>%
    select(-id)
  
  yeast_df <- yeast_df %>%
    filter(id %in% shared_doms) %>%
    select(-id)
  
  # combine both dataframes into one
  microsp_yeast_df = data.frame(
    species = rep(sp, nrow(microsp_df)),
    orthogroup = microsp_df$orthogroup,
    microsp_ortho = microsp_df$microsp,
    yeast_ortho = yeast_df$yeast,
    microsp_dom = microsp_df$microsp_DA,
    yeast_dom = yeast_df$yeast_DA,
    microsp_clan = microsp_df$microsp_DA_clan,
    yeast_clan = yeast_df$yeast_DA_clan,
    microsp_dom_len = microsp_df$microsp_merged_dom_length,
    yeast_dom_len = yeast_df$yeast_merged_dom_length
  )
  
  # calculate p-value for differences in microsporidia-yeast domain lengths
  # for this microsporidia species
  microsp_yeast_df$p_species <- wilcox.test(microsp_yeast_df$microsp_dom_len,
                                            microsp_yeast_df$yeast_dom_len,
                                            paired = TRUE,
                                            alternative = 'two.sided')$p.value
  
  # initialize columns for corrected p values, and p value for domain length
  # differences amongst all microsporidia species
  microsp_yeast_df$p_species_corrected <- NA
  microsp_yeast_df$p_overall <- NA
  
  return(microsp_yeast_df)
}

################################################################################

main()