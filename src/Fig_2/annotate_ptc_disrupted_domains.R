# -----------------------------------------------------------------------------
#
# Find dispensible C-terminal yeast domains from PTC data
#
# Jason Jiang - Created: 2022/03/07
#               Last edited: 2022/03/09
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Using supp. table 11 from X et al's 2018 paper, infer dispensible
#       C-terminal domains in yeast proteins using premature stop codon effects
#       on yeast viability.
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))

################################################################################

# Column names for cath-resolve-hits output files
CRH_HEADER = c('query_id', 'match_id', 'score', 'boundaries',
               'resolved', 'cond_evalue', 'indp_evalue')

main <- function() {
  # ----------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = filepath to supp. table 11 from X et al's paper (PTC effects on yeast
  #        viability)
  #
  #   $2 = filepath to hashtable mapping SGD gene names to Uniprot names
  #        (from map_sgd_names_to_uniprot.R)
  #
  #   $3 = filepath to directory with domain architectures for each
  #        microsporidia-yeast ortholog pair
  #
  #   $4 = filepath to csv of yeast gene intron boundaries
  #        (from parse_gene_introns.R)
  #
  #   $5 = filepath to csv with all microsporidia - yeast ortholog pairs
  #
  #   $6 = filepath to save output files from this script into
  # ----------------------------------------------------------------------------
  # Parse in command line arguments
  args <- commandArgs(trailingOnly = T)
  
  ptc <- readxl::read_xls(args[1])
  sgd_to_uniprot <- readRDS(args[2])
  domain_archs_dir <- args[3]
  yeast_introns <- read_csv(args[4], show_col_types = F)
  ortho_pairs <- read_csv(args[5], show_col_types = F)
  out <- args[6]
  
  # make hashtable mapping each yeast gene conserved in microsporidia to its
  # domain architecture
  yeast_domain_archs <- map_yeast_genes_to_archs(ortho_pairs, domain_archs_dir)
  
  # annotate supp fig 11 dataframe with what protein domain are disrupted by
  # each ptc
  ptc_domain_disruptions <- get_ptc_domain_disruptions(ptc, yeast_domain_archs,
                                                       sgd_to_uniprot)
  
  # write annotated ptc dataframe to specified out directory
  write_csv(ptc_domain_disruptions, out)
}

################################################################################

### Helper functions

map_yeast_genes_to_archs <- function(ortho_pairs, domain_archs_dir) {
  # ----------------------------------------------------------------------------
  # Returns a hashtable (environment) mapping yeast uniprot gene names to
  # dataframes of their domain architectures (from cath-resolve-hits)
  #
  # Arguments:
  #   ortho_pairs: dataframe w/ info on all microsporidia - yeast gene pairs,
  #                from compare_sco_lengths snakefile
  #
  #   domain_archs_dir: filepath to directory with all cath-resolve-hits
  #                     domain architectures for each microsporidia - yeast
  #                     ortholog pair
  # ----------------------------------------------------------------------------
  # Extract each unique ortholog conserved between yeast and microsporidia
  # Use the species, orthogroup and ortholog name information later to load
  # in domain architecture information for each yeast ortholog
  yeast_orthos <- ortho_pairs %>%
    select(species, orthogroup, yeast_ortho) %>%
    distinct(yeast_ortho, .keep_all = T)
  
  # use above dataframe to load in domain arch data for all yeast orthologs to
  # microsporidia
  yeast_domain_archs <- new.env()
  for (i in 1 : nrow(yeast_orthos)) {
    # filepath to cath-resolve-hits domain architecture file for this yeast
    # ortholog to microsporidia
    crh_file <- str_c(domain_archs_dir, str_c(yeast_orthos$species[i], '_S_cere'),
                      yeast_orthos$orthogroup[i], 'crh',
                      yeast_orthos$yeast_ortho[i], sep = '/')
    
    # add entry to yeast_domain_archs, with key as uniprot name of yeast
    # ortholog and value as dataframe of the cath-resolve-hits domain
    # architecture file
    yeast_domain_archs[[yeast_orthos$yeast_ortho[i]]] <-
      read_delim(file = crh_file, comment = '#', col_names = CRH_HEADER,
                 delim = ' ', show_col_types = FALSE) %>%
      # add a column for the end boundaries of each domain in the ortholog
      mutate(dom_end = get_domain_end(resolved))
  }
  
  return(yeast_domain_archs)
}


get_domain_end <- Vectorize(function(resolved_boundaries) {
  # ----------------------------------------------------------------------------
  # Returns position of the last residue of a domain in a protein
  # Vectorized to with dplyr mutate
  #
  # Arguments:
  #   resolved_boundaries: 'resolved' column entry from a cath-resolve-hits
  #                         domain architecture file
  # ----------------------------------------------------------------------------
  domain_segments <- str_split(resolved_boundaries, ',')[[1]]
  last_segment <- domain_segments[length(domain_segments)]
  end <- str_split(last_segment, '-')[[1]][2]
  
  return(end)
})


get_ptc_disrupted_domains <- function(gene_uniprot, CDS_length, dist_from_CDS_end,
                                      yeast_domain_archs) {
  # ----------------------------------------------------------------------------
  # Returns a comma-separated string of all domains disrupted by a PTC in a 
  # yeast gene, with domains going from N -> C terminal.
  #
  # Helper function for get_ptc_domain_disruptions
  # ----------------------------------------------------------------------------
  doms_df <- yeast_domain_archs[[gene_uniprot]]
  ptc_position <- CDS_length - dist_from_CDS_end + 1
  
  doms_df <- doms_df %>%
    mutate(disrupted_by_ptc = (ptc_position <= as.integer(dom_end))) %>%
    filter(disrupted_by_ptc)
  
  disrupted_doms <- str_c(doms_df$match_id, collapse = ', ')
  if (disrupted_doms == '') {
    return(NA)
  }
  
  return(disrupted_doms)
}


get_ptc_domain_disruptions <- function(ptc, yeast_domain_archs, sgd_to_uniprot) {
  # ----------------------------------------------------------------------------
  # Annotates supplementary figure 11 from X et al's paper with what protein
  # domains are disrupted by each PTC
  #
  # Arguments:
  #   ptc: dataframe of supp figure 11
  #`  yeast_domain_archs: hashtable of yeast genes and their domain architectures
  #   sgd_to_uniprot: hashtable mapping sgd yeast gene names to uniprot names
  # ----------------------------------------------------------------------------
  ptc <- ptc %>%
    rename(gene = Gene, survival = `HMM PTC classification`) %>%
    # add this column just to account for disrepancy in CDS_length found in
    # supp fig 11 and protein lengths reported on uniprot
    mutate(prt_length = CDS_length - 1) %>%  # CDS length = number of amino acids + stop codon
    rowwise() %>%
    mutate(gene_uniprot = sgd_to_uniprot[[gene]]) %>%
    select(gene, gene_uniprot, prt_length, CDS_length, dist_from_CDS_end, survival) %>%
    # remove rows w/ missing survival data and filter dataframe to only have
    # PTCs for yeast genes that are in microsporidia + have domain architectures
    # assigned
    filter(!is.na(survival), !is.null(yeast_domain_archs[[gene_uniprot]])) %>%
    mutate(domains_disrupted = get_ptc_disrupted_domains(gene_uniprot, CDS_length,
                                                         dist_from_CDS_end,
                                                         yeast_domain_archs),
           num_domains_disrupted = ifelse(!is.na(domains_disrupted),
                                          length(str_split(domains_disrupted, ', ')[[1]]),
                                          0))
  
  return(ptc)
}

################################################################################

main()