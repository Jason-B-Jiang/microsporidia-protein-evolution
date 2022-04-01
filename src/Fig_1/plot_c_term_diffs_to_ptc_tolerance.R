# -----------------------------------------------------------------------------
#
# Get microsporidia-yeast ortholog c-terminal differences
#
# Jason Jiang - Created: 2022/02/28
#               Last edited: 2022/03/01
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Align microsporidia and yeast ortholog pairs with Needleman-Wunsch
#       algorithm, and determine extent of c-terminal reduction in microsporidia
#       orthologs to yeast
#
# Thanks to Brandon Murareanu for this RScript layout
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(seqinr))
suppressMessages(library(glue))

################################################################################

### Global constants

# Parameters for Needleman-Wunsch alignment
GAP_OPEN = 10.0
GAP_EXTEND = 0.5
MATRIX = 'EBLOSUM45'
OUT_FORMAT = 'fasta'

# Regex pattern for microsporidia species names
MICROSP_REGEX <- '[:upper:]_[:lower:]{2,4}\\.?'

# Regex pattern for yeast uniprot gene names
UNIPROT_REGEX <-
  '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'

################################################################################

main <- function() {
  # ----------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = directory with orthofinder results for each microsporidia species
  #        with yeast
  #
  #   $2 = filepath to supp table 12 from X et al's 2018 paper
  #
  #   $3 = filepath to 'yeast.txt' from Uniprot, which maps Saccharomyces
  #        Genome Database gene names to Uniprot gene names
  #
  #   $4 = filepath to target directory for writing output files into
  # ----------------------------------------------------------------------------
  # Parse command line arguments
  args <- commandArgs(trailingOnly = T)
  orthofinder_dir <- args[1]
  ptc_tolerance <- readxl::read_xls(args[2])
  out <- args[4]
  
  # Create hashtable mapping SGD yeast gene names to their Uniprot gene names,
  # using yeast.txt
  yeast_names <- as.list(readLines(args[3])) %>%
    lapply(function(x) {str_split(substr(x[1], 76, str_length(x[1])), "\\s+")[[1]][c(1, 2)]})
  
  yeast_hash <- new.env()
  
  for (pair in yeast_names) {
    yeast_hash[[pair[1]]] <- pair[2]
  }
  
  # Construct dataframe of extent of c-terminus truncation for each microsporidia
  # ortholog in each single-copy orthogroup, for all microsporidia species
  c_term_diffs <- get_c_term_diffs(orthofinder_dir, out)
  
  # Plot graphs comparing number of PTCs tolerated to extent of microsporidia
  # ortholog C-terminus for all ortholog pairs, in each microsporidia species.
  plot_c_term_vs_ptc_tolerance(c_term_diffs, ptc_tolerance, yeast_hash, out)
  
  # Write c_term_diffs dataframe as csv in specified out folder
  write_csv(c_term_diffs, str_c(out, '/c_term_diffs.csv'))
}

################################################################################

### Helper functions

get_c_term_diffs <- function(orthofinder_dir, out) {
  # ----------------------------------------------------------------------------
  # For each single-copy orthogroup for each microsporidia species and yeast,
  # align the ortholog pairs together and calculate extent of c-terminus
  # truncation in the microsporidia ortholog.
  #
  # Return a dataframe with microsporidia ortholog c-terminus truncation info
  # for all single-copy orthogroups, for all microsporidia species.
  # ----------------------------------------------------------------------------
  # Get all orthofinder results folders for each microsporidia species with
  # yeast
  species_results <- Sys.glob(str_c(orthofinder_dir, '/*'))
  sp_names <- str_extract(species_results, pattern = MICROSP_REGEX)
  
  # dataframe to hold c-terminal truncation information for each microsporidia
  # yeast ortholog pair
  c_term_diffs <- data.frame(species = character(),
                             yeast_gene = character(),
                             microsp_aligned = character(),
                             yeast_aligned = character(),
                             microsp_len = integer(),
                             yeast_len = integer(),
                             c_term_diff = integer())
  
  for (i in 1 : length(species_results)) {
    # Get filepaths to single-copy orthogroups + ortholog sequences for each
    # microsporidia species and yeast
    orthogroup_seqs <-
      Sys.glob(str_c(species_results[i],
                     '/Ortho*/Res*/Work*/Ortho*/Res*/SCO_seqs_split/OG*'))
    
    for (og in orthogroup_seqs) {
      # Load in microsporidia and yeast ortholog sequences in each orthogroup
      ortholog_seqs <- get_ortholog_seqs(og)
      
      cat(str_c('Aligning ', basename(og), ' for ', sp_names[i]), sep = '\n')
      
      # Align the ortholog sequences using NW algorithm, using specified
      # alignment parameters at top of this script
      aligned_orthos <- read.fasta(
        NW_align_orthologs(ortholog_seqs, out, sp_names[i]),
        seqtype = 'AA',
        as.string = T)
      
      # add row for this ortholog pair to c_term_diffs dataframe
      c_term_diffs <- add_to_c_term_diffs(aligned_orthos,
                                          c_term_diffs,
                                          sp_names[i])
    }
  }
  
  # Normalize each c-terminus difference by the length of the yeast ortholog
  c_term_diffs <- c_term_diffs %>%
    mutate(c_term_diff_normalized = c_term_diff / yeast_len)
  
  return(c_term_diffs)
}


get_ortholog_seqs <- function(og) {
  # ----------------------------------------------------------------------------
  # og: folder containing fasta files of ortholog sequences for a single copy
  #     orthogroup
  # ----------------------------------------------------------------------------
  seqs <- as.list(Sys.glob(str_c(og, '/*.fa')))
  
  # !any(str_detect(x, c('KKO', 'LWDP', 'ELQ') condition added to prevent any
  # microsporidia sequences that match the uniprot regex pattern from being
  # classified as yeast orthologs
  names(seqs) <- lapply(seqs,
                        function(x)
                        {ifelse(str_detect(x, pattern = UNIPROT_REGEX) & !any(str_detect(x, c('KKO', 'LWDP', 'ELQ'))),
                                  'yeast',
                                'microsp')}
                        )
  
  return(seqs)
}


NW_align_orthologs <- function(ortholog_seqs, out, sp_name) {
  # ----------------------------------------------------------------------------
  # Arguments:
  #   ortholog_seqs: named list of microsporidia and yeast ortholog pair
  #
  #   out: directory to save aligned orthologs into
  #
  #   sp_name: name of microsporidia species in the ortholog pair
  # ----------------------------------------------------------------------------
  dir.create(str_c(out, 'ortho_alignments', sp_name, sep = '/'), showWarnings = F, recursive = T)
  outfile = str_c(out, 'ortho_alignments', sp_name, basename(ortholog_seqs$yeast), sep = '/')
  
  # generate shell command to perform NW alignment on these ortholog sequences
  cmd = glue(str_c('needle -asequence {ortholog_seqs$microsp} ',
                   '-bsequence {ortholog_seqs$yeast} ',
                   '-gapopen {GAP_OPEN} -gapextend {GAP_EXTEND} ',
                   '-aformat3 {OUT_FORMAT} ',
                   '-datafile {MATRIX} -outfile {outfile}'))
  
  system(cmd)  # run the shell command for doing the needleman wunsch alignment
  
  # return filepath to the alignment file generated
  return(outfile)
}


add_to_c_term_diffs <- function(aligned_orthos, c_term_diffs, sp_name) {
  # ----------------------------------------------------------------------------
  # Pass
  # ----------------------------------------------------------------------------
  microsp_aln <- toString(aligned_orthos[[1]])
  yeast_aln <- toString(aligned_orthos[[2]])
  c_term_diff <- get_c_terminal_gap(microsp_aln, yeast_aln)
  
  return(c_term_diffs %>%
           add_row(data.frame(species = sp_name,
                              yeast_gene = getName(aligned_orthos[[2]]),
                              microsp_aligned = microsp_aln,
                              yeast_aligned = yeast_aln,
                              microsp_len = nchar(str_remove_all(aligned_orthos[[1]], '-')),
                              yeast_len = nchar(str_remove_all(aligned_orthos[[2]], '-')),
                              c_term_diff = c_term_diff)))
}


get_c_terminal_gap <- function(microsp_seq_aligned, yeast_seq_aligned) {
  # ----------------------------------------------------------------------------
  # Returns the c-terminal length difference between an aligned pair of microsporidia
  # and yeast orthologs.
  # ----------------------------------------------------------------------------
  i = nchar(microsp_seq_aligned)  # aligned seqs are of same length anyways
  microsp_c_term_diff = 0
  yeast_c_term_diff = 0
  
  if (substr(microsp_seq_aligned, i, i) == '-') {
    # microsporidia ortholog has a c-terminal truncation
    microsp_c_term_diff = microsp_c_term_diff + 1
    i = i - 1
    
    # keep moving back up the sequence and get length of continuous c-terminal
    # gap
    while (i >= 1 & substr(microsp_seq_aligned, i, i) == '-') {
      microsp_c_term_diff = microsp_c_term_diff + 1
      i = i - 1
    }
    
  } else if (substr(yeast_seq_aligned, i, i) == '-') {
    # yeast ortholog has a c-terminal truncation
    yeast_c_term_diff = yeast_c_term_diff + 1
    i = i - 1
    
    while (i >= 1 & substr(yeast_seq_aligned, i, i) == '-') {
      yeast_c_term_diff = yeast_c_term_diff + 1
      i = i - 1
    }
  }
  
  return(microsp_c_term_diff - yeast_c_term_diff)
}


###

plot_c_term_vs_ptc_tolerance <- function(c_term_diffs, ptc_tolerance,
                                         yeast_hash, out) {
  # ----------------------------------------------------------------------------
  # Using c_term_diffs dataframe produced by get_c_term_diffs, plot c-terminus
  # difference of microsporidia orthologs against number of PTCs tolerated
  # by their yeast orthologs (using data from supp. table 12 in X et al's 2018
  # paper)
  # ----------------------------------------------------------------------------
  # Format supp table 12 from X et al's paper
  ptc_tolerance <- ptc_tolerance %>%
    rename(gene = GENEID, ptc_tolerated = `HMM count of PTCs tolerated per gene`) %>%
    select(gene, ptc_tolerated) %>%
    rowwise() %>%
    # filter PTC tolerance data to only yeast genes that are conserved in
    # microsporidia
    mutate(gene = yeast_hash[[gene]]) %>%
    filter(gene %in% unique(c_term_diffs$yeast_gene))
  
  # Create hashtable from supp table 12, mapping each yeast gene to the number
  # of PTCs they tolerate
  ptc_tolerance_hash <- new.env()
  for (i in 1 : nrow(ptc_tolerance)) {
    ptc_tolerance_hash[[ptc_tolerance$gene[i]]] <- ptc_tolerance$ptc_tolerated[i]
  }
  
  # For each row / ortholog pair in c_term_diffs, add a column for the number
  # of PTCs tolerated by the yeast ortholog
  c_term_diffs_ptc <- c_term_diffs %>%
    rowwise() %>%
    mutate(yeast_ptc_tolerance = ifelse(
      yeast_gene %in% ptc_tolerance$gene,
      ptc_tolerance_hash[[yeast_gene]],
      NA
    )) %>%
    # only keep ortholog pairs w/ PTC data available for the yeast ortholog
    filter(!is.na(yeast_ptc_tolerance))
  
  # Plot c-terminus diffs vs ptc tolerance for all microsporidia species
  # (except R. allo)
  microsp_ptc_data <- filter(c_term_diffs_ptc, species != 'R_allo')
  microsp_cor <- cor.test(x = microsp_ptc_data$yeast_ptc_tolerance,
                          y = microsp_ptc_data$c_term_diff_normalized,
                          method = 'spearman')
  
  microsp_ptc_plot <- ggplot(data = microsp_ptc_data, aes(x = yeast_ptc_tolerance,
                                                          y = c_term_diff_normalized)) +
    geom_point(alpha = 0.2) +
    labs(title = str_c('All microsporidia\n', nrow(microsp_ptc_data), ' ortholog pairs\n',
                       'spearman rho = ', round(microsp_cor$estimate, 2), ', p = ',
                       round(microsp_cor$p.value)),
         x = 'Number of PTCs tolerated by yeast ortholog',
         y = 'C-terminus difference / yeast ortholog length') +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 14, color = 'black')) +
    theme_bw()
  
  ggsave(plot = microsp_ptc_plot, filename = str_c(out, '/All_species.tiff'),
         units = 'in', width = 4.9, height = 4.9, dpi = 600)
  
  # Plot c-terminus diffs vs ptc tolerance for each individual microsporidia
  # species (including R. allo outgroup)
  for (sp in unique(c_term_diffs_ptc$species)) {
    sp_df <- filter(c_term_diffs_ptc, species == sp)
    
    sp_cor <- cor.test(x = sp_df$yeast_ptc_tolerance,
                       y = sp_df$c_term_diff,
                       method = 'spearman')
    
    sp_ptc_plot <- ggplot(data = sp_df, aes(x = yeast_ptc_tolerance,
                                            y = c_term_diff)) +
      geom_point(alpha = 0.2) +
      labs(title = str_c(sp, '\n', nrow(sp_df), ' ortholog pairs\n',
                         'spearman rho = ', round(sp_cor$estimate, 2), ', p = ',
                         round(sp_cor$p.value)),
           x = 'Number of PTCs tolerated by yeast ortholog',
           y = 'C-terminus difference / yeast ortholog length') +
      theme(axis.title = element_text(size = 18),
            axis.text = element_text(size = 14, color = 'black')) +
      theme_bw()
      theme_bw()
    
    ggsave(plot = sp_ptc_plot, filename = str_c(sp, '.tiff'), path = out,
           units = 'in', width = 4.9, height = 4.9, dpi = 600)
  }
}

################################################################################

main()