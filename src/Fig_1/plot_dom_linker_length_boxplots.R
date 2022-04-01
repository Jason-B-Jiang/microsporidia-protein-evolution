# -----------------------------------------------------------------------------
#
# Plot linker and domain length difference boxplots
#
# Jason Jiang - Created: 2022/01/30
#               Last edited: 2022/03/30
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Create boxplots comparing differences in linker and domain lengths in
#       essential and non-essential microsporidia - yeast ortholog pairs.
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(ggsignif))

################################################################################

main <- function() {
  # ----------------------------------------------------------------------------
  # $1 = filepath to csv of microsporidia and yeast ortholog linker lengths
  # $2 = filepath to csv of microsporidia and yeast ortholog domain lengths
  # $3 = filepath to txt list of essential yeast genes, in UniProt format
  # $4 = filepath to directory for saving output plots to
  # ----------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  linkers <- read_csv(args[1], show_col_types = F)
  domains <- read_csv(args[2], show_col_types = F)
  essential_genes <- readLines(args[3])
  out_dir <- args[4]

  plot_and_save_boxplots(linkers, domains, essential_genes, out_dir)
}

################################################################################

plot_and_save_boxplots <- function(linkers, domains, essential_genes, out_dir) {
  # ----------------------------------------------------------------------------
  # Creates boxplots comparing domain/linker length differences between essential
  # and non-essential ortholog pairs, and saves to directory of dir_name in
  # out_dir.
  #
  # Arguments:
  #   linkers: dataframe of linker lengths
  #
  #   domains: dataframe of domain lengths
  #
  #   essential_genes: vector of essential gene names, for annotating ortholog
  #                  essentiality
  # 
  #   out_dir: base directory for storing results
  # ----------------------------------------------------------------------------
  domains <- domains %>%
    rename(microsp_len = microsp_dom_len, yeast_len = yeast_dom_len) %>%
    select(species, microsp_ortho, yeast_ortho, microsp_len, yeast_len) %>%
    filter(yeast_len > 0) %>%  # to avoid division by zero, just in case
    mutate(length_ratio = microsp_len / yeast_len,
           essential = yeast_ortho %in% essential_genes,
           type = 'domain')
  
  linkers <- linkers %>%
    select(-microsp_len, -yeast_len) %>%
    rename(microsp_len = microsp_linker_len, yeast_len = yeast_linker_len) %>%
    select(species, microsp_ortho, yeast_ortho, microsp_len, yeast_len) %>%
    filter(yeast_len > 0) %>%
    mutate(length_ratio = microsp_len / yeast_len,
           essential = yeast_ortho %in% essential_genes,
           type = 'linker')
  
  # calculate p-values for differences in domain/linker lengths between the
  # microsporidia species and yeast, and the difference in domain/linker lengths
  # between essential and non-essential ortholog pairs
  p_domains <- wilcox.test(domains$microsp_len, domains$yeast_len,
                           paired = TRUE)$p.value # bonferroni corrected
  
  p_domains_essential <- wilcox.test(filter(domains, essential)$length_ratio,
                                     filter(domains, !essential)$length_ratio)$p.value
  
  p_linkers <- wilcox.test(linkers$microsp_len, linkers$yeast_len,
                           paired = TRUE)$p.value
  
  p_linkers_essential <- wilcox.test(filter(linkers, essential)$length_ratio,
                                     filter(linkers, !essential)$length_ratio)$p.value
  
  # make boxplot comparing relative domain and linker lengths between the
  # microsporidia species and yeast, as well as between essential and
  # non-essential ortholog pairs
  dom_linkers <- rbind(domains, linkers)  # combine dataframes for plotting
  
  boxplt <- ggplot(data = dom_linkers, aes(x = essential, y = log(length_ratio, base = 2))) +
    geom_boxplot() +
    geom_signif(comparisons = list(c('TRUE', 'FALSE'))) +
    facet_wrap(~type) +
    labs(x = 'Essential?', y = 'log2(microsporidia / yeast)') +
    theme(axis.title = element_text(color = 'black', size = 18),
          axis.text = element_text(color = 'black', size = 18)) +
    theme_bw()
  
  # helper function to round p values to n significant digits
  # taken from https://stackoverflow.com/questions/43050903/round-to-significant-digits-only-with-decimal-portion-of-number-in-r
  my_signif = function(x, digits) floor(x) + signif(x %% 1, digits)
  
  # dataframe holding overall N. parisii vs yeast domain and linker length diffs
  # p values
  f_labels <- data.frame(type = c('domain', 'linker'),
                         label = c(str_c('overall: p = ', my_signif(p_domains, 2)),
                                   str_c('overall: p = ', my_signif(p_linkers, 2))))
  
  # add overall p-values for microsporidia domains vs yeast domains lengths,
  # microsporidia linkers vs yeast linkers lengths to the boxplot
  boxplt <- boxplt +
    geom_text(x = 1.9, y = -10.7, aes(label = label), data = f_labels)
  
  ggsave(plot = boxplt, filename = str_c(out_dir, '/Fig_1C.tiff'), units = 'in',
         width = 5.05, height = 5.44, dpi = 600)
}

################################################################################

main()

domains <- read_csv('~/Desktop/smsk_microsporidia_domains/results/sco_summary_stats/SCO_domain_lengths_summary.csv')
linkers <- read_csv('~/Desktop/smsk_microsporidia_domains/results/sco_summary_stats/SCO_linker_lengths_summary.csv')
essential_genes <- readLines('~/Desktop/smsk_microsporidia_domains/data/essential_yeast_genes/essential.txt')