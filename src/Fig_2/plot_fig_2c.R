# -----------------------------------------------------------------------------
#
# Plot Figure 2C
#
# Jason Jiang - Created: 2022/03/08
#               Last edited: 2022/03/08
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Plot fig 2c, too lazy to properly document
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(ggsignif))

################################################################################

### Global variables
SP_ORDER <- c('A_alge', 'A_locu', 'C_dike', 'D_muel', 'D_roes', 'E_aedi',
              'E_bien', 'E_brev', 'E_canc', 'E_cuni', 'E_hell', 'E_hepa', 'E_inte', 'E_roma',
              'H_erio', 'H_magn', 'H_tvae', 'N_apis', 'N_ausu', 'N_bomb', 'N_cera', 'N_cide',
              'N_disp', 'N_ferr', 'N_gran', 'N_homo', 'N_iron', 'N_majo', 'N_mino', 'N_okwa',
              'N_pari', 'O_coll', 'P_epip', 'P_neur', 'P_phil', 'S_loph', 'T_cont', 'T_homi',
              'T_rati', 'V_corn', 'V_culi', 'A_sp.', 'M_incu', 'M_daph', 'P_sacc', 'R_allo')

################################################################################

main <- function() {
  # ----------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = filepath to csv of PTC data annotated w/ domain disruptions
  #        (from annotate_ptc_disrupted_domains.R)
  #
  #   $2 = filepath to csv of aligned ortholog domain architectures for all
  #        microsporidia - yeast ortholog pairs
  #
  #   $3 = filepath + filename of fig 2C plot to save
  # ----------------------------------------------------------------------------
  # Parse in command line arguments
  args <- commandArgs(trailingOnly = T)
  ptc_domain_disruptions <- read_csv(args[1], show_col_types = F)
  aligned_domain_archs <- read_csv(args[2], show_col_types = F)
  out <- args[3]
  
  # get hashmap mapping yeast genes to dispensible c-terminal domains
  # if gene doesn't have entry in hashmap, then it has no dispensible
  # c-terminal domains
  dispensable_doms <- get_dispensable_doms(ptc_domain_disruptions)
  
  # add columns to aligned_domain_archs indicating lost/swapped domains
  # in the microsporidia ortholog that are dispensable c-terminal domains
  # in yeast
  dispensable_c_term_losses <- annotate_dispensable_c_term_changes(aligned_domain_archs,
                                                                   dispensable_doms[[1]],
                                                                   dispensable_doms[[2]])
  
  plot_proportion_dispensable_losses(dispensable_c_term_losses)
  plot_proportion_dispensable_losses_per_species(dispensable_c_term_losses, out)
  
  write_csv(dispensable_c_term_losses, stc_c(dirname(out),
                                             'dispensable_c_term_losses', sep = '/'))
  
  cat('Figure 2C plotted in results/Fig_2/Fig_2C', sep = '\n')
}

################################################################################

### Helper functions

get_dispensable_doms <- function(ptc_domain_disruptions) {
  # ----------------------------------------------------------------------------
  # Returns list -> [[1]] = dataframe, [[2]] = hashtable
  # ----------------------------------------------------------------------------
  dispensable_domains <- ptc_domain_disruptions %>%
    filter(survival == 'A', !is.na(domains_disrupted)) %>%
    group_by(gene_uniprot) %>%
    filter(num_domains_disrupted == max(num_domains_disrupted)) %>%
    distinct(gene_uniprot, num_domains_disrupted, .keep_all = T) %>%
    select(gene, gene_uniprot, survival, domains_disrupted, num_domains_disrupted)
  
  # turn above dataframe into a hashtable, with uniprot gene names as keys and
  # a vector of dispensable c-terminal domains as values
  dispensable_domains_hash <- new.env()
  for (i in 1 : nrow(dispensable_domains)) {
    dispensable_domains_hash[[dispensable_domains$gene_uniprot[i]]] <-
      str_split(dispensable_domains$domains_disrupted[i], ', ')[[1]]
  }
  
  return(list(dispensable_domains, dispensable_domains_hash))
}


annotate_dispensable_c_term_changes <- function(aligned_domain_archs,
                                                dispensable_doms,
                                                dispensable_doms_hash) {
  # ----------------------------------------------------------------------------
  # Write docstring later
  # ----------------------------------------------------------------------------
  losses_c_term <- aligned_domain_archs %>%
    filter(!is.na(lost_doms), yeast %in% dispensable_doms$gene_uniprot)
  
  # Define helper function for counting number of lost domains in each
  # microsporidia ortholog that is a dispensable c-terminal domain
  get_num_dispensable_losses <- function(gene, lost_doms, dispensable_doms_hash) {
    lost_doms <- str_split(lost_doms, ', ')[[1]]
    return(sum(lost_doms %in% dispensable_doms_hash[[gene]]))
  }
  
  losses_c_term <- losses_c_term %>%
    rowwise() %>%
    mutate(num_dispensable_losses = get_num_dispensable_losses(yeast,
                                                               lost_doms,
                                                               dispensable_doms_hash),
           total_dom_losses = length(str_split(lost_doms, ', ')[[1]]))
  
  return(losses_c_term)
}


plot_proportion_dispensable_losses <- function(dispensable_c_term_losses, out) {
  # ----------------------------------------------------------------------------
  # Write docstring later
  # ----------------------------------------------------------------------------
  total_losses = sum(dispensable_c_term_losses$total_dom_losses)
  total_dispensable_losses = sum(dispensable_c_term_losses$num_dispensable_losses)
  
  p_val <- prop.test(x = total_dispensable_losses, n = total_losses,
                     alternative = 'two.sided')$p.value
  
  losses_count_table <- matrix(c(total_dispensable_losses / total_losses,
                                 (total_losses - total_dispensable_losses) / total_losses),
                               ncol = 2)
  colnames(losses_count_table) <- c('Dispensable C-terminal domain', 'Other')
  losses_count_table <- as.table(losses_count_table)
  
  svg(filename = out, width = 5.7, height = 5.4)
  
  barplot(losses_count_table,
          ylab = 'proportion',
          border = 'black',
          col = 'black')
  
  dev.off()
}


plot_proportion_dispensable_losses_per_species <-
  function(dispensable_c_term_losses, out) {
    # ----------------------------------------------------------------------------
    # Write docstring later
    # ----------------------------------------------------------------------------
    domain_loss_summary <- dispensable_c_term_losses %>%
      group_by(species) %>%
      mutate(species_num_orthos = n(),
             species_dom_losses = sum(total_dom_losses),
             species_dom_losses_dispensable = sum(num_dispensable_losses),
             percent_dispensable_losses = species_dom_losses_dispensable / species_dom_losses) %>%
      select(species, species_num_orthos, species_dom_losses,
             species_dom_losses_dispensable, percent_dispensable_losses) %>%
      distinct()
    
    # order rows by phylogenetic arrangement of species
    domain_loss_summary <- left_join(data.frame(species = SP_ORDER),
                                     domain_loss_summary,
                                     by = 'species') %>%
      mutate(order = row_number())  # for ordering rows in plotting
    
    species_loss_plot <- ggplot(domain_loss_summary, aes(x = reorder(species, -order),
                                                         y = percent_dispensable_losses,
                                                         label = species_num_orthos)) +
      geom_segment(aes(x = reorder(species, -order), xend = species, y = 0,
                       yend = percent_dispensable_losses), color = "skyblue") +
      geom_point(color = "blue", size = 3, alpha = 0.6) +
      theme_light() +
      coord_flip() +
      geom_text(size = 3.5, hjust = -1) +
      labs(y = '% lost domains that are C-terminal dispensable') +
      theme(
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14, color = 'black')
      )
    
    ggsave(species_loss_plot, units = 'in', dpi = 300, width = 6.4, height = 7.2,
           filename = 'Fig_2C_B.tiff', path = out)
  }

################################################################################

main()


