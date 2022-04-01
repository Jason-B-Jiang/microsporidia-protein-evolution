# -----------------------------------------------------------------------------
#
# Make microsporidia domain architecture heatmap
#
# Jason Jiang - Created: 2022/01/28
#               Last edited: 2022/01/28
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Create a heatmap of domain architecture conservation of all microsporidia
#       orthologs with yeast.
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(cluster))
suppressMessages(library(ggdendro))

################################################################################

main <- function() {
  # ----------------------------------------------------------------------------
  # $1 = filepath to csv of domain architecture alignments for all microsporidia
  #      species
  # $2 = filepath to hashtable mapping microsporidia species to their clades
  # $3 = filepath + name of where to save data for generating the heatmap
  # $4 = filepath + name of the final heatmap image
  # ----------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  aligned_DA <- read_csv(args[1], show_col_types = F)
  clades <- readRDS(args[2])
  heatmap_data_out <- args[3]
  heatmap_out <- args[4]
  
  heatmap_data <- get_ortholog_DA_conservation(aligned_DA, clades)
  heatmap <- create_ortholog_heatmap(heatmap_data)
  
  write_csv(heatmap_data, heatmap_data_out)
  # save the heatmap image somehow
}

################################################################################

# Helper functions

get_ortholog_DA_conservation <- function(aligned_DA, clades) {
  # ----------------------------------------------------------------------------
  # Formats a dataframe of aligned domain architectures for making a heatmap of
  # domain architecture conservation with different yeast orthologs.
  #
  # aligned_DA: dataframe of aligned domain architectures of all microsporidia
  #             species with their yeast orthologs
  #
  # clades: hashtable mapping microsporidia species to their clades
  # ----------------------------------------------------------------------------
  ortholog_DA_conservation <- aligned_DA %>%
    rowwise() %>%
    # add column for the evolutionary clade of each microsporidia species
    mutate(clade = clades[[species]]) %>%
    ungroup() %>%
    # group species alphabetically by species names, clades amd yeast orthologs
    mutate(DA_conservation = get_DA_conservation_class(lost_doms,
                                                       gained_doms,
                                                       swapped_doms)) %>%
    select(species, clade, yeast, DA_conservation)
  
  # add rows for missing orthologs in each species
  return(fill_missing_orthos(ortholog_DA_conservation))
}


get_DA_conservation_class <- function(lost, gained, swapped) {
  # perfectly conserved DA
  if (all(is.na(c(lost, gained, swapped)))) {
    return('Conserved DA')
  } else if (!is.na(lost)) {
    return('Lost domain(s)')
  } else {
    return('Gained/swapped domain(s)')
  }
}

get_DA_conservation_class <- Vectorize(get_DA_conservation_class,
                                       vectorize.args = c('lost', 'gained',
                                                          'swapped'))


fill_missing_orthos <- function(ortholog_DA_conservation) {
  # ----------------------------------------------------------------------------
  # For each species, add rows of zero for missing orthologs
  # ----------------------------------------------------------------------------
  # get all unique orthologs that microsporidia share with yeast
  yeast_orthos <- unique(ortholog_DA_conservation$yeast)
  
  # separate lost_clans_list dataframe into separate dataframes for each species
  sp_DA_conservation <- ortholog_DA_conservation %>%
    group_by(species) %>%
    group_split()
  
  for (i in 1 : length(sp_DA_conservation)) {
    # get orthos that this microsporidia doesn't share with yeast
    missing_orthos <- setdiff(yeast_orthos, sp_DA_conservation[[i]]$yeast)
    
    # create dataframe of all missing domains for this species, with a count
    # of zero
    missing_orthos <- data.frame(
      species = unique(sp_DA_conservation[[i]]$species),
      clade = unique(sp_DA_conservation[[i]]$clade),
      yeast = missing_orthos,
      DA_conservation = 'Ortholog not conserved'
    )
    
    # add these missing domains to the original species dataframe
    sp_DA_conservation[[i]] <- rbind(sp_DA_conservation[[i]], missing_orthos)
  }
  
  # join all the split dataframes back together, now with rows added for missing
  # clans in each species
  combined <- bind_rows(sp_DA_conservation) %>%
    arrange(clade, species, yeast)
  
  return(combined)
}


get_agnes_diana_cluster <- function(ortholog_DA_conservation) {
  spread_DA <- ortholog_DA_conservation %>%
    select(-clade) %>%
    spread(species, DA_conservation)
  
  # convert all character columns to factors, for calculating gower distance
  spread_DA[sapply(spread_DA, is.character)] <-
    lapply(spread_DA[sapply(spread_DA, is.character)], as.factor)
  
  # gower dist
  gower_dist <- daisy(spread_DA, metric = 'gower')
  
  # diana and agnes clustering
  divisive_clust <- diana(as.matrix(gower_dist), diss = TRUE, keep.diss = TRUE)
  
  agnes_clust <- agnes(as.matrix(gower_dist), diss = TRUE, keep.diss = TRUE,
                       method = 'ward')
  
  # dendrogrames
  divisive_dendro <- as.dendrogram(divisive_clust)
  agnes_dendro <- as.dendrogram(agnes_clust)
  
  # get dendrogram orders
  divisive_order <- order.dendrogram(divisive_dendro)
  agnes_order <- order.dendrogram(agnes_dendro)
}


create_ortholog_heatmap <- function(ortholog_DA_conservation) {
  sp_order <- c('A_alge', 'A_locu', 'C_dike', 'D_muel', 'D_roes', 'E_aedi',
  'E_bien', 'E_brev', 'E_canc', 'E_cuni', 'E_hell', 'E_hepa', 'E_inte', 'E_roma',
  'H_erio', 'H_magn', 'H_tvae', 'N_apis', 'N_ausu', 'N_bomb', 'N_cera', 'N_cide',
  'N_disp', 'N_ferr', 'N_gran', 'N_homo', 'N_iron', 'N_majo', 'N_mino', 'N_okwa',
  'N_pari', 'O_coll', 'P_epip', 'P_neur', 'P_phil', 'S_loph', 'T_cont', 'T_homi',
  'T_rati', 'V_corn', 'V_culi', 'A_sp.', 'M_incu', 'M_daph', 'P_sacc', 'R_allo')
  
  agnes_diana_cluster <- get_agnes_diana_cluster(ortholog_DA_conservation)
  
  heatmap <- ggplot(ortholog_DA_conservation, aes(x = yeast, y = species)) +
    geom_tile(aes(fill = DA_conservation)) +
    scale_fill_manual(values=c("cadetblue", "yellow", "red", "black")) +
    scale_y_discrete(limits = rev(sp_order)) +
    labs(x = 'S. cerevisiae ortholog', title = 'Agglomerative clustering') +
    theme(axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  return(heatmap)
}


# Testing
set.seed(42)  # for consistency in random sampling
aligned_DA <- read_csv('~/Desktop/smsk_microsporidia_domains/results/domain_arch_conservation/merged_DA_alignments.csv')
clades <- readRDS('~/Desktop/smsk_microsporidia_domains/data/microsp_clades/clades.rds')
ortholog_DA_conservation <- get_ortholog_DA_conservation(aligned_DA, clades)

spread_DA <- ortholog_DA_conservation %>%
  select(-clade) %>%
  spread(species, DA_conservation)

# convert all character columns to factors, for calculating gower distance
spread_DA[sapply(spread_DA, is.character)] <-
  lapply(spread_DA[sapply(spread_DA, is.character)], as.factor)

gower_dist <- daisy(spread_DA, metric = 'gower')

# divisive clustering -> start with one cluster and gradually separate
# good for finding larger clusters
divisive_clust <- diana(as.matrix(gower_dist), diss = TRUE, keep.diss = TRUE)
plot(divisive_clust, main = "Divisive")

# agglomerative clustering -> start with small clusters, and gradually merge
# into larger clusters
aggl_clust <- hclust(gower_dist, method = "complete")
plot(aggl_clust, main = "Agglomerative, complete linkages")