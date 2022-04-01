# -----------------------------------------------------------------------------
#
# Cluster microsporidia by ortholog domain architectures
#
# Jason Jiang - Created: 2022/03/20
#               Last edited: 2022/03/20
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Use k-medoids clustering to cluster microsporidia species by their
#       domain architecture conservation of orthologs to yeast
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(cluster))
suppressMessages(library(Rtsne))

################################################################################

main <- function() {
  # ----------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = filepath to csv of aligned microsporidia - yeast ortholog domain
  #        architectures, produced by assign_and_align_domain_archs Snakefile
  #
  #   $2 = filepath to rds file of hashtable mapping microsporidia species to
  #        their evolutionary clades
  #
  #   $3 = filepath to directory for saving outputs in
  # ----------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  aligned_DA <- read_csv(args[1], show_col_types = F)
  microsp_clades <- readRDS(args[2])
  out <- args[3]
  
  DA_clusters <- cluster_microsp_by_DA(aligned_DA)
  cluster_tsne_plot <- make_tsne_cluster_plot(DA_clusters[[1]],
                                              DA_clusters[[2]],
                                              microsp_clades)
  
  saveRDS(DA_clusters[[1]], str_c(out, '/clustered_aligned_DA.rds'))
  saveRDS(DA_clusters[[2]], str_c(out, '/gower_dissimilarity.rds'))
  ggsave(cluster_tsne_plot)  # fill out arguments later
}

################################################################################

### Helper functions

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
                                       vectorize.args = c('lost',
                                                          'gained',
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
      yeast = missing_orthos,
      DA_conservation = 'Ortholog not conserved'
    )
    
    # add these missing domains to the original species dataframe
    sp_DA_conservation[[i]] <- rbind(sp_DA_conservation[[i]], missing_orthos)
  }
  
  # join all the split dataframes back together, now with rows added for missing
  # clans in each species
  combined <- bind_rows(sp_DA_conservation) %>%
    arrange(species, yeast)
  
  return(combined)
}


get_ortholog_DA_conservation <- function(aligned_DA) {
  # ----------------------------------------------------------------------------
  # Formats a dataframe of aligned domain architectures for making a heatmap of
  # domain architecture conservation with different yeast orthologs.
  #
  # aligned_DA: dataframe of aligned domain architectures of all microsporidia
  #             species with their yeast orthologs
  # ----------------------------------------------------------------------------
  ortholog_DA_conservation <- aligned_DA %>%
    # group species alphabetically by species names, clades amd yeast orthologs
    mutate(DA_conservation = get_DA_conservation_class(lost_doms,
                                                       gained_doms,
                                                       swapped_doms)) %>%
    select(species, yeast, DA_conservation)
  
  ortholog_DA_conservation <- fill_missing_orthos(ortholog_DA_conservation)
  
  ortholog_DA_conservation$DA_conservation <-
    as.factor(ortholog_DA_conservation$DA_conservation)
  
  ortholog_DA_conservation <- ortholog_DA_conservation %>%
    arrange(species) %>%
    spread(yeast, DA_conservation, fill = 'Ortholog not conserved')
  
  row.names(ortholog_DA_conservation) <-
    ortholog_DA_conservation$species
  
  return(ortholog_DA_conservation)
}


cluster_microsp_by_DA <- function(aligned_DA) {
  # Format aligned_DA for clustering + label by ortholog domain architecture
  # conservation
  ortholog_DA_conservation <- get_ortholog_DA_conservation(aligned_DA)
  
  # Get dissimilarity matrix between species and their ortholog domain
  # architectures using Gower's distance
  # Exclude first column (species names) from dissimilarity calculations
  gower_dist <- daisy(ortholog_DA_conservation[, -1], metric = 'gower')
  
  # Cluster colleges w/ k-medoids clustering (PAM)
  # First, choose ideal number of clusters to maximize silhouette width
  sil_width <- c(NA)
  
  for (i in 2 : 100) {
    
    pam_fit <- pam(gower_dist,
                   diss = TRUE,
                   k = i)
    
    sil_width[i] <- pam_fit$silinfo$avg.width
    
  }
  
  # Number of clusters maximizing silhouette width is the "optimal" number of
  # clusters
  max_sil_width <- which.max(sil_width)
  pam_fit <- pam(gower_dist, diss = TRUE, k = max_sil_width)
  
  return(pam_fit)
}


make_tsne_cluster_plot <- function(DA_clusters, gower_dist, microsp_clades) {
  tsne_obj <- Rtsne(gower_dist, is_distance = TRUE)
  
  tsne_data <- tsne_obj$Y %>%
    data.frame() %>%
    # x and y coordinates after projecting onto 2d plane
    setNames(c("X", "Y")) %>%
    mutate(cluster = factor(pam_fit$clustering),
           species = ortholog_DA_conservation$species)
  
  ggplot(aes(x = X, y = Y, label = species), data = tsne_data) +
    geom_point(aes(color = cluster)) +
    geom_text(angle = 45)
  
  
}
