library(tidyverse)

ortho_lengths <- read_csv('~/Desktop/smsk_microsporidia_domains/results/sco_summary_stats/SCO_protein_lengths_summary.csv')
ptc_tolerance <- readxl::read_xls('~/Desktop/smsk_microsporidia_domains/data/PTC_data/supp_12.xls')

yeast_names <- readLines('~/Desktop/smsk_microsporidia_domains/data/PTC_data/yeast.txt')

yeast_names <- as.list(yeast_names)
yeast_names <- lapply(yeast_names, function(x) {str_split(substr(x[1], 76, str_length(x[1])), "\\s+")[[1]][c(1, 2)]})

yeast_hash <- new.env()

for (pair in yeast_names) {
  yeast_hash[[pair[1]]] <- pair[2]
}

################################################################################

ortho_lengths <- ortho_lengths %>%
  select(species, yeast_ortho, microsp_len, yeast_len, microsp_yeast_ratio)

ptc_tolerance <- ptc_tolerance %>%
  rename(gene = GENEID, ptc_tolerated = `HMM count of PTCs tolerated per gene`) %>%
  select(gene, ptc_tolerated) %>%
  rowwise() %>%
  mutate(gene = yeast_hash[[gene]]) %>%
  # 777 of 2163 orthologs shared between microsporidia/rozella and yeast are
  # essential AND have PTC data available
  filter(gene %in% unique(ortho_lengths$yeast_ortho))

################################################################################

ptc_tolerance_hash <- new.env()
for (i in 1 : nrow(ptc_tolerance)) {
  ptc_tolerance_hash[[ptc_tolerance$gene[i]]] <- ptc_tolerance$ptc_tolerated[i]
}

ortho_lengths_ptc <- ortho_lengths %>%
  rowwise() %>%
  mutate(yeast_ptc_tolerance = ifelse(
    yeast_ortho %in% ptc_tolerance$gene,
    ptc_tolerance_hash[[yeast_ortho]],
    NA
  )) %>%
  # only keep ortholog pairs w/ PTC data available for the yeast ortholog
  filter(!is.na(yeast_ptc_tolerance))

################################################################################

# expectation: as PTC tolerance increases, average length ratio between ortholog
# pairs should decrease to reflect more PTC tolerant genes getting reduced
microsp_ptc_data <- filter(ortho_lengths_ptc, species != 'R_allo')
microsp_cor <- cor.test(x = microsp_ptc_data$yeast_ptc_tolerance,
                        y = microsp_ptc_data$microsp_yeast_ratio,
                        method = 'spearman')

microsp_ptc_plot <- ggplot(data = microsp_ptc_data, aes(x = yeast_ptc_tolerance,
                                                        y = microsp_yeast_ratio)) +
  geom_point(alpha = 0.2) +
  labs(title = str_c('All microsporidia\n', nrow(microsp_ptc_data), ' ortholog pairs\n',
                     'spearman rho = ', round(microsp_cor$estimate, 2), ', p = ',
                     round(microsp_cor$p.value))) +
  theme_bw()

ggsave(plot = microsp_ptc_plot, filename = 'All species.png',
       path = '~/Desktop/smsk_microsporidia_domains/results/Fig_1/Fig_1F')

################################################################################

# Make same plot but for each individual species
for (sp in unique(ortho_lengths_ptc$species)) {
  sp_ptc_data <- filter(ortho_lengths_ptc, species == sp)
  sp_cor <- cor.test(x = sp_ptc_data$yeast_ptc_tolerance,
                     y = sp_ptc_data$microsp_yeast_ratio,
                     method = 'spearman')
  
  sp_ptc_plot <- ggplot(data = sp_ptc_data, aes(x = yeast_ptc_tolerance,
                                                y = microsp_yeast_ratio)) +
    geom_point(alpha = 0.2) +
    labs(title = str_c(sp, '\n', nrow(sp_ptc_data), ' ortholog pairs\n',
                       'spearman rho = ', round(sp_cor$estimate, 2), ', p = ',
                       round(sp_cor$p.value))) +
    theme_bw()
  
  ggsave(plot = sp_ptc_plot, filename = str_c(sp, '.png'),
         path = '~/Desktop/smsk_microsporidia_domains/results/Fig_1/Fig_1F')
}