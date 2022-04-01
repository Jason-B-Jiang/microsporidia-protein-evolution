library(tidyverse)
setwd('~/Desktop/smsk_microsporidia_domains/results/Fig_2/Fig_2E/PANTHER_enrichment')


go_enrichment <- read_tsv('merged.txt') %>%
  mutate(negative_log_fdr = -log(fdr, base = 2)) %>%
  group_by(class) %>%
  arrange(-negative_log_fdr, .by_group = T) %>%
  ungroup() %>%
  mutate(order = row_number()) %>%  # create ordering for graphing
  rename(`PANTHER GO Slim Category` = class)


ggplot(data = go_enrichment, aes(x = reorder(term, -order),
                                 y = negative_log_fdr,
                                 fill = `PANTHER GO Slim Category`)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  labs(y = '-log₂(FDR)') +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, color = 'black'),
        axis.title.x = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black'))