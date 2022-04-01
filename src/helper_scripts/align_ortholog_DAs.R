# -----------------------------------------------------------------------------
#
# Align and compare ortholog domain architectures
#
# Jason Jiang - Created: 2022/01/17
#               Last edited: 2022/01/17
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Align domain architectures of single copy ortholog pairs, and detect
#       domain architectural differences between orthologs
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

library(tidyverse)
library(NameNeedle)  # For Needleman-Wunsch string alignments
library(stringdist)  # For calculating number of differences between domain archs


main <- function() {
  # $1 = filepath to csv of domain architectures
  args <- commandArgs(trailingOnly = T)
  
  DA_df <- read_csv(args[1])
  
  aligned_DA_df <- align_and_compare_DAs(DA_df)
  
  write_csv(aligned_DA_df,
            file = str_c(str_remove(args[1], '\\.csv'), '_aligned.csv'))
}


align_and_compare_DAs <- function(DA_df) {
  # Adds new column to dataframe of domain architectures for aligned domain
  # architectures + domain architecture differences between ortholog pairs
  DA_df <- DA_df %>%
    mutate(microsp_clan_DA = get_clan_DA(microsp_DA, microsp_DA_clan),
           yeast_clan_DA = get_clan_DA(yeast_DA, yeast_DA_clan)) %>%
    rowwise() %>%
    # Switch to rowwise operations as functions using environments/hashmaps
    # don't support vectorized operations
    mutate(aligned_DA = align_DA(microsp_clan_DA, yeast_clan_DA),
           DA_score = get_DA_score(aligned_DA),
           DA_diffs = ifelse(DA_score < 1,
                             get_DA_diffs(aligned_DA, microsp_DA, yeast_DA),
                             NA)) %>%
    # Switch back to vectorized operations
    ungroup() %>%
    mutate(lost_doms = get_lost_doms(DA_diffs),
           gained_doms = get_gained_doms(DA_diffs),
           swapped_doms = get_swapped_doms(DA_diffs))

  return(DA_df)
}

################################################################################

# Helper functions

get_clan_DA <- function(DA, DA_clan) {
  DA <- str_split(DA, ', ')[[1]]
  DA_clan <- str_split(DA_clan, ', ')[[1]]
  DA_clan[DA_clan == 'NA'] <- NA  # convert NA strings to real NAs
  
  return(str_c(coalesce(DA_clan, DA), collapse = ', '))
}

get_clan_DA <- Vectorize(get_clan_DA, vectorize.args = c('DA', 'DA_clan'))


get_DA_letter_map <- function(DA_1, DA_2) {
  # Get all domains from a pair of domain architectures
  doms_1 <- str_split(DA_1, ', ')[[1]]
  doms_2 <- str_split(DA_2, ', ')[[1]]
  all_doms <- unique(c(doms_1, doms_2))
  
  # Create a mapping of domains to unique letter representations
  dom_to_letters = new.env()
  
  for (i in 1:length(all_doms)) {
    dom_to_letters[[all_doms[i]]] <- toupper(letters[i])
  }
  
  return(dom_to_letters)
}

  
align_DA <- function(DA_1, DA_2) {
  dom_to_letters <- get_DA_letter_map(DA_1, DA_2)
  
  DA_1_letters <- str_c(
    sapply(str_split(DA_1, ', ')[[1]],
                         function(x) {dom_to_letters[[x]]}),
    collapse = '')
  
  DA_2_letters <- str_c(
    sapply(str_split(DA_2, ', ')[[1]],
                         function(x) {dom_to_letters[[x]]}),
    collapse = '')
  
  # do NW alignment using letter representations for domain architectures
  needle_params = list(MATCH = 0, MISMATCH = -3, GAP = -10, GAPCHAR = '*')
  aligned = needles(DA_1_letters, DA_2_letters, params = needle_params)
  
  # Return '|' separated string of aligned domain architectures
  return(str_c(aligned$align1, aligned$align2, sep = ' | '))
}


get_DA_score <- function(DA_aligned) {
  DA = str_split(DA_aligned, ' \\| ')[[1]]
  
  # DA score = number of matching domains / total alignment length
  return((nchar(DA[1]) - stringdist(DA[1], DA[2])) / nchar(DA[1]))
}

get_DA_score <- Vectorize(get_DA_score)


get_DA_diffs <- function(aligned_DA, DA_1, DA_2) {
  # Identifies lost, gained and swapped domains between two aligned DAs
  
  # keep track of domain architectural changes
  DA_changes <- new.env()
  DA_changes[['lost']] <- c()
  DA_changes[['gained']] <- c()
  DA_changes[['swapped']] <- c()
  
  # turn aligned letter DAs into separate vectors of letters
  split_DA_letters = str_split(aligned_DA, ' \\| ')[[1]]
  DA_1_letters = str_split(split_DA_letters[1], '')[[1]]
  DA_2_letters = str_split(split_DA_letters[2], '')[[1]]
  
  # create mapping between letters and domains to map letters back to domains
  # can push this into a helper function to make things clearer
  doms_to_letters <- get_DA_letter_map(DA_1, DA_2)
  letters_to_doms <- new.env()
  
  for (dom in names(doms_to_letters)) {
    letters_to_doms[[doms_to_letters[[dom]]]] <- dom
  }
  
  # iterate over DA alignment to find DA changes
  for (i in 1 : length(DA_1_letters)) {
    if (DA_1_letters[i] != DA_2_letters[i]) {
      if (DA_1_letters[i] == '*') {
        # domain was lost in first DA relative to the second DA
        DA_changes[['lost']] <- append(DA_changes[['lost']],
                                       letters_to_doms[[DA_2_letters[i]]])
      } else if (DA_2_letters[i] == '*') {
        # domain was gained in first DA relative to the second
        DA_changes[['gained']] <- append(DA_changes[['gained']],
                                         letters_to_doms[[DA_1_letters[i]]])
      } else {
        # domain swap/change between the two DAs
        DA_changes[['swapped']] <-
          append(DA_changes[['swapped']], str_c(letters_to_doms[[DA_2_letters[i]]],
                                                ' -> ',
                                                letters_to_doms[[DA_1_letters[i]]]))
      }
    }
  }
  
  return(format_DA_changes_string(DA_changes))
}

format_DA_changes_string <- function(DA_changes) {
  return(str_c('Lost: ', str_c(DA_changes[['lost']], collapse = ', '), ' | ',
               'Gained: ', str_c(DA_changes[['gained']], collapse = ', '), ' | ',
               'Swapped: ', str_c(DA_changes[['swapped']], collapse = ', ')))
}


get_lost_doms <- function(DA_diffs) {
  if (is.na(DA_diffs)) {
    return(NA)
  }
  
  lost_doms = str_remove(str_split(DA_diffs, ' \\| ')[[1]][1], 'Lost: ')
  
  if (lost_doms == '') {
    return(NA)
  } else {
    return(lost_doms)
  }
}

get_gained_doms <- function(DA_diffs) {
  if (is.na(DA_diffs)) {
    return(NA)
  }
  
  gained_doms = str_remove(str_split(DA_diffs, ' \\| ')[[1]][2], 'Gained: ')
  
  if (gained_doms == '') {
    return(NA)
  } else {
    return(gained_doms)
  }
}

get_swapped_doms <- function(DA_diffs) {
  if (is.na(DA_diffs)) {
    return(NA)
  }
  
  swapped_doms = str_remove(str_split(DA_diffs, ' \\| ')[[1]][3], 'Swapped: ')
  
  if (swapped_doms == '') {
    return(NA)
  } else {
    return(swapped_doms)
  }
}

get_lost_doms <- Vectorize(get_lost_doms)
get_gained_doms <- Vectorize(get_gained_doms)
get_swapped_doms <- Vectorize(get_swapped_doms)

################################################################################

main()  # run the script