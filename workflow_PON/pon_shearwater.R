#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("deepSNV")
library("abind")
library("tidyverse")

piles_to_counts <- function(file, regions, extended = F) {
  #INPUT:
  #files is a list of pileup files generated and named as by pysamstat.
  #regions is a data.frame with three first columns giving chrom, start, end
  #OUTPUT:
  # Array with dim pts x positions x bases where bases are A,T,C,G,D,a,t,c,g,d or
  #
  # Requires pileup_reader, subset_ranges, dplyr

  #require(dplyr)

  #Make skeleton of positions using the regions argument.
  owarn <- getOption("warn")
  options(warn=-1)
  skeleton <- apply(regions, 1, function(x){
      data.frame('chrom' = x[1], 'pos' = seq.int(as.numeric(x[2]), as.numeric(x[3])), stringsAsFactors = F)
    }) %>%
    bind_rows() %>%
    unique()

  options(warn=owarn)
  pileup_reader <-function(x) {
    data <-read.table(x,header = T, stringsAsFactors = F)
    return (data)
    }

  nucleotide_freq <-
    function(x) {
      if(extended) {
        tmp <-
          pileup_reader(x) %>%
          dplyr::select(chrom, pos, ref,
                        A = A_pp_fwd, `T` = T_pp_fwd, C = C_pp_fwd, G = G_pp_fwd, D = deletions_pp_fwd, I = insertions_pp_fwd, N = N_pp_fwd,
                        a = A_pp_rev, t = T_pp_rev, c = C_pp_rev, g = G_pp_rev, d = deletions_pp_rev,  i = insertions_pp_rev, n = N_pp_rev)
        
      } else {
        tmp <-
          pileup_reader(x) %>%
          dplyr::select(chrom, pos, ref,
                        A = A_pp_fwd, `T` = T_pp_fwd, C = C_pp_fwd, G = G_pp_fwd, D = deletions_pp_fwd,
                        a = A_pp_rev, t = T_pp_rev, c = C_pp_rev, g = G_pp_rev, d = deletions_pp_rev)
      }
      
      # Note subset_ranges fails when nrow(tmp) == 0
      res <-
        tmp %>%
        dplyr::left_join(skeleton, ., by = c('chrom' = 'chrom', 'pos' = 'pos')) %>%
        dplyr::select(-c(1,2,3)) %>% #Just get the counts, delete columns of chrom,pos,ref
        as.matrix()

      #Region positions not sequenced (NA) is set to 0 for all bases
      res[ is.na(res) ] <- 0
      
      return(res)
    }
  res <-nucleotide_freq(file)
  # abind::abind(res, along = 0) #old abind dependent way
  #res <- aperm(array(unlist(freq), c(dim(freq[[1]]), length(freq) ) ), c(3,1,2) )
  res3d <-array(res, dim=c(1, dim(res)[1], dim(res)[2]))
  if(extended) {
    dimnames(res3d) <- list(file, NULL, c('A','T','C','G','D', 'I', 'N', 
                                        'a','t','c','g','d', 'i', 'n') )
    return(res3d)  
  }
  dimnames(res3d) <- list(file, NULL, c('A','T','C','G','-','a','t','c','g','_') )
  return(res3d)
}

pon_pileups <-args[1]
panel_bed <-read.table(args[2])
pon_lists <-lapply(pon_pileups, piles_to_counts, regions=panel_bed)
pon_counts <-abind(pon_lists, along=1)
saveRDS(pon_counts, file = args[3])