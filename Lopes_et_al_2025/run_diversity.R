rm(list = ls())

library(adegenet)
library(hierfstat)
library(dplyr)
# import read_genepop function that keeps population assignments
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/refs/heads/master/scripts/adegent_functions.R")


#############################
# calculate diversity metrics
#############################

species <- "Glyphorynchus spirurus"
project <- "glyphorynchus_"
postfix <- "discosnp_all_final"
fltr <- "_01_00"

# root directory where your data are
data_path <- "/media/tomas/Data2/ddRAD/glyphorynchus/discosnp_glyphorynchus/"
# directory where individual data are
res_path <- paste0(data_path, postfix, "/")

if (file.exists(paste0(res_path, project, postfix, fltr, '_genepop.gen'))) { # execute loop only if taxon exists
  taxon <- read_genepop(file = paste0(res_path, project, postfix, fltr, '_genepop.gen'))
  
  # convert to data frame
  tmp <- genind2df(taxon) %>%
    mutate(across(2:ncol(.), ~as.numeric(.)))
  
  stats <- basic.stats(tmp)
  ho <- c(mean(stats$Ho[,1]), mean(stats$Ho[,2]), mean(stats$Ho[,3]), 
          mean(stats$Ho[,4]), mean(stats$Ho[,5]), mean(stats$Ho[,6]))
  he <- c(mean(stats$Hs[,1]), mean(stats$Hs[,2]), mean(stats$Hs[,3]),
          mean(stats$Hs[,4]), mean(stats$Hs[,5]), mean(stats$Hs[,6]))
  
  richness <- allelic.richness(tmp)
  ri <- richness$Ar %>%  
    summarise(across(everything(),  ~mean(.))) %>%  
    pull()
  
  tibble(pops = colnames(stats$Ho),
         Ho = ho,
         He = he,
         Ri = ri,
         FIS = 1 - (ho/he))
}
