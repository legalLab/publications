rm(list = ls())

library(smartsnp)
library(dplyr)
library(ggplot2)
library(svglite)
library(data.table)

# dataset type
species <- "Glyphorynchus spirurus"
project <- "glyphorynchus_"
postfix <- "discosnp_all_final"
fltr <- "_01_00"
strt <- "strata"
grps <- "groups"

# root directory where your data are
data_path <- "/media/tomas/Data2/ddRAD/glyphorynchus/discosnp_glyphorynchus/"
# directory of smartpca
analysis_path <- "/media/tomas/Data2/ddRAD/glyphorynchus/analyses/smartpca/"

# directory where individual data are
res_path <- paste0(data_path, postfix, "/")

# run smartPCA
if (file.exists(paste0(res_path, project, postfix, fltr, '_smartsnp.txt'))) { # execute loop only if taxon exists
  taxon <- read.table(paste0(res_path, project, postfix, fltr, '_smartsnp.txt'), header = TRUE, check.names = FALSE)
  # get samples from input file
  samples <- colnames(taxon) %>%
    as_tibble() %>%
    rename(id = 1)
  
  # read sample to group assignment
  strata <- read.table(paste0(data_path, strt), header = TRUE) %>%
    as_tibble() %>%
    mutate(id = as.character(id))
  
  # assign samples in table to groups
  sample_groups <- samples %>%
    left_join(strata)
  
  # write genotypes without header
  write.table(taxon, file = paste0(analysis_path, project, postfix, fltr), quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE, append = FALSE)
  
  # simple index vectors that determines which samples to remove and which to use for PCA ordination or projection:
  SR <- c() # remove species
  removed_sample <- which(samples$id %in% SR) # Column numbers of samples to remove
  SA <- c() # project species
  ancient_sample <- which(samples$id %in% SA) # Column numbers of samples to project (i.e. aDNA)
  
  # run PCA, PERMANOVA and PERMDISP in one run and assign results to object mvaR
  # without projected or removed taxa 
  mvaR <- smart_mva(snp_data = paste0(analysis_path, project, postfix, fltr), 
                    sample_group = sample_groups$pop, missing_value = 9, 
                    missing_impute = "mean", scaling = "drift", program_svd = "RSpectra", 
                    pc_axes = 3, sample_remove = removed_sample) #, sample_project = ancient_sample, pc_project = c(1:3)) # sample_remove = removed_sample
  
  # show PERMANOVA table
  print(mvaR$test$permanova.global_test)
  
  #show PERMDISP table
  print(mvaR$test$permdisp.global_test) # extract PERMDISP table
  
  # Plot PCA 1 x PCA 2 x PCA 3
  bob <- data.table(mvaR$pca$pca.sample_coordinates)
  bob[, c("PC1mean", "PC2mean", "PC3mean") := .(mean(na.omit(PC1)), mean(na.omit(PC2)), mean(na.omit(PC3))), Group]
  bob <- cbind(samples, bob)
  
  # rearrange bob based on order in groups
  group_order <- read.table(paste0(data_path, grps), header = TRUE, sep = "\t") %>%
    as_tibble()
  bob <- bob %>%
    mutate(Group = factor(Group, levels = group_order$pop)) %>%
    arrange(Group)
  
  PC1 <- paste0('PC1 (', round(mvaR$pca$pca.eigenvalues[2,1], 2), '%)')
  PC2 <- paste0('PC2 (', round(mvaR$pca$pca.eigenvalues[2,2], 2), '%)')
  PC3 <- paste0('PC3 (', round(mvaR$pca$pca.eigenvalues[2,3], 2), '%)')
  
  # plot pca by locality - without projected samples
  p <- ggplot(bob, aes(x=PC1, y=PC2, color=forcats::fct_inorder(Group))) + #color=as.factor(Group)
    geom_point(size=3.5) +
    labs(title=bquote(paste(italic(.(species)), " ", "(", .(postfix), .(fltr), ")")), x=PC1, y=PC2, color="regions")
  
  # output to pdf and svg
  ggsave(plot=p, filename=paste0(analysis_path, project, postfix, fltr, "_pc1_pc2.pdf"), width=9, height=6, bg="transparent", limitsize=FALSE)
  ggsave(plot=p, filename=paste0(analysis_path, project, postfix, fltr, "_pc1_pc2.svg"), device="svg", width=9, height=6, bg="transparent", limitsize=FALSE)
  ggsave(plot=p, filename=paste0(analysis_path, project, postfix, fltr, "_pc1_pc2.png"), device="png", width=9, height=6, bg="transparent", limitsize=FALSE)
  
  # plot pca by locality - without projected samples
  p <- ggplot(bob, aes(x=PC1, y=PC3, color=forcats::fct_inorder(Group))) + 
    geom_point(size=3.5) +
    labs(title=bquote(paste(italic(.(species)), " ", "(", .(postfix), .(fltr), ")")), x=PC1, y=PC3, color="regions")
  
  # output to pdf and svg
  ggsave(plot=p, filename=paste0(analysis_path, project, postfix, fltr, "_pc1_pc3.pdf"), width=9, height=6, bg="transparent", limitsize=FALSE)
  ggsave(plot=p, filename=paste0(analysis_path, project, postfix, fltr, "_pc1_pc3.svg"), device="svg", width=9, height=6, bg="transparent", limitsize=FALSE)
  ggsave(plot=p, filename=paste0(analysis_path, project, postfix, fltr, "_pc1_pc3.png"), device="png", width=9, height=6, bg="transparent", limitsize=FALSE)
  
  # plot pca by locality - without projected samples
  p <- ggplot(bob, aes(x=PC2, y=PC3, color=forcats::fct_inorder(Group))) + 
    geom_point(size=3.5) +
    labs(title=bquote(paste(italic(.(species)), " ", "(", .(postfix), .(fltr), ")")), x=PC2, y=PC3, color="regions")
  
  # output to pdf and svg
  ggsave(plot=p, filename=paste0(analysis_path, project, postfix, fltr, "_pc2_pc3.pdf"), width=9, height=6, bg="transparent", limitsize=FALSE)
  ggsave(plot=p, filename=paste0(analysis_path, project, postfix, fltr, "_pc2_pc3.svg"), device="svg", width=9, height=6, bg="transparent", limitsize=FALSE)
  ggsave(plot=p, filename=paste0(analysis_path, project, postfix, fltr, "_pc2_pc3.png"), device="png", width=9, height=6, bg="transparent", limitsize=FALSE)
  
  # save associated data
  write.table(bob, file = paste0(analysis_path, project, postfix, fltr, "_data.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)
}
