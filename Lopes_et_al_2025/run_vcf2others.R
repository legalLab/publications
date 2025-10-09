rm(list = ls())

library(vcfR)
library(vcf2others)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(colorRamps)
library(svglite)

#############################
# final filter - prepare VCF for analyses
#############################

# directory where your data are
data_path <- "/media/tomas/Data2/ddRAD/glyphorynchus/discosnp_glyphorynchus/"
species <- "Glyphorynchus spirurus"
project <- "glyphorynchus_"
postfix <- "discosnp_all"
postfix_sub <- "discosnp_all_final"
fltr <- "_01_00"

vcf <- read.vcfR(paste0(data_path, project, postfix, ".vcf.gz"))

# read individuals to include/exclude
indivs <- read.table(paste0(data_path, "indivs_b"), header = TRUE)$id %>%
  as.character()
# check if all indivs are in vcf_names
if (any(!(indivs %in% colnames(vcf@gt)[-1]))) stop(paste("Some individuals in list not in VCF"))

# directory where your results will be written
res_path_sub <- paste0(data_path, postfix_sub, "/")
# create results directory if does not exist
if (!file.exists(res_path_sub)) dir.create(res_path_sub)

# remove offending individuals from VCF and filter
vcf1 <- vcf_extract_indivs(vcf, indivs, whitelist = FALSE) %>%
  vcf_filter_oneSNP() %>%
  vcf_filter_maf(.02) %>% 
  vcf_filter_coverage(6) %>%
  vcf_filter_missingness(0.01)

print(vcf1)
write.vcf(vcf1, file = paste0(res_path_sub, project, postfix_sub, fltr, "_oneSNP.vcf.gz"))
assess_vcf_missing_data(vcf1, data_path, res_path_sub, project, postfix_sub, fltr, species)
vcf_stats(vcf1, res_path_sub, paste0(project, postfix_sub, fltr))
project_info <- get_vcf_group_info(vcf1, data_path)
vcf2smartsnp(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_smartsnp.txt'))
vcf2migrate(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_migrate.txt'))
vcf2arlequin(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_arlequin.arp'))
vcf2structure(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_structure.str'))
vcf2structure(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_fstructure.str'), method = "F")
vcf2genepop(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_genepop.gen'))


#############################
# subset final filtered data
#############################

postfix_sub <- "discosnp_all_final_sub"
# directory where your results will be written
res_path_sub <- paste0(data_path, postfix_sub, "/")
# create results directory if does not exist
if (!file.exists(res_path_sub)) dir.create(res_path_sub)

# subsample SNPs (default 10000 SNPs)
vcf <- vcf1
vcf1 <- vcf_sub_loci(vcf, n_loci = 10000)
print(vcf1)
write.vcf(vcf1, file = paste0(res_path_sub, project, postfix_sub, fltr, "_oneSNP.vcf.gz"))
assess_vcf_missing_data(vcf1, data_path, res_path_sub, project, postfix_sub, fltr, species)
vcf_stats(vcf1, res_path_sub, paste0(project, postfix_sub, fltr))
project_info <- get_vcf_group_info(vcf1, data_path)
vcf2smartsnp(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_smartsnp.txt'))
vcf2migrate(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_migrate.txt'))
vcf2arlequin(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_arlequin.arp'))
vcf2structure(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_structure.str'))
vcf2structure(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_fstructure.str'), method = "F")
vcf2genepop(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_genepop.gen'))

# subsample SNPs (default 10000 SNPs) - each filter 5x
fltrr <- fltr
for (i in 1:5) {
  vcf1 <- vcf_sub_loci(vcf, n_loci = 10000)
  print(vcf1)
  fltr <- paste0(fltrr, "_r", i)
  write.vcf(vcf1, file = paste0(res_path_sub, project, postfix_sub, fltr, "_oneSNP.vcf.gz"))
  assess_vcf_missing_data(vcf1, data_path, res_path_sub, project, postfix_sub, fltr, species)
  vcf_stats(vcf1, res_path_sub, paste0(project, postfix_sub, fltr))
  project_info <- get_vcf_group_info(vcf1, data_path)
  vcf2smartsnp(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_smartsnp.txt'))
  vcf2migrate(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_migrate.txt'))
  vcf2arlequin(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_arlequin.arp'))
  vcf2structure(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_structure.str'))
  vcf2structure(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_fstructure.str'), method = "F")
  vcf2genepop(vcf1, project_info$indiv_group, project_info$groups, out_file = paste0(res_path_sub, project, postfix_sub, fltr, '_genepop.gen'))
}
