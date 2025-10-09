rm(list = ls())

library(vcfR)
library(vcf2others)

# classify relatedness based on coefficients
classify_relatedness <- function(R0, R1, KING) {
  case_when(
    KING > 0.354 ~ "Twins",
    KING > 0.177 & R0 < 0.1 ~ "Parent_Offspring",
    KING > 0.177 ~ "Full_Sibs",
    KING > 0.0884 ~ "Half_Sibs", # also grandparents_grandchildren and aunts/uncles_nephews
    KING > 0.0442 ~ "1_degree_cousins",
    KING > 0.0221 ~ "2_degree_cousins",
    TRUE ~ "Unrelated"
  )
}

# dataset type
species <- "Glyphorynchus spirurus"
project <- "glyphorynchus_"
postfix <- "discosnp_all"
fltr <- "_01_00"

# directory where your data are
data_path <- "/media/tomas/Data2/ddRAD/glyphorynchus/discosnp_glyphorynchus/"
res_path <- paste0(data_path, postfix, "/")
# directory where your data are
result_path <- "/media/tomas/Data2/ddRAD/glyphorynchus/analyses/ngsrelate/"


### extract families from multi VCF ###
vcf <- read.vcfR(paste0(res_path, project, postfix, fltr, "_oneSNP.vcf.gz"))
project_info <- get_vcf_group_info(vcf, data_path)

for(group in project_info$groups){
  vcf1 <- vcf_extract_pops(vcf, project_info$indiv_group, group)
  write.vcf(vcf1, file = paste0(result_path, group, "_oneSNP.vcf.gz"))
}

### run NGSRelate - external command ###
for(group in project_info$groups){
  cmd <- paste0("/home/tomas/bin/NgsRelate/ngsRelate -h ", result_path, group, "_oneSNP.vcf.gz -O ", result_path, group, ".res")
  system(cmd)
}

### rename NGSRelated samples ###
for(group in project_info$groups){
  vcf <- read.vcfR(paste0(result_path, group, "_oneSNP.vcf.gz"))
  ngsrelate_res <- read.table(paste0(result_path, group, ".res"), header = TRUE, sep = "\t")
  
  # get names of samples from VCF column names
  individuals <- colnames(vcf@gt)[-1]
  
  # function to generate unique pairs and return them as a data frame
  generate_unique_pairs <- function(values) {
    # initialize vectors to store the pairs
    col1 <- c()
    col2 <- c()
    
    # loop through the values to create pairs
    for (i in 1:(length(values) - 1)) {
      for (j in (i + 1):length(values)) {
        # Append values to the vectors
        col1 <- c(col1, values[i])
        col2 <- c(col2, values[j])
      }
    }
    
    # create a data frame from the vectors
    df <- data.frame(a = col1, b = col2)
    return(df)
  }
  
  # generate unique pairs and return as a data frame
  unique_pairs <- generate_unique_pairs(individuals)
  
  # substitute columns in NGSrelate table
  ngsrelate_res$a <- unique_pairs$a
  ngsrelate_res$b <- unique_pairs$b
  
  # classify pairs of individuals
  ngsrelate_res$relationship <- classify_relatedness(ngsrelate_res$R0, ngsrelate_res$R1, ngsrelate_res$KING)
  
  # output new NGSrelate table
  write.table(ngsrelate_res, file = paste0(result_path, group, ".res"), row.names = FALSE, quote = FALSE, sep = "\t")
}

### run all individuals together ###
vcf <- read.vcfR(paste0(res_path, project, postfix, fltr, "_oneSNP.vcf.gz"))

### run NGSRelate - external command ###
cmd <- paste0("/home/tomas/bin/NgsRelate/ngsRelate -h ", res_path, project, postfix, fltr, "_oneSNP.vcf.gz -O ", result_path, project, postfix, fltr, ".res")
system(cmd)

### rename NGSRelated samples ###
ngsrelate_res <- read.table(paste0(result_path, project, postfix, fltr, ".res"), header = TRUE, sep = "\t")
  
# get names of samples from VCF column names
individuals <- colnames(vcf@gt)[-1]

# function to generate unique pairs and return them as a data frame
generate_unique_pairs <- function(values) {
  # initialize vectors to store the pairs
  col1 <- c()
  col2 <- c()
  
  # loop through the values to create pairs
  for (i in 1:(length(values) - 1)) {
    for (j in (i + 1):length(values)) {
      # Append values to the vectors
      col1 <- c(col1, values[i])
      col2 <- c(col2, values[j])
    }
  }
  
  # create a data frame from the vectors
  df <- data.frame(a = col1, b = col2)
  return(df)
}

# generate unique pairs and return as a data frame
unique_pairs <- generate_unique_pairs(individuals)

# substitute columns in NGSrelate table
ngsrelate_res$a <- unique_pairs$a
ngsrelate_res$b <- unique_pairs$b

# classify pairs of individuals
ngsrelate_res$relationship <- classify_relatedness(ngsrelate_res$R0, ngsrelate_res$R1, ngsrelate_res$KING)

# output new NGSrelate table
write.table(ngsrelate_res, file = paste0(result_path, project, postfix, fltr, ".res"), row.names = FALSE, quote = FALSE, sep = "\t")


