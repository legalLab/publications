#----------------------------------------
##### Ecological Niche Models of #####
###### Saguinus midas #####
#### National Institute of Amazonian Research - Ecology Program #################

#-------------------------------------
##### THIAGO CAVALCANTE_2020_08_04 #####

#Based and adapted from 
#Luisa Maria Diele-Viegas course_Ecological Niche Modeling.
# For more information, visit:
# https://github.com/thiago-cav/Disciplina_ENMs
rm(list = ls())


#######################################################################
##### Loading packages
#######################################################################
library(rgbif)
library(mapview)
library(scrubr)
library(sp)
library(dplyr)
library(dismo)
library(raster)
library(ggplot2)
library(maptools)
library(scrubr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(sf)
library(jsonlite)
library(mapdata)
library(rgdal)
library(rJava)
library(spThin)
library(psych)
library(beepr)
library(kernlab)
library(tidyverse)
library(wesanderson)
library(vegan)
library(rgeos)
library(sp)

# Using stack function to read variables in all bands
predictors <- 
  stack(paste0("predictors"))
predictors@layers

# reading occurrence data
occ <- readr::read_csv("occ.csv")
occ <- dplyr::select(occ, species, lat, lon)


##### Modeling procedures #####

# definindo os parametros a priori
replica <- 100 #numeros de replica??es que serao usadas no modelo
partition <- .7 #definindo j? aqui que serao 70% treino, 30% teste!

# algoritmos = vamos usar um loop para rodar todos de uma vez!
for(i in occ$species %>% unique){ # para cada especie
  
  # diretorio
  dir.create(i)
  setwd(i)
  
  # informacao importante!
  paste0("Preparing data for modeling ", i, " in ", getwd()) %>% print
  
  # objeto para avaliacao 
  eval_species <- tibble::tibble()
  
  # selecionando dados de presenca e (pseudo)ausencia 
  # dados de presenca
  pr_specie <- occ %>% 
    dplyr::filter(species == i) %>% 
    dplyr::select(lon, lat) %>% 
    dplyr::mutate(id = seq(nrow(.)))
  
  # dados de pseudoausencia 
  pa_specie <- dismo::randomPoints(mask = var, n = 1000) %>% 
    tibble::as_tibble() %>%
    dplyr::rename(lon = x, lat = y) %>% 
    dplyr::mutate(id = seq(nrow(.)))
  
  # mudando o diretorio de trabalho novamente
  dir.create("00_replicas")
  setwd("00_replicas")
  
  # replicas
  for(r in replica %>% seq){	# numero de replicas do modelo 
    
    # objeto para a avaliacao 
    eval_algorithm <- tibble::tibble()
    
    # particionando os dados com base na nossa selecao (70% treino, 30% teste)
    # dados de presenca
    pr_sample_train <- pr_specie %>% 
      dplyr::sample_frac(partition) %>% 
      dplyr::select(id) %>% 
      dplyr::pull()
    
    # dados de pseudo ausencia
    pa_sample_train <- pa_specie %>% 
      dplyr::sample_frac(partition) %>% 
      dplyr::select(id) %>% 
      dplyr::pull()
    
    # dados de treino e teste !
    # treino
    train <- dismo::prepareData(x = var, 
                                p = pr_specie %>% 
                                  dplyr::filter(id %in% pr_sample_train) %>% 
                                  dplyr::select(lon, lat), 
                                b = pa_specie %>% 
                                  dplyr::filter(id %in% pa_sample_train) %>% 
                                  dplyr::select(lon, lat)) %>% na.omit
    
    # teste
    test <- dismo::prepareData(x = var, 
                               p = pr_specie %>% 
                                 dplyr::filter(!id %in% pr_sample_train) %>% 
                                 dplyr::select(lon, lat), 
                               b = pa_specie %>% 
                                 dplyr::filter(!id %in% pa_sample_train) %>% 
                                 dplyr::select(lon, lat)) %>% na.omit
    
    
    # Ajuste do modelo
    # informacao
    print(paste("Models fitting to", i, "replica", r, "of", replica))
    
    # SVM (machine learning) - Dados de presen?a e plano de fundo
    SVM <- kernlab::ksvm(x = pb ~ ., data = train)
    
    # Lista com todos os algoritmos :)
    fit <- list(svm = SVM)
    
    # Previsoes (mais uma vez no loop)
    for(a in seq(fit)){
      
      # informacao 
      print(paste("Model predict algorithm", fit[a] %>% names))
      
      # previsao do modelo 
      model_predict <- dismo::predict(var2, fit[[a]], progress = "text")
      
      # exporta os valores da previsao do modelo 
      raster::writeRaster(x = model_predict, 
                          filename = paste0("enm_", i, "_", fit[a] %>% names, "_r", ifelse(r < 10, paste0("0", r), r)), 
                          format = "GTiff", 
                          options = c("COMPRESS=DEFLATE"), 
                          overwrite = TRUE)
      
      # avaliacao do modelo
      eval <- dismo::evaluate(p = test %>% dplyr::filter(pb == 1) %>% dplyr::select(-pb), 
                              a = test %>% dplyr::filter(pb == 0) %>% dplyr::select(-pb), 
                              model = fit[[a]])
      
      # indices de avaliacao 
      id_eval_spec_sens <- which(eval@t == dismo::threshold(eval, "spec_sens"))
      tss_spec_sens <- eval@TPR[id_eval_spec_sens] + eval@TNR[id_eval_spec_sens] - 1
      
      # dados da avaliacao 
      eval_data <- tibble::tibble(species = i, 
                                  replica = r, 
                                  algorithm = fit[a] %>% names, 
                                  thr_max_spec_sens = dismo::threshold(eval, "spec_sens"),
                                  tss_spec_sens = tss_spec_sens,
                                  auc = eval@auc, 
                                  file = paste0("enm_", i, "_", fit[a] %>% names, "_r", ifelse(r < 10, paste0("0", r), r), ".tif"))
      
      # combina a avaliacao dos modelos 
      eval_algorithm <- dplyr::bind_rows(eval_algorithm, eval_data)
      
    } 
    
    # combina as avaliacoes
    eval_species <- dplyr::bind_rows(eval_species, eval_algorithm)
    
  }
  
  # exporta as avaliacoes numa nova pasta 
  setwd("..")
  
  dir.create("01_evaluation")
  setwd("01_evaluation")
  dir.create("00_raw")
  setwd("00_raw")
  
  readr::write_csv(eval_species, paste0("eval_", i, ".csv"))
  
  # diretorio de trabalho
  setwd(".."); setwd(".."); setwd("..") 
  
  # notifica que a analise acabou porque nao somos obrigados a ficar conferindo toda hora :)
  beepr::beep(8)
  
} 


##### Models evaluations #####

for(i in eva$species %>% unique){
  
  # diretorio de trabalho
  setwd(i); setwd("01_evaluation")
  
  # seleciona a especie 
  setwd("00_raw")
  eva_sp <- eva %>% 
    dplyr::filter(species == i)
  
  # tabelas
  setwd("..")
  dir.create("01_tables")
  setwd("01_tables")
  
  # tabela para avaliar os modelos por TSS e AUC
  eva_table <- eva_sp %>% 
    dplyr::mutate(species = species %>% stringr::str_to_title() %>% stringr::str_replace("_", " ")) %>% 
    dplyr::group_by(species, algorithm) %>% 
    dplyr::summarise(tss_mean = mean(tss_spec_sens) %>% round(3), 
                     tss_sd = sd(tss_spec_sens) %>% round(3),
                     auc_mean = mean(auc) %>% round(3), 
                     auc_sd = sd(auc) %>% round(3))
  eva_table
  
  # exportando a avaliacao dos modelos
  readr::write_csv(eva_table, paste0("evaluation_summary_table_", i, ".csv"))
  
  # boxplots
  # definindo diretorio de trabalho
  setwd("..")
  dir.create("02_boxplot")
  setwd("02_boxplot")
  
  for(j in c("tss_spec_sens", "auc")){
    
    # informacao
    print(paste(i, j))
    
    # plot dos boxplots referentes aos diferentes algoritmos  
    ggplot(data = eva_sp) + 
      aes_string(x = "algorithm", y = j, color = "algorithm") +
      geom_boxplot(size = .5, fill = "gray90", color = "black") +
      geom_jitter(width = 0.2, size = 4, alpha = .7) +
      scale_color_manual(values = wesanderson::wes_palette(name = "Darjeeling1", n = eva$algorithm %>% unique %>% length, 
                                                           type = "continuous")) +
      labs(x = "Algorithms", 
           y = stringr::str_to_upper(j) %>% stringr::str_replace("_", " "), 
           title = i %>% stringr::str_to_title() %>% stringr::str_replace("_", " ")) + 
      ylim(c(-.01, 1.05)) + 
      theme_bw() +
      geom_hline(yintercept = ifelse(j == "tss_spec_sens", .5, .75), color = "red") +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold.italic", size = 20), 
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 15), 
            axis.title = element_text(size = 17))
    ggsave(paste0("boxplot_jitter_an_", j, "_", i, ".tiff"), he = 20, wi = 30, un = "cm", dpi = 300)
    
  }
  
}



##### Creating the consensus model of the 100 replicates #####
##### using only the models above the mean AUC of the 100 replicates #####
auc_limit <- .898

# algoritmos
alg <- eva$algorithm %>% unique
alg

# fazendo o consenso
for(i in eva$species %>% unique){
  
  # informacao
  print(paste("Consenso de", i))
  
  # selecao de modelos = somente aqueles com AUC maior ou igual a 0.75
  eva_i <- eva %>% 
    dplyr::filter(species == i, 
                  auc >= auc_limit, 
                  algorithm %in% alg)
  
  # importando os modelos
  enm <- eva_i %>% 
    dplyr::select(file) %>% 
    dplyr::mutate(file = paste0(i, "/00_replicas/", file)) %>% 
    dplyr::pull() %>% 
    raster::stack()
  
  # AUC
  auc <- eva_i %>% 
    dplyr::select(auc) %>% 
    dplyr::mutate(auc = (auc - .5) ^ 2) %>% 
    dplyr::pull()
  
  # padronizacao 
  print("Pode demorar... mas vai dar bom!")
  enm_st <- enm %>% 
    values %>% 
    vegan::decostand("range", na.rm = TRUE)
  print("N?o disse? Sucesso!")
  
  # consenso da media ponderada 
  ens <- enm[[1]]
  ens[] <- apply(enm_st, 1, function(x){sum(x * auc) / sum(auc)})
  
  # diretorio de trabalho 
  setwd(path)
  dir.create("S_midas_344NRPfinal_consensus")
  setwd("S_midas_344NRPfinal_consensus")
  
  # exporta o ensemble 
  raster::writeRaster(x = ens, 
                      filename = paste0("consenso_media_ponderada_", i), 
                      format = "GTiff", 
                      options = c("COMPRESS=DEFLATE"), 
                      overwrite = TRUE)
  
} 

