#setwd('/Users/prade/OneDrive - Indiana University/Documents/AgePredictionProject/')
setwd('/Users/ppugale/OneDrive - Indiana University/Documents/AgePredictionProject/')
#install.packages("rempsyc")
#devtools::install_github("PradoVarathan/AssociatoR")

library(AssociatoR)
library(rempsyc)
library(dplyr)
library(ggpubr)
prep_factors = function(data, targets, add_on){

  if(add_on == ""){
    data$Sex = as.factor(data$Sex)
  }
  data$APOE4 = as.factor(data$APOE4)
  colnames(data)[which(colnames(data) == 'VISCODE2')] = 'timepoint'
  data$delta_age = data$Corrected_Predicted_Age - data$Actual_Age

  for(target in targets){
    data = data[which(data$RID %in% unique(data$RID[which(data$timepoint == 'bl')])),]
    for(sub in unique(data[,'RID'])){
      data[which(data$RID == sub),target] = data[which(data$RID == sub & data$timepoint == 'bl'), target]
    }
  }
  return(data)

}

#Long
add_on = ""
targets = c("Corrected_Predicted_Age")
thresh_year = 2
#folder = 'LongitudinalDeltaAge'
folder = 'LongitudinalPredictedAge'
for(add_on in c("","__FEMALE","__MALE")){
  mri_cross = read.csv(paste0('Data/CrossSectionalDataPrepared/MRI_ALL_ADNI',add_on,'.csv'))
  csf_cross = read.csv(paste0('Data/CrossSectionalDataPrepared/CSF_ALL_ADNI',add_on,'.csv'))
  plasma_cross = read.csv(paste0('Data/CrossSectionalDataPrepared/Plasma_ALL_ADNI',add_on,'.csv'))
  fdg_cross = read.csv(paste0('Data/CrossSectionalDataPrepared/FDG_ALL_ADNI',add_on,'.csv'))
  cogn_cross = readRDS(paste0('Data/CrossSectionalDataPrepared/Cognitivie_Scores',add_on,'_ALL.RDS'))





  # MRI ROIs ----------------------------------------------------------------


    mri_cross = prep_factors(mri_cross,c('HippVol', 'entorhinal_thickness'), add_on)
    if(add_on == ""){
      formula = '1 + Actual_Age + Sex + APOE4 + BMI + Phenotype * year + ( 1 + year | RID )'
    }else{
      formula = '1 + Actual_Age + APOE4 + BMI + Phenotype * year + ( 1 + year | RID )'
    }
    for(tar in c('HippVol', 'entorhinal_thickness')){

      formula_tar = gsub("Phenotype", tar, formula, fixed = TRUE)
      mri_out = run_longitudinal_interaction_association(data = mri_cross,
                                                         targets = targets,
                                                         threshold_years  = thresh_year, plot_residuals = T,
                                                         formula = formula_tar)
      sep_idx = grep("[+]",strsplit(formula_tar,"")[[1]])[floor(length(grep("[+]",strsplit(formula_tar,"")[[1]]))/2)]
      formula_tar = paste0(paste(strsplit(formula_tar,"")[[1]][1:sep_idx],collapse = "")," \n",paste(strsplit(formula_tar,"")[[1]][(sep_idx+1):length(strsplit(formula_tar,"")[[1]])],collapse = ""))
      ggsave(paste0('LongitudinalPlots/',folder,'/MRI_BMI_Long_',tar,'_',add_on,thresh_year+1,'.png'),plot = annotate_figure(ggarrange(mri_out$effect_plots[[1]],
                                                                                                                  mri_out$residual_plots[[1]],
                                                                                                                  ncol = 2, nrow = 1), top = text_grob(paste0("MRI (",thresh_year+1," years)\n",targets[1]," ~ ",formula_tar),
                                                                                                                                                       color = "red", face = "bold", size = 14)))
    }






  # CSF ---------------------------------------------------------------------



    csf_cross = prep_factors(csf_cross, c('ABETA','TAU','PTAU'), add_on)
    csf_cross$ABETA = log(csf_cross$ABETA)
    csf_cross$TAU = log(csf_cross$TAU)
    csf_cross$PTAU = log(csf_cross$PTAU)
    #removing subject with only baseline timepoint
    csf_bl_delta_age_subjects_more_than_one = names(table(csf_cross$RID)[table(csf_cross$RID) != 1])
    csf_cross = csf_cross %>% filter(RID %in% csf_bl_delta_age_subjects_more_than_one)

    formula = '1 + Actual_Age + APOE4 + Phenotype * year + BMI + ( 1 + year | RID )'
    for(tar in  c('ABETA','TAU','PTAU')){

      formula_tar = gsub("Phenotype", tar, formula, fixed = TRUE)
      csf_out = run_longitudinal_interaction_association(data = csf_cross,
                                                         targets = targets,
                                                         threshold_years  = thresh_year, plot_residuals = T,
                                                         formula = formula_tar)

      sep_idx = grep("[+]",strsplit(formula_tar,"")[[1]])[floor(length(grep("[+]",strsplit(formula_tar,"")[[1]]))/2)]
      formula_tar = paste0(paste(strsplit(formula_tar,"")[[1]][1:sep_idx],collapse = "")," \n",paste(strsplit(formula_tar,"")[[1]][(sep_idx+1):length(strsplit(formula_tar,"")[[1]])],collapse = ""))

      ggsave(paste0('LongitudinalPlots/',folder,'/CSF_Long_',tar,'_',add_on,thresh_year+1,'.png'),plot = annotate_figure(ggarrange(csf_out$effect_plots[[1]],
                                                                                                              csf_out$residual_plots[[1]],
                                                                                                              ncol = 2, nrow = 1), top = text_grob(paste0("CSF (",thresh_year+1," years)\n ",targets[1]," ~ ",formula_tar),
                                                                                                                                                   color = "red", face = "bold", size = 14)))
    }




  # PLASMA ------------------------------------------------------------------


    plasma_cross = prep_factors(plasma_cross, c('AB40','AB42'),add_on)
    formula = '1 + Actual_Age + APOE4 + BMI + Phenotype * year + ( 1 + year | RID )'
    for(tar in  c('AB40','AB42')){

      formula_tar = gsub("Phenotype", tar, formula, fixed = TRUE)
      plasma_out = run_longitudinal_interaction_association(data = plasma_cross,
                                                            targets =targets,
                                                            threshold_years  = thresh_year,plot_residuals = T,
                                                            formula = formula_tar)
      sep_idx = grep("[+]",strsplit(formula_tar,"")[[1]])[floor(length(grep("[+]",strsplit(formula_tar,"")[[1]]))/2)]
      formula_tar = paste0(paste(strsplit(formula_tar,"")[[1]][1:sep_idx],collapse = "")," \n",paste(strsplit(formula_tar,"")[[1]][(sep_idx+1):length(strsplit(formula_tar,"")[[1]])],collapse = ""))
      ggsave(paste0('LongitudinalPlots/',folder,'/PLASMA_Long_',tar,'_',add_on,thresh_year+1,'.png'),plot = annotate_figure(ggarrange(plasma_out$effect_plots[[1]],
                                                                                                                 plasma_out$residual_plots[[1]],
                                                                                                                 ncol = 2, nrow = 1), top = text_grob(paste0("Plasma (",thresh_year+1," years)\n ",targets[1]," ~ ", formula_tar),
                                                                                                                                                      color = "red", face = "bold", size = 14)))
    }


  # FDG MEASURES ------------------------------------------------------------


    fdg_cross = prep_factors(fdg_cross,c('Cing_Mean','Ang_Mean','LAng_Mean','RAng_Mean','Temp_Mean','RTemp_Mean','LTemp_Mean'), add_on)
    fdg_cross = fdg_cross[which(fdg_cross$RID %in% names(table(fdg_cross$RID))[table(fdg_cross$RID) > 1]),]
    if(add_on == ""){
      formula = '1 + Sex + Actual_Age + BMI + APOE4 + Phenotype * year + ( 1 + year | RID )'
    }else{
      formula = '1 + Actual_Age + APOE4 + BMI + Phenotype * year + ( 1 + year | RID )'
    }
    for(tar in  c('Cing_Mean','Ang_Mean','LAng_Mean','RAng_Mean','Temp_Mean','RTemp_Mean','LTemp_Mean')){

      formula_tar = gsub("Phenotype", tar, formula, fixed = TRUE)
      fdg_out = run_longitudinal_interaction_association(data = fdg_cross,
                                                         targets = targets,
                                                         threshold_years  = thresh_year, plot_residuals = T,
                                                         formula = formula_tar)
      sep_idx = grep("[+]",strsplit(formula_tar,"")[[1]])[floor(length(grep("[+]",strsplit(formula_tar,"")[[1]]))/2)]
      formula_tar = paste0(paste(strsplit(formula_tar,"")[[1]][1:sep_idx],collapse = "")," \n",paste(strsplit(formula_tar,"")[[1]][(sep_idx+1):length(strsplit(formula_tar,"")[[1]])],collapse = ""))

      ggsave(paste0('LongitudinalPlots/',folder,'/FDG_Long_',tar,'_',add_on,thresh_year+1,'.png'),plot = annotate_figure(ggarrange(fdg_out$effect_plots[[1]],
                                                                                                              fdg_out$residual_plots[[1]],
                                                                                                              ncol = 2, nrow = 1), top = text_grob(paste0("FDG Scores (",thresh_year+1," years)\n ",targets[1]," ~ ",formula_tar),
                                                                                                                                                   color = "red", face = "bold", size = 14)))
    }

  # COGNITION COMPOSITE SCORES ----------------------------------------------


    cogn_cross_main = merge(prep_factors(cogn_cross$`_ADNI_MEM`,'ADNI_MEM', add_on)[,c('RID','timepoint','ADNI_MEM')],
                            prep_factors(cogn_cross$`_ADNI_EF`,'ADNI_EF', add_on)[,c('RID','timepoint','ADNI_EF')],
                            by = c('RID','timepoint'))
    cogn_cross_main = merge(cogn_cross_main, prep_factors(cogn_cross$`_ADNI_LAN`,'ADNI_LAN', add_on),by = c('RID','timepoint'))
    if(add_on == ""){
      formula = '1 + Sex + Actual_Age + BMI + APOE4 + Phenotype * year + ( 1 + year | RID )'
    }else{
      formula = '1 + Actual_Age + BMI + APOE4 + Phenotype * year + ( 1 + year | RID )'
    }
    for(tar in  c('ADNI_MEM','ADNI_EF','ADNI_LAN')){

      formula_tar = gsub("Phenotype", tar, formula, fixed = TRUE)
      cogn_p_vals = run_longitudinal_interaction_association(data = cogn_cross_main,
                                                             targets = targets,
                                                             threshold_years  = thresh_year, plot_residuals = T,
                                                             formula = formula_tar)

      sep_idx = grep("[+]",strsplit(formula_tar,"")[[1]])[floor(length(grep("[+]",strsplit(formula_tar,"")[[1]]))/2)]
      formula_tar = paste0(paste(strsplit(formula_tar,"")[[1]][1:sep_idx],collapse = "")," \n",paste(strsplit(formula_tar,"")[[1]][(sep_idx+1):length(strsplit(formula_tar,"")[[1]])],collapse = ""))

      ggsave(paste0('LongitudinalPlots/',folder,'/COGN_Long_',tar,'_',add_on,thresh_year+1,'.png'),plot = annotate_figure(ggarrange(cogn_p_vals$effect_plots[[1]],
                                                                                                               cogn_p_vals$residual_plots[[1]],
                                                                                                               ncol = 2, nrow = 1), top = text_grob(paste0("Cognitive Scores (",thresh_year+1," years)\n ",targets[1]," ~ ",formula_tar),
                                                                                                                                                    color = "red", face = "bold", size = 14)))

  }
}


