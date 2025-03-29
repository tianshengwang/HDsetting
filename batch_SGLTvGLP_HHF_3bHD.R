rm(list=ls())
setwd("/local/projects/medicare/sglt_cvs/data/full/analysis")
library(MASS)
library(grf)
library(tidyverse)
library(rlang)
library(rlist)
library(plyr)
library(caret)
library(caTools)
library(randomForest)
library(data.table)
library(grid)
library(broom)
library(rstatix)
library(knitr)
library(ggplot2)
library(ggridges)
library(glmnet)
library(pROC)

source("/local/projects/medicare/DPP4i_HTE/programs/macros/best_tree_MSegar.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/iCF_GG_toolbox.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/iCF_TREE_build.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/iCF_PARENT_node.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/iCF_PRE_majority.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/iCF_MAJORITY_VOTE.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/iCF_SUBGROUP_DECISION.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/iCF_SUBGROUP_PIPELINE.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/iCF_CV.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/iCF_SUBGROUP_ANALYSIS.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/sim_Truth_tree.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/sim_GenSimData.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/iCF_SUBGROUP_MODEL.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/iCF_HD.R")

#####################################################
#####################################################
#
# STEP 1: HD feature identificaiton (done in SAS)   #
#
#####################################################
#####################################################
#-----------------------#
analysis ="sen7"        #
#-----------------------#
drug1="sglt"; drug2="glp"; pct_inter=0;

if (analysis=="primary"){
c=0.01; n=200; dxgroup <<- 3; atcgroup <<- 3; Lcutoff=20; pcttop=0.95; mls_D5=115; mls_D4=95;  mls_D3=70;  mls_D2=50;

} else if (analysis=="sen1"){
c=0.01; n=200; dxgroup <<- 3; atcgroup <<- 3; Lcutoff=20; pcttop=0.99; mls_D5=300;  mls_D4=190;  mls_D3=140;  mls_D2=110;
  
} else if (analysis=="sen2"){
c=0.01; n=200; dxgroup <<- 3; atcgroup <<- 3; Lcutoff=35; pcttop=0.9;  mls_D5=90;  mls_D4=70;  mls_D3=50;  mls_D2=50;
  
} else if (analysis=="sen3"){
c=0.01; n=100; dxgroup <<- 3; atcgroup <<- 3; Lcutoff=25; pcttop=0.95; mls_D5=115; mls_D4=90;  mls_D3=70;  mls_D2=50;

} else if (analysis=="sen4"){
c=0.02; n=200; dxgroup <<- 3; atcgroup <<- 3; Lcutoff=20; pcttop=0.95; mls_D5=115;  mls_D4=90;  mls_D3=70;  mls_D2=50;

} else if (analysis=="sen5"){
c=0.05; n=200; dxgroup <<- 3; atcgroup <<- 3; Lcutoff=20; pcttop=0.95; mls_D5=110;  mls_D4=90;  mls_D3=70;  mls_D2=50;

} else if (analysis=="sen6"){
c=0.01; n=200; dxgroup <<- 4; atcgroup <<- 3; Lcutoff=30; pcttop=0.95; mls_D5=120;  mls_D4=95;  mls_D3=75;  mls_D2=50;

} else if (analysis=="sen7"){
c=0.01; n=200; dxgroup <<- 3; atcgroup <<- 4; Lcutoff=25; pcttop=0.95; mls_D5=115;  mls_D4=90;  mls_D3=75;  mls_D2=50;

} 


#4/27/2024 add senstivity analysis 7 for ATC4
Train_0 = haven::read_sas(paste0("hdps_" #"hdps_"
                                 , 
                                 drug1, "v", drug2, "_p",c*100, 
                                #"n",  
                                "_n",
                                n,
                                # "i",
                                "_i",
                                 dxgroup, 
                                #"a",
                                "_a",
                                 atcgroup,
                                #"_v2",
                                 ".sas7bdat"))
str(Train_0)
nrow(Train_0)
Train_0_date <- Train_0 %>% 
            dplyr::filter(FillDate2 <=  as.Date("2019-12-31") & #21914 &  #31DEC2019
                          IndexDate >= (as.Date("2015-10-15") + 365) & #(20376 + 365) &#20376 & #15OCT2015 + 365 daysBL 
                          FillDate2 <= (as.Date("2019-12-31") - 365*2) & #(21914-365*2) & #& #31DEC2019 due to limited sample size
                          is.na(FillDate2)==F &
                          excludeFlag_preFill2Initiator ==0 &
                          excludeFlag_sameDayInitiator==0 &
                          excludeFlag_prevalentUser==0 
                          ) %>%
            dplyr::rename(HFPRIMARY_ICD10DX_date=outcome_date)#for 3bHD setting data generated from SAS

table(Train_0_date$HFPRIMARY_ICD10DX_date, useNA = "always")
is.na(Train_0_date$HFPRIMARY_ICD10DX_date)

Train1 <- Train_0_date %>% 
          transform(min_CensorF2_730= ifelse(is.na(HFPRIMARY_ICD10DX_date)==T, 
                                             pmin(censorDate_ITT, FillDate2 + 365*2),
                                             pmin(censorDate_ITT, FillDate2 + 365*2, HFPRIMARY_ICD10DX_date)
                                                  )) %>%
          dplyr::mutate(
            HHFfu_2yr = as.Date(min_CensorF2_730)     - FillDate2,
            HHF   = ifelse( is.na( HFPRIMARY_ICD10DX_date) != T & #has event date
                                   HFPRIMARY_ICD10DX_date  >= FillDate2, #event date is later than Filldate2
                            1, 0),
            age = as.numeric( cut(age, c(65,70,75,80,85,Inf) ,
                                  labels=c("65<age<=70 ","70<age<=75","75<age<=80","80<age<=85", "age>85")
            )) 
            
          ) %>%
          dplyr::mutate(race = case_when(race == "1" ~ 1,
                                         race == "2" ~ 2,
                                         !(race %in% c("1", "2")) ~  3)
          )

nrow(Train1) #15388
table(Train1$HHF, useNA = "always")

             
Train2 <- Train1 %>%
          dplyr::mutate(HHF_2yr_2yr =  ifelse(  HHF==1 & # event date is later than the start of follow up 
                                                    (#not dead by 12/31/2019
                                                    (is.na(death_dt) == T & 
                                                     HFPRIMARY_ICD10DX_date <=  FillDate2+730 &
                                                     HFPRIMARY_ICD10DX_date <=  censorDate_ITT    ) |
                                                      #dead by 12/31/2019
                                                    (is.na(death_dt) != T & 
                                                     HFPRIMARY_ICD10DX_date <=  FillDate2+730 &
                                                     HFPRIMARY_ICD10DX_date <=  censorDate_ITT & 
                                                     HFPRIMARY_ICD10DX_date <=  death_dt
                                                     ) 
                                                    )
                                               , 
                                               1, 
                                               0)                                                      

                        )



table(Train2$HHF_2yr_2yr, useNA = "always")
#######################################################################
#  STEP 1: READ DATA from Medicare or MarketScan: Fx to get CF-ready  #
#######################################################################
PREPARE_HD <-function(train, dxgroup, atcgroup){
    train0 <- train %>%  dplyr::mutate(Y = ifelse(HHF_2yr_2yr==1, 1, 0), 
                                       W = ifelse(SGLT==1,1,0),     
                                       sex=as.numeric(sex),
                                       race=as.numeric(race))

  train00 <- train0 %>% 
            dplyr::select(BENE_ID,  IndexDate, Y, W, age, sex , race,
                          starts_with(c(paste0("dx",dxgroup) , 
                                        "cpt5", 
                                       paste0("atc", atcgroup)
                                        ),
                                      )

                  ) 
              # %>%  as.data.frame.matrix() #may need to delete this line
  #remove columns with only one level
  train00 <- train00[, sapply(train00, function(col) length(unique(col))) > 1] 
  
  
  return(train00)
}

#---------------------------------------------------------------------------------------------------------------
Train_BENEID_all <<- PREPARE_HD(Train2, dxgroup, atcgroup) 

if (atcgroup == 4) {
  columns_with_pattern <- grep("A10BJ|A10BK", names(Train_BENEID_all), value = TRUE) #atc4_outpt_once_A10BJ, atc4_outpt_once_A10BK
  print(columns_with_pattern)
  Train_BENEID_all <- Train_BENEID_all %>% select(-all_of(columns_with_pattern))
}
# Group by BENE_ID and count how many rows have the same ID
duplicate_counts <- Train_BENEID_all %>%
                    dplyr::group_by(BENE_ID) %>%
                    dplyr::summarise(count = n()) %>%
                    dplyr::filter(count > 1)
# View the number of duplicate rows
nrow(duplicate_counts)

table(Train_BENEID_all$Y,useNA = "always")
table(Train_BENEID_all$W,useNA = "always")


list_vars_with_missing_data <- function(df) {
  vars_with_missing <- sapply(df, function(x) anyNA(x))
  names(df)[vars_with_missing]
}


#---------------------------------------------------------------------------------------------------------------
vars_with_missing <- list_vars_with_missing_data(Train_BENEID_all)
print(vars_with_missing)

table(Train_BENEID_all$Y, useNA = "always") #for sensitivity analysis 7

Train <- Train_BENEID_all %>% select(-c("BENE_ID", "IndexDate"))
vars_with_missing <- list_vars_with_missing_data(Train)
print(vars_with_missing)


vars_forest = colnames( Train %>% dplyr::select(-c("Y", "W"))  ) 

XYW <- function(TrainDat){
  X <- TrainDat[,vars_forest]
  Y <- as.vector( as.numeric( TrainDat[,"Y"] ) )
  W <- as.vector( as.numeric( TrainDat[,"W"] ) )
  return(list(x=X, y=Y, w=W))
}
#all patients
X <<- XYW(Train)$x
Y <<- XYW(Train)$y
W <<- XYW(Train)$w

ncol(X); nrow(X); length(Y); length(W)
#Z<-Train[,vars_IV]
cf_raw_key.tr <- CF_RAW_key(Train, 1, "hd", hdPctTop=pct_inter) 
#==============================================#==============================================
Y.hat  <<- cf_raw_key.tr$Y.hat                 #
W.hat  <<- cf_raw_key.tr$W.hat  
HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw    # run for overall population or each subgroup
varimp_cf  <- cf_raw_key.tr$varimp_cf          #
#==============================================#==============================================

#W.hat    <- predict(grf::regression_forest(X, W))$predictions
summary(varimp_cf)
PSplot_allV <- GG_PS(Train, W.hat, "Propensity Score", "PS_allV")
pROC::roc(Y, Y.hat )

VI_lab_priortrim <-  PlotVI(varimp_cf, paste0(ncol(X), ' HD variables'), colnames(X))

VI_heat_priortrim <- GG_VI(varimp_cf#[which(varimp_cf > quantile(varimp_cf, pct_inter) )], 
                          ,
                          paste0( '', 
                                  #pct_inter*100, 
                                  paste0(#'% percentile (',
                                    ncol(X),
                                    #length(which(varimp_cf > quantile(varimp_cf, pct_inter) )),
                                    ' HD variables')
                                  ),
                          colnames( X#[which(varimp_cf > quantile(varimp_cf, pct_inter) )]
                          ) )

##############################################################
##############################################################
#   STEP 2.1 PS TRIM & HD feature preparation                #
##############################################################
##############################################################
library("plyr")
detach("package:plyr", unload = TRUE)

PS_trim_results  <- PS_trim( Train_BENEID_all, W.hat, "commonrange", NA)


Train= PS_trim_results[[1]]
ID_post_trim = PS_trim_results[[2]]

write.csv(ID_post_trim , paste0("c", c*100, "_n", n, "_I", dxgroup, "_A", atcgroup, "_L",Lcutoff,"_ID_posttrim.csv"), row.names=FALSE)

#---------------------------------------------------------------------------------------------

nrow(Train)
X <<- XYW(Train)$x
Y <<- XYW(Train)$y
W <<- XYW(Train)$w


cf_raw_key.tr <- CF_RAW_key(Train, 1, "hd", hdPctTop=pct_inter) 
#==============================================#==============================================
Y.hat  <<- cf_raw_key.tr$Y.hat                 #
W.hat  <<- cf_raw_key.tr$W.hat  
HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw    # run for overall population or each subgroup
varimp_cf  <- cf_raw_key.tr$varimp_cf          #
#==============================================#==============================================
PSplot_allV_trim_reesti <- GG_PS(Train, W.hat, "Propensity Score", "PS_allV_trim")


VI_lab_posttrim <- PlotVI(varimp_cf, paste0(ncol(X), ' HD variables'), colnames(X))

VI_heat_posttrim <- GG_VI(varimp_cf ,
                          paste0( '', 
                                  paste0(ncol(X), ' HD variables')
                                  ),
                          colnames( X ) )

####################################################################
# STEP 2.2 reassign extreme low frequent levels to the next levels #
####################################################################
#redefine X
# STEP 3: redefine dataset, selected_cf.idx and vars_catover2
X <<- FIX_LOW_FREQ (X, Lcutoff)


#redefine training set
Train <<- Train[,c("Y", "W", colnames(X)) ]

#ncol(Train)
vars_forest = colnames( Train %>% dplyr::select(-c("Y", "W" ))  ) #extract varaible names by index, #excluded IV


table(Train$Y, useNA = "always")
nrow(Train)
ncol(X)
################################################
################################################
# STEP 3: running iCF                         #
################################################
################################################


####################################################################################
######################
  cf_raw_key.tr <- CF_RAW_key(Train, 1, 
                              "hd", 
                               hdPctTop=pcttop) #use all selected variable in the 1st step
#==============================================#==============================================
Y.hat  <<- cf_raw_key.tr$Y.hat                 #
W.hat  <<- cf_raw_key.tr$W.hat  
HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw    # run for overall population or each subgroup
varimp_cf  <- cf_raw_key.tr$varimp_cf          #
#==============================================#==============================================

  #redefine selected X from shrinked dataset
  selected_cf.idx <<- cf_raw_key.tr$selected_cf.idx #MUST reselect important covaraites to run CF!!!!
  
  
  #if < 5 variables selected, then reselect top 5 variables
  if (  length(selected_cf.idx) <5) {
    cf_raw_key.tr <- CF_RAW_key(Train, 1, 
                                "hd", 
                                hdPctTop=5) #use all selected variable in the 1st step
    Y.hat  <<- cf_raw_key.tr$Y.hat
    W.hat  <<- cf_raw_key.tr$W.hat
    HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw
    HTE_P_cf.raw
    varimp_cf  <- cf_raw_key.tr$varimp_cf
    selected_cf.idx <<- cf_raw_key.tr$selected_cf.idx #MUST reselect important covaraites to run CF!!!!
    pcttop <<-5
  }
  

colnames(X[,c(selected_cf.idx)]) #sex not involved
  
write.csv(colnames(X[,c(selected_cf.idx)]) , "hdiCF_top10pct.csv", row.names=FALSE)


PSplot_allV_trim_reesti_fixL <- GG_PS(Train, W.hat, "Propensity Score", "PS_allV_trim_fixL")


length(selected_cf.idx)
time_rawCF <- cf_raw_key.tr$time_rawCF

VI_lab_posttrim_fixL <- PlotVI(varimp_cf, paste0(ncol(X), ' HD variables'), colnames(X))

VI_heat_posttrim_fixL <- GG_VI(varimp_cf , paste0( '', 
                                  paste0( ncol(X), ' HD variables') ),
                          colnames( X
                          ) )



cowplot::plot_grid( PSplot_allV, PSplot_allV_trim_reesti, PSplot_allV_trim_reesti_fixL, 
                    ncol  = 3, nrow=1,
                    labels = c("D)", "E)", "F)"), 
                    label_size = 15)

cowplot::plot_grid( VI_lab_priortrim, VI_lab_posttrim, VI_lab_posttrim_fixL, 
                    ncol  = 3, nrow=1,
                    labels = c("D)", "E)", "F)"),
                    label_size = 15) 

cowplot::plot_grid( VI_heat_priortrim, VI_heat_posttrim, VI_heat_posttrim_fixL, 
                    ncol  = 3, nrow=1,
                    labels = c("D)", "E)", "F)"),
                    label_size = 15) 


VI_heat_top5_sen3  <- GG_Xs(0.95) 
VI_heat_top5_sen4  <- GG_Xs(0.95) 
VI_heat_top5_sen5  <- GG_Xs(0.95) 
VI_heat_top5_sen6  <- GG_Xs(0.95) 
VI_heat_top5_sen7  <- GG_Xs(0.95) 
VI_heat_top5  <- GG_Xs(0.95) 
VI_heat_top1  <- GG_Xs(0.99) 
VI_heat_top10 <- GG_Xs(0.9) 


cowplot::plot_grid( VI_heat_top1, VI_heat_top5, VI_heat_top10, 
                    ncol  = 3, nrow=1,
                    labels = c("D)", "E)", "F)"),
                    label_size = 15) 

cowplot::plot_grid( VI_heat_top5_sen3, VI_heat_top5_sen4, VI_heat_top5_sen5, 
                    ncol  = 3, nrow=1,
                    labels = c("A)", "B)", "C)"), 
                    label_size = 15) 

cowplot::plot_grid( VI_heat_top5_sen6, VI_heat_top5_sen7,
                    ncol  = 2, nrow=1,
                    labels = c("A)", "B)"), 
                    label_size = 15) 

#Specify the decimal position for continuous variables in the subgroup definition.
split_val_round_posi=0
#Define categorical variables with more than two levels:
vars_catover2 <<- find_level_over2(X) 


P_threshold <<- 0.1
HTE_P_cf.raw <<- 0.1

#p1n200, if leaf size too large, can't form 2 leavies, i.e., one leaf only, will have the error:
#: parent_sign not found!!!
D2_MLS=MinLeafSizeTune(dat=Train, denominator=50, treeNo = 1000, iterationNo=50, split_val_round_posi=0, "D2", "firebrick")
D2_MLS$depth_gg

#1/29/2023 debug at /local/projects/medicare/DPP4i_HTE/programs/macros/iCF_SG_PIPELINE.R#19: iCF_D5 <- iCF(leafsize$D5, treeNo, iterationNo, Ntrain, "D5", split_val_round_posi)
D3_MLS=MinLeafSizeTune(dat=Train, denominator=70, treeNo = 1000, iterationNo=50, split_val_round_posi=0, "D3", "firebrick")
D3_MLS$depth_gg

#depth 4: 
D4_MLS=MinLeafSizeTune(dat=Train, denominator=95, treeNo = 1000, iterationNo=50, split_val_round_posi=0, "D4", "firebrick")
D4_MLS$depth_gg

#depth 5
D5_MLS=MinLeafSizeTune(dat=Train, denominator=115, treeNo = 1000, iterationNo=50, split_val_round_posi=0, "D5", "firebrick")
D5_MLS$depth_gg


leafsize <<- list(D5= mls_D5, #D5_MLS$denominator,
                  D4= mls_D4, #D4_MLS$denominator,
                  D3= mls_D3, #D3_MLS$denominator, 
                  D2= mls_D2 #D2_MLS$denominator
                  ) 

rm(list=ls(pattern="vote_D"))
rm(list=ls(pattern="stability_D"))

        hdiCF <- iCFCV(dat=Train,K=5,treeNo=1000,iterationNo=100,
                                                            min.split.var=4, split_val_round_posi=0, P_threshold=0.1, 
                                                            variable_type = "hd", hdPctTop= pcttop, HTE_P_cf.raw = 0.1)
saveRDS(hdiCF, 
file = paste0("RRcd_",analysis,"_c", c*100, "_n", n, 
              "_I", dxgroup, "_A", atcgroup, 
              "_K5_B1000_i100_Tc_L", Lcutoff, "_V", pcttop*100, ".rds") 
)
