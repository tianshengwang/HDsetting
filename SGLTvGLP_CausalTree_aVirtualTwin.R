rm(list=ls())
setwd("/local/projects/medicare/sglt_cvs/data/full")
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
# Data preparation                          
#
#####################################################
#####################################################
#-----------------------#
analysis ="primary"        #
#-----------------------#
drug1="sglt"; drug2="glp"; pct_inter=0

if (analysis=="primary"){
  c=0.01; n=200; dxgroup <<- 3; atcgroup <<- 3; Lcutoff=20; pcttop=0.95; 
  
} else if (analysis=="sen1"){
  c=0.01; n=200; dxgroup <<- 3; atcgroup <<- 3; Lcutoff=20; pcttop=0.99; 
  
} else if (analysis=="sen2"){
  c=0.01; n=200; dxgroup <<- 3; atcgroup <<- 3; Lcutoff=35; pcttop=0.9;
  
} else if (analysis=="sen3"){
  c=0.01; n=100; dxgroup <<- 3; atcgroup <<- 3; Lcutoff=25; pcttop=0.95; 
  
} else if (analysis=="sen4"){
  c=0.02; n=200; dxgroup <<- 3; atcgroup <<- 3; Lcutoff=20; pcttop=0.95;
  
} else if (analysis=="sen5"){
  c=0.05; n=200; dxgroup <<- 3; atcgroup <<- 3 ; Lcutoff=20; pcttop=0.95;
  
} else if (analysis=="sen6"){
  c=0.01; n=200; dxgroup <<- 4; atcgroup <<- 3 ; Lcutoff=30; pcttop=0.95;
  
} else if (analysis=="sen7"){
  c=0.01; n=200; dxgroup <<- 3; atcgroup <<- 4 ; Lcutoff=25; pcttop=0.95;
  
} else if (analysis=="sen8"){
  c=0.02; n=200; dxgroup <<- 3; atcgroup <<- 4 ; Lcutoff=20; pcttop=0.95;
}


Train_0 = haven::read_sas(paste0("hdicf_", drug1, "v", drug2, "_p",c*100, "_n",n, "_i",dxgroup, "_a", atcgroup,".sas7bdat"))

Train_0_date <- Train_0 %>% 
  dplyr::filter(FillDate2 <= as.Date("2019-12-31") & #21914 &  #31DEC2019
                IndexDate >= (as.Date("2015-10-15") + 365) & #(20376 + 365) &#20376 & #15OCT2015 + 365 daysBL 
                FillDate2 <= (as.Date("2019-12-31")-365*2) & #(21914-365*2) & #& #31DEC2019 due to limited sample size
                is.na(FillDate2)==F &
                excludeFlag_preFill2Initiator ==0 &
                excludeFlag_sameDayInitiator==0 &
                excludeFlag_prevalentUser==0 
            ) 

is.na(Train_0_date$HFPRIMARY_ICD10DX_date)
Train1 <- Train_0_date %>% 
  transform(min_CensorF2_730 = pmin(censorDate_ITT, FillDate2+ 365*2),
            min_CensorF2_730_HFF = pmin(censorDate_ITT, HFPRIMARY_ICD10DX_date, FillDate2 + 365*2))%>%
  dplyr::mutate(
    HHFfu_2yr = ifelse(is.na(HFPRIMARY_ICD10DX_date)==T,
                       min_CensorF2_730 - FillDate2,
                       min_CensorF2_730_HFF - FillDate2),
    HHF   = ifelse( is.na( HFPRIMARY_ICD10DX_date) != T & #has event date
                      HFPRIMARY_ICD10DX_date >= FillDate2, 
                    1, 0),
    
    age = as.numeric( cut(age, c(65,70,75,80,85,Inf) ,
                          labels=c("65<age<=70 ","70<age<=75","75<age<=80","80<age<=85", "age>85")
    )) 
    
  ) %>%
  dplyr::mutate(race = case_when(race == "1" ~ 1,
                                 race == "2" ~ 2,
                                 !(race %in% c("1", "2")) ~  3)
  )

#rm(Train_0)
nrow(Train1) #15388
table(Train1$HHF)


is.na(Train1$death_dt)
Train2 <- Train1 %>%
  dplyr::mutate(HHF_3yr_3yr = ifelse(HHF==1 &
                                       HFPRIMARY_ICD10DX_date <=censorDate_ITT &
                                       HFPRIMARY_ICD10DX_date <= IndexDate + 365*3 ,
                                     1, 0),
                HHF_2yr_2yr =  ifelse(  HHF==1 & # event date is later than the start of follow up 
                                          
                                          (#not dead by 12/31/2019
                                            (is.na(death_dt) == T & 
                                               #HFPRIMARY_ICD10DX_date <=  min(FillDate2+730, censorDate_ITT) this min Fx doesn't work!need to figure out why
                                               HFPRIMARY_ICD10DX_date <=  FillDate2+730 &
                                               HFPRIMARY_ICD10DX_date <=  censorDate_ITT    ) |
                                              #dead by 12/31/2019
                                              (is.na(death_dt) != T & 
                                                 #HFPRIMARY_ICD10DX_date <=  min(FillDate2+730, censorDate_ITT, death_dt) this min Fx doesn't work!need to figure out why
                                                 HFPRIMARY_ICD10DX_date <=  FillDate2+730 &
                                                 HFPRIMARY_ICD10DX_date <=  censorDate_ITT & 
                                                 HFPRIMARY_ICD10DX_date <=  death_dt
                                              ) 
                                          )
                                        , 
                                        1, 
                                        0) ,                                                       
                HHF_1yr_1yr = ifelse(HHF==1 &
                                       HFPRIMARY_ICD10DX_date <=censorDate_ITT &
                                       HFPRIMARY_ICD10DX_date <= IndexDate + 365*1 ,
                                     1, 0),
                HHF_05yr_05yr = ifelse(HHF==1 &
                                         HFPRIMARY_ICD10DX_date <=censorDate_ITT &
                                         HFPRIMARY_ICD10DX_date <= IndexDate + 365*0.5 ,
                                       1, 0),
                
  )

nrow(Train2) #15388
table(Train2$HHF)
table(Train2$HHF_2yr_2yr)

#######################################################################
#  STEP 1: READ DATA from Medicare or MarketScan: Fx to get CF-ready  #
#######################################################################
PREPARE_HD <-function(train, dxgroup, atcgroup){  train0 <- train %>%  
                                      dplyr::mutate(Y = ifelse(HHF_2yr_2yr==1, 1, 0), 
                                     W = ifelse(SGLT==1,1,0),     
                                     sex=as.numeric(sex),
                                     race=as.numeric(race))
  
  train00 <- train0 %>% 
    
    dplyr::select( BENE_ID,  IndexDate, Y, W, age, sex , race,
      starts_with(c(paste0("dx",dxgroup) , 
                    "cpt5", 
                    paste0("atc", atcgroup)
      ),
      
      )
      
    ) %>%  as.data.frame.matrix() 
  train00 <- train00[, sapply(train00, function(col) length(unique(col))) > 1] 
  
  
  return(train00)
}

#---------------------------------------------------------------------------------------------------------------
Train_BENEID_all <<- PREPARE_HD(Train2, dxgroup, atcgroup) 

#ATC4th level will specify treatment (SGLT2i) and comparator (GLP1RA) as baseline period include indexdate
if(atcgroup==4){
  Train_BENEID_all <- Train_BENEID_all %>% select(-atc4_outpt_A10BK, -atc4_outpt_A10BJ)
}

table(Train_BENEID_all$atc4_outpt_A10BK, useNA = "always")
table(Train_BENEID_all$atc4_outpt_A10BJ, useNA = "always")


list_vars_with_missing_data <- function(df) {
  vars_with_missing <- sapply(df, function(x) anyNA(x))
  names(df)[vars_with_missing]
}

#---------------------------------------------------------------------------------------------------------------

vars_with_missing <- list_vars_with_missing_data(Train_BENEID_all)
print(vars_with_missing)

table(Train_BENEID_all$Y, useNA = "always")
Train_BENEID_all$Y[is.na(Train_BENEID_all$Y)] <- 0
table(Train_BENEID_all$Y, useNA = "always") #for sensitivity analysis 7


Train <- Train_BENEID_all %>% select(-c("BENE_ID", "IndexDate"))
vars_with_missing <- list_vars_with_missing_data(Train)
print(vars_with_missing)



vars_forest = colnames( Train %>% dplyr::select(-c("Y", "W" ))  ) 

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


# PS matched cohort for non-CF methods
#------------------------------------------
#traditional logistic regression for PS
#------------------------------------------
psform1 <<- as.formula(paste0("W ~ ", paste0(colnames(X), collapse = " + ")))

(aGlm_tr <- glm(psform1, family=binomial(), data=Train))
ps_tr <- predict(object = aGlm_tr, type = "response")

ps_m_tr <- optmatch::match_on(aGlm_tr, caliper = 0.25*sd(ps_tr))

fm_tr   <- optmatch::fullmatch(ps_m_tr, data = Train)
Train_m <- na.omit(cbind(Train, fm_tr)) %>% select(-c(fm_tr))



#####################################################
#####################################################
#
#  aVirtualTwin ANALYSIS: only works for binary outcome                          
#
#####################################################
#####################################################
library("aVirtualTwins")
split_val_round_posi=0

s_VT =  Sys.time()

if (length(unique(Y)) ==2){
  vt.o <- vt.data(Train_m %>% dplyr::mutate(Y= ifelse(Y==1,0,1)), 
                  "Y", 
                  "W", 
                  interactions = TRUE)
  set.seed(123)
  model.rf <- randomForest::randomForest(x = vt.o$getX(interactions = T),
                                         y = vt.o$getY(),
                                         ntree = 1000)
  vt.f.rf <- vt.forest("one", vt.data = vt.o, 
                       model = model.rf, 
                       interactions = T)
  # grow RF for T = 1
  model.rf.trt1 <- randomForest(x = vt.o$getX(trt = 1), 
                                y = vt.o$getY(trt = 1))
  # grow RF for T = 0
  model.rf.trt0 <- randomForest(x = vt.o$getX(trt = 0), y = vt.o$getY(trt = 0))
  # initialize VT.forest.double()
  vt.doublef.rf <- vt.forest("double",
                             vt.data = vt.o, 
                             model_trt1 = model.rf.trt1, 
                             model_trt0 = model.rf.trt0)
  #model.fold <- vt.forest("fold", vt.data = vt.o, fold = 5, ratio = 1, interactions = T, ntree = 200)
  
  
  # initialize classification tree
  tr.class <- vt.tree("class",
                      #vt.difft = vt.f.rf,
                      vt.difft = vt.doublef.rf,
                      sens = ">",
                      threshold = quantile(vt.f.rf$difft, seq(.5, .8, .1)),
                      maxdepth = 3,
                      cp = 0,
                      maxcompete = 2) 
  # tr.class is a list if threshold is a vectoor
  class(tr.class)
  class(tr.class$tree1)
  
  tr.reg <- vt.tree("reg",
                    # vt.difft = vt.f.rf,
                    vt.difft = vt.doublef.rf,
                    sens = ">",
                    threshold = quantile(vt.f.rf$difft, seq(.5, .8, .1)))
  # tr.class is a list if threshold is a vectoor
  class(tr.reg)
  class(tr.reg$tree1)
  
  vt.sbgrps <- vt.subgroups(tr.class)
  
  # print tables with knitr package, subgroup decision from aVirtualTwin!!!
  knitr::kable(vt.sbgrps)
  
  Deci_VT0 <- vt.sbgrps %>%
    dplyr::mutate(Subgroup_size = as.numeric(vt.sbgrps$`Subgroup size`)) %>%
    dplyr::filter(Subgroup_size > 500 ) %>% #require sample size > 50
    `rownames<-`( NULL ) %>% #remove original long rowname from VirtualTwin
    dplyr::select(Subgroup) %>%
    `colnames<-`( NULL ) 
  #modity the way to present binary (0,1) split value
  #remove "&" to prepare for sorting
  #dplyr::mutate(subgroup2 = stringr::str_replace_all(Subgroup, " & ", " "))  %>%
  #sort conditions in each subgroup definition
  #convert df to list 
  Deci_VT_L <-  split(Deci_VT0, seq(nrow(Deci_VT0)))
  #apply function that convert VirtualTwin format to Causal forest format
  Deci_VT_cf.format <- lapply(Deci_VT_L, VT2CF_format )
  
  Deci_VT_con_split2 <- data.frame(matrix(unlist(Deci_VT_cf.format), nrow=length(Deci_VT_cf.format), byrow=TRUE)) %>%
    `colnames<-`("Subgroup") %>% 
    #sort conditions in each subgroup definition
    rowwise() %>% 
    mutate(subgroup = paste(sort(unlist(strsplit(as.character(Subgroup), " & ", fixed = TRUE)), decreasing = TRUE), collapse = " & ")) 
  #sort subgroup definitions
  Deci_VT <-  Deci_VT_con_split2[order(Deci_VT_con_split2$subgroup),] %>%
    tibble::rowid_to_column() %>%
    dplyr::select(subgroup, rowid)%>%
    dplyr::rename(subgroupID =rowid) %>%
    dplyr::select(subgroupID, everything())
  
  
} else if (length(unique(Y)) >8){
  Deci_VT = "NA" 
} 
SG2INT_VT <- SG2INT(Deci_VT)

e_VT =  Sys.time()
time_VT =  e_VT - s_VT

saveRDS(Deci_VT, 
        file = paste0("hdPS_c", c*100, "_n", n, 
                      "_I", dxgroup, "_A", atcgroup, 
                      "VT", ".rds") 
)


#####################################################
#####################################################
#
#  causal tree                          
#
#####################################################
#####################################################
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/causalForest.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/causalTee.control.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/causalTree.anova.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/causalTree.branch.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/causalTree.matrix.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/causalTree.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/causalTreecallback.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/causalTreeco.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/est.causalTree.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/estimate.causalTree.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/formatg.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/honest.causalTree.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/honest.est.rparttree.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/honest.rparttree.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/importance.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/model.frame.causalTree.R")
source("/local/projects/medicare/DPP4i_HTE/programs/macros/causalTree/na.causalTree.R")


#library(causalTree) installation on R.4.4.1 doesn't work
yform1 <<- as.formula(paste0("Y ~ ", paste0(colnames(X), collapse = " + ")))
s_CT =  Sys.time()
tree <- causalTree(yform1, 
                   data = Train_m, 
                   treatment = Train_m$W,
                   split.Rule = "CT", 
                   cv.option = "CT", 
                   split.Honest = T, 
                   cv.Honest = T, 
                   split.Bucket = F, 
                   xval = 5, #CV fold
                   cp = 0, #complexity parameter
                   minsize = 500, #min.leaf.size
                   propensity = 0.5)

opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]

opfit <- prune(tree, opcp)

rpart.plot(opfit)

saveRDS(tree, 
        file = paste0("hdPS_c", c*100, "_n", n, 
                      "_I", dxgroup, "_A", atcgroup, 
                      "CT", ".rds") 
)

e_CT =  Sys.time()
time_CT =  e_CT - s_CT
