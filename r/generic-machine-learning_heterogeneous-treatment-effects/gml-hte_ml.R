#------------------------------------------------------------------------------#
# This programme applies the Generic Machine Learning for Heterogeneous Treatment Effects
# methodology from Chernozhukhov et al. 2019. The code has been modified from previous code written by Demirer.
# Several modifications have been made to make the code more flexible, including the introduction of
# k-fold cross-validation and separating the machine learning code and analysis code into two files.
#------------------------------------------------------------------------------#
# Authors:  Paul Brimble - Mind and Behaviour Research Group,
#           Centre for the Study of African Economies (CSAE), University of Oxford
# Email:    paul.brimble@bsg.ox.ac.uk
#------------------------------------------------------------------------------#
# Preamble for original code:
# This program obtains empirical results for the paper "Generic Machine Learning Discovery
# and Classification Analysis of Heterogenous Treatment Effects in Randomized Experiments"
# by V. CHERNOZHUKOV, M. DEMIRER, E. DUFLO, I. FERNANDEZ-VAL
# Authors:  V. CHERNOZHUKOV, M. DEMIRER, E. DUFLO, I. FERNANDEZ-VAL
#------------------------------------------------------------------------------#

# The code is split into two files:
# 1) "gml-hte_ml" - machine learning estimation.
# 2) "gml-hte_an" - analysis using inputs from "gml-hte_ml".

# The "gml-hte_ml" code returns an RData file to be used in the "gml_hte_an" code

# The "gml_hte_an" code returns five outputs:
# 1) PLOTS: Two plots for each outcome variable reporting ATE and GATES (first plot is top 2 methods, second plot is all methods)
# 2) BLP: Latex table reporting ATE and heterogeneity loading factor (HET).
# 3) GATES: Latex table reporting the GATE for the most and least affected groups, and the difference.
# 4) CLAN: Latex table reporting the most affected and least affected averages for the CLAN variables and the difference.
# 5) BEST: Latex table reporting best machine learning methods.

# ML_Functions.R code is also required.

#------------------------------------------------------------------------------#
## A) INPUT DECISIONS
#------------------------------------------------------------------------------#
## There are 15 INPUT blocks to assist with customising this programme.

## A.1) LOAD PACKAGES ##########################################################
rm(list=ls(all=TRUE))
list_packages= c("foreign", "quantreg", "gbm", "glmnet",
           "MASS", "rpart", "doParallel", "sandwich", "randomForest",
           "nnet", "matrixStats", "xtable", "readstata13", "car", "lfe",
           "caret", "foreach", "multcomp","cowplot")
lapply(list_packages, require, character.only = TRUE)

time_start_ml <- proc.time()
print(noquote(paste0("START"," (date time: ", Sys.time(),"; elapsed time: ", (proc.time() - time_start_ml)[3],")")))
set.seed(1211);

## A.2) PRELIMINARIES ##########################################################

## Preliminaries [INPUT 1/15]
#setwd()                    # set working directory to save output
source("ML_Functions.R")    # source ML_Functions file

## Set Clusters [INPUT 2/15]
num_clusters <- detectCores(all.tests = FALSE, logical = TRUE) - 1           # number of clusters for parallel procesing
cl   <- makeCluster(num_clusters, outfile="")
registerDoParallel(cl)

## Name [INPUT 3/15]
name        <- "name"       # choose name

## A.3) DATA PREPARATION #######################################################

## Load Data [INPUT 4/15]
data        <- read.dta13("datafile.dta") # load data
# If multiple treatment arms, reduce to one:
data        <- data[data$treatment == 1 | data$control == 1,]

## Key Parameters [INPUT 5/15]
num_reps    <- 100  # number of repetitions
num_K       <- 2    # number of folds
tog_kfold   <- 0    # K-fold cross-validation (1) or use default split sampling (0)
tog_holarge <- 0    # If split sampling is selected, then holdout large (1) means
                    # that the holdout sample is larger than the training sample, otherwise set it to 0.
                    # To adjust split proportions, choose K > 2, and then toggle the holarge option accordinly.
tog_strat   <- 0    # stratified sampling (1) or random sampling (0)

##  Outcome Variables [INPUT 6/15]
names_outcomes  <- c("Var1_name", "Var2_name", "Var3_name")  # vector of labels for outcome variables
var_Y           <- c("var1", "var2", "var3")                 # vector of outcome variables

## Treatment Binary Variable [INPUT 7/15]
var_D           <- rep("treatment", length(var_Y))      # vector of treatment variables

## Cluster Variable [INPUT 8/15]
tog_cluster     <- 0    # set to 1 if there is a cluster variable for standard errors
var_cluster     <- "cluster"  # replace with cluster variable name

if(tog_cluster == 0){ var_cluster <- ""}

## Fixed Effects Variables [INPUT 9/15]
tog_fe1  <- 0        # set to 1 if there is a fixed effect, otherwise 0
tog_fe2  <- 0        # set to 1 if there is a second fixed effect, otherwise 0
var_fe1 <- "fixedeffect1"   # replace with fixed effect variable name
var_fe2 <- "fixedeffect2"   # replace with fixed effect variable name

## Fixed Effects Factorisation
if(tog_fe1 == 0){ var_fe1 <- "" }
if(tog_fe2 == 0){ var_fe2 <- "" }
if(tog_fe1 == 1){
    data[,var_fe1] <- factor(data[,var_fe1])
    a           <- as.data.frame(model.matrix(~data[,var_fe1]-1))
    colnames(a) <- paste0(var_fe1,"_",substring( names(a), nchar("data[,var_fe1]") + 2, 50) )
    data <- cbind(data, a)
    rm(a)
}
if(tog_fe2 == 1){
    data[,var_fe2] <- factor(data[,var_fe2])
    a           <- as.data.frame(model.matrix(~data[,var_fe2]-1))
    colnames(a) <- paste0(var_fe2,"_",substring( names(a), nchar("data[,var_fe2]") + 2, 50) )
    data <- cbind(data, a)
    rm(a)
}

## Covariates [INPUT 10/15]
var_covariates     <- c("cov1", "cov2", "cov3", "cov4", "cov5")    # covariates for estimation
tog_covariates_fe  <- 0         # set to 1 if covariates should include fixed effects, otherwise 0

if(tog_covariates_fe == 1) {
    var_covariates     <- c(var_covariates,                            # final covariates with fixed effects
                            names(data)[(substring( names(data), 1, nchar(var_fe1) + 1)  == paste0(var_fe1,"_") )],
                            names(data)[(substring( names(data), 1, nchar(var_fe2) + 1)  == paste0(var_fe2,"_") )] )
}

## Additional Control Variables [INPUT 11/15]
tog_var_add    <- 0                          # set to 1 if additional control variables for regressions
var_add        <- c("addvar1","addvar2")     # easier to make this list exhaustive (can reduce additional variables in analysis file)

## Additional Variable Form
form_var_add   <- ""
if(tog_var_add == 1){
    for(i in 1:length(var_add)){ form_var_add <- (paste0(form_var_add,"+",var_add[i])) }
}

# CLAN Variables [INPUT 12/15]
var_affected   <- c("cov1","cov2","cov3")    # easier to make this list exhaustive (can reduce affected variables in analysis file)
names_affected <- c("cov1_name","cov2_name","cov3_name")

## A.4) MACHINE LEARNING PREPARATION ###########################################

## svmPoly    : Support Vector Machines with Polynomial Kernel , package: kernlab, tuning parameters: degree (Polynomial Degree), scale (Scale), C (Cost)
## svmLinear  : Support Vector Machines with Linear Kernel , package: kernlab , C (Cost)
## svmLinear2 : Support Vector Machines with Linear Kernel , package: e1071 , cost (Cost)
## gbm        : Stochastic Gradient Boosting , package: gbm  , n.trees (# Boosting Iterations), interaction.depth (Max Tree Depth), shrinkage (Shrinkage), n.minobsinnode (Min. Terminal Node Size)
## glmnet     : Regularized Generalized Linear Models, package: glmnet , alpha (Mixing Percentage), lambda (Regularization Parameter)
## blackboost : Boosted Tree , package : mboost , mstop (#Trees), maxdepth (Max Tree Depth)
## nnet       : Neural Network , package : nnet  , size (#Hidden Units) , decay (Weight Decay)
## pcaNNet    : Neural Networks with Feature Extraction, package:nnet , size (#Hidden Units) , decay (Weight Decay)
## rpart      : CART, package:rpart, cp (Complexity Parameter)
## rf         : Random Forest,  package:randomForest, mtry (#Randomly Selected Predictors)

## Model names. For a list of available model names in caret package see: http://topepo.github.io/caret/available-models.html

## Machine Learning Methods [INPUT 13/15]
ml_methods       <- c("glmnet", "gbm", "pcaNNet", "rf")     # list of machine learning methods, please select at least 2.
names_methods    <- c("Elastic Net", "Boosting", "Nnet", "Random Forest")   # names for list

## A list of arguments for models used in the estimation
ml_args         <- list(svmLinear2=list(type='eps-regression'), svmLinear=list(type='nu-svr'), svmPoly=list(type='nu-svr'), gbm=list(verbose=FALSE), rf=list(ntree=1000), gamboost=list(baselearner='btree'), avNNet=list(verbose = 0, linout = TRUE, trace = FALSE), pcaNNet=list(linout = TRUE, trace = FALSE, MaxNWts=100000, maxit=10000), nnet=list(linout = TRUE, trace = FALSE, MaxNWts=100000, maxit=10000))

## Machine Learning Parameters [INPUT 14/15]
ml_resample    <- c("repeatedcv", "repeatedcv", "repeatedcv", "none")  # resampling method for chosing tuning parameters. available options: boot, boot632, cv, LOOCV, LGOCV, repeatedcv, oob, none
ml_tune       <- c(100, 20, 20, NA)                                    # number of elements per parameter in the grid. the grid size is tune^{number of tuning parameters}.

ml_proces     <- c("range", "range","range", "range")                  # pre-processing method
ml_select     <- c("best", "best","best", NA)                          # optimality criteria for choosing tuning parameter in cross validation. available options: best, oneSE, tolerance
ml_cv         <- c(2, 2,2, 2)                                          # the number of folds in cross-validation
ml_cvrep        <- c(2, 2,2, NA)                                       # number of iteration in repeated cross-validations

## Machine Learning Tuning Parameters [INPUT 15/15]

## If there is a parameter of the model that user doesn't want to choose with cross validation,
## it should be set using tune_param variable. Below mtry of random forest is set to 5
## for glmnet we want to choose both tuning parameters using cross validation so it is set to NULL

ml_tune_param       <- list(0)
ml_tune_param[[1]]  <- 0
ml_tune_param[[2]]  <- 0
ml_tune_param[[3]]  <- 0
ml_tune_param[[4]]  <- data.frame(mtry=5)

## Output Name
name_output <- paste("range","-","best", "-", 2, "-" ,2, "-",num_reps,sep="")

#------------------------------------------------------------------------------#
## B) ESTIMATION PROCEDURE
#------------------------------------------------------------------------------#

## B1) DATA FILTERING ##########################################################

## generate formula for x, xl is for linear models
var_X <- ""
for(i in 1:length(var_covariates)){ var_X <- paste(var_X, var_covariates[i], "+", sep = "") }
var_X  <- substr(var_X, 1, nchar(var_X)-1)
var_XL <- paste("(", var_X , ")", sep="")
var_XL <- var_X
rm(i)

## Data Preparation
data_temp     <- data[,c(var_Y, var_D, var_covariates)]
if(tog_var_add == 1){ data_temp <- cbind(data[,c(var_add)],     data_temp) ; names(data_temp)[1:length(var_add)] <- var_add }
if(tog_fe1 == 1){     data_temp <- cbind(data[,c(var_fe1)],     data_temp) ; names(data_temp)[1] <- var_fe1 }
if(tog_fe2 == 1){     data_temp <- cbind(data[,c(var_fe2)],     data_temp) ; names(data_temp)[1] <- var_fe2 }
if(tog_cluster == 1){ data_temp <- cbind(data[,c(var_cluster)], data_temp) ; names(data_temp)[1] <- var_cluster }
data <- data_temp
rm(data_temp)

## B2) MACHINE LEARNING ESTIMATION #############################################

## Start of Repetition Loop
results_ml <- foreach(t = 1:num_reps, .combine = 'cbind', .inorder = FALSE, .packages = list_packages) %dopar% {

    set.seed(t);

    ## B2.1) ESTIMATION PREPARATION ############################################

    ## Start of Randomisation (Create Splits and Obtain results_cv)
    if(tog_strat == 0){
        split             <- runif(nrow(data))
        results_cv        <- as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/num_K)),include.lowest = TRUE))
        rm(split)
    }
    if(tog_strat == 1){
        split_coni        <- which(data[,var_D[1]] == 1)
        split_con         <- runif(nrow(data[split_coni,]))
        cvgroup_con       <- as.numeric(cut(split_con,quantile(split_con,probs = seq(0, 1, 1/num_K)),include.lowest = TRUE))
        split_trei        <- which(data[,var_D[1]] == 0)
        split_tre         <- runif(nrow(data[split_trei,]))
        cvgroup_tre       <- as.numeric(cut(split_tre,quantile(split_tre,probs = seq(0, 1, 1/num_K)),include.lowest = TRUE))
        results_cv        <- matrix(NA,nrow(data),1)
        results_cv[split_coni,]  <- cvgroup_con
        results_cv[split_trei,]  <- cvgroup_tre
        rm(split_coni,split_con,cvgroup_con,split_trei,split_tre,cvgroup_tre)
    }

    ## Holdout
    if(tog_kfold == 0){ num_kfold <- 1 }

    ## K-Fold Cross-Validation
    if(tog_kfold == 1){
        num_kfold    <- num_K
        results_B    <- matrix(NA,nrow(data)*length(var_Y),length(ml_methods))
        results_S    <- matrix(NA,nrow(data)*length(var_Y),length(ml_methods))
        results_mdx  <- matrix(NA,nrow(data)*length(var_Y),length(ml_methods))
    }

    ## B2.2) CROSS-VALIDATION PREPARATION ######################################

    ## K-Fold Cross-Validation Loop
    for(k in 1:num_kfold){
        if(tog_kfold == 0){
            if(tog_holarge == 1){
                datause_raw       <- as.data.frame(data[results_cv == k,])
                dataout_raw       <- as.data.frame(data[results_cv != k,])
            }
            if(tog_holarge == 0){
                datause_raw       <- as.data.frame(data[results_cv != k,])
                dataout_raw       <- as.data.frame(data[results_cv == k,])
                }
        }
        if(tog_kfold == 1){
            datause_raw       <- as.data.frame(data[results_cv != k,])
            dataout_raw       <- as.data.frame(data[results_cv == k,])
        }

        ## Results to Store
        results_B1    <- matrix(NA,nrow(dataout_raw)*length(var_Y),length(ml_methods))
        results_S1    <- matrix(NA,nrow(dataout_raw)*length(var_Y),length(ml_methods))
        results_mdx1  <- matrix(NA,nrow(dataout_raw)*length(var_Y),length(ml_methods))

        ## B2.3) OUTCOME PREPARATION ###########################################

        ## Outcome Loop
        for(i in 1:length(var_Y)){
            y      <- var_Y[i]
            d      <- var_D[i]

            ## Clean Data
            datause   <- data.frame(datause_raw[complete.cases(datause_raw[, c(var_covariates, y, d )]),])
            dataout   <- data.frame(dataout_raw[complete.cases(dataout_raw[, c(var_covariates, y, d )]),])

            ## Treatment Indicator
            index_treat <- which(datause[,d] == 1)

            ## B2.4) MACHINE LEARNING PREPARATION ##############################

            ## Machine Learning Methods Loop
            for(l in 1:length(ml_methods)){

                if(ml_methods[l]=="glmnet"){ x <- var_XL }
                if(ml_methods[l]!="glmnet"){ x <- var_X  }
                if(ml_tune_param[[l]] == 0){ f = NULL }
                if(ml_tune_param[[l]] != 0){ f = ml_tune_param[[l]] }

                form_ml        <- as.formula(paste(y,"~",x,sep=""));

                ## B2.5) MACHINE LEARNING SCORE ESTIMATION #####################

                ml_mdx        <- rep((nrow(datause[datause[,d] == 1,]) + nrow(dataout[dataout[,d] == 1,]))/(nrow(datause) + nrow(dataout)), nrow(dataout))

                ml_fitControl <- trainControl(method = ml_resample[l], number = ml_cv[l], repeats = ml_cvrep[l], allowParallel = FALSE, verboseIter=FALSE, search="random", selectionFunction=ml_select[l])
                ml_arg        <- c(list(form = form_ml, data = datause[index_treat,],  method = ml_methods[l],  tuneGrid = f, trControl = ml_fitControl, preProcess=ml_proces[l], tuneLength=ml_tune[l]), ml_args[[ml_methods[l]]])
                ml_fit.yz1    <- suppressWarnings(do.call(caret::train, ml_arg))
                ml_z1x        <- predict(ml_fit.yz1, newdata = dataout, type = "raw")

                ml_fitControl <- trainControl(method = ml_resample[l], number = ml_cv[l], repeats = ml_cvrep[l], allowParallel = FALSE, verboseIter=FALSE, search="random", selectionFunction=ml_select[l])
                ml_arg        <- c(list(form = form_ml, data = datause[-index_treat,],  method = ml_methods[l], tuneGrid = f, trControl = ml_fitControl, preProcess=ml_proces[l], tuneLength=ml_tune[l]), ml_args[[ml_methods[l]]])
                ml_fit.yz0    <- suppressWarnings(do.call(caret::train, ml_arg))
                ml_z0x        <- predict(ml_fit.yz0, newdata = dataout, type = "raw")

                index_mdx     <- (ml_mdx>0.01 & ml_mdx<0.99)
                dataout       <- dataout[index_mdx, ]
                ml_z1x        <- ml_z1x[index_mdx]
                ml_z0x        <- ml_z0x[index_mdx]
                ml_mdx        <- ml_mdx[index_mdx]

                ## Basline Conditional Average (BCA)
                results_B1[  (nrow(dataout_raw)*(i-1)+1):(nrow(dataout_raw)*(i-1)+nrow(dataout)),l] <- ml_z0x

                ## Conditional Average Treatment Effect (CATE)
                results_S1[  (nrow(dataout_raw)*(i-1)+1):(nrow(dataout_raw)*(i-1)+nrow(dataout)),l] <- ml_z1x - ml_z0x

                ## Treatment Propensity
                results_mdx1[(nrow(dataout_raw)*(i-1)+1):(nrow(dataout_raw)*(i-1)+nrow(dataout)),l] <- ml_mdx

            }   ## End of Machine Learning Method Loop
        }   ## End of Outcome Loop

        ## B2.6) STORE RESULTS #################################################

        if(tog_kfold == 0){
            results_B   <- results_B1
            results_S   <- results_S1
            results_mdx <- results_mdx1
        }
        if(tog_kfold == 1){
            results_B[results_cv == k]   <- results_B1
            results_S[results_cv == k]   <- results_S1
            results_mdx[results_cv == k] <- results_mdx1
        }
        print(noquote(paste0(t," (",k,"/",num_K,")")))

    }   ## End of K-Fold Cross-Validation Loop

    result_all_vector <- c(as.vector(results_cv), as.vector(results_B), as.vector(results_S), as.vector(results_mdx))
    print(noquote(paste0(t," (date time: ", Sys.time(),"; elapsed time: ", (proc.time() - time_start_ml)[3],")")))
    results_ml <- data.frame(result_all_vector)

}   ## End of Repetition Loop

## B3) DATA SAVING #############################################################

## Data Timing
time_end_ml <- proc.time() - time_start_ml
print(noquote(paste0("END"," (date time: ", Sys.time(),"; elapsed time: ", (time_end_ml)[3],")")))

## Save File
save.image(file=paste0(name, "_ml.RData"))

## Stop Cluster
stopCluster(cl)
