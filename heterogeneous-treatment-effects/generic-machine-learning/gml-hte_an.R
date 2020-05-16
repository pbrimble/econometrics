#--------------------------------------------------------------------------------------------------------------------
# This program obtains empirical results for the paper "Generic Machine Learning Discovery
# and Classification Analysis of Heterogenous Treatment Effects in Randomized Experiments"
# by V. CHERNOZHUKOV, M. DEMIRER, E. DUFLO, I. FERNANDEZ-VAL
#--------------------------------------------------------------------------------------------------------------------
# Authors:  V. CHERNOZHUKOV, M. DEMIRER, E. DUFLO, I. FERNANDEZ-VAL
#--------------------------------------------------------------------------------------------------------------------
# Several modifications have been made to make the code more flexible, including the introduction of
# k-fold cross-validation and separating the machine learning code and analysis code into two files.
#--------------------------------------------------------------------------------------------------------------------
# Authors:  Paul Brimble - Mind and Behaviour Research Group,
#           Centre for the Study of African Economies (CSAE), University of Oxford
# Email:    paul.brimble@bsg.ox.ac.uk

# The code is split into two files:
# 1) "ml" - machine learning estimation.
# 2) "an" - analysis using inputs from ml.

# The "ml" code returns an RData file to be used in the "an" code

# The "an" code returns five outputs
# 1) Two plots for each outcome variable reporting ATE by quantile groups (first plot is top 2 ML methods, second plot is all ML methods)
# 2) het_results: Latex table reporting ATE and heterogeneity loading coefficient
# 3) test_results: A latex table reporting estimated ATE_group1 - ATE_group5 and its p-value
# 4) index_ml_best: A latex table reporting index_ml_best ML methods
# 5) RData file.

# ML_Functions.R code is also required.

####################################### Input Decisions  #######################################
# There are 4 INPUT blocks to assist with customising this programme.

####################################### Load Packages #######################################
rm(list=ls(all=TRUE))

list_packages= c("foreign", "quantreg", "gbm", "glmnet",
           "MASS", "rpart", "doParallel", "sandwich", "randomForest",
           "nnet", "matrixStats", "xtable", "readstata13", "car", "lfe",
           "caret", "foreach", "multcomp","cowplot")
lapply(list_packages, require, character.only = TRUE)
list_packages_an = c("lfe", "multcomp","matrixStats","xtable")

time_start_an <- proc.time()

set.seed(1211);

####################################### Load and Process Data  #######################################
# Name [INPUT 1/4]
name        <- "name"

# Load Data
load(paste0(name, "_ml.RData"))

####################################### Parameters  #######################################
# Key Parameters [INPUT 2/4]
num_groups     <- 5        # number of quantile groups (3, 4 or 5)
num_thres      <- 1/num_groups     # quantile for most/least var_affected group
num_alpha      <- 0.05    # signifigance level
tog_mono       <- 0      # rearrange for GATES monotonicity ("1"), or no rearrangement ("0")
num_kfold_spec <- 0      # if kfold=="1", choose specific fold to run analysis on to return to split sampling method
                     # or run on full sample ("0")
tog_endostrat  <- 0     # choose "1" for endogenous stratification, otherwise will sort by predicted treatment effect ("0")

tog_silence    <- 0

# Additional Control Variables (Change if Necessary) [INPUT 3/4]
tog_var_add_change  <- 0                         # if no additional variables <- "0"
if(tog_var_add_change == 1){
    var_add        <- c("addvar1","addvar2")     # easier to make this list exhaustive (can reduce additional variables in analysis file)
    form_var_add   <- ""
    for(i in 1:length(var_add)){ form_var_add <- (paste0(form_var_add,"+",var_add[i])) }
}

# CLANs (Change if Necessary) [INPUT 4/4]
tog_var_affected_change <- 0
if(tog_var_affected_change == 1){
    var_affected       <- c("cov1","cov2")
    names_affected <- c("cov1_name","cov2_name")
}

# Group Variable Add Form
form_var_group = ""
for(g in 1:num_groups){
    form_var_group <- paste0(form_var_group,"+G",g)
}
####################################### Estimation  #######################################

results_an <- foreach(t = 1:num_reps, .combine='cbind', .inorder=FALSE, .packages=list_packages_an)  %dopar% {

    set.seed(t);

    results_blp_ate      <- matrix(NA,5*length(var_Y), length(ml_methods))
    results_blp_het      <- matrix(NA,5*length(var_Y), length(ml_methods))
    results_gates_tests  <- matrix(NA,15*length(var_Y), length(ml_methods))
    results_gates        <- matrix(NA,3*num_groups*length(var_Y), length(ml_methods))
    results_clan        <- matrix(NA, (length(var_affected)*15)*length(var_Y), length(ml_methods))
    results_ml_best      <- matrix(NA, 2*length(var_Y), length(ml_methods))

  ## Call ML Results for All Methods
    if(tog_kfold == 0){
        results_cv        <- array(c(as.vector(results_ml[1:nrow(data),t])))
        if(tog_holarge == 1){ dataout_raw <- as.data.frame(data[results_cv != 1,]) }
        if(tog_holarge == 0){ dataout_raw <- as.data.frame(data[results_cv == 1,]) }
        B_all             <- array(c(as.matrix(results_ml[(nrow(data)+1):(nrow(data)+nrow(dataout_raw)*length(var_Y)*length(ml_methods)),t])),c(nrow(dataout_raw)*length(var_Y),length(ml_methods)), )
        S_all             <- array(c(as.matrix(results_ml[(nrow(data)+1+nrow(dataout_raw)*length(var_Y)*length(ml_methods)):(nrow(data)+2*nrow(dataout_raw)*length(var_Y)*length(ml_methods)),t])),c(nrow(dataout_raw)*length(var_Y),length(ml_methods)), )
        mdx_all           <- array(c(as.matrix(results_ml[(nrow(data)+1+2*nrow(dataout_raw)*length(var_Y)*length(ml_methods)):(nrow(data)+3*nrow(dataout_raw)*length(var_Y)*length(ml_methods)),t])),c(nrow(dataout_raw)*length(var_Y),length(ml_methods)), )
    }
    if(tog_kfold == 1){
        B_all             <- array(c(as.matrix(results_ml[(nrow(data)+1):(nrow(data)+nrow(data)*length(var_Y)*length(ml_methods)),t])),c(nrow(data)*length(var_Y),length(ml_methods)), )
        S_all             <- array(c(as.matrix(results_ml[(nrow(data)+1+nrow(data)*length(var_Y)*length(ml_methods)):(nrow(data)+2*nrow(data)*length(var_Y)*length(ml_methods)),t])),c(nrow(data)*length(var_Y),length(ml_methods)), )
        mdx_all           <- array(c(as.matrix(results_ml[(nrow(data)+1+2*nrow(data)*length(var_Y)*length(ml_methods)):(nrow(data)+3*nrow(data)*length(var_Y)*length(ml_methods)),t])),c(nrow(data)*length(var_Y),length(ml_methods)), )
        if(num_kfold_spec>0){
            results_cv    <- array(c(as.vector(results_ml[1:nrow(data),t])))
            dataout_raw   <- as.data.frame(data[results_cv == num_kfold_spec,])
            B_all         <- B_all[  results_cv == num_kfold_spec,]
            S_all         <- S_all[  results_cv == num_kfold_spec,]
            mdx_all       <- mdx_all[results_cv == num_kfold_spec,]
        }
    }

    # Prepare Data
    for(i in 1:length(var_Y)){
        y      <- var_Y[i]
        d      <- var_D[i]

        if((tog_kfold == 0) | (tog_kfold == 1 & num_kfold_spec > 0)){
            dataest   <- data.frame(dataout_raw[complete.cases(dataout_raw[, c(var_covariates, y, d)]),])
            mdx_raw   <- mdx_all[(nrow(dataout_raw)*(i-1)+1):(nrow(dataout_raw)*i),]
            B_raw     <- B_all[(  nrow(dataout_raw)*(i-1)+1):(nrow(dataout_raw)*i),]
            S_raw     <- S_all[(  nrow(dataout_raw)*(i-1)+1):(nrow(dataout_raw)*i),]
            mdx_est   <- mdx_raw[complete.cases(mdx_raw),]
            B_est     <- B_raw[  complete.cases(B_raw),]
            S_est     <- S_raw[  complete.cases(S_raw),]
        }
        if(tog_kfold == 1 & num_kfold_spec == 0){
            dataest   <- data.frame(data[complete.cases(data[, c(var_covariates, y, d)]),])
            mdx_raw   <- mdx_all[(nrow(data)*(i-1)+1):(nrow(data)*i),]
            B_raw     <- B_all[(  nrow(data)*(i-1)+1):(nrow(data)*i),]
            S_raw     <- S_all[(  nrow(data)*(i-1)+1):(nrow(data)*i),]
            mdx_est   <- mdx_raw[complete.cases(mdx_raw),]
            B_est     <- B_raw[  complete.cases(B_raw),]
            S_est     <- S_raw[  complete.cases(S_raw),]
        }

        for(l in 1:length(ml_methods)){
            ############ Load Scores from ML ############
            md_x         <- mdx_est[,l]
            B            <- B_est[,l]
            S            <- S_est[,l]
            ############################################# GATES #############################################

            if(tog_endostrat == 0){ S2 <- S+runif(length(S), 0, 0.00001) }
            if(tog_endostrat == 1){ S2 <- B+runif(length(B), 0, 0.00001) }

            breaks    <- quantile(S2, seq(0,1, num_thres),  include.lowest =T)
            breaks[1] <- breaks[1] - 0.001
            breaks[num_groups+1] <- breaks[num_groups+1] + 0.001
            SG        <- cut(S2, breaks = breaks)
            SGX       <- model.matrix(~-1+SG)
            DSG       <- data.frame(as.numeric(I(as.numeric(dataest[,d])-md_x))*SGX)

            # Add Variables to Dataest
            dataest[, c("B", "S")] <- cbind(B, S)
            for(g in 1:num_groups){ dataest[, c(paste0("G",g))] <- cbind(DSG[,g]) }
            dataest[, c("weight")] <- cbind(as.numeric((1/(md_x * (1 - md_x)))))

            if(var(dataest$B) == 0){ dataest$B <- dataest$B + rnorm(length(dataest$B),  mean=0, sd=0.1) }
            if(var(dataest$S) == 0){ dataest$S <- dataest$S + rnorm(length(dataest$S),  mean=0, sd=0.1) }

            form1 <- as.formula(paste(y, "~", "B+S", form_var_group, form_var_add, "|",  var_fe1, "+", var_fe2, "| 0 |", var_cluster, sep=""))

            a <- tryCatch({
            a <- felm(form1, data=dataest, weights=dataest$weight)
            },error=function(e){
                if(tog_silence == 0){cat("ERROR :",ml_methods[l], t, i, "\n")}
                form1  <- as.formula(paste(y, "~", form_var_group, form_var_add, "|", var_fe1, "+", var_fe2, "| 0 |", var_cluster, sep=""))
                reg    <- felm(form1, data=dataest, weights=dataest$weight)
                return(reg)
            }, warning = function(war) {
                if(tog_silence == 0){cat("WARNING :",ml_methods[l], t, i, "\n")}
                form1  <- as.formula(paste(y, "~", form_var_group, " | ", var_fe1, "+", var_fe2, "| 0 |", var_cluster, sep=""))

                reg    <- felm(form1, data=dataest, weights=dataest$weight)
                return(reg)
            })
            reg   <- a

            crit <- qnorm(1-num_alpha/2)

            mean <- numeric(0)
            sd   <- numeric(0)
            for(g in 1:num_groups){
                mean[g] <- summary(reg)$coef[c(paste0("G",g)),1]
                sd[g]   <- summary(reg)$coef[c(paste0("G",g)),2]
            }

            if(tog_mono == 1){
                results_gates[((i-1)*3*num_groups+1):((i-1)*3*num_groups+num_groups),l]        <- sort(mean)
                results_gates[((i-1)*3*num_groups+num_groups+1):((i-1)*3*num_groups+2*num_groups),l]   <- sort(mean +crit*sd)
                results_gates[((i-1)*3*num_groups+2*num_groups+1):((i-1)*3*num_groups+3*num_groups),l] <- sort(mean -crit*sd)

                Gmax <- paste("G",toString(which.max(mean)),sep="")
                Gmin <- paste("G",toString(which.min(mean)),sep="")
            }
            if(tog_mono == 0){
                results_gates[((i-1)*3*num_groups+1):((i-1)*3*num_groups+num_groups),l]        <- mean
                results_gates[((i-1)*3*num_groups+num_groups+1):((i-1)*3*num_groups+2*num_groups),l]   <- mean +crit*sd
                results_gates[((i-1)*3*num_groups+2*num_groups+1):((i-1)*3*num_groups+3*num_groups),l] <- mean -crit*sd

                Gmax <- paste("G",num_groups,sep="")
                Gmin <- "G1"
            }

            coef <- (summary(reg)$coefficients[Gmax,1])
            pval <- (summary(reg)$coefficients[Gmax,4])
            results_gates_tests[(1 + (i - 1) * 15):(5 + ((i - 1) * 15)), l] <-
                    c(coef, confint(reg, Gmax, level = 1 - num_alpha)[1:2],
                    (as.numeric(coef < 0) * (pval/2) + as.numeric(coef > 0) * (1 - pval/2)),
                    (as.numeric(coef > 0) * (pval/2) + as.numeric(coef < 0) * (1 - pval/2)) )

            coef <- (summary(reg)$coefficients[Gmin,1])
            pval <- (summary(reg)$coefficients[Gmin,4])
            results_gates_tests[(6+(i-1)*15):(10+((i-1)*15)),l]  <-
                    c(coef, confint(reg, Gmin, level = 1 - num_alpha)[1:2],
                    (as.numeric(coef < 0) * (pval/2) + as.numeric(coef > 0) * (1 - pval/2)),
                    (as.numeric(coef > 0) * (pval/2) + as.numeric(coef < 0) * (1 - pval/2)) )

            Gdif <- paste(Gmax,"-",Gmin," == 0",sep="")
            test <- glht(reg, linfct = c(Gdif))
            coef <- (summary(reg)$coefficients[Gmax,1]) - (summary(reg)$coefficients[Gmin,1])
            pval <- summary(test)$test$pvalues[1]
            results_gates_tests[(11 + (i - 1) * 15):(15 + ((i - 1) * 15)), l] <-
                    c((confint(test, level = 1 - num_alpha))$confint[1:3],
                    (as.numeric(coef < 0) * (pval/2) + as.numeric(coef > 0) * (1 - pval/2)),
                    (as.numeric(coef > 0) * (pval/2) + as.numeric(coef < 0) * (1 - pval/2)) )

            results_ml_best[(1+(i-1)*2),l]  <- (sum(mean^2)/num_groups)

            ################################### Best Linear Prediction Regression  ###################################

            Sd            <- dataest$S- mean(dataest$S)
            dataest$S_ort <- I((as.numeric(dataest[,d])-md_x)*Sd)
            dataest$d_ort <- I((as.numeric(dataest[,d])-md_x))

            form1 <- as.formula(paste(y, "~", "B+S+d_ort+S_ort", form_var_add, "|", var_fe1, "+", var_fe2, "| 0 |", var_cluster, sep=""))

            a  <- tryCatch({
                a  <- felm(form1, data=dataest, weights=dataest$weight)
            },error=function(e){
                if(tog_silence == 0){cat("ERROR2 :",ml_methods[l], t, i, "\n")}
                form1 <- as.formula(paste(y, "~", "d_ort+S_ort", form_var_add, "|", var_fe1, "+", var_fe2, "| 0 |", var_cluster, sep=""))
                reg   <- felm(form1, data=dataest, weights=dataest$weight)
                return(reg)
            }, warning = function(war) {
                if(tog_silence == 0){cat("WARNING2 :",ml_methods[l], t, i, "\n")}
                form1 <- as.formula(paste(y, "~", "d_ort+S_ort| ", var_fe1, "+", var_fe2, "| 0 |", var_cluster, sep=""))
                reg   <- felm(form1, data=dataest, weights=dataest$weight)
                return(reg)
            })
            reg <- a

            coef <- (summary(reg)$coefficients['d_ort',1])
            pval <- (summary(reg)$coefficients['d_ort',4])
            results_blp_ate[(1+(i-1)*5):(i*5),l]      <-
                    c(coef, confint(reg, 'd_ort', level = 1-num_alpha)[1:2],
                    (as.numeric(coef < 0) * (pval/2) + as.numeric(coef > 0) * (1 - pval/2)),
                    (as.numeric(coef > 0) * (pval/2) + as.numeric(coef < 0) * (1 - pval/2)) )

            coef <- (summary(reg)$coefficients['S_ort',1])
            pval <- (summary(reg)$coefficients['S_ort',4])
            results_blp_het[(1 + (i - 1) * 5):(i * 5), l] <-
                    c(coef, confint(reg, "S_ort", level = 1 - num_alpha)[1:2],
                    (as.numeric(coef < 0) * (pval/2) + as.numeric(coef > 0) * (1 - pval/2)),
                    (as.numeric(coef > 0) * (pval/2) + as.numeric(coef < 0) * (1 - pval/2)) )

            results_ml_best[(2+(i-1)*2),l]      <- abs(summary(reg)$coefficients['S_ort',1])*sqrt(var(dataest$S))


            ################################################ CLANS  ################################################
            if(tog_endostrat == 0){
                high.effect     <- quantile(dataest$S, 1-num_thres);
                low.effect      <- quantile(dataest$S, num_thres);
                dataest$h       <- as.numeric(dataest$S>high.effect)
                dataest$l       <- as.numeric(dataest$S<low.effect)
            }
            if(tog_endostrat == 1){
                high.effect     <- quantile(dataest$B, 1-num_thres);
                low.effect      <- quantile(dataest$B, num_thres);
                dataest$h       <- as.numeric(dataest$B>high.effect)
                dataest$l       <- as.numeric(dataest$B<low.effect)
            }

            if(var(dataest$h) == 0){ dataest$h <- as.numeric(runif(length(dataest$h))<0.1) }
            if(var(dataest$l) == 0){ dataest$l <- as.numeric(runif(length(dataest$l))<0.1) }

            for(m in 1:length(var_affected)){
                a  <- tryCatch({
                    form <-  paste(var_affected[m],"~h+l-1", sep="")
                    reg  <-  lm(form, data=dataest[(dataest$h == 1)| (dataest$l == 1),])
                    coef <-  reg$coefficients['h'] - reg$coefficients['l']
                    test <-  glht(reg, linfct = c("h-l == 0"))

                    coef <-  (summary(reg)$coefficients['h',1])
                    pval <-  (summary(reg)$coefficients['h',4])
                    res1 <- c(coef, confint(reg, "h", level = 1 - num_alpha)[1:2],
                             (as.numeric(coef < 0) * (pval/2) + as.numeric(coef > 0) * (1 - pval/2)),
                             (as.numeric(coef > 0) * (pval/2) + as.numeric(coef < 0) * (1 - pval/2)) )

                    coef <-  (summary(reg)$coefficients['l',1])
                    pval <-  (summary(reg)$coefficients['l',4])
                    res2 <- c(coef, confint(reg, "l", level = 1 - num_alpha)[1:2],
                             (as.numeric(coef < 0) * (pval/2) + as.numeric(coef > 0) * (1 - pval/2)),
                             (as.numeric(coef > 0) * (pval/2) + as.numeric(coef < 0) * (1 - pval/2))  )

                    coef <- (summary(reg)$coefficients['h',1]) - (summary(reg)$coefficients['l',1])
                    pval <- summary(test)$test$pvalues[1]
                    res3 <- c((confint(test, level = 1 - num_alpha))$confint[1:3],
                             (as.numeric(coef < 0) * (pval/2) + as.numeric(coef > 0) * (1 - pval/2)),
                             (as.numeric(coef > 0) * (pval/2) + as.numeric(coef < 0) * (1 - pval/2)) )

                    a    <- c(res1, res2, res3)

                },error=function(e){
                    if(tog_silence == 0){cat("ERROR3 :",ml_methods[l], t, i, "\n")}
                    res1 <- c(mean(dataest[(dataest$h == 1), var_affected[m]]), mean(dataest[(dataest$h == 1), var_affected[m]]),
                              mean(dataest[(dataest$h == 1), var_affected[m]]), 0.5, 0.5 )
                    res2 <- c(mean(dataest[(dataest$l == 1), var_affected[m]]), mean(dataest[(dataest$l == 1), var_affected[m]]),
                              mean(dataest[(dataest$l == 1), var_affected[m]]), 0.5, 0.5 )
                    res3 <- c((res1[1] - res2[1]), (res1[1] - res2[1]), (res1[1] - res2[1]), 0.5, 0.5)
                    a    <- c(res1, res2, res3)
                    return(a)
                })
                results_clan[((i-1)*length(var_affected)*15+(m-1)*15+1):((i-1)*length(var_affected)*15+(m)*15),l]   <- a
            }
        }
    }
    results_all_vector <- c(as.vector(results_gates_tests), as.vector(results_blp_ate), as.vector(results_blp_het), as.vector(results_gates), as.vector(results_ml_best), as.vector(results_clan))
    print(t)
    results_an <- data.frame(results_all_vector)
}

####################################### Prepare Latex Tables for BLP and Most/Least Affected  #######################################

results_gates_tests  <- array(c(as.matrix(results_an[1:(15*length(var_Y)*length(ml_methods)),])), c(15*length(var_Y),length(ml_methods), num_reps))
results_blp_ate      <- array(c(as.matrix(results_an[((15*length(var_Y)*length(ml_methods))+1):((15+5)*length(var_Y)*length(ml_methods)),])), c(5*length(var_Y),length(ml_methods), num_reps))
results_blp_het      <- array(c(as.matrix(results_an[(((20)*length(var_Y)*length(ml_methods))+1):((20+5)*length(var_Y)*length(ml_methods)),])), c(5*length(var_Y),length(ml_methods), num_reps))
results_gates        <- array(c(as.matrix(results_an[(((25)*length(var_Y)*length(ml_methods))+1):((25+3*num_groups)*length(var_Y)*length(ml_methods)),])), c(3*num_groups*length(var_Y),length(ml_methods), num_reps))
results_ml_best      <- array(c(as.matrix(results_an[(((25+3*num_groups)*length(var_Y)*length(ml_methods))+1):((25+3*num_groups+2)*length(var_Y)*length(ml_methods)),])), c(2*length(var_Y),length(ml_methods), num_reps))
results_clan         <- array(c(as.matrix(results_an[(((25+3*num_groups+2)*length(var_Y)*length(ml_methods))+1):((25+3*num_groups+2+length(var_affected)*15)*length(var_Y)*length(ml_methods)),])), c(length(var_affected)*15*length(var_Y),length(ml_methods), num_reps))

results_blp_ate_all       <- t(sapply(seq(1:nrow(results_blp_ate[,,1])), function(x) colMedians(t(results_blp_ate[x,,]))))
results_blp_het_all       <- t(sapply(seq(1:nrow(results_blp_het[,,1])), function(x) colMedians(t(results_blp_het[x,,]))))
results_gates_tests_all   <- t(sapply(seq(1:nrow(results_gates_tests[,,1])), function(x) colMedians(t(results_gates_tests[x,,]))))
results_clan_all          <- t(sapply(seq(1:nrow(results_clan[,,1])), function(x) colMedians(t(results_clan[x,,]))))
results_gates_all         <- t(sapply(seq(1:nrow(results_gates[,,1])), function(x) colMedians(t(results_gates[x,,]))))
results_ml_best_all       <- t(sapply(seq(1:nrow(results_ml_best[,,1])), function(x) colMedians(t(results_ml_best[x,,]))))

index_ml_best1 <- order(-results_ml_best_all[1,])[1:2]
index_ml_best2 <- order(-results_ml_best_all[2,])[1:2]

if(tog_endostrat == 0){
    if(index_ml_best1[1] != index_ml_best2[1]){ index_ml_best <- c(index_ml_best1[1],index_ml_best2[1]) }
    if(index_ml_best1[1] == index_ml_best2[1]){ index_ml_best <- c(index_ml_best2[1],index_ml_best2[2]) }
}

if(tog_endostrat == 1){ index_ml_best <- c(index_ml_best1[1],index_ml_best1[2]) }

results_clan_all2         <- matrix(0,length(var_affected)*12*length(var_Y), length(ml_methods))
results_blp_ate_all2      <- matrix(NA,4*length(var_Y), length(ml_methods))
results_blp_het_all2      <- matrix(NA,4*length(var_Y), length(ml_methods))
results_gates_tests_all2  <- matrix(NA,12*length(var_Y), length(ml_methods))

l <- 1
for(i in seq(1, nrow(results_clan_all), 5)){
    results_clan_all2[l:(l+2),]  <- results_clan_all[i:(i+2),]
    results_clan_all2[l+3,]      <- sapply(seq(1:length(ml_methods)), function(x) min(1,4*min(results_clan_all[i+3,x], results_clan_all[i+4,x])))
    if(l < nrow(results_blp_ate_all2)){
        results_blp_ate_all2[l:(l+2),] <- results_blp_ate_all[i:(i+2),]
        results_blp_het_all2[l:(l+2),] <- results_blp_het_all[i:(i+2),]
        results_blp_ate_all2[l+3,]     <- sapply(seq(1:length(ml_methods)), function(x) min(1,4*min(results_blp_ate_all[i+3,x], results_blp_ate_all[i+4,x])))
        results_blp_het_all2[l+3,]     <- sapply(seq(1:length(ml_methods)), function(x) min(1,4*min(results_blp_het_all[i+3,x], results_blp_het_all[i+4,x])))
    }
    if(l < nrow(results_gates_tests_all2)){
        results_gates_tests_all2[l:(l+2),] <- results_gates_tests_all[i:(i+2),]
        results_gates_tests_all2[l+3,]     <- sapply(seq(1:length(ml_methods)), function(x) min(1,4*min(results_gates_tests_all[i+3,x], results_gates_tests_all[i+4,x])))
    }
    l <- l+4
}

results_blp_ate_all      <- round(results_blp_ate_all2, digits = 3)
results_blp_het_all      <- round(results_blp_het_all2, digits = 3)
results_gates_tests_all  <- round(results_gates_tests_all2, digits = 3)
results_ml_best_all      <- format(round(results_ml_best_all, pmax(0,4-nchar(floor(abs(results_ml_best_all))))), nsmall= pmax(0,4-nchar(floor(abs(results_ml_best_all)))))
results_clan_all         <- round(results_clan_all2, digits = 3)

results_gates_tests2     <- matrix(0,9*length(var_Y), length(ml_methods))
results_blp_ate2         <- matrix(0,3*length(var_Y), length(ml_methods))
results_blp_het2         <- matrix(0,3*length(var_Y), length(ml_methods))
results_clan_all2        <- matrix(0,9*length(var_Y)*length(var_affected), length(ml_methods))

seq3 <- seq(1, nrow(results_clan_all), 4)
l    <- 1
for(i in seq(1, nrow(results_clan_all2), 3)){
    k <- seq3[l]
    if(i < nrow(results_blp_ate2)){
        results_blp_ate2[i,]       <- format(round(results_blp_ate_all[k,], pmax(0,4-nchar(floor(abs(results_blp_ate_all[k,]))))), nsmall= pmax(0,4-nchar(floor(abs(results_blp_ate_all[k,])))))
        results_blp_ate2[i+1,]     <- sapply(seq(1:ncol(results_blp_ate_all)), function(x) paste("(", format(round(results_blp_ate_all[k+1,x],pmax(0,4-nchar(floor(abs(results_blp_ate_all[k+1,x]))))), nsmall=pmax(0,4-nchar(floor(abs(results_blp_ate_all[k+1,x]))))), ",", format(round(results_blp_ate_all[k+2,x],pmax(0,4-nchar(floor(abs(results_blp_ate_all[k+2,x]))))) , nsmall=pmax(0,4-nchar(floor(abs(results_blp_ate_all[k+2,x]))))), ")", sep=""))
        results_blp_ate2[i+2,]     <- paste("[", format(results_blp_ate_all[k+3,], nsmall = pmax(0,4-nchar(floor(abs(results_blp_ate_all[k+3,]))))), "]", sep="")

        results_blp_het2[i,]      <- format(round(results_blp_het_all[k,],max(0,4-nchar(floor(abs(results_blp_het_all[k,]))))) , nsmall=pmax(0,4-nchar(floor(results_blp_het_all[k,]))))
        results_blp_het2[i+1,]    <- sapply(seq(1:ncol(results_blp_het_all)), function(x) paste("(", format(round(results_blp_het_all[k+1,x], pmax(0,4-nchar(floor(abs(results_blp_het_all[k+1,x]))))) , nsmall=pmax(0,4-nchar(floor(abs(results_blp_het_all[k+1,x]))))), ",", format(round(results_blp_het_all[k+2,x],pmax(0,4-nchar(floor(abs(results_blp_het_all[k+2,x]))))) , nsmall=pmax(0,4-nchar(floor(abs(results_blp_het_all[k+2,x]))))), ")", sep=""))
        results_blp_het2[i+2,]    <- paste("[", format(results_blp_het_all[k+3,], nsmall=max(0,4-nchar(floor(abs(results_blp_het_all[k+3,]))))), "]", sep="")
    }
    if(i < nrow(results_gates_tests2)){
        results_gates_tests2[i,]    <- format(round(results_gates_tests_all[k,],pmax(0,4-nchar(floor(abs(results_gates_tests_all[k,]))))) , nsmall=pmax(0,4-nchar(floor(results_gates_tests_all[k,]))))
        results_gates_tests2[i+1,]  <- sapply(seq(1:ncol(results_gates_tests_all)), function(x) paste("(", format(round(results_gates_tests_all[k+1,x], pmax(0,4-nchar(floor(abs(results_gates_tests_all[k+1,x]))))),nsmall=pmax(0,4-nchar(floor(abs(results_gates_tests_all[k+1,x]))))), ",", format(round(results_gates_tests_all[k+2,x],pmax(0,4-nchar(abs(floor(results_gates_tests_all[k+2,x]))))),  nsmall=pmax(0,4-nchar(floor(abs(results_gates_tests_all[k+2,x]))))), ")", sep=""))
        results_gates_tests2[i+2,]  <- paste("[", format(results_gates_tests_all[k+3,], nsmall=pmax(0,4-nchar(floor(abs(results_gates_tests_all[k+3,]))))), "]", sep="")
    }
    results_clan_all2[i,]       <- format(round(results_clan_all[k,],pmax(0,4-nchar(floor(abs(results_clan_all[k,]))))) ,nsmall=pmax(0,4-nchar(floor(abs(results_clan_all[k,])))))
    results_clan_all2[i+1,]     <- sapply(seq(1:ncol(results_clan_all)), function(x) paste("(", format(round(results_clan_all[k+1,x], pmax(0,4-nchar(floor(abs(results_clan_all[k+1,x]))))), nsmall=pmax(0,4-nchar(floor(abs(results_clan_all[k+1,x]))))), ",", format(round(results_clan_all[k+2,x],pmax(0,4-nchar(floor(abs(results_clan_all[k+2,x]))))) , nsmall=pmax(0,4-nchar(floor(abs(results_clan_all[k+2,x]))))), ")", sep=""))
    if(i%%9==7){  results_clan_all2[i+2,]     <- paste("[", format(results_clan_all[k+3,], nsmall=pmax(0,4-nchar(floor(abs(results_clan_all[k+3,]))))), "]", sep="") }
    if(i%%9!=7){  results_clan_all2[i+2,]     <- "-" }

    l <- l+1
}

results_clan_final    <- matrix(NA, length(var_Y)*(length(var_affected)*3+1), length(index_ml_best)*3)
results_gates_final   <- matrix(NA, length(var_Y)*3, length(index_ml_best)*3)
results_blp_final     <- matrix(NA, length(var_Y)*3, length(index_ml_best)*2)
results_ml_best_final <- results_ml_best_all

rownames_CLAN    <- matrix(NA, nrow(results_clan_final),1)
rownames_GATES   <- matrix(NA, nrow(results_gates_final),1)
rownames_BEST    <- matrix(NA, nrow(results_ml_best_all),1)

a  <- 1
b  <- 1
c  <- 1
c2 <- 1

for(l in 1:length(var_Y)){
    rownames_CLAN[a] <- names_outcomes[l]
    a <- a+1
    for(i in 1:length(var_affected)){
        for(j in 1:length(index_ml_best)){
            k <- index_ml_best[j]
            results_clan_final[(a):(a+2),((j-1)*3+1):(j*3)] <- matrix(results_clan_all2[(b):(b+8),k], 3, 3)
            if(i == 1){
                results_gates_final[(c):(c+2),((j-1)*3+1):(j*3)] <- matrix(results_gates_tests2[(c2):(c2+8),k], 3, 3)
                rownames_GATES[c]   <- names_outcomes[l]
                results_blp_final[(c):(c+2),((j-1)*2+1):(j*2)] <- cbind(results_blp_ate2[(c):(c+2),k], results_blp_het2[(c):(c+2),k])
            }
            rownames_CLAN[a]   <- names_affected[i]
        }
        a <- a+3
        b <- b+9
    }
    c  <- c+3
    c2 <- c2+9
    rownames_BEST[((l-1)*2+1):((l-1)*2+2)] <- c(names_outcomes[l], names_outcomes[l])
}

rownames(results_clan_final)    <- rownames_CLAN
rownames(results_gates_final)   <- rownames_GATES
rownames(results_blp_final)     <- rownames_GATES
rownames(results_ml_best_final) <- rownames_BEST

colnames(results_clan_final)    <- rep(c("Most Affected", 	"Least Affected",	"Difference"), length(index_ml_best))
colnames(results_gates_final)   <- rep(c("Most Affected", 	"Least Affected",	"Difference"), length(index_ml_best))
colnames(results_blp_final)     <- rep(c("ATE", 	"HET"), length(index_ml_best))
colnames(results_ml_best_final) <- names_methods

print(xtable(cbind(rownames(results_blp_final),results_blp_final)),     include.rownames=FALSE,file=paste(name,"_BLP"  ,"-",name_output,".txt",sep=""), digits=3)
print(xtable(cbind(rownames(results_gates_final),results_gates_final)), include.rownames=FALSE,file=paste(name,"_GATES","-",name_output,".txt",sep=""), digits=3)
print(xtable(cbind(rownames(results_ml_best_final),results_ml_best_final)),   include.rownames=FALSE,file=paste(name,"_BEST" ,"-",name_output,".txt",sep=""), digits=3)
print(xtable(cbind(rownames(results_clan_final),results_clan_final)),   include.rownames=FALSE,file=paste(name,"_CLAN" ,"-",name_output,".txt",sep=""), digits=3)


for(i in 1:length(var_Y)){
    if(length(ml_methods)>1){ par(mfrow=c(2,2)) }

    y_range     <- 1*range(results_gates_all[(3*num_groups*(i-1)+num_groups+1):(3*num_groups*(i-1)+2*num_groups),],results_gates_all[(3*num_groups*(i-1)+2*num_groups+1):(3*num_groups*(i-1)+3*num_groups),])
    y_range2    <- y_range
    y_range2[1] <- y_range[1]- (y_range[2] -  y_range[1])*0.1
    y_range2[2] <- y_range[2]+ (y_range[2] -  y_range[1])*0.1

    result=list(0)

    for(j in 1:length(ml_methods)){
        ATE <- data.frame( x = c(-Inf, Inf), y = results_blp_ate_all[(4*(i-1)+1),j] , cutoff = factor(50))
        U   <- data.frame( x = c(-Inf, Inf), y = results_blp_ate_all[(4*(i-1)+3),j] , cutoff = factor(50))
        L   <- data.frame( x = c(-Inf, Inf), y = results_blp_ate_all[(4*(i-1)+2),j] , cutoff = factor(50))

        group_factor = "2"
        for(g in 2:num_groups){
            group_factor <- paste0(group_factor,",","2")
        }
        group_factor = factor(c(group_factor))

        label_crit <- (1-2*num_alpha)*100
        label_ci_gates <- paste0(label_crit,"% CI (GATES)")
        label_ci_ate   <- paste0(label_crit,"% CI (ATE)")

        result[[j]] <- ggplot() +
            theme_gray(base_size = 14) +
            geom_point(data=df,aes(y = F, x = x, colour=label_ci_gates), size = 3) +
            geom_errorbar(data=df, aes(ymax = U, ymin = L ,x=x, width=0.7, colour="GATES"), show.legend = TRUE) +
            geom_line(aes( x, y, linetype = cutoff, colour='ATE' ),ATE, linetype = 2) +
            geom_line(aes( x, y, linetype = cutoff, colour=label_ci_ate ), U, linetype = 2) +
            geom_line(aes( x, y, linetype = cutoff ), L, linetype = 2, color="red") +
            scale_colour_manual(values = c("red", "black", "blue", "black"),
                                breaks=c('ATE',label_ci_ate,"GATES",label_ci_gates),
                                guide = guide_legend(override.aes = list(
                                linetype = c("dashed", "dashed"  ,"blank", "solid"),
                                shape = c(NA,NA, 16, NA)), ncol =2,byrow=TRUE)) +
        theme(plot.title = element_text(hjust = 0.5,size = 11, face = "bold"), axis.title=element_text(size=10), legend.text=element_text(size=7), legend.key = element_rect(colour = NA, fill = NA), legend.key.size = unit(1, 'lines'), legend.title=element_blank(),legend.justification=c(0,1), legend.position=c(0,1), legend.background=element_rect(fill=alpha('blue', 0)))  +
        ylim(y_range) +
        labs(title=names_methods[j], y = "Treatment Effect", x = "Group by Het Score")
    }
    print(var_Y[i])
    if(length(ml_methods) >= 4){
        p      <- plot_grid(result[[1]], result[[2]], result[[3]], result[[4]], ncol=2)
        ggsave(paste(name,"_plot","-",name_output,"-",var_Y[i],".pdf",sep=""), p, width = 10, height =10)
    }
    p_best <- plot_grid(result[[index_ml_best[1]]], result[[index_ml_best[2]]], ncol=2)
    ggsave(paste(name,"_plot_best","-",name_output,"-",var_Y[i],".pdf",sep=""), p_best, width = 10, height = 5)
}
rm(df,L,p,p_best,result,rownames_BEST,rownames_CLAN,rownames_GATES,U)
rm(a,b,c,c2,i,j,k,l,label_ci_ate,label_ci_gates,label_crit,seq3,y_range,y_range2,group_factor)
rm(results_blp_ate_all,results_blp_het_all,results_clan_all,results_gates_tests_all)

time_end_an <- proc.time() - time_start_an
print(time_end_an)

# Save File
save.image(file=paste(name, "_an.RData", sep=""))

# Stop Cluster
stopCluster(cl)
