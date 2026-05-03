#! /usr/bin/env Rscript

################
### set path ###
################
path = "/gs/home/kangyf/lily/"
setwd(path)


##################
### load data ###
##################
load("febama/M4.RData")
load("febama/RData/M1392_lpd_features.RData")

index <- c()
for (i in 1:48000) {
  if (length(M4[[47000+i]]$x)>500){
    index <- c(index, 47000+i)
  }
}

ind = c(1:length(lpd_features))
ind = ind[-c(501, 1201)]
lpd_features3 = lpd_features
lpd_features3[[501]] = NULL
lpd_features3[[1200]] = NULL
index3 = index[-c(501, 1201)]

################
### packages ###
################
library("tsfeatures")
library("M4metalearning")
library("forecast")
library("tseries")
library("purrr")
library("ggplot2")
library("numDeriv")
library("mvtnorm")
library("base")
library("MASS")
library(parallel)


#################
### sourceDir ###
#################
source("febama/R/features.R")
source("febama/R/models.R")
source("febama/R/mcmc.R")
source("febama/R/priors.R")
source("febama/R/posterior.R")
source("febama/R/logscore.R")
source("febama/R/febama.R")

### Model config template ###
# shrinkage = 1
# nIter = 100
# max_batchSize = 108,
# nEpoch = 10,
# stepsize = 0.1
num_models = 3
model_conf_default = list(
  frequency = 12
  , ets_model = "ANN" 
  , forecast_h = 18 
  , train_h = 1 
  , history_burn = 60 
  , PI_level = 90 
  , roll = NULL 
  , feature_window = 60 
  , features_used = rep(list(c("entropy", "arch_acf", "alpha", "beta", "unitroot_kpss")), num_models - 1)
  , fore_model = c("ets_fore",  "naive_fore", "rw_drift_fore")
  , varSelArgs = rep(list(list(cand = "2:end", init = "all-in")), num_models - 1)
  
  , priArgs = rep(list(list("beta" = list(type = "cond-mvnorm",
                                          mean = 0, covariance = "identity", shrinkage = 1),
                            "betaIdx" = list(type = "beta", alpha0 = 1, beta0 = 1))), num_models - 1)
  , algArgs = list(initOptim = TRUE, 
                   algName = "sgld", 
                   nIter = 20, 
                   "sgld" = list(max_batchSize = 108,
                                 nEpoch = 10,
                                 burninProp = 0.4, 
                                 stepsize = 0.1,
                                 gama = 0.55,
                                 a = 0.4,
                                 b = 10)
  )
)


## Model without variable selection
model_conf_NoVS = model_conf_default
model_conf_NoVS[["varSelArgs"]] = rep(list(list(cand = NULL, init = "all-in")), num_models - 1)

## Model with only intercept (Bayesian optimal pool)
model_conf_NoFeat = model_conf_default
model_conf_NoFeat[["features_used"]] = rep(list(NULL), num_models - 1)
model_conf_NoFeat[["varSelArgs"]] = rep(list(list(cand = NULL, init = "all-in")), num_models - 1)

## Experiments ##
model_conf_curr = model_conf_default

t1 = proc.time()
OUT = mclapply(lpd_features3[1:36], febama_mcmc, model_conf = model_conf_curr, detail_out = T,
               mc.cores = 36)
# out = lapply(lpd_features3[1:36], febama_mcmc, model_conf = model_conf_curr, detail_out = T)
cat("The time of febama_mcmc is")
proc.time()-t1

save(model_conf_curr, OUT, file = "febama/RData/M1390_par3.RData")


t2 = proc.time()
RES = mcmapply(forecast_feature_results_multi, ts = M4[index3][1:36],
               data = lpd_features3[1:36],
               beta_out = OUT, mc.cores = 36,
               SIMPLIFY = FALSE, USE.NAMES = FALSE)
# res = mapply(forecast_feature_results_multi, ts = M4[index3][1:36], 
#                data = lpd_features3[1:36],
#                beta_out = out, SIMPLIFY = FALSE, USE.NAMES = FALSE)

cat("The time of forecast is")
proc.time()-t2

forecast_feature_performance(data = RES)

save(model_conf_curr, OUT, RES, file = "febama/RData/M1390_par3.RData")


