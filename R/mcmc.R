#' Inference procedure of FEBAMA framework
#' 
#' In the inference procedure, we take the MAP estimation with the standard BFGS algorithm.
#' If variable selection is considered (\code{model_conf$algArgs$nIter > 1}) simultaneously with MAP, 
#' we utilize the Gibbs sampler to perform variable selection over all forecasting models.  
#'
#' @param data A list with \code{lpd} and \code{feat} 
#' (the output of function \code{lpd_features_multi}).
#' @param model_conf Parameter settings for the FEBAMA framework. Defaults can be created by \code{model_conf_default()}.
#' 
#' @return \code{febama_mcmc} returns a list with the entries:
#' \describe{
#'   \item{beta}{A list of (number of models -1) matrices of beta in every iteration. }
#'   \item{betaIdx}{A list of (number of models -1) matrices of betaIdx in every iteration.}
#'   \item{accept_prob}{A list of (number of models -1) matrices of accept probabilities in every iteration.}
#' }
#' @export
febama_mcmc <- function(data, model_conf)
{
    ## Extract arguments
    algArgs = model_conf$algArgs
    varSelArgs = model_conf$varSelArgs
    priArgs = model_conf$priArgs

    nIter = algArgs$nIter
    num_models_updated = ncol(data$lpd) - 1

    ## Reserve space for OUT beta and initialize variable selection indicators
    OUT = list()
    for(iComp in 1:num_models_updated)
    {
        ## browser()
        ## Determine number of features used.
        nFeat = length(model_conf$features_used[[iComp]]) # default number of features = features_used

        ## nFeat = 0 if only intercept is used.
        if(length(model_conf$features_used[[iComp]]) == 0)
        {
            nFeat = 0
        }

        OUT[["beta"]][[iComp]] = matrix(NA, nIter, nFeat + 1)
        OUT[["betaIdx"]][[iComp]] = matrix(NA, nIter, nFeat + 1)

        ## Initialize variable selection indicators
        OUT[["betaIdx"]][[iComp]][1, ] = initialize_beta_idx(nFeat + 1, varSelArgs[[iComp]])

        ## Reserve space for acceptance probabilities with MH corrections
        OUT[["accept_prob"]][[iComp]] = matrix(1, nIter, 1)
    }

    betaIdx_curr = lapply(OUT[["betaIdx"]], function(x) x[1,])
    beta_curr = lapply(betaIdx_curr, function(x) stats::rnorm(length(x)))

    ## Numeric optimization to obtain MAP (Maximum a Posteriori)
    if(algArgs$initOptim == TRUE)
    {
        beta_optim = optimize_map_coefficients(data = data,
                                               beta = beta_curr,
                                               betaIdx = betaIdx_curr,
                                               priArgs = priArgs,
                                               varSelArgs = varSelArgs,
                                               features_used = model_conf$features_used,
                                               model_update = seq_along(betaIdx_curr))
        if( beta_optim$convergence != 0){
            stop("The optimization of initial values is not convergent")
        }
        beta_curr = beta_optim$beta
    }

    ## Assign initial values (conditional of variable selections)
    OUT[["beta"]] = mapply(function(x, y){
        x[1, ] = y
        return(x)
    }, x = OUT[["beta"]], y = beta_curr, SIMPLIFY = FALSE)
                       
    if(nIter > 1){
    for (iIter in 2:nIter)
    { # Loop start with the second iteration. 
      # The first iteration is considered as initial values.
      
      # SGLD did not seem to work well in our experiments, so we use simple MAP here. 
      if(algArgs$algName == "sgld"){
        beta_betaIdx <-  SGLD_gibbs(data = data,  beta_curr = beta_curr,
                                    betaIdx_curr = betaIdx_curr, model_conf = model_conf)
      }else if(algArgs$algName == "MAP"){
        beta_betaIdx <-  MAP_gibbs(data = data,  beta_curr = beta_curr,
                                    betaIdx_curr = betaIdx_curr, model_conf = model_conf)
      }

        ## Extract parameters for loop use
        beta_curr = beta_betaIdx[["beta"]]
        betaIdx_curr = beta_betaIdx[["betaIdx"]]
        accept_prob_curr = beta_betaIdx[["accept_prob"]]

        ## Assign the final output
        OUT[["beta"]] = mapply(function(x, y){
            x[iIter, ] = y
            return(x)
        }, x = OUT[["beta"]], y = beta_betaIdx[["beta"]], SIMPLIFY = FALSE)

        OUT[["betaIdx"]] = mapply(function(x, y){
            x[iIter, ] = y
            return(x)
        }, x = OUT[["betaIdx"]], y = beta_betaIdx[["betaIdx"]], SIMPLIFY = FALSE)

        OUT[["accept_prob"]] = mapply(function(x, y){
            x[iIter, ] = y
            return(x)
        }, x = OUT[["accept_prob"]], y = beta_betaIdx[["accept_prob"]], SIMPLIFY = FALSE)
    
    }
    }
    return(OUT)
}


# Variable selection using MAP with Gibbs
MAP_gibbs <- function(data, beta_curr, betaIdx_curr, model_conf)
{

    priArgs = model_conf$priArgs
    varSelArgs = model_conf$ varSelArgs

    features_used = model_conf$features_used

    num_models_updated = ncol(data$lpd) - 1

    beta_prop = beta_curr
    betaIdx_prop = betaIdx_curr
    accept_prob = list()
    for (iComp in 1:num_models_updated)
    {
        nPar_full = length(betaIdx_curr[[iComp]])

        candIdx = resolve_varsel_candidates(nPar_full, model_conf$varSelArgs[[iComp]]$cand)
        if(length(candIdx) > 0)
        {
            betaIdx_prop[[iComp]][candIdx] = stats::rbinom(length(candIdx), 1, prob = 0.5)
            betaIdx_prop[[iComp]][1] = 1
        }

        beta_optim = optimize_map_coefficients(data = data,
                                               beta = beta_prop,
                                               betaIdx = betaIdx_prop,
                                               priArgs = priArgs,
                                               varSelArgs = varSelArgs,
                                               features_used = model_conf$features_used,
                                               model_update = iComp)
        if( beta_optim$convergence != 0){
          stop("The optimization of MAP is not convergent")
        }
        
        beta_prop = beta_optim$beta

        ## Metropolis-Hasting accept/reject for variable selection
        if(length(candIdx) > 0)
        {
            log_post_prop = log_posterior(data = data,
                                          beta = beta_prop,
                                          betaIdx = betaIdx_prop,
                                          priArgs = priArgs,
                                          varSelArgs = varSelArgs,
                                          features_used = features_used,
                                          model_update = iComp)
            log_post_curr = log_posterior(data = data,
                                          beta = beta_curr,
                                          betaIdx = betaIdx_curr,
                                          priArgs = priArgs,
                                          varSelArgs = varSelArgs,
                                          features_used = features_used,
                                          model_update = iComp)


            ## The jump density for the variable selection indicators. TODO: Add adaptive scheme
            logJump.Idx.currATprop <- 1
            logJump.Idx.propATcurr <- 1

            logMHRatio <- (log_post_prop - log_post_curr +
                           logJump.Idx.currATprop - logJump.Idx.propATcurr)


            if(is.na(logMHRatio))
            { ## bad proposal, i.e logJump.currATpropRev = -Inf, or logJump.propATprop = -Inf
                accept_prob_curr <- 0
            } else{
                accept_prob_curr <- exp(min(0, logMHRatio))
            }

            if(stats::runif(1) < accept_prob_curr) #!is.na(accept.prob.curr)
            {  ## keep the proposal
                beta_curr[[iComp]] <- beta_prop[[iComp]]
                betaIdx_curr[[iComp]] <- betaIdx_prop[[iComp]]
            }
            else
            { ## keep the current
                beta_prop[[iComp]] = beta_curr[[iComp]]
                betaIdx_prop[[iComp]] = betaIdx_curr[[iComp]]
            }
        }
        else
        {
            accept_prob_curr = 1
            beta_curr[[iComp]] <- beta_prop[[iComp]] ## Accepted
            ## betaIdx_curr unchanged.
        }
        accept_prob[[iComp]] <- accept_prob_curr
        
    }

    out = list(beta = beta_curr, betaIdx = betaIdx_curr, 
               accept_prob = accept_prob)
    return(out)
}



SGLD_gibbs <- function(data, beta_curr, betaIdx_curr, model_conf)
{
  nObs = nrow(data$lpd)
  sgldArgs = validate_sgld_settings(nObs, model_conf$algArgs$sgld)
  stepsize = sgldArgs$stepsize
  burninProp = sgldArgs$burninProp
  
  priArgs = model_conf$priArgs
  varSelArgs = model_conf$varSelArgs
  
  features_used = model_conf$features_used
  
  batchSize = min(nObs, sgldArgs$max_batchSize)
  nBatch = ceiling(nObs / batchSize)
  nEpoch = sgldArgs$nEpoch
  nSgldIter = nEpoch * nBatch
  
  num_models_updated = ncol(data$lpd) - 1
  
  beta_prop = beta_curr
  betaIdx_prop = betaIdx_curr
  
  beta_sgld = list()
  accept_prob = vector("list", num_models_updated)
  for (iComp in 1:num_models_updated)
  {
    nPar_full = length(betaIdx_curr[[iComp]])
    
    ## 1. propose an update of variable selection indicators. Random proposal when
    ## variable selection is enabled: NOT NULL.
    candIdx = resolve_varsel_candidates(nPar_full, model_conf$varSelArgs[[iComp]]$cand)
    if(length(candIdx) > 0)
    {
      betaIdx_prop[[iComp]][candIdx] = stats::rbinom(length(candIdx), 1, prob = 0.5)
      betaIdx_prop[[iComp]][1] = 1
    }
    
    ## 2. conditional on this variable selection indicators, update beta via SGLD
    beta_iComp_sgld = matrix(0, nSgldIter, nPar_full)
    for (iEpoch in 1:nEpoch)
    {
      ## Re-split the data into small batches after finish one complete epoch.
      dataIdxLst = make_sgld_batches(nObs, batchSize)
      for (iBatch in seq_along(dataIdxLst))
      {
        iIter = (iEpoch - 1) * nBatch + iBatch
        stepSizeCurr = sgld_step_size(stepsize, iIter, sgldArgs)

        data_curr = lapply(data[c("lpd","feat")], function(x) x[dataIdxLst[[iBatch]], ,drop=FALSE])
        batchRatio = length(dataIdxLst[[iBatch]]) / nObs # n/N

        grad_iComp = log_posterior_grad(data = data_curr,
                                        beta = beta_prop,
                                        betaIdx = betaIdx_prop,
                                        priArgs = priArgs,
                                        varSelArgs = varSelArgs,
                                        features_used = features_used,
                                        model_update = iComp,
                                        batchRatio = batchRatio)[[1]]

        ## SGLD
        betaIdxActive = betaIdx_prop[[iComp]] == 1
        nPar1 = sum(betaIdxActive) # length of non-zero parameters
        beta_new <- (as.vector(beta_prop[[iComp]][betaIdxActive] + stepSizeCurr / 2 * grad_iComp) +
                       as.vector(mvtnorm::rmvnorm(1, rep(0, nPar1), stepSizeCurr * diag(nPar1))))

        beta_iComp_sgld[iIter, betaIdxActive] = beta_new

        beta_prop[[iComp]][betaIdxActive] = beta_new
        beta_prop[[iComp]][!betaIdxActive] = 0
      }
    }
    
    ## Polyak-Ruppert averaging improve the efficiency of SGLD
    nDrop = floor(burninProp * nSgldIter)
    keepRows = seq_len(nSgldIter)
    if(nDrop > 0)
    {
      keepRows = keepRows[-seq_len(nDrop)]
    }
    beta_prop[[iComp]] = colMeans(beta_iComp_sgld[keepRows,, drop = FALSE])
    beta_prop[[iComp]][betaIdx_prop[[iComp]] == 0] = 0
    beta_sgld[[iComp]] = beta_iComp_sgld
    
    ## Metropolis-Hasting accept/reject for variable selection
    if(length(candIdx) > 0)
    {
      log_post_prop = log_posterior(data = data,
                                    beta = beta_prop,
                                    betaIdx = betaIdx_prop,
                                    priArgs = priArgs,
                                    varSelArgs = varSelArgs,
                                    features_used = features_used,
                                    model_update = iComp)
      log_post_curr = log_posterior(data = data,
                                    beta = beta_curr,
                                    betaIdx = betaIdx_curr,
                                    priArgs = priArgs,
                                    varSelArgs = varSelArgs,
                                    features_used = features_used,
                                    model_update = iComp)
      
      
      ## The jump density for the variable selection indicators. TODO: Add adaptive scheme
      logJump.Idx.currATprop <- 1
      logJump.Idx.propATcurr <- 1
      
      logMHRatio <- (log_post_prop - log_post_curr +
                       logJump.Idx.currATprop - logJump.Idx.propATcurr)
      
      
      if(is.na(logMHRatio))
      { ## bad proposal, i.e logJump.currATpropRev = -Inf, or logJump.propATprop = -Inf
        accept_prob_curr <- 0
      }
      else
      {
        accept_prob_curr <- exp(min(0, logMHRatio))
      }
      
      if(stats::runif(1) < accept_prob_curr) #!is.na(accept.prob.curr)
      {  ## keep the proposal
        beta_curr[[iComp]] <- beta_prop[[iComp]]
        betaIdx_curr[[iComp]] <- betaIdx_prop[[iComp]]
      }
      else
      { ## keep the current
        beta_prop[[iComp]] = beta_curr[[iComp]]
        betaIdx_prop[[iComp]] = betaIdx_curr[[iComp]]
      }
    }
    else
    {
      accept_prob_curr = 1 # SGLD always accepts
      beta_curr[[iComp]] <- beta_prop[[iComp]] ## Accepted
      ## betaIdx_curr unchanged.
    }
    accept_prob[[iComp]] <- accept_prob_curr

  }

  out = list(beta = beta_curr, betaIdx = betaIdx_curr, 
             accept_prob = accept_prob, beta_sgld = beta_sgld)
  return(out)
}

initialize_beta_idx <- function(nPar, varSelArg)
{
  idx = rep(1, nPar)
  candIdx = resolve_varsel_candidates(nPar, varSelArg$cand)
  if(length(candIdx) == 0)
  {
    return(idx)
  }

  init = varSelArg$init
  if(length(init) == 0 || is.na(init))
  {
    init = "all-in"
  }

  if(init == "all-in")
  {
    idx[candIdx] = 1
  }
  else if(init == "all-out")
  {
    idx[candIdx] = 0
  }
  else if(init == "random")
  {
    idx[candIdx] = stats::rbinom(length(candIdx), 1, 0.5)
  }
  else
  {
    stop("No such init for betaIdx!", call. = FALSE)
  }

  idx[1] = 1
  idx
}

resolve_varsel_candidates <- function(nPar, cand)
{
  if(is.null(cand) || length(cand) == 0)
  {
    return(integer())
  }

  if(length(cand) == 1 && is.character(cand) && tolower(cand) == "2:end")
  {
    if(nPar < 2)
    {
      return(integer())
    }
    return(seq.int(2, nPar))
  }

  if(!is.numeric(cand) && !is.integer(cand))
  {
    stop("Variable-selection candidates must be numeric indices or \"2:end\".",
         call. = FALSE)
  }

  candIdx = as.integer(cand)
  if(any(is.na(candIdx)) || any(candIdx != cand) ||
     any(candIdx < 1) || any(candIdx > nPar))
  {
    stop("Variable-selection candidate indices are out of range.", call. = FALSE)
  }

  sort(unique(setdiff(candIdx, 1L)))
}

active_beta_vector <- function(beta, betaIdx, model_update = seq_along(betaIdx))
{
  unlist(lapply(model_update, function(iComp) {
    beta[[iComp]][betaIdx[[iComp]] == 1]
  }), use.names = FALSE)
}

replace_active_beta <- function(beta, betaIdx, par, model_update = seq_along(betaIdx))
{
  start = 1L
  for(iComp in model_update)
  {
    active = betaIdx[[iComp]] == 1
    nActive = sum(active)
    if(nActive > 0)
    {
      end = start + nActive - 1L
      beta[[iComp]][active] = par[start:end]
      start = end + 1L
    }
    beta[[iComp]][!active] = 0
  }
  beta
}

log_posterior_active <- function(par,
                                 data,
                                 beta,
                                 betaIdx,
                                 priArgs,
                                 varSelArgs,
                                 features_used,
                                 model_update)
{
  beta_full = replace_active_beta(beta, betaIdx, par, model_update)
  log_posterior(data = data,
                beta = beta_full,
                betaIdx = betaIdx,
                priArgs = priArgs,
                varSelArgs = varSelArgs,
                features_used = features_used,
                model_update = model_update)
}

log_posterior_active_grad <- function(par,
                                      data,
                                      beta,
                                      betaIdx,
                                      priArgs,
                                      varSelArgs,
                                      features_used,
                                      model_update)
{
  beta_full = replace_active_beta(beta, betaIdx, par, model_update)
  grad = log_posterior_grad(data = data,
                            beta = beta_full,
                            betaIdx = betaIdx,
                            priArgs = priArgs,
                            varSelArgs = varSelArgs,
                            features_used = features_used,
                            model_update = model_update)
  as.numeric(unlist(grad, use.names = FALSE))
}

optimize_map_coefficients <- function(data,
                                      beta,
                                      betaIdx,
                                      priArgs,
                                      varSelArgs,
                                      features_used,
                                      model_update,
                                      maxit = 1000)
{
  par = active_beta_vector(beta, betaIdx, model_update)
  beta_template = beta
  out = stats::optim(par = par,
                     fn = log_posterior_active,
                     gr = log_posterior_active_grad,
                     data = data,
                     beta = beta_template,
                     betaIdx = betaIdx,
                     priArgs = priArgs,
                     varSelArgs = varSelArgs,
                     features_used = features_used,
                     model_update = model_update,
                     method = "BFGS",
                     control = list(fnscale = -1, maxit = maxit))
  out$beta = replace_active_beta(beta_template, betaIdx, out$par, model_update)
  out
}

validate_sgld_settings <- function(nObs, sgldArgs)
{
  max_batchSize = as.integer(sgldArgs$max_batchSize)
  if(length(max_batchSize) != 1 || is.na(max_batchSize) || max_batchSize < 1)
  {
    stop("`max_batchSize` must be a positive integer.", call. = FALSE)
  }

  nEpoch = as.integer(sgldArgs$nEpoch)
  if(length(nEpoch) != 1 || is.na(nEpoch) || nEpoch < 1)
  {
    stop("`nEpoch` must be a positive integer.", call. = FALSE)
  }

  burninProp = sgldArgs$burninProp
  if(!is.numeric(burninProp) || length(burninProp) != 1 ||
     is.na(burninProp) || burninProp < 0 || burninProp >= 1)
  {
    stop("`burninProp` must be a number in [0, 1).", call. = FALSE)
  }

  stepsize = sgldArgs$stepsize
  if(!is.numeric(stepsize) || length(stepsize) != 1 ||
     (!is.na(stepsize) && stepsize <= 0))
  {
    stop("`stepsize` must be a positive number or NA.", call. = FALSE)
  }

  if(is.na(stepsize))
  {
    for (arg in c("a", "b", "gama"))
    {
      value = sgldArgs[[arg]]
      if(!is.numeric(value) || length(value) != 1 || is.na(value) || value <= 0)
      {
        stop("SGLD decay parameters `a`, `b`, and `gama` must be positive numbers.",
             call. = FALSE)
      }
    }
  }

  list(
    max_batchSize = min(nObs, max_batchSize),
    nEpoch = nEpoch,
    burninProp = burninProp,
    stepsize = stepsize,
    a = sgldArgs$a,
    b = sgldArgs$b,
    gama = sgldArgs$gama
  )
}

make_sgld_batches <- function(nObs, batchSize)
{
  idx = sample.int(nObs)
  split(idx, ceiling(seq_along(idx) / batchSize))
}

sgld_step_size <- function(stepsize, iIter, sgldArgs)
{
  if(is.na(stepsize))
  {
    return(sgldArgs$a * (sgldArgs$b + iIter) ^ (-sgldArgs$gama))
  }
  stepsize
}
