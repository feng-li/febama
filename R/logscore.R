#' Calculate the log predictive score
#'
#' Calculate the log predictive score of a time series for FEBAMA framework. 
#' The weights in forecast combination are related to time series features.
#' 
#' @param data A list with \code{lpd} and \code{feat} 
#' (the output of function \code{lpd_features_multi}).
#' @param beta A list of coefficient vectors of the features.
#' @param betaIdx A list of indicator vectors. 
#' @param features_used The features used for forecast combination. 
#' See the parameter settings for more information.
#' @param sum If TRUE, return the sum of log predictive densities.
#' 
#' @return \code{logscore} returns the value of log predictive score.
#' @export

logscore <- function(data, beta, betaIdx, features_used, sum = TRUE)
{
   if(!is.list(beta))
    {
        beta = betaVec2Lst(beta, betaIdx)
    }

    eta = febama_linear_predictors(data, beta, betaIdx, features_used)
    log_weights = sweep(eta, 1, row_logsumexp(eta), "-")
    out = row_logsumexp(log_weights + data$lpd)

    if(sum == TRUE)
    {
        out = sum(out)
    }

    return(out)
}

# When only optimize the coefficients of a particular model. 
logscore_comp <- function(data, beta_comp, beta, betaIdx, features_used, sum = TRUE, model_update)
{
    if(!is.list(beta))
    {
        beta = betaVec2Lst(beta, betaIdx)
    }

    beta[[model_update]] = beta_comp
    eta = febama_linear_predictors(data, beta, betaIdx, features_used)
    log_weights = sweep(eta, 1, row_logsumexp(eta), "-")
    out = row_logsumexp(log_weights + data$lpd)
    
    if(sum == TRUE)
    {
        out = sum(out)
    }
    
    return(out)
}

# Gradient of the log score with respect to given models

logscore_grad <- function(data, beta, betaIdx, features_used, model_update = 1:length(betaIdx))
{
    nObs = nrow(data$lpd)
    eta = febama_linear_predictors(data, beta, betaIdx, features_used)
    log_weights = sweep(eta, 1, row_logsumexp(eta), "-")
    log_mix = row_logsumexp(log_weights + data$lpd)

    ## The gradient wrt me=x'beta
    grad0 = exp(log_weights[, model_update, drop = FALSE] +
                data$lpd[, model_update, drop = FALSE] -
                log_mix) -
            exp(log_weights[, model_update, drop = FALSE])

    ## The gradient wrt beta
    out = list()
    iCompdx = 0
    for(iComp in model_update)
    {
        iCompdx = iCompdx + 1
        features_used_curr = features_used[[iComp]]
        betaIdxCurr = betaIdx[[iComp]]
        features0 = cbind(rep(1, nObs), data$feat[, features_used_curr, drop = FALSE])
        out[[iCompdx]] = colSums(grad0[, iCompdx] * features0[, betaIdxCurr == 1, drop = FALSE])
    }
    return(out)
}

febama_linear_predictors <- function(data, beta, betaIdx, features_used)
{
    num_models_updated <- length(betaIdx)
    nObs = nrow(data$lpd)
    eta = matrix(0, nObs, num_models_updated + 1)

    for(iComp in 1:num_models_updated)
    {
        betaCurr = beta[[iComp]]
        betaIdxCurr = betaIdx[[iComp]]
        features_used_curr = features_used[[iComp]]
        features0 = cbind(rep(1, nObs), data$feat[, features_used_curr, drop = FALSE])

        me <- features0[, betaIdxCurr == 1, drop = FALSE] %*% matrix(betaCurr[betaIdxCurr == 1])
        me[me > 709] <- 709 # avoid overflow in downstream exp calls
        eta[, iComp] = as.numeric(me)
    }

    eta
}

row_logsumexp <- function(x)
{
    maxes = apply(x, 1, max)
    shifted = sweep(x, 1, maxes, "-")
    out = maxes + log(rowSums(exp(shifted)))
    out[!is.finite(maxes)] = maxes[!is.finite(maxes)]
    out
}
