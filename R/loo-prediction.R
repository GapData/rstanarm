#' Compute weighted expectations using LOO
#' 
#' @export
#' @aliases loo_predict loo_linpred
#' 
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param lw An optional matrix of (smoothed) log-weights. If \code{lw} is 
#'   missing then \code{\link[loo]{psislw}} is executed internally, which may be
#'   time consuming for large amounts of data.
#' @param ... Optional arguments passed to \code{\link[loo]{psislw}} if 
#'   \code{lw} is not specified.
#'   
#' @return A vector with \code{nrow(newdata)} elements.
#' 
#' @examples 
#' ### Logistic regression
#' head(wells)
#' wells$dist100 <- wells$dist / 100
#' fit <- stan_glm(
#'   switch ~ dist100 + arsenic, 
#'   data = wells, 
#'   family = binomial(link = "logit"), 
#'   QR = TRUE,
#'   chains = 2, iter = 200 # for speed of example only
#' )
#' pred <- loo_predict(fit)
#' 
loo_predict.stanreg <- function(object, lw, ...) {
  wts <- make_weights(object, lw, ...)
  preds <- posterior_predict(object)
  colSums(preds * wts)
}

#' @rdname loo_predict.stanreg
#' @export
#' @param transform Passed to \code{\link{posterior_linpred}}.
#'    
loo_linpred.stanreg <- function(object, lw, transform = FALSE, ...) {
  wts <- make_weights(object, lw, ...)
  linpreds <- posterior_linpred(object, transform = transform)
  colSums(linpreds * weight_matrix)
}


# internal ----------------------------------------------------------------

# @param object,lw,... Same as above.
# @return A matrix.
make_weights <- function(object, lw, ...) {
  if (!missing(lw)) {
    stopifnot(is.matrix(lw))
    return(exp(lw))
  }
  psis <- loo::psislw(llfun = ll_fun(object), llargs = ll_args(object), ...)
  exp(psis[["lw_smooth"]]) 
}

