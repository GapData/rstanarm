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
#' @param ... For \code{loo_predict}, arguments passed to 
#'   \code{\link{posterior_predict}} (e.g. \code{newdata}). For
#'   \code{loo_linpred}, arguments passed to \code{\link{posterior_linpred}} 
#'   (e.g., \code{transform}, \code{newdata}).
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
  wts <- make_weights(object, lw)
  preds <- posterior_predict.stanreg(object, ...)
  colSums(preds * wts)
}

#' @rdname loo_predict.stanreg
#' @export
#'   
loo_linpred.stanreg <- function(object, lw, ...) {
  wts <- make_weights(object, lw)
  linpreds <- posterior_linpred.stanreg(object, ...)
  colSums(linpreds * weight_matrix)
}


# internal ----------------------------------------------------------------

# @param object,lw Same as above.
make_weights <- function(object, lw) {
  if (!missing(lw)) {
    stopifnot(is.matrix(lw))
    return(exp(lw))
  }
  psis <- loo::psislw(llfun = ll_fun(object), llargs = ll_args(object))
  exp(psis[["lw_smooth"]]) 
}

