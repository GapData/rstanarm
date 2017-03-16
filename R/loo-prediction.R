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
#' @param type The type of expectation to compute. Currently the options are
#'   \code{"mean"} and \code{"var"} (variance).
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
loo_predict.stanreg <-
  function(object, 
           lw, 
           type = c("mean", "var"), 
           ...) {
    
    type <- match.arg(type)
    wts <- loo_weights(object, lw, ...)
    preds <- posterior_predict(object)
    loo_expectation(preds, wts, type)
  }

#' @rdname loo_predict.stanreg
#' @export
#' @param transform Passed to \code{\link{posterior_linpred}}.
#'    
loo_linpred.stanreg <-
  function(object,
           lw,
           type = c("mean", "var"),
           transform = FALSE,
           ...) {
    
    type <- match.arg(type)
    wts <- loo_weights(object, lw, ...)
    linpreds <- posterior_linpred(object, transform = transform)
    loo_expectation(linpreds, wts, type)
  }


# internal ----------------------------------------------------------------

# @param object,lw,... Same as above.
# @return A matrix.
loo_weights <- function(object, lw, ...) {
  if (!missing(lw)) {
    stopifnot(is.matrix(lw))
    return(exp(lw))
  }
  message("'lw' argument was not specified. Running PSIS to compute weights...")
  psis <- loo::psislw(llfun = ll_fun(object), llargs = ll_args(object), ...)
  exp(psis[["lw_smooth"]]) 
}


# @param X,W matrices of the same size
# @param type either "mean" or "var"
loo_expectation <- function(X, W, type) {
  fun <- switch(
    type,
    "mean" = function(x, w) sum(w * x),
    "var" = function(x, w) sum(w * (x - sum(w * x)) ^ 2)
  )
  vapply(seq_len(ncol(X)), function(j) {
    fun(X[, j], W[, j])
  }, FUN.VALUE = numeric(1))
}
