#' Compute weighted expectations using LOO
#' 
#' @export
#' @aliases loo_predict loo_linpred loo_pit
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
#' @return A vector with one element per observation.
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
    if (is_polr(object) && !is_scobit(object))
      preds <- polr_yrep_to_numeric(preds)
    
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

#' @rdname loo_predict.stanreg
#' @export
#' 
loo_pit.stanreg <- function(object, lw, ...) {
  lw <- loo_weights(object, lw, log = TRUE, ...)
  yrep <- posterior_predict(object)
  y <- get_y(object)
  if (is_polr(object) && !is_scobit(object)) {
    yrep <- polr_yrep_to_numeric(yrep)
    y <- as.integer(y)
  } else if (is_binomial_ppc(object) && NCOL(y) == 2) {
    y <- y[, 1]
  }
  rstantools::loo_pit(object = yrep, y = y, lw = lw)
}

# internal ----------------------------------------------------------------

# @param object,lw,... Same as above.
# @return A matrix.
loo_weights <- function(object, lw, log = FALSE, ...) {
  if (!missing(lw)) {
    stopifnot(is.matrix(lw))
  } else {
    message("'lw' argument was not specified. Running PSIS to compute weights...")
    psis <- loo::psislw(llfun = ll_fun(object), llargs = ll_args(object), ...)
    lw <- psis[["lw_smooth"]]
  }
  if (log) 
    return(lw) 
  else 
    return(exp(lw))
}


# @param X,W matrices of the same size
# @param type either "mean" or "var"
loo_expectation <- function(X, W, type) {
  fun <- switch(
    type,
    "mean" = function(x, w)
      sum(w * x),
    "var" = function(x, w)
      sum(w * (x - sum(w * x)) ^ 2)
  )
  vapply(seq_len(ncol(X)), function(j) {
    fun(X[, j], W[, j])
  }, FUN.VALUE = numeric(1))
}

# @param x A numeric vector
exp_log_sum_exp <- function(x) {
  max_x <- max(x, na.rm = TRUE)
  exp(max_x + log(sum(exp(x - max_x))))
}
