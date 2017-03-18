#' Compute weighted expectations using LOO
#' 
#' @export
#' @aliases loo_predict loo_linpred loo_pit loo_predictive_interval
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
#' @return \code{loo_predictive_interval} returns a matrix with one row per
#'   observation and two columns. The other functions return a vector with one
#'   element per observation.
#' 
#' @examples
#' # data from help("lm")
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' d <- data.frame(
#'   weight = c(ctl, trt), 
#'   group = gl(2, 10, 20, labels = c("Ctl","Trt"))
#' ) 
#' fit <- stan_glm(weight ~ group, data = d)
#' loo_predictive_interval(fit, prob = 0.8)
#' 
#' 
loo_predict.stanreg <-
  function(object, 
           lw, 
           type = c("mean", "var", "quantile"), 
           probs = 0.5,
           ...) {
    
    type <- match.arg(type)
    wts <- loo_weights(object, lw, ...)
    preds <- posterior_predict(object)
    if (is_polr(object) && !is_scobit(object))
      preds <- polr_yrep_to_numeric(preds)
    
    loo_expectation(preds, wts, type, probs)
  }

#' @rdname loo_predict.stanreg
#' @export
#' @param prob For \code{loo_predictive_interval}, a scalar in \eqn{(0,1)}
#'   indicating the desired probability mass to include in the intervals. The
#'   default is \code{prob=0.9} (\eqn{90}\% intervals).
loo_predictive_interval.stanreg <- function(object, lw, prob = 0.9, ...) {
  stopifnot(length(prob) == 1)
  alpha <- (1 - prob) / 2
  probs <- c(alpha, 1 - alpha)
  labs <- paste0(100 * probs, "%")
  intervals <- loo_predict(object, lw, type = "quantile", probs)
  rownames(intervals) <- labs
  t(intervals)
}

#' @rdname loo_predict.stanreg
#' @export
#' @param transform Passed to \code{\link{posterior_linpred}}.
#' @param probs If \code{type} is \code{"quantile"}, a vector of probabilities
#'   in \eqn{(0,1)}. Ignored if \code{type} is \code{"mean"} or \code{"var"}.
#'    
loo_linpred.stanreg <-
  function(object,
           lw,
           type = c("mean", "var", "quantile"),
           probs = 0.5,
           transform = FALSE,
           ...) {
    
    type <- match.arg(type)
    wts <- loo_weights(object, lw, ...)
    linpreds <- posterior_linpred(object, transform = transform)
    loo_expectation(linpreds, wts, type, probs)
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
# @param log If FALSE (default) the weights are exponentiated before returning
# @return A matrix.
loo_weights <- function(object, lw, log = FALSE, ...) {
  if (!missing(lw)) {
    stopifnot(is.matrix(lw))
  } else {
    message("'lw' argument not specified. Running PSIS to compute weights...")
    psis <- loo::psislw(llfun = ll_fun(object), llargs = ll_args(object), ...)
    lw <- psis[["lw_smooth"]]
  }
  if (log) 
    return(lw) 
  
  exp(lw)
}

# Compute weighted expectations
#
# @param X,W matrices of the same size
# @param type either "mean", "var", or "quantile"
# @param probs a vector of probabilities if type="quantile"
loo_expectation <- function(X, W, type, probs) {
  stopifnot(identical(dim(X), dim(W)))
  fun <- switch(type,
                "mean" = .wmean,
                "var" = .wvar,
                "quantile" = .wquant)
  if (type == "quantile") {
    stopifnot(all(probs > 0 | probs < 1))
    formals(fun)[["probs"]] <- probs
    funval <- numeric(length(probs))
  } else {
    funval <- numeric(1)
  }
  vapply(seq_len(ncol(X)), function(j) {
    fun(X[, j], W[, j])
  }, FUN.VALUE = funval)
}


# loo-weighted mean, variance, and quantiles
#
# @param x,w vectors of the same length. this should be checked inside
#   loo_expectation() before calling these functions.
# @param probs vector of probabilities.
#
.wmean <- function(x, w) {
  sum(w * x)
}
.wvar <- function(x, w) {
  r <- (x - .wmean(x, w))^2
  sum(w * r)
}
.wquant <- function(x, w, probs) {
  x <- sort(x)  
  ww <- cumsum(w)
  ww <- ww / ww[length(ww)]

  y <- numeric(length(probs))
  for (j in seq_along(probs)) {
    ids <- which(ww >= probs[j])
    wi <- min(ids)
    if (wi == 1) {
      y[j] <- x[1]
    } else {
      w1 <- ww[wi - 1]
      x1 <- x[wi - 1]
      y[j] <- x1 + (x[wi] - x1) * (probs[j] - w1) / (ww[wi] - w1)
    }
  }
  return(y)
}
