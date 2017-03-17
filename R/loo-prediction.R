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
#' fit <- stan_glm(mpg ~ wt, data = mtcars)
#' loo_preds <- loo_predict(fit)
#' preds <- colMeans(posterior_predict(fit))
#' 
#' library("ggplot2")
#' ggplot() + 
#'  geom_line(aes(x = seq_along(preds), y = preds), color = "black", size = 2) + 
#'  geom_line(aes(x = seq_along(loo_preds), y = loo_preds), color = "gray", size = 1) + 
#'  labs(x = "Data point", y = "Mean") +
#'  bayesplot::theme_default()
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
  intervals <- t(loo_predict(object, lw, type = "quantile", probs))
  colnames(intervals) <- labs
  return(intervals)
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
loo_expectation <- function(X, W, type, probs) {
  fun <- switch(
    type,
    "mean" = function(x, w)
      sum(w * x),
    "var" = function(x, w)
      sum(w * (x - sum(w * x)) ^ 2), 
    "quantile" = wprctile
  )
  if (type == "quantile") {
    formals(fun)[["probs"]] <- probs
    funval <- numeric(length(probs))
  } else {
    funval <- numeric(1)
  }
  vapply(seq_len(ncol(X)), function(j) {
    fun(X[, j], W[, j])
  }, FUN.VALUE = funval)
}


# @param x,w vectors of the same length
# @param p vector of probabilities
wprctile <- function(x, w, probs = 0.5) {
  stopifnot(length(x) == length(w), all(probs > 0 | probs < 1))
  
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
