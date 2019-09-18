#Test driven dev!

#' Compress extremes non-linearly
#'
#' Takes a Gradient Forest (GF) object and a
#' data.frame of environmental predictors
#' that are ready to be passed to
#' predict.gradientForest.
#'
#' This function finds the upper and lower
#' limits sampled by the GF object. Any prection
#' sites with a value outside of the range is compressed
#' to sit between capping (extreme prediction is identical
#' to max sampled prediction) and linear extrapolation.
#'
#' @param gf A gradientForest or combinedGradientForest object
#' @param env_grid Data.frame of environmental predictors to predict
#' @param pow numeric in range [0, 1]. Extremes will be raised to this power, eg 1/4 takes the 4th root. 0 corresponds to capping, 1 corresponds to linear extrapolation.
#' @param gf_weight Character, how should gradient forest objects be weighted. See ?cumimp.combinedGradientForest
#' @param ... paramaters passed to predict.gradientForest or predict.combinedGradientForest
#'
#' @return Data.frame of transformed environmental predictors
#'
#' @importFrom gradientForest importance.gradientForest importance.combinedGradientForest predict.gradientForest predict.combinedGradientForest
#'
#' @export
#'
#' @examples
#'
#' if (requireNamespace("gradientForest", quietly = TRUE)) {
#' library(gradientForest) #required to attach extendedForest
#' #8 species,
#' #9 env vars
#' #100 sample sites ranging from 1 to 2
#' #1000 grid sites ranging from 0 to 3
#' #poisson dependency on env vars
#'
#' set.seed(1000)
#'
#' species_dep <- matrix(runif(72, -10, 20), 9, 8)
#'
#' env_samp <- matrix(runif(900, 1, 2), 100, 9)
#'
#' species_response <- env_samp %*% species_dep
#' species_abundance <- data.frame(matrix(rpois(length(as.vector(species_response)), as.vector(species_response)), 100, 8))
#' names(species_abundance) <- LETTERS[1:8]
#'
#' env_samp <- as.data.frame(env_samp)
#' names(env_samp) <- letters[1:9]
#'
#' gf <- gradientForest::gradientForest(
#'     cbind(env_samp, species_abundance),
#'     letters[1:9],
#'     LETTERS[1:8]
#' )
#'
#' env_grid_upper <- matrix(runif(9000, 1.1, 3), 1000, 9)
#' env_grid_upper <- as.data.frame(env_grid_upper)
#' names(env_grid_upper) <- letters[1:9]
#' env_grid_lower <- matrix(runif(9000, 0, 1.8), 1000, 9)
#' env_grid_lower <- as.data.frame(env_grid_lower)
#' names(env_grid_lower) <- letters[1:9]
#'
#' #Testing upper extremes
#'
#' ##Extrapolating
#' pred_extrap <- predict(gf, env_grid_upper, extrap = TRUE)
#' compress_extrap <- gf_extrap_compress(gf, env_grid_upper, pow = 1)
#' testthat::expect_equal(as.data.frame(pred_extrap), as.data.frame(compress_extrap))
#'
#' ##Capping
#' pred_cap <- predict(gf, env_grid_upper, extrap = FALSE)
#' compress_cap <- gf_extrap_compress(gf, env_grid_upper, 0)
#'
#' testthat::expect_equal(as.data.frame(pred_cap), as.data.frame(compress_cap))
#'
#' ##Compression
#' compress_compress <- gf_extrap_compress(gf, env_grid_upper, 0.25)
#' #Here, extremes should lie between pred_cap and pred_extrap
#' testthat::expect_true(all(compress_compress >= pred_cap[,names(compress_cap)] & compress_compress <= pred_extrap[,names(compress_cap)]))
#'
#' hist(compress_compress$a, breaks = seq(-0.05, 0.05, 0.005), ylim = c(0, 600))
#' hist(pred_extrap$a, breaks = seq(-0.05, 0.05, 0.005), ylim = c(0,600))
#'
#' #' #Testing lower extremes
#' ##Capping
#' pred_cap <- predict(gf, env_grid_lower, extrap = FALSE)
#' compress_cap <- gf_extrap_compress(gf, env_grid_lower, 0)
#' testthat::expect_equal(as.data.frame(pred_cap[,names(compress_cap)]), as.data.frame(compress_cap), tolerance = 1e-4)
#'
#' ##Extrapolating
#' pred_extrap <- predict(gf, env_grid_lower, extrap = TRUE)
#' compress_extrap <- gf_extrap_compress(gf, env_grid_lower, 1)
#' testthat::expect_equal(as.data.frame(pred_extrap[,names(compress_cap)]), as.data.frame(compress_extrap))
#'
#' ##Compression
#' compress_compress <- gf_extrap_compress(gf, env_grid_lower, 0.25)
#' #Here, extremes should lie between pred_cap and pred_extrap
#' testthat::expect_true(all(compress_compress <= pred_cap & compress_compress >= pred_extrap))
#'
#' testthat::expect_error(gf_extrap_compress(gf, env_grid_lower, -1), "pow not greater than or equal to 0")
#' testthat::expect_error(gf_extrap_compress(gf, env_grid_lower, 1.1), "pow not less than or equal to 1")
#'
#' }
#'
gf_extrap_compress <- function(gf,
                               env_grid,
                               pow,
                               gf_importance_type = c("Weighted", "Raw", "Species")[1],
                               gf_weight = c("uniform", "species", "rsq.total",
                                             "rsq.mean", "site", "site.species", "site.rsq.total", "site.rsq.mean")[3],
                               ...){

  assertthat::assert_that(pow >= 0)
  assertthat::assert_that(pow <= 1)
  assertthat::assert_that(class(gf)[1] %in% c("gradientForest", "combinedGradientForest"))

  env_extrap <- predict(gf, env_grid, extrap = TRUE, ...)

  env_compress <- sapply(names(env_grid), env_extrap = env_extrap, function(varX, env_extrap) {
    tmp_y <- env_extrap[ , varX]
    tmp_x <- env_grid[, varX]
    #to compress extremes, I need x cap, y cap, extrapolation gradient
    #Adapted from gradientForest::predict.gradientForest()
    ci <- cumimp(gf, varX, ...)
    x_range <- range(ci$x)
    y_range <- range(ci$y)
    grad <- diff(y_range) / diff(x_range)
print(grad)
    upper_extremes <- tmp_x > x_range[2]
    if(length(upper_extremes) > 0){
      tmp_y[upper_extremes] <- rphildyerphd:::compress_extrap_z(tmp_x[upper_extremes] - x_range[2], pow, grad, y_range[2])
    }


    lower_extremes <- tmp_x < x_range[1]
    if(length(lower_extremes) > 0){
      tmp_y[lower_extremes] <- rphildyerphd:::compress_extrap_z( x_range[1] - tmp_x[lower_extremes], pow, grad, -y_range[1])
      tmp_y[lower_extremes] <- -tmp_y[lower_extremes]
    }
    return(tmp_y)
  })

  env_compress <- as.data.frame(env_compress)
  return(env_compress)

}

#' Compression Function
#'
#' Helper function to apply power compression to
#' GF predictions.
#'
#' The compressed value lies between a linear extrapolation (ext)
#' and capping at the maximum known value (cap)
#'
#' @description
#'
#' The logic is:
#'
#' We have two lines:
#' \deqn{y_{cap} = b}{y_cap = b}
#' \deqn{y_{ext} = ax + b}{y_ext = ax + b}
#'
#' where we set \eqn{x = 0} to be the point where extrapolation begins, so \eqn{y_{cap} = y_{ext}} at \eqn{x = 0}, and
#' \eqn{b} is the capped value.
#'
#' We want to define a third line, \eqn{z(x)}, st.
#'
#' \deqn{y_{cap} \le z(x) \le y_{ext} \all x \ge 0}{y_cap <= z(x) <= y_ext AA x >= 0}
#'
#' Two clear choices are sigmoid functions that asymptote, and power functions that do not asymptote.
#'
#' Here I will apply a power function to allow compositional turnover to grow indefinitely.
#'
#' Let \eqn{0 \le p \le 1}{0 <= p <= 1} be the power x is raised to. Then:
#'
#' \deqn{z = (x + c)^p + d}
#'
#' where \eqn{c,d} are constants.
#'
#' \eqn{x+c} for \eqn{p} between 0 and 1 grows faster than linear when \eqn{x+c} is close to 0, but it is also convex and monotonically increasing.
#' Therefore choosing \eqn{c,d} st. \eqn{y_{cap}(x)}{y_cap} is tangent to \eqn{z(x)} at \eqn{x = 0} will satisfy
#' \eqn{y_{cap} \le z(x) \le y_{ext} \all x \ge 0}{y_cap <= z(x) <= y_ext AA x >= 0}.
#'
#' At the tangent we know:
#'
#' \deqn{x = 0}
#' \deqn{y_{ext}(0) = z(0)}{y_ext(0) = z(0)}
#' \deqn{y_{ext}'(0) = z'(0)}{y'_ext(0) = z'(0)}
#'
#' \deqn{y_{ext}'(0) = z'(0)}{y'_ext(0) = z'(0)}
#' \deqn{a = p(x+c)^{p-1}}{a = p(x+c)^(p-1)}
#' \deqn{a = p(c)^{p-1} by x = 0}{a = p(c)^(p-1) by x = 0}
#' \deqn{c = \frac{a}{p}^\frac{1}{p-1}}{c = (a/p)^(1/(p-1))}
#'
#' \deqn{y_{ext}(0) = z(0)}{y_ext(0) = z(0)}
#' \deqn{b = c^p + d}{b = c^p + d}
#' \deqn{d = b - c^p}{d = b - c^p}
#' \deqn{d = b - \frac{a}{p}^\frac{p}{p-1}}{d = b - (a/p)^(p/(p-1))}
#'
#' therefore, to compress the extrapolation by power \eqn{p}, use:
#'
#' \deqn{z(x) = (x + \frac{a}{p}^\frac{1}{p-1})^p + b - \frac{a}{p}^\frac{p}{p-1}}{z(x) = (x + (a/p)^(1/(p-1)))^p + b - (a/p)^(p/(p-1))}
#'
#' When extrapolating into the negative domain, just take the negative of both y and x, then call this function, then negate x and y again.
#'
#' @param x numeric vector, all >= 0
#' @param p power, in range [0, 1]
#' @param a gradient of linear extrapolation
#' @param b cap value
#'
#' @return numeric vector of compressed values \eqn{b \le z(x) \le ax + b}{b <= z(x) <= ax + b}
#'
#' @examples
#'
#' a <- 1
#' p <- 0.25
#' b <- 1
#' x <- seq(0, 5, 0.1)
#'
#' testthat::expect_true(all(rphildyerphd:::compress_extrap_z(x, p, a, b) <= a*x+b))
#' testthat::expect_true(all(rphildyerphd:::compress_extrap_z(x, p, a, b) >= b))
#'
#' testthat::expect_error(rphildyerphd:::compress_extrap_z(x, p = 2, a, b), "p not less than or equal to 1")
#' testthat::expect_error(rphildyerphd:::compress_extrap_z(x, p, a = -1, b), "a not greater than 0")
#' testthat::expect_error(rphildyerphd:::compress_extrap_z(x, p, a = 0, b), "a not greater than 0")
#' testthat::expect_error(rphildyerphd:::compress_extrap_z(x, p = -1, a, b), "p not greater than or equal to 0")
#' testthat::expect_error(rphildyerphd:::compress_extrap_z(x = seq(-1, 1, 0.1), p, a, b), "Elements 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 of x >= 0 are not true")
#'
#' z <- rphildyerphd:::compress_extrap_z(x, p, a, b)
#' pl_all <- data.frame(x = x, y = a*x+b, z = z, cap = b)
#' matplot(pl_all$x, pl_all[,c("y","z", "cap")])
#'
compress_extrap_z <- function(x, p, a, b){

  assertthat::assert_that(p >=0, p<=1)
  assertthat::assert_that(a >0)
  assertthat::assert_that(all(x >= 0))

  if(p == 1){
    z <- a * x + b
  } else {
    z <- (x + (a / p) ^ (1/(p-1)) ) ^ p + b - (a / p) ^ (p / (p-1))
  }

  return(z)

}

