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
#' env_grid_upper <- matrix(runif(9000, 1, 3), 1000, 9)
#' env_grid_upper <- as.data.frame(env_grid_upper)
#' names(env_grid_upper) <- letters[1:9]
#' env_grid_lower <- matrix(runif(9000, 0, 2), 1000, 9)
#' env_grid_lower <- as.data.frame(env_grid_lower)
#' names(env_grid_lower) <- letters[1:9]
#'
#' #Testing upper extremes
#' ##Capping
#' pred_cap <- predict(gf, env_grid_upper, extrap = FALSE)
#' compress_cap <- gf_extrap_compress(gf, env_grid_upper, 0)
#' testthat::expect_equal(as.data.frame(pred_cap), as.data.frame(compress_cap))
#'
#' ##Extrapolating
#' pred_extrap <- predict(gf, env_grid_upper, extrap = FALSE)
#' compress_extrap <- gf_extrap_compress(gf, env_grid_upper, 1)
#' testthat::expect_equal(as.data.frame(pred_cap), as.data.frame(compress_cap))
#'
#' ##Compression
#' compress_compress <- gf_extrap_compress(gf, env_grid_upper, 0.25)
#' #Here, extremes should lie between pred_cap and pred_extrap
#' testthat::expect_true(all(compress_compress >= pred_cap & compress_compress <= pred_extrap))
#'
#' hist(compress_compress$a, breaks = seq(-0.05, 0.05, 0.005), ylim = c(0, 250))
#' hist(pred_extrap$a, breaks = seq(-0.05, 0.05, 0.005), ylim = c(0,250))
#'
#' #' #Testing lower extremes
#' ##Capping
#' pred_cap <- predict(gf, env_grid_lower, extrap = FALSE)
#' compress_cap <- gf_extrap_compress(gf, env_grid_lower, 0)
#' testthat::expect_equal(as.data.frame(pred_cap), as.data.frame(compress_cap))
#'
#' ##Extrapolating
#' pred_extrap <- predict(gf, env_grid_lower, extrap = FALSE)
#' compress_extrap <- gf_extrap_compress(gf, env_grid_lower, 1)
#' testthat::expect_equal(as.data.frame(pred_cap), as.data.frame(compress_cap))
#'
#' ##Compression
#' compress_compress <- gf_extrap_compress(gf, env_grid_lower, 0.25)
#' #Here, extremes should lie between pred_cap and pred_extrap
#' testthat::expect_true(all(compress_compress <= pred_cap & compress_compress >= pred_extrap))
#'
#' testthat::expect_error(gf_extrap_compress(gf, env_grid_lower, -1), "gf_extrap_compress: power must lie between 0 and 1")
#' testthat::expect_error(gf_extrap_compress(gf, env_grid_lower, 1.1), "gf_extrap_compress: power must lie between 0 and 1")
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
  #Check power is in [0,1]
  if (pow < 0 | pow > 1){
    stop("gf_extrap_compress: power must lie between 0 and 1")
  }

  env_extrap <- predict(gf, env_grid, extrap = TRUE, ...)

  #Rolands code
  # xrange <- range(Phys_grid[ ,varX]) #find range before transformation# Not used for curbing, just plotting
  # tmp <- Trns_grid[,varX] #get vector after transformation
  # ex_i <- tmp > max_cum_imp[varX] #max_cum_imp is already calculated. find values that exceed the max
  # tmp[ex_i] <- (tmp[ex_i]/max_cum_imp[varX])^0.25 * max_cum_imp[varX] # ratio is preferable to #tmp[ex_i] <- (tmp[ex_i]+1)^0.25 - (max_cum_imp[varX]+1)^0.25 + max_cum_imp[varX]

  if(class(gf)[1] == "combinedGradientForest" ) {
    max_cum_imp <- importance(gf, type = gf_importance_type, weight = gf_weight, sort = TRUE)
  } else if (class(gf)[1] == "gradientForest") {
    max_cum_imp <- importance(gf, type = gf_importance_type, sort = TRUE)
  }
  min_cum_imp <- as.numeric(predict(gf , as.data.frame(t(apply(gf$X, 2, min))-1), extrap = FALSE, ...)[names(max_cum_imp)])
  names(min_cum_imp) <- names(max_cum_imp)



  env_compress <- sapply(names(max_cum_imp), env_extrap = env_extrap, function(varX, env_extrap) {
    tmp <- env_extrap[ , varX]
    upper_extremes <- tmp > max_cum_imp[varX]
    #tmp[upper_extremes] <- (tmp[upper_extremes]/max_cum_imp[varX])^pow * max_cum_imp[varX] # ratio is preferable to #tmp[ex_i] <- (tmp[ex_i]+1)^0.25 - (max_cum_imp[varX]+1)^0.25 + max_cum_imp[varX]
    diff <- abs(tmp[upper_extremes] - max_cum_imp[varX])
    tmp[upper_extremes] <- diff^pow + max_cum_imp[varX]

    lower_extremes <- tmp < min_cum_imp[varX]
    #Roland has done:
    #negate lower extreme (which is likely negative)
    #add 1, so that "ratio" becomes 1 at min, and larger at smaller values
    #power of ratio, which compresses towards 1
    #remove 1, to undo offset
    #negate again.

    # I want the negative equivalent of the positive direction
    # which should be
    #diff <- |linear extrap - capped |
    #diff^pow
    #compress = capped + diff
    diff <- (min_cum_imp[varX] - tmp[lower_extremes])
    tmp[upper_extremes] <- diff^pow + max_cum_imp[varX]

    tmp[lower_extremes] <- -((-tmp[lower_extremes]+1)^pow - 1)
    return(tmp)
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

#' Maximum cumulative importance
#'
#' Internal Helper function
#'
#' Calculates the maximum cumulative importance for all variables
#'
#' @param gf
#' @param importance_type Character, one of c("Weighted", "Raw", "Species"), see importance.gradientForest and importance.combinedGradientForest.
#' @param gf_weight Character, how should gradient forest objects be weighted. See ?cumimp.combinedGradientForest
#' @param sort logical, sort variables by importance
#'
#' @return numeric vector of maximum cumuliative importances, named by variable
#'
#' @examples
#'
#' if (requireNamespace("gradientForest", quietly = TRUE)) {
#' library(gradientForest) #required to attach extendedForest
#'
#' set.seed(1000)
#' species_dep <- matrix(runif(72, -10, 20), 9, 8)
#'
#' env_samp <- matrix(runif(900, 1, 2), 100, 9)
#'
#' species_response <- env_samp %*% species_dep
#' species_abundance <- data.frame(matrix(rpois(length(as.vector(species_response)), as.vector(species_response)), 100, 8))
#' names(species_abundance) <- LETTERS[1:8]
#' env_samp <- as.data.frame(env_samp)
#' names(env_samp) <- letters[1:9]
#'
#'
#' gf1 <- gradientForest::gradientForest(
#'     cbind(env_samp, species_abundance),
#'     letters[1:9],
#'     LETTERS[1:4]
#' )
#'
#' gf2 <- gradientForest::gradientForest(
#'     cbind(env_samp, species_abundance),
#'     letters[1:9],
#'     LETTERS[1:4+4]
#' )
#'
#' gf1_2 <- combinedGradientForest(west = gf1, east = gf2)
#'
#' max_cum_imp <- get_max_cum_imp(gf1_2, importance_type = "Weighted", sort = TRUE)
#' expectation_max <-  c(1,2,3,4,5,6,7,8,9)
#' names(expectation_max) <- letters[1:9]
#' testthat::expect_equal(max_cum_imp, expectation_max, tolerance = 1e-4)
#' }
get_max_cum_imp <- function(x) UseMethod("get_max_cum_imp")


get_max_cum_imp.default <- function(x, ...){

  warning(paste("'get_max_cum_imp' does not know how to handle object of class ",
                class(x),
                "and can only be used on classes gradientForest and combinedGradientForest"))

}

#'
get_max_cum_imp.combinedGradientForest <- function(gf,
                                           importance_type = c("Weighted", "Raw", "Species")[1],
                                           gf_weight = c("uniform", "species", "rsq.total",
                                                         "rsq.mean", "site", "site.species", "site.rsq.total", "site.rsq.mean")[3],
                                           sort_vars = TRUE){
  comb_imp <- importance.combinedGradientForest(gf, type = importance_type, weight = gf_weight, sort = sort_vars)
  cumimp(gf1_2, "a")
  max_cum_imp <- sapply(names(comb_imp), function(var_x){
    tmp <- cumimp.combinedGradientForest(gf, var_x)
    return( max(tmp$y))
  })
  max_cum_imp <- sort(max_cum_imp, decreasing = T )
}


comb_imp <- importance.combinedGradientForest(NWS.combinedGf, type = c("Weighted", "Raw", "Species")[1], sort = T)
max_cum_imp <- comb_imp
cum_imps <- list()
for (varX in names(comb_imp)) {
  cum_imps[[varX]] <- cumimp.combinedGradientForest(NWS.combinedGf, varX, weight = "site.rsq.total")
  max_cum_imp[varX] <- max(cum_imps[[varX]]$y)
}
max_cum_imp <- sort(max_cum_imp, decreasing = T )
imp.vars <- names(max_cum_imp)

#' importance.gradientForest already gives the maximumum cumulative imp, sorted.
get_max_cum_imp.gradientForest <- function(gf,
                            importance_type = c("Accuracy", "Impurity", "Weighted", "Raw", "Species")[3],
                            sort_vars = TRUE){
  comb_imp <- importance.gradientForest(gf, type = importance_type, sort = sort_vars)
  return(comb_imp)
}


