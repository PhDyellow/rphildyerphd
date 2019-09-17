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

  env_compress <- sapply(names(max_cum_imp), env_extrap = env_extrap, function(varX, env_extrap) {
    tmp <- env_extrap[ , varX]
    upper_extremes <- tmp > max_cum_imp[varX]
    tmp[upper_extremes] <- (tmp[upper_extremes]/max_cum_imp[varX])^pow * max_cum_imp[varX] # ratio is preferable to #tmp[ex_i] <- (tmp[ex_i]+1)^0.25 - (max_cum_imp[varX]+1)^0.25 + max_cum_imp[varX]

    lower_extremes <- tmp < 0
    tmp[lower_extremes] <- -((-tmp[lower_extremes]+1)^pow - 1)
    return(tmp)
  })
  env_compress <- as.data.frame(env_compress)


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


