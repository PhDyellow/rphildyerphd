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
                               ...){
  env_grid <- as.data.frame(env_grid)
  assertthat::assert_that(pow >= 0)
  assertthat::assert_that(pow <= 1)
  assertthat::assert_that(class(gf)[1] %in% c("gradientForest", "combinedGradientForest"))

  env_extrap <- predict(gf, env_grid, extrap = TRUE, ...)

  env_compress <- sapply(names(env_grid), env_extrap = env_extrap, function(varX, env_extrap) {
    tmp_y <- env_extrap[ , varX]
    tmp_x <- env_grid[, varX]
    #to compress extremes, I need x cap, y cap, extrapolation gradient
    #Adapted from gradientForest::predict.gradientForest()
    ci <- gradientForest::cumimp(gf, varX, ...)
    x_range <- range(ci$x)
    y_range <- range(ci$y)
    grad <- diff(y_range) / diff(x_range)

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
#' testthat::expect_true(all(rphildyerphd:::compress_extrap_z(x, p = 0, a, b) == b))
#'
#' testthat::expect_error(rphildyerphd:::compress_extrap_z(x, p = 2, a, b), "p not less than or equal to 1")
#' testthat::expect_error(rphildyerphd:::compress_extrap_z(x, p, a = -1, b), "a not greater than or equal to 0")
#' #testthat::expect_true(all(rphildyerphd:::compress_extrap_z(x, p, a = 0, b) == b))
#' testthat::expect_error(rphildyerphd:::compress_extrap_z(x, p = -1, a, b), "p not greater than or equal to 0")
#' testthat::expect_error(rphildyerphd:::compress_extrap_z(x = seq(-1, 1, 0.1), p, a, b), "Elements 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 of x >= 0 are not true")
#'
#' z <- rphildyerphd:::compress_extrap_z(x, p, a, b)
#' pl_all <- data.frame(x = x, y = a*x+b, z = z, cap = b)
#' matplot(pl_all$x, pl_all[,c("y","z", "cap")])
#'
compress_extrap_z <- function(x, p, a, b){

  assertthat::assert_that(p >=0, p<=1)
  assertthat::assert_that(a >=0)
  assertthat::assert_that(all(x >= 0))

  if(p == 1){
    z <- a * x + b
  } else if (a == 0 | p == 0){
    z <- b
  } else {
    z <- (x + (a / p) ^ (1/(p-1)) ) ^ p + b - (a / p) ^ (p / (p-1))
  }

  return(z)

}


#I need to cluster with a combination of MVPART, F-ratio and confusion matrix.
#That sounds like at least 4 functions, and MVPART feeds into F-ratio and confusion matrix. The 4th function is an integrator function

#' Gradient Forest MVPART
#'
#' Apply an MVPART analysis to a dataset that is also being passed to Gradient Forest
#'
#' @param gf A GF model, not a combined GF model though
#'
#' @return factor mapping between leaf nodes and samples sites. Each entry is a sample site, each level is a leaf node.
#'
#' @export
#'
#' @examples
#'
#' if (requireNamespace("gradientForest", quietly = TRUE)) {
#' library(gradientForest) #required to attach extendedForest
#' set.seed(1001)
#'
#' species_dep <- matrix(runif(144, -5, 20), 9, 16)
#'
#' env_samp <- matrix(runif(900, 1, 2), 100, 9)
#'
#' species_response <- env_samp %*% species_dep
#' species_abundance <- data.frame(matrix(rpois(length(as.vector(species_response)), as.vector(species_response)), 100, 16))
#' names(species_abundance) <- LETTERS[1:16]
#'
#' env_samp <- as.data.frame(env_samp)
#' names(env_samp) <- letters[1:9]
#'
#' gf1 <- gradientForest::gradientForest(
#'     cbind(env_samp, species_abundance),
#'     letters[1:9],
#'     LETTERS[1:8]
#' )
#' gf2 <- gradientForest::gradientForest(
#'     cbind(env_samp, species_abundance),
#'     letters[1:9],
#'     LETTERS[1:8 + 8]
#' )
#'
#' gf3 <- gradientForest::combinedGradientForest(west = gf1, east = gf2)
#'
#' #The returned object is not trivial or easy to construct in other ways, so testing with hash functions.
#' testthat::expect_known_hash(gf_mvpart(gf1), "0b7d05ed")
#' testthat::expect_known_hash(gf_mvpart(gf2), "02c9ba96")
#' testthat::expect_error(gf_mvpart(gf3), "class(gf)[1] not equal to \"gradientForest\"", fixed = TRUE)
#' }

gf_mvpart <- function(gf){
  assertthat::assert_that(assertthat::are_equal(class(gf)[1], "gradientForest"))

  x <- gf$X
  y <- gf$Y

  y <- sapply(y, function(y_col){
    unclass(y_col)
  })

  spp <- names(gf$result)
  y <- y[,spp]
  y <- as.matrix(y)

  mvp <- mvpart::mvpart(y ~ .,
                x,
                xv=c("1se", "min")[1],
                xval=10,
                xvmult=10,
                xvse=1,
                plot.add=FALSE,
                text.add=FALSE,
                pretty=FALSE)

  #only terminal nodes appear in mvp$where, but other nodes occupy id slots leading to "missing" id numbers
  #relevel to avoid confusion over "missing" nodes by end users
  fac <- factor(mvp$where)
  levels(fac) <- 1:nlevels(fac)
  return(fac)

}


#' F-ratio over many clusterings
#'
#' This function ties together a number of helper functions
#' to get the f-ratio for a gf fit.
#'
#' It takes a gf model, an environmental grid, and a set of parameters for passing down into
#' helper functions.
#'
#' You get bace a long table of F-ratio scores for each k.
#'
#' @param gf a gf model. Currently only single GF model, not combinedGradientForest
#' @param gf_grid_sites data.frame, same number of rows as gf$Y. Columns are keys from env_grid that map each gf sample to a grid cell. These columns will be excluded from transformation.
#' Usually lat,lon or a grid id. See examples if rows are already aligned.
#' @param env_grid environmental data in original units. Must have columns with the same names as gf_grid_sites for matching grid cells to biological samples.
#'
#'
#' @param k_range integer vector, cluster values to fit
#' @param reps integer, number of fitting repetitions for each entry in k_range
#' @param pow compression power for extrapolation, between 0 (capping) and 1 (linear). For example, 0.5 giver square root, 0.25 gives 4th root.
#' @param parallel TRUE/FALSE Run in parallel. Assumes a foreach registerDo* backend has been defined.
#'
#' @param gf_predict_args list of additional arguments to predict.gradientForest
#'
#' @return returns list containing a list of all clusterings as clara objects. Each clara object has been extended with an anova element from gf_anova() and a confusion score from opt_confusion().
#' The return list also contains a factor vector of the mvpart assignments to sample sites, generated by gf_mvpart.
#'
#' TODO
#'
#'
#' @export
#'
#' @examples
#'
#' if (requireNamespace("gradientForest", quietly = TRUE)) {
#' library(gradientForest) #required to attach extendedForest
#'
#' cluster_tests <- 3:5
#'
#' data(CoMLsimulation)
#' preds <- colnames(Xsimulation)
#' specs <- colnames(Ysimulation)
#' f1 <- gradientForest(data.frame(Ysimulation,Xsimulation), preds, specs[1:6], ntree=10)
#' f2 <- gradientForest(data.frame(Ysimulation,Xsimulation), preds, specs[1:6+6], ntree=10)
#' f12 <- combinedGradientForest(west=f1,east=f2)
#'
#' env_grid <- f1$X
#' env_grid$site_id <- 1:nrow(env_grid)
#' gf_grid_sites <- data.frame(site_id = 1:nrow(f1$Y))
#' k_range <- 2:10
#' reps <- 1
#' is_parallel <- FALSE
#' pow <- 0.25
#' set.seed(1000)
#' test <- gf_clust_f_ratio(gf = f1, gf_grid_sites = gf_grid_sites, env_grid = env_grid, k_range = k_range, reps =  reps, is_parallel = is_parallel, pow = pow)
#' testthat::expect_named(test, c("mvpart", "clust_list"))
#'
#' testthat::expect_equal(length(test$clust_list), length(k_range)*reps)
#'
#' testthat::expect_true(all(sapply(test$clust_list, function(x){class(x) == c("clara", "partition")} )))
#' testthat::expect_equal(class(test$mvpart), c("factor") )
#' set.seed(1000)
#' test2 <- gf_clust_f_ratio(gf = f1, gf_grid_sites = gf_grid_sites, env_grid = env_grid, k_range = k_range, reps =  reps, is_parallel = is_parallel, pow = pow*0.5)
#' #testthat::expect_true(test$clust_list[[1]]$anova$f_ratio != test2$clust_list[[1]]$anova$f_ratio)
#'
#' #Adjust clara fittings
#' set.seed(1000)
#' test3 <- gf_clust_f_ratio(gf = f1, gf_grid_sites = gf_grid_sites, env_grid = env_grid, k_range = k_range, reps =  reps+2, is_parallel = is_parallel, pow = pow,
#' clara_args = list(samples = 20, sampsize = 50, trace = 0, rngR = TRUE, pamLike = TRUE, correct.d = TRUE))
#'
#' #Extract f-ratios per k into long form
#' f_ratio_test <- do.call(rbind, lapply(test3$clust_list, function(clust){
#'   return(data.frame(k = clust$anova$cluster, f_ratio = clust$anova$f_ratio))
#' }))
#'
#' plot(f_ratio_test)
#'
#' #Check confusion scores
#' testthat::expect_equal(sapply(test3$clust_list, function(clust){clust$conf$score}),
#'                  c(29, 40, 48, 55, 66, 76, 87, 85, 86, 29, 40, 48, 55, 67, 75, 86, 86, 87, 29, 40, 48, 55, 66, 75, 86, 87, 87))
#'
#' #chech tolerance to mapping vars
#' env_grid <- f1$X
#' env_grid$site_idd <- 1:nrow(env_grid)
#' gf_grid_sites <- data.frame(site_idd = 1:nrow(f1$Y))
#' set.seed(1000)
#' test4 <- gf_clust_f_ratio(gf = f1, gf_grid_sites = gf_grid_sites, env_grid = env_grid, k_range = k_range, reps =  reps+2, is_parallel = is_parallel, pow = pow,
#' clara_args = list(samples = 20, sampsize = 50, trace = 0, rngR = TRUE, pamLike = TRUE, correct.d = TRUE))
#'
#' testthat::expect_equal( f_ratio_test,  do.call(rbind, lapply(test4$clust_list, function(clust){
#'   return(data.frame(k = clust$anova$cluster, f_ratio = clust$anova$f_ratio))
#' })))
#'
#' }
#'
#'
#'
gf_clust_f_ratio <- function(gf,
                            gf_grid_sites,
                            env_grid = gf$X,
                            k_range,
                            reps = 1,
                            is_parallel = TRUE,
                            pow = 0.25,
                            gf_predict_args = list(),
                            clara_args = list()){

  spatial_vars <- names(gf_grid_sites)
  #Predict + compress

  env_trans <- do.call("gf_extrap_compress", c(list(gf = gf, env_grid = env_grid[!names(env_grid) %in% spatial_vars], pow = pow), gf_predict_args))

  #reattach spatial data
  env_trans <- cbind(env_grid[names(env_grid) %in% spatial_vars], env_trans)

  #Cluster
  cluster_list <- do.call("cluster_range", c(list(x = env_trans[, !names(env_trans) %in% spatial_vars], k = k_range, reps = reps, is_parallel = is_parallel), clara_args))

  #MVPART
  pdf(file = NULL)
  mvpart_result <- gf_mvpart(gf)
  dev.off()
  #fratio

  cluster_fratio <- lapply(cluster_list, spatial_vars = spatial_vars, function(clust, spatial_vars){
    #get f ratio and a
    ret <- clust
    #need to align clusters and samples somehow

    #cluster of samples
    sample_clust <- merge(data.frame(env_trans[, spatial_vars, drop = FALSE], clust = clust$clustering), gf_grid_sites, by = spatial_vars)
    node_clust <- merge(data.frame(gf_grid_sites, leaf = mvpart_result), sample_clust, by = spatial_vars)
    ret$conf <- rphildyerphd:::opt_confusion(node_clust$leaf, node_clust$clust)
    ret$anova <- gf_anova(gf  , k = length(clust$i.med), sample_clust$clust)

    return(ret)
  })



  return(list(mvpart = mvpart_result, clust_list = cluster_fratio))

}


#' Gradient Forest Clustering F-ratio
#'
#' Calculate the F-ratio for predicting the species response using only the clustering.
#'
#' @param gf A gradient forest model
#' @param k number of clusters in fitting
#' @param clust integer vector assigning sample sites to clusters, must align with gf_spp
#'
#' @return list containing var_model, var_resid, f_ratio, p_value, cluster, inertia_exp
#'
#' @importFrom vegan capscale
#'
#' @examples
#'
#' set.seed(1000)
#'
#' if (requireNamespace("gradientForest", quietly = TRUE)) {
#' library(gradientForest) #required to attach extendedForest
#'
#' cluster_tests <- 3:5
#'
#' data(CoMLsimulation)
#' preds <- colnames(Xsimulation)
#' specs <- colnames(Ysimulation)
#' f1 <- gradientForest(data.frame(Ysimulation,Xsimulation), preds, specs[1:6], ntree=10)
#' f2 <- gradientForest(data.frame(Ysimulation,Xsimulation), preds, specs[1:6+6], ntree=10)
#' f12 <- combinedGradientForest(west=f1,east=f2)
#' f1.pred<-predict(f1) # defaults to predicting the samples
#'
#' clust_list <- rphildyerphd:::cluster_range(f1.pred, k = cluster_tests) #no parameters passed to cluster::clara, not important
#'
#'
#' expected_f <- c(23.328, 18.5, 19.1)
#'
#' for(i in 1:length(cluster_tests)){
#'   results <- rphildyerphd:::gf_anova(gf = f1, k = cluster_tests[i], clust = clust_list[[i]]$clustering)
#'   testthat::expect_named(results, c("var_model", "var_resid", "f_ratio", "p_value", "cluster", "inertia_exp"))
#'   testthat::expect_equal(results$cluster, cluster_tests[i])
#'   testthat::expect_equal(results$f_ratio, expected_f[i], tolerance = 1e-3)
#' }
#'
#' results <- rphildyerphd:::gf_anova(gf = f1, k = cluster_tests[3], clust = clust_list[[3]]$clustering)
#' testthat::expect_equal(results, list(var_model = 16.1579,
#'    var_resid = 20.10891,
#'    f_ratio = 19.08359,
#'    p_value = 0.001,
#'    cluster = 5,
#'    inertia_exp = 0.5776277
#'    ),
#'   tolerance = 1e-3)
#'
#' }
gf_anova <- function(gf,
                     k,
                     clust) {
  gf_spp <- gf$Y[names(gf$result)] #species with R^2 > 0
  xy <- cbind(gf_spp, cluster = clust)
  xy <- xy[rowSums(gf_spp) > 0, ]
  lhs <- "xy[names(gf_spp)]"
  rhs <- "factor(cluster)"
  form <- as.formula(paste(lhs, "~", rhs))

  fit_rda <- vegan::capscale(form, xy, dist="bray")
  varexp <- if (is.null(fit_rda$CCA)) {
    NA
  } else {
    fit_rda$CCA$tot.chi / fit_rda$tot.chi
  }

  rda_anova <- anova(fit_rda)



  ret <- list(var_model = ifelse(!is.null(rda_anova$SumOfSqs[1]),
                            rda_anova$SumOfSqs[1], NA),
              var_resid = ifelse(!is.null(rda_anova$SumOfSqs[2]),
                            rda_anova$SumOfSqs[2], NA),
              f_ratio = rda_anova$F[1],
              p_value = rda_anova[["Pr(>F)"]][1],
              cluster = k,
              inertia_exp = varexp)

  return(ret)

}


#' Calculate a range of k-medoid clusterings
#'
#' Uses the clara method to find clusterings for a dataset, but fits a range of k and will run in parallel if either a foreach registerDo*() backend is set or
#' future::plan() is set.
#'
#'
#' @param x dataset to cluster, passed to cluster::clara
#' @param k integer vector, each element will create a fitting with k[i] clusters
#' @param reps number of replicates for each element of k (redundant for clara?)
#' @param parallel boolean. Allow parallel excecution. Requires a foreach registerDo*() backend to be set. doFuture::registerDoFuture() is recommended for flexibility.
#' @param ... arguments passed to clara
#'
#' @return list of cluster::clara objects, NOT ordered with respect to k vector or rep.
#'
#'
#' @export
#' @importFrom foreach foreach %dopar% %do%
#' @examples
#' set.seed(1000)
#'
#' samples <- cluster::xclara
#'
#' test_one  <- cluster_range(samples, k = c(2), reps = 1, is_parallel = FALSE)
#' testthat::expect_equal(class(test_one), "list")
#' testthat::expect_equal(class(test_one[[1]]), c("clara", "partition"))
#' plot(test_one[[1]], which.plot = 1)
#'
#' testthat::expect_error(cluster_range(x, k = c(2), reps = 1, is_parallel = FALSE, not_a_param = TRUE), "unused argument (not_a_param = TRUE)", fixed = TRUE)
#'
#' test_many <- cluster_range(samples, k = 2:10, reps = 1, is_parallel = FALSE) #this is also a test, because it can fail
#' testthat::expect_equal(class(test_many), "list")
#' testthat::expect_true(all(sapply(test_many, function(x){class(x) == c("clara", "partition")})))
#' testthat::expect_equal(length(test_many), length(2:10))
#'
#' test_many_reps <- cluster_range(samples, k = 2:10, reps = 5, is_parallel = FALSE)
#'
#' testthat::expect_equal(class(test_many_reps), "list")
#' testthat::expect_true(all(sapply(test_many_reps, function(x){class(x) == c("clara", "partition")})))
#' testthat::expect_equal(length(test_many_reps), length(2:10)*5)
#'
#' k_set <- sapply(test_many_reps, function(x){length(x$i.med)})
#' testthat::expect_equal(sort(k_set), sort(rep(2:10, 5)))
#'
#' #Parallel aware
#'
#' if (requireNamespace("future", quietly = TRUE) & requireNamespace("doFuture", quietly = TRUE)) {
#' library(future)
#' library(doFuture)
#' doFuture::registerDoFuture()
#'
#' future::plan(multisession, workers = 2)
#'
#' testthat::expect_gt(system.time( cluster_range(samples, k = 4, reps = 4, samples = 500, is_parallel = FALSE))[1],
#'               system.time(cluster_range(samples, k = 4, reps = 4, samples = 500, is_parallel = TRUE))[1])
#'}
cluster_range <- function(x, k, reps = 1, is_parallel = TRUE, ...) {

  # two choices if parallel is false:
  #   1. if statement, %do% vs %dopar%
  #   2. backup the foreach env, run registerDoSEQ(), reinstate the origenal foreach setup.
  #
  # going with 1, for simplicity. Could be improved with rlang call2 or base call, then eval

  #I could also just registerDoSeq() if either parallel==FALSE or no backend exists.
  #Not so good, I don't want to clobber a backend just because is_parallel==fALSE
  if(!foreach::getDoParRegistered()) {
    #No parallel backend defined. Prevent warning
    foreach::registerDoSEQ()
  }

  if(is_parallel){
    ret <- foreach(k = rep(k, reps), .inorder = FALSE) %dopar% {
      cluster::clara(x = x, k = k, ...)
    }
    return(ret)
  } else {
    ret <- foreach(k = rep(k, reps), .inorder = FALSE) %do% {
      cluster::clara(x = x, k = k, ...)
    }
    return(ret)

  }
}


#' Find optimal mvpart to cluster mapping
#'
#' Given mvpart classifications and another clustering classification,
#' return a mapping that maximises the diagonal of the confusion matrix.
#'
#' The mvpart leaf node to cluster id mapping is 1 to 1. Any excess nodes or id's will be assigned NA.
#'
#' @param mvpart_map
#' @param clust_map
#'
#' @return list of: best confusion matrix, sum of diagonal of best confusion matrix, data.frame mapping mvpart_map entries to clust_map entries
#'
#' @examples
#'
#' if (requireNamespace("gradientForest", quietly = TRUE)) {
#' library(gradientForest) #required to attach extendedForest
#'
#'
#' data(CoMLsimulation)
#' preds <- colnames(Xsimulation)
#' specs <- colnames(Ysimulation)
#' f1 <- gradientForest(data.frame(Ysimulation,Xsimulation), preds, specs[1:6], ntree=10)
#' f2 <- gradientForest(data.frame(Ysimulation,Xsimulation), preds, specs[1:6+6], ntree=10)
#' f12 <- combinedGradientForest(west=f1,east=f2)
#'
#' env_grid <- f1$X
#' env_grid$site_id <- 1:nrow(env_grid)
#' gf_grid_sites <- data.frame(site_id = 1:nrow(f1$Y))
#' k_range <- 2:10
#' reps <- 1
#' is_parallel <- FALSE
#' pow <- 0.25
#' set.seed(1000)
#'
#' #Adjust clara fittings
#'
#' test3 <- gf_clust_f_ratio(gf = f1, gf_grid_sites = gf_grid_sites, env_grid = env_grid, k_range = k_range, reps =  reps, is_parallel = is_parallel, pow = pow,
#' clara_args = list(samples = 20, sampsize = 50, trace = 0, rngR = TRUE, pamLike = TRUE, correct.d = TRUE))
#'
#'
#' spatial_vars <- names(gf_grid_sites)
#' #Predict + compress#'
#'
#' env_trans <- do.call("gf_extrap_compress", c(list(gf = f1, env_grid = env_grid[!names(env_grid) %in% spatial_vars], pow = pow)))
#'
#' #reattach spatial data
#' env_trans <- cbind(env_grid[names(env_grid) %in% spatial_vars], env_trans)
#'
#' clara_clust <- merge(data.frame(env_trans[, spatial_vars, drop = FALSE], clust = test3$clust_list[[3]]$clustering), gf_grid_sites, by = "site_id")
#' node_clara <- merge(data.frame(site_id = gf_grid_sites, leaf = test3$mvpart), clara_clust, by = "site_id")
#' opt_test <- rphildyerphd:::opt_confusion(as.numeric(node_clara$leaf), node_clara$clust)
#' #Extract f-ratios per k into long form
#' set.seed(1000)
#' opt_test2 <- lapply(test3$clust_list, mvpart = test3$mvpart, function(clust, mvpart){
#'   clara_clust <- merge(data.frame(env_trans[, spatial_vars, drop = FALSE], clust = clust$clustering), gf_grid_sites, by = "site_id")
#'   node_clara <- merge(data.frame(site_id = gf_grid_sites, leaf = mvpart), clara_clust, by = "site_id")
#'   return(rphildyerphd:::opt_confusion(as.numeric(node_clara$leaf), node_clara$clust))
#' })
#'
#' scores <- sapply(opt_test2, function(x){x$score})
#' testthat::expect_equal(scores, c(29, 40, 48, 58, 71, 80, 91, 91, 89))
#'
#' }
#'
opt_confusion <- function(x, y){


  x <- factor(x)
  y <- factor(y)
  #confusion matrix
  ifelse(nlevels(x) <= nlevels(y), transpose <- FALSE, transpose <- TRUE)
  conf <- table(x, y)
  if(transpose) conf <- t(conf)
  mapping <- clue::solve_LSAP(conf, maximum = TRUE)
  conf_max <- conf[,mapping]
  if(transpose) conf_max <- t(conf_max)
  sum(diag(conf_max))


  if(transpose) {
    map_id <- data.frame(x_id = c(row.names(conf_max), levels(x)[!levels(x) %in% row.names(conf_max)]), y_id = c(colnames(conf_max), rep(NA, nlevels(x) - nlevels(y))))
  } else {
    map_id <- data.frame(x_id = c(row.names(conf_max), rep(NA, nlevels(y) - nlevels(x))), y_id = c(colnames(conf_max), levels(y)[!levels(y) %in% colnames(conf_max)]))
  }

  return(list(conf_max = conf_max, score = sum(diag(conf_max)), mapping = map_id))

}



