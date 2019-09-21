#' Align and join spatial data to a regular grid
#'
#' Takes two data.frames, each with lat and lon (or an equivalent spatial grid),
#' and a target grid.
#'
#' Every lat lon point is assigned to the nearest grid point.
#' Points that fall on the same grid cell are aggregated by fun.
#'
#' The target grid will align with the origin unless an offset is specified.
#'
#' The two data.frames are then merged, so that only rows appearing in both datasets are returned.
#'
#' This function is suitable for sparse datasets that are not necessarily on grid, and are in data.frames.
#' The raster package should be used if both datasets are on a regular grid.
#' sp and sf packages should be used if you have want more complex spatial manipulations.
#'
#' @param x first data.frame to merge
#' @param y second data.frame to merge
#' @param spatial_cols character vector length 2, columns in x and y that specify spatial position
#' @param res numeric. target grid resolution
#' @param offset numeric. target grid offset
#' @param fun function, applied to data that falls in the same grid cell
#' @param ... additional arguments to fun
#'
#' @return A data.frame. Each row has a unique lat and lon, duplicates are aggregated by fun. All other columns from x and y are included.
#'
#'
#' @export
#'
#' @examples
#'
#' set.seed(1000)
#'
#' target_res <- 1/10
#'
#' target_offset <- 0
#'
#' #Parameters are checked
#' invalid_data <- cluster::clara(cluster::xclara, 5)
#' valid_data <- data.frame(lat = 1:5, lon = 1:5, val = 1:5)
#' testthat::expect_error(align_merge_sp(invalid_data,  valid_data), "cannot coerce class ‘c(\"clara\", \"partition\")’ to a data.frame", fixed = TRUE)
#' testthat::expect_error(align_merge_sp(valid_data,  invalid_data), "cannot coerce class ‘c(\"clara\", \"partition\")’ to a data.frame", fixed = TRUE)
#' testthat::expect_error(align_merge_sp(valid_data[,2:3],  valid_data)  , "x does not have all of these name(s): 'lat', 'lon'",  fixed = TRUE)
#' testthat::expect_error(align_merge_sp(valid_data,  valid_data[,2:3])  , "y does not have all of these name(s): 'lat', 'lon'",  fixed = TRUE)
#' testthat::expect_error(align_merge_sp(valid_data,  valid_data, res = "a")  , "res is not a numeric or integer vector",  fixed = TRUE)
#' testthat::expect_error(align_merge_sp(valid_data,  valid_data, res = c(1,2))  , "length(res) not equal to 1",  fixed = TRUE)
#' testthat::expect_error(align_merge_sp(valid_data,  valid_data, offset = "a")  , "offset is not a numeric or integer vector",  fixed = TRUE)
#' testthat::expect_error(align_merge_sp(valid_data,  valid_data, offset = c(1,2))  , "length(offset) not equal to 1",  fixed = TRUE)
#' testthat::expect_error(align_merge_sp(valid_data,  valid_data, fun = notafun)  , "object 'notafun' not found",  fixed = TRUE)
#' testthat::expect_error(align_merge_sp(valid_data,  valid_data, fun = setdiff)  , "argument \"y\" is missing, with no default", fixed = TRUE) #fun is not an aggregator function
#'
#'
#'
#' #Align to grid
#' x <- data.frame(lat = c(0.1, .14, .26, 0.2, 0.35), lon = c(0.1, 0.2, 0.26, .14, 0.45), val_x = 1:5)
#'
#' testthat::expect_equal(align_merge_sp(x,x,spatial_cols =  c("lat", "lon"), res = target_res, offset = target_offset, fun = mean),
#'    data.frame(lat = c(0.1, .1, 0.2, .3, 0.3), lon = c(0.1, 0.2, .1, 0.3, 0.4), val_x.x = c(1,2,4,3,5), val_x.y =  c(1,2,4,3,5)))
#'
#' #Alternative spatial_cols
#' sp_test <- data.frame(latitude = c(0.1, .14, .26, 0.2, 0.35), longitude = c(0.1, 0.2, 0.26, .14, 0.45), val_x = 1:5)
#'
#' testthat::expect_equal(align_merge_sp(sp_test,sp_test,spatial_cols =  c("latitude", "longitude"), res = target_res, offset = target_offset, fun = mean),
#'    data.frame(latitude = c(0.1, .1, 0.2, .3, 0.3), longitude = c(0.1, 0.2, .1, 0.3, 0.4), val_x.x = c(1,2,4,3,5), val_x.y = c(1,2,4,3,5)))
#'
#'
#' #Aggregate function
#' y <- data.frame(lat = c(0.1, .14, .16, 0.2), lon = c(0.1, 0.14, 0.16, .2), val_y = 6:9)
#'
#' testthat::expect_equal(align_merge_sp(y,y,spatial_cols =  c("lat", "lon"), res = target_res, offset = target_offset, fun = mean),
#'    data.frame(lat = c(0.1, .2), lon = c(0.1, .2), val_y.x = c(6.5, 8.5), val_y.y =  c(6.5, 8.5)))
#'
#' #Merge data.frames
#'
#' z <- data.frame(lat = c(0.1, .14, .2, 0.26, 0.35), lon = c(0.1, 0.2, 0.14, .26, 0.45), val_z = 11:15)
#'
#' testthat::expect_equal(align_merge_sp(x,z,spatial_cols =  c("lat", "lon"), res = target_res, offset = target_offset, fun = mean),
#'    data.frame(lat = c(0.1, .1, 0.2, .3, 0.3), lon = c(0.1, 0.2, .1, 0.3, 0.4), val_x = c(1,2,4,3,5), val_z =  11:15))
#'
#' #MAggregate and merge
#'
#' y2 <- data.frame(lat = c(0.1, .14, .16, 0.2), lon = c(0.1, 0.14, 0.16, .2), val_y2 = 16:19)
#'
#' testthat::expect_equal(align_merge_sp(y,y2,spatial_cols =  c("lat", "lon"), res = target_res, offset = target_offset, fun = mean),
#'    data.frame(lat = c(0.1, .2), lon = c(0.1, .2), val_y = c(6.5, 8.5), val_y2 =  c(16.5, 18.5)))
#'
align_merge_sp <- function(x, y,
                           spatial_cols = c("lat", "lon"),
                           res = 1,
                           offset = 0,
                           fun = mean,
                           ...){


  x <- as.data.frame(x)
  y <- as.data.frame(y)

  assertthat::assert_that(assertthat::has_name(x, spatial_cols))
  assertthat::assert_that(assertthat::has_name(y, spatial_cols))

  assertthat::assert_that(is.numeric(res))
  assertthat::assert_that(assertthat::are_equal(length(res), 1))

  assertthat::assert_that(is.numeric(offset))
  assertthat::assert_that(assertthat::are_equal(length(offset), 1))


  x[,spatial_cols] <- round(x[,spatial_cols]/res + offset)
  x_agg <- aggregate(x, by = list(x[,spatial_cols[1]], x[,spatial_cols[2]]), fun, ..., simplify = TRUE)
  x_agg <- x_agg[,names(x)]
  y[,spatial_cols] <- round(y[,spatial_cols]/res + offset)
  y_agg <- aggregate(y, by = list(y[,spatial_cols[1]], y[,spatial_cols[2]]), fun, ..., simplify = TRUE)
  y_agg <- y_agg[,names(y)]
  xy <- merge(x_agg, y_agg, by = spatial_cols)
  xy[, spatial_cols] <- (xy[, spatial_cols] - offset)*res
  return(xy)
}
