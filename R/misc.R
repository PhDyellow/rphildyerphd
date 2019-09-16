#' load a knitr cache that is NOT lazy
#'
#' qwraps2::lazyload_cache_dir() works if you have set lazy = TRUE in knitr
#'
#' For very large files this can crash, so I may need to turn off lazy caching.
#'
#' @param path path to cache dir
#' @param cache_env an env that has already been created to capture the objects
#' @param ext cache file extension
load_cache_dir <- function(path, cache_env, ext = "RData"){
  files <- do.call(list.files, list(path = path, pattern = paste0("\\.", ext, "$"), full.names = full.names))
  lapply(files, load, cache_env)
}
