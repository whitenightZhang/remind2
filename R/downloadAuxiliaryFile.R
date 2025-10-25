#' Helper to download auxiliary file needed for reporting from RSE server
#' @param fname name of the file to be downloaded from RSE server folder with extra data
#' @author Falk Benke
downloadAuxiliaryFile <- function(fname) {

  dir <- tempdir()

  url <- paste0("https://rse.pik-potsdam.de/data/example/remind2_extraData/", fname)
  message("Downloading auxiliary file ", fname, " from ", url)
  f <- file.path(dir, fname)

  utils::download.file(url = url, destfile = f, mode = "wb", quiet = TRUE)

  return(f)

}
