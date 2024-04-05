#' Load all .Rdata or .rdata files in a directory
#'
#' This function scans a specified directory for all files with .Rdata or .rdata
#' extensions and loads them into the R environment. It provides error messages
#' for non-existent directories and directories with no relevant files.
#'
#' @param directory A character string of the directory path to scan for .Rdata or .rdata files.
#'
#' @examples
#' load_all_Rdata(here::here("inst","report","backend))
#'
#' @export
load_all_Rdata <- function(directory) {
  # Check if the directory exists
  if (!dir.exists(directory)) {
    stop("Directory does not exist: ", directory)
  }
  
  # List all .Rdata or .rdata files in the directory
  file_paths <- list.files(directory, pattern = "\\.Rdata$|\\.rdata$", full.names = TRUE)
  
  # Check if there are any files to load
  if (length(file_paths) == 0) {
    warning("No .Rdata or .rdata files found in the directory: ", directory)
    return(invisible(NULL))
  }
  
  # Load each file with error handling
  for (file_path in file_paths) {
    tryCatch({
      load(file_path, envir = .GlobalEnv)
    }, error = function(e) {
      warning("Failed to load file: ", file_path, "\nError: ", e$message)
    })
  }
}
