#' martys_format_p: Formatting p-values
#'
#' This function takes a p-value and formats it based on specific conditions.
#' The function aims to make p-values easier to interpret and present.
#'
#' @param p A p-value, which should ideally be of type "double".
#'
#' @return Returns the p-value formatted as a character string according
#' to specific rules:
#' - If p is NA, returns "-".
#' - If p is not of type "double", returns p as is.
#' - If p < 0.0001, returns "<0.0001".
#' - If p > 0.1, returns p rounded to 3 decimal places.
#' - Otherwise, returns p rounded to 4 decimal places.
#'
#' @examples
#' martys_format_p(0.00009) # should return "<0.0001"
#' martys_format_p(0.5)     # should return "0.5"
#' martys_format_p(NA)      # should return "-"
#' 
#' @export
martys_format_p <- function (p) {
  # Check for NA values first
  if (is.na(p)) {
    return("-")
  }
  
  # Check for type "double"
  if (typeof(p) != "double") {
    return(as.character(p))
  }
  
  # Handle values less than 0.0001
  if (p < 0.0001) {
    return("<0.0001")
  } 
  
  else if (p == 0.0001) {
    return("0.0001")
  }
  
  # Handle values greater than 0.1
  else if (p > 0.1) {
    return(format(round(p, digits = 3), digits = 3))
  }
  
  # Handle all other cases
  else {
    return(format(round(p, digits = 4)))
  }
}
