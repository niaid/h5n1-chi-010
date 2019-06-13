factor.date <- function(x, frmt="%m/%d/%y") {
  factor(x, 
         levels = unique(x[
           order(as.Date(as.character(x),format=frmt))]))
  
}
