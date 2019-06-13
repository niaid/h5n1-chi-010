repeated_elements_2vec <- function(x,y) {
  xi = which(c(diff(x)!=0,TRUE))
  yi = which(c(diff(y)!=0,TRUE))
  xi = sort(unique(c(xi,yi)))
  df = data.frame(value=x[xi], length=diff(c(0,xi)), iend=xi)
  return(df)
}
