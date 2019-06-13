hours2days <- function(time, format.h="d%03d_%02dh", format.d="d%03d") {
  t = time/24
  d = as.integer(t)
  h = time - (d*24)
  days.by.hours = unique(d[h>0])
  string = ifelse(d %in% days.by.hours,
    sprintf(format.h, d, h), sprintf(format.d, d) )
}
