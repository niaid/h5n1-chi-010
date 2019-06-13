gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

darken <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col/factor
	col[col>255]=255
    col <- rgb(t(col), maxColorValue=255)
    col
}


lighten <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col*factor
	col[col>255]=255
    col <- rgb(t(col), maxColorValue=255)
    col
}
