"deflate.f2" <-
function(theta)
  {
    t <- recomb(0.5*genetic.dist(theta))
    num <- 6*t^4 - 12*t^3 +10*t^2 - 4*t + 1
    den <- 8*t^4 - 16*t^3 +12*t^2 - 4*t + 1
    add.factor <- deflate.bc(theta)
    dom.factor <- (add.factor^2)*num/den
    list(add=add.factor,dom=dom.factor)
  }

