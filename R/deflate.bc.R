"deflate.bc" <-
function(theta)
  {
    theta1 <- recomb(0.5*genetic.dist(theta))
    q <- ((1-theta1)^2)/(theta1^2 + (1-theta1)^2)
    A <- (1-4*q*(1-q))*(1-theta)
    A
  }

