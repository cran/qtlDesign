"nullinfo" <-
function(alpha)
  {
    z <- -qnorm(alpha/2)
    2*( z*dnorm(z) + alpha/2 )
  }

"info" <- function(alpha,theta=0,cross)
  {
    if(cross=="bc")
      info.bc(alpha,theta)
    else if(cross=="f2")
      info.f2(alpha,theta)
    else
      stop("Unknown cross ", cross, ".")
  }

"info.bc" <-
function(alpha,theta=0)
  {
    nullinfo(alpha)*deflate.bc(theta)
  }

"info.f2" <-
function(alpha,theta=0)
  {
    defl <- deflate.f2(theta)
    list( add=nullinfo(alpha)*defl$add, dom=nullinfo(alpha)*defl$dom )
  }

"deflate" <- function(theta,cross)
  {
    if(cross=="bc")
      deflate.bc(theta)
    else if(cross=="f2")
      deflate.f2(theta)
    else
      stop("Unknown cross ", cross, ".")
  }

"deflate.bc" <-
function(theta)
  {
    theta1 <- recomb(0.5*genetic.dist(theta))
    q <- ((1-theta1)^2)/(theta1^2 + (1-theta1)^2)
    A <- (1-4*q*(1-q))*(1-theta)
    A
  }

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


