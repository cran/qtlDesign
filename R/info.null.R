"info.null" <-
function(alpha)
  {
    z <- -qnorm(alpha/2)
    2*( z*dnorm(z) + alpha/2 )
  }

"info.bc.null" <-
function(alpha,theta=0)
  {
    info.null(alpha)*deflate.bc(theta)
  }

"info.f2.null" <-
function(alpha,theta=0)
  {
    defl <- deflate.f2(theta)
    list( add=info.null(alpha)*defl$add, dom=info.null(alpha)*defl$dom )
  }

