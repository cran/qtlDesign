"fracmiss0" <-
function(y,delta)
{
  density <- 0.5 * ( dnorm(y,mean=delta) + dnorm(y,mean=-delta) )
  A <- exp(delta*y)
  B <- exp(-delta*y)
  genovar <- (A*B)/((A+B)^2)
  # print(c(density,genovar))
  4*y*y*density*genovar
}

