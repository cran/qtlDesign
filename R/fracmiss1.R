"fracmiss1" <-
function(y,m,delta,theta)
{
  if(m==0)
    {
      density <- 0.5 * ( theta * dnorm(y,mean=delta)
                        + (1-theta) * dnorm(y,mean=-delta) )
      A <- exp(delta*y) * theta
      B <- exp(-delta*y) * (1-theta)
    }
  if(m==1)
    {
      density <- 0.5 * ( (1-theta) * dnorm(y,mean=delta)
                        + theta * dnorm(y,mean=-delta) )
      A <- exp(delta*y) * (1-theta)
      B <- exp(-delta*y) * theta
    }

  genovar <- (A*B)/((A+B)^2)
  # print(c(density,genovar))
  4*y*y*density*genovar
}

