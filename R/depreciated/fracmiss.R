"fracmiss" <-
function(y,delta,m=NULL,theta=NULL)
  {
    if(is.null(m)&is.null(theta))
      {
        ans <- fracmiss0(y,delta)
      }
    else
      {
        if( (length(m)==1) & (length(theta)==1) )
          {
            ans <- fracmiss1(y,m,delta,theta)
          }
        else
          {
            if( (length(m)==2) & (length(theta)==2) )
              {
                ans <- fracmiss2(y,m[1],m[2],delta,theta[1],theta[2])
              }
            else
              {
                stop("Incompatible m and theta lengths.")
              }
          }
      }
    ans
  }


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

"fracmiss2" <-
function(y,m1,m2,delta,theta1,theta2)
{
  # y = phenotype
  # m1 = left marker genotype
  # m2 = right marker genotype
  # delta = QTL effect; means are -delta and +delta
  # theta1 = recombination fraction to left marker
  # theta2 = recombination fraction to right marker

  lik0.pheno <- dnorm(y,mean=-delta)
  lik0.m1 <- theta1*m1 + (1-theta1)*(1-m1)
  lik0.m2 <- theta2*m2 + (1-theta2)*(1-m2)      
  lik0 <- lik0.pheno * lik0.m1 * lik0.m2 * 0.5

  lik1.pheno <- dnorm(y,mean=delta)
  lik1.m1 <- (1-theta1)*m1 + theta1*(1-m1)
  lik1.m2 <- (1-theta2)*m2 + theta2*(1-m2)      
  lik1 <- lik1.pheno * lik1.m1 * lik1.m2 * 0.5

  density <- lik0 + lik1
  genovar <- (lik0*lik1)/((lik0+lik1)^2)
  4*y*y*density*genovar
}

