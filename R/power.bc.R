"power.bc" <- function(n,prop,thresh=3,alpha=1,theta=0,effective.n=FALSE)
  {
    # convert to deviation from overall mean
    delta <- prop2delta.bc(prop)
    # effective sample size
    m <- n * info.bc.null(alpha,theta)
    # non-centrality parameter
    ncp <- m*delta^2
    if( m<=30 )
      {
        stop("Approximation not reliable as effective sample size < 30.")
      }
    # threshold in 2*loglikelihood units
    T <- 2*log(10)*thresh
    # power using non-central chi-square
    pow <- 1-pchisq(T,df=1,ncp=ncp)
    if(!effective.n)
     {
       return(pow)
     }
    else
      {
        return(list(power=pow,effective.n=m))
      }
  }

"detectable.bc" <- function (n, power = 0.8, thresh = 3,
                           alpha = 1, theta = 0, delta = FALSE) 
{
  # proportion of variance explained for given sample size,
  # power, threshold, selection fraction, and size of marker interval
    prop <- uniroot(function(x) {
        power.bc(n, x, thresh, alpha, theta, effective.n = FALSE) - 
            power
    }, interval = c(0, (1000/n)/(1+1000/n)))$root
    
  # decide what to return depending on delta flag
    if (!delta) {
        return(prop)
    }
    else {
        return(prop2delta.bc(prop))
    }
}


"prop2delta.bc" <- function(prop)
  {
    sqrt(1/(1-prop)-1)
  }

"delta2prop.bc" <- function(delta)
  {
    delta^2/(1+delta^2)
  }
