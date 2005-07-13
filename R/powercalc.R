"powercalc" <- function(cross,n,effect,sigma2,
                        thresh=3,alpha=1,theta=0)
  {
    if(cross=="bc")
      powercalc.bc(n,effect,sigma2,thresh,alpha,theta)
    else if(cross=="f2")
      powercalc.f2(n,effectsigma2,thresh,alpha,theta)      
    else if(cross=="ri")
      powercalc.ri(n,effect,sigma2,thresh)      
    else
      stop("Unknown cross ", cross, ".")
  }

"detectable" <- function(cross,n,effect=NULL,
                         sigma2,power=0.8,thresh=3,alpha=1,theta=0) 
  {
    if(cross=="bc")
      detectable.bc(n,sigma2,power,thresh,alpha,theta)
    else if(cross=="f2")
      {
        if(is.null(effect))
          {
            detectable.f2(n,effect="add",sigma2,power,thresh,alpha,theta)
          }
        else
          {
            detectable.f2(n,effect=effect,sigma2,power,thresh,alpha,theta)
          }
      }
    else if(cross=="ri")
      detectable.ri(n,sigma2,power,thresh)
    else
      stop("Unknown cross ", cross, ".")
  }

"samplesize" <- function(cross,effect,sigma2,power=0.8,thresh=3,alpha=1,
                         theta=0) 
  {
    if(cross=="bc")
      samplesize.bc(effect,sigma2,power,thresh,alpha,theta)
    else if(cross=="f2")
      samplesize.f2(effect,sigma2,power,thresh,alpha,theta)
    else if(cross=="ri")
      samplesize.ri(effect,sigma2,power,thresh)
    else
      stop("Unknown cross ", cross, ".")
  }



"powercalc.bc" <- function(n,effect,sigma2,thresh=3,alpha=1,theta=0)
  {
    delta <- (effect/2)/sqrt(sigma2)
    # effective sample size
    m <- n * info.bc(alpha,theta)
    # non-centrality parameter
    ncp <- m*delta^2
    if( m<30 )
      {
        if( (alpha<1) | (theta>0) )
          {
            warning("Approximation not reliable as effective sample size < 30.")
          }
      }
    # threshold in 2*loglikelihood units
    T <- 2*log(10)*thresh
    # power using non-central chi-square
    pow <- 1-pchisq(T,df=1,ncp=ncp)
    return(pow)
  }


"detectable.bc" <- function (n, sigma2,
                             power = 0.8, thresh = 3, alpha = 1, theta = 0) 
{
  # proportion of variance explained for given sample size,
  # power, threshold, selection fraction, and size of marker interval
    effect <- uniroot(function(x) {
      powercalc.bc(n, x, sigma2, thresh, alpha, theta) -
        power}, interval = c(0,10*sqrt(sigma2/n)))$root
    effect
}

"samplesize.bc" <- function (effect, sigma2, power = 0.8, thresh = 3,
                             alpha = 1, theta = 0)
{

  # find an interval for the sample size to search
  # search in powers of 2
  p <- 0
  m <- 0
  while(p<power)
    {
      m <- m+1
      p <- suppressWarnings(powercalc.bc(2^m, effect, sigma2, thresh,
                                         alpha, theta))
    }

  # refine the number solving the power equation
    n <- uniroot(function(n) {
        powercalc.bc(n, effect, sigma2, thresh, alpha, theta) - power
    }, interval = c(2^(m-1), 2^m))$root

  # return the nearest largest integer
    return(ceiling(n))

}


"powercalc.f2" <-  function (n, effect, sigma2, thresh = 3, alpha = 1,
                             theta = 0) 
{
  # get info per individual
  iii <- info.f2(alpha, theta)
  # additive and dominance components
  a <- effect[1]/sqrt(sigma2)
  d <- effect[2]/sqrt(sigma2)
  # calculate non-centrality parameter
  ncp <- n * ( iii$add*a^2/2 + iii$dom*d^2/4 )
  m <- n*min(iii$add,iii$dom)
  # if effective sample size not big enough, stop
  if (m < 30)
    {
      if( (alpha<1) | (theta>0) )
        {
          warning("Approximation not reliable as effective sample size < 30.")
        }
    }
  # calculate threshold in chi-square scale
  T <- 2 * log(10) * thresh
  # calculate power
  pow <- 1 - pchisq(T, df = 2, ncp = ncp)
  pow

}


"detectable.f2" <- function (n, effect="add", sigma2, power = 0.8, thresh = 3,
                             alpha = 1, theta = 0)
{
  # effect
  if(effect=="add")
    {
      a <- 1
      d <- 0
      effect <- c(a,d)
    }
  else if(effect=="dom")
    {
      a <- 1
      d <- 1
      effect <- c(a,d)
    }
  else if( is.numeric(effect) && (length(effect)==2) )
    {
      effect <- effect
    }
  else
    {
      stop("Cannot understand effect argument.")
    }

  # proportion of variance explained for given sample size,
  # power, threshold, selection fraction, and size of marker interval
  del <- uniroot(function(x) {
    powercalc.f2(n, x*effect, sigma2, thresh, alpha, theta) - power  },
                 interval = c(0,10*sqrt(sigma2/n)))$root

    # decide what to return depending on delta flag
  return(del*effect)
}

"samplesize.f2" <- function (effect, sigma2, power = 0.8, thresh = 3,
                             alpha = 1, theta = 0)
{

  if(length(effect)!=2)
    {
      warning("Assuming additive effect.")
      effect <- c(effect[1],0)
    }
  # find an interval for the sample size to search
  # search in powers of 2
  p <- 0
  m <- 0
  while(p<power)
    {
      m <- m+1
      p <- suppressWarnings(powercalc.f2(2^m, effect, sigma2, thresh,
                                         alpha, theta))
    }

  # refine the number solving the power equation
    n <- uniroot(function(n) {
        powercalc.f2(n, effect, sigma2, thresh, alpha, theta) - 
          power }, interval = c(2^(m-1), 2^m))$root

  # return the nearest largest integer
    return(ceiling(n))

}



"powercalc.ri" <- function(n,effect,sigma2,thresh=3)
{
  powercalc.bc(n,effect*2,sigma2,thresh,alpha=1,theta=0)
}

"detectable.ri" <- function (n, sigma2, power = 0.8, thresh = 3) 
{
  # proportion of variance explained for given sample size,
  # power, threshold, selection fraction, and size of marker interval
    effect <- uniroot(function(x) {
      powercalc.ri(n, x, sigma2, thresh) -  power},
                     interval = c(0,10*sqrt(sigma2/n)))$root
    effect
}


"samplesize.ri" <- function (effect,sigma2,power = 0.8, thresh = 3)
{

   # find an interval for the sample size to search
  # search in powers of 2
  p <- 0
  m <- 0
  while(p<power)
    {
      m <- m+1
      p <- suppressWarnings(powercalc.ri(2^m, effect, sigma2, thresh))
    }

  # refine the number solving the power equation
    n <- uniroot(function(n) {
        powercalc.ri(n, effect, sigma2, thresh) - power
    }, interval = c(2^(m-1), 2^m))$root

  # return the nearest largest integer
    return(ceiling(n))
}



"error.var" <- function(cross,bio.var=1,gen.var=0,bio.reps=1)
  {
    # get the genetic variance multiplier 
    if(cross=="bc")
      CC <- 1/4
    else if (cross=="f2")
      CC <- 1/2
    else if (cross=="ri")
      CC <- 1
    else
      error("Cross type not recognized.")

    bio.var/bio.reps + CC*gen.var
    
  }



# function to convert the genotype means to additive and dominant
# effects segegating in a cross
"gmeans2effect" <- function(cross,means)
{
aa <- means[1]
ab <- means[2]
bb <- means[3]

if( cross == "f2" )
  {
    a <- (aa-bb)/2
    d <- ab - (bb+aa)/2
    effect <- c(a,d)
  }
else if( cross == "bc" )
  {
    effect <- c(aa-ab,ab-bb)
  }
else if (cross == "ri" )
  {
    effect <- (aa-bb)/2
  }

effect
}
