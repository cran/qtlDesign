"power.f2" <-  function (n, model, thresh = 3, alpha = 1,
                       theta = 0, effective.n = FALSE) 
{
  # get info per individual
  iii <- info.f2.null(alpha, theta)
  # additive and dominance components
  a <- model[1]
  d <- model[2]
  # calculate non-centrality parameter
  ncp <- n * ( iii$add*a^2/2 + iii$dom*d^2/4 )
  m <- n*min(iii$add,iii$dom)
  # if effective sample size not big enough, stop
  if (m <= 30) {
    stop("Approximation not reliable as effective sample size < 30.")
  }
  # calculate threshold in chi-square scale
  T <- 2 * log(10) * thresh
  # calculate power
  pow <- 1 - pchisq(T, df = 2, ncp = ncp)
  # decide what to return depending on effective.n flag
  if (!effective.n) {
    return(pow)
  }
  else {
    return(list(power = pow, effective.n = m))
  }
}


"detectable.f2" <- function (n, model="add", power = 0.8, thresh = 3,
                           alpha = 1, theta = 0, delta=FALSE)
{
  # model
  if(model=="add")
    {
      a <- 1
      d <- 0
      model <- c(a,d)
    }
  else if(model=="dom")
    {
      a <- 1
      d <- 1
      model <- c(a,d)
    }
  else if( is.numeric(model) && (length(model)==2) )
    {
      model <- model
    }
  else
    {
      stop("Cannot understand model argument.")
    }

  # proportion of variance explained for given sample size,
  # power, threshold, selection fraction, and size of marker interval
  del <- uniroot(function(x) {
    power.f2(n, x*model, thresh, alpha, theta, effective.n = FALSE) - 
      power  }, interval = c(0, (1000/n)/(1+1000/n)))$root

    # decide what to return depending on delta flag
    if (!delta)
      {
        return(delta2prop.f2(del,model))
      }
    else
      {
        return(del*model)
      }
  
}

# from proportion of variance explained to a delta relative to a model
"prop2delta.f2" <- function(prop,model)
  {
    # additive and dominance components
    a <- model[1]
    d <- model[2]
    # get genetic variance
    gv <- prop/(1-prop)
    # convert to delta
    delta <- sqrt(gv/(d^2/4+a^2/2))
    delta <- delta*model
    delta
  }

"delta2prop.f2" <- function(delta,model)
  {
    # additive and dominance components
    a <- delta*model[1]
    d <- delta*model[2]
    gv <- (d^2/4+a^2/2)
    gv/(1+gv)
  }

# proportion of variance explained to genetic variance 
"prop2gv" <- function (prop) 
{
    prop/(1 - prop)
}
"gv2prop" <- function (gv) 
{
    gv/(1 + gv)
}

