"missinfo" <-
function(delta,alpha,theta=NULL)
  {
    if(is.null(theta))
      {
        ans <- missinfo0(delta,alpha)
      }
    else
      if(length(theta)==1)
        {
          ans <- missinfo1(delta,alpha,theta)
        }
      else
        {
          ans <- missinfo2(delta,alpha,theta[1],theta[2])
        }
    mean(ans$mi)
  }

"missinfo0" <-
function(delta,alpha)
{
  limit <- uniroot(pmixnorm,interval=c(0,delta+5),mean=c(-delta,delta),
               level=1-alpha/2)$root
  mi <- 2*integrate(fracmiss0,rel.tol=1e-7,lower=0,upper=limit,
                    delta=delta)$value
  list( mi=mi, lim=limit )
}

"missinfo1" <-
function(delta,alpha,theta)
{
  limit <- uniroot(pmixnorm,interval=c(0,delta+5),mean=c(-delta,delta),
               level=1-alpha/2)$root
  mi.mid <- 2*integrate(fracmiss0,sub=10000,rel.tol=1e-7,lower=0,upper=limit,
                         delta=delta)$value
  mi0.upp <- 2*integrate(fracmiss1,sub=10000,rel.tol=1e-7,lower=limit,
                         upper=delta+10,
                         m=0,delta=delta,theta=theta)$value
  mi1.upp <- 2*integrate(fracmiss1,sub=10000,rel.tol=1e-7,lower=limit,
                         upper=delta+10,
                         m=1,delta=delta,theta=theta)$value

  list( mi=mi.mid+mi0.upp+mi1.upp, lim=limit )
}

"missinfo2" <-
function(delta,alpha,theta1,theta2)
{
  limit <- uniroot(pmixnorm,interval=c(0,delta+5),mean=c(-delta,delta),
               level=1-alpha/2)$root
  mi.mid <- 2*integrate(fracmiss0,sub=10000,rel.tol=1e-7,lower=0,upper=limit,
                         delta=delta)$value

  mi00.upp <- 2*integrate(fracmiss2,sub=10000,rel.tol=1e-7,lower=limit,
                          upper=delta+10,m1=0,m2=0,delta=delta,
                          theta1=theta1,theta2=theta2)$value

  mi01.upp <- 2*integrate(fracmiss2,sub=10000,rel.tol=1e-7,lower=limit,
                          upper=delta+10,m1=0,m2=1,delta=delta,
                          theta1=theta1,theta2=theta2)$value

  mi10.upp <- 2*integrate(fracmiss2,sub=10000,rel.tol=1e-7,lower=limit,
                          upper=delta+10,m1=1,m2=0,delta=delta,
                          theta1=theta1,theta2=theta2)$value

  mi11.upp <- 2*integrate(fracmiss2,sub=10000,rel.tol=1e-7,lower=limit,
                          upper=delta+10,m1=1,m2=1,delta=delta,
                          theta1=theta1,theta2=theta2)$value


  list( mi=mi.mid+mi00.upp+mi01.upp+mi10.upp+mi11.upp, lim=limit )
}

