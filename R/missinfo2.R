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

