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

