"missinfo0" <-
function(delta,alpha)
{
  limit <- uniroot(pmixnorm,interval=c(0,delta+5),mean=c(-delta,delta),
               level=1-alpha/2)$root
  mi <- 2*integrate(fracmiss0,rel.tol=1e-7,lower=0,upper=limit,
                    delta=delta)$value
  list( mi=mi, lim=limit )
}

