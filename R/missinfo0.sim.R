"missinfo0.sim" <-
function(delta,n,alpha)
{
y <- rnorm(n=n)
g <- rbinom(n=n,size=1,prob=0.5)
y <- y+(2*g-1)*delta

srt <- order(y)

y <- y[srt]
g <- g[srt]

m <- round(n*alpha/2)

if(m>=1)
	qq <- c(g[1:m], rep(0.5,n-2*m), g[(n-m+1):n])
else
	qq <- rep(0.5,n)

a <- qq*dnorm(y,mean=delta)
b <- (1-qq)*dnorm(y,mean=-delta)
qstar <- a/(a+b)

ans <- 4*y*y*qstar*(1-qstar)
list(mi=ans,y=y,g=g,delta=delta,n=n,alpha=alpha,qstar=qstar,
     ulim=y[n-m+1],llim=y[m])
}

