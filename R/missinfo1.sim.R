"missinfo1.sim" <-
function(delta,n,alpha,theta)
{
y <- rnorm(n=n)
# simulate the QTL
g <- rbinom(n=n,size=1,prob=0.5)
# simulate the phenotype
y <- y+(2*g-1)*delta
# simulate a marker theta recombination fraction away
g1 <- rbinom(n=n,size=1,prob=theta)
g1 <- (g+g1) %% 2

srt <- order(y)

y <- y[srt]
g <- g[srt]
g1 <- g1[srt]

m <- round(n*alpha/2)

if(m>=1)
	qq <- c(g1[1:m]*(1-theta)+theta*(1-g1[1:m]),
                rep(0.5,n-2*m),
                g1[(n-m+1):n]*(1-theta)+theta*(1-g1[(n-m+1):n]))
else
	qq <- rep(0.5,n)

a <- qq*dnorm(y,mean=delta)
b <- (1-qq)*dnorm(y,mean=-delta)
qstar <- a/(a+b)

ans <- 4*y*y*qstar*(1-qstar)
list(mi=ans,y=y,g=g,delta=delta,n=n,alpha=alpha,theta=theta,qstar=qstar)
}

