"missinfo2.sim" <-
function(delta,n,alpha,theta1,theta2)
{
y <- rnorm(n=n)
# simulate the QTL
g <- rbinom(n=n,size=1,prob=0.5)
# simulate the phenotype
y <- y+(2*g-1)*delta
# simulate a marker theta1 recombination fraction away
g1 <- rbinom(n=n,size=1,prob=theta1)
g1 <- (g+g1) %% 2
# simulate a marker theta2 recombination fraction away
g2 <- rbinom(n=n,size=1,prob=theta2)
g2 <- (g+g2) %% 2

# find order of phenotypes
srt <- order(y)
# sort the phenotypes, QTL genotypes, and flanking marker genotypes
y <- y[srt]
g <- g[srt]
g1 <- g1[srt]
g2 <- g2[srt]

# calculate prior probabilities of QTL=1 given the flanking markers
# prob QTL=1
qq.1 <- g1*(1-theta1) + (1-g1)*theta1 + g2*(1-theta2) + (1-g2)*theta2
# prob QTL=0
qq.0 <- (1-g1)*(1-theta1) + g1*theta1 + (1-g2)*(1-theta2) + g2*theta2
# prob QTL=1 given flanking markers
qq <- qq.1/(qq.0+qq.1)

# number genotyped at each flank
m <- round(n*alpha/2)

# if number genotyped is greater than 2, get the prior probs
# else replace it by the equilibrium distribution
if(m>=1)
	qq <- c(qq[1:m], rep(0.5,n-2*m), qq[(n-m+1):n])
else
	qq <- rep(0.5,n)

# likelihood for y given QTL=1
a <- qq*dnorm(y,mean=delta)
# likelihood for y given QTL=0
b <- (1-qq)*dnorm(y,mean=-delta)
# posterior prob of QTL=1
qstar <- a/(a+b)

ans <- 4*y*y*qstar*(1-qstar)
list(mi=ans,y=y,g=g,delta=delta,n=n,alpha=alpha,theta=c(theta1,theta2),
     qstar=qstar)
}

