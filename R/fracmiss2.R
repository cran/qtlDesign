"fracmiss2" <-
function(y,m1,m2,delta,theta1,theta2)
{
  # y = phenotype
  # m1 = left marker genotype
  # m2 = right marker genotype
  # delta = QTL effect; means are -delta and +delta
  # theta1 = recombination fraction to left marker
  # theta2 = recombination fraction to right marker

  lik0.pheno <- dnorm(y,mean=-delta)
  lik0.m1 <- theta1*m1 + (1-theta1)*(1-m1)
  lik0.m2 <- theta2*m2 + (1-theta2)*(1-m2)      
  lik0 <- lik0.pheno * lik0.m1 * lik0.m2 * 0.5

  lik1.pheno <- dnorm(y,mean=delta)
  lik1.m1 <- (1-theta1)*m1 + theta1*(1-m1)
  lik1.m2 <- (1-theta2)*m2 + theta2*(1-m2)      
  lik1 <- lik1.pheno * lik1.m1 * lik1.m2 * 0.5

  density <- lik0 + lik1
  genovar <- (lik0*lik1)/((lik0+lik1)^2)
  4*y*y*density*genovar
}

