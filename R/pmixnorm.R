"pmixnorm" <-
function(x,mean=c(0,0),sd=c(1,1),mix.prop=0.5,level=0)
{
  ans <- mix.prop * pnorm(x,mean[1],sd[1]) +
    (1-mix.prop) * pnorm(x,mean[2],sd[2]) - level
  ans
}

