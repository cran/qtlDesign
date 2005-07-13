"optspacing" <- function(cost,G=NULL,alpha=NULL,cross)
  {
    if(cross=="bc")
      optspacing.bc(cost,G,alpha)
    else if(cross=="f2")
      optspacing.f2(cost,G,alpha)
    else
      stop("Unknown cross ", cross, ".")
  }

"optspacing.bc" <-
function(cost,G=NULL,alpha=NULL)
{
  if(is.null(alpha))
      {
        tmp <- optim(par=c(50,0.5),fn=function(x,cost,G)
              {-info2cost.bc(x[2],cost,x[1],G)},
              method="L-BFGS-B",
              lower=c(.Machine$double.eps,.Machine$double.eps),
              upper=c(Inf,1),             
              G=G,cost=cost)$par
        tmp
        # list(d=tmp[1],alpha=tmp[2])
      }
    else
      {
        optim(par=50,fn=function(d,alpha,cost,G)
              {-info2cost.bc(alpha,cost,d,G)},
              method="L-BFGS-B",lower=.Machine$double.eps,
              upper=Inf,             
              G=G,alpha=alpha,cost=cost)$par
      }
  }

"optspacing.f2" <-
function(cost,G=NULL,alpha=NULL)
  {
    if(is.null(alpha))
      {
        tmp <- optim(par=c(50,0.5),fn=function(x,cost,G)
              {-info2cost.f2(x[2],cost,x[1],G)},
              method="L-BFGS-B",
              lower=c(.Machine$double.eps,.Machine$double.eps),
              upper=c(Inf,1),             
              G=G,cost=cost)$par
        tmp
        # list(d=tmp[1],alpha=tmp[2])
      }
    else
      {
        optim(par=50,fn=function(d,alpha,cost,G)
              {-info2cost.f2(alpha,cost,d,G)},
              method="L-BFGS-B",lower=.Machine$double.eps,
              upper=Inf,             
              G=G,alpha=alpha,cost=cost)$par
      }
  }


