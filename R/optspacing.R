"optspacing" <- function(cost,G=NULL,sel.frac=NULL,cross)
  {
    if(cross=="bc")
      optspacing.bc(cost,G,sel.frac)
    else if(cross=="f2")
      optspacing.f2(cost,G,sel.frac)
    else
      stop("Unknown cross ", cross, ".")
  }

"optspacing.bc" <-
function(cost,G=NULL,sel.frac=NULL)
{
  if(is.null(sel.frac))
      {
        tmp <- optim(par=c(50,0.5),fn=function(x,cost,G)
              {-info2cost.bc(x[2],cost,x[1],G)},
              method="L-BFGS-B",
              lower=c(.Machine$double.eps,.Machine$double.eps),
              upper=c(Inf,1),             
              G=G,cost=cost)$par
        names(tmp) <- c("Marker spacing (cM)", "Selection fraction")
        tmp
      }
    else
      {
        tmp <- optim(par=50,fn=function(d,sel.frac,cost,G)
              {-info2cost.bc(sel.frac,cost,d,G)},
              method="L-BFGS-B",lower=.Machine$double.eps,
              upper=Inf,             
              G=G,sel.frac=sel.frac,cost=cost)$par
        tmp <- c(tmp,sel.frac)
        names(tmp) <- c("Marker spacing (cM)", "Selection fraction")
        tmp
      }
  }

"optspacing.f2" <-
function(cost,G=NULL,sel.frac=NULL)
  {
    if(is.null(sel.frac))
      {
        tmp <- optim(par=c(50,0.5),fn=function(x,cost,G)
              {-info2cost.f2(x[2],cost,x[1],G)},
              method="L-BFGS-B",
              lower=c(.Machine$double.eps,.Machine$double.eps),
              upper=c(Inf,1),             
              G=G,cost=cost)$par
        names(tmp) <- c("Marker spacing (cM)", "Selection fraction")
        tmp
        # list(d=tmp[1],sel.frac=tmp[2])
      }
    else
      {
        tmp <- optim(par=50,fn=function(d,sel.frac,cost,G)
              {-info2cost.f2(sel.frac,cost,d,G)},
              method="L-BFGS-B",lower=.Machine$double.eps,
              upper=Inf,             
              G=G,sel.frac=sel.frac,cost=cost)$par
        tmp <- c(tmp,sel.frac)
        names(tmp) <- c("Marker spacing (cM)", "Selection fraction")
        tmp
      }
  }


