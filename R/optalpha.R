"optalpha" <- function(cost,d=0,G=NULL,cross)
  {
    if(cross=="bc")
      optalpha.bc(cost,d,G)
    else if(cross=="f2")
      optalpha.f2(cost,d,G)
    else
      stop("Unknown cross ", cross, ".")
  }

"optalpha.bc" <-
function(cost,d=0,G=NULL)
  {
    optimize(f=info2cost.bc,interval=c(0.0001,0.9999),max=TRUE,
             G=G,d=d,cost=cost)$maximum
  }

"optalpha.f2" <-
function(cost,d=0,G=NULL)
  {
    optimize(f=info2cost.f2,interval=c(0.0001,0.9999),max=TRUE,
             G=G,d=d,cost=cost)$maximum
  }

