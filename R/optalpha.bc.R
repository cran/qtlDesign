"optalpha.bc" <-
function(cost,d=0,G=NULL)
  {
    optimize(f=info2cost.bc.null,interval=c(0.0001,0.9999),max=TRUE,
             G=G,d=d,cost=cost)$maximum
  }

