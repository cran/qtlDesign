"info2cost.bc.null" <-
function(alpha,cost,d=0,G=NULL)
  {
    if((d==0) & is.null(G))
      {
        ans <- info.bc.null(alpha,theta=0)/(1+cost*alpha)
      }
    else
      {
        if((d==0)|is.null(G)|(G<=0))
          {
            stop("Cannot compute with given d and G.")
          }
        else
          {
            theta <- recomb(d)
            ans <- info.bc.null(alpha,theta)/(1+cost*alpha*G/d)
          }
      }
    ans
  }

