"info2cost" <- function(alpha,cost,d=0,G=NULL,cross)
  {
    if(cross=="bc")
      info2cost.bc(alpha,cost,d,G)
    else if(cross=="f2")
      info2cost.f2(alpha,cost,d,G)
    else
      stop("Unknown cross ", cross, ".")
  }


"info2cost.bc" <-
function(alpha,cost,d=0,G=NULL)
  {
    if((d==0) & is.null(G))
      {
        ans <- info.bc(alpha,theta=0)/(1+cost*alpha)
      }
    else
      {
        if((d==0)|is.null(G)|(G<=0))
          {
            stop("Cannot compute with given d and G.")
          }
        else
          {
            theta <- recomb(d/100)
            ans <- info.bc(alpha,theta)/(1+cost*alpha*G/d)
          }
      }
    ans
  }


"info2cost.f2" <-
function(alpha,cost,d=0,G=NULL)
  {
    if((d==0) & is.null(G))
      {
        ans <- info.f2(alpha,theta=0)/(1+cost*alpha)
      }
    else
      {
        if((d==0)|is.null(G)|(G<=0))
          {
            stop("Cannot compute with given d and G.")
          }
        else
          {
            theta <- recomb(d/100)
            ans <- info.f2(alpha,theta)$add/(1+cost*alpha*G/d)
          }
      }
    ans
  }
