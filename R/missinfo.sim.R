"missinfo.sim" <-
function(delta,n,alpha,theta=NULL)
  {
    if(is.null(theta))
      {
        ans <- missinfo0.sim(delta,n,alpha)
      }
    else
      if(length(theta)==1)
        {
          ans <- missinfo1.sim(delta,n,alpha,theta)
        }
      else
        {
          ans <- missinfo2.sim(delta,n,alpha,theta[1],theta[2])
        }
    mean(ans$mi)
  }

