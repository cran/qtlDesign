"missinfo" <-
function(delta,alpha,theta=NULL)
  {
    if(is.null(theta))
      {
        ans <- missinfo0(delta,alpha)
      }
    else
      if(length(theta)==1)
        {
          ans <- missinfo1(delta,alpha,theta)
        }
      else
        {
          ans <- missinfo2(delta,alpha,theta[1],theta[2])
        }
    mean(ans$mi)
  }

