"fracmiss" <-
function(y,delta,m=NULL,theta=NULL)
  {
    if(is.null(m)&is.null(theta))
      {
        ans <- fracmiss0(y,delta)
      }
    else
      {
        if( (length(m)==1) & (length(theta)==1) )
          {
            ans <- fracmiss1(y,m,delta,theta)
          }
        else
          {
            if( (length(m)==2) & (length(theta)==2) )
              {
                ans <- fracmiss2(y,m[1],m[2],delta,theta[1],theta[2])
              }
            else
              {
                stop("Incompatible m and theta lengths.")
              }
          }
      }
    ans
  }

