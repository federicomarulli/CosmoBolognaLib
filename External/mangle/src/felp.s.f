c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      real*10 function felp(epoch)
      real*10 epoch
c
c        parameters
      include 'frames.par'
c        local (automatic) variables
      real*10 t
c *
c * Ecliptic latitude of NCP = Dec of ecliptic NP
c * as a function of epoch (e.g. 1950, 2000).
c *
c        RA & Dec epoch in centuries since 1900
      t=(epoch-1900._10)/100._10
c        ecliptic latitude of NCP = Dec of ecliptic NP
      felp=90._10-(E1+t*(E2+t*(E3+t*E4)))
      return
      end
c
