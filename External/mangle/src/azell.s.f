c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      subroutine azell(rag,decg,l2p,raz,elp,azp,azg,elg,l2z)
      real*10 rag,decg,l2p,raz,elp,azp,azg,elg,l2z
c
c        parameters
      real*10 CIRCLE,PI,RADIAN
      parameter (CIRCLE = 360._10,
     *           PI = 3.1415926535897932384626_10,
     *           RADIAN = 180._10/PI)
c        local (automatic) variables
      integer iz
      real*10 cazm,cdecg,celg,celp,cl2m,cra
      real*10 sazm,sdecg,selg,selp,sl2m,sra
c *
c * Given transformations g <-> p between and z <-> p between spherical
c * frames, determine transformation g <-> z
c * (I've imagined g = galactic coordinates, p = celestial coordinates,
c * z = dome coordinates, but of course the transformation is quite
c * general).
c *
c  Input: rag = RA of NGP in deg [192.25]
c         decg = Dec of NGP in deg = latitude of NCP [27.4]
c         l2p = longitude of NCP in deg [123]
c         raz = RA of zenith in deg
c         elp = elevation of NCP in deg = Dec of zenith
c         azp = azimuth of NCP in deg
c Output: azg = azimuth of NGP in deg (0.le.az.lt.360)
c         elg = elevation NGP in deg (-90.le.el.le.90)
c         l2z = longitude of zenith in deg
c
c        sines and cosines of input angles
      sdecg=sin(decg/RADIAN)
      cdecg=cos(decg/RADIAN)
      selp=sin(elp/RADIAN)
      celp=cos(elp/RADIAN)
      sra=sin((rag-raz)/RADIAN)
      cra=cos((rag-raz)/RADIAN)
c        sine and cosine of elevation of NGP
      selg=cdecg*celp*cra+sdecg*selp
      if (selg.gt.1._10) then
        selg=1._10
      elseif (selg.lt.-1._10) then
        selg=-1._10
      endif
      celg=sqrt(1._10-selg*selg)
c        elevation of NGP in deg
      elg=asin(selg)*RADIAN
c        at NGP el +- 90 deg, set NGP az & zenith long consistently
      if (celg.eq.0._10) then
        azg=azp
        l2z=l2p+CIRCLE/2._10
      elseif (celg.ne.0._10) then
c        sine and cosine of azimuth of NGP relative to NCP azimuth
        sazm=-cdecg*sra/celg
        cazm=(sdecg*celp-cdecg*selp*cra)/celg
c        azimuth in deg
        azg=atan2(sazm,cazm)*RADIAN+azp
c        sine and cosine of longitude of zenith relative to NCP longitude
        sl2m=celp*sra/celg
        cl2m=(selp*cdecg-celp*sdecg*cra)/celg
c        longitude of zenith in deg
        l2z=atan2(sl2m,cl2m)*RADIAN+l2p
      endif
c        ensure azimuthal angles are in interval [0,360)
      iz=azg/CIRCLE
      if (azg.lt.0._10) iz=iz-1
      azg=azg-iz*CIRCLE
      iz=l2z/CIRCLE
      if (l2z.lt.0._10) iz=iz-1
      l2z=l2z-iz*CIRCLE
      return
      end
c
