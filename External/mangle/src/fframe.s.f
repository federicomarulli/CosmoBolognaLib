c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      subroutine fframe(framei,azi,eli,framef,azf,elf)
      integer framei,framef
      real*10 azi,eli,azf,elf
c
c        parameters
      real*10 BEPOCH,JEPOCH
      parameter (BEPOCH=1950._10,JEPOCH=2000._10)
      include 'frames.par'
      include 'radian.par'
c        externals
      real*10 felp,SLA_epj2d
c        data variables
      logical init
c        saved local variables
      real*10 azg,elg,elp,l2z
      save azg,elg,elp,l2z
c        local (automatic) variables
      integer iaz
      real*10 date,dd,dec2k,dr,ra2k
c *
c * Transform azimuth (phi) and elevation (90-theta) in degrees
c * from one frame to another frame.
c *
c * Subroutines beginning SLA_ are from the STARLINK SLA library.
c * Please respect their copyright.
c *
c * WARNING:
c * The 2000 transformations are done consistently with the SLA
c * library, while the 1950 transformations are done consistently
c * using the transformations actually used by the catalogue makers
c * (notably, IRAS mask is in 1950 ecliptic coordinates).
c * The 1950 and 2000 transformations do not commute exactly,
c * but the residuals are < 1 arcsec.
c *
c * WARNING:
c * The SDSS eta, lambda system is unconventional
c * (latitude <-> longitude, and longitude goes backwards).
c * This subroutine returns conventional quantities
c *    azimuth = -eta
c *    elevation = lambda.
c *
c  Input: framei = initial frame, as defined in frames.par .
c         azi = azimuth (longitude) in degrees wrt framei .
c         eli = elevation (latitude) in degrees wrt framei .
c Output: framef = final frame, as defined in frames.par .
c         azf = azimuth (longitude) in degrees wrt framef
c               in interval [0,360).
c         elf = elevation (latitude) in degrees wrt framef
c               in interval [-90,90].
c
      data init /.true./
c
c--------initialize
      if (init) then
c        ecliptic latitude of North Celestial Pole (1950 FK4)
        elp=felp(BEPOCH)
c        angles to transform between ecliptic and galactic coordinates
        call azell(RAG,DECG,L2P,RAEZ,elp,EAZP,azg,elg,l2z)
        init=.false.
      endif

c--------initial and final frames identical
      if (framei.eq.framef) then
        azf=azi
        elf=eli

c--------1950 frames ->
      elseif (framei.eq.EQUATORIAL.or.framei.eq.ECLIPTIC) then
c........RA & Dec B1950.0 FK4 ->
        if (framei.eq.EQUATORIAL) then
          if (framef.eq.ECLIPTIC) then
            call azel(azi,eli,RAEZ,elp,EAZP,azf,elf)
          elseif (framef.eq.GALACTIC) then
            call azel(azi,eli,RAG,DECG,L2P,azf,elf)
          else
            call SLA_fk45z(azi/RADIAN,eli/RADIAN,BEPOCH,ra2k,dec2k)
            ra2k=ra2k*RADIAN
            dec2k=dec2k*RADIAN
          endif
c........Ecliptic 1950 ->
        elseif (framei.eq.ECLIPTIC) then
          if (framef.eq.EQUATORIAL) then
            call azel(azi,eli,EAZP,elp,RAEZ,azf,elf)
          elseif (framef.eq.GALACTIC) then
            call azel(azi,eli,azg,elg,l2z,azf,elf)
          else
            date=SLA_epj2d(BEPOCH)
            call SLA_ecleq(azi/RADIAN,eli/RADIAN,date,ra2k,dec2k)
            ra2k=ra2k*RADIAN
            dec2k=dec2k*RADIAN
          endif
        endif
c........-> 2000 frames
        if (framef.eq.EQUATORIAL2K) then
          azf=ra2k
          elf=dec2k
        elseif (framef.eq.ECLIPTIC2K) then
          date=SLA_epj2d(JEPOCH)
          call SLA_eqecl(ra2k/RADIAN,dec2k/RADIAN,date,azf,elf)
          azf=azf*RADIAN
          elf=elf*RADIAN
        elseif (framef.eq.SDSS) then
          call azel(ra2k,dec2k,RASDNP,DECSDNP,-ETANCP,azf,elf)
        endif

c--------2000 frames ->
      elseif (framei.eq.EQUATORIAL2K.or.framei.eq.ECLIPTIC2K
     *  .or.framei.eq.SDSS) then
c........initial frame -> RA & Dec 2000 FK5
        if (framei.eq.EQUATORIAL2K) then
          ra2k=azi
          dec2k=eli
        elseif (framei.eq.ECLIPTIC2K) then
          date=SLA_epj2d(JEPOCH)
          call SLA_ecleq(azi/RADIAN,eli/RADIAN,date,ra2k,dec2k)
          ra2k=ra2k*RADIAN
          dec2k=dec2k*RADIAN
        elseif (framei.eq.SDSS) then
          call azel(azi,eli,-ETANCP,DECSDNP,RASDNP,ra2k,dec2k)
        endif
c........RA & Dec 2000 FK5 -> final frame
        if (framef.eq.EQUATORIAL) then
          call SLA_fk54z(ra2k/RADIAN,dec2k/RADIAN,BEPOCH,azf,elf,dr,dd)
          azf=azf*RADIAN
          elf=elf*RADIAN
        elseif (framef.eq.ECLIPTIC) then
          date=SLA_epj2d(BEPOCH)
          call SLA_eqecl(ra2k/RADIAN,dec2k/RADIAN,date,azf,elf)
          azf=azf*RADIAN
          elf=elf*RADIAN
        elseif (framef.eq.GALACTIC) then
          call SLA_eqgal(ra2k/RADIAN,dec2k/RADIAN,azf,elf)
          azf=azf*RADIAN
          elf=elf*RADIAN
        elseif (framef.eq.EQUATORIAL2K) then
          azf=ra2k
          elf=dec2k
        elseif (framef.eq.ECLIPTIC2K) then
          date=SLA_epj2d(JEPOCH)
          call SLA_eqecl(ra2k/RADIAN,dec2k/RADIAN,date,azf,elf)
          azf=azf*RADIAN
          elf=elf*RADIAN
        elseif (framef.eq.SDSS) then
          call azel(ra2k,dec2k,RASDNP,DECSDNP,-ETANCP,azf,elf)
        endif

c--------Galactic frame ->
      elseif (framei.eq.GALACTIC) then
        if (framef.eq.EQUATORIAL) then
c         call SLA_ge50(azi/RADIAN,eli/RADIAN,azf,elf)
c         azf=azf*RADIAN
c         elf=elf*RADIAN
          call azel(azi,eli,L2P,DECG,RAG,azf,elf)
        elseif (framef.eq.ECLIPTIC) then
          call azel(azi,eli,l2z,elg,azg,azf,elf)
        else
          call SLA_galeq(azi/RADIAN,eli/RADIAN,ra2k,dec2k)
          if (framef.eq.EQUATORIAL2K) then
            azf=ra2k*RADIAN
            elf=dec2k*RADIAN
          elseif (framef.eq.ECLIPTIC2K) then
            date=SLA_epj2d(JEPOCH)
            call SLA_eqecl(ra2k,dec2k,date,azf,elf)
            azf=azf*RADIAN
            elf=elf*RADIAN
          elseif (framef.eq.SDSS) then
            ra2k=ra2k*RADIAN
            dec2k=dec2k*RADIAN
            call azel(ra2k,dec2k,RASDNP,DECSDNP,-ETANCP,azf,elf)
          endif
        endif

      endif
c--------put azf in interval [0,360)

      iaz=azf/360._10
      if (azf.lt.0._10) iaz=iaz-1
      azf=azf-iaz*360._10

      return
      end
c
