c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      real*10 function twodf230k(ra,dec)
      real*10 ra,dec
c
c        parameters
      include 'pi.par'
      real*10 RADIAN
      parameter (RADIAN=180._10/PI)
c        data variables
C     logical init
C     real magcut
c        local (automatic) variables
      real rar,decr,compl,maglim
c *
c * Completeness of 30 June 2001 2dF 100k galaxy survey.
c * This is an interface to
c * the 2dFGRS mask software by Peder Norberg and Shaun Cole.
c 
c * WARNING:
c * There is a hard-wired magnitude cut.
c * Comment it out if you don't want it.
c *
c  Input: ra, dec = RA & Dec (B1950) in degrees.
c Output: twodf230k: 0 = empty to 1 = complete.
c
C     data init /.true./
c This magnitude cut gives the maximum number of survivors
C     data magcut /19.27/

c convert from real*10 degrees to real*4 radians
      rar=ra
      decr=dec
      if (rar.lt.0.) rar=rar+360.
      if (rar.lt.0.) rar=rar+360.
      rar=rar/RADIAN
      decr=decr/RADIAN
c
c the 2dFGRS mask software by Peder Norberg and Shaun Cole
      call in_2df_mask2(rar,decr,compl,maglim)

c compl = -1. means the position is outside the main 2dF boundary
      if (compl.eq.-1.) then
        twodf230k=0._10

c compl = -2. means the position is inside a drill hole
      elseif (compl.eq.-2.) then
        twodf230k=0._10

c discard fields with limiting magnitude brighter than magcut
C COMMENT THIS OUT IF YOU DON'T WANT IT
C     elseif (maglim.lt.magcut) then
C       if (init) then
C         write (*,'(" WARNING from twodf230k: USING HARD-WIRED MAGNITUDE CU
C    *T OF ",f5.2)') magcut
C         init=.false.
C       endif
C       twodf230k=0._10

c standard
      else

c completeness
        twodf230k=compl

c magnitude limit
c       twodf230k=maglim

      endif

c completeness: round to 8 decimal places
      twodf230k=dble(nint(twodf230k*1.e8_10))/1.e8_10

c     print *,ra,dec,compl,maglim,twodf230k
c
      return
      end
c
c=======================================================================
c
c The remaining code below is (with minor modifications)
c the mask_compl.f subroutine provided by the 2dFGRS team
c in the 30 June 2001 public release of data at
c http://www.mso.anu.edu.au/2dFGRS/Public/Release/index.html
c
c The minor modifications are options to read/write
c formatted as well as fortran unformatted data.
c
c The documentation at
c http://www.mso.anu.edu.au/2dFGRS/Public/Release/Masks/index.html
c states the following:
c
c 2. Copyright and use of the mask codes
c    As a user of the mask codes, we request that you observe the following guidelines:
c    The code provided here has been written by several members of the 2dFGRS team.
c    Peder Norberg and Shaun Cole (University of Durham, UK) wrote the
c    software used to create the survey masks.
c    If you use any part of this 2dFGRS software, please acknowledge
c    "the 2dFGRS mask software by Peder Norberg and Shaun Cole".
c    This code is supplied as-is (i.e. we do not support any modifications
c    to the file mask_compl.f containing the subroutines) and without guarantees
c    - use it with caution and report any bugs.
c    Read all of this documentation before trying to use the code!
c
c-----------------------------------------------------------------------------
c     For a given position in ra & dec, this subroutine tells if this
c     position is inside the 2df_mask or not, by returning the actual
c     completeness. It returns in the same time the corresponding
c     magnitude limit.

      subroutine in_2df_mask2(ra,dec,compl,maglim)
*************************************variables********************************
      implicit none

      real ra,dec,compl,maglim

      integer np_x,np_y
      parameter (np_x=2*1800,np_y=2*600)
      integer np_xx,np_yy,ifirst,ix,iy,np,sgn
      real rpix(np_x,np_y),magpix(np_x,np_y),mupix(np_x,np_y)
      real x,y,z,rx,ry
      character name_mask*20,last_reg*3,reg*3
      save ifirst,np_xx,np_yy,np,rpix,magpix,last_reg
      data ifirst/1/
******************************************************************************

c     Everytime check which region one points to: either ngp, sgp or ran
      call which_reg2(ra,dec,reg)

      if (reg.eq.'ran') then
         call in_2df_ran2(ra,dec,compl,maglim)
         return
      endif

c     On the first call, the program reads the mask
      if ((ifirst.eq.1).or.(last_reg.ne.reg)) then
         ifirst=0
         last_reg=reg
         name_mask='mask.'//reg//'.dat'
         call read_mask2(rpix,np_xx,np_yy,name_mask,np_x,np_y,np)
         name_mask='maglim.'//reg//'.dat'
         call read_mask2(magpix,np_xx,np_yy,name_mask,np_x,np_y,np)
         name_mask='mumask.'//reg//'.dat'
         call read_mask2(mupix,np_xx,np_yy,name_mask,np_x,np_y,np)
c         write(0,*) 'Read masks for ',reg,' region!'
      endif

      call radec_xyz2(ra,dec,x,y,z)
      call eq_2dfrx2(np,x,y,z,rx,ry,sgn)
      call rx_ix_map2(rx,ry,ix,iy)

      if ((ix.ge.1).and.(ix.le.np_xx).and.
     :    (iy.ge.1).and.(iy.le.np_yy)) then
         compl=rpix(ix,iy)
         maglim=magpix(ix,iy)
      else
         compl= -1.0  !default value (i.e. outside 2dF boundary)
         maglim=-2.0  !default value (i.e. outside 2dF boundary)
      endif

      return
      end

c-----------------------------------------------------------------------------
c     For a given position in ra & dec, this subroutine tells if this
c     position is inside one of the used random fields or not, by returning
c     the actual completeness of the field.

      subroutine in_2df_ran2(ra,dec,compl,maglim)
*************************************variables********************************
      implicit none

      real ra,dec,compl,maglim

      integer nmax,NP,NPX,NPY
      real theta_min
      parameter (nmax=98,NP=2400,NPX=100,NPY=100,theta_min=1.014)
      logical in
      integer i,ifirst,nb_field,io
      real costheta,costheta_r,costheta_min,xc(nmax),yc(nmax),zc(nmax)
     &,    comp(nmax),r,x,y,z,magpix(nmax,NPX,NPY),mupix(nmax,NPX,NPY)
      character name_mask*20
      data ifirst/1/
      save ifirst,nb_field,xc,yc,zc,comp,costheta_min,magpix
******************************************************************************

c     On the first call, reads the file containing the random fields
      if (ifirst.eq.1) then
         ifirst=0
         name_mask='mask.ran.dat'
         call ran2df_cen_used2(nb_field,xc,yc,zc,comp,name_mask)
         name_mask='maglim.ran.dat'
         call readranmask2(name_mask,NPX,NPY,nmax,magpix)
         name_mask='mumask.ran.dat'
         call readranmask2(name_mask,NPX,NPY,nmax,mupix)
         costheta_min=cos(theta_min*atan(1.)/45.)
c         write(0,*) 'Read mask for ran region!'
      endif

      i=0
      in=.false.
      call radec_xyz2(ra,dec,x,y,z)
      r=sqrt(x**2+y**2+z**2)
      costheta_r=costheta_min*r
      do while ((.not.in).and.(i.lt.nb_field))
         i=i+1
         costheta=(x*xc(i)+y*yc(i)+z*zc(i))
         if (costheta.gt.costheta_r) then   ! We use here the fact the random
            in=.true.                       ! fields doesn't overlapp with
         endif                              ! each other...
      enddo

      if (in) then
         call holes_2df_xyz2(x,y,z,io)
c         compl=comp(i)*real(io)       ! io = 0 if in a hole; 1 otherwise
         if (io.eq.0) then
            compl=-2.           ! default value in holes
            maglim=-2.          ! default value when not specified
         else
            compl=comp(i)
            call get_maglimran2(NP,NPX,NPY,x,y,z,xc(i),yc(i),zc(i)
     &,                        maglim,nmax,magpix,i)
         endif


      else
         compl=-1.              ! default value (ie. outside 2dF boundary)
         maglim=-2.             ! default value (when not specified)
      endif

      return
      end

c----------------------------------------------------------------------------
c     Given x,y,z position and centre of the corresponding random field,
c     returns the magnitude limit at that position
      subroutine get_maglimran2(NP,NPX,NPY,x,y,z,xc,yc,zc,maglim
     &,                        nmax,magpix,ifield)

      implicit none

      integer NP,NPX,NPY,nmax,ifield
      real x,y,z,xc,yc,zc,maglim,magpix(nmax,NPX,NPY)

      intrinsic max,min

      integer sgn,ixc_min,iyc_min,ix,iy
      real rx,ry

c      integer ifirst,max_ix,max_iy,min_ix,min_iy,inum
c      data ifirst/-1/
c      save ifirst,max_ix,max_iy,min_ix,min_iy,inum

c      if (ifirst.eq.-1) then
c         ifirst=0
c         inum=0
c         max_ix=0
c         min_ix=NPX
c         max_iy=0
c         min_iy=NPY
c      endif

c      ifirst=ifirst+1
      call eq_2dfrx2(NP,xc,yc,zc,rx,ry,sgn)
      call bound_map_ran2(NPX,NPY,rx,ry,ixc_min,iyc_min)

      call eq_2dfrx2(NP,x,y,z,rx,ry,sgn)
      call rx_ix_map2(rx,ry,ix,iy)

      ix=min(max(1,ix-ixc_min),NPX)
      iy=min(max(1,iy-iyc_min),NPY)
c      ix = ix-ixc_min
c      iy = iy-iyc_min
c      if ((ix.ge.1).and.(ix.le.NPX).and.(iy.ge.1).and.(iy.le.NPY)) then
         maglim=magpix(ifield,ix,iy)
c      else
c         write(*,*) 'ix=',ix,'iy=',iy,ifield
c         max_ix=max(ix,max_ix)
c         min_ix=min(ix,min_ix)
c         max_iy=max(iy,max_iy)
c         min_iy=min(iy,min_iy)
c         maglim=1.0
c         inum=inum+1
c      endif
c
c      if (ifirst.ge.57012) then
c         write(*,*) ifirst,inum,max_ix,max_iy,min_ix,min_iy
c      endif

      return
      end

c----------------------------------------------------------------------------
c     Given ra & dec determines which region to look at
      subroutine which_reg2(ra,dec,reg)
*************************************variables********************************
      implicit none

      real ra,dec
      character reg*(*)

      real EPS
      parameter (EPS=1.e-5)
      integer i,ifirst,nb_strip,nb_stripsgp,nb_stripngp
      real ra_min(4),ra_max(4),dec_min(4),dec_max(4)
      data ifirst/1/
      save ra_min,ra_max,dec_min,dec_max,nb_strip,nb_stripsgp
     &,    nb_stripngp,ifirst
******************************************************************************

      if (ifirst.eq.1) then
         ifirst=0
         nb_stripsgp=3
         nb_stripngp=4
         nb_strip=5
c        SGP strips
         ra_min(1)= 5.711590    ! 21h49m
         ra_max(1)= 0.911935    !  3h29m
         dec_min(1)= -0.479966  ! -27.5 deg
         dec_max(1)= -0.392699  ! -22.5 deg
         ra_min(2)= 5.670138    ! 21h39.5m
         ra_max(2)= 0.975203    !  3h43.5m
         dec_min(2)= -0.567232  ! -32.5 deg
         dec_max(2)= -0.479966  ! -27.5 deg
         ra_min(3)= 5.707227    ! 21h48m
         ra_max(3)= 0.890118    !  3h24m
         dec_min(3)= -0.654498  ! -37.5 deg
         dec_max(3)= -0.567232  ! -32.5 deg
c        NGP strips
         ra_min(4)= 2.574361    !  9h50m
         ra_max(4)= 3.883358    ! 14h50m
         dec_min(4)= -0.130900  !  -7.5 deg
         dec_max(4)= +0.043633  !  +2.5 deg
c        Make boundaries more secure
         do i=1,4
            ra_min(i)=ra_min(i)-EPS
            ra_max(i)=ra_max(i)+EPS
            dec_min(i)=dec_min(i)-EPS
            dec_max(i)=dec_max(i)+EPS
         enddo
      endif

      i=0
      do while (i.lt.nb_stripsgp)
         i=i+1
         if (((ra_min(i).le.ra).or.(ra_max(i).ge.ra)).and.
     &        (dec_min(i).le.dec).and.(dec_max(i).ge.dec)) then
            i=nb_strip
            reg='sgp'
         endif
      enddo

      do while (i.lt.nb_stripngp)
         i=i+1
         if ((ra_min(i).le.ra).and.(ra_max(i).ge.ra).and.
     &        (dec_min(i).le.dec).and.(dec_max(i).ge.dec)) then
            i=nb_strip
            reg='ngp'
         endif
      enddo

      if (i.ne.nb_strip) reg='ran'

      return
      end

c----------------------------------------------------------------------------
c  Subroutine to read the 2dF mask
c
      subroutine  read_mask2(mask,np_xx,np_yy,name_mask,np_x,np_y,np)
*************************************variables********************************
      implicit none

      integer np_x,np_y,np_xx,np_yy,np
      real mask(np_x,np_y)
      character name_mask*(*)

      integer ix,iy

      integer access,lnblnk
      character*1 go
      character*128 dat
      logical np_set2,ok
      logical ex
      character*3 re
******************************************************************************

      dat=name_mask(1:lnblnk(name_mask))//'.fmt'

      ok=.false.

c try unformatted
      inquire(FILE=name_mask, EXIST=ex, READ=re)
      if (ex.and.(re.eq.'YES')) then
c      if (access(name_mask,' r').eq.0) then
        open(11,file=name_mask,status='old',form='unformatted')
        rewind(11)
        read(11,end=200,err=200) np_xx,np_yy
        ok=np_set2(np_xx,np_yy,np,0)
        if (.not.ok) then
          print *,'try formatted version of file instead ...'
          close(11)
          goto 200
        endif
        ok=.false.
        if ((np_xx.gt.np_x).or.(np_yy.gt.np_y)) then
           write(0,*) 'Dimension of the mask too big!'
           write(0,*) ' Mask: NPX= ',np_xx,' NPY= ',np_yy
           print *,'try formatted version of file instead ...'
           close(11)
           goto 200
        endif
        do ix=1,np_xx
          read(11,end=200,err=200) (mask(ix,iy),iy=1,np_yy)
        enddo
        print *,np_xx,' x',np_yy,' values read from ',
     *    name_mask(1:lnblnk(name_mask))
        close(11)
        ok=.true.

c COMMENT THIS OUT IF YOU DON'T WANT IT
C        inquire(FILE=dat, EXIST=ex)
C        if (.not.ex) then
Cc       if (access(dat,' ').ne.0) then
C         print *,'write FORMATTED values? [CR,n=no, y=yes]'
C         read (*,'(a1)',end=150,err=150) go
C         if (go.eq.'y'.or.go.eq.'Y') then
C           open(11,file=dat)
C           write(11,'(2i8)') np_xx,np_yy
C           do ix=1,np_xx
C             write (11,'(5g16.8)') (mask(ix,iy),iy=1,np_yy)
C           enddo
C           print *,np_xx,' x',np_yy,' values written to ',
C    *        dat(1:lnblnk(dat))
C           close(11)
C         endif
C 150     continue
C       endif

      endif

c try formatted
  200 if (.not.ok) then
         inquire(FILE=dat, EXIST=ex)
         if (.not.ex) goto 300
c        if (access(dat,' ').ne.0) goto 300
        open(11,file=dat)
        rewind(11)
        read(11,'(2i8)') np_xx,np_yy
        ok=np_set2(np_xx,np_yy,np,0)
        if (.not.ok) goto 300
        ok=.false.
        if ((np_xx.gt.np_x).or.(np_yy.gt.np_y)) then
          write(0,*) 'Dimension of the mask too big!'
          write(0,*) ' Mask: NPX= ',np_xx,' NPY= ',np_yy
          close(11)
          goto 300
        endif
        do ix=1,np_xx
          read (11,'(5g16.8)') (mask(ix,iy),iy=1,np_yy)
        enddo
        print *,np_xx,' x',np_yy,' values read from ',
     *    dat(1:lnblnk(dat))
        close(11)
        ok=.true.

c COMMENT THIS OUT IF YOU DON'T WANT IT
        inquire(FILE=name_mask, EXIST=ex)
        if (.not.ex) then
c       if (access(name_mask,' ').ne.0) then
         print *,'write UNFORMATTED values? [CR,n=no, y=yes]'
         read (*,'(a1)',end=250,err=250) go
         if (go.eq.'y'.or.go.eq.'Y') then
           open(11,file=name_mask,form='unformatted')
           write(11) np_xx,np_yy
           do ix=1,np_xx
             write(11) (mask(ix,iy),iy=1,np_yy)
           enddo
           print *,np_xx,' x',np_yy,' values written to from ',
     *        name_mask(1:lnblnk(name_mask))
           close(11)
         endif
 250     continue
       endif

      endif
c
      return
c
c error
  300 print *,'failed to read unformatted data from ',
     *  name_mask(1:lnblnk(name_mask))
      print *,'or formatted data from ',
     *  dat(1:lnblnk(dat))
      inquire(FILE=name_mask, EXIST=ex, READ=re)
      if (ex.and.(re.eq.'YES')) then
c      if (access(name_mask,' r').eq.0) then
        print *,'you may have a problem with endianness;'
        print *,'please read HELP.unformatted in the mangle directory'
      endif
      stop
c
      end

c----------------------------------------------------------------------------
      subroutine readranmask2(namemask,NPX,NPY,nbf,rpix)

c     This subroutine reads the mask stored in the file name_mask.
c     N.B.: This subroutine contains the same stuff as mask_2df, but has the
c     advantage of giving back also sgn!

      implicit none

      integer nbf,NPX,NPY
      real rpix(nbf,NPX,NPY)
      character namemask*(*)

      integer i,j,k,npxx,npyy,nbff

      integer access,lnblnk
      character*1 go
      character*128 dat
      logical ok
      logical ex
      character*3 re

      ok=.false.

      dat=namemask(1:lnblnk(namemask))//'.fmt'

c try unformatted
      inquire(FILE=namemask, EXIST=ex, READ=re)
      if (ex.and.(re.eq.'YES')) then
c     if (access(namemask,' r').eq.0) then
         open(unit=14,file=namemask,status='old',form='unformatted')
         rewind(14)
         read(14,end=200,err=200) nbff,npxx,npyy
         if ((npxx.ne.NPX).or.(npyy.ne.NPY).or.(nbff.ne.nbf)) then
            write(*,*) 'Wrong dimensions for ran mask. Should have:'
     &           ,             npxx,npyy,nbff
            print *,'try formatted version of file instead ...'
            close(14)
            goto 200
         endif
         do k=1,nbff
            do i=1,npxx
               read(14,end=200,err=200) (rpix(k,i,j),j=1,npyy)
            enddo
         enddo
         print *,nbff,' x',npxx,' x',npyy,
     *        ' values read from ',namemask(1:lnblnk(namemask))
         ok=.true.
         close(14)

c     COMMENT THIS OUT IF YOU DON'T WANT IT
C     inquire(FILE=dat, EXIST=ex)
C     if (.not.ex) then
Cc       if (access(dat,' ').ne.0) then
C         print *,'write FORMATTED values? [CR,n=no, y=yes]'
C         read (*,'(a1)',end=150,err=150) go
C         if (go.eq.'y'.or.go.eq.'Y') then
C           open(14,file=dat)
C           write(14,'(3i8)') nbff,npxx,npyy
C           do k=1,nbff
C             do i=1,npxx
C               write (14,'(5g16.8)') (rpix(k,i,j),j=1,npyy)
C             enddo
C           enddo
C           print *,nbff,' x',npxx,' x',npyy,
C    *        ' values written to ',dat(1:lnblnk(dat))
C           close(14)
C         endif
C 150     continue
C       endif

      endif

c try formatted
 200  if (.not.ok) then
         inquire(FILE=dat, EXIST=ex)
         if (.not.ex) goto 300
c        if (access(dat,' ').ne.0) goto 300
        open(unit=14,file=dat)
        rewind(14)
        read(14,'(3i8)') nbff,npxx,npyy
        if ((npxx.ne.NPX).or.(npyy.ne.NPY).or.(nbff.ne.nbf)) then
          write(*,*) 'Wrong dimensions for ran mask. Should have:'
     &  ,             npxx,npyy,nbff
          goto 300
        endif
        do k=1,nbff
          do i=1,npxx
            read (14,'(5g16.8)') (rpix(k,i,j),j=1,npyy)
          enddo
        enddo
        print *,nbff,' x',npxx,' x',npyy,
     *    ' values read from ',dat(1:lnblnk(dat))
        close(14)
        ok=.true.

c COMMENT THIS OUT IF YOU DON'T WANT IT
         inquire(FILE=namemask, EXIST=ex)
         if (.not.ex) then
c       if (access(namemask,' ').ne.0) then
         print *,'write UNFORMATTED values? [CR,n=no, y=yes]'
         read (*,'(a1)',end=250,err=250) go
         if (go.eq.'y'.or.go.eq.'Y') then
           open(14,file=namemask,form='unformatted')
           write(14) nbff,npxx,npyy
           do k=1,nbff
             do i=1,npxx
               write (14) (rpix(k,i,j),j=1,npyy)
             enddo
           enddo
           print *,nbff,' x',npxx,' x',npyy,
     *        ' values written to ',namemask(1:lnblnk(namemask))
           close(14)
         endif
 250     continue
       endif

      endif
c
      return
c
c error
  300 print *,'failed to read unformatted data from ',
     *  namemask(1:lnblnk(namemask))
      print *,'or formatted data from ',
     *  dat(1:lnblnk(dat))
      inquire(FILE=namemask, EXIST=ex, READ=re)
      if (ex.and.(re.eq.'YES')) then
c      if (access(namemask,' r').eq.0) then
        print *,'you may have a problem with endianness;'
        print *,'please read HELP.unformatted in the mangle directory'
      endif
      stop
c
      end

c-----------------------------------------------------------------------------
c     This subroutine reads the used 2dF random field centres from the file
c     name, which contains also the overall completeness of the field (with
c     respect to the underlying density field). Read the information written
c     by write_ran2df_cen.
      subroutine ran2df_cen_used2(nb_field,xc,yc,zc,comp,name)
*************************************variables********************************
      implicit none

      integer nb_field
      real xc(*),yc(*),zc(*),comp(*)
      character name*(*)

      integer i

      integer access,lnblnk
      character*1 go
      character*128 dat
      logical ok
      logical ex
      character*3 re

******************************************************************************

      dat=name(1:lnblnk(name))//'.fmt'

      ok=.false.

c try unformatted
      inquire(FILE=name, EXIST=ex, READ=re)
      if (ex.and.(re.eq.'YES')) then
c      if (access(name,' r').eq.0) then
        open(unit=11,file=name,status='old',form='unformatted')
        rewind(11)
        read(11,end=200,err=200) nb_field
        read(11,end=200,err=200) (xc(i),i=1,nb_field)
        read(11,end=200,err=200) (yc(i),i=1,nb_field)
        read(11,end=200,err=200) (zc(i),i=1,nb_field)
        read(11,end=200,err=200) (comp(i),i=1,nb_field)
        print *,nb_field,' x y z w values read from ',
     *    name(1:lnblnk(name))
        close(11)
        ok=.true.

c COMMENT THIS OUT IF YOU DON'T WANT IT
C         inquire(FILE=dat, EXIST=ex)
C         if (.not.ex) then
Cc       if (access(dat,' ').ne.0) then
C         print *,'write FORMATTED values? [CR,n=no, y=yes]'
C         read (*,'(a1)',end=150,err=150) go
C         if (go.eq.'y'.or.go.eq.'Y') then
C           open(11,file=dat)
C           write(11,'(i8," fields")') nb_field
C           do i=1,nb_field
C             write (11,'(4g16.8)') xc(i),yc(i),zc(i),comp(i)
C           enddo
C           print *,nb_field,' lines written to ',dat(1:lnblnk(dat))
C           close(11)
C         endif
C 150     continue
C       endif

      endif

c try formatted
  200 if (.not.ok) then
         inquire(FILE=dat, EXIST=ex)
         if (.not.ex) goto 300
c        if (access(dat,' ').ne.0) goto 300
        open(11,file=dat)
        rewind(11)
        read(11,'(i8)') nb_field
        do i=1,nb_field
          read (11,'(4g16.8)') xc(i),yc(i),zc(i),comp(i)
        enddo
        print *,nb_field,' lines read from ',dat(1:lnblnk(dat))
        close(11)
        ok=.true.

c COMMENT THIS OUT IF YOU DON'T WANT IT
         inquire(FILE=name, EXIST=ex)
         if (.not.ex) then
c       if (access(name,' ').ne.0) then
         print *,'write UNFORMATTED values? [CR,n=no, y=yes]'
         read (*,'(a1)',end=250,err=250) go
         if (go.eq.'y'.or.go.eq.'Y') then
           open(unit=11,file=name,form='unformatted')
           write(11) nb_field
           write(11) (xc(i),i=1,nb_field)
           write(11) (yc(i),i=1,nb_field)
           write(11) (zc(i),i=1,nb_field)
           write(11) (comp(i),i=1,nb_field)
           close(11)
           print *,nb_field,' x y z w values written to ',
     *        name(1:lnblnk(name))
         endif
 250     continue
       endif

      endif

      return
c
c error
  300 print *,'failed to read unformatted data from ',
     *  name(1:lnblnk(name))
      print *,'or formatted data from ',
     *  dat(1:lnblnk(dat))
      inquire(FILE=name, EXIST=ex, READ=re)
      if (ex.and.(re.eq.'YES')) then
c      if (access(name,' r').eq.0) then
        print *,'you may have a problem with endianness;'
        print *,'please read HELP.unformatted in the mangle directory'
      endif
      stop
c
      end

c-----------------------------------------------------------------------------
c     Convert an equatorial cartesian coordinate to the rx,ry pixel
c     coordinate used for the 2df mask. This subroutine is valid for both
c     sgp and ngp elements.
c
c
      subroutine  eq_2dfrx2(NP,x,y,z,rx,ry,sgn)
*************************************variables********************************
      implicit none
      intrinsic abs
      real x,y,z,r,PI,EPS
      parameter (EPS=3.e-13)
      integer ifirst,sgn,NP,ix_min,iy_min
      real racen,sc,deccen,rx,ry
      double precision costheta
      real xc,yc,zc,phi,dec_gc,ra_gc,xgc,ygc,zgc,sintheta_2
     & ,xgp,ygp,zgp,xx,yy,zz,phi0
      save ifirst,sc,xc,yc,zc,xgc,ygc,zgc,xgp,ygp,zgp,phi0,PI !ix_min,iy_min
      data ifirst/1/
******************************************************************************

c     On the first call define the values which define the
c     centre, orientation and scale of the transformation
      if (ifirst.eq.1) then
         ifirst=0
         PI= atan(1.)*4.
c        This is the direction of the Galactic Centre as adopted
c        by Steve Maddox
         racen  = 12.3*PI/180.0
         deccen = -27.5*PI/180.0
         sc = 120.0*PI/180.0
         xc=cos(deccen)*cos(racen)
         yc=cos(deccen)*sin(racen)
         zc=sin(deccen)
c        This is a direction perpendicular to above but otherwise
c        arbitrary provided that the phi0 has then been correctly
c        set to be the offset between this arbitrary direction and
c        that used for the projection of the APM SGP region
         ra_gc   = -94.40593*PI/180.0
         dec_gc  = -28.90771*PI/180.0
         phi0=236.97694*PI/180.0
         xgc=cos(dec_gc)*cos(ra_gc)
         ygc=cos(dec_gc)*sin(ra_gc)
         zgc=sin(dec_gc)
c        Generate the mutually perpendicular unit vector via
c        the cross product
         xgp = yc*zgc-zc*ygc
         ygp = zc*xgc-xc*zgc
         zgp = xc*ygc-yc*xgc
      endif

c     Section utilized on each call
c     Find components in the new 3D Cartesian system
      r=sqrt(x**2+y**2+z**2)
      zz=(x*xc+y*yc+z*zc)
      xx=(x*xgc+y*ygc+z*zgc)
      yy=(x*xgp+y*ygp+z*zgp)

c     Compute corresponding spherical polar angles
      costheta= abs(dble(zz)/dble(r))
      sgn=int(sign(1.,zz/r)) ! tells if x,y,z is in the ngp or sgp region
      phi=phi0-atan2(yy,xx)
      if (phi.gt.2.0*PI) phi=phi-2.0*PI

c     Apply the Zenithal Equal Area Projection
      if (sngl(dble(1.0)-costheta).ge.EPS) then
         sintheta_2=sngl(dsqrt(dble(0.5)*(dble(1.0)-costheta)))
      else
         sintheta_2=0.
      endif
      r=2.*sintheta_2 !If sc changes, change this like:r=sintheta_2/sin(sc/4.)

c     convert to x,y coordinate
      xx=-r*sin(phi)*real(NP)*0.5
      yy= r*cos(phi)*real(NP)*0.5

      call bound_map_2df2(NP,sgn,ix_min,iy_min)

      rx=xx+0.5*real(NP)-1.*real(ix_min)
      ry=-yy+0.5*real(NP)-1.*real(iy_min)

      return
      end

c-----------------------------------------------------------------------------
c     Convert the pixel postion rx,ry to ra and dec.
c     This subroutine is the inverse of eq_2dfrx2. Works for both ngp and sgp,
c     as long as the sgn is given (respectively by -1 and 1).

      subroutine inv_eq_2dfrx2(NP,rx,ry,sgn,x,y,z)
************************************variables*********************************
      implicit none
      intrinsic abs
      real x,y,z,r,PI,EPS,EPS2,a(3,3),PI_hf,PI_3hf
      parameter (EPS=3.e-13)
      integer ifirst,sgn,NP,iy_min,ix_min
      real racen,sc,deccen,costheta,rx,ry,delta_phi
     & ,xc,yc,zc,phi,dec_gc,ra_gc,xgc,ygc,zgc,sintheta_2
     & ,xgp,ygp,zgp,xx,yy,zz,phi0
      save ifirst,sc,phi0,PI,EPS2,a,PI_hf,PI_3hf,xc,yc,zc   !ix_min,iy_min
      data ifirst/1/
******************************************************************************

c     On the first call define the values which define the
c     centre, orientation and scale of the transformation
      if (ifirst.eq.1) then
         ifirst=0
         EPS2=sqrt(0.5*EPS)
         PI= atan(1.)*4.
         PI_hf=PI/2.
         PI_3hf=3.*PI/2.
c        This is the direction of the Galactic Centre as adopted
c        by Steve Maddox
         racen  = 12.3*PI/180.0
         deccen = -27.5*PI/180.0
         sc = 120.0*PI/180.0
         xc=cos(deccen)*cos(racen)
         yc=cos(deccen)*sin(racen)
         zc=sin(deccen)
c        This is a direction perpendicular to above but otherwise
c        arbitrary provided that the phi0 has then been correctly
c        set to be the offset between this arbitrary direction and
c        that used for the projection of the APM SGP region
         ra_gc   = -94.40593*PI/180.0
         dec_gc  = -28.90771*PI/180.0
         phi0=236.97694*PI/180.0
         xgc=cos(dec_gc)*cos(ra_gc)
         ygc=cos(dec_gc)*sin(ra_gc)
         zgc=sin(dec_gc)
c        Generate the mutually perpendicular unit vector via
c        the cross product
         xgp = yc*zgc-zc*ygc
         ygp = zc*xgc-xc*zgc
         zgp = xc*ygc-yc*xgc
c        Define Matrix a which is the inverse (transpose) of the base matrix
         a(1,1)=xc
         a(1,2)=xgc
         a(1,3)=xgp
         a(2,1)=yc
         a(2,2)=ygc
         a(2,3)=ygp
         a(3,1)=zc
         a(3,2)=zgc
         a(3,3)=zgp
      endif

c     Section utilized on each call
      call bound_map_2df2(NP,sgn,ix_min,iy_min)

      xx=rx+1.*real(ix_min)-0.5*real(NP)
      yy=-ry-1.*real(iy_min)+0.5*real(NP)

      if ((xx.eq.0.).and.(yy.eq.0.)) then
         x=xc*real(sgn)
         y=yc*real(sgn)
         z=zc*real(sgn)
         return
      endif

      phi=atan2(-xx,yy)+2.*PI
      if (xx.eq.0.) then
         r=yy/(cos(phi)*real(NP)*0.5)
      else if (yy.eq.0.) then
         r=-xx/(sin(phi)*real(NP)*0.5)
      else if (    ((sgn.eq. 1).and.(abs(xx).gt.50.))
     &         .or.((sgn.eq.-1).and.(abs(yy).lt.50.))) then
         r=-xx/(sin(phi)*real(NP)*0.5)
      else
        r=yy/(cos(phi)*real(NP)*0.5)
      endif

      sintheta_2=r*0.5!If sc changes, change like this:sintheta_2=r*sin(sc/4.)
      if (sintheta_2.gt.EPS2) then
         costheta=sngl((dble(1.)-dble(2.*sintheta_2**2)))*real(sgn)
      else
         costheta=1.*real(sgn)
      endif

      zz=costheta
      if (zz.gt.1.)  zz= 1. ! rounding error corrections
      if (zz.lt.-1.) zz=-1.

      delta_phi=phi0-phi
      if (delta_phi.lt.0.) delta_phi=delta_phi+2.*PI
      if ((delta_phi.le.PI_hf).or.(delta_phi.ge.PI_3hf)) then ! xx is positive
         xx=sqrt((1.-zz**2)/(1.+tan(delta_phi)**2))
         yy=xx*tan(delta_phi)
      else                                                    ! xx is negative
         xx=-sqrt((1.-zz**2)/(1.+tan(delta_phi)**2))
         yy=xx*tan(delta_phi)
      endif

      x=a(1,1)*zz+a(1,2)*xx+a(1,3)*yy
      y=a(2,1)*zz+a(2,2)*xx+a(2,3)*yy
      z=a(3,1)*zz+a(3,2)*xx+a(3,3)*yy

      return
      end

c-----------------------------------------------------------------------------
c     This subroutine contains the boundaries of the sgp and ngp pixel map.
c     It returns the ix_min and iy_min which are the offsets of the two maps.

      subroutine bound_map_2df2(NP,sgn,ix_min,iy_min)
******************************************************************************
      implicit none

      integer NP,sgn,ix_min,iy_min
******************************************************************************

      if (sgn.eq.1) then        ! use the sgp map
         ix_min=NP*15/100
         iy_min=NP*44/100
      else if (sgn.eq.-1) then  ! use the ngp map
         ix_min=NP*7/100
         iy_min=NP*60/100
      endif

      return
      end

c----------------------------------------------------------------------------
c This subroutine gives the right boundaries for each random field mask, ie.
c calculates the offset needed such that the pixel center of each random
c field, given by (rxc,ryc), is located at (NPX_ran/2, NPY_ran/2) .
c
      subroutine bound_map_ran2(NPX_ran,NPY_ran,rxc,ryc,ixc_min,iyc_min)

      implicit none

      integer NPX_ran,NPY_ran,ixc_min,iyc_min
      real rxc,ryc

      integer ixc,iyc


      call rx_ix_map2(rxc,ryc,ixc,iyc)
      ixc_min=ixc-int(NPX_ran/2)
      iyc_min=iyc-int(NPY_ran/2)

      return
      end

c-----------------------------------------------------------------------------
c     This subroutine converts real pixel coordinates rx,ry to integer
c     ix,iy pixel coordinates.

      subroutine rx_ix_map2(rx,ry,ix,iy)
******************************************************************************
      implicit none

      integer ix,iy
      real rx,ry
******************************************************************************

      ix=int(rx+0.5)
      iy=int(ry+0.5)

      return
      end

c-----------------------------------------------------------------------------
c     This subroutine converts integer pixel coordinates ix,iy to real
c     rx,ry pixel coordinates which corresponds to the center of the pixel.

      subroutine ix_rx_map2(ix,iy,rx,ry)
******************************************************************************
      implicit none

      integer ix,iy
      real rx,ry
******************************************************************************

      rx=real(ix)
      ry=real(iy)

      return
      end

c-----------------------------------------------------------------------------
c     This subroutine transforms x,y & z cartesian coordinates to ra & dec.

      subroutine xyz_radec2(x,y,z,ra,dec)
******************************************************************************
      implicit none

      real x,y,z,ra,dec

      real PI2
      parameter (PI2=6.2831853072)
      real zz
******************************************************************************

      zz=sqrt(x**2+y**2+z**2)
      dec= asin(z/zz)
      if (y.lt.0) then
         ra= atan2(y,x)+PI2
      else
         ra= atan2(y,x)
      endif

      return
      end

c-----------------------------------------------------------------------------
c This subroutine transforms ra & dec into x,y & z cartesian coordinates.

      subroutine radec_xyz2(ra,dec,x,y,z)
******************************************************************************
      implicit none

      real x,y,z,ra,dec
******************************************************************************

      x=cos(dec)*cos(ra)
      y=cos(dec)*sin(ra)
      z=sin(dec)

      return
      end

c-----------------------------------------------------------------------------
c     This subroutine sets and checks that NP and NPX_NPY are correct and that
c     they correspond to the value used for NPX and NPY, the dimensions of the
c     pixel grid.
c     N.B.: NPX_NPY is also accepted if equal to 0 (default value if not
c     used!)

      logical function np_set2(NPX,NPY,NP,NPX_NPY)
******************************************************************************
      implicit none

      integer NPX,NPY,NP,NPX_NPY

      integer num
******************************************************************************

      np_set2=.true.

      NP=NPX*4/3
      if (NP.ne.NPY*4) then
         write(0,*) 'Distortion of the pixel grid: NPX= ',NPX,' NPY='
     &,              NPY,' do not obey NPX=3*NPY!'
         np_set2=.false.
      endif

c     We check also that NPX_NPY is NPX*NPY or NPX_NPY = 0
      num=NPX*NPY
      if ((num.ne.NPX_NPY).and.(NPX_NPY.ne.0)) then
         write(0,*) 'Wrong dimension on NPX_NPY; should be ',num
         np_set2=.false.
      endif

      return
      end

c-----------------------------------------------------------------------------
c     This subroutine holes_2df_xyz2 tells if (x,y,z) is in a hole or not,
c     depending on the value of io (io=1 not in a hole; io=0 in a hole). The
c     value of sgn tells if it is a position in the ngp or sgp/ran.

      subroutine holes_2df_xyz2(x,y,z,io)
*************************************variables********************************
      implicit none

      integer io
      real x,y,z

      intrinsic abs

      integer jf,ifirst,in,sgn
      real ra,dec,racen,deccen,xc,yc,zc,PI,r,zz
      save ifirst,xc,yc,zc

      real pival,degs_rads,ramins_rads,radfac,tpi
      parameter (pival=3.1415926535897)
      parameter (degs_rads=pival/180.,ramins_rads=pival/720.
     &,          radfac=pival/10800.,tpi=2.*pival)
      data ifirst/1/
******************************************************************************

      if (ifirst.eq.1) then
         ifirst=0
         PI= atan(1.)*4.
c        This is the direction of the Galactic Centre as adopted
c        by Steve Maddox
         racen  = 12.3*PI/180.0
         deccen = -27.5*PI/180.0
         xc=cos(deccen)*cos(racen)
         yc=cos(deccen)*sin(racen)
         zc=sin(deccen)
      endif

      in=0
      r=sqrt(x**2+y**2+z**2)
      zz=(x*xc+y*yc+z*zc)
      sgn=int(sign(1.,zz/r))

      call initialisation2(sgn)
      call xyz_radec2(x,y,z,ra,dec)

      if ((ra.gt.pival).and.(sgn.eq.1))then
         ra=ra-tpi
      endif

      call fnumber2(ra,dec,jf)

      if (sgn.eq.1) then
         call test_sgp_holes2(ra,dec,jf,in)
      else
         call test_ngp_holes2(ra,dec,jf,in)
      endif

      io=abs(in-1) ! N.B.: in = 1 if in a hole and 0 otherwise(opposite to io)

      return
      end

c----------------------------------------------------------------------------
      subroutine initialisation2(sgn)
******************************************************************************
      implicit none

      integer sgn

      integer sgn_ifirst
      save sgn_ifirst
      data sgn_ifirst/0/
******************************************************************************

      if (sgn.ne.sgn_ifirst) then
         sgn_ifirst=sgn
         if (sgn.eq.1) then
            call sgp_mask_init2
         else
            call ngp_mask_init2
         endif
      endif

      return
      end

c-----------------------------------------------------------------------------
c     Given ra,dec in radians returns corresponding field no.

      subroutine fnumber2(rar,decr,ifield)
******************************************************************************
      implicit none

      intrinsic abs

      integer ifield
      real decr,rar

      integer iygrid,ixgrid,i
      real ra,dec,raminus,yspace,ygrid,xspace,xgrid,xoff,yoff

      integer*4 field0h(22)
      real*4 ragaps(22), decbands(22)

      real pival,degs_rads,ramins_rads,radfac,tpi
      parameter (pival=3.1415926535897)
      parameter (degs_rads=pival/180.,ramins_rads=pival/720.
     &,          radfac=pival/10800.,tpi=2.*pival)

      data (ragaps(i),i=1,22)/ 0.,144.,90.,66.,52.,44.,38.,33.,30.,28.,
     +                 26.,24.,23.,22.,21.,20.,20.,20.,20.,20.,20.,20./
      data (decbands(i),i=1,22)/-90.,-85.,-80.,-75.,-70.,-65.,-60.,-55.,
     +   -50.,-45.,-40.,-35.,-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15./
      data (field0h(i),i=1,22)/1,2,12,28,50,78,111,149,193,241,293,349,
     +   409,472,538,607,679,751,823,895,967,1039 /
******************************************************************************

c     find the nearest grid centre
      ra = rar*180./pival
      dec = decr*180./pival

      if (dec.gt.15.) then
         ifield = 0
         return
      end if

      if (ra.lt.0.) ra = ra + 360.
      if (ra.gt.180.) then
         raminus = ra - 360.
      else
         raminus = ra
      end if

      yspace = 5.0
      iygrid = nint(dec/5.)+19
      ygrid = decbands(iygrid)

      xspace = ragaps(iygrid)
      if (abs(raminus).lt.xspace/4.0/2.0) then ! 4 for degs, 2 for half field
c     it's a 0hr field
         xgrid = 0.
         ixgrid = 0
      else
         xspace = ragaps(iygrid) / 4.0 ! in degs
         ixgrid = int(ra/xspace + 0.5)
         xgrid = ixgrid * xspace
      end if

      ifield = field0h(iygrid) + ixgrid
      xoff = ra - xgrid
      if (xoff.gt.180.) xoff = xoff - 360.
      yoff = dec - ygrid

      return
      end

c-----------------------------------------------------------------------------
c     Initialisation of the sgp_holes positions

      subroutine sgp_mask_init2
******************************************************************************
      implicit none

      integer i,j,jhole(2000),nshols,ilast,jfld
      real sxhole(3,2000),syhole(3,2000)

      real pival,degs_rads,ramins_rads,radfac,tpi
      parameter (pival=3.1415926535897)
      parameter (degs_rads=pival/180.,ramins_rads=pival/720.
     &,          radfac=pival/10800.,tpi=2.*pival)
      common/sgphole/sxhole,syhole,jhole,nshols
      common/last_field/ilast,jfld
******************************************************************************

      ilast = 0
      open (unit=30,file='sgpholes.lis',status='old',form='formatted')
      do i=1,2000
         read(30,*,end=999) (sxhole(j,i),syhole(j,i),j=1,3),jhole(i)
      enddo
 999  nshols = i-1
      close(30)

      return
      end

c-----------------------------------------------------------------------------
c     Initialisation of the ngp_holes positions

      subroutine ngp_mask_init2
******************************************************************************
      implicit none

      integer i,j,jhole(300),nnhols,ilast,jfld
      real xnhole(3,300),ynhole(3,300)

      real pival,degs_rads,ramins_rads,radfac,tpi
      parameter (pival=3.1415926535897)
      parameter (degs_rads=pival/180.,ramins_rads=pival/720.
     &,          radfac=pival/10800.,tpi=2.*pival)
      common/ngphole/xnhole,ynhole,jhole,nnhols
      common/last_field/ilast,jfld
******************************************************************************

      open (unit=30,file='ngpholes.lis',status='old',form='formatted')
      do i=1,300
         read(30,*,end=999) (xnhole(j,i),ynhole(j,i),j=1,3),jhole(i)
      enddo
 999  nnhols = i-1
      close(30)

      return
      end

c----------------------------------------------------------------------------
c     Returns inhol = 1 if point (xx,yy) lies in a drilled region

      subroutine test_sgp_holes2(xx,yy,jf,inhol)
******************************************************************************
      implicit none

      integer jf,inhol
      real xx,yy

      intrinsic abs,max,min

      integer j,jhole(2000),nholes
      real xhole(3,2000),yhole(3,2000),dx,dy,dx1,dy1,dx2,dy2,x0,y0,a1,a2

      real pival,degs_rads,ramins_rads,radfac,tpi
      parameter (pival=3.1415926535897)
      parameter (degs_rads=pival/180.,ramins_rads=pival/720.
     &,          radfac=pival/10800.,tpi=2.*pival)
      common/sgphole/xhole,yhole,jhole,nholes
******************************************************************************

      inhol = 0
      do j = 1,nholes
         if (jhole(j).ne.jf) goto 124
         x0 = xhole(1,j)
         y0 = yhole(1,j)

         if(y0.eq.yhole(3,j)) then ! It's a square hole
            if(xx.le.min(x0,xhole(3,j)).or.xx.ge.max(x0,xhole(3,j)))
     &         goto 124 ! speed-up
         end if
         dx = xx-x0
         dy = yy-y0
         dx1 = xhole(2,j)-x0
         dy1 = yhole(2,j)-y0
         dx2 = xhole(3,j)-x0
         dy2 = yhole(3,j)-y0

         if (abs(dx2).gt.0) then
            a1 = (dx*dy2 - dy*dx2)/(dx1*dy2 - dy1*dx2)
            a2 = (dx - a1*dx1)/dx2
            if (0.le.a1.and.a1.le.1.and.0.le.a2.and.a2.le.1) then
               inhol = 1
               return
            end if
         end if
 124     continue
      end do                    ! loop over holes

      return
      end

c-----------------------------------------------------------------------------
c     Returns inhol = 1 if point (xx,yy) lies in a drilled region

      subroutine test_ngp_holes2(xx,yy,jf,inhol)
******************************************************************************
      implicit none

      integer jf,inhol
      real xx,yy

      intrinsic abs,max,min

      integer j,jhole(300),nholes
      real xhole(3,300),yhole(3,300),dx,dy,dx1,dy1,dx2,dy2,x0,y0,a1,a2

      real pival,degs_rads,ramins_rads,radfac,tpi
      parameter (pival=3.1415926535897)
      parameter (degs_rads=pival/180.,ramins_rads=pival/720.
     &,          radfac=pival/10800.,tpi=2.*pival)
      common/ngphole/xhole,yhole,jhole,nholes
******************************************************************************

      inhol = 0
      do j = 1,nholes
         if (jhole(j).ne.jf) goto 123
         x0 = xhole(1,j)
         y0 = yhole(1,j)
         if(y0.eq.yhole(3,j)) then ! It's a square hole
            if(xx.le.min(x0,xhole(3,j)).or.xx.ge.max(x0,xhole(3,j)))
     &         goto 123 ! speed-up
         end if
         dx = xx-x0
         dy = yy-y0
         dx1 = xhole(2,j)-x0
         dy1 = yhole(2,j)-y0
         dx2 = xhole(3,j)-x0
         dy2 = yhole(3,j)-y0

         if (abs(dx2).gt.0) then
            a1 = (dx*dy2 - dy*dx2)/(dx1*dy2 - dy1*dx2)
            a2 = (dx - a1*dx1)/dx2
            if (0.le.a1.and.a1.le.1.and.0.le.a2.and.a2.le.1) then
               inhol = 1
               return
            end if
         end if
 123     continue
      end do  ! loop over holes

      return
      end

c-----------------------------------------------------------------------------
