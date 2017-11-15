c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      real*10 function twoqz(ra,dec,verbose)
      real*10 ra,dec
      integer verbose
c
c        parameters
      integer NHEM
      parameter (NHEM=2)
c        externals
      integer access,lnblnk
c        data variables
      include 'mangdir.inc'
      character*24 subdir
      integer unit
      logical init
c        local (automatic) variables
      character*1 go
      character*64 dat
      character*128 tmpnam
      integer idec,idir,ihem,ifmt,ira,found
      real*4 decmn,decmx,r,ramn,ramx
c        local variables to be saved
      integer ndec(NHEM),nra(NHEM)
      real*4 decmin(NHEM),decstep(NHEM),ramin(NHEM),rastep(NHEM)
      real*4 array(4510,330,2)
      save ndec,nra
      save decmin,decstep,ramin,rastep
      save array
      logical ex
      character*3 re
c *
c * Completeness of 2QZ, the April 2001 2dF 10k quasar redshift survey.
c *
c  Input: ra, dec = RA & Dec in degrees.
c         verbose = 0 quiet
c                 = 1 normal verbosity
c                 = 3 details (useful if twoqz is run interactively)
c Output: twoqz: 0 = empty to 1 = complete.
c
c        mangle directory
      include 'mangdir.data'
c        subdirectory of mangle directory containing 2QZ mask data
      data subdir /'masks/2qz10k'/
c        initialize on first call
      data init /.true./
c        fortran unit
      data unit /7/
c
c--------read arrays
      if (init) then
c        each hemisphere
        do 170 ihem=1,NHEM
          found=0
c        try current directory, then subdirectory, then mangle directory
          do 160 idir=1,3
c        try unformatted then formatted files
            do 150 ifmt=1,2
c        unformatted
              if (ifmt.eq.1) then
                if (ihem.eq.1) then
                  dat='ngc_obscomp.dat'
                elseif (ihem.eq.2) then
                  dat='sgc_obscomp.dat'
                endif
c        formatted
              elseif (ifmt.eq.2) then
                if (ihem.eq.1) then
                  dat='ngc_obscomp.txt'
                elseif (ihem.eq.2) then
                  dat='sgc_obscomp.txt'
                endif
              endif
c        current directory
              if (idir.eq.1) then
                tmpnam=dat
c        subdirectory
              elseif (idir.eq.2) then
                tmpnam=subdir(1:lnblnk(subdir))
     *            //'/'//dat(1:lnblnk(dat))
c        mangle directory
              elseif (idir.eq.3) then
                call getenv(mangenv,tmpnam)
                if (found.eq.1) goto 310
                if (tmpnam.eq.' ') goto 320
                mangdir=tmpnam
                tmpnam=mangdir(1:lnblnk(mangdir))
     *            //'/'//subdir(1:lnblnk(subdir))
     *            //'/'//dat(1:lnblnk(dat))
              endif
c        found the file
              inquire(FILE=tmpnam(1:lnblnk(tmpnam)), EXIST=ex, READ=re)
              if (ex.and.(re.eq.'YES')) then
c              if (access(tmpnam(1:lnblnk(tmpnam)),' r').eq.0) then
                if (ifmt.eq.1) then
                  open (unit,file=tmpnam,form='unformatted',
     *              status='old',err=340)
                else
                  open (unit,file=tmpnam,status='old',err=340)
                endif
                found=1
c        did not find the file
              else
c        try formatted file
                if (ifmt.eq.1) goto 150
c        try another directory
                if (idir.le.2) goto 160
                if (found.eq.1) goto 310
c        ran out of options
                goto 330
              endif
c        found the file: now read it
              rewind (unit)
              if (verbose.ge.1)
     *          print *,'reading ',tmpnam(1:lnblnk(tmpnam)),' ...'
c        unformatted
              if (ifmt.eq.1) then
                read (unit,end=120,err=120) ramin(ihem)
                read (unit,end=120,err=120) decmin(ihem)
                read (unit,end=120,err=120) nra(ihem)
                read (unit,end=120,err=120) ndec(ihem)
                read (unit,end=120,err=120) rastep(ihem)
                read (unit,end=120,err=120) decstep(ihem)
                read (unit,end=120,err=120) ((array(ira,idec,ihem),
     *            ira=1,nra(ihem)),idec=1,ndec(ihem))
                goto 130
c        error reading unformatted file
  120           if (verbose.ge.1)
     *            print *,'error reading ',tmpnam(1:lnblnk(tmpnam)),
     *            '; retry:'
c        try formatted file
                goto 150
c        successful read
  130           continue
c        formatted
              elseif (ifmt.eq.2) then
                read (unit,*,end=350,err=360) ramin(ihem)
                read (unit,*,end=350,err=360) decmin(ihem)
                read (unit,*,end=350,err=360) nra(ihem)
                read (unit,*,end=350,err=360) ndec(ihem)
                read (unit,*,end=350,err=360) rastep(ihem)
                read (unit,*,end=350,err=360) decstep(ihem)
                read (unit,*,end=350,err=360) ((array(ira,idec,ihem),
     *            ira=1,nra(ihem)),idec=1,ndec(ihem))
              endif
c        close file
              close (unit)
c             write (*,'(2x,i8,9i12)') (ira,ira=1,10)
c             do idec=1,10
c               write (*,'(i2,10g12.4)')
c    *            idec,(array(ira,idec,ihem),ira=1,10)
c             enddo

c        option to write formatted copy
              if (verbose.ge.2.and.ifmt.eq.1) then
                write (*,'(" write formatted copy of ",a,
     *            "? [CR,n=no, y=yes]: ",$)') tmpnam(1:lnblnk(tmpnam))
                read (*,'(a1)',end=135,err=135) go
                if (go.eq.' '.or.go.eq.'n'.or.go.eq.'N'
     *            .or.go.eq.'x'.or.go.eq.'X'
     *            .or.go.eq.'q'.or.go.eq.'Q') then
                  continue
                else
                  if (ihem.eq.1) then
                    dat='ngc_obscomp.txt'
                  elseif (ihem.eq.2) then
                    dat='sgc_obscomp.txt'
                  endif
c        current directory
                  if (idir.eq.1) then
                    tmpnam=dat
c        subdirectory
                  elseif (idir.eq.2) then
                    tmpnam=subdir(1:lnblnk(subdir))
     *                //'/'//dat(1:lnblnk(dat))
c        mangle directory
                  elseif (idir.eq.3) then
                    tmpnam=mangdir(1:lnblnk(mangdir))
     *                //'/'//subdir(1:lnblnk(subdir))
     *                //'/'//dat(1:lnblnk(dat))
                  endif
                  open (unit,file=tmpnam,err=370)
                  rewind (unit)
                  print *,'writing ',tmpnam(1:lnblnk(tmpnam)),' ...'
                  write (unit,'(g16.8)',err=380) ramin(ihem)
                  write (unit,'(g16.8)',err=380) decmin(ihem)
                  write (unit,'(i8)',err=380) nra(ihem)
                  write (unit,'(i8)',err=380) ndec(ihem)
                  write (unit,'(g16.8)',err=380) rastep(ihem)
                  write (unit,'(g16.8)',err=380) decstep(ihem)
                  write (unit,'(6g16.8)',err=380)((array(ira,idec,ihem),
     *              ira=1,nra(ihem)),idec=1,ndec(ihem))
                  close (unit)
                endif
 135            continue
              endif

c        option to write unformatted copy
c              if (verbose.ge.2.and.ifmt.eq.2) then
              if (verbose.ge.1.and.ifmt.eq.2) then
                write (*,'(" write unformatted copy of ",a,
     *            "? [CR,n=no, y=yes]: ",$)') tmpnam(1:lnblnk(tmpnam))
                read (*,'(a1)',end=140,err=140) go
                if (go.eq.' '.or.go.eq.'n'.or.go.eq.'N'
     *            .or.go.eq.'x'.or.go.eq.'X'
     *            .or.go.eq.'q'.or.go.eq.'Q') then
                  continue
                else
                  if (ihem.eq.1) then
                    dat='ngc_obscomp.dat'
                  elseif (ihem.eq.2) then
                    dat='sgc_obscomp.dat'
                  endif
c        current directory
                  if (idir.eq.1) then
                    tmpnam=dat
c        subdirectory
                  elseif (idir.eq.2) then
                    tmpnam=subdir(1:lnblnk(subdir))
     *                //'/'//dat(1:lnblnk(dat))
c        mangle directory
                  elseif (idir.eq.3) then
                    tmpnam=mangdir(1:lnblnk(mangdir))
     *                //'/'//subdir(1:lnblnk(subdir))
     *                //'/'//dat(1:lnblnk(dat))
                  endif
                  open (unit,file=tmpnam,form='unformatted',err=370)
                  rewind (unit)
                  print *,'writing ',tmpnam(1:lnblnk(tmpnam)),' ...'
                  write (unit,err=380) ramin(ihem)
                  write (unit,err=380) decmin(ihem)
                  write (unit,err=380) nra(ihem)
                  write (unit,err=380) ndec(ihem)
                  write (unit,err=380) rastep(ihem)
                  write (unit,err=380) decstep(ihem)
                  write (unit,err=380) ((array(ira,idec,ihem),
     *              ira=1,nra(ihem)),idec=1,ndec(ihem))
                  close (unit)
                endif
  140           continue
              endif
c        successful read
              goto 170
  150       continue
  160     continue
  170   continue
        init=.false.
      endif
c--------the routine
      twoqz=0._10
      do ihem=1,NHEM
        decmn=decmin(ihem)-decstep(ihem)*.5_10
        decmx=decmn+decstep(ihem)*ndec(ihem)
        if (dec.ge.decmn.and.dec.le.decmx) then
          ramn=ramin(ihem)-rastep(ihem)*.5_10
          ramx=ramn+rastep(ihem)*nra(ihem)
          ira=(ra-ramn)/360._10
          r=ra-ira*360._10
          if (r.lt.ramn) r=r+360._10
          if (verbose.ge.3) then
            print *,'ra, ramin, ramax =',r,ramn,ramx
            print *,'dec, decmin, decmax =',dec,decmn,decmx
          endif
          if (r.le.ramx) then
            ira=1+int((r-ramn)/rastep(ihem))
            idec=1+int((dec-decmn)/decstep(ihem))
            if (verbose.ge.3) then
              print *,'ira, idec, ihem, ra, dec =',ira,idec,ihem,
     *          ramin(ihem)+(ira-1)*rastep(ihem),
     *          decmin(ihem)+(idec-1)*decstep(ihem)
            endif
            twoqz=array(ira,idec,ihem)
c        round to 6 sig fig
            twoqz=dble(nint(twoqz*1.e6_10))/1.e6_10
            goto 200
          endif
        endif
      enddo
  200 continue
      return
c
c--------errors
  310 print *,'twoqz: you may have a problem with endianness;'
      print *,'please read HELP.unformatted in the mangle directory'
      stop
  320 print *,'twoqz: can''t find ',dat(1:lnblnk(dat))
      print *,'maybe I''d find it if environment variable ',
     *  mangenv(1:lnblnk(mangenv)),' were set to'
      print *,mangdir(1:lnblnk(mangdir))
      stop
  330 print *,'twoqz: can''t find ',dat(1:lnblnk(dat)),
     *  ' in directories . or '
      print *,mangdir(1:lnblnk(mangdir))//'/'//subdir(1:lnblnk(subdir))
      stop
  340 print *,'twoqz: error opening ',tmpnam(1:lnblnk(tmpnam))
      stop
  350 print *,'twoqz: premature EOF on ',tmpnam(1:lnblnk(tmpnam))
      stop
  360 print *,'twoqz: error reading ',tmpnam(1:lnblnk(tmpnam))
      stop
  370 print *,'twoqz: error opening ',tmpnam(1:lnblnk(tmpnam))
      stop
  380 print *,'twoqz: error writing to ',tmpnam(1:lnblnk(tmpnam))
      stop
c
      end
c
