c ****************************************************************
c
      program read_pk_sum
c
c ****************************************************************
c
      call load_pkrenorm
c
      call load_pkcorr
c
      call save_pkred
c
      end
c
c ****************************************************************
c
      subroutine load_pkrenorm
c
c ****************************************************************
c
      implicit none
      integer ik, ikmax, ik_max, izred
      parameter(ikmax=5000)
      real*8 ak(ikmax), pk0(ikmax), pk2(ikmax), pk4(ikmax)
      real*8 dummy, pk0_EH(ikmax), pk2_EH(ikmax), pk4_EH(ikmax)
      common /pkrenorm/ ak,pk0_EH,pk0,pk2_EH,pk2,pk4_EH,pk4,ik_max
c
      write(6,*) 'redshift: z=0.5[0], z=1[1], z=2[2], z=3[3]'
      read(5,*) izred
      if(izred.eq.0) then 
         open(10,file='pkrenorm_Closure_red_2loop_zred05.dat',
     &     status='unknown')
      elseif(izred.eq.1) then 
         open(10,file='pkrenorm_Closure_red_2loop_zred1.dat',
     &     status='unknown')
      elseif(izred.eq.2) then 
         open(10,file='pkrenorm_Closure_red_2loop_zred2.dat',
     &     status='unknown')
      elseif(izred.eq.3) then 
         open(10,file='pkrenorm_Closure_red_2loop_zred3.dat',
     &     status='unknown')
      else
         stop
      endif
c
      do ik=1, ikmax
         read(10,*,END=10) dummy, pk0_EH(ik), dummy, pk0(ik), 
     &        pk2_EH(ik), dummy, pk2(ik), pk4_EH(ik), dummy, pk4(ik)
      enddo
c
 10   ik_max = ik
c
      end
c
c ****************************************************************
c
      subroutine load_pkcorr
c
c ****************************************************************
c
      implicit none
      integer ik, ikmax, ik_max, ifile
      parameter(ikmax=5000)
      real*8 ak(ikmax), pk0(ikmax), pk2(ikmax), pk4(ikmax)
      real*8 pk0_EH(ikmax), pk2_EH(ikmax), pk4_EH(ikmax)
      real*8  pk0_corr(ikmax), pk2_corr(ikmax), pk4_corr(ikmax)
      real*8 boost, dummy
      common /pkrenorm/ ak,pk0_EH,pk0,pk2_EH,pk2,pk4_EH,pk4,ik_max
c
 5    write(6,*) 'load file [1--3] or exit [0]:  '
      write(6,*) '   corr_pkred.dat  [1]'
      write(6,*) '   corr_pkred2.dat [2]'
      write(6,*) '   corr_pkred3.dat [3]'
      write(6,*) '   '
      read(5,*) ifile
      if(ifile.eq.0) then 
         goto 100
      elseif(ifile.eq.1) then
         open(9,file='corr_pkred.dat',status='unknown')
      elseif(ifile.eq.2) then
         open(9,file='corr_pkred2.dat',status='unknown')
      elseif(ifile.eq.3) then
         open(9,file='corr_pkred3.dat',status='unknown')
      else
         stop
      endif
c
      do ik=1, ikmax
         read(9,*,END=20) ak(ik), dummy, dummy, pk0_corr(ik), dummy, 
     &        pk2_corr(ik), dummy, pk4_corr(ik)
      enddo
c
 20   ik_max = min(ik, ik_max)
      close(9)
c
      write(6,*) 'input boost factor'
      read(5,*) boost
c
      do ik=1, ik_max
         pk0(ik) = pk0(ik) + boost*pk0_corr(ik)
         pk2(ik) = pk2(ik) + boost*pk2_corr(ik)
         pk4(ik) = pk4(ik) + boost*pk4_corr(ik)
      enddo
c
      goto 5
c
 100  continue
c
      end
c
c ****************************************************************
c
      subroutine save_pkred
c
c ****************************************************************
c     
      implicit none
      integer ik, ikmax, ik_max
      parameter(ikmax=5000)
      real*8  ak(ikmax), pk0_EH(ikmax), pk2_EH(ikmax), pk4_EH(ikmax)
      real*8  pk0(ikmax), pk2(ikmax), pk4(ikmax)
      common /pkrenorm/ ak,pk0_EH,pk0,pk2_EH,pk2,pk4_EH,pk4,ik_max
c
c
      open(11,file='pkrenorm_Closure_red_2loop_corr.dat',
     &     status='unknown')
c
      do ik=1, ik_max
         write(11,'(1p7e16.8)') ak(ik), pk0_EH(ik), pk0(ik), 
     &        pk2_EH(ik), pk2(ik), pk4_EH(ik), pk4(ik)
      enddo
      close(11)
c
      end
