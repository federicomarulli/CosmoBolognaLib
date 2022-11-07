!     Standalone HEALPix subroutine converted to f77 by J C Hill 07/13/06
!     Get the original at http://healpix.jpl.nasa.gov
!     a f 'f88 pix2vecnest.f'

      !=========================================================================
      subroutine pix2vec_nest(nside, ipix, vector_x, vector_y, 
     & vector_z, vertex_n_x, vertex_n_y, vertex_n_z, vertex_e_x, 
     & vertex_e_y, vertex_e_z, vertex_s_x, vertex_s_y, vertex_s_z, 
     & vertex_w_x, vertex_w_y, vertex_w_z)
      !=========================================================================
      !     renders vector (x,y,z) coordinates of the nominal pixel center
      !     for the pixel number ipix (NESTED scheme)
      !     given the map resolution parameter nside
      !     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
      !     in the order N,E,S,W
      !=========================================================================
      IMPLICIT none
      INTEGER ipix, nside
      REAL*10 vector_x, vector_y, vector_z, vertex_n_x, 
     & vertex_n_y, vertex_n_z, vertex_e_x, vertex_e_y, 
     & vertex_e_z, vertex_s_x, vertex_s_y, vertex_s_z, 
     & vertex_w_x, vertex_w_y, vertex_w_z

      INTEGER npix, npface, ipf, ip_low, ip_trunc, ip_med, 
     & ip_hi, jrt, jr, nr, jpt, jp, kshift, nl4, ns_max
      parameter(ns_max=8192) ! 2^13 : largest nside allowed
      REAL*10 z, fn, fact1, fact2, sth, phi, PI
      parameter(PI=3.141592653589793238462643383279502884197_10)
      INTEGER pix2x(0:1023), pix2y(0:1023)
      INTEGER ix, iy, face_num
      ! common /xy_nest/ ix, iy, face_num ! can be useful to calling routine

      INTEGER jrll(1:12)
      INTEGER jpll(1:12)

      REAL*10 phi_nv, phi_wv, phi_sv, phi_ev, phi_up, phi_dn, 
     & z_nv, z_sv, sth_nv, sth_sv, hdelta_phi
      INTEGER iphi_mod, iphi_rat, kpix, jpix, ix_mk, 
     & iy_mk, ip_mk, id_mk, i
      ! LOGICAL do_vertex    ! this is unnecessary, since we will always want the vertices
      !------------------------------------------------------------------------
      if (nside.lt.1 .or. nside.gt.ns_max) stop 'nside out of range'
      npix = 12*nside**2     ! total number of pixels
      if (ipix.lt.0 .or. ipix.gt.npix-1) stop 'ipix out of range'

      do 10 i = 0, 1023
         pix2x(i) = 0
         pix2y(i) = 0
 10   continue

      ! coordinate of the lowest corner of each face, in unit of nside
      jrll(1) = 2
      jrll(2) = 2
      jrll(3) = 2
      jrll(4) = 2
      jrll(5) = 3
      jrll(6) = 3
      jrll(7) = 3
      jrll(8) = 3
      jrll(9) = 4
      jrll(10) = 4
      jrll(11) = 4
      jrll(12) = 4

      ! coordinate of the lowest corner of each face, in unit of nside/2
      jpll(1) = 1
      jpll(2) = 3
      jpll(3) = 5
      jpll(4) = 7
      jpll(5) = 0
      jpll(6) = 2
      jpll(7) = 4
      jpll(8) = 6
      jpll(9) = 1
      jpll(10) = 3
      jpll(11) = 5
      jpll(12) = 7

      !     initializes the array for the pixel number -> (x,y) mapping
      if (pix2x(1023).le.0) then
         ! constructs the array giving x and y in the face from pixel number
         ! for the nested (quad-cube like) ordering of pixels
         !
         ! the bits corresponding to x and y are interleaved in the pixel number
         ! one breaks up the pixel number by even and odd bits

         ! cc cf block data      data      pix2x(1023) /0/
         !
         !      print *, 'initiate pix2xy'
         do kpix = 0,1023               ! pixel number
            jpix = kpix
            ix_mk = 0
            iy_mk = 0
            ip_mk = 1                 ! bit position (in x and y)
         ! do while (jpix/=0)         ! go through all the bits
            do
               if (jpix.eq.0) exit    ! go through all the bits
!               id_mk = jpix           ! this and the next 4 lines are equivalent to id_mk=MODULO(jpix,2)
! 10            continue
!                  if (id_mk.ge.2) id_mk = id_mk - 2
!                  if (id_mk.lt.0) id_mk = id_mk + 2
!               if (id_mk.ge.2 .or. id_mk.lt.0) goto 10  ! bit value (in kpix), goes in ix_mk
               id_mk = MOD(jpix,2)
               jpix = jpix/2
               ix_mk = id_mk*ip_mk + ix_mk

!               id_mk = jpix           ! this and the next 4 lines are equivalent to id_mk=MODULO(jpix,2)
! 20            continue
!                  if (id_mk.ge.2) id_mk = id_mk - 2
!                  if (id_mk.lt.0) id_mk = id_mk + 2
!               if (id_mk.ge.2 .or. id_mk.lt.0) goto 20  ! bit value (in kpix), goes in iy_mk
               id_mk = MOD(jpix,2)
               jpix = jpix/2
               iy_mk = id_mk*ip_mk + iy_mk

               ip_mk = 2*ip_mk        ! next bit (in x and y)
            enddo
            pix2x(kpix) = ix_mk       ! in 0,31
            pix2y(kpix) = iy_mk       ! in 0,31
         enddo
      endif

      fn = DBLE(nside)
      fact1 = 1.0_10/(3.0_10*fn*fn)
      fact2 = 2.0_10/(3.0_10*fn)
      nl4   = 4*nside

      !    finds the face, and the number in the face
      npface = nside**2

      face_num = ipix/npface  ! face number in {0,11}
      ipf = MOD(ipix,npface)  ! pixel number in the face {0,npface-1}

      !     finds the x,y on the face (starting from the lowest corner)
      !     from the pixel number
      ip_low = MOD(ipf,1024)          ! content of the last 10 bits
      ip_trunc = ipf/1024             ! truncation of the last 10 bits
      ip_med = MOD(ip_trunc,1024)     ! content of the next 10 bits
      ip_hi = ip_trunc/1024           ! content of the high weight 10 bits
      
      ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
      iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)

      !     transforms this in (horizontal, vertical) coordinates
      jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
      jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}
      
      !     computes the z coordinate on the sphere
      jr = jrll(face_num+1)*nside - jrt - 1  ! ring number in {1,4*nside-1}

      nr = nside                 ! equatorial region (the most frequent)
      z = (2*nside-jr)*fact2
      kshift = MOD(jr - nside, 2)
      
      z_nv = (2*nside-jr+1)*fact2
      z_sv = (2*nside-jr-1)*fact2
      if (jr.eq.nside) then ! northern transition
         z_nv = 1.0_10 - (nside-1)**2 * fact1
      elseif (jr.eq.3*nside) then  ! southern transition
         z_sv = -1.0_10 + (nside-1)**2 * fact1
      endif

      if (jr.lt.nside) then     ! north pole region
         nr = jr
         z = 1.0_10 - nr*nr*fact1
         kshift = 0
         z_nv = 1.0_10 - (nr-1)**2 * fact1
         z_sv = 1.0_10 - (nr+1)**2 * fact1
      elseif (jr.gt.3*nside) then ! south pole region
         nr = nl4 - jr
         z = -1.0_10 + nr*nr*fact1
         kshift = 0
         z_nv = -1.0_10 + (nr+1)**2 * fact1
         z_sv = -1.0_10 + (nr-1)**2 * fact1
      endif

      !     computes the phi coordinate on the sphere, in [0,2Pi]
      jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
      if (jp.gt.nl4) jp = jp - nl4
      if (jp.lt.1)   jp = jp + nl4

      phi = (jp - (kshift+1)*0.5_10) * (PI/(2.0_10*nr))

      sth = SQRT((1.0_10-z)*(1.0_10+z))
      vector_x = sth * COS(phi)
      vector_y = sth * SIN(phi)
      vector_z = z

      phi_nv = phi
      phi_sv = phi

      phi_up = 0.0_10
      iphi_mod = MOD(jp-1,nr) ! in {0,1,... nr-1}
      iphi_rat = (jp-1) / nr  ! in {0,1,2,3}
      if (nr.gt.1) phi_up=(PI/2.0_10)*(iphi_rat+iphi_mod/(DBLE(nr-1)))
      phi_dn         =(PI/2.0_10)*(iphi_rat+(iphi_mod+1)/(DBLE(nr+1)))
      if (jr.lt.nside) then          ! North polar cap
         phi_nv = phi_up
         phi_sv = phi_dn
      elseif (jr.gt.3*nside) then    ! South polar cap
         phi_nv = phi_dn
         phi_sv = phi_up
      elseif (jr.eq.nside) then      ! North transition
         phi_nv = phi_up
      elseif (jr.eq.3*nside) then    ! South transition
         phi_sv = phi_up
      endif

      hdelta_phi = PI / (4.0_10*nr)

      ! west vertex
      phi_wv     = phi - hdelta_phi
      vertex_w_x = sth * COS(phi_wv)
      vertex_w_y = sth * SIN(phi_wv)
      vertex_w_z = z

      ! east vertex
      phi_ev     = phi + hdelta_phi
      vertex_e_x = sth * COS(phi_ev)
      vertex_e_y = sth * SIN(phi_ev)
      vertex_e_z = z

      ! north vertex
      sth_nv = SQRT((1.0_10-z_nv)*(1.0_10+z_nv))
      vertex_n_x = sth_nv * COS(phi_nv)
      vertex_n_y = sth_nv * SIN(phi_nv)
      vertex_n_z = z_nv

      ! south vertex
      sth_sv = SQRT((1.0_10-z_sv)*(1.0_10+z_sv))
      vertex_s_x = sth_sv * COS(phi_sv)
      vertex_s_y = sth_sv * SIN(phi_sv)
      vertex_s_z = z_sv

      return
      end subroutine pix2vec_nest ! pix2vec_nest
