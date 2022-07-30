! *****************************************************************
! Copyright (C) 2022 by Federico Marulli and Giorgio Lesci        *
! federico.marulli3@unibo.it                                      *
!                                                                 *
! This program is free software; you can redistribute it and/or   *
! modify it under the terms of the GNU General Public License as  *
! published by the Free Software Foundation; either version 2 of  *
! the License, or (at your option) any later version.             *
!                                                                 *
! This program is distributed in the hope that it will be useful, *
! but WITHOUT ANY WARRANTY; without even the implied warranty of  *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   *
! GNU General Public License for more details.                    *
!                                                                 *
! You should have received a copy of the GNU General Public       *
! License along with this program; if not, write to the Free      *
! Software Foundation, Inc.,                                      *
! 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.       *
! *****************************************************************

! CAMB interface for C++
! Author: Giorgio Lesci

module CAMBinterface
  use camb
  use DarkEnergyPPF
  use Quintessence
  use Recombination
  implicit none

contains

  ! »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  function GetCAMBparams () RESULT(params) BIND(C, NAME='GetCAMBparams')

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
    type(c_ptr) :: params
    type(CAMBparams), pointer :: p

    allocate(p)

    ! Use the C address as an opaque handle
    params = c_loc(p)

  end function GetCAMBparams

  ! »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  function GetCAMBdata () RESULT(data) BIND(C, NAME='GetCAMBdata')

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
    type(c_ptr) :: data
    type(CAMBdata), pointer :: pp

    allocate(pp)

    ! Use the C address as an opaque handle
    data = c_loc(pp)

  end function GetCAMBdata

  ! »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»  

  subroutine SetCAMBparams (params, ombh2, omch2, omnuh2, massless_nu, massive_nu, neutrino_hierarchy, omk, H0, dark_energy_model, w, wa, tau, cs2_lam, T_cmb, helium_fraction) BIND(C, NAME='SetCAMBparams')

    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_double, c_int
    type(c_ptr), intent(in), value :: params
    real(c_double), intent(in), value :: ombh2, omch2, omnuh2, massless_nu, omk, H0, w, wa, tau, T_cmb, helium_fraction, cs2_lam
    integer(c_int), intent(in), value :: massive_nu, neutrino_hierarchy, dark_energy_model
    real*8 :: omnuh2_sterile
    type(CAMBparams), pointer :: p

    call c_f_pointer(params, p)

    call CAMB_SetDefParams(p)

    p%ombh2 = ombh2
    p%omch2 = omch2
    p%omk = omk
    p%H0 = H0
    p%TCMB = T_cmb
    p%Yhe = helium_fraction

    select type(Reion=>p%Reion)
    class is (TTanhReionization)
       Reion%Reionization = .false.
       Reion%use_optical_depth = .true.
       Reion%optical_depth = tau
       Reion%delta_redshift = 1.5
       Reion%fraction = -1
       Reion%helium_redshift = 3.5
       Reion%helium_delta_redshift = 0.4
       Reion%helium_redshiftstart = Reion%helium_redshift + 5*Reion%helium_delta_redshift
    end select
    
    select type(Recomb=>p%Recomb)
    class is (TRecfast)
       Recomb%RECFAST_fudge_He = 0.86
       Recomb%RECFAST_Heswitch = 6
       Recomb%RECFAST_Hswitch = .true.
       Recomb%RECFAST_fudge = 1.14
       if (Recomb%RECFAST_Hswitch) then
          Recomb%RECFAST_fudge = Recomb%RECFAST_fudge - (1.14_dl - (1.105d0 + 0.02d0))
       end if
    end select
    
    do_bispectrum = .false.
    
    p%Scalar_initial_condition = 1 ! adiabatic initial perturbations
    
    ! Neutrinos settings
    omnuh2_sterile = 0.0
    p%share_delta_neff = .false.
    p%MassiveNuMethod = 1
    call CAMBparams_SetNeutrinoHierarchy(p, omnuh2, omnuh2_sterile, massless_nu+massive_nu, neutrino_hierarchy, massive_nu)

    ! Dark Energy settings
    if (allocated(p%DarkEnergy)) deallocate(p%DarkEnergy)
    if (dark_energy_model == 0) then
       allocate (TDarkEnergyFluid::p%DarkEnergy)
    else if (dark_energy_model == 1) then
       allocate (TDarkEnergyPPF::p%DarkEnergy)
    else if (dark_energy_model == 2) then
       allocate (TAxionEffectiveFluid::p%DarkEnergy)
    else if (dark_energy_model == 3) then
       allocate (TEarlyQuintessence::p%DarkEnergy)
    else
       error stop 'Unknown dark energy model!'
    end if

    select type(DarkEnergy=>p%DarkEnergy)
    class is (TDarkEnergyEqnOfState)
       DarkEnergy%no_perturbations = .false. ! don't change this: no perturbations is unphysical
       DarkEnergy%use_tabulated_w = .false. ! avoid tabulated w
       DarkEnergy%w_lam = w
       DarkEnergy%wa = wa
       DarkEnergy%cs2_lam = cs2_lam
    end select
    
    ! Additional settings to enhance performance  
    p%Accuracy%AccuratePolarization = .false.
    p%Accuracy%AccurateBB = .false.
    p%Accuracy%AccurateReionization = .false.
    p%Accuracy%AccuracyBoost = 1 ! Increase to decrease time steps, use more k values, etc.; decrease to speed up at cost of worse accuracy. Suggest 0.8 to 3.
    p%Accuracy%lAccuracyBoost = 1 ! Larger to keep more terms in the hierarchy evolution
    p%Accuracy%lSampleBoost = 1
      
    p%Max_l = 2200
    p%Min_l = 2
    p%Max_l_tensor = 1500
    p%Max_eta_k_tensor = 3000
    
    p%WantDerivedParameters = .false.

    p%WantScalars = .false.
    p%WantTensors = .false.
    p%WantVectors = .false.
    p%WantTransfer = .true.

    p%Reion%Reionization = .false.
    p%Want_CMB = .false.
    p%Want_CMB_lensing = .false.
    p%DoLensing = .false.
    p%WantCls = .false.

    p%DoLateRadTruncation = .true.
    p%Evolve_baryon_cs = .false.
    p%Evolve_delta_xe = .false.
    p%Evolve_delta_Ts =.false.
    p%Do21cm = .false.
    p%transfer_21cm_cl = .false.
    p%Log_lvalues  = .false.
    p%use_cl_spline_template = .true.

  end subroutine SetCAMBparams

  ! »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  subroutine SetCAMBPk (params, redshift, ns, As, pivot_scalar, accurate_massive_nu, kmax, nonlinear) BIND(C, NAME='SetCAMBPk')

    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_double, c_bool
    type(c_ptr), intent(in), value :: params
    real(c_double), intent(in), value :: redshift, ns, As, pivot_scalar, kmax
    logical(c_bool), intent(in), value :: accurate_massive_nu, nonlinear
    type(CAMBparams), pointer :: p

    call c_f_pointer(params, p)
    
    if (nonlinear) then
       p%NonLinear = 1
    end if
    
    ! Set InitialPower and Transfer parameters
    select type(InitPower=>p%InitPower)
    class is (TInitialPowerLaw)
       InitPower%As = As
       InitPower%ns = ns
       InitPower%pivot_scalar = pivot_scalar
       InitPower%r = 0
       InitPower%At = 0
    end select

    p%Transfer%PK_num_redshifts = 1
    p%Transfer%PK_redshifts(1) = redshift
    p%Transfer%kmax = kmax / (p%H0 / 100.)
    
    ! Additional settings for managing performances
    p%Transfer%high_precision = .false.
    p%Transfer%accurate_massive_neutrinos = accurate_massive_nu
    p%Transfer%k_per_logint = 0

  end subroutine SetCAMBPk

  ! »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  subroutine GetCAMBresults (params, data) BIND(C, NAME='GetCAMBresults')

    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_bool
    type(c_ptr), intent(in), value :: params, data
    type(CAMBparams), pointer :: p
    type(CAMBdata), pointer :: pp
    integer :: error

    call c_f_pointer(params, p)
    call c_f_pointer(data, pp)

    if (CAMBparams_Validate(p) .eqv. .false.) then
       error stop 'CAMB stopped'
    end if

    error = 0
    ThreadNum = 1
    transfer_power_var = 7
    pp%num_transfer_redshifts = 1
    
    !feedbacklevel = 1
    !debugmsgs = .true.
    call CAMB_GetResults(pp, p, error, .true., .true.)

  end subroutine GetCAMBresults

  ! »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  subroutine GetCAMBPk (data, Pk_out, minkh, dlnkh, npoints) BIND(C, NAME='GetCAMBPk')

    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_double, c_int
    type(c_ptr), intent(in), value ::  data
    real(c_double), intent(in), value :: minkh, dlnkh
    real(c_double) :: Pk_out(npoints)
    integer(c_int), intent(in), value :: npoints
    real(dl) :: Pk(npoints)
    integer :: var, i
    type(CAMBdata), pointer :: pp

    call c_f_pointer(data, pp)

    var = transfer_power_var
    call Transfer_GetMatterPowerD(pp, pp%MT, Pk, 1, minkh, dlnkh, npoints, var, var)

    do i=1, npoints+1, 1
       Pk_out(i) = Pk(i)
    end do
    return

  end subroutine GetCAMBPk

  ! »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  subroutine ReleaseCAMBparams (params) BIND(C, name='ReleaseCAMBparams')

    use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
    type(c_ptr), intent(in), value :: params
    type(CAMBparams), pointer :: p

    call c_f_pointer(params, p)
    deallocate(p)

  end subroutine ReleaseCAMBparams

  ! »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  subroutine ReleaseCAMBdata (data) BIND(C, name='ReleaseCAMBdata')

    use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
    type(c_ptr), intent(in), value :: data
    type(CAMBdata), pointer :: pp

    call c_f_pointer(data, pp)
    deallocate(pp)

  end subroutine ReleaseCAMBdata

end module CAMBinterface
