!     Modules used by cmbmain and other routines.
!
!     Code for Anisotropies in the Microwave Background
!     by Antony Lewis (http://cosmologist.info) and Anthony Challinor
!     See readme.html for documentation. 
!
!     Based on CMBFAST  by  Uros Seljak and Matias Zaldarriaga, itself based
!     on Boltzmann code written by Edmund Bertschinger, Chung-Pei Ma and Paul Bode.
!     Original CMBFAST copyright and disclaimer:
!
!     Copyright 1996 by Harvard-Smithsonian Center for Astrophysics and
!     the Massachusetts Institute of Technology.  All rights reserved.
!
!     THIS SOFTWARE IS PROVIDED "AS IS", AND M.I.T. OR C.f.A. MAKE NO
!     REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
!     By way of example, but not limitation,
!     M.I.T. AND C.f.A MAKE NO REPRESENTATIONS OR WARRANTIES OF
!     MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
!     THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
!     ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
!
!     portions of this software are based on the COSMICS package of
!     E. Bertschinger.  See the LICENSE file of the COSMICS distribution
!     for restrictions on the modification and distribution of this software.

!     Modified by Shun Saito
!     Every modification is noted by !===== MODIFIED BY S.S. =====


!==================================================================================
module ModelParams
  use precision
  use Ranges
  use InitialPower
  use Reionization
  use Recombination
  
  implicit none    
  public
  
  character(LEN=*), parameter :: version = 'July_11'
  
  integer :: FeedbackLevel = 0 !if >0 print out useful information about the model
  
  logical, parameter :: DebugMsgs=.false. !Set to true to view progress and timing
  
  logical, parameter :: DebugEvolution = .false. 
  !Set to true to do all the evolution for all k
  
  logical ::  do_bispectrum  = .false. 
  logical, parameter :: hard_bispectrum = .false. 
  ! e.g. warm inflation where delicate cancellations
  
  logical, parameter :: full_bessel_integration = .false. 
  !(go into the tails when calculating the sources)
  
  integer, parameter :: Nu_int = 0, Nu_trunc=1, Nu_approx = 2, Nu_best = 3
  !For CAMBparams%MassiveNuMethod
  !Nu_int: always integrate distribution function
  !Nu_trunc: switch to expansion in velocity once non-relativistic
  !Nu_approx: approximate scheme - good for CMB, but not formally 
  !correct and no good for matter power
  !Nu_best: automatically use mixture which is fastest and most accurate
  
  integer, parameter :: max_Nu = 5 !Maximum number of neutrino species    
  integer, parameter :: max_transfer_redshifts = 50    
  integer, parameter :: fileio_unit = 13 !Any number not used elsewhere will do
  integer, parameter :: outCOBE=0, outNone=1
  
  integer :: max_bessels_l_index  = 1000000
  real(dl) :: max_bessels_etak = 1000000*2
  
  
  real(dl), parameter ::  OutputDenominator =twopi
  !When using outNone the output is l(l+1)Cl/OutputDenominator
  
  Type(Regions) :: TimeSteps
  
  type TransferParams
     logical     ::  high_precision
     integer     ::  num_redshifts
     real(dl)    ::  kmax         
     !these are acutally q values, but same as k for CP%flat
     integer     ::  k_per_logint ! ..
     real(dl)    ::  redshifts(max_transfer_redshifts)         
  end type TransferParams
  
  !other variables, options, derived variables, etc.
  
  !===== MODIFIED BY S.S. =====
  integer, parameter :: NonLinear_none=0, NonLinear_PkHF = 1, NonLinear_Lens=2
  integer, parameter :: NonLinear_PkSPT1loop = 3

  integer, parameter :: LinKaiser=1, NLKaiser=2
  integer, parameter :: NLKaiser_FoG=3, NLKaiser_AB_FoG=4
  !NonLinear RSD model:
  ! 1 :: Linear Kaiser
  ! 2 :: Nonlinear Kaiser
  ! 3 :: Nonlinear Kaiser * FoG
  ! 4 :: (Nonlinear Kaiser + A&B corrections) * FoG
  
  
  ! Main parameters type
  type CAMBparams
     
     logical   :: WantCls, WantTransfer
     logical   :: WantScalars, WantTensors, WantVectors
     logical   :: DoLensing
     logical   :: want_zstar, want_zdrag     !!JH for updated BAO likelihood.
     integer   :: NonLinear
     !added by S.S.
     integer   :: NonLinear_RSDmodel
     logical   :: DoNonLocalBias, DoNonLocalBias_A01
     
     integer   :: Max_l, Max_l_tensor
     real(dl)  :: Max_eta_k, Max_eta_k_tensor
     ! _tensor settings only used in initialization, 
     !Max_l and Max_eta_k are set to the tensor variables if only tensors requested
     
     real(dl)  :: omegab, omegac, omegav, omegan
     !Omega baryon, CDM, Lambda and massive neutrino
     real(dl)  :: H0,TCMB,yhe,Num_Nu_massless
     integer   :: Num_Nu_massive
     
     logical :: Nu_mass_splittings
     integer   :: Nu_mass_eigenstates  !1 for degenerate masses
     real(dl)  :: Nu_mass_degeneracies(max_nu)
     real(dl)  :: Nu_mass_fractions(max_nu)
     !The ratios of the masses
     
     integer   :: Scalar_initial_condition 
     !must be one of the initial_xxx values defined in GaugeInterface
     
     integer   :: OutputNormalization  
     !outNone, outCOBE, or C_OutputNormalization=1 if > 1
     
     logical   :: AccuratePolarization
     !Do you care about the accuracy of the polarization Cls?
     
     logical   :: AccurateBB
     !Do you care about BB accuracy (e.g. in lensing)
     
     !Reionization settings - used if Reion%Reionization=.true.
     logical   :: AccurateReionization
     !Do you care about pecent level accuracy on EE signal from reionization?
     
     integer   :: MassiveNuMethod
     
     type(InitialPowerParams) :: InitPower  
     !see power_tilt.f90 - you can change this
     type(ReionizationParams) :: Reion
     type(RecombinationParams):: Recomb
     type(TransferParams)     :: Transfer 
     
     real(dl) ::  InitialConditionVector(1:10) 
     !Allow up to 10 for future extensions
     !ignored unless Scalar_initial_condition == initial_vector
     
     logical OnlyTransfers !Don't use initial power spectrum data, 
     !instead get Delta_q_l array
     !If trye, sigma_8 is not calculated either
     
     !Derived parameters, not set initially
     type(ReionizationHistory) :: ReionHist
     
     logical flat,closed,open
     real(dl) omegak
     real(dl) curv,r, Ksign !CP%r = 1/sqrt(|CP%curv|), CP%Ksign = 1,0 or -1
     real(dl) tau0,chi0 !time today and rofChi(CP%tau0/CP%r) 
     
     
  end type CAMBparams
  
  type(CAMBparams) CP  !Global collection of parameters
  
  
  real(dl) scale !relative to CP%flat. e.g. for scaling lSamp%l sampling.
  
  logical ::call_again = .false.
  !if being called again with same parameters to get different thing
  
  
  !     grhom =kappa*a^2*rho_m0
  !     grhornomass=grhor*number of massless neutrino species
  !     taurst,taurend - time at start/end of recombination
  !     dtaurec - dtau during recombination
  !     adotrad - a(tau) in radiation era
  
  real(dl) grhom,grhog,grhor,grhob,grhoc,grhov,grhornomass,grhok
  real(dl) taurst,dtaurec,taurend, tau_maxvis,adotrad
  
  !Neutrinos
  real(dl) grhormass(max_nu)
  
  !     nu_masses=m_nu*c**2/(k_B*T_nu0)      
  real(dl) :: nu_masses(max_nu) 
  
  real(dl) akthom !sigma_T * (number density of protons now)
  real(dl) fHe !n_He_tot / n_H_tot
  real(dl) Nnow
  
  
  integer :: ThreadNum = 0 
  !If zero assigned automatically, obviously only used if parallelised
  
  !Parameters for checking/changing overall accuracy
  !If HighAccuracyDefault=.false., the other parameters equal to 1 
  !corresponds to ~0.3% scalar C_l accuracy
  !If HighAccuracyDefault=.true., the other parameters equal to 1 
  !corresponds to ~0.1% scalar C_l accuracy (at L>600)
  logical :: HighAccuracyDefault = .false.
  
  real(dl) :: lSampleBoost=1._dl
  !Increase lSampleBoost to increase sampling in lSamp%l for Cl interpolation
  
  real(dl) :: AccuracyBoost =1._dl  
  
  !Decrease step sizes, etc. by this parameter. Useful for checking accuracy.
  !Can also be used to improve speed significantly if less accuracy is required.
  !or improving accuracy for extreme models. 
  !Note this does not increase lSamp%l sampling or massive neutrino q-sampling
  
  real(sp) :: lAccuracyBoost=1. 
  !Boost number of multipoles integrated in Boltzman heirarchy
  
  integer, parameter :: lmin = 2  
  !must be either 1 or 2       
  
  real(dl), parameter :: OmegaKFlat = 5e-7_dl !Value at which to use flat code
  
  real(dl),parameter :: tol=1.0d-4 !Base tolerance for integrations
  
  !     used as parameter for spline - tells it to use 'natural' end values
  real(dl), parameter :: spl_large=1.e40_dl
  
  integer, parameter:: l0max=4000
  
  !     lmax is max possible number of l's evaluated
  integer, parameter :: lmax_arr = l0max 
  
contains
  !////////////////////////////////////////////////////////////////////////
  subroutine CAMBParams_Set(P, error, DoReion)
    use constants
    type(CAMBparams), intent(in) :: P
    real(dl) GetOmegak, fractional_number
    integer, optional :: error !Zero if OK
    logical, optional :: DoReion
    logical WantReion
    integer nu_i,actual_massless
    external GetOmegak
    real(dl), save :: last_tau0
    !Constants in SI units

    if ((P%WantTensors .or. P%WantVectors).and. & 
         P%WantTransfer .and. .not. P%WantScalars) then
       write (*,*) 'Cannot generate tensor C_l and transfer without scalar C_l'
       if (present(error)) then
          error = 1
          return
       else
          stop
       end if
    end if
    
    if (present(DoReion)) then
       WantReion = DoReion
    else
       WantReion = .true.
    end if
    
    CP=P
    
    CP%Max_eta_k = max(CP%Max_eta_k,CP%Max_eta_k_tensor)
    
    if (CP%WantTransfer) then
       CP%WantScalars=.true.
       if (.not. CP%WantCls) then
          CP%AccuratePolarization = .false.
          CP%Reion%Reionization = .false.
       end if
    else
       CP%transfer%num_redshifts=0
    end if
    
    if (CP%Omegan == 0 .and. CP%Num_Nu_Massive /=0) then
       CP%Num_Nu_Massless = CP%Num_Nu_Massless + CP%Num_Nu_Massive
       CP%Num_Nu_Massive  = 0
    end if
    
    if (CP%Num_nu_massive > 0) then
       if (.not. CP%Nu_mass_splittings) then
          !Default totally degenerate masses
          CP%Nu_mass_eigenstates = 1
          CP%Nu_mass_degeneracies(1) = CP%Num_Nu_Massive 
          CP%Nu_mass_fractions(1) = 1
       else
          if (CP%Nu_mass_degeneracies(1)==0) & 
               CP%Nu_mass_degeneracies(1) = CP%Num_Nu_Massive 
          if (abs(sum(CP%Nu_mass_fractions(1:CP%Nu_mass_eigenstates))-1) > 1e-4) &
               stop 'Nu_mass_fractions do not add up to 1'
          
          if (abs(sum(CP%Nu_mass_degeneracies(1:CP%Nu_mass_eigenstates)) & 
               -CP%Num_nu_massive) >1e-4 ) &
               stop 'nu_mass_eigenstates do not add up to num_nu_massive'
          if (CP%Nu_mass_eigenstates==0) stop  & 
               'Have Num_nu_massive>0 but no nu_mass_eigenstates'
          
       end if
    else
       CP%Nu_mass_eigenstates = 0
    end if
    
    if ((CP%WantTransfer).and. CP%MassiveNuMethod==Nu_approx) then
       CP%MassiveNuMethod = Nu_trunc
    end if
    
    CP%omegak = GetOmegak()
    
    CP%flat = (abs(CP%omegak) <= OmegaKFlat)
    CP%closed = CP%omegak < -OmegaKFlat
    
    CP%open = .not.CP%flat.and..not.CP%closed
    if (CP%flat) then
       CP%curv=0
       CP%Ksign=0
       CP%r=1._dl !so we can use tau/CP%r, etc, where CP%r's cancel
    else   
       CP%curv=-CP%omegak/((c/1000)/CP%h0)**2
       CP%Ksign =sign(1._dl,CP%curv)
       CP%r=1._dl/sqrt(abs(CP%curv))
    end if
    !  grho gives the contribution to the expansion rate from: (g) photons,
    !  (r) one flavor of relativistic neutrino (2 degrees of freedom),
    !  (m) nonrelativistic matter (for Omega=1).  grho is actually
    !  8*pi*G*rho/c^2 at a=1, with units of Mpc**(-2).
    !  a=tau(Mpc)*adotrad, with a=1 today, assuming 3 neutrinos.
    !  (Used only to set the initial conformal time.)
    
    !H0 is in km/s/Mpc
    
    grhom = 3*CP%h0**2/c**2*1000**2 !3*h0^2/c^2 (=8*pi*G*rho_crit/c^2)
    
    !grhom=3.3379d-11*h0*h0 
    grhog = kappa/c**2*4*sigma_boltz/c**3*CP%tcmb**4*Mpc**2 
    !8*pi*G/c^2*4*sigma_B/c^3 T^4
    ! grhog=1.4952d-13*tcmb**4
    grhor = 7._dl/8*(4._dl/11)**(4._dl/3)*grhog 
    !7/8*(4/11)^(4/3)*grhog (per neutrino species)
    !grhor=3.3957d-14*tcmb**4
    !correction for fractional number of neutrinos, e.g. 3.04 to give 
    !slightly higher T_nu hence rhor
    !Num_Nu_massive is already integer, Num_Nu_massless can contain fraction
    !We assume all eigenstates affected the same way
    fractional_number  = CP%Num_Nu_massless + CP%Num_Nu_massive
    actual_massless = int(CP%Num_Nu_massless + 1e-6)
    grhor = grhor * fractional_number/(actual_massless + CP%Num_Nu_massive)
    
    grhornomass=grhor*actual_massless
    grhormass=0
    do nu_i = 1, CP%Nu_mass_eigenstates
       grhormass(nu_i)=grhor*CP%Nu_mass_degeneracies(nu_i)
    end do
    grhoc=grhom*CP%omegac
    grhob=grhom*CP%omegab
    grhov=grhom*CP%omegav
    grhok=grhom*CP%omegak
    !  adotrad gives the relation a(tau) in the radiation era:
    adotrad = sqrt((grhog+grhornomass+sum(grhormass(1:CP%Nu_mass_eigenstates)))/3)
    
    
    Nnow = CP%omegab*(1-CP%yhe)*grhom*c**2/kappa/m_H/Mpc**2
    
    akthom = sigma_thomson*Nnow*Mpc
    !sigma_T * (number density of protons now)
    
    fHe = CP%YHe/(mass_ratio_He_H*(1.d0-CP%YHe))  !n_He_tot / n_H_tot
    
    if (CP%omegan==0) then
       CP%Num_nu_massless = CP%Num_nu_massless + CP%Num_nu_massive
       CP%Num_nu_massive = 0
    end if
    
    if (.not.call_again) then
       
       call init_massive_nu(CP%omegan /=0)
       call init_background
       CP%tau0=TimeOfz(0._dl)
       ! print *, 'chi = ',  (CP%tau0 - TimeOfz(0.15_dl)) * CP%h0/100   
       last_tau0=CP%tau0
       if (WantReion) call Reionization_Init(CP%Reion,CP%ReionHist, & 
            CP%YHe, akthom, CP%tau0, FeedbackLevel)
    else
       CP%tau0=last_tau0
    end if
    
    if ( CP%NonLinear==NonLinear_Lens) then
       CP%Transfer%kmax = max(CP%Transfer%kmax, CP%Max_eta_k/CP%tau0) 
       if (FeedbackLevel > 0 .and. CP%Transfer%kmax== CP%Max_eta_k/CP%tau0) &
            write (*,*) 'max_eta_k changed to ', CP%Max_eta_k
    end if
    
    
    if (CP%closed .and. CP%tau0/CP%r >3.14) then
       if (present(error)) then
          error = 2
          return
       else
          stop 'chi >= pi in closed model not supported'
       end if
    end if
    
    if (present(error)) then
       error = 0
    else if (FeedbackLevel > 0 .and. .not. call_again) then
       write(*,'("Om_b h^2             = ",f9.6)') CP%omegab*(CP%H0/100)**2
       write(*,'("Om_c h^2             = ",f9.6)') CP%omegac*(CP%H0/100)**2
       write(*,'("Om_nu h^2            = ",f9.6)') CP%omegan*(CP%H0/100)**2
       write(*,'("Om_Lambda            = ",f9.6)') CP%omegav
       write(*,'("Om_K                 = ",f9.6)') CP%omegak
       write(*,'("Om_m (1-Om_K-Om_L)   = ",f9.6)') 1-CP%omegak-CP%omegav
       write(*,'("100 theta (CosmoMC)  = ",f9.6)') 100*CosmomcTheta()
       if (CP%Num_Nu_Massive > 0) then
          do nu_i=1, CP%Nu_mass_eigenstates 
             write(*,'(f5.2, " nu, m_nu*c^2/k_B/T_nu0   = ",f8.2, & 
                  " (m_nu = ",f6.3," eV)")') &
                  CP%nu_mass_degeneracies(nu_i), nu_masses(nu_i), & 
                  1.68e-4*nu_masses(nu_i)
          end do
       end if
    end if
    CP%chi0=rofChi(CP%tau0/CP%r)
    scale= CP%chi0*CP%r/CP%tau0  
    !e.g. changel sampling depending on approx peak spacing      
    
  end subroutine CAMBParams_Set
  !////////////////////////////////////////////////////////////////////////
  function GetTestTime()
    real(sp) GetTestTime
    real(sp) atime
    
    !GetTestTime = etime(tarray)
    !Can replace this if etime gives problems
    !Or just comment out - only used if DebugMsgs = .true.
    call cpu_time(atime)
    GetTestTime = atime
    
  end function GetTestTime
  !////////////////////////////////////////////////////////////////////////
  function rofChi(Chi) !sinh(chi) for open, sin(chi) for closed.
    real(dl) Chi,rofChi
    
    if (CP%closed) then
       rofChi=sin(chi)
    else if (CP%open) then
       rofChi=sinh(chi)
    else
       rofChi=chi
    endif
  end function rofChi  
  !////////////////////////////////////////////////////////////////////////
  function cosfunc (Chi)
    real(dl) Chi,cosfunc
    
    if (CP%closed) then
       cosfunc= cos(chi)
    else if (CP%open) then
       cosfunc=cosh(chi)
    else
       cosfunc = 1._dl
    endif
  end function cosfunc
  !////////////////////////////////////////////////////////////////////////
  function tanfunc(Chi)
    real(dl) Chi,tanfunc
    if (CP%closed) then
       tanfunc=tan(Chi)
    else if (CP%open) then
       tanfunc=tanh(Chi)
    else
       tanfunc=Chi
    end if
    
  end  function tanfunc
  !////////////////////////////////////////////////////////////////////////
  function invsinfunc(x)
    real(dl) invsinfunc,x
    
    if (CP%closed) then
       invsinfunc=asin(x)
    else if (CP%open) then
       invsinfunc=log((x+sqrt(1._dl+x**2)))  
    else
       invsinfunc = x
    endif
  end function invsinfunc
  !////////////////////////////////////////////////////////////////////////
  function f_K(x)
    real(dl) :: f_K
    real(dl), intent(in) :: x
    f_K = CP%r*rofChi(x/CP%r)
    
  end function f_K
  !////////////////////////////////////////////////////////////////////////
  function DeltaTime(a1,a2, in_tol)
    implicit none
    real(dl) DeltaTime, atol
    real(dl), intent(IN) :: a1,a2
    real(dl), optional, intent(in) :: in_tol
    real(dl) dtauda, rombint !diff of tau w.CP%r.t a and integration
    external dtauda, rombint
    
    if (present(in_tol)) then
       atol = in_tol
    else
       atol = tol/1000/exp(AccuracyBoost-1)
    end if
    DeltaTime=rombint(dtauda,a1,a2,atol)
    
  end function DeltaTime
  !////////////////////////////////////////////////////////////////////////
  function TimeOfz(z)
    implicit none
    real(dl) TimeOfz
    real(dl), intent(IN) :: z
    
    TimeOfz=DeltaTime(0._dl,1._dl/(z+1._dl))
  end function TimeOfz
  !////////////////////////////////////////////////////////////////////////
  function AngularDiameterDistance(z)
    real(dl) AngularDiameterDistance
    real(dl), intent(in) :: z
    
    AngularDiameterDistance = CP%r/(1+z)*rofchi(DeltaTime(1/(1+z),1._dl)/CP%r)
    
  end function AngularDiameterDistance
  !////////////////////////////////////////////////////////////////////////
  function dsound_da(a)
    implicit none
    real(dl) dsound_da,dtauda,a,R,cs
    external dtauda
    
    R=3.0d4*a*CP%omegab*(CP%h0/100.0d0)**2
    !R = 3*grhob*a / (4*grhog) //above is mostly within 0.2% and 
    !used for previous consistency
    cs=1.0d0/sqrt(3*(1+R))
    dsound_da=dtauda(a)*cs
    
  end function dsound_da
  !////////////////////////////////////////////////////////////////////////
  function CosmomcTheta()
    real(dl) zstar, astar, atol, rs, DA
    real(dl) CosmomcTheta
    real(dl) ombh2, omdmh2
    real(dl) rombint
    external rombint
    
    ombh2 = CP%omegab*(CP%h0/100.0d0)**2
    omdmh2 = (CP%omegac+CP%omegan)*(CP%h0/100.0d0)**2
    
    !!From Hu & Sugiyama
    zstar =  1048*(1+0.00124*ombh2**(-0.738))*(1+ &
         (0.0783*ombh2**(-0.238)/(1+39.5*ombh2**0.763)) * &
         (omdmh2+ombh2)**(0.560/(1+21.1*ombh2**1.81)))
    
    astar = 1/(1+zstar)
    atol = 1e-6
    rs = rombint(dsound_da,1d-8,astar,atol)
    DA = AngularDiameterDistance(zstar)/astar
    CosmomcTheta = rs/DA
    !       print *,'z* = ',zstar, 'r_s = ',rs, 'DA = ',DA, rs/DA
    
  end function CosmomcTheta
  !////////////////////////////////////////////////////////////////////////
end module ModelParams
!==================================================================================
module LambdaGeneral
  use precision
  implicit none
  
  real(dl)  :: w_lam = -1 !p/rho for the dark energy (assumed constant) 
  real(dl) :: cs2_lam = 1_dl 
  !comoving sound speed. Always exactly 1 for quintessence 
  !(otherwise assumed constant, though this is almost certainly unrealistic)
  real(dl) :: w0=-1,wa=0
  logical :: w_perturb = .true.
  
  
contains
  
  function Funcofw(a,deriv)
    implicit none
    
    
    real(dl) a,z,Funcofw, a0
    integer deriv,i,j
    
    
    if (a .lt. 1.d-8) then
       a0 = 1.d-8
    else
       a0 = a
    end if
    
    z=1.d0/a0-1.d0 
    

    if (deriv==0) then
       Funcofw = w0+wa*(1.d0-a0)
    else if (deriv ==1) then
       Funcofw = -wa*a0
    else if (deriv ==2) then
       Funcofw = a0**(2.d0 - 3.d0*(1.d0+w0+wa))*exp(3.d0*wa*(a0-1.d0))
    end if
    
  end  function Funcofw
   
end module LambdaGeneral
!==================================================================================
module lvalues
  use precision
  use ModelParams
  implicit none
  public
  
  Type lSamples
     integer l0
     integer l(lmax_arr)
  end Type lSamples
  
  Type(lSamples) :: lSamp
  
contains
  !////////////////////////////////////////////////////////////////////////
  subroutine initlval(lSet,max_l)
    
    ! This subroutines initializes lSet%l arrays. 
    !Other values will be interpolated.
    
    implicit none
    type(lSamples) :: lSet
    
    integer, intent(IN) :: max_l
    integer lind, lvar, step,top,bot,ls(lmax_arr)
    real(dl) AScale
    
    Ascale=scale/lSampleBoost       
    
    if (lSampleBoost >=50) then
       !just do all of them
       lind=0
       do lvar=lmin, max_l
          lind=lind+1
          ls(lind)=lvar 
       end do
       lSet%l0=lind
       lSet%l(1:lind) = ls(1:lind)
       return       
    end if
    
    lind=0
    do lvar=lmin, 10
       lind=lind+1
       ls(lind)=lvar 
    end do
    
    if (CP%AccurateReionization) then
       if (lSampleBoost > 1) then
          do lvar=11, 37,1
             lind=lind+1
             ls(lind)=lvar 
          end do
       else
          do lvar=11, 37,2
             lind=lind+1
             ls(lind)=lvar 
          end do
       end if
       
       step = max(nint(5*Ascale),2)           
       bot=40
       top=bot + step*10
    else
       
       if (lSampleBoost >1) then
          do lvar=11, 15
             lind=lind+1
             ls(lind)=lvar 
          end do
       else
          lind=lind+1
          ls(lind)=12
          lind=lind+1
          ls(lind)=15
       end if
       step = max(nint(10*Ascale),3)           
       bot=15+max(step/2,2)
       top=bot + step*7
    end if
    
    do lvar=bot, top, step
       lind=lind+1
       ls(lind)=lvar          
    end do
    
    step=max(nint(20*Ascale),4)
    bot=ls(lind)+step
    top=bot+step*2
    
    do lvar = bot,top,step 
       lind=lind+1
       ls(lind)=lvar
    end do
    
    if (ls(lind)>=max_l) then
       do lvar=lind,1,-1
          if (ls(lvar)<=max_l) exit  
       end do
       lind=lvar
       if (ls(lind)<max_l) then
          lind=lind+1
          ls(lind)=max_l
       end if
    else
       
       step=max(nint(25*Ascale),4)
       !Get EE right around l=200 by putting extra point at 175
       bot=ls(lind)+step
       top=bot+step
       
       do lvar = bot,top,step 
          lind=lind+1
          ls(lind)=lvar
       end do
       
       
       if (ls(lind)>=max_l) then
          do lvar=lind,1,-1
             if (ls(lvar)<=max_l) exit  
          end do
          lind=lvar
          if (ls(lind)<max_l) then
             lind=lind+1
             ls(lind)=max_l
          end if
       else
          
          step=max(nint(50*Ascale),7)
          bot=ls(lind)+step
          top=min(5000,max_l)
          
          do lvar = bot,top,step
             lind=lind+1
             ls(lind)=lvar
          end do
          
          if (max_l > 5000) then
             !Should be pretty smooth or tiny out here   
             step=max(nint(400*Ascale),50)
             lvar = ls(lind)
            
             do
              lvar = lvar + step
              if (lvar > max_l) exit
              lind=lind+1
              ls(lind)=lvar
              step = nint(step*1.5) !log spacing
             end do

         end if

         if (ls(lind) /=max_l) then          
           lind=lind+1
           ls(lind)=max_l
         end if
         if (.not. CP%flat) ls(lind-1)=int(max_l+ls(lind-2))/2
         !Not in CP%flat case so interpolation table is the same 
         !when using lower l_max
      end if
   end if
   lSet%l0=lind
   lSet%l(1:lind) = ls(1:lind)
   
 end subroutine initlval
 !////////////////////////////////////////////////////////////////////////
 subroutine InterpolateClArr(lSet,iCl, all_Cl, max_ind)
   type (lSamples), intent(in) :: lSet        
   real(dl), intent(in) :: iCl(*)
   real(dl), intent(out):: all_Cl(lmin:*)
   integer, intent(in) :: max_ind
   integer il,llo,lhi, xi
   real(dl) ddCl(lSet%l0)
   real(dl) xl(lSet%l0)
   
   real(dl) a0,b0,ho
   real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
   
   if (max_ind > lSet%l0) stop 'Wrong max_ind in InterpolateClArr'
   
   xl = real(lSet%l(1:lSet%l0),dl)
   call spline(xl,iCl(1),max_ind,cllo,clhi,ddCl(1))
   
   llo=1
   do il=lmin,lSet%l(max_ind)
      xi=il
      if ((xi > lSet%l(llo+1)).and.(llo < max_ind)) then
         llo=llo+1
      end if
      lhi=llo+1
      ho=lSet%l(lhi)-lSet%l(llo)
      a0=(lSet%l(lhi)-xi)/ho
      b0=(xi-lSet%l(llo))/ho
      
      all_Cl(il) = a0*iCl(llo)+ b0*iCl(lhi)+((a0**3-a0)* ddCl(llo) &
           +(b0**3-b0)*ddCl(lhi))*ho**2/6
      
   end do
   
 end subroutine InterpolateClArr
 !////////////////////////////////////////////////////////////////////////
end module lvalues
!==================================================================================
MODULE ModelData
  use precision
  use ModelParams
  use InitialPower
  use lValues
  use Ranges
  use AMlUtils
  implicit none
  public
  
  Type ClTransferData
     !Cl transfer function variables
     !values of q for integration over q to get C_ls
     Type (lSamples) :: ls ! scalar and tensor l that are computed
     integer :: NumSources 
     !Changes -scalars:  2 for just CMB, 3 for lensing
     !- tensors: T and E and phi (for lensing), and T, E, B respectively
     
     Type (Regions) :: q
     
     real(dl), dimension(:,:,:), pointer :: Delta_p_l_k => NULL() 
     
  end Type ClTransferData
  
  Type(ClTransferData), save, target :: CTransScal, CTransTens, CTransVec
  
  !Computed output power spectra data
  
  integer, parameter :: C_Temp = 1, C_E = 2, C_Cross =3, C_Phi = 4
  integer, parameter :: C_PhiTemp = 5, C_PhiE=6
  integer :: C_last = C_PhiE
  integer, parameter :: CT_Temp =1, CT_E = 2, CT_B = 3, CT_Cross=  4

  real(dl), dimension (:,:,:), allocatable :: Cl_scalar, Cl_tensor, Cl_vector
  !Indices are Cl_xxx( l , intial_power_index, Cl_type)
  !where Cl_type is one of the above constants
  
  !The following are set only if doing lensing
  integer lmax_lensed !Only accurate to rather less than this 
  real(dl) , dimension (:,:,:), allocatable :: Cl_lensed
  !Cl_lensed(l, power_index, Cl_type) are the interpolated Cls
  
  real(dl), dimension (:), allocatable ::  COBElikelihoods,COBE_scales
  !Set by COBEnormalize if using outCOBE
contains
  !////////////////////////////////////////////////////////////////////////
  subroutine Init_ClTransfer(CTrans)
    !Need to set the Ranges array q before calling this
    Type(ClTransferData) :: CTrans
    integer st
    
    deallocate(CTrans%Delta_p_l_k, STAT = st)
    call Ranges_getArray(CTrans%q, .true.)
    
    allocate(CTrans%Delta_p_l_k(CTrans%NumSources,min(max_bessels_l_index, & 
         CTrans%ls%l0), CTrans%q%npoints))
    CTrans%Delta_p_l_k = 0
    
  end subroutine Init_ClTransfer
  !////////////////////////////////////////////////////////////////////////
  subroutine Free_ClTransfer(CTrans)
    Type(ClTransferData) :: CTrans
    integer st
    
    deallocate(CTrans%Delta_p_l_k, STAT = st)
    nullify(CTrans%Delta_p_l_k)
    call Ranges_Free(CTrans%q)
    
  end subroutine Free_ClTransfer
  !////////////////////////////////////////////////////////////////////////
  subroutine Init_Cls
    
    if (CP%WantScalars) then
       if (allocated(Cl_scalar)) deallocate(Cl_scalar)
       allocate(Cl_scalar(lmin:CP%Max_l, CP%InitPower%nn, C_Temp:C_last))
       Cl_scalar = 0
    end if
    
    if (CP%WantVectors) then
       if (allocated(Cl_vector)) deallocate(Cl_vector)
       allocate(Cl_vector(lmin:CP%Max_l, CP%InitPower%nn, CT_Temp:CT_Cross))
       Cl_vector = 0
    end if
    
    
    if (CP%WantTensors) then
       if (allocated(Cl_tensor)) deallocate(Cl_tensor)
       allocate(Cl_tensor(lmin:CP%Max_l_tensor, CP%InitPower%nn, CT_Temp:CT_Cross))
       Cl_tensor = 0
    end if
    
  end subroutine Init_Cls
  !////////////////////////////////////////////////////////////////////////
  subroutine output_cl_files(ScalFile,TensFile,TotFile,LensFile,LensTotFile,factor)
    implicit none
    integer in,il
    character(LEN=*) ScalFile, TensFile, TotFile, LensFile, LensTotFile
    real(dl), intent(in), optional :: factor
    real(dl) fact
    integer last_C
    
    if (present(factor)) then
       fact = factor
    else
       fact =1
    end if
    
    if (CP%WantScalars .and. ScalFile /= '') then
       last_C=min(C_PhiTemp,C_last)
       open(unit=fileio_unit,file=ScalFile,form='formatted',status='replace')
       do in=1,CP%InitPower%nn
          do il=lmin,min(10000,CP%Max_l)
             write(fileio_unit,trim(numcat('(1I6,',last_C))//'E15.5)')il, & 
                  fact*Cl_scalar(il,in,C_Temp:last_C)
          end do
          do il=10100,CP%Max_l, 100
             write(fileio_unit,trim(numcat('(1E15.5,',last_C))//'E15.5)') real(il),&
                  fact*Cl_scalar(il,in,C_Temp:last_C)
          end do
       end do
       close(fileio_unit)
    end if
    
    if (CP%WantTensors .and. TensFile /= '') then
       open(unit=fileio_unit,file=TensFile,form='formatted',status='replace')
       do in=1,CP%InitPower%nn
          do il=lmin,CP%Max_l_tensor
             write(fileio_unit,'(1I6,4E15.5)')il, & 
                  fact*Cl_tensor(il, in, CT_Temp:CT_Cross)
          end do
       end do
       close(fileio_unit)
    end if
    
    if (CP%WantTensors .and. CP%WantScalars .and. TotFile /= '') then
       open(unit=fileio_unit,file=TotFile,form='formatted',status='replace')
       do in=1,CP%InitPower%nn
          do il=lmin,CP%Max_l_tensor
             
             write(fileio_unit,'(1I6,4E15.5)')il, & 
                  fact*(Cl_scalar(il, in, C_Temp:C_E)+ Cl_tensor(il,in, C_Temp:C_E)),&
                  fact*Cl_tensor(il,in, CT_B), & 
                  fact*(Cl_scalar(il, in, C_Cross) + Cl_tensor(il, in, CT_Cross))
          end do
          do il=CP%Max_l_tensor+1,CP%Max_l
             write(fileio_unit,'(1I6,4E15.5)')il ,fact*Cl_scalar(il,in,C_Temp:C_E), & 
                  0._dl, fact*Cl_scalar(il,in,C_Cross)
          end do
       end do
       close(fileio_unit)
    end if
    
    if (CP%WantScalars .and. CP%DoLensing .and. LensFile /= '') then
       open(unit=fileio_unit,file=LensFile,form='formatted',status='replace')
       do in=1,CP%InitPower%nn
          do il=lmin, lmax_lensed
             write(fileio_unit,'(1I6,4E15.5)')il, & 
                  fact*Cl_lensed(il, in, CT_Temp:CT_Cross)
          end do
       end do
       close(fileio_unit)      
    end if
    
    
    if (CP%WantScalars .and. CP%WantTensors .and. CP%DoLensing & 
         .and. LensTotFile /= '') then
       open(unit=fileio_unit,file=LensTotFile,form='formatted',status='replace')
       do in=1,CP%InitPower%nn
          do il=lmin,min(CP%Max_l_tensor,lmax_lensed)
             write(fileio_unit,'(1I6,4E15.5)')il, & 
                  fact*(Cl_lensed(il, in, CT_Temp:CT_Cross) & 
                  + Cl_tensor(il,in, CT_Temp:CT_Cross))
          end do
          do il=min(CP%Max_l_tensor,lmax_lensed)+1,lmax_lensed
             write(fileio_unit,'(1I6,4E15.5)')il, & 
                  fact*Cl_lensed(il, in, CT_Temp:CT_Cross)
          end do
       end do
       
    end if
  end subroutine output_cl_files
  !////////////////////////////////////////////////////////////////////////
  subroutine output_lens_pot_files(LensPotFile, factor)
    !Write out L TT EE BB TE PP PT PE where P is the lensing potential, all unlensed  
    !This input supported by LensPix from 2010
    implicit none
    integer in,il
    real(dl), intent(in), optional :: factor
    real(dl) fact, scale, BB, TT, TE, EE
    character(LEN=*) LensPotFile
    !output file of dimensionless [l(l+1)]^2 C_phi_phi/2pi 
    !and [l(l+1)]^(3/2) C_phi_T/2pi 
    !This is the format used by Planck_like but original LensPix uses 
    !scalar_output_file.
    
    !(Cl_scalar and scalar_output_file numbers are instead l^4 C_phi and l^3 C_phi 
    ! - for historical reasons) 
    
    if (present(factor)) then
       fact = factor
    else
       fact =1
    end if
    
    if (CP%WantScalars .and. CP%DoLensing .and. LensPotFile/='') then
       
       open(unit=fileio_unit,file=LensPotFile,form='formatted',status='replace')
       do in=1,CP%InitPower%nn
          do il=lmin,min(10000,CP%Max_l)
             
             TT = Cl_scalar(il, in, C_Temp)
             EE = Cl_scalar(il, in, C_E)
             TE = Cl_scalar(il, in, C_Cross)              
             if (CP%WantTensors .and. il <= CP%Max_l_tensor) then
                TT= TT+Cl_tensor(il,in, CT_Temp)
                EE= EE+Cl_tensor(il,in, CT_E)
                TE= TE+Cl_tensor(il,in, CT_Cross)
                BB= Cl_tensor(il,in, CT_B)               
             else
                BB=0
             end if
             scale = (real(il+1)/il)**2/OutputDenominator 
             !Factor to go from old l^4 factor to new
               
             write(fileio_unit,'(1I6,7E15.5)') il , fact*TT, fact*EE, & 
                  fact*BB, fact*TE, scale*Cl_scalar(il,in,C_Phi),&
                  (real(il+1)/il)**1.5/OutputDenominator*sqrt(fact) & 
                  *Cl_scalar(il,in,C_PhiTemp:C_PhiE)
                   
          end do
          do il=10100,CP%Max_l, 100
             scale = (real(il+1)/il)**2/OutputDenominator
             write(fileio_unit,'(1E15.5,7E15.5)') real(il), & 
                  fact*Cl_scalar(il,in,C_Temp:C_E),0., & 
                  fact*Cl_scalar(il,in,C_Cross), &
                  scale*Cl_scalar(il,in,C_Phi),&
                  (real(il+1)/il)**1.5/OutputDenominator*sqrt(fact) & 
                  *Cl_scalar(il,in,C_PhiTemp:C_PhiE)
          end do
       end do
       close(fileio_unit)
    end if
  end subroutine output_lens_pot_files
  !////////////////////////////////////////////////////////////////////////
  subroutine output_veccl_files(VecFile, factor)
    implicit none
    integer in,il
    character(LEN=*) VecFile
    real(dl), intent(in), optional :: factor
    real(dl) fact
    
    if (present(factor)) then
       fact = factor
    else
       fact =1
    end if
    
    if (CP%WantVectors .and. VecFile /= '') then
       open(unit=fileio_unit,file=VecFile,form='formatted',status='replace')
       do in=1,CP%InitPower%nn
          do il=lmin,CP%Max_l
             write(fileio_unit,'(1I5,4E15.5)')il, & 
                  fact*Cl_vector(il, in, CT_Temp:CT_Cross)
          end do
       end do
       
       close(fileio_unit)
    end if
    
  end subroutine output_veccl_files
  !////////////////////////////////////////////////////////////////////////
  subroutine output_COBElikelihood
    integer in
    do in=1, CP%InitPower%nn
       write(*,*)'COBE Likelihood relative to CP%flat=',COBElikelihoods(in)
    end do
  end  subroutine output_COBElikelihood
  !////////////////////////////////////////////////////////////////////////
  subroutine NormalizeClsAtL(lnorm)
    implicit none
    integer, intent(IN) :: lnorm
    integer in
    real(dl) Norm
    
    do in=1,CP%InitPower%nn
       
       if (CP%WantScalars) then
          Norm=1/Cl_scalar(lnorm,in, C_Temp)
          Cl_scalar(lmin:CP%Max_l, in, C_Temp:C_Cross) = Cl_scalar(lmin:CP%Max_l, & 
               in, C_Temp:C_Cross) * Norm
       end if
       
       if (CP%WantTensors) then
          if (.not.CP%WantScalars) Norm = 1/Cl_tensor(lnorm,in, C_Temp)
          !Otherwise Norm already set correctly
          Cl_tensor(lmin:CP%Max_l_tensor, in, CT_Temp:CT_Cross) =  &
               Cl_tensor(lmin:CP%Max_l_tensor, in, CT_Temp:CT_Cross) * Norm
       end if
    end do
    
  end  subroutine NormalizeClsAtL
  !////////////////////////////////////////////////////////////////////////
  subroutine COBEnormalize
    use precision
    use ModelParams
    
    integer in
    real(dl) xlog10
    real(dl) c10, d1,d2,d3,d4,d5,d6,d7, xlogl, COBE_scale
    real(dl) x1, x2,x3,x4,x5,x6,x7,sy,s,sx,sxy,sxx,delt,d1pr,d1ppr
    real(dl) Ctot(lmin:20)
    
    if (allocated(COBElikelihoods)) deallocate(COBElikelihoods)
    if (allocated(COBE_scales)) deallocate(COBE_scales)
    allocate(COBElikelihoods(CP%InitPower%nn))
    allocate(COBE_scales(CP%InitPower%nn))
    
    xlog10=log(10._dl)

    ! COBE normalization
    ! fit the spectrum to a quadratic around C_10 with equal weights in logl
    
    do in=1,CP%InitPower%nn
       
       if (CP%WantTensors) then
          Ctot =  Cl_tensor(lmin:20, in, C_Temp)
       else
          Ctot = 0
       end if
       if (CP%WantScalars) then
          Ctot=Ctot + Cl_scalar(lmin:20, in, C_Temp)
          
       end if
       c10=Ctot(10)
       
       d1=(Ctot(3))/c10-1._dl
       d2=(Ctot(4))/c10-1._dl
       d3=(Ctot(6))/c10-1._dl
       d4=(Ctot(8))/c10-1._dl
       d5=(Ctot(12))/c10-1._dl
       d6=(Ctot(15))/c10-1._dl
       d7=(Ctot(20))/c10-1._dl
       
       
       x1=log(3._dl)/xlog10-1._dl
       x2=log(4._dl)/xlog10-1._dl
       x3=log(6._dl)/xlog10-1._dl
       x4=log(8._dl)/xlog10-1._dl
       x5=log(12._dl)/xlog10-1._dl
       x6=log(15._dl)/xlog10-1._dl
       x7=log(20._dl)/xlog10-1._dl
       sy=x1*d1+x2*d2+x3*d3+x4*d4+x5*d5+x6*d6+x7*d7
       s=x1*x1+x2*x2+x3*x3+x4*x4+x5*x5+x6*x6+x7*x7
       sx=x1**3+x2**3+x3**3+x4**3+x5**3+x6**3+x7**3
       sxy=x1**2*d1+x2**2*d2+x3**2*d3+x4**2*d4+ &
            x5**2*d5+x6**2*d6+x7**2*d7
       sxx=x1**4+x2**4+x3**4+x4**4+x5**4+x6**4+x7**4
       delt=s*sxx-sx*sx
       d1pr=(sxx*sy-sx*sxy)/delt
       d1ppr=2._dl*(s*sxy-sx*sy)/delt
       
       ! Bunn and White fitting formula
       c10=(0.64575d0+0.02282d0*d1pr+0.01391d0*d1pr*d1pr &
            -0.01819d0*d1ppr-0.00646d0*d1pr*d1ppr &
            +0.00103d0*d1ppr*d1ppr)/c10
       ! logl
       xlogl=-0.01669d0+1.19895d0*d1pr-0.83527d0*d1pr*d1pr &
            -0.43541d0*d1ppr-0.03421d0*d1pr*d1ppr &
            +0.01049d0*d1ppr*d1ppr
       ! write(*,*)'COBE Likelihood relative to CP%flat=',exp(xlogl)
       COBElikelihoods(in) = exp(xlogl)
       
       ! density power spectrum normalization;
       
       COBE_scale=c10/OutputDenominator*1.1d-9
       COBE_scales(in)=COBE_scale
       
!!$!delta^2 = k^4*(tf)^2*ScalarPower(k,in)*COBE_scale where (tf) is output 
       !in the transfer function file
!!$!delta^2 = 4*pi*k^3 P(k)
       
       ! C_l normalization; output l(l+1)C_l/twopi
       c10=c10*2.2d-9/fourpi
       
       if (CP%WantScalars) Cl_scalar(lmin:CP%Max_l, in, C_Temp:C_last) = &
            Cl_scalar(lmin:CP%Max_l, in, C_Temp:C_last)*c10
       if (CP%WantTensors) Cl_tensor(lmin:CP%Max_l_tensor, in, CT_Temp:CT_Cross) = &
            Cl_tensor(lmin:CP%Max_l_tensor, in, CT_Temp:CT_Cross)*c10
       
    end do !in
  end subroutine COBEnormalize
  !////////////////////////////////////////////////////////////////////////
  subroutine ModelData_Free

    call Free_ClTransfer(CTransScal)
    call Free_ClTransfer(CTransVec)
    call Free_ClTransfer(CTransTens)
    if (allocated(Cl_vector)) deallocate(Cl_vector)
    if (allocated(Cl_tensor)) deallocate(Cl_tensor)
    if (allocated(Cl_scalar)) deallocate(Cl_scalar)
    if (allocated(Cl_lensed)) deallocate(Cl_lensed)
    if (allocated(COBElikelihoods)) deallocate(COBElikelihoods)
    if (allocated(COBE_scales)) deallocate(COBE_scales) 
    
  end subroutine ModelData_Free
  !////////////////////////////////////////////////////////////////////////
end module ModelData
!==================================================================================
module MassiveNu
  use precision
  use ModelParams
  implicit none
  private 
  
  real(dl), parameter  :: const  = 7._dl/120*pi**4 ! 5.68219698_dl
  !const = int q^3 F(q) dq = 7/120*pi^4
  real(dl), parameter  :: const2 = 5._dl/7/pi**2   !0.072372274_dl
  real(dl), parameter  :: zeta3  = 1.2020569031595942853997_dl
  real(dl), parameter  :: zeta5  = 1.0369277551433699263313_dl
  real(dl), parameter  :: zeta7  = 1.0083492773819228268397_dl
  
  integer, parameter  :: nrhopn=2000  
  real(dl), parameter :: am_min = 0.01_dl  !0.02_dl
  !smallest a*m_nu to integrate distribution function rather than using series
  real(dl), parameter :: am_max = 600._dl 
  !max a*m_nu to integrate
  
  real(dl),parameter  :: am_minp=am_min*1.1
  real(dl), parameter :: am_maxp=am_max*0.9
  
  real(dl) dlnam
  
  real(dl), dimension(:), allocatable ::  r1,p1,dr1,dp1,ddr1
  
  !Sample for massive neutrino momentum
  !These settings appear to be OK for P_k accuate at 1e-3 level
  integer, parameter :: nqmax0=80 !maximum array size of q momentum samples 
  real(dl) :: nu_q(nqmax0), nu_int_kernel(nqmax0)
  
  integer nqmax !actual number of q modes evolves
  
  public const,Nu_Init,Nu_background, Nu_rho, Nu_drho,  nqmax0, nqmax, &
       nu_int_kernel, nu_q
contains
  !////////////////////////////////////////////////////////////////////////
  subroutine Nu_init
    
    !  Initialize interpolation tables for massive neutrinos.
    !  Use cubic splines interpolation of log rhonu and pnu vs. log a*m.
    
    integer i
    real(dl) dq,dlfdlq, q, am, rhonu,pnu
    real(dl) spline_data(nrhopn)
    
    !  nu_masses=m_nu(i)*c**2/(k_B*T_nu0).
    !  Get number density n of neutrinos from
    !  rho_massless/n = int q^3/(1+e^q) / int q^2/(1+e^q)=7/180 pi^4/Zeta(3)
    !  then m = Omega_nu/N_nu rho_crit /n
    !  Error due to velocity < 1e-5
    
    do i=1, CP%Nu_mass_eigenstates 
       nu_masses(i)=const/(1.5d0*zeta3)*grhom/grhor & 
            *CP%omegan*CP%Nu_mass_fractions(i)/CP%Nu_mass_degeneracies(i)
    end do
    
    if (allocated(r1)) return
    allocate(r1(nrhopn),p1(nrhopn),dr1(nrhopn),dp1(nrhopn),ddr1(nrhopn))
    
    
    nqmax=3
    if (AccuracyBoost >1) nqmax=4
    if (AccuracyBoost >2) nqmax=5
    if (AccuracyBoost >3) nqmax=nint(AccuracyBoost*10) 
    !note this may well be worse than the 5 optimized points
    
    if (nqmax > nqmax0) call MpiStop('Nu_Init: qmax > nqmax0')
    
    !We evolve evolve 4F_l/dlfdlq(i), so kernel includes dlfdlnq factor
    !Integration scheme gets (Fermi-Dirac thing)*q^n exact,for n=-4, -2..2
    !see CAMB notes
    if (nqmax==3) then
       !Accurate at 2e-4 level
       nu_q(1:3) = (/0.913201, 3.37517, 7.79184/)
       nu_int_kernel(1:3) = (/0.0687359, 3.31435, 2.29911/)
       
    else if (nqmax==4) then
       !This seems to be very accurate (limited by other numerics)
       nu_q(1:4) = (/0.7, 2.62814, 5.90428, 12.0/)
       nu_int_kernel(1:4) = (/0.0200251, 1.84539, 3.52736, 0.289427/)
       
    else if (nqmax==5) then
       !exact for n=-4,-2..3 
       !This seems to be very accurate (limited by other numerics)
       nu_q(1:5) = (/0.583165, 2.0, 4.0, 7.26582, 13.0/)  
       nu_int_kernel(1:5) = (/0.0081201, 0.689407, 2.8063, 2.05156, 0.126817/) 
       
    else
       dq = (12 + nqmax/5)/real(nqmax)
       do i=1,nqmax
          q=(i-0.5d0)*dq
          nu_q(i) = q 
          dlfdlq=-q/(1._dl+exp(-q))
          nu_int_kernel(i)=dq*q**3/(exp(q)+1._dl) * (-0.25_dl*dlfdlq) 
          !now evolve 4F_l/dlfdlq(i)
          
       end do
    end if
    nu_int_kernel=nu_int_kernel/const
    
    dlnam=-(log(am_min/am_max))/(nrhopn-1)
    
    
    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC) &
    !$OMP & PRIVATE(am, rhonu,pnu) 
    do i=1,nrhopn
       am=am_min*exp((i-1)*dlnam)
       call nuRhoPres(am,rhonu,pnu)
       r1(i)=log(rhonu)
       p1(i)=log(pnu)
    end do
    !$OMP END PARALLEL DO
    
    
    call splini(spline_data,nrhopn)
    call splder(r1,dr1,nrhopn,spline_data)
    call splder(p1,dp1,nrhopn,spline_data)
    call splder(dr1,ddr1,nrhopn,spline_data)       
      
  end subroutine Nu_init
  !////////////////////////////////////////////////////////////////////////
  subroutine nuRhoPres(am,rhonu,pnu)
    !  Compute the density and pressure of one eigenstate of massive neutrinos,
    !  in units of the mean density of one flavor of massless neutrinos.
    
    real(dl),  parameter :: qmax=30._dl
    integer, parameter :: nq=100
    real(dl) dum1(nq+1),dum2(nq+1)
    real(dl), intent(in) :: am
    real(dl), intent(out) ::  rhonu,pnu
    integer i
    real(dl) q,aq,v,aqdn,adq
    
    
    !  q is the comoving momentum in units of k_B*T_nu0/c.
    !  Integrate up to qmax and then use asymptotic expansion for remainder.
    adq=qmax/nq
    dum1(1)=0._dl
    dum2(1)=0._dl
    do  i=1,nq
       q=i*adq
       aq=am/q
       v=1._dl/sqrt(1._dl+aq*aq)
       aqdn=adq*q*q*q/(exp(q)+1._dl)
       dum1(i+1)=aqdn/v
       dum2(i+1)=aqdn*v
    end do
    call splint(dum1,rhonu,nq+1)
    call splint(dum2,pnu,nq+1)
    !  Apply asymptotic corrrection for q>qmax and normalize by relativistic
    !  energy density.
    rhonu=(rhonu+dum1(nq+1)/adq)/const
    pnu=(pnu+dum2(nq+1)/adq)/const/3._dl
    
  end subroutine nuRhoPres
  !////////////////////////////////////////////////////////////////////////
  subroutine Nu_background(am,rhonu,pnu)
    use precision
    use ModelParams
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu, pnu
    
    !  Compute massive neutrino density and pressure in units of the mean
    !  density of one eigenstate of massless neutrinos.  Use cubic splines to
    !  interpolate from a table.
    
    real(dl) d
    integer i
    
    if (am <= am_minp) then
       rhonu=1._dl + const2*am**2  
       pnu=(2-rhonu)/3._dl
       return
    else if (am >= am_maxp) then
       rhonu = 3/(2*const)*(zeta3*am + (15*zeta5)/2/am)
       pnu = 900._dl/120._dl/const*(zeta5-63._dl/4*Zeta7/am**2)/am
       return
    end if
    
    d=log(am/am_min)/dlnam+1._dl
    i=int(d)
    d=d-i
    
    !  Cubic spline interpolation.
    rhonu=r1(i)+d*(dr1(i)+d*(3._dl*(r1(i+1)-r1(i))-2._dl*dr1(i) &
         -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2._dl*(r1(i)-r1(i+1)))))
    pnu=p1(i)+d*(dp1(i)+d*(3._dl*(p1(i+1)-p1(i))-2._dl*dp1(i) &
         -dp1(i+1)+d*(dp1(i)+dp1(i+1)+2._dl*(p1(i)-p1(i+1)))))
    rhonu=exp(rhonu)
    pnu=exp(pnu)
    
  end subroutine Nu_background
  !////////////////////////////////////////////////////////////////////////
  subroutine Nu_rho(am,rhonu)
    use precision
    use ModelParams
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu
    
    !  Compute massive neutrino density in units of the mean
    !  density of one eigenstate of massless neutrinos.  Use cubic splines to
    !  interpolate from a table.
    
    real(dl) d
    integer i
    
    if (am <= am_minp) then
       rhonu=1._dl + const2*am**2  
       return
    else if (am >= am_maxp) then
       rhonu = 3/(2*const)*(zeta3*am + (15*zeta5)/2/am)
       return
    end if
    
    d=log(am/am_min)/dlnam+1._dl
    i=int(d)
    d=d-i
    
    !  Cubic spline interpolation.
    rhonu=r1(i)+d*(dr1(i)+d*(3._dl*(r1(i+1)-r1(i))-2._dl*dr1(i) &
         -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2._dl*(r1(i)-r1(i+1)))))
    rhonu=exp(rhonu)
  end subroutine Nu_rho
  !////////////////////////////////////////////////////////////////////////
  function Nu_drho(am,adotoa,rhonu) result (rhonudot)
    use precision
    use ModelParams
    
    !  Compute the time derivative of the mean density in massive neutrinos
    !  and the shear perturbation.
    real(dl) adotoa,rhonu,rhonudot
    real(dl) d
    real(dl), intent(IN) :: am
    integer i
    
    if (am< am_minp) then
       
       rhonudot = 2*const2*am**2*adotoa
       
    else if (am>am_maxp) then
       
       rhonudot = 3/(2*const)*(zeta3*am - (15*zeta5)/2/am)*adotoa
       
    else
       
       d=log(am/am_min)/dlnam+1._dl
       i=int(d)
       d=d-i
       !  Cubic spline interpolation for rhonudot.
       rhonudot=dr1(i)+d*(ddr1(i)+d*(3._dl*(dr1(i+1)-dr1(i)) &
            -2._dl*ddr1(i)-ddr1(i+1)+d*(ddr1(i)+ddr1(i+1) &
            +2._dl*(dr1(i)-dr1(i+1)))))
       
       rhonudot=rhonu*adotoa*rhonudot/dlnam
    end if
    
  end function Nu_drho
  !////////////////////////////////////////////////////////////////////////
end module MassiveNu
!=====================================================================================
! wrapper function to avoid cirular module references
subroutine init_massive_nu(has_massive_nu)
  use MassiveNu
  use ModelParams
  implicit none
  logical, intent(IN) :: has_massive_nu
  
  if (has_massive_nu) then
     call Nu_Init  
  else
     nu_masses = 0
  end if
end subroutine init_massive_nu
!====================================================================================
!===== LARGELY MODIFIED BY S.S. =====
module Transfer
  use ModelData
  use constants
  use PersonalUtils
  use LambdaGeneral
  implicit none
  public
  integer, parameter :: Transfer_kh =1, Transfer_cdm=2, Transfer_b=3
  integer, parameter :: Transfer_g=4, Transfer_r=5, Transfer_nu = 6
  !massless and massive neutrino
  integer, parameter :: Transfer_tot=7
  
  integer, parameter :: Transfer_max = Transfer_tot

  !For 1-loop Standard Perturbation Theory corrections
  integer, parameter :: id_22dd=1, id_22dv=2, id_22vv=3
  integer, parameter :: id_b2dd=4, id_b2dv=5, id_b22 =6
  integer, parameter :: id_13dd=7, id_13dv=8, id_13vv=9
  
  integer, parameter :: id_kh=1, id_PkLin=2
  integer, parameter :: id_PkNL_dd=3, id_PkNL_dv=4, id_PkNL_vv=5 
  integer, parameter :: id_Pkb2_dd=6, id_Pkb2_dv=7, id_Pkb22=8

  !For RSD A & B correction terms
  integer, parameter :: id_A11=1, id_A12=2, id_A22=3, id_A23=4, id_A33=5
  integer, parameter :: id_B111=1, id_B112=2, id_B121=3, id_B122=4
  integer, parameter :: id_B211=5, id_B212=6, id_B221=7, id_B222=8
  integer, parameter :: id_B312=9, id_B321=10,id_B322=11,id_B422=12
  integer, parameter :: numid_RSDTNS = id_A33+id_B422

  !For local & non-local bias corrections
  integer, parameter :: id_bs2d=1, id_bs2v=2, id_b3=5
  integer, parameter :: id_b2s2=3, id_bs22=4
  integer, parameter :: numid_bias=id_b3

  !For A01 corrections
  integer, parameter :: id_Bb1=1, id_Bb2=2, id_Bbs2=3
  !integer, parameter :: id_tBb2=4, id_tBbs2=5
  integer, parameter :: numid_A01=id_Bbs2

  !For output monopole
  integer, parameter :: id_mu0dd=2,id_mu0b2d=3,id_mu0b22=4
  integer, parameter :: id_mu2dv=5,id_mu2b2v=6,id_mu2AB=7
  integer, parameter :: id_mu4vv=8,id_mu4AB=9
  integer, parameter :: id_mu6AB=10,id_mu8AB=11
  
  logical :: transfer_interp_matterpower  = .true. !output regular grid in log k
  !set to false to output calculated values for later interpolation
  
  integer :: transfer_power_var = Transfer_tot 
  !What to use to calulcate the output matter power spectrum and sigma_8
  !Transfer_tot uses total matter perturbation
  
  Type MatterTransferData
     !Computed data
     integer   ::  num_q_trans   !    number of steps in k for transfer calculation
     real(dl), dimension (:), pointer :: q_trans => NULL() 
     real(dl), dimension (:,:), pointer ::  sigma_8 => NULL() 
     real, dimension(:,:,:), pointer :: TransferData => NULL() 
     !TransferData(entry,k_index,z_index) for entry=Tranfer_kh.. Transfer_tot
  end Type MatterTransferData
  
  Type MatterPowerData
     !everything is a function of k/h
     integer   ::  num_k, num_z          
     real(dl), dimension(:), pointer :: log_kh, redshifts => NULL() 
     !matpower is log(Pk)
     real(dl), dimension(:,:), pointer :: matpower, ddmat => NULL() 
     real(dl), dimension(:,:), pointer :: matpower_cb, ddmat_cb => NULL()
     !real(dl), dimension(:,:), allocatable :: matpower_cb, ddmat_cb
     !Pk_total = fcb**2*Pk_cb + Pk_rest
     !where Pk_rest = 2*fcb*fnu*P_cbnu + fnu**2*Pnu stays linear 
     !because of linear neutrino approximation
     real(dl), dimension(:,:), pointer :: matpower_rest => NULL() 
     !real(dl), dimension(:,:,:), pointer :: outputPk => NULL()
     !if NonLinear, nonlin_ratio =  sqrt(P_nonlinear/P_linear)
     !function of k and redshift NonLinearScaling(k_index,z_index)         
     real(dl), dimension(:,:), pointer :: nonlin_ratio => NULL() 
  end Type MatterPowerData
  
  Type(MatterTransferData), save :: MT
  
contains
  !////////////////////////////////////////////////////////////////////////
  subroutine Transfer_GetMatterPowerData_SPT(MTrans, PK_data, in, itf_only)
    !Does *NOT* include non-linear corrections
    !Get total matter power spectrum in units of 
    !(h Mpc^{-1})^3 ready for interpolation.
    !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    !We are assuming that Cls are generated so any baryonic wiggles are 
    !well sampled and that matter power
    !sepctrum is generated to beyond the CMB k_max

    Type(MatterTransferData), intent(in) :: MTrans
    Type(MatterPowerData) :: PK_data
    integer, intent(in) :: in
    integer, intent(in), optional :: itf_only
    real(dl) h, kh, k, fcb, fnu, Tk_cb, Tk_nu, Tk_tot
    integer ik
    integer nz,itf, itf_start, itf_end
          
    if (present(itf_only)) then
       itf_start=itf_only
       itf_end = itf_only
       nz = 1
    else
       itf_start=1
       nz= size(MTrans%TransferData,3)
       itf_end = nz
    end if

    PK_data%num_k = MTrans%num_q_trans
    PK_Data%num_z = nz
   
    allocate(PK_data%matpower(PK_data%num_k,nz))
    allocate(PK_data%ddmat(PK_data%num_k,nz))
    allocate(PK_data%matpower_cb(PK_data%num_k,nz))
    allocate(PK_data%ddmat_cb(PK_data%num_k,nz))
    allocate(PK_data%matpower_rest(PK_data%num_k,nz))
    allocate(PK_data%nonlin_ratio(PK_data%num_k,nz))
    allocate(PK_data%log_kh(PK_data%num_k))
    allocate(PK_data%redshifts(nz))

    PK_data%redshifts = CP%Transfer%Redshifts(itf_start:itf_end)

    h = CP%H0/100

    fcb = (CP%omegab+CP%omegac)/(CP%omegab+CP%omegac+CP%omegan)
    fnu = 1._dl - fcb

    
    do ik=1, MTrans%num_q_trans
       kh = MTrans%TransferData(Transfer_kh,ik,1)
       k = kh*h
       PK_data%log_kh(ik) = log(kh)

       do itf = 1, nz
          Tk_cb  = (CP%omegab)/(CP%omegab+CP%omegac) &
               *MTrans%TransferData(Transfer_b,ik,itf_start+itf-1) &
               +(CP%omegac)/(CP%omegab+CP%omegac) &
               *MTrans%TransferData(Transfer_cdm,ik,itf_start+itf-1)
          
          PK_data%matpower_cb(ik,itf) = &
               log(Tk_cb**2*k*pi*twopi*h**3*ScalarPower(k,in))
          
          if (fnu>0) then
             Tk_nu = MTrans%TransferData(Transfer_nu,ik,itf_start+itf-1)
             Tk_tot = fcb*Tk_cb + fnu*Tk_nu
             PK_data%matpower_rest(ik,itf) = &
                  log((Tk_tot**2-fcb**2*Tk_cb**2) & 
                  *k*pi*twopi*h**3*ScalarPower(k,in))
          else
             Tk_nu = 0.
             Tk_tot = fcb*Tk_cb + fnu*Tk_nu
             PK_data%matpower_rest(ik,itf) = -30.
          end if
          
          PK_data%matpower(ik,itf_start+itf-1) = &
               log(Tk_tot**2*k*pi*twopi*h**3*ScalarPower(k,in))
                    
       end do
    end do
    
    call MatterPowerdata_getsplines(PK_data)
    
  end subroutine Transfer_GetMatterPowerData_SPT
  
  subroutine Transfer_GetMatterPowerData(MTrans, PK_data, in, itf_only)
    !Does *NOT* include non-linear corrections
    !Get total matter P(k) in units of (h Mpc^{-1})^3 ready for interpolation.
    !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    !We are assuming that Cls are generated so any baryonic wiggles are well sampled 
    !and that matter power sepctrum is generated to beyond the CMB k_max
    Type(MatterTransferData), intent(in) :: MTrans
    Type(MatterPowerData) :: PK_data
    integer, intent(in) :: in
    integer, intent(in), optional :: itf_only
    real(dl) h, kh, k
    integer ik
    integer nz,itf, itf_start, itf_end
    
    if (present(itf_only)) then
       itf_start=itf_only
       itf_end = itf_only
       nz = 1
    else
       itf_start=1
       nz= size(MTrans%TransferData,3)
       itf_end = nz
    end if
    PK_data%num_k = MTrans%num_q_trans
    PK_Data%num_z = nz
    
    allocate(PK_data%matpower(PK_data%num_k,nz))
    allocate(PK_data%ddmat(PK_data%num_k,nz))
    allocate(PK_data%nonlin_ratio(PK_data%num_k,nz))
    allocate(PK_data%log_kh(PK_data%num_k))
    allocate(PK_data%redshifts(nz))
    PK_data%redshifts = CP%Transfer%Redshifts(itf_start:itf_end)
    
    h = CP%H0/100
    
    do ik=1,MTrans%num_q_trans
       kh = MTrans%TransferData(Transfer_kh,ik,1)
       k = kh*h
       PK_data%log_kh(ik) = log(kh)
       do itf = 1, nz
          PK_data%matpower(ik,itf) = &
               log(MTrans%TransferData(transfer_power_var,ik,itf_start+itf-1)**2*k & 
               *pi*twopi*h**3*ScalarPower(k,in))
       end do
    end do
    
    call MatterPowerdata_getsplines(PK_data)
    
  end subroutine Transfer_GetMatterPowerData
  !////////////////////////////////////////////////////////////////////////
  subroutine MatterPowerData_Load(PK_data,fname)
    !Loads in kh, P_k from file for one redshiftr and one initial power spectrum
    !Not redshift is not stored in file, so not set correctly
    !Also note that output _matterpower file is already interpolated, 
    !so re-interpolating is probs not a good idea
    
    !Get total matter power spectrum in units of (h Mpc^{-1})^3 
    !ready for interpolation.
    !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    use AmlUtils
    character(LEN=*) :: fname
    Type(MatterPowerData) :: PK_data
    real(dl)kh, Pk
    integer ik
    integer nz
    
    
    nz = 1
    call openTxtFile(fname, fileio_unit)
    
    PK_data%num_k = FileLines(fileio_unit)
    PK_Data%num_z = 1
    
    allocate(PK_data%matpower(PK_data%num_k,nz))
    allocate(PK_data%ddmat(PK_data%num_k,nz))
    allocate(PK_data%matpower_cb(PK_data%num_k,nz))
    allocate(PK_data%ddmat_cb(PK_data%num_k,nz))
    allocate(PK_data%nonlin_ratio(PK_data%num_k,nz))
    allocate(PK_data%log_kh(PK_data%num_k))
    
    allocate(PK_data%redshifts(nz))
    PK_data%redshifts = 0
    
    do ik=1,PK_data%num_k
       read (fileio_unit,*) kh, Pk
       PK_data%matpower(ik,1) = log(Pk) 
       PK_data%log_kh(ik) = log(kh)
    end do
    
    call MatterPowerdata_getsplines(PK_data)
    
  end subroutine MatterPowerData_Load
  !////////////////////////////////////////////////////////////////////////
  subroutine MatterPowerdata_getsplines(PK_data)
    Type(MatterPowerData) :: PK_data
    integer i
    real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
    
    do i = 1,PK_Data%num_z
       call spline(PK_data%log_kh,PK_data%matpower(1,i),PK_data%num_k,&
            cllo,clhi,PK_data%ddmat(1,i))
       call spline(PK_data%log_kh,PK_data%matpower_cb(1,i),PK_data%num_k,&
            cllo,clhi,PK_data%ddmat_cb(1,i))
    end do
    
  end subroutine MatterPowerdata_getsplines
  !////////////////////////////////////////////////////////////////////////
  subroutine MatterPowerdata_MakeNonlinear(PK_data)
    Type(MatterPowerData) :: PK_data
    
    call NonLinear_GetRatios(PK_data)
    PK_data%matpower = PK_data%matpower +  2*log(PK_data%nonlin_ratio)
    call MatterPowerdata_getsplines(PK_data)
    
  end subroutine MatterPowerdata_MakeNonlinear
  !////////////////////////////////////////////////////////////////////////
  subroutine MatterPowerdata_Free(PK_data)
    Type(MatterPowerData) :: PK_data
    integer i

    
    deallocate(PK_data%log_kh,stat=i)
    deallocate(PK_data%matpower,stat=i)
    deallocate(PK_data%ddmat,stat=i)
    deallocate(PK_data%matpower_cb,stat=i)
    deallocate(PK_data%ddmat_cb,stat=i)
    deallocate(PK_data%matpower_rest,stat=i)
    deallocate(PK_data%nonlin_ratio,stat=i)
    deallocate(PK_data%redshifts,stat=i)
    
    
  end subroutine MatterPowerdata_Free
  !////////////////////////////////////////////////////////////////////////
  !!===== ADDED BY S.S. =====
  function FuncPk(kh, Alogkh, Amatp, Addmatp, num_k) result(outpower)
    !Get power spectrum at particular k/h by interpolation
    !precomputed arrays, Alogkh, Amatp are prepared in logkh-logPk space
    real(dl), intent(in)  :: kh
    integer, intent(in)   :: num_k
    real(dl), intent(in)  :: Alogkh(num_k), Amatp(num_k), Addmatp(num_k)
    real(dl) :: logk
    integer llo,lhi
    real(dl) outpower, dp
    real(dl) ho, a0, b0
    integer, save :: i_last = 1  
    
    logk = log(kh)
    if (logk < Alogkh(1)) then
       dp = (Amatp(2)-Amatp(1))/(Alogkh(2)-Alogkh(1))
       outpower = Amatp(1) + dp*(logk-Alogkh(1))

    else if (logk > Alogkh(num_k)) then 
       dp = (Amatp(num_k)-Amatp(num_k-1))/(Alogkh(num_k)-Alogkh(num_k-1))
       outpower = Amatp(num_k) + dp*(logk-Alogkh(num_k))
    else 
       
       llo=min(i_last,num_k)
       do while (Alogkh(llo) > logk)
          llo=llo-1
       end do
       do while (Alogkh(llo+1) < logk)
          llo=llo+1
       end do
       i_last =llo  
       lhi=llo+1
       ho=Alogkh(lhi)-Alogkh(llo) 
       a0=(Alogkh(lhi)-logk)/ho
       b0=1-a0
       
       outpower = a0*Amatp(llo)+b0*Amatp(lhi) + &
            ((a0**3-a0)*Addmatp(llo) &
            +(b0**3-b0)*Addmatp(lhi))*ho**2/6
       
    end if
    
    outpower = exp(max(-30._dl,outpower))

  end function FuncPk
  !////////////////////////////////////////////////////////////////////////
  !!===== ADDED BY S.S. =====
  function FuncPk_cb(PK, kh, itf) result(outpower)
    !Get 'total' matter power spectrum at particular k/h by interpolation
    Type(MatterPowerData) :: PK
    integer, intent(in)   :: itf
    real(dl), intent(in)  :: kh
    real(dl) :: logk
    integer llo,lhi
    real(dl) outpower, dp
    real(dl) ho,a0,b0
    integer, save :: i_last = 1  
    
    logk = log(kh)
    if (logk < PK%log_kh(1)) then
       dp = (PK%matpower_cb(2,itf) -  PK%matpower_cb(1,itf)) / &
            ( PK%log_kh(2)-PK%log_kh(1) )
       outpower = PK%matpower_cb(1,itf) + dp*(logk - PK%log_kh(1))
    else if (logk > PK%log_kh(PK%num_k)) then
       !Do dodgy linear extrapolation on assumption accuracy of result won't matter
       
       dp = (PK%matpower_cb(PK%num_k,itf) -  PK%matpower_cb(PK%num_k-1,itf)) / &
            ( PK%log_kh(PK%num_k)-PK%log_kh(PK%num_k-1) )
       outpower = PK%matpower_cb(PK%num_k,itf) + dp*(logk - PK%log_kh(PK%num_k))
    else 
       
       llo=min(i_last,PK%num_k)
       do while (PK%log_kh(llo) > logk)
          llo=llo-1
       end do
       do while (PK%log_kh(llo+1)< logk)
          llo=llo+1
       end do
       i_last =llo  
       lhi=llo+1
       ho=PK%log_kh(lhi)-PK%log_kh(llo) 
       a0=(PK%log_kh(lhi)-logk)/ho
       b0=1-a0
       
       outpower = a0*PK%matpower_cb(llo,itf)+ b0*PK%matpower_cb(lhi,itf)+&
            ((a0**3-a0)*PK%ddmat_cb(llo,itf) &
            +(b0**3-b0)*PK%ddmat_cb(lhi,itf))*ho**2/6
       
    end if
    
    outpower = exp(max(-30._dl,outpower))
    
  end function FuncPk_cb
  !////////////////////////////////////////////////////////////////////////
  subroutine Transfer_GetMatterPower(MTrans,outpower, itf, in, minkh, dlnkh, npoints)
    !Allows for non-smooth priordial spectra
    !if CP%Nonlinear/ = NonLinear_none includes non-linear evolution
    !Get total matter power spectrum at logarithmically equal intervals dlnkh 
    !of k/h starting at minkh
    !in units of (h Mpc^{-1})^3.   
    !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    !We are assuming that Cls are generated so any baryonic wiggles are well 
    !sampled and that matter power
    !sepctrum is generated to beyond the CMB k_max
    Type(MatterTransferData), intent(in) :: MTrans
    Type(MatterPowerData) :: PK
    
    integer, intent(in) :: itf, in, npoints
    real, intent(out) :: outpower(npoints)
    real, intent(in) :: minkh, dlnkh
    real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
    integer ik, llo,il,lhi,lastix
    real(dl) :: matpower(MTrans%num_q_trans), kh
    real(dl) :: kvals(MTrans%num_q_trans), ddmat(MTrans%num_q_trans)
    real(dl) atransfer,xi, a0, b0, ho, logminkh,k, h, fcb, fnu
          

    if (npoints < 2) stop 'Need at least 2 points in Transfer_GetMatterPower'
    
    !         if (minkh < MTrans%TransferData(Transfer_kh,1,itf)) then
    !            stop 'Transfer_GetMatterPower: kh out of computed region'
    !          end if
    if (minkh*exp((npoints-1)*dlnkh) > & 
         MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf) &
         .and. FeedbackLevel > 0 ) &
         write(*,*) 'Warning: extrapolating matter power in Transfer_GetMatterPower'
    
    
    if (CP%NonLinear/=NonLinear_None) then
       call Transfer_GetMatterPowerData(MTrans, PK, in, itf)
       call NonLinear_GetRatios(PK)
      
    end if
    
    h = CP%H0/100
    logminkh = log(minkh)
    
    do ik=1,MTrans%num_q_trans
       kh = MTrans%TransferData(Transfer_kh,ik,itf)
       k = kh*h
       kvals(ik) = log(kh)
       atransfer=MTrans%TransferData(transfer_power_var,ik,itf)
       !In the case of Halofit
       !if (CP%NonLinear/=NonLinear_None) &
       if (CP%NonLinear==NonLinear_PkHF) & 
            atransfer = atransfer*PK%nonlin_ratio(ik,1) !only one element, this itf
   
            
       matpower(ik) = log(atransfer**2*k*pi*twopi*h**3)
       !Put in power spectrum later: transfer functions should be smooth, 
       !initial power may not be                
    end do
    
    call spline(kvals,matpower,MTrans%num_q_trans,cllo,clhi,ddmat)
    
    llo=1
    lastix = npoints + 1
    
    do il=1, npoints
       xi=logminkh + dlnkh*(il-1)
       if (xi < kvals(1)) then
          outpower(il)=-30.
          cycle
       end if
       do while ((xi > kvals(llo+1)).and.(llo < MTrans%num_q_trans))
          llo=llo+1
          if (llo >= MTrans%num_q_trans) exit
       end do
       if (llo == MTrans%num_q_trans) then
          lastix = il
          exit
       end if
       lhi=llo+1
       ho=kvals(lhi)-kvals(llo) 
       a0=(kvals(lhi)-xi)/ho
       b0=(xi-kvals(llo))/ho
       
       outpower(il) = a0*matpower(llo)+ b0*matpower(lhi)+((a0**3-a0)*ddmat(llo) &
            +(b0**3-b0)*ddmat(lhi))*ho**2/6
       
    end do
    
    do while (lastix <= npoints)
       !Do linear extrapolation in the log
       !Obviouly inaccurate, non-linear etc, but OK 
       !if only using in tails of window functions
       outpower(lastix) = 2*outpower(lastix-1) - outpower(lastix-2)
       lastix = lastix+1
    end do

    outpower = exp(max(-30.,outpower))
    
    do il = 1, npoints
       k = exp(logminkh + dlnkh*(il-1))*h
       outpower(il) = outpower(il)*ScalarPower(k,in) 
    end do
        
    if (CP%NonLinear/=NonLinear_None) call MatterPowerdata_Free(PK)
    
  end subroutine Transfer_GetMatterPower
  !////////////////////////////////////////////////////////////////////////
  subroutine Transfer_GetMatterPower_SPT(MT,outputPk, & 
       itf, in, minkh, dlnkh, npoints, num_id)
    Type(MatterTransferData), intent(in) :: MT
    integer, intent(in)  :: in, itf, npoints, num_id
    real, intent(in) :: minkh, dlnkh
    real, intent(out) :: outputPk(npoints,num_id)
   
    real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
    real(dl) :: kh, k, h, fcb, fnu, Tk_cb, Tk_nu, Tk_tot, Pk_rest_tmp
    real(dl) :: khvals(npoints), Pk_tot(npoints), Pk_cb(npoints), Pk_rest(npoints)
    real(dl), dimension(:), allocatable :: logkh, matp_cb, matp_rest, matp_tot
    real(dl), dimension(:), allocatable :: ddmatp_cb, ddmatp_rest, ddmatp_tot
    real(dl) :: Pk1loop(npoints,id_13vv) 
    real(dl) :: PkBias(npoints,numid_bias), A01(npoints,numid_A01)
    real(dl) :: Ak(npoints,id_A33), Bk(npoints,id_B422) 
    integer :: num_k, ik, id_tmp

    !if (CP%NonLinear.NE.NonLinear_PkSPT1loop) then 
    !   stop "inproper calling "
    !end if

    !First get Linear Pk
    h = CP%H0/100
    fcb = (CP%omegab+CP%omegac)/(CP%omegab+CP%omegac+CP%omegan)
    fnu = 1._dl - fcb
    num_k = MT%num_q_trans

    allocate(logkh(num_k))
    allocate(matp_cb(num_k))
    allocate(matp_rest(num_k))
    allocate(matp_tot(num_k))
    allocate(ddmatp_cb(num_k))
    allocate(ddmatp_rest(num_k))
    allocate(ddmatp_tot(num_k))
    
    do ik=1, num_k
       kh = MT%TransferData(Transfer_kh,ik,itf)
       k = kh*h
       logkh(ik) = dlog(kh)

       Tk_cb  = (CP%omegab)/(CP%omegab+CP%omegac) &
               *MT%TransferData(Transfer_b,ik,itf) &
               +(CP%omegac)/(CP%omegab+CP%omegac) &
               *MT%TransferData(Transfer_cdm,ik,itf)
       
       matp_cb(ik) = dlog(Tk_cb**2*k*pi*twopi*h**3*ScalarPower(k,in))
       
       if (fnu>0) then
          Tk_nu = MT%TransferData(Transfer_nu,ik,itf)
          Tk_tot = fcb*Tk_cb + fnu*Tk_nu

          Pk_rest_tmp = (2.*fcb*fnu*Tk_cb*Tk_nu+fnu**2.*Tk_nu**2.) & 
               *k*pi*twopi*h**3*ScalarPower(k,in)

          !Added by S.S. on Oct 21st
          !this avoids Pk_rest returns NAN in the case of tiny fnu
          if (Pk_rest_tmp.LE.1.e-11) then 
             matp_rest(ik) = -25.
          else
             matp_rest(ik) = dlog(Pk_rest_tmp)
             !dlog((Tk_tot**2.-fcb**2.*Tk_cb**2.)*k*pi*twopi*h**3*ScalarPower(k,in))
             !dlog((2.*fcb*fnu*Tk_cb*Tk_nu+fnu**2.*Tk_nu**2.)*k*pi*twopi*h**3*ScalarPower(k,in))
          end if
       else
          Tk_nu = 0.
          Tk_tot = fcb*Tk_cb + fnu*Tk_nu
          matp_rest(ik) = -30.
       end if
       
       matp_tot(ik) = log(Tk_tot**2*k*pi*twopi*h**3*ScalarPower(k,in))
      
    end do

    call spline(logkh,matp_cb,num_k,cllo,clhi,ddmatp_cb)
    call spline(logkh,matp_rest,num_k,cllo,clhi,ddmatp_rest)
    call spline(logkh,matp_tot,num_k,cllo,clhi,ddmatp_tot)

    !Evaluate power spectrum at each kh
    outputPk = 0.
    !kh & total Linear Pk
    do ik = 1, npoints
       khvals(ik)  = minkh*exp(dlnkh*(ik-1))
       outputPk(ik,id_kh) = real(khvals(ik))
       
       !Note Pk_tot = fcb**2*Pk_cb + Pk_rest
       Pk_tot(ik)  = FuncPk(khvals(ik),logkh,matp_tot,ddmatp_tot,num_k)
       Pk_cb(ik)   = FuncPk(khvals(ik),logkh,matp_cb ,ddmatp_cb,num_k)
       Pk_rest(ik) = FuncPk(khvals(ik),logkh,matp_rest,ddmatp_rest,num_k)
       outputPk(ik,id_PkLin) = real(Pk_tot(ik))
    end do
    
    if (CP%NonLinear==NonLinear_PkSPT1loop) then
       !Compute Nonlinear correction

       call GetOneLoopCorrections(khvals,npoints,Pk1loop, & 
            logkh,matp_cb,ddmatp_cb,matp_tot,ddmatp_tot,num_k)
       
       !modified by SS on Oct/7/2014
       outputPk(:,id_PkNL_dd) = real( & 
            fcb**2.*(Pk_cb(:)+Pk1loop(:,id_22dd)+Pk1loop(:,id_13dd)) & 
            + Pk_rest(:))
       !outputPk(:,id_PkNL_dd) = real(fcb**2.*(Pk1loop(:,id_22dd)))
       
       outputPk(:,id_PkNL_dv) = real( & 
            fcb**2.*(Pk_cb(:)+Pk1loop(:,id_22dv)+Pk1loop(:,id_13dv)) & 
            + Pk_rest(:))
       
       outputPk(:,id_PkNL_vv) = real( & 
            fcb**2.*(Pk_cb(:)+Pk1loop(:,id_22vv)+Pk1loop(:,id_13vv)) & 
            + Pk_rest(:))
       
       outputPk(:,id_Pkb2_dd) = real(Pk1loop(:,id_b2dd))
       outputPk(:,id_Pkb2_dv) = real(Pk1loop(:,id_b2dv))
       outputPk(:,id_Pkb22)   = real(Pk1loop(:,id_b22))

       id_tmp = id_Pkb22

       if (CP%DoNonLocalBias) then
          
          call GetBiasCorrections(khvals,npoints,PkBias, & 
               logkh,matp_tot,ddmatp_tot,num_k)
          
          
          outputPk(:,id_tmp+1) = real(PkBias(:,id_bs2d))
          outputPk(:,id_tmp+2) = real(PkBias(:,id_bs2v))
          outputPk(:,id_tmp+3) = real(PkBias(:,id_b2s2))
          outputPk(:,id_tmp+4) = real(PkBias(:,id_bs22))
          outputPk(:,id_tmp+5) = real(PkBias(:,id_b3))
          
          call GetA01Corrections(khvals,npoints,A01, & 
               logkh,matp_tot,ddmatp_tot,num_k)
          
          outputPk(:,id_tmp+6:) = real(A01(:,:))
          id_tmp = id_tmp + numid_bias + numid_A01
       end if
       
       if (CP%NonLinear_RSDmodel==NLKaiser_AB_FoG) then
          call GetRSDCorrections(khvals,npoints,Ak,Bk, &
               logkh,matp_cb,ddmatp_cb,num_k)
       
          !assuming neutrino component is negligibly small
          outputPk(:,id_tmp+1:id_tmp+id_A33) = fcb**3.*real(Ak(:,:))
          outputPk(:,id_tmp+id_A33+1:id_tmp+id_A33+id_B422) & 
               = fcb**4.*real(Bk(:,:))
       end if
    else
       outputPk(:,id_PkLin+1:) = 0.
    end if
    
    deallocate(logkh,matp_cb,matp_rest,matp_tot)
    deallocate(ddmatp_cb,ddmatp_rest,ddmatp_tot)
    
  end subroutine Transfer_GetMatterPower_SPT
  !////////////////////////////////////////////////////////////////////////
  subroutine Transfer_Get_sigma8(MTrans, sigr8)
    use MassiveNu
    Type(MatterTransferData) :: MTrans
    integer ik, itf, in
    real(dl) kh, k, h, x, win, delta
    real(dl) lnk, dlnk, lnko
    real(dl) dsig8, dsig8o, sig8, sig8o, powers
    real(dl), intent(IN) :: sigr8
    
    !Calculate MTrans%sigma_8^2 = int dk/k win**2 T_k**2 P(k), 
    !where win is the FT of a spherical top hat
    !of radius sigr8 h^{-1} Mpc
    
    H=CP%h0/100._dl
    do in = 1, CP%InitPower%nn
       do itf=1,CP%Transfer%num_redshifts
          lnko=0
          dsig8o=0
          sig8=0
          sig8o=0
          do ik=1, MTrans%num_q_trans
             kh = MTrans%TransferData(Transfer_kh,ik,itf)
             if (kh==0) cycle
             k = kh*H
             
             delta = k**2*MTrans%TransferData(transfer_power_var,ik,itf)
             !if (CP%NonLinear/=NonLinear_None) delta= delta* 
             !MTrans%NonLinearScaling(ik,itf)
             !sigma_8 defined "as though it were linear"
             
             x= kh *sigr8
             win =3*(sin(x)-x*cos(x))/x**3
             lnk=log(k)
             if (ik==1) then
                dlnk=0.5_dl 
                !Approx for 2._dl/(CP%InitPower%an(in)+3)  
                ![From int_0^k_1 dk/k k^4 P(k)]
                !Contribution should be very small in any case 
             else
                dlnk=lnk-lnko
             end if
             powers = ScalarPower(k,in)
             dsig8=(win*delta)**2*powers
             sig8=sig8+(dsig8+dsig8o)*dlnk/2
             dsig8o=dsig8
             lnko=lnk
             
          end do
          
          MTrans%sigma_8(itf,in) = sqrt(sig8)
       end do
    end do
    
  end subroutine Transfer_Get_sigma8
  !////////////////////////////////////////////////////////////////////////
  subroutine Transfer_output_Sig8(MTrans)
    Type(MatterTransferData), intent(in) :: MTrans
    real(dl) :: zin, dlnDdlna
    integer in, j
    
    do in=1, CP%InitPower%nn
       if (CP%InitPower%nn>1)  write(*,*) 'Power spectrum : ', in
       do j = 1, CP%Transfer%num_redshifts
          
          write(*,*) 'at z = ',real(CP%Transfer%redshifts(j)), & 
               ' sigma8 (all matter)=', real(MTrans%sigma_8(j,in))

          zin = CP%Transfer%redshifts(j)
          call ComputeGrowth(zin,dlnDdlna)
       end do
    end do
    
  end subroutine Transfer_output_Sig8
  !////////////////////////////////////////////////////////////////////////
  subroutine Transfer_output_Sig8AndNorm(MTrans)
    Type(MatterTransferData), intent(in) :: MTrans
    integer in, j
    
    do in=1, CP%InitPower%nn
       write(*,*) 'Power spectrum ',in, ' COBE_scale = ',real(COBE_scales(in))
       do j = 1, CP%Transfer%num_redshifts
          write(*,*) 'at z = ',real(CP%Transfer%redshifts(j)), &
               ' sigma8(all matter) = ', &
               real(MTrans%sigma_8(j,in)*sqrt(COBE_scales(in)))
       end do
    end do
    
  end subroutine Transfer_output_Sig8AndNorm
  !////////////////////////////////////////////////////////////////////////
  subroutine Transfer_Allocate(MTrans)
    Type(MatterTransferData) :: MTrans
    integer st
    
    deallocate(MTrans%q_trans, STAT = st)
    deallocate(MTrans%TransferData, STAT = st)
    deallocate(MTrans%sigma_8, STAT = st)
    allocate(MTrans%q_trans(MTrans%num_q_trans))            
    allocate(MTrans%TransferData(Transfer_max,MTrans%num_q_trans, & 
         CP%Transfer%num_redshifts))  
    allocate(MTrans%sigma_8(CP%Transfer%num_redshifts, CP%InitPower%nn))
    
  end  subroutine Transfer_Allocate
  !////////////////////////////////////////////////////////////////////////
  subroutine Transfer_Free(MTrans)
    Type(MatterTransferData):: MTrans
    integer st
    
    deallocate(MTrans%q_trans, STAT = st)
    deallocate(MTrans%TransferData, STAT = st)
    deallocate(MTrans%sigma_8, STAT = st)
    nullify(MTrans%q_trans)
    nullify(MTrans%TransferData)
    nullify(MTrans%sigma_8)
    
  end subroutine Transfer_Free
  !////////////////////////////////////////////////////////////////////////
  subroutine Transfer_SetForNonlinearLensing(P)
    Type(TransferParams) :: P
    integer i
    
    P%kmax = 5*AccuracyBoost
    P%k_per_logint  = 0
    P%num_redshifts =  nint(10*AccuracyBoost)
    if (P%num_redshifts > max_transfer_redshifts) &
         stop 'Transfer_SetForNonlinearLensing: Too many redshifts'
    do i=1,P%num_redshifts
       P%redshifts(i) = real(P%num_redshifts-i)/(P%num_redshifts/10)
    end do
    
  end subroutine Transfer_SetForNonlinearLensing
  !////////////////////////////////////////////////////////////////////////
  subroutine Transfer_SaveToFiles(MTrans,FileNames)
    use IniFile
    Type(MatterTransferData), intent(in) :: MTrans
    integer i,ik
    character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
    
    do i=1, CP%Transfer%num_redshifts
       if (FileNames(i) /= '') then
          open(unit=fileio_unit,file=FileNames(i),form='formatted',status='replace')
          do ik=1,MTrans%num_q_trans
             if (MTrans%TransferData(Transfer_kh,ik,i)/=0) then
                write(fileio_unit,'(7E14.6)')  & 
                     MTrans%TransferData(Transfer_kh:Transfer_max,ik,i)
             end if
          end do
          close(fileio_unit)
       end if
    end do
    
    
  end subroutine Transfer_SaveToFiles
  !////////////////////////////////////////////////////////////////////////
  subroutine Transfer_SaveMatterPower(MTrans, FileNames, Pk_NLB_FileNames, Bk_NLB_FileNames)
    use IniFile
    !Export files of total  matter power spectra in h^{-1} Mpc units, 
    !against k/h.
    Type(MatterTransferData), intent(in) :: MTrans
    character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
    character(LEN=Ini_max_string_len), intent(IN), OPTIONAL :: Pk_NLB_FileNames(*), Bk_NLB_FileNames(*)
    integer itf,in,i, num_id
    integer points
    real, dimension(:,:), allocatable :: outpower
    real, dimension(:,:,:), allocatable :: outputPk
    real(dl), dimension(:,:), allocatable :: Ps0
    character(LEN=80) fmt, fmt3
    real minkh,dlnkh
    real(dl) :: dlnDdlna, zin, f, sigma_v
    Type(MatterPowerData) :: PK_data
    
    write (fmt,*) CP%InitPower%nn+1
    fmt  = '('//trim(adjustl(fmt))//'E15.5)'

    do itf=1, CP%Transfer%num_redshifts
       if (FileNames(itf) /= '') then
          if (.not. transfer_interp_matterpower ) then
             
             points = MTrans%num_q_trans
             allocate(outpower(points,CP%InitPower%nn))
             
             do in = 1, CP%InitPower%nn
                
                call Transfer_GetMatterPowerData(MTrans, PK_data, in, itf)
                
                if (CP%NonLinear/=NonLinear_None) & 
                     call MatterPowerdata_MakeNonlinear(PK_Data)
                
                outpower(:,in) = exp(PK_data%matpower(:,1))
                call MatterPowerdata_Free(PK_Data)
             end do
             
             open(unit=fileio_unit,file=FileNames(itf),form='formatted', & 
                  status='replace')
             do i=1,points
                write (fileio_unit, fmt) MTrans%TransferData(Transfer_kh,i,1), & 
                     outpower(i,1:CP%InitPower%nn)
             end do
             close(fileio_unit)
             
          else
             !minkh = 1e-4
             !dlnkh = 0.02
             !points = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf) & 
             !     /minkh)/dlnkh+1

             minkh = real(MTrans%TransferData(Transfer_kh,1,itf))
             points = MTrans%num_q_trans
             dlnkh = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf) & 
             /minkh)/(points-0.999)

             if (CP%NonLinear==NonLinear_PkSPT1loop) then     
                if (CP%DoNonLocalBias) then 
                   num_id = id_Pkb22 + numid_bias + numid_A01
                else
                   num_id = id_Pkb22
                end if
             else
                num_id = id_PkLin
             end if

             if (CP%NonLinear_RSDmodel==NLKaiser_AB_FoG) then
                num_id = num_id + numid_RSDTNS
             end if

             write (fmt3,*) num_id
             
             fmt3 = '('//trim(adjustl(fmt3))//'E15.5)'

             allocate(outputPk(points,num_id,CP%InitPower%nn))

             allocate(outpower(points,CP%InitPower%nn))

             !zin = 0.55
             !call ComputeGrowth(zin,dlnDdlna)
             !stop

             do in = 1, CP%InitPower%nn
                !call Transfer_GetMatterPower(MTrans,outpower(1,in),itf,in, & 
                !     minkh,dlnkh, points)

                call Transfer_GetMatterPower_SPT(MTrans, &
                     outputPk(:,:,in),itf,in,minkh,dlnkh,points,num_id)

                !allocate(Ps0(points,num_id))
                !zin = 0.
                !call ComputeGrowth(zin,f)
                !sigma_v = 3.0
                !call CalcMonopole(dble(outputPk(:,:,in)),points,num_id,f,sigma_v,Ps0)
                !stop

                if (CP%OutputNormalization == outCOBE) then
                   if (allocated(COBE_scales)) then
                      outpower(:,in) = outpower(:,in)*COBE_scales(in)
                   else
                      if (FeedbackLevel>0) then 
                         write (*,*) 'Cannot COBE normalize - no Cls generated'
                      end if
                   end if
                end if
             end do
             
             open(unit=fileio_unit,file=FileNames(itf),form='formatted', & 
                  status='replace')
             if (CP%NonLinear==NonLinear_PkSPT1loop) then
                !zin = PK_data%redshifts(itf)
                !call ComputeGrowth(zin,dlnDdlna)
                write(*,*) "  # of output Pk: ", num_id
                write(*,*) "  ouput: k P_lin P_MC ..."
                do i=1,points
                   write(fileio_unit, fmt3) outputPk(i,:,1:CP%InitPower%nn)
                end do
             else if (CP%NonLinear==NonLinear_None) then
                do i=1,points
                   !write (fileio_unit, fmt) minkh*exp((i-1)*dlnkh), & 
                   !     outpower(i,1:CP%InitPower%nn)
                   write(fileio_unit, fmt) outputPk(i,1:2,1:CP%InitPower%nn)
                end do
             else
                !zin = PK_data%redshifts(itf)
                !call ComputeGrowth(zin,dlnDdlna)
                write(*,*) "  # of outpout Pk: ", num_id
                do i = 1, points
                   write(fileio_unit, fmt3) outputPk(i,:,1:CP%InitPower%nn)
                end do
             end if
             close(fileio_unit)
             
          end if

          !if (CP%NonLinear==NonLinear_PkSPT1loop) then
          deallocate(outputPk)
          !end if
          
          deallocate(outpower) 
       end if
    end do
    
  end subroutine Transfer_SaveMatterPower
  !////////////////////////////////////////////////////////////////////////
  !===== ADDED BY S.S. =====
  subroutine GetOneLoopCorrections(khvals,npoints,Pk1loop, &
       alogkh,amatp_cb,addmatp_cb,amatp_tot,addmatp_tot,num_k)
    !******************************************************************
    !  compute the nonlinear correction for P_cb part 
    !  use total P(k) for nonlinear bias terms 
    !  based on the Standard Perturbation Theory at 1-loop level
    !
    !
    !  khvals(npoints): precomputed array for kh to be evaluated
    !  Pk1loop(npoints): output 1loop nonlinear corrections
    !    
    !  alogkh,amatp_cb,addmatp_cb,amatp_tot,addmatp_tot:
    !     precompured array with dimension of (num_k) 
    !     for logkh, logPk_cb, logPk_tot
    !
    !******************************************************************
   
    integer, intent(in)   :: npoints, num_k
    real(dl), intent(in)  :: khvals(npoints), alogkh(num_k)
    real(dl), intent(in)  :: amatp_cb(num_k), addmatp_cb(num_k) 
    real(dl), intent(in)  :: amatp_tot(num_k), addmatp_tot(num_k)

    real(dl), intent(OUT) :: Pk1loop(npoints,id_13vv)

    integer, parameter :: numx=500, nummu=30
    real(dl) :: tmin = -1._dl, tmax = 1._dl
    real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
    real(dl) :: aXorg(numx), aWXorg(numx), aMUorg(nummu), aWMUorg(nummu)
    real(dl) :: aX(numx), aWX(numx), aMU(nummu), aWMU(nummu)
    real(dl) :: muintegrand(id_b22), Pkt, P2, P3
    integer  :: ik, ik_start, ix, imu, id
    real(dl) :: kh_val, kh_tmp, logxmin, logxmax, mumin, mumax


    call GauLeg(tmin,tmax,numx,aXorg,aWXorg)
    call GauLeg(tmin,tmax,nummu,aMUorg,aWMUorg)

    Pk1loop = 0

    !need not to compute kh<0.005 where nonlinear corrections 
    !are negligible
    ik = 1
    ik_start = 1
    do while (khvals(ik)<0.005)
       ik = ik +1
       ik_start = ik
    end do

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC) &
    !$OMP & PRIVATE(kh_val,logxmin,logxmax,aX,aWX,P2,P3, & 
    !$OMP &         mumin,mumax,aMU,aWMU,muintegrand,kh_tmp)
    !do ik = ik_start, npoints
    do ik = 1, npoints
       kh_val = khvals(ik)
       
       logxmin = dlog(dexp(alogkh(1))/kh_val)
       logxmax = dlog(20._dl/kh_val)

       !rescale Gaussi-Legendre roots and weights
       aX  = ((logxmax-logxmin)*aXorg+(logxmax+logxmin))/2._dl
       aWX = (logxmax-logxmin)/2._dl*aWXorg

       !x-integration
       do ix = 1, numx
          mumin = (1._dl + dexp(aX(ix))**2 - dexp(logxmax)**2)/(2._dl*dexp(aX(ix)))
          mumax = (1._dl + dexp(aX(ix))**2 - dexp(logxmin)**2)/(2._dl*dexp(aX(ix)))

          !integral region
          if (mumin < -1._dl) mumin = -1._dl
          
          if (dexp(aX(ix)) <= 0.5_dl .and. mumax >= 1._dl) then
             mumax = 1._dl
          else if (dexp(aX(ix))>0.5_dl .and. & 
                   mumax>1._dl/(2._dl*dexp(aX(ix)))) then
             mumax = 1._dl/(2._dl*dexp(aX(ix)))
          end if
          !rescale Gauss-Legendre roots and weights
          aMU  = ((mumax-mumin)*aMUorg+(mumax+mumin))/2._dl
          aWMU = (mumax-mumin)/2._dl*aWMUorg

          P3 = FuncPk(kh_val*dexp(aX(ix)),alogkh,amatp_tot,addmatp_tot,num_k)

          !mu-integration
          muintegrand = 0
          do imu = 1, nummu
             kh_tmp = kh_val & 
                  *dsqrt(1._dl+dexp(aX(ix))**2.-2._dl*dexp(aX(ix))*aMU(imu))

             P2 = FuncPk(kh_tmp,alogkh,amatp_cb,addmatp_cb,num_k)

             do id = id_22dd, id_b22
                if (id<=id_22vv) then
                   muintegrand(id) = muintegrand(id) + &
                        aWMU(imu)*Kernel22(dexp(aX(ix)),aMU(imu),id)*P2
                        !*FuncPk(kh_tmp,alogkh,amatp_cb,addmatp_cb,num_k)
               
                else if (id<=id_b2dv) then
                   muintegrand(id) = muintegrand(id) + &
                        aWMU(imu)*Kernel22(dexp(aX(ix)),aMU(imu),id)*P2
                        !*FuncPk(kh_tmp,alogkh,amatp_tot,addmatp_tot,num_k)
                else
                   !muintegrand(id) = muintegrand(id) + aWMU(imu) & 
                   !     *( FuncPk(kh_tmp,alogkh,amatp_tot,addmatp_tot,num_k) & 
                   !       -FuncPk(kh_val*dexp(aX(ix)),alogkh,amatp_tot, & 
                   !       addmatp_tot,num_k))/8._dl
                   muintegrand(id) = muintegrand(id) + aWMU(imu) &
                        *(2.*P2*P3-P2**2.-P3**2.)/8.
                end if
             end do
          end do
          
          do id = id_22dd, id_13vv
             if (id<id_b22) then
                Pk1loop(ik,id) = Pk1loop(ik,id) + &
                     aWX(ix)*muintegrand(id)*dexp(aX(ix))**3.*P3
                     !*FuncPk(kh_val*dexp(aX(ix)),alogkh,amatp_cb,addmatp_cb,num_k)
             else if (id==id_b22) then
                Pk1loop(ik,id) = Pk1loop(ik,id) + &
                     aWX(ix)*muintegrand(id)*dexp(aX(ix))**3.
             else
                Pk1loop(ik,id) = Pk1loop(ik,id) + &
                     aWX(ix)*Kernel13(dexp(aX(ix)),id)*dexp(aX(ix))*P3
                     !*FuncPk(kh_val*dexp(aX(ix)),alogkh,amatp_cb,addmatp_cb,num_k)
             end if
          end do
       end do

       !Be careful: this notation is based on Saito et al. (2010).
       !Pk1loop(ik,:id_b22) = Pk1loop(ik,:id_b22)*kh_val**3/const_pi2
       !Pk1loop(ik,id_13dd:id_13vv) = Pk1loop(ik,id_13dd:id_13vv) &
       !     *FuncPk(kh_val,alogkh,amatp_cb,addmatp_cb,num_k) & 
       !     *kh_val**3/(2._dl*const_pi2)

       !corrected for nonlinear bias
       !note not corrected for P_b22
       Pk1loop(ik,:id_22vv) = Pk1loop(ik,:id_22vv)*kh_val**3/const_pi2
       Pk1loop(ik,id_b2dd:id_b2dv) = Pk1loop(ik,id_b2dd:id_b2dv)*kh_val**3/const_pi2/2._dl
       Pk1loop(ik,id_b22)   = Pk1loop(ik,id_b22)*kh_val**3/const_pi2
       Pk1loop(ik,id_13dd:id_13vv) = Pk1loop(ik,id_13dd:id_13vv) &
            *FuncPk(kh_val,alogkh,amatp_cb,addmatp_cb,num_k) & 
            *kh_val**3/(2._dl*const_pi2)
 
    end do
    !$OMP END PARALLEL DO
  
  end subroutine GetOneLoopCorrections
  !////////////////////////////////////////////////////////////////////////
  ! Standard Perturbation Theory Kernels
  !////////////////////////////////////////////////////////////////////////
  function Kernel22(r,mu,id)
    real(dl), intent(IN) :: r, mu
    integer, intent(IN)  :: id
    real(dl) :: Kernel22

    
    integer, parameter :: id_13dd=7, id_13dv=8, id_13vv=9
    
    if (id==id_22dd) then
       Kernel22 = (3._dl*r+7._dl*mu-10._dl*r*mu**2)**2
       Kernel22 = Kernel22/(1._dl+r*r-2._dl*r*mu)**2/196._dl/r**2

    else if (id==id_22dv) then
       Kernel22 = (-r+7._dl*mu-6._dl*r*mu**2)*(3._dl*r+7._dl*mu-10._dl*r*mu**2)
       Kernel22 = Kernel22/(1._dl+r*r-2._dl*r*mu)**2/196._dl/r**2

    else if (id==id_22vv) then
       Kernel22 = (-r+7._dl*mu-6._dl*r*mu**2)**2
       Kernel22 = Kernel22/(1._dl+r*r-2._dl*r*mu)**2/196._dl/r**2

    else if (id==id_b2dd) then
       Kernel22 = 3._dl*r+7._dl*mu-10._dl*r*mu**2
       Kernel22 = Kernel22/(1._dl+r*r-2._dl*r*mu)/14._dl/r

    else if (id==id_b2dv) then
       Kernel22 = -r+7._dl*mu-6._dl*r*mu**2
       Kernel22 = Kernel22/(1._dl+r*r-2._dl*r*mu)/14._dl/r

    else if (id==id_b22 ) then
       Kernel22 = 1._dl

    else
       write(*,*) "stop: an error occurs in function Kernel22"
       write(*,*) "index is out of supported range"
       stop

    end if
      
  end function Kernel22
  !////////////////////////////////////////////////////////////////////////
  function Kernel13(r,id)
    real(dl), intent(IN) :: r
    integer,  intent(IN) :: id
    real(dl) :: r2, Kernel13
    
    r2 = r * r
    if (id==id_13dd) then
       if (r<200._dl) then
          Kernel13 = 12._dl/r2-158._dl+100._dl*r2-42._dl*r2**2 & 
               +(3._dl/(r**3))*((r2-1._dl)**3) &
               *(7._dl*r2+2._dl)*log((1._dl+r)/abs(1._dl-r))
          Kernel13 = Kernel13/504._dl
       else
          Kernel13 = -97.6_dl+19.2_dl/r2
          Kernel13 = Kernel13/504._dl
       end if

    else if (id==id_13dv) then
       if (r<200._dl) then
          Kernel13 = 24._dl/r2-202._dl+56._dl*r2-30._dl*r2**2 & 
               +(3._dl/(r**3))*((r2-1._dl)**3) &
               *(5._dl*r2+4._dl)*log((1._dl+r)/abs(1._dl-r))
          Kernel13 = Kernel13/504._dl
       else
          Kernel13 = -200._dl+2208._dl/35._dl/r2
          Kernel13 = Kernel13/504._dl
       end if
    else if (id==id_13vv) then
       if (r<200._dl) then
          Kernel13 = 12._dl/r2-82._dl+4._dl*r2-6._dl*r2**2 & 
               +(3._dl/(r**3))*((r2-1._dl)**3) &
               *(r2+2._dl)*log((1._dl+r)/abs(1._dl-r))
          Kernel13 = Kernel13/168._dl
       else
          Kernel13 = -504._dl/5._dl +1248._dl/35._dl/r2
          Kernel13 = Kernel13/168._dl
       end if
    
    else
       write(*,*) "stop: an error occurs in function Kernel13"
       write(*,*) "index is out of supported range"
       stop
       
    end if
  end function Kernel13
  !////////////////////////////////////////////////////////////////////////
  !===== ADDED BY S.S. =====
  subroutine GetRSDCorrections(khvals,npoints,Ak,Bk, &
       alogkh,amatp_cb,addmatp_cb,num_k)
    !******************************************************************
    !  compute A & B correction terms in redshift-space power spectrum
    !  
    !  see Taruya, Nishimichi, Saito (2010)
    !
    !  Note that neutrino contributions to the corrections are assumed 
    !  to be neglected. 
    !
    !  khvals(npoints): precomputed array for kh to be evaluated
    !  Ak(npoints,): A correction terms
    !  Bk(npoints): B correction terms 
    !    
    !  alogkh,amatp_cb,addmatp_cb:
    !     precompured array with dimension of (num_k) 
    !     for logkh, logPk_cb
    !
    !******************************************************************
   
    integer, intent(in)   :: npoints, num_k
    real(dl), intent(in)  :: khvals(npoints), alogkh(num_k)
    real(dl), intent(in)  :: amatp_cb(num_k), addmatp_cb(num_k) 

    real(dl), intent(OUT) :: Ak(npoints,id_A33), Bk(npoints,id_B422)

    integer, parameter :: numx=500, nummu=30
    real(dl) :: tmin = -1._dl, tmax = 1._dl
    real(dl) :: aXorg(numx), aWXorg(numx), aMUorg(nummu), aWMUorg(nummu)
    real(dl) :: aX(numx), aWX(numx), aMU(nummu), aWMU(nummu)
   
    real(dl) :: P1, P2, P3
    integer  :: ik, ik_start, ix, imu, id
    real(dl) :: kh_val, xmin, xmax, mumin, mumax, x, mu
    real(dl), dimension(id_A33) :: A23, A31, A12
    real(dl), dimension(id_B422) :: B23

    call GauLeg(tmin,tmax,numx,aXorg,aWXorg)
    call GauLeg(tmin,tmax,nummu,aMUorg,aWMUorg)

    !need not to compute kh<0.005 where nonlinear corrections 
    !are negligible
    ik = 1
    ik_start = 1
    do while (khvals(ik)<0.001)
       ik = ik +1
       ik_start = ik
    end do

    Ak = 0.
    Bk = 0.
    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC) &
    !$OMP & PRIVATE(kh_val,x,xmin,xmax,aX,aWX, & 
    !$OMP &         mu,mumin,mumax,aMU,aWMU,P1,P2,P3,A23,A31,A12,B23)
    do ik = ik_start, npoints
       kh_val = khvals(ik)
       
       xmin = dlog(dexp(alogkh(1))/kh_val)
       xmax = dlog(20._dl/kh_val)

       !rescale Gaussi-Legendre roots and weights
       aX  = ((xmax-xmin)*aXorg+(xmax+xmin))/2._dl
       aWX = (xmax-xmin)/2._dl*aWXorg

       P1 = FuncPk(kh_val,alogkh,amatp_cb,addmatp_cb,num_k)

       !x-integration
       do ix = 1, numx
          x = dexp(aX(ix))

          P3 = FuncPk(kh_val*x,alogkh,amatp_cb,addmatp_cb,num_k)

          mumin = max(-1., (1.+x**2-dexp(xmax)**2)/2./x)
          mumax = min( 1., (1.+x**2-dexp(xmin)**2)/2./x)
          if (x>0.5_dl) mumax = 0.5_dl/x

          !rescale Gauss-Legendre roots and weights
          aMU  = ((mumax-mumin)*aMUorg+(mumax+mumin))/2._dl
          aWMU = (mumax-mumin)/2._dl*aWMUorg
          
          A23 = 0.
          A31 = 0.
          A12 = 0.
          B23 = 0.
          !mu-integration
          do imu = 1, nummu
             mu = aMU(imu)
             
             P2 = FuncPk(kh_val*dsqrt(1.+x**2.-2.*x*mu), & 
                  alogkh,amatp_cb,addmatp_cb,num_k)

             !A terms
             do id = id_A11, id_A33
                A23(id) = A23(id) + aWMU(imu)*KernelA23(id,mu,x)*P2*P3
                A31(id) = A31(id) + aWMU(imu)*KernelA31(id,mu,x)*P3*P1
                A12(id) = A12(id) + aWMU(imu)*KernelA12(id,mu,x)*P1*P2
             end do

             !B terms
             do id = id_B111, id_B422
                B23(id) = B23(id) + aWMU(imu)*KernelB23(id,mu,x)*P2*P3
             end do
          end do

          Ak(ik,:) = Ak(ik,:) + aWX(ix)*x*(A23(:)+A31(:)+A12(:))
          Bk(ik,:) = Bk(ik,:) + aWX(ix)*x*B23(:)
       end do

       !factor 2 comes from the half-plane 
       Ak(ik,:) = 2.* Ak(ik,:) * kh_val**3/(2.*pi)**2.
       Bk(ik,:) = 2.* Bk(ik,:) * kh_val**3/(2.*pi)**2.
      
    end do
    !$OMP END PARALLEL DO
  
  end subroutine GetRSDCorrections
  !////////////////////////////////////////////////////////////////////////
  !This kernels correspond to tilde{A}_mn in TNS(2010)
  function KernelA23(id,mu,x)
    real(dl), intent(IN) :: mu, x
    integer,  intent(IN) :: id
    real(dl) :: KernelA23
    
    if (id==id_A11) then
       KernelA23 = (-mu+x*(2.*mu*mu-1.))*(-3.*x-7.*mu+10.*x*mu*mu) & 
            /7./(1.+x*x-2.*mu*x)**2.
    else if (id==id_A12) then
       KernelA23 = x*(mu*mu-1.)*(3.*x+7.*mu-10.*x*mu*mu) & 
            /14./(1.+x*x-2.*mu*x)**2.
    else if (id==id_A22) then
       KernelA23 = (28.*mu*mu+x*mu*(25.-81.*mu*mu)+ & 
            x**2*(1.-27.*mu*mu+54.*mu**4)) & 
            /14./(1.+x*x-2.*mu*x)**2.
    else if (id==id_A23) then
       KernelA23 = x*(1.-mu*mu)*(x-7.*mu+6.*x*mu*mu) & 
            /14./(1.+x*x-2.*mu*x)**2.
    else if (id==id_A33) then
       KernelA23 = (x-7.*mu+6.*x*mu*mu)*(-2.*mu-x+3.*x*mu*mu) & 
            /14./(1.+x*x-2.*mu*x)**2.
    end if
    
  end function KernelA23
  !////////////////////////////////////////////////////////////////////////
  !This kernels correspond to a_mn in TNS(2010)
  function KernelA31(id,mu,x)
    real(dl), intent(IN) :: mu, x
    integer,  intent(IN) :: id
    real(dl) :: KernelA31

    if (id==id_A11) then
       KernelA31 =  (-7.*mu*mu+x**3*mu*(-3.+10.*mu*mu) & 
            +3.*x*(mu+6.*mu**3)+x*x*(6.-19.*mu*mu-8.*mu**4)) & 
            /7./(1.+x*x-2.*mu*x)
    else if (id==id_A12.or.id==id_A23) then
       KernelA31 =  x*(-1.+mu*mu)*(6.*x-7.*(1.+x*x)*mu+8.*x*mu*mu) & 
            /14./(1.+x*x-2.*mu*x)
    else if (id==id_A22) then
       KernelA31 = (-28.*mu*mu+x**3*mu*(-13.+41.*mu*mu)+ & 
            x*mu*(11.+73.*mu*mu)-2.*x*x*(-9.+31.*mu*mu+20.*mu**4)) & 
            /14./(1.+x*x-2.*mu*x)
    else if (id==id_A33) then
       KernelA31 = (7.*mu+x*(-6.+7.*x*mu-8.*mu*mu))*(-2.*mu+x*(-1.+3.*mu*mu)) & 
            /14./(1.+x*x-2.*mu*x)
    end if
  
  end function KernelA31
  !////////////////////////////////////////////////////////////////////////
  !This kernels correspond to A_mn in TNS(2010)
  function KernelA12(id,mu,x)
    real(dl), intent(IN) :: mu, x
    integer,  intent(IN) :: id
    real(dl) :: KernelA12

    if (id==id_A11) then
       KernelA12 = -x**3.*(mu+6.*mu**3+x*x*mu*(-3.+10.*mu*mu)+ & 
            x*(-3.+mu*mu-12.*mu**4))/7./(1.+x*x-2.*mu*x)**2.
    else if (id==id_A12.or.id==id_A23) then
       KernelA12 = x**4.*(mu*mu-1.)*(-1.+7.*x*mu-6.*mu*mu) & 
            /14./(1.+x*x-2.*mu*x)**2.
    else if (id==id_A22) then
       KernelA12 = x**3.*(x*x*mu*(13.-41.*mu*mu)-4.*(mu+6.*mu**3) & 
            +x*(5.+9.*mu*mu+42.*mu**4.)) & 
            /14./(1.+x*x-2.*mu*x)**2.
    else if (id==id_A33) then
       KernelA12 = x**3.*(1.-7.*x*mu+6.*mu*mu)*(-2.*mu+x*(-1.+3.*mu*mu)) & 
            /14./(1.+x*x-2.*mu*x)**2.
    end if
    
  end function KernelA12
  !////////////////////////////////////////////////////////////////////////
  !This kernels correspond to B^n_ab in TNS(2010)
  function KernelB23(id,mu,x)
    real(dl), intent(IN) :: mu, x
    integer,  intent(IN) :: id
    real(dl) :: KernelB23

    if (id==id_B111) then
       KernelB23 = x**2.*(mu*mu-1.)/2./(1.+x*x-2.*mu*x)
    else if (id==id_B112) then
       KernelB23 = 3.*x**2*(mu*mu-1.)**2./8./(1.+x*x-2.*mu*x)
    else if (id==id_B121) then
       KernelB23 = 3.*x**4*(mu*mu-1.)**2./8./(1.+x*x-2.*mu*x)**2.
    else if (id==id_B122) then
       KernelB23 = 5.*x**4*(mu*mu-1.)**3./16./(1.+x*x-2.*mu*x)**2.
    else if (id==id_B211) then
       KernelB23 = x*(x+2.*mu-3.*x*mu*mu)/2./(1.+x*x-2.*mu*x)
    else if (id==id_B212) then
       KernelB23 = -3.*x*(mu*mu-1.)*(-x-2.*mu+5.*x*mu*mu)/4./(1.+x*x-2.*mu*x)
    else if (id==id_B221) then
       KernelB23 = 3.*x**2.*(mu*mu-1.)*(-2.+x*x+6.*x*mu-5.*x*x*mu*mu) & 
            /4./(1.+x*x-2.*mu*x)**2.
    else if (id==id_B222) then
       KernelB23 = -3.*x**2.*(mu*mu-1.)**2.*(6.-5.*x*x-30.*x*mu+35.*x*x*mu*mu) & 
            /16./(1.+x*x-2.*mu*x)**2.
    else if (id==id_B312) then
       KernelB23 = x*(4.*mu*(3.-5.*mu*mu)+x*(3.-30.*mu*mu+35.*mu**4)) & 
            /8./(1.+x*x-2.*mu*x)
    else if (id==id_B321) then
       KernelB23 = x*(-8.*mu+x*(-12.+36.*mu*mu+12.*x*mu*(3.-5.*mu*mu)+ & 
            x**2*(3.-30.*mu*mu+35.*mu**4))) & 
            /8./(1.+x*x-2.*mu*x)**2.
    else if (id==id_B322) then
       KernelB23 = 3.*x*(mu*mu-1.)*(-8.*mu+x*(-12.+60.*mu*mu+ & 
            20.*x*mu*(3.-7.*mu*mu)+5.*x*x*(1.-14.*mu*mu+21.*mu**4))) & 
            /16./(1.+x*x-2.*mu*x)**2.
    else if (id==id_B422) then
       KernelB23 = x*(8.*mu*(-3.+5.*mu*mu)-6.*x*(3.-30.*mu*mu+35.*mu**4) & 
            +6.*x*x*mu*(15.-70.*mu*mu+63*mu**4)+ & 
            x**3*(5.-21.*mu*mu*(5.-15.*mu*mu+11.*mu**4))) &
            /16./(1.+x*x-2.*mu*x)**2.
    else
       KernelB23 = 0.
    end if
    
  end function KernelB23
  !////////////////////////////////////////////////////////////////////////
  subroutine GetBiasCorrections(khvals,npoints,Pkbias, &
       alogkh,amatp_tot,addmatp_tot,num_k)
    integer, intent(in)   :: npoints, num_k
    real(dl), intent(in)  :: khvals(npoints), alogkh(num_k) 
    real(dl), intent(in)  :: amatp_tot(num_k), addmatp_tot(num_k)

    real(dl), intent(OUT) :: Pkbias(npoints,numid_bias)

    integer, parameter :: numx=500, nummu=30
    real(dl) :: tmin = -1._dl, tmax = 1._dl
    real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
    real(dl) :: aXorg(numx), aWXorg(numx), aMUorg(nummu), aWMUorg(nummu)
    real(dl) :: aX(numx), aWX(numx), aMU(nummu), aWMU(nummu)

    real(dl) :: P1, P2, P3, Pktmp(4)
    integer  :: ik, ik_start, ix, imu, id
    real(dl) :: kh_val, xmin, xmax, mumin, mumax, x, mu
    integer, parameter:: numr=501
    real(dl) :: alogr(numr), aIR(numr), addIR(numr)

    call ComputeIR(numr,alogr,aIR)
    call spline(alogr,aIR,numr,cllo,clhi,addIR) 
    !open(11,file='testIR.dat')
    !do ix = 1, numr
    !   write(11,'(2E15.5)') ar(ix), aIR(ix)
    !end do
    !stop "stop at GetOneLoopCorrections"

    call GauLeg(tmin,tmax,numx,aXorg,aWXorg)
    call GauLeg(tmin,tmax,nummu,aMUorg,aWMUorg)

    Pkbias = 0.
    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC) &
    !$OMP & PRIVATE(kh_val,x,xmin,xmax,aX,aWX, & 
    !$OMP &         mu,mumin,mumax,aMU,aWMU,P1,P2,P3,Pktmp)
    do ik = 1, npoints
       kh_val = khvals(ik)
       
       xmin = dlog(dexp(alogkh(1))/kh_val)
       xmax = dlog(20._dl/kh_val)

       !rescale Gaussi-Legendre roots and weights
       aX  = ((xmax-xmin)*aXorg+(xmax+xmin))/2._dl
       aWX = (xmax-xmin)/2._dl*aWXorg

       P1 = FuncPk(kh_val,alogkh,amatp_tot,addmatp_tot,num_k)

       !x-integration
       do ix = 1, numx
          x = dexp(aX(ix))
          
          P3 = FuncPk(kh_val*x,alogkh,amatp_tot,addmatp_tot,num_k)

          mumin = max(-1., (1.+x**2-dexp(xmax)**2)/2./x)
          mumax = min( 1., (1.+x**2-dexp(xmin)**2)/2./x)
          
          if (x>0.5_dl) mumax = 0.5_dl/x

          !rescale Gauss-Legendre roots and weights
          aMU  = ((mumax-mumin)*aMUorg+(mumax+mumin))/2._dl
          aWMU = (mumax-mumin)/2._dl*aWMUorg
          
          Pktmp = 0.
          !mu-integration
          do imu = 1, nummu
             mu = aMU(imu)
             
             P2 = FuncPk(kh_val*dsqrt(1.+x**2.-2.*x*mu), & 
                  alogkh,amatp_tot,addmatp_tot,num_k)

             do id = id_bs2d, id_bs2v
                Pktmp(id) = Pktmp(id) + aWMU(imu)*KernelBias(id,mu,x)*P2*P3
             end do

             Pktmp(id_b2s2) = Pktmp(id_b2s2) + aWMU(imu) & 
                  *x**2.*(2.*P2*P3*((mu-x)**2./(1.+x*x-2.*x*mu)-1./3.) & 
                  -2.*P2**2./3.-2.*P3**2./3.)/4.
             
             !b22 check
             !Pktmp(id_b2s2) = Pktmp(id_b2s2) + aWMU(imu) & 
             !     *x**2.*(2.*P2*P3-P2**2.-P3**2.)/4.

             Pktmp(id_bs22) = Pktmp(id_bs22) + aWMU(imu) & 
                  *x**2.*(2.*P2*P3*((mu-x)**2./(1.+x*x-2.*x*mu)-1./3.)**2. & 
                  -4.*P2**2./9.-4.*P3**2./9.)/4.

          end do

          PkBias(ik,id_bs2d:id_bs22) = PkBias(ik,id_bs2d:id_bs22) + aWX(ix)*x*Pktmp(id_bs2d:id_bs22)
      
          PkBias(ik,id_b3) = PkBias(ik,id_b3) & 
               + aWX(ix)*x**3.*P3*FuncIR(x,alogr,aIR,addIR,numr)
          
       end do

       !factor 2 comes from the half-plane 
       PkBias(ik,id_bs2d:id_bs22) = 2.* PkBias(ik,id_bs2d:id_bs22) & 
            * kh_val**3/(2.*pi)**2.
       
       !no factor 2 because of no mu-integration
       PkBias(ik,id_b3) = PkBias(ik,id_b3)*kh_val**3/(2.*pi**2.)*P1
      
    end do
    !$OMP END PARALLEL DO

  end subroutine GetBiasCorrections
  !////////////////////////////////////////////////////////////////////////
  FUNCTION KernelBias(id,mu,x)
    real(dl), intent(IN) :: mu, x
    integer,  intent(IN) :: id
    real(dl) :: KernelBias

    if (id==id_bs2d) then
       KernelBias = (3.*x+7.*mu-10.*x*mu**2.)/(1.+x*x-2.*x*mu)/14./x & 
            *((mu-x)**2./(1.+x*x-2.*x*mu)-1./3.)*x**2.
    else if (id==id_bs2v) then
       KernelBias = (-x+7.*mu-6.*x*mu**2)/(1.+x*x-2.*x*mu)/14./x & 
            *((mu-x)**2./(1.+x*x-2.*x*mu)-1./3.)*x**2.
    end if
  END FUNCTION KernelBias
  !////////////////////////////////////////////////////////////////////////
  FUNCTION FuncIR(r_eval,alogr,aIR,addIR,numr)
    real(dl), intent(IN) :: r_eval
    integer, intent(IN)  :: numr
    real(dl), intent(IN) :: alogr(numr), aIR(numr), addIR(numr)
    real(dl) :: logr, FuncIR
    integer llo,lhi
    real(dl) dp, ho, a0, b0
    integer, save :: i_last = 1

    logr = dlog(r_eval)

    if (logr<alogr(1)) then
       FuncIR = 1.
    else if  (logr>alogr(numr)) then
       FuncIR = 0.
    else 
       llo=min(i_last,numr)
       do while (alogr(llo) > logr)
          llo=llo-1
       end do
       do while (alogr(llo+1) < logr)
          llo=llo+1
       end do
       i_last =llo  
       lhi=llo+1
       ho=alogr(lhi)-alogr(llo) 
       a0=(alogr(lhi)-logr)/ho
       b0=1-a0
       
       FuncIR = a0*aIR(llo)+b0*aIR(lhi) + &
            ((a0**3-a0)*addIR(llo) &
            +(b0**3-b0)*addIR(lhi))*ho**2./6.
    end if

  END FUNCTION FuncIR
  !////////////////////////////////////////////////////////////////////////
  SUBROUTINE ComputeIR(numr,alogr,aIR)
    integer, intent(IN) :: numr
    real(dl), intent(OUT) :: alogr(numr), aIR(numr)
    integer, parameter  :: nummu=30
    real(dl) :: r, mu, aMU(nummu), aWMU(nummu)
    real(dl) :: mumin = -1._dl, mumax = 1._dl
    integer  :: ir, imu

    call GauLeg(mumin,mumax,nummu,aMU,aWMU)

    alogr = 0.
    aIR = 0.
    do ir = 1, numr
       r = 10.**(-2.+4./(numr-1)*(ir-1))
       alogr(ir) = dlog(r)
 
       do imu = 1, nummu
          aIR(ir) = aIR(ir) + aWMU(imu)*KernelIR(r,aMU(imu))
       end do

       aIR(ir) = -21./8.*(aIR(ir)/2.-34./63.)
    end do

  END SUBROUTINE ComputeIR
  !////////////////////////////////////////////////////////////////////////
  FUNCTION KernelIR(r,mu)
    real(dl), intent(IN) :: r, mu
    real(dl) :: KernelIR

    KernelIR = (5./7.-mu/2.*(r+1./r)+2.*mu**2./7.) & 
         *((mu-r)**2./(1.+r**2.-2.*r*mu)-1./3.)
  END FUNCTION KernelIR
  !////////////////////////////////////////////////////////////////////////
  !added by S.S. on Dec/19/2012
  subroutine GetA01Corrections(khvals,npoints,A01, &
       alogkh,amatp_tot,addmatp_tot,num_k)
    integer, intent(in)   :: npoints, num_k
    real(dl), intent(in)  :: khvals(npoints), alogkh(num_k) 
    real(dl), intent(in)  :: amatp_tot(num_k), addmatp_tot(num_k)

    real(dl), intent(OUT) :: A01(npoints,numid_A01)

    integer, parameter :: numx=500, nummu=30
    real(dl) :: tmin = -1._dl, tmax = 1._dl
    real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
    real(dl) :: aXorg(numx), aWXorg(numx), aMUorg(nummu), aWMUorg(nummu)
    real(dl) :: aX(numx), aWX(numx), aMU(nummu), aWMU(nummu)

    real(dl) :: P1, P2, P3, B12(numid_A01), B31(numid_A01), B23(numid_A01)
    integer  :: ik, ik_start, ix, imu, id
    real(dl) :: kh_val, xmin, xmax, mumin, mumax, x, mu, zin, f

    !zin = 0.
    !call ComputeGrowth(zin, f)
   
    call GauLeg(tmin,tmax,numx,aXorg,aWXorg)
    call GauLeg(tmin,tmax,nummu,aMUorg,aWMUorg)

    A01 = 0.
    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC) &
    !$OMP & PRIVATE(kh_val,x,xmin,xmax,aX,aWX, & 
    !$OMP &         mu,mumin,mumax,aMU,aWMU,P1,P2,P3,B12,B31,B23)
    do ik = 1, npoints
       kh_val = khvals(ik)
       
       xmin = dlog(dexp(alogkh(1))/kh_val)
       xmax = dlog(20._dl/kh_val)

       !rescale Gaussi-Legendre roots and weights
       aX  = ((xmax-xmin)*aXorg+(xmax+xmin))/2._dl
       aWX = (xmax-xmin)/2._dl*aWXorg

       P1 = FuncPk(kh_val,alogkh,amatp_tot,addmatp_tot,num_k)

       !x-integration
       do ix = 1, numx
          x = dexp(aX(ix))
          
          P3 = FuncPk(kh_val*x,alogkh,amatp_tot,addmatp_tot,num_k)

          mumin = max(-1., (1.+x**2-dexp(xmax)**2)/2./x)
          mumax = min( 1., (1.+x**2-dexp(xmin)**2)/2./x)
          if (x>0.5_dl) mumax = 0.5_dl/x

          !rescale Gauss-Legendre roots and weights
          aMU  = ((mumax-mumin)*aMUorg+(mumax+mumin))/2._dl
          aWMU = (mumax-mumin)/2._dl*aWMUorg
          
          B12 = 0.
          B31 = 0.
          B23 = 0.
          !mu-integration
          do imu = 1, nummu
             mu = aMU(imu)
             
             P2 = FuncPk(kh_val*dsqrt(1.+x**2.-2.*x*mu), & 
                  alogkh,amatp_tot,addmatp_tot,num_k)

             do id = id_Bb1, numid_A01
                B12(id) = B12(id) + aWMU(imu)*KernelA01_B12(id,mu,x)*P1*P2
                B31(id) = B31(id) + aWMU(imu)*KernelA01_B31(id,mu,x)*P3*P1
                B23(id) = B23(id) + aWMU(imu)*KernelA01_B23(id,mu,x)*P2*P3
             end do

          end do

          do id = id_Bb1, numid_A01
             A01(ik,id)  = A01(ik,id)  + aWX(ix)*x*(B12(id)+B31(id)+B23(id))
          end do
         
       end do

       !factor 2 comes from the half-plane 
       A01(ik,:) = 2.* A01(ik,:) * kh_val**3/(2.*pi)**2.
      
    end do
    !$OMP END PARALLEL DO

  end subroutine GetA01Corrections
  !////////////////////////////////////////////////////////////////////////
  function KernelA01_B12(id,mu,x)
    real(dl), intent(IN) :: mu, x
    integer, intent(IN) :: id
    real(dl) :: KernelA01_B12
    
    if (id==id_Bb1) then
       KernelA01_B12 = (-x+3.*mu+11.*x**2.*mu-5.*x*mu**2.-10.*mu**3. & 
            -18.*x**2.*mu**3.+20.*x*mu**4.) & 
            *x**3./14./(1.+x*x-2.*mu*x)**2.
    !else if (id==id_tBb2) then
    !   KernelA01_B12 = x**2.*(1.-x*mu)/2./(1.+x*x-2.*mu*x)
    !else if (id==id_tBbs2) then
    !   KernelA01_B12 = x**2.*(1.-x*mu)*(2.-4.*x*mu-x**2.+3.*x**2.*mu**2.) & 
    !        /6./(1.+x*x-2.*mu*x)
    else
       KernelA01_B12 = 0.
    end if
       
  end function KernelA01_B12
  !////////////////////////////////////////////////////////////////////////
  function KernelA01_B31(id,mu,x)
    real(dl), intent(IN) :: mu, x
    integer, intent(IN) :: id
    real(dl) :: KernelA01_B31
    
    if (id==id_Bb1) then
       KernelA01_B31 = (-7.*mu**2.-mu*x+22.*mu**3.*x+10.*x**2.-15.*mu**2.*x**2. & 
            -16.*mu**4.*x**2.-11.*mu*x**3.+18.*mu**3.*x**3.)& 
            /14./(1.+x*x-2.*mu*x)
    !else if (id==id_tBb2) then
    !   KernelA01_B31 = x*mu/2.
    !else if (id==id_tBbs2) then
    !   KernelA01_B31 = x*mu*(-1.+3*mu**2.)/6.
    else
       KernelA01_B31 = 0.
    end if
    
  end function KernelA01_B31
  !////////////////////////////////////////////////////////////////////////
  function KernelA01_B23(id,mu,x)
    real(dl), intent(IN) :: mu, x
    integer, intent(IN) :: id
    real(dl) :: KernelA01_B23
    
    if (id==id_Bb1) then
       KernelA01_B23 = (-x-mu+2.*x*mu**2.)*(-3.*x-7.*mu+10.*x*mu**2.) & 
            /14./(1.+x*x-2.*mu*x)**2.
    else if (id==id_Bb2) then
       KernelA01_B23 = x*(x+mu-2.*x*mu**2.)/(1.+x*x-2.*mu*x)/2.
    else if (id==id_Bbs2) then
       KernelA01_B23 = x*(-1.+2.*x**2.-4.*x*mu+3.*mu**2.) & 
            *(x+mu-2.*x*mu**2.)/3./(1.+x*x-2.*mu*x)**2./2.
    else
       KernelA01_B23 = 0.
    end if
    
  end function KernelA01_B23
  !////////////////////////////////////////////////////////////////////////
  FUNCTION D1(zred)
    !Linear growth function
    real(dl), intent(IN)  :: zred
    real(dl) :: D1
    real(dl) :: zred1, zi, zz, a, b, c, g0, gz

    if (wa.NE.0) then
       stop "D1 supports only in the case of constant w"
    end if

    a = -1./(3.*w0)
    b = (w0-1.)/(2.*w0)
    c = 1.-5./(6.*w0)
 
    zred1 = 1. + zred
    zi = - CP%omegav /(CP%omegab+CP%omegac+CP%omegan)
    zz = zi * zred1**(3.*w0)
    
    call HYGFX(a,b,c,zi, g0)  
    call HYGFX(a,b,c,zz, gz)  
    
    D1 = (gz/g0) / zred1

  END FUNCTION D1
  !////////////////////////////////////////////////////////////////////////
  SUBROUTINE ComputeGrowth(zin, f)
    real(dl), intent(IN) :: zin
    real(dl), intent(OUT) :: f
    real(dl) :: zred1, zi, zz, a, b, c, g1, g2

    !write(*,*) " "
    !write(*,*) "computing dlnD/dlna..."

    if (wa.NE.0) then
       stop "ComputeGrowth supports only in the case of constant w"
    end if

    zred1 = 1. + zin
    zi = - CP%omegav /(CP%omegab+CP%omegac+CP%omegan)
    zz = zi * zred1**(3.*w0)
    
    a = 1.5 - 1./(2.*w0)
    b = 1.  - 1./(3.*w0)
    c = 2.  - 5./(6.*w0)

    call HYGFX(a,b,c,zz,g1)

    b = -1./(3.*w0)
    a = (-1.+w0)/(2.*w0)
    c = 1.-5./(6.*w0)
    
    call HYGFX(a,b,c,zz,g2)
    
    f = 1. + 3.*(-1.+w0)/(6.*w0-5.)*zz*g1/g2

    !write(*,'(4X,A,F8.5)') "input redshift = ", zin
    !write(*,'(4X,A,F8.5)') "f = dlnD/dlna  = ", f

    write(*,*) '   z = ', real(zin), & 
         " f =   dlnD/dlna    =", real(f)

  END SUBROUTINE ComputeGrowth
  !////////////////////////////////////////////////////////////////////////
  SUBROUTINE CalcLinearSigmaV(alogkh,amatp,addmatp,num_k,sigma_v)
    !sigma_v^2 = 1/6pi^2 * Integral dq P^L(q)
    integer, intent(IN)   :: num_k 
    real(dl), intent(IN)  :: alogkh(num_k), amatp(num_k), addmatp(num_k)
    real(dl), intent(OUT) :: sigma_v
    
    real(dl) :: kmin = 1.e-5_dl, kmax = 20._dl
    real(dl) :: a = -1._dl, b = 1._dl
    integer, parameter :: ny = 300
    integer :: i, iy
    real(dl) :: y(ny), w(ny)
    real(dl) :: ymin, ymax, scale, sigmav2
    
    ymin = dlog(kmin)
    ymax = dlog(kmax)
    
    call GauLeg(a,b,ny,y,w)
    !rescale
    do i = 1, ny
       y(i) = ((ymax-ymin)*y(i)+ymax+ymin)/2._dl
    end do
    
    scale = (ymax-ymin)/2._dl
    
    do iy = 1, ny
       sigmav2 = sigmav2 + w(iy)*dexp(y(iy)) & 
            *FuncPk(dexp(y(iy)),alogkh,amatp,addmatp,num_k)
    end do
    
    sigmav2 = sigmav2*scale/(6._dl*pi*pi)

    sigma_v = dsqrt(sigmav2)
    
  end subroutine CalcLinearSigmaV
  !////////////////////////////////////////////////////////////////////////
  SUBROUTINE CalcMonopole(Pk,num_k,num_id,f,sigmav,Ps0)
    integer, intent(IN)   :: num_k, num_id
    real(dl), intent(IN)  :: f, sigmav
    real(dl), intent(IN)  :: Pk(num_k,num_id)
    real(dl), intent(OUT) :: Ps0(num_k,id_mu8AB)
    !Ps0 is output monopole power spectrum
    
    integer  :: ik
    real(dl) :: kh, x
    real(dl) :: A0k, A2k, A4k, A6k, A8k

    !check data consistensy
    if (CP%NonLinear_RSDmodel==NLKaiser_AB_FoG) then
       if (num_id.NE.id_Pkb22+id_A33+id_B422) then 
          stop "inconsistent data in CalcMonopole"
       end if
    else
       if (num_id.NE.id_Pkb22) then 
          stop "inconsistent data in CalcMonopole"
       end if
    end if

    Ps0 = 0.
    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC) &
    !$OMP & PRIVATE(kh,x,A0k,A2k,A4k,A6k,A8k)
    do ik = 1, num_k
       kh = Pk(ik,id_kh)
       
       if (CP%NonLinear_RSDmodel==3.or.CP%NonLinear_RSDmodel==4) then
          !first compute prefactor assuming gaussian FoG factor
          !1/2 Int dmu exp[-x*mu^2] mu^2i = 1/2 x^-i-1/2 LImcompleteGamma(i+1/2,x)
          x = f**2. * sigmav**2. * kh**2. 

          A0k = x**(-0.5)*LImcompleteGamma(0.5_dl,x)/2._dl
          A2k = x**(-1.5)*LImcompleteGamma(1.5_dl,x)/2._dl
          A4k = x**(-2.5)*LImcompleteGamma(2.5_dl,x)/2._dl
          A6k = x**(-3.5)*LImcompleteGamma(3.5_dl,x)/2._dl
          A8k = x**(-4.5)*LImcompleteGamma(4.5_dl,x)/2._dl
       else
          A0k = 1.
          A2k = 1./3.
          A4k = 1./5.
          A6k = 0.
          A8k = 0.
       end if

       !kh
       Ps0(ik,id_kh) = kh
       !mu^0 terms 
       Ps0(ik,id_mu0dd)  = A0k*Pk(ik,id_PkNL_dd)
       Ps0(ik,id_mu0b2d) = A0k*Pk(ik,id_Pkb2_dd)
       Ps0(ik,id_mu0b22) = A0k*Pk(ik,id_Pkb22)
       !mu^2 terms
       Ps0(ik,id_mu2dv)  = A2k*2._dl*f*Pk(ik,id_PkNL_dv)
       Ps0(ik,id_mu2b2v) = A2k*2._dl*f*Pk(ik,id_Pkb2_dv)
       !mu^4 terms
       Ps0(ik,id_mu4vv)  = A4k*f**2.*Pk(ik,id_PkNL_vv)

       if (CP%NonLinear_RSDmodel==NLKaiser_AB_FoG) then
          !mu^2 terms
          Ps0(ik,id_mu2AB)  = A2k*( f*Pk(ik,id_Pkb22+id_A11) & 
               +f**2.*Pk(ik,id_Pkb22+id_A12) & 
               +f**2.*Pk(ik,id_Pkb22+id_A33+id_B111) & 
               -f**3.*Pk(ik,id_Pkb22+id_A33+id_B112) &
               -f**3.*Pk(ik,id_Pkb22+id_A33+id_B121) &
               +f**4.*Pk(ik,id_Pkb22+id_A33+id_B122) )
          !mu^4 terms
          Ps0(ik,id_mu4AB)  = A4k*( f**2.*Pk(ik,id_Pkb22+id_A22) & 
               +f**3.*Pk(ik,id_Pkb22+id_A23) & 
               +f**2.*Pk(ik,id_Pkb22+id_A33+id_B211) & 
               -f**3.*Pk(ik,id_Pkb22+id_A33+id_B212) &
               -f**3.*Pk(ik,id_Pkb22+id_A33+id_B221) &
               +f**4.*Pk(ik,id_Pkb22+id_A33+id_B222) )
          !mu^6 terms 
          Ps0(ik,id_mu6AB)  = A6k*( f**3.*Pk(ik,id_Pkb22+id_A33) & 
               -f**3.*Pk(ik,id_Pkb22+id_A33+id_B312) &
               -f**3.*Pk(ik,id_Pkb22+id_A33+id_B321) &
               +f**4.*Pk(ik,id_Pkb22+id_A33+id_B322) )
          !mu^8 terms 
          Ps0(ik,id_mu8AB)  = A8k*f**4.*Pk(ik,id_Pkb22+id_A33+id_B422)
       end if
       
    end do
    !$OMP END PARALLEL DO

  END SUBROUTINE CalcMonopole
  !////////////////////////////////////////////////////////////////////////
end module Transfer
!==========================================================================
module ThermoData
  use ModelData
  implicit none
  private
  integer,parameter :: nthermo=20000
  
  real(dl) tb(nthermo),cs2(nthermo),xe(nthermo)
  real(dl) dcs2(nthermo)
  real(dl) dotmu(nthermo), ddotmu(nthermo)
  real(dl) sdotmu(nthermo),emmu(nthermo)
  real(dl) demmu(nthermo)
  real(dl) dddotmu(nthermo),ddddotmu(nthermo)
  real(dl) tauminn,dlntau,Maxtau
  real(dl), dimension(:), allocatable :: vis,dvis,ddvis,expmmu,dopac, opac
  
  real(dl) :: tight_tau, actual_opt_depth
  !Times when 1/(opacity*tau) = 0.01, for use switching tight coupling approximation
  real(dl) :: matter_verydom_tau
  real(dl) :: r_drag0, z_star, z_drag  !!JH for updated BAO likelihood.
  public thermo,inithermo,vis,opac,expmmu,dvis,dopac,ddvis, tight_tau,&
       Thermo_OpacityToTime,matter_verydom_tau, ThermoData_Free,&
       z_star, z_drag !!JH for updated BAO likelihood.
contains
  !////////////////////////////////////////////////////////////////////////
  subroutine thermo(tau,cs2b,opacity, dopacity)
    !Compute unperturbed sound speed squared,
    !and ionization fraction by interpolating pre-computed tables.
    !If requested also get time derivative of opacity
    implicit none
    real(dl) tau,cs2b,opacity
    real(dl), intent(out), optional :: dopacity
    
    integer i
    real(dl) d
    
    d=log(tau/tauminn)/dlntau+1._dl
    i=int(d)
    d=d-i
    if (i < 1) then
       !Linear interpolation if out of bounds (should not occur).
       cs2b=cs2(1)+(d+i-1)*dcs2(1)
       opacity=dotmu(1)+(d-1)*ddotmu(1)
       stop 'thermo out of bounds'
    else if (i >= nthermo) then
       cs2b=cs2(nthermo)+(d+i-nthermo)*dcs2(nthermo)
       opacity=dotmu(nthermo)+(d-nthermo)*ddotmu(nthermo)
       if (present(dopacity)) then
          dopacity = 0
          stop 'thermo: shouldn''t happen'
       end if
    else
       !Cubic spline interpolation.
       cs2b=cs2(i)+d*(dcs2(i)+d*(3*(cs2(i+1)-cs2(i))  &
            -2*dcs2(i)-dcs2(i+1)+d*(dcs2(i)+dcs2(i+1)  &
            +2*(cs2(i)-cs2(i+1)))))
       opacity=dotmu(i)+d*(ddotmu(i)+d*(3*(dotmu(i+1)-dotmu(i)) &
            -2*ddotmu(i)-ddotmu(i+1)+d*(ddotmu(i)+ddotmu(i+1) &
            +2*(dotmu(i)-dotmu(i+1)))))
       
       if (present(dopacity)) then
          
          dopacity=(ddotmu(i)+d*(dddotmu(i)+d*(3*(ddotmu(i+1)  &
               -ddotmu(i))-2*dddotmu(i)-dddotmu(i+1)+d*(dddotmu(i) &
               +dddotmu(i+1)+2*(ddotmu(i)-ddotmu(i+1))))))/(tau*dlntau)
          
       end if
    end if
  end subroutine thermo
  !////////////////////////////////////////////////////////////////////////
  function Thermo_OpacityToTime(opacity)
    real(dl), intent(in) :: opacity
    integer j
    real(dl) Thermo_OpacityToTime
    !Do this the bad slow way for now..
    !The answer is approximate
    j =1
    do while(dotmu(j)> opacity)
       j=j+1
    end do
    
    Thermo_OpacityToTime = exp((j-1)*dlntau)*tauminn
    
  end function Thermo_OpacityToTime
  !////////////////////////////////////////////////////////////////////////
  subroutine inithermo(taumin,taumax)
    !  Compute and save unperturbed baryon temperature and ionization fraction
    !  as a function of time.  With nthermo=10000, xe(tau) has a relative 
    ! accuracy (numerical integration precision) better than 1.e-5.
    use constants
    use precision
    use ModelParams
    use MassiveNu
    real(dl) taumin,taumax
    
    
    real(dl) tau01,adot0,a0,a02,x1,x2,barssc,dtau
    real(dl) xe0,tau,a,a2
    real(dl) adot,tg0,ahalf,adothalf,fe,thomc,thomc0,etc,a2t
    real(dl) dtbdla,vfi,cf1,maxvis, vis
    integer ncount,i,j1,j2,iv,ns
    real(dl) spline_data(nthermo)
    real(dl) last_dotmu
    real(dl) dtauda  !diff of tau w.CP%r.t a and integration
    external dtauda
    real(dl) a_verydom
     
    call Recombination_Init(CP%Recomb, CP%omegac, CP%omegab,CP%Omegan, & 
         CP%Omegav, CP%h0,CP%tcmb,CP%yhe)
    !almost all the time spent here
    
    Maxtau=taumax
    tight_tau = 0
    actual_opt_depth = 0
    ncount=0
    thomc0= Compton_CT * CP%tcmb**4 
    !thomc0=5.0577d-8*CP%tcmb**4
    
    tauminn=0.05d0*taumin
    dlntau=log(CP%tau0/tauminn)/(nthermo-1)
    last_dotmu = 0
    
    matter_verydom_tau = 0
    a_verydom = AccuracyBoost*5*(grhog+grhornomass)/(grhoc+grhob)
    
    !  Initial conditions: assume radiation-dominated universe.
    tau01=tauminn
    adot0=adotrad
    a0=adotrad*tauminn
    a02=a0*a0
    !  Assume that any entropy generation occurs before tauminn.
    !  This gives wrong temperature before pair annihilation, but
    !  the error is harmless.
    tb(1)=CP%tcmb/a0
    xe0=1._dl
    x1=0._dl
    x2=1._dl
    xe(1)=xe0+0.25d0*CP%yhe/(1._dl-CP%yhe)*(x1+2*x2)
    barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*xe(1))
    cs2(1)=4._dl/3._dl*barssc*tb(1)
    dotmu(1)=xe(1)*akthom/a02
    sdotmu(1)=0
    
    do i=2,nthermo
       tau=tauminn*exp((i-1)*dlntau)
       dtau=tau-tau01
       !  Integrate Friedmann equation using inverse trapezoidal rule.
       
       a=a0+adot0*dtau
       a2=a*a
       
       adot=1/dtauda(a)
       
       if (matter_verydom_tau ==0 .and. a > a_verydom) then
          matter_verydom_tau = tau  
       end if
       
       a=a0+2._dl*dtau/(1._dl/adot0+1._dl/adot)         
       !  Baryon temperature evolution: adiabatic except for Thomson cooling.
       !  Use  quadrature solution.
       ! This is redundant as also calculated in REFCAST, but agrees 
       !well before reionization
       tg0=CP%tcmb/a0
       ahalf=0.5d0*(a0+a)
       adothalf=0.5d0*(adot0+adot)
       !  fe=number of free electrons divided by total number of free baryon
       !  particles (e+p+H+He).  Evaluate at timstep i-1 for convenience; if
       !  more accuracy is required (unlikely) then this can be iterated with
       !  the solution of the ionization equation.
       fe=(1._dl-CP%yhe)*xe(i-1)/(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*xe(i-1))
       thomc=thomc0*fe/adothalf/ahalf**3
       etc=exp(-thomc*(a-a0))
       a2t=a0*a0*(tb(i-1)-tg0)*etc-CP%tcmb/thomc*(1._dl-etc)
       tb(i)=CP%tcmb/a+a2t/(a*a)
       
       ! If there is re-ionization, smoothly increase xe to the 
       ! requested value.
       if (CP%Reion%Reionization .and. tau > CP%ReionHist%tau_start) then
          if(ncount == 0) then
             ncount=i-1
          end if
          
          xe(i) = Reionization_xe(a, tau, xe(ncount))
          !print *,1/a-1,xe(i)
          if (CP%AccurateReionization .and. FeedbackLevel > 0) then
             dotmu(i)=(Recombination_xe(a) - xe(i))*akthom/a2
             
             if (last_dotmu /=0) then
                actual_opt_depth = actual_opt_depth & 
                     - 2._dl*dtau/(1._dl/dotmu(i)+1._dl/last_dotmu)
             end if
             last_dotmu = dotmu(i) 
          end if
          
       else
          xe(i)=Recombination_xe(a)
       end if
       
       
       !  Baryon sound speed squared (over c**2).
       dtbdla=-2._dl*tb(i)-thomc*adothalf/adot*(a*tb(i)-CP%tcmb)
       barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*xe(i))
       cs2(i)=barssc*tb(i)*(1-dtbdla/tb(i)/3._dl)
       
       
       ! Calculation of the visibility function
       dotmu(i)=xe(i)*akthom/a2
       
       if (tight_tau==0 .and. 1/(tau*dotmu(i)) > 0.005) tight_tau = tau !0.005
       !Tight coupling switch time when k/opacity is smaller than 1/(tau*opacity)
       
       if (tau < 0.001) then
          sdotmu(i)=0
       else
          sdotmu(i)=sdotmu(i-1)+2._dl*dtau/(1._dl/dotmu(i)+1._dl/dotmu(i-1))
       end if
       
       a0=a
       tau01=tau
       adot0=adot
    end do !i
    
    if (CP%Reion%Reionization .and. (xe(nthermo) < 0.999d0)) then
       write(*,*)'Warning: xe at redshift zero is < 1'
       write(*,*) 'Check input parameters an Reionization_xe'
       write(*,*) 'function in the Reionization module'
    end if
    
    do j1=1,nthermo
       if (sdotmu(j1) - sdotmu(nthermo)< -69) then
          emmu(j1)=1.d-30
       else
          emmu(j1)=exp(sdotmu(j1)-sdotmu(nthermo))
          if (.not. CP%AccurateReionization .and. &
               actual_opt_depth==0 .and. xe(j1) < 1e-3) then
             actual_opt_depth = -sdotmu(j1)+sdotmu(nthermo) 
          end if
       end if
    end do
    
    if (CP%AccurateReionization .and. FeedbackLevel > 0) then                         
       write(*,'("Reion opt depth      = ",f7.4)') actual_opt_depth
    end if
    
    iv=0
    vfi=0._dl
    ! Getting the starting and finishing times for decoupling and time 
    !of maximum visibility
    if (ncount == 0) then
       cf1=1._dl
       ns=nthermo
    else
       cf1=exp(sdotmu(nthermo)-sdotmu(ncount))
       ns=ncount
    end if
    maxvis = 0
    do j1=1,ns
       vis = emmu(j1)*dotmu(j1)
       tau = tauminn*exp((j1-1)*dlntau)
       vfi=vfi+vis*cf1*dlntau*tau
       if ((iv == 0).and.(vfi > 1.0d-7/AccuracyBoost)) then  
          taurst=9._dl/10._dl*tau
          iv=1
       elseif (iv == 1) then 
          if (vis > maxvis) then
             maxvis=vis
             tau_maxvis = tau
          end if
          if (vfi > 0.995) then 
             taurend=tau 
             iv=2
             exit
          end if
       end if
    end do
    
    if (iv /= 2) then
       stop 'inithermo: failed to find end of recombination'
       !              taurend=1.5d0*(tauminn*exp((ncount-1)*dlntau))
    end if
    
    ! Calculating the timesteps during recombination.
    
    if (CP%WantTensors) then
       dtaurec=min(dtaurec,taurst/160)/AccuracyBoost 
    else
       dtaurec=min(dtaurec,taurst/40)/AccuracyBoost  
       if (do_bispectrum .and. hard_bispectrum) dtaurec = dtaurec / 4
    end if
    
    if (CP%Reion%Reionization) taurend=min(taurend,CP%ReionHist%tau_start)
    
    if (DebugMsgs) then
       write (*,*) 'taurst, taurend = ', taurst, taurend
    end if
    
    call splini(spline_data,nthermo)
    call splder(cs2,dcs2,nthermo,spline_data)
    call splder(dotmu,ddotmu,nthermo,spline_data)
    call splder(ddotmu,dddotmu,nthermo,spline_data)  
    call splder(dddotmu,ddddotmu,nthermo,spline_data)
    call splder(emmu,demmu,nthermo,spline_data)
    
    call SetTimeSteps
    
    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC) 
    do j2=1,TimeSteps%npoints
       call DoThermoSpline(j2,TimeSteps%points(j2))
    end do
    !$OMP END PARALLEL DO 
    
    if (CP%want_zdrag .or. CP%want_zstar) then !JH: calculate exact zstar and/or zdrag
       
       r_drag0 = 3.d0/4.d0*CP%omegab*grhom/grhog
       
       if (CP%want_zstar) call find_z(optdepth,z_star)
       if (CP%want_zdrag) call find_z(dragoptdepth,z_drag)
       
    end if
  end subroutine inithermo
  !////////////////////////////////////////////////////////////////////////
  subroutine SetTimeSteps
    real(dl) dtau0
    integer nri0, nstep
    
    call Ranges_Init(TimeSteps)
    
    call Ranges_Add_delta(TimeSteps, taurst, taurend, dtaurec)
    
    ! Calculating the timesteps after recombination
    if (CP%WantTensors) then
       dtau0=max(taurst/40,Maxtau/2000._dl/AccuracyBoost)
    else       
       dtau0=Maxtau/500._dl/AccuracyBoost 
       if (do_bispectrum) dtau0 = dtau0/3 
       !Don't need this since adding in Limber on small scales
       !  if (CP%DoLensing) dtau0=dtau0/2 
       !  if (CP%AccurateBB) dtau0=dtau0/3 !Need to get C_Phi accurate on small scales
    end if
    
    call Ranges_Add_delta(TimeSteps,taurend, CP%tau0, dtau0)
    
    if (CP%Reion%Reionization) then
       
       nri0=int(Reionization_timesteps(CP%ReionHist)*AccuracyBoost) 
       !Steps while reionization going from zero to maximum
       call Ranges_Add(TimeSteps,CP%ReionHist%tau_start,CP%ReionHist%tau_complete, & 
            nri0) 
       
    end if
    
    !Create arrays out of the region information.
    call Ranges_GetArray(TimeSteps)
    nstep = TimeSteps%npoints
    
    if (allocated(vis)) then
       deallocate(vis,dvis,ddvis,expmmu,dopac, opac)
    end if
    allocate(vis(nstep),dvis(nstep),ddvis(nstep),expmmu(nstep), & 
         dopac(nstep),opac(nstep))
    
    if (DebugMsgs .and. FeedbackLevel > 0) write(*,*) 'Set ',nstep, ' time steps'
    
  end subroutine SetTimeSteps
  !////////////////////////////////////////////////////////////////////////
  subroutine ThermoData_Free
    if (allocated(vis)) then
       deallocate(vis,dvis,ddvis,expmmu,dopac, opac)
    end if
    call Ranges_Free(TimeSteps)
    
  end subroutine ThermoData_Free
  !////////////////////////////////////////////////////////////////////////
  subroutine DoThermoSpline(j2,tau)
    integer j2,i
    real(dl) d,ddopac,tau
    
    !     Cubic-spline interpolation.
    d=log(tau/tauminn)/dlntau+1._dl
    i=int(d)
    
    d=d-i
    if (i < nthermo) then
       opac(j2)=dotmu(i)+d*(ddotmu(i)+d*(3._dl*(dotmu(i+1)-dotmu(i)) &
            -2._dl*ddotmu(i)-ddotmu(i+1)+d*(ddotmu(i)+ddotmu(i+1) &
            +2._dl*(dotmu(i)-dotmu(i+1)))))
       dopac(j2)=(ddotmu(i)+d*(dddotmu(i)+d*(3._dl*(ddotmu(i+1)  &
            -ddotmu(i))-2._dl*dddotmu(i)-dddotmu(i+1)+d*(dddotmu(i) &
            +dddotmu(i+1)+2._dl*(ddotmu(i)-ddotmu(i+1))))))/(tau &
            *dlntau)
       ddopac=(dddotmu(i)+d*(ddddotmu(i)+d*(3._dl*(dddotmu(i+1) &
            -dddotmu(i))-2._dl*ddddotmu(i)-ddddotmu(i+1)  &
            +d*(ddddotmu(i)+ddddotmu(i+1)+2._dl*(dddotmu(i) &
            -dddotmu(i+1)))))-(dlntau**2)*tau*dopac(j2)) &
            /(tau*dlntau)**2
       expmmu(j2)=emmu(i)+d*(demmu(i)+d*(3._dl*(emmu(i+1)-emmu(i)) &
            -2._dl*demmu(i)-demmu(i+1)+d*(demmu(i)+demmu(i+1) &
            +2._dl*(emmu(i)-emmu(i+1)))))
       
       vis(j2)=opac(j2)*expmmu(j2)
       dvis(j2)=expmmu(j2)*(opac(j2)**2+dopac(j2))
       ddvis(j2)=expmmu(j2)*(opac(j2)**3+3*opac(j2)*dopac(j2)+ddopac)
    else
       opac(j2)=dotmu(nthermo)
       dopac(j2)=ddotmu(nthermo)
       ddopac=dddotmu(nthermo)
       expmmu(j2)=emmu(nthermo)
       vis(j2)=opac(j2)*expmmu(j2)
       dvis(j2)=expmmu(j2)*(opac(j2)**2+dopac(j2))
       ddvis(j2)=expmmu(j2)*(opac(j2)**3+3._dl*opac(j2)*dopac(j2)+ddopac)
       
    end if
  end subroutine DoThermoSpline
  !////////////////////////////////////////////////////////////////////////
  !JH: functions and subroutines for calculating z_star and z_drag
  function doptdepth_dz(z)
    real(dl) :: doptdepth_dz
    real(dl), intent(in) :: z
    real(dl) :: a
    real(dl) :: dtauda
    external dtauda
    
    a = 1._dl/(1._dl+z)
    
    !ignoring reionisation, not relevant for distance measures
    doptdepth_dz = Recombination_xe(a)*akthom*dtauda(a)
    
  end function doptdepth_dz
  !////////////////////////////////////////////////////////////////////////
  function optdepth(z)
    real(dl) :: rombint2
    external rombint2
    real(dl) optdepth
    real(dl),intent(in) :: z
    
    optdepth = rombint2(doptdepth_dz, 0.d0, z, 1d-5, 20, 100)
    
  end function optdepth
  !////////////////////////////////////////////////////////////////////////
  function ddragoptdepth_dz(z)
    real(dl) :: ddragoptdepth_dz
    real(dl), intent(in) :: z
    real(dl) :: a
    real(dl) :: dtauda
    external dtauda
    
    a = 1._dl/(1._dl+z)
    ddragoptdepth_dz = doptdepth_dz(z)/r_drag0/a
    
  end function ddragoptdepth_dz
  !////////////////////////////////////////////////////////////////////////
  function dragoptdepth(z)
    real(dl) :: rombint2
    external rombint2
    real(dl) dragoptdepth
    real(dl),intent(in) :: z
    
    dragoptdepth =  rombint2(ddragoptdepth_dz, 0.d0, z, 1d-5, 20, 100)
    
  end function dragoptdepth
  !////////////////////////////////////////////////////////////////////////
  subroutine find_z(func,zout)  
    !find redshift at which (photon/drag) optical depth = 1
    real(dl), external :: func
    real(dl), intent(out) :: zout
    real(dl) :: try1,try2,diff,avg
    integer :: i
    
    try1 = 0.d0
    try2 = 10000.d0
    
    i=0
    diff = 10.d0
    do while (diff .gt. 1d-3)
       i=i+1
       if (i .eq. 100) stop 'optical depth redshift finder did not converge'
       
       diff = func(try2)-func(try1)
       avg = 0.5d0*(try2+try1)
       if (func(avg) .gt. 1.d0) then
          try2 = avg
       else
          try1 = avg
       end if
    end do
    
    zout = avg
    
  end subroutine find_z
  !////////////////////////////////////////////////////////////////////////
  !!!!!!!!!!!!!!!!!!! end JH
end module ThermoData
!====================================================================================
