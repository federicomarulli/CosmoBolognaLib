module MGCAMB
    use precision

    ! new model selection flags
    integer :: MG_flag
    integer :: pure_MG_flag
    integer :: alt_MG_flag
    integer :: QSA_flag
    integer :: mugamma_par
    integer :: muSigma_par
    integer :: QR_par

    ! DE model flag
    integer :: DE_model

    real(dl) :: GRtrans                     !< scale factor at which MG is switched on

    ! BZ parametrization (and QS f(R))
    real(dl) :: B1
    real(dl) :: B2
    real(dl) :: lambda1_2
    real(dl) :: lambda2_2
    real(dl) :: ss

    ! Planck Parametrization
    real(dl) :: E11
    real(dl) :: E22

    ! Q-R parametrization 1
    real(dl) :: MGQfix
    real(dl) :: MGRfix

    ! Q-R parametrization 2
    real(dl) :: Qnot
    real(dl) :: Rnot
    real(dl) :: sss

    ! Growth rate gamma
    real(dl) :: Linder_gamma

    ! Symmetron
    real(dl) :: beta_star
    real(dl) :: a_star
    real(dl) :: xi_star

    ! Dilaton
    real(dl) :: beta0
    real(dl) :: xi0
    real(dl) :: DilR
    real(dl) :: DilS

    ! Hu-Sawicki f(R) gravity
    real(dl) :: F_R0
    real(dl) :: FRn

    ! DES parametrization
    real(dl) :: mu0
    real(dl) :: sigma0


    ! effective Newton's constant
    real(dl) :: ga
    real(dl) :: nn

    ! DE model parameters
    real(dl) :: w0DE              !< w0 parameters for DE
    real(dl) :: waDE              !< waDE parameters for DE

    character(len=(10)) :: MGCAMB_version = 'v 3.0'


    ! define the type MGCAMB_par_cache
    type :: MGCAMB_parameter_cache
        real(dl) :: omegab
        real(dl) :: omegac
        real(dl) :: omegav
        real(dl) :: h0
        real(dl) :: h0_Mpc
        character(len=30) :: output_root
    end type MGCAMB_parameter_cache

    type(MGCAMB_parameter_cache) :: mgcamb_par_cache

    ! define the tyoe MGCAMB_timestep_cache
    type :: MGCAMB_timestep_cache

        ! 1. Background quantities
        real(dl) :: adotoa
        real(dl) :: Hdot
        real(dl) :: grho
        real(dl) :: gpres
        real(dl) :: grhob_t
        real(dl) :: grhoc_t
        real(dl) :: grhog_t
        real(dl) :: grhor_t
        real(dl) :: grhov_t
        real(dl) :: gpresv_t
        real(dl) :: grhonu_t
        real(dl) :: gpresnu_t

        ! 2. Perturbation quantities
        real(dl) :: k
        real(dl) :: k2
        real(dl) :: dgrho
        real(dl) :: dgq
        real(dl) :: pidot_sum
        real(dl) :: dgpi_w_sum
        real(dl) :: dgpi
        real(dl) :: dgpi_diff
        real(dl) :: dgpidot
        real(dl) :: rhoDelta
        real(dl) :: rhoDeltadot

        ! 3. MG functions
        real(dl) :: mu
        real(dl) :: mudot
        real(dl) :: gamma
        real(dl) :: gammadot
        real(dl) :: q
        real(dl) :: qdot
        real(dl) :: r
        real(dl) :: rdot

        !> 4. Perturbations evolution variables
        real(dl) :: z
        real(dl) :: sigma
        real(dl) :: sigmadot
        real(dl) :: etak
        real(dl) :: etadot

        !> 5. ISW and lensing realted quantities
        real(dl) :: MG_alpha
        real(dl) :: MG_alphadot
        real(dl) :: MG_phi
        real(dl) :: MG_phidot
        real(dl) :: MG_psi
        real(dl) :: MG_psidot
        real(dl) :: MG_ISW
        real(dl) :: MG_lensing
        real(dl) :: source1
        real(dl) :: source3

    end type MGCAMB_timestep_cache

#ifdef DEBUG
    logical , parameter :: DebugMGCAMB = .true.              !< MGCAMB debug flag.This will turn on printing of many things to aid debugging the code.
#else
    logical , parameter :: DebugMGCAMB = .false.             !< MGCAMB debug flag.This will turn on printing of many things to aid debugging the code.
#endif


contains

    !---------------------------------------------------------------------------
    !> this subroutine computes the MG functions at a time-step
    subroutine MGCAMB_compute_MG_functions( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        ! Divide the cases here
        if (( MG_flag == 1 .and. pure_MG_flag /= 3 ) & ! all the mu, gamma parametrizations
            .or. MG_flag == 2 &
            .or. MG_flag == 3 ) then

            mg_cache%mu         = MGCAMB_Mu( a, mg_par_cache, mg_cache )
            mg_cache%mudot      = MGCAMB_MuDot( a, mg_par_cache, mg_cache )
            mg_cache%gamma      = MGCAMB_Gamma( a, mg_par_cache, mg_cache )
            mg_cache%gammadot   = MGCAMB_GammaDot( a, mg_par_cache, mg_cache )

            !write(*,*) 'a, k, mu, mudot, gamma, gammadot', a, mg_cache%k, mg_cache%mu,&
            !            mg_cache%mudot, mg_cache%gamma, mg_cache%gammadot

            ! other EFT functions are zero
            mg_cache%q      = 0._dl
            mg_cache%qdot   = 0._dl
            mg_cache%r      = 0._dl
            mg_cache%rdot   = 0._dl

        else if (  MG_flag == 1 .and. pure_MG_flag == 3  ) then ! the Q,R parametrization

            mg_cache%q      = MGCAMB_Q( a, mg_par_cache, mg_cache )
            mg_cache%qdot   = MGCAMB_Qdot( a, mg_par_cache, mg_cache )
            mg_cache%r      = MGCAMB_R( a, mg_par_cache, mg_cache )
            mg_cache%rdot   = MGCAMB_Rdot( a, mg_par_cache, mg_cache )

            ! other MG functions are zero
            mg_cache%mu         = 0._dl
            mg_cache%mudot      = 0._dl
            mg_cache%gamma      = 0._dl
            mg_cache%gammadot   = 0._dl

        end if

    end subroutine MGCAMB_compute_MG_functions


    !---------------------------------------------------------------------------
    !> this subroutine computes the shear sigma in MG
    subroutine MGCAMB_compute_sigma( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        if (( MG_flag == 1 .and. pure_MG_flag /= 3 ) & ! all the mu, gamma parametrizations
            .or. MG_flag == 2 &
            .or. MG_flag == 3 ) then

            ! first calculate MG_alpha
            mg_cache%MG_alpha = ( mg_cache%etak/mg_cache%k + mg_cache%mu * ( mg_cache%gamma*mg_cache%rhoDelta+ &
                                ( mg_cache%gamma -1._dl )*2._dl* mg_cache%dgpi)/(2._dl*mg_cache%k2)) / mg_cache%adotoa

            ! then calculate sigma
            mg_cache%sigma = mg_cache%k * mg_cache%MG_alpha

        else if (  MG_flag == 1 .and. pure_MG_flag == 3  ) then

            mg_cache%MG_phi      = - mg_cache%rhoDelta * mg_cache%q/(2._dl*mg_cache%k2)
            mg_cache%sigma       = (mg_cache%etak - mg_cache%k * mg_cache%MG_phi)/mg_cache%adotoa
            mg_cache%MG_alpha    = mg_cache%sigma/mg_cache%k

        end if

    end subroutine MGCAMB_compute_sigma

    !---------------------------------------------------------------------------
    !> this subroutine computes the perturbation Z in MG
    subroutine MGCAMB_compute_z( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        !> other parameters
        real(dl) :: fmu
        real(dl) :: f1
        real(dl) :: fQ
        real(dl) :: term1
        real(dl) :: term2
        real(dl) :: term3
        real(dl) :: term4
        real(dl) :: term5
        real(dl) :: term6
        real(dl) :: k2alpha

        if (( MG_flag == 1 .and. pure_MG_flag /= 3 ) & ! all the mu, gamma parametrizations
            .or. MG_flag == 2 &
            .or. MG_flag == 3 ) then

            !> adding the massive neutrinos contibutions
            fmu = mg_cache%k2+0.5d0*mg_cache%gamma*mg_cache%mu*(3._dl*(mg_cache%grhoc_t+mg_cache%grhob_t) &
                & + 4._dl*(mg_cache%grhog_t+mg_cache%grhor_t) +3._dl * (mg_cache%grhonu_t + mg_cache%gpresnu_t ))

            !> adding massive neutrinos contributions

            f1 = mg_cache%k2+3._dl*( mg_cache%adotoa**2 - mg_cache%Hdot )
            !f1 = mg_cache%k2+0.5d0*(3._dl*(mg_cache%grhoc_t+mg_cache%grhob_t) &
            !    & + 4._dl*(mg_cache%grhog_t+mg_cache%grhor_t) + 3._dl*(mg_cache%grhonu_t+mg_cache%gpresnu_t) &
            !    & + 3._dl*(mg_cache%grhov_t+mg_cache%gpresv_t))

            term1 = mg_cache%gamma*mg_cache%mu* f1 * mg_cache%dgq/mg_cache%k

            !> adding massive neutrinos contribution, if w_DE /= -1 this has to be changed
            !term2 = mg_cache%k2*mg_cache%MG_alpha* ((mg_cache%mu* mg_cache%gamma- 1._dl)*(mg_cache%grhoc_t+mg_cache%grhob_t&
            !        & +(4._dl/3._dl)*(mg_cache%grhog_t+mg_cache%grhor_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t)) &
            !        & - (mg_cache%grhov_t+ mg_cache%gpresv_t))

            term2 = mg_cache%k2*mg_cache%MG_alpha* (mg_cache%mu* mg_cache%gamma*( mg_cache%grhoc_t+mg_cache%grhob_t   &
                    & +(4._dl/3._dl)*(mg_cache%grhog_t+mg_cache%grhor_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t) ) &
                    & - 2._dl*(mg_cache%adotoa**2 - mg_cache%Hdot))

            term3= (mg_cache%mu * ( mg_cache%gamma -1._dl)* mg_cache%adotoa - mg_cache%gamma*mg_cache%mudot &
                    & - mg_cache%gammadot*mg_cache%mu )*mg_cache%rhoDelta

            ! typo corrected here
            term4 = 2._dl*mg_cache%mu*(mg_cache%gamma - 1._dl)*mg_cache%adotoa*mg_cache%dgpi_w_sum

            ! separated fromt the previous term
            term5 = -2._dl*((mg_cache%gamma-1._dl)*mg_cache%mudot -mg_cache%gammadot*mg_cache%mu)*mg_cache%dgpi

            !> adding massive neutrinos contribution
            term6= 2._dl * mg_cache%mu*(1._dl - mg_cache%gamma)* mg_cache%pidot_sum

            !> calculate etadot
            mg_cache%etadot = (term1 + term2 + term3 + term4 + term5 + term6)/( 2._dl * fmu)

            !> finally calculate Z
            mg_cache%z = mg_cache%sigma - 3._dl * mg_cache%etadot/mg_cache%k

            !> Calculate the Newtonian potential
            mg_cache%MG_psi = - mg_cache%mu * ( mg_cache%rhoDelta + 2._dl* mg_cache%dgpi)/(2._dl*mg_cache%k2)

            !> calculate the curvature perturbation potential
            mg_cache%MG_phi = mg_cache%gamma * mg_cache%MG_psi + mg_cache%mu* 1._dl*mg_cache%dgpi/mg_cache%k2

            mg_cache%MG_phidot = mg_cache%etadot - mg_cache%adotoa * (mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha) &
                                & - mg_cache%Hdot * mg_cache%MG_alpha

        else if (  MG_flag == 1 .and. pure_MG_flag == 3  ) then

            ! adding massive neutrinos contributions
            fQ = mg_cache%k2 + 0.5d0*mg_cache%q * (3._dl*(mg_cache%grhob_t+mg_cache%grhoc_t)+&
                & 4._dl*(mg_cache%grhor_t+mg_cache%grhog_t)+3._dl*(mg_cache%grhonu_t + mg_cache%gpresnu_t))

            ! fixed for w_DE /= -1
            !f1=mg_cache%k2+3._dl*( mg_cache%adotoa**2 - mg_cache%Hdot )
            f1 = mg_cache%k2+0.5d0*(3._dl*(mg_cache%grhoc_t+mg_cache%grhob_t) &
                & + 4._dl*(mg_cache%grhog_t+mg_cache%grhor_t) + 3._dl*(mg_cache%grhonu_t+mg_cache%gpresnu_t) &
                & + 3._dl*(mg_cache%grhov_t+mg_cache%gpresv_t))

            k2alpha= mg_cache%k * mg_cache%sigma

            term1 = mg_cache%q * f1 * mg_cache%dgq/mg_cache%k

            term2 = k2alpha * ((mg_cache%q - 1._dl) * ( mg_cache%grhob_t+mg_cache%grhoc_t+(4._dl/3._dl) &
                    & *(mg_cache%grhor_t+mg_cache%grhog_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t) &
                    & ) -mg_cache%grhov_t - mg_cache%gpresv_t)

            !term2 = k2alpha * ((mg_cache%q) * ( mg_cache%grhob_t+mg_cache%grhoc_t+(4._dl/3._dl) &
            !    & *(mg_cache%grhor_t+mg_cache%grhog_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t)) &
            !    & - 2._dl *(mg_cache%adotoa**2 - mg_cache%Hdot))

            term3 = -( mg_cache%qdot + (mg_cache%r-1._dl) * mg_cache%q * mg_cache%adotoa ) * mg_cache%rhoDelta

            mg_cache%etadot = (term1 + term2 + term3)/( 2._dl * fQ )

            mg_cache%z = mg_cache%sigma - 3._dl * mg_cache%etadot/mg_cache%k

            !calculating also ISW related quantities
            mg_cache%MG_psi     = mg_cache%r * mg_cache%MG_phi - mg_cache%q * 1._dl * mg_cache%dgpi/mg_cache%k2
            mg_cache%MG_phidot  = mg_cache%etadot - mg_cache%adotoa * (mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha) &
                                & - mg_cache%Hdot * mg_cache%MG_alpha

        end if

        ! calculate sigmadot
        mg_cache%sigmadot = mg_cache%k * (mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha)

    end subroutine MGCAMB_compute_z

    !---------------------------------------------------------------------------
    !> this subroutine computes the ISW term in MG
    subroutine MGCAMB_compute_ISW( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        !local variables
        real(dl) :: term0

        term0 = mg_cache%k2 + 3._dl* (mg_cache%adotoa**2._dl - mg_cache%Hdot)

        !adding MG_rhoDeltadot
        mg_cache%rhoDeltadot = -term0 * mg_cache%dgq/mg_cache%k - (mg_cache%grho + mg_cache%gpres)* mg_cache%k*mg_cache%z &
                            & - mg_cache%adotoa * mg_cache%rhoDelta - 2._dl * mg_cache%adotoa * mg_cache%dgpi

        !adding dgpidot
        mg_cache%dgpidot = mg_cache%pidot_sum - (2._dl*mg_cache%dgpi+ mg_cache%dgpi_diff )*mg_cache%adotoa

        if (( MG_flag == 1 .and. pure_MG_flag /= 3 ) & ! all the mu, gamma parametrizations
            .or. MG_flag == 2 &
            .or. MG_flag == 3 ) then

            mg_cache%MG_psidot = - 0.5d0*mg_cache%mu/mg_cache%k2*(mg_cache%rhoDeltadot+2._dl*mg_cache%dgpidot) &
                                & - 0.5d0*mg_cache%mudot/mg_cache%k2*(mg_cache%rhoDelta+2._dl*mg_cache%dgpi)

        else if (  MG_flag == 1 .and. pure_MG_flag == 3  ) then

            mg_cache%MG_psidot = mg_cache%R * mg_cache%MG_phidot + mg_cache%Rdot * mg_cache%MG_phi - &
                            & mg_cache%Qdot*mg_cache%dgpi/mg_cache%k2 - mg_cache%Q * mg_cache%dgpidot /mg_cache%k2

        end if

        mg_cache%MG_ISW = mg_cache%MG_phidot+mg_cache%MG_psidot

        mg_cache%MG_alphadot = mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha

    end subroutine MGCAMB_compute_ISW

    !---------------------------------------------------------------------------
    !> this subroutine computes the lensing term in MG
    subroutine MGCAMB_compute_lensing( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        mg_cache%MG_lensing = mg_cache%MG_phi + mg_cache%MG_psi

    end subroutine MGCAMB_compute_lensing


    !-----------------------------------------------
    !> mu(a,k) function
    function MGCAMB_Mu( a, mg_par_cache, mg_cache )
        !use ModelParams
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Mu                                       !< MG mu function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        ! local variables
        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: LKA2 ! \lambda_1^2 k^2 a^s
        real(dl) :: t1, t2, t1dot, t2dot
        real(dl) :: omm, ommdot

        real(dl) :: omegaDE_t

        ! beta, m parametrization
        real(dl) :: beta, m

        !> pure MG models
        if ( MG_flag == 1 .and. pure_MG_flag /= 3 ) then

            if ( pure_MG_flag == 1 ) then ! mu-gamma
                if ( mugamma_par == 1 ) then ! BZ parametrization
                    LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                    LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                    MGCAMB_Mu = (1._dl + B1 * LKA1)/(1._dl + LKA1)

                else if ( mugamma_par == 2 ) then ! Planck parametrization

                    ! changing the following
                    !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2

                    omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                    MGCAMB_Mu = 1._dl + E11*omegaDE_t

                else if ( mugamma_par == 3 ) then ! effective Newton constant
                    MGCAMB_Mu = 1._dl+ga*(1._dl)**nn - ga*(1._dl)**(2._dl*nn)

                else if ( mugamma_par == 4 ) then
                    MGCAMB_Mu = 1._dl

                end if

            else if ( pure_MG_flag == 2 ) then ! mu-Sigma

                if ( muSigma_par == 1 ) then ! DES parametrization

                    !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    !MGCAMB_Mu = 1._dl + mu0 * omegaDE_t/mg_par_cache%omegav

                    ! this is being changed
                    omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                    MGCAMB_Mu = 1._dl + mu0 * omegaDE_t/mg_par_cache%omegav


                else if ( muSigma_par == 2 ) then
                    MGCAMB_Mu = 1._dl

                end if

            end if

        !> alternative MG
        else if ( MG_flag == 2 ) then

            if (alt_MG_flag == 1) then !(Linder Gamma)
                omm=(mg_par_cache%omegab+mg_par_cache%omegac)/((mg_par_cache%omegab+mg_par_cache%omegac) &
                & + (1-mg_par_cache%omegab-mg_par_cache%omegac)*a**3)
                ommdot=-3._dl*omm**2*a**3*mg_cache%adotoa*(1-mg_par_cache%omegab-mg_par_cache%omegac) &
                & /(mg_par_cache%omegab+mg_par_cache%omegac)

                MGCAMB_Mu=2._dl/3._dl*omm**(Linder_gamma-1._dl)*&
                (omm**Linder_gamma+2-3._dl*Linder_gamma+3._dl*(Linder_gamma-0.5d0)*omm)

            else if ( alt_MG_flag == 2 ) then
                MGCAMB_Mu = 1._dl
            end if


        !> QSA models
        else if ( MG_flag == 3 ) then

            if ( QSA_flag == 1 ) then ! f(R)
                LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                LKA2 = lambda2_2 * mg_cache%k2 * a**ss
                MGCAMB_Mu = (1._dl + B1 * LKA1)/(1._dl + LKA1)
                MGCAMB_Mu = MGCAMB_Mu/(1._dl - 1.4d-8 * lambda1_2 * a**3)

            else if ( QSA_flag == 2 .or. &  ! beta, m parametrization
                      QSA_flag == 3 .or. &
                      QSA_flag == 4 ) then
                beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                m       = MGCAMB_M( a, mg_par_cache, mg_cache )
                t1      = (2._dl*beta**2._dl)*mg_cache%k2
                t2      = (m**2._dl)*a**2._dl

                MGCAMB_Mu = (mg_cache%k2 + t1 + t2)/(mg_cache%k2 + t2)
                


            else if ( QSA_flag == 5 )  then
                MGCAMB_Mu = 1._dl

            end if

        end if

    end function MGCAMB_Mu

    !-----------------------------------------------
    !> \dot{mu}(a,k) function
    function MGCAMB_Mudot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Mudot                                    !< MG mudot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        ! local variables
        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: LKA2 ! \lambda_1^2 k^2 a^s
        real(dl) :: t1,t2,t1dot,t2dot
        real(dl) :: omm, ommdot

        ! mapping beta,m into mu,gamma
        real(dl) :: beta, betadot, m, mdot
        real(dl) :: mu

        real(dl) :: omegaDEdot

        !> pure MG models
        if ( MG_flag == 1 .and. pure_MG_flag /= 3 ) then

            if ( pure_MG_flag == 1 ) then ! mu-gamma
                if ( mugamma_par == 1 ) then ! BZ parametrization
                    LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                    LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                    MGCAMB_Mudot = ((B1 - 1._dl) * mg_cache%adotoa * ss * LKA1) / ((1._dl+LKA1)**2._dl)

                else if ( mugamma_par == 2 ) then ! Planck parametrization

                    ! changingh the following quantity
                    !omegaDEdot = - 3._dl * mg_cache%adotoa * (mg_cache%grhov_t + mg_cache%gpresv_t) &
                    !            & / a**2 / 3._dl / mg_par_cache%h0_Mpc**2

                    omegaDEdot=-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                            & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t

                    MGCAMB_Mudot = E11*omegaDEdot

                else if ( mugamma_par == 3 ) then
                    MGCAMB_Mudot = mg_cache%adotoa*a*ga*nn*(-1._dl+2._dl*(1._dl-a)**nn)*(1._dl-a)**(nn-1._dl)

                else if ( mugamma_par == 4 ) then
                    MGCAMB_Mudot = 0._dl

                end if

            else if ( pure_MG_flag == 2 ) then ! mu-Sigma

                if ( muSigma_par == 1 ) then ! DES parametrization
                    ! changing the following
                    !omegaDEdot = - 3._dl * mg_cache%adotoa * (mg_cache%grhov_t + mg_cache%gpresv_t) &
                    !            & / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    omegaDEdot=-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                                & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t

                    MGCAMB_Mudot =  mu0 * omegaDEdot/mg_par_cache%omegav

                else if ( muSigma_par == 2 ) then
                    MGCAMB_Mudot = 0._dl

                end if

            end if

        !> alternative MG
        else if ( MG_flag == 2 ) then

            if (alt_MG_flag == 1) then !(Linder Gamma)
                mu = MGCAMB_Mu( a, mg_par_cache, mg_cache )

                omm=(mg_par_cache%omegab+mg_par_cache%omegac)/((mg_par_cache%omegab+mg_par_cache%omegac) &
                    & +(1-mg_par_cache%omegab-mg_par_cache%omegac)*a**3)
                ommdot=-3._dl*omm**2*a**3*mg_cache%adotoa*(1-mg_par_cache%omegab-mg_par_cache%omegac) &
                    & /(mg_par_cache%omegab+mg_par_cache%omegac)

                MGCAMB_Mudot = mu/omm*(Linder_gamma-1._dl)*ommdot+&
                    2._dl/3._dl*omm**(Linder_gamma-1._dl)*ommdot*&
                    (Linder_gamma*omm**(Linder_gamma-1._dl)+3._dl*(Linder_gamma-0.5d0))

            else if ( alt_MG_flag == 2 ) then
                MGCAMB_Mudot = 0._dl
            end if


        !> QSA models
        else if ( MG_flag == 3 ) then

            if ( QSA_flag == 1 ) then ! f(R)
                LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                MGCAMB_Mudot = ((B1 - 1._dl) * mg_cache%adotoa * ss * LKA1) / ((1._dl+LKA1)**2._dl)
                mu = MGCAMB_Mu( a, mg_par_cache, mg_cache )
                MGCAMB_Mudot = MGCAMB_Mudot/(1._dl - 1.4d-8 * lambda1_2 * a**3) + 3._dl * &
                                mu* mg_cache%adotoa *a**3 *(1.4d-8 * lambda1_2 ) &
                                /(1._dl - 1.4d-8 * lambda1_2 * a**3)

            else if ( QSA_flag == 2 .or. &  ! beta, m parametrization
                    QSA_flag == 3 .or. &
                    QSA_flag == 4 ) then

                beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                m       = MGCAMB_M( a, mg_par_cache, mg_cache )
                betadot = MGCAMB_Betadot( a, mg_par_cache, mg_cache )
                mdot    = MGCAMB_Mdot( a, mg_par_cache, mg_cache )

                t1 = (2._dl*beta**2._dl)*mg_cache%k2
                t2 = (m**2._dl)*a**2._dl
                t1dot = 4._dl*beta*betadot*mg_cache%k2
                t2dot = (2._dl*a**2._dl)*(m*mdot+ (m**2._dl) *mg_cache%adotoa)

                MGCAMB_Mudot = (t1dot*(mg_cache%k2 + t2) - t1*t2dot)/((mg_cache%k2 + t2)**2._dl)


            else if ( QSA_flag == 5 )  then
                MGCAMB_Mudot = 0._dl

            end if

        end if

    end function MGCAMB_Mudot

    !-----------------------------------------------
    ! gamma(a,k) function
    function MGCAMB_Gamma( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Gamma                                    !< MG gamma function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: LKA2 ! \lambda_1^2 k^2 a^s
        real(dl) :: t1,t2, t1dot, t2dot

        real(dl) :: beta, m
        real(dl) :: omegaDE_t

        real(dl) :: sigma_t
        real(dl) :: mu_t

        !> pure MG models
        if ( MG_flag == 1 .and. pure_MG_flag /= 3 ) then

            if ( pure_MG_flag == 1 ) then ! mu-gamma
                if ( mugamma_par == 1 ) then ! BZ parametrization
                    LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                    LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                    MGCAMB_Gamma = (1._dl + B2 * LKA2)/(1._dl +LKA2)

                else if ( mugamma_par == 2 ) then ! Planck parametrization
                    ! changing the following
                    !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                    MGCAMB_Gamma = 1._dl+E22*omegaDE_t

                else if ( mugamma_par == 3 ) then
                    MGCAMB_Gamma = 1._dl

                else if ( mugamma_par == 4 ) then
                    MGCAMB_Gamma = 1._dl
                end if

            else if ( pure_MG_flag == 2 ) then ! mu-Sigma

                if ( muSigma_par == 1 ) then ! DES parametrization
                    ! changing the following
                    !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                    sigma_t = 1._dl + sigma0 * omegaDE_t / mg_par_cache%omegav
                    mu_t    = 1._dl + mu0 * omegaDE_t / mg_par_cache%omegav
                    MGCAMB_Gamma = 2._dl * sigma_t / mu_t - 1._dl

                else if ( muSigma_par == 2 ) then
                    MGCAMB_Gamma = 1._dl

                end if

            end if

        !> alternative MG
        else if ( MG_flag == 2 ) then

            if (alt_MG_flag == 1) then !(Linder Gamma)
                MGCAMB_Gamma = 1._dl

            else if ( alt_MG_flag == 2 ) then
                MGCAMB_Gamma = 1._dl
            end if


        !> QSA models
        else if ( MG_flag == 3 ) then

            if ( QSA_flag == 1 ) then ! f(R)
                LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                MGCAMB_Gamma = (1._dl + B2 * LKA2)/(1._dl +LKA2)

            else if ( QSA_flag == 2 .or. &  ! beta, m parametrization
                    QSA_flag == 3 .or. &
                    QSA_flag == 4 ) then

                beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                m       = MGCAMB_M( a, mg_par_cache, mg_cache )

                t1 = (2._dl*beta**2._dl)*mg_cache%k2
                t2 = (m**2._dl)*a**2._dl

                MGCAMB_Gamma = (mg_cache%k2 - t1 + t2)/(mg_cache%k2 + t1 + t2)


            else if ( QSA_flag == 5 )  then
                MGCAMB_Gamma = 1._dl

            end if

        end if


    end function MGCAMB_Gamma


    !-----------------------------------------------
    ! \dot{gamma}(a,k) function
    function MGCAMB_Gammadot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Gammadot                                 !< MG gammadot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters


        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: LKA2 ! \lambda_1^2 k^2 a^s
        real(dl) :: t1,t2,t1dot,t2dot

        real(dl) :: beta, betadot, m, mdot
        real(dl) :: omegaDE_t, omegaDEdot
        real(dl) :: sigma_t, sigmadot_t
        real(dl) :: mu_t, mudot_t



        !> pure MG models
        if ( MG_flag == 1 .and. pure_MG_flag /= 3 ) then

            if ( pure_MG_flag == 1 ) then ! mu-gamma
                if ( mugamma_par == 1 ) then ! BZ parametrization
                    LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                    LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                    MGCAMB_Gammadot = ((B2 -1._dl)*mg_cache%adotoa * ss* LKA2)/((1._dl+LKA2)**2._dl)

                else if ( mugamma_par == 2 ) then ! Planck parametrization
                    ! changing the following
                    !omegaDEdot = - 3._dl * mg_cache%adotoa * (mg_cache%grhov_t + mg_cache%gpresv_t) &
                    !            & / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    omegaDEdot=-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                                & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t

                    MGCAMB_Gammadot = E22*omegaDEdot

                else if ( mugamma_par == 3 ) then
                    MGCAMB_Gammadot = 0._dl

                else if ( mugamma_par == 4 ) then
                    MGCAMB_Gammadot = 0._dl

                end if

            else if ( pure_MG_flag == 2 ) then ! mu-Sigma

                if ( muSigma_par == 1 ) then ! DES parametrization

                ! changing the following
                omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                omegaDEdot  =-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                                & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t
                sigma_t     = 1._dl + sigma0 * omegaDE_t / mg_par_cache%omegav
                sigmadot_t  = sigma0 * omegaDEdot / mg_par_cache%omegav
                mu_t        = 1._dl + mu0 * omegaDE_t / mg_par_cache%omegav
                mudot_t     = mu0 * omegaDEdot / mg_par_cache%omegav
                MGCAMB_Gammadot = 2._dl * sigmadot_t / mu_t - 2._dl *sigma_t*mudot_t/mu_t**2

                else if ( muSigma_par == 2 ) then
                    MGCAMB_Gammadot = 0._dl

            end if

            end if

        !> alternative MG
        else if ( MG_flag == 2 ) then

            if (alt_MG_flag == 1) then !(Linder Gamma)
                MGCAMB_Gammadot = 0._dl

            else if ( alt_MG_flag == 2 ) then
                MGCAMB_Gammadot = 0._dl
            end if


        !> QSA models
        else if ( MG_flag == 3 ) then

            if ( QSA_flag == 1 ) then ! f(R)
                LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                MGCAMB_Gammadot = ((B2 -1._dl)*mg_cache%adotoa * ss* LKA2)/((1._dl+LKA2)**2._dl)

            else if ( QSA_flag == 2 .or. &  ! beta, m parametrization
                    QSA_flag == 3  .or. &
                    QSA_flag == 4 ) then
                beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                m       = MGCAMB_M( a, mg_par_cache, mg_cache )
                betadot = MGCAMB_Betadot( a, mg_par_cache, mg_cache )
                mdot    = MGCAMB_Mdot( a, mg_par_cache, mg_cache )

                t1      = (2._dl*beta**2._dl)*mg_cache%k2
                t2      = (m**2._dl)*a**2._dl
                t1dot   = 4._dl*beta*betadot*mg_cache%k2
                t2dot   = (2._dl*a**2._dl)*(m*mdot + (m**2._dl) *mg_cache%adotoa)

                MGCAMB_Gammadot = 2._dl*(t1*t2dot-t1dot*(mg_cache%k2 + t2))/((mg_cache%k2 + t1 + t2)**2._dl)

            else if ( QSA_flag == 5 )  then
                MGCAMB_Gammadot = 0._dl

            end if

        end if

    end function MGCAMB_Gammadot


!----------------------------------------------------------------------------------------------
!> MGCAMB (beta, m) parametrization, QSA for scalar-tensor models

    !-----------------------------------------------
    !> m(a) function
    function MGCAMB_M( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_M                                        !< MG m function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        real(dl) :: FRm0

        ! SYMMETRON
        if( QSA_flag ==  2 ) then
            MGCAMB_M = (mg_par_cache%H0/3.0D05) / (xi_star) * sqrt(1._dl-(a_star/a)**3._dl)

        ! DILATON: based on 1206.3568
        else if ( QSA_flag ==  3 ) then
            MGCAMB_M = (mg_par_cache%H0/3.0D05) /(xi0) * a**(- DilR)

        ! Hu-Sawicki f(R) model: m, beta parametrization as in 1305.5647
        else if ( QSA_flag ==  4 )then
            FRm0 = (mg_par_cache%h0/3.0D05)*sqrt((4._dl*mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac) &
                    & /((FRn+1._dl)*F_R0))!note factor of c here
            MGCAMB_M = FRm0 * ((4._dl * mg_par_cache%omegav + (mg_par_cache%omegab + mg_par_cache%omegac)*a**(-3._dl)) &
                    & /(4._dl * mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac))**(FRn/2._dl+1._dl)

        end if

    end function MGCAMB_M


    !-----------------------------------------------
    !> \dot{m}(a) function
    function MGCAMB_Mdot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Mdot                                     !< MG mdot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        real(dl) :: FRm0
        real(dl) :: m

        m = MGCAMB_M( a, mg_par_cache, mg_cache )

        ! SYMMETRON
        if( QSA_flag ==  2 ) then
            MGCAMB_Mdot = 1.5d0*(mg_par_cache%H0/3.0D05)/(xi_star)*((a_star/a)**3._dl*mg_cache%adotoa)/&
                        & (sqrt(1._dl-(a_star/a)**3._dl))

        ! DILATON
        else if ( QSA_flag ==  3 ) then
            MGCAMB_Mdot = - DilR * m * mg_cache%adotoa

        ! Hu-Sawicki f(R) model
        else if ( QSA_flag ==  4 )then

            FRm0 = (mg_par_cache%h0/3.0D05)*sqrt((4._dl*mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac)/ &
                    & ((FRn+1._dl)*F_R0))
            MGCAMB_Mdot = m / (4._dl * mg_par_cache%omegav + (mg_par_cache%omegab + mg_par_cache%omegac)*a**(-3._dl)) &
                    & * (-3._dl*FRn/2._dl-3._dl)*((mg_par_cache%omegab + mg_par_cache%omegac)* a**(-3._dl)*mg_cache%adotoa)!/(4._dl * mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac)) ! complete this

        end if

    end function MGCAMB_Mdot

    !-----------------------------------------------
    !> beta(a) function
    function MGCAMB_Beta( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Beta                                     !< MG beta function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        ! SYMMETRON
        if( QSA_flag == 2 ) then
            MGCAMB_Beta =  beta_star * sqrt(1._dl-(a_star/a)**3._dl)

        ! DILATON
        else if ( QSA_flag == 3 ) then
            MGCAMB_Beta = beta0 * exp((DilS)/(2._dl* DilR - 3._dl)*(a**(2._dl* DilR - 3._dl)-1._dl))

        ! Hu-Sawicki f(R) model
        else if ( QSA_flag == 4 )then
            MGCAMB_Beta = beta0

        end if

    end function MGCAMB_Beta

    !-----------------------------------------------
    !> \dot{beta}(a) function
    function MGCAMB_Betadot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Betadot                                  !< MG betadot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        real(dl) :: beta

        beta = MGCAMB_Beta( a, mg_par_cache, mg_cache )

        ! SYMMETRON
        if( QSA_flag == 2 ) then
            MGCAMB_Betadot = 1.5d0 * (beta_star * (a_star/a)**3._dl * mg_cache%adotoa) /( sqrt(1._dl-(a_star/a)**3._dl))

        ! DILATON
        else if ( QSA_flag == 3 ) then
            MGCAMB_Betadot = beta * (DilS * a**(2._dl* DilR - 3._dl) *  mg_cache%adotoa)

        ! Hu-Sawicki f(R) model
        else if ( QSA_flag == 4 )then
            MGCAMB_Betadot = 0._dl


        end if

    end function MGCAMB_Betadot

    !----------------------------------------------------------------------------------------------
    !> MGCAMB (Q,R) parametrization, QSA for scalar-tensor models

    !-----------------------------------------------
    !> Q(a,k) function
    function MGCAMB_Q( a, mg_par_cache, mg_cache )
    implicit none
    real(dl) :: a                                               !< scale factor
    real(dl) :: MGCAMB_Q                                        !< MG Q function
    type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
    type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        if ( QR_par == 1) then
            MGCAMB_Q = MGQfix

        else if ( QR_par == 2 ) then
            MGCAMB_Q = 1._dl + (Qnot - 1._dl)* a**sss

        end if

    end function MGCAMB_Q

    !-----------------------------------------------
    !> \dot{Q}(a,k) function
    function MGCAMB_Qdot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Qdot                                     !< MG Qdot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        if ( QR_par == 1 ) then
            MGCAMB_Qdot = 0._dl

        else if ( QR_par == 2 ) then
            MGCAMB_Qdot = (Qnot - 1._dl)*mg_cache%adotoa* sss* a**(sss)
        end if

    end function MGCAMB_Qdot

    !-----------------------------------------------
    ! R(a,k) function
    function MGCAMB_R( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_R                                        !< MG R function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters


        if ( QR_par == 1 ) then
            MGCAMB_R=MGRfix

        else if ( QR_par == 2 ) then
            MGCAMB_R = 1._dl + (Rnot - 1._dl)* a**sss

        end if

    end function MGCAMB_R

    !-----------------------------------------------
    ! \dot{R}(a,k) function
    function MGCAMB_Rdot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Rdot                                     !< MG Rdot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        if ( QR_par == 1 ) then
            MGCAMB_Rdot = 0._dl

        else if ( QR_par == 2 ) then
            MGCAMB_Rdot = (Rnot - 1._dl)*mg_cache%adotoa* sss* a**(sss)

        end if

    end function MGCAMB_Rdot

    ! ---------------------------------------------------------------------------------------------
    !> Modifying the background
    subroutine MGCAMB_DarkEnergy( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache),  intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        real(dl) :: wnow

        if ( DE_model == 0 ) then
            mg_cache%grhov_t = 3._dl*mg_par_cache%h0_Mpc**2 * mg_par_cache%omegav *a**2
            mg_cache%gpresv_t = - mg_cache%grhov_t
        else if ( DE_model == 1 ) then
            mg_cache%grhov_t = 3._dl*mg_par_cache%h0_Mpc**2*mg_par_cache%omegav*a**(-1._dl-3._dl*w0DE)
            mg_cache%gpresv_t = mg_cache%grhov_t * w0DE
        else if (DE_model == 2 ) then
            wnow = w0DE+(1._dl-a)*waDE
            !mg_cache%grhov_t = 3._dl*mg_par_cache%h0_Mpc**2*mg_par_cache%omegav*a**(-1._dl-3._dl*wnow)
            mg_cache%grhov_t = 3._dl*mg_par_cache%h0_Mpc**2*mg_par_cache%omegav*a**(-1.d0-3.d0*w0DE-3.d0*waDE)*Exp(3.d0*waDE*(a-1.d0))
            mg_cache%gpresv_t = mg_cache%grhov_t * wnow
        else if ( DE_model == 3 ) then
            write(*,*) 'This will contain the reconstruction of w_DE(a)'
            write(*,*) 'Not implemented yet'
            stop
        else if ( DE_model == 4 ) then
            write(*,*) 'This will contain the reconstruction of rho_DE(a)'
            write(*,*) 'Not implemented yet'
            stop
        else
            write(*,*) 'choose a DE model'
            stop
        end if

    end subroutine MGCAMB_DarkEnergy


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the MGCAMB model parameters
    subroutine MGCAMB_read_model_params( mg_par_cache )
        use IniFile
        Type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        ! 1. MG_flag
        MG_flag = Ini_Read_Int('MG_flag', 0)

        if ( MG_flag /= 0 ) then
            call print_MGCAMB_header
            write(*,*)
            write(*,*) 'MG_flag:', MG_flag

            write(*,*) 'Debug:', DebugMGCAMB


            ! read GRtrans
            GRtrans = Ini_Read_Double('GRtrans',0.01_dl)
            write(*,*) '    GRtrans:', GRtrans

            ! 1. pure MG models
            if ( MG_flag == 1 ) then

                pure_MG_flag = Ini_Read_Int('pure_MG_flag', 1)

                if ( pure_MG_flag == 1 ) then ! mu-gamma
                    write(*,*) '    MGCAMB: mu-gamma parametrization'
                    mugamma_par = Ini_Read_Int('mugamma_par' , 1)
                    if ( mugamma_par == 1 ) then
                        write(*,*) '        BZ parametrization'
                        B1= Ini_Read_Double('B1',0._dl)
                        B2= Ini_Read_Double('B2',0._dl)
                        lambda1_2= Ini_Read_Double('lambda1_2',0._dl)
                        lambda2_2= Ini_Read_Double('lambda2_2',0._dl)
                        ss= Ini_Read_Double('ss',0._dl)
                    else if ( mugamma_par == 2 ) then
                        write(*,*) '        Planck parametrization'
                        E11     = Ini_Read_Double('E11', 0._dl)
                        E22     = Ini_Read_Double('E22', 0._dl)
                        write(*,*) 'E11, E22', E11, E22
                    else if ( mugamma_par == 3 ) then
                        write(*,*) '        Effective Newton constant'
                        ga      = Ini_Read_Double('ga', 0._dl)
                        nn      = Ini_Read_Double('nn', 0._dl)
                        write(*,*) 'ga, nn:', ga, nn
                    else
                        write(*,*) ' write your own mu-gamma parametrization in mgcamb.f90'
                        stop
                    end if


                else if ( pure_MG_flag == 2 ) then ! mu-Sigma
                    write(*,*) '    MGCAMB: mu-Sigma parametrization'
                    muSigma_par = Ini_Read_Int('musigma_par', 1)
                    if ( muSigma_par == 1 ) then
                        write(*,*) '        DES parametrization'
                        mu0     = Ini_Read_Double('mu0', 0._dl)
                        sigma0  = Ini_Read_Double('sigma0', 0._dl)
                        write(*,*) 'mu0, sigma0:', mu0, sigma0
                    else if ( muSigma_par == 2 ) then
                        write(*,*) 'write you own mu-sigma parametrization in mgcamb.f90'
                        stop
                    else
                        write(*,*) 'Please choose a model in params_MG.ini'
                        stop
                    end if

                else if ( pure_MG_flag == 3 ) then ! Q-R
                    write(*,*) '    MGCAMB: Q-R parametrization'
                    QR_par = Ini_Read_Int('QR_par', 1)
                    if ( QR_par == 1 ) then
                        MGQfix=Ini_Read_Double('MGQfix', 0._dl)
                        MGRfix=Ini_Read_Double('MGRfix', 0._dl)
                    else if ( QR_par == 2 ) then
                        Qnot=Ini_Read_Double('Qnot', 0._dl)
                        Rnot=Ini_Read_Double('Rnot', 0._dl)
                        sss=Ini_Read_Double('sss', 0._dl)
                    else if ( QR_par == 3 ) then
                        write(*,*) 'write your own QR parametrization in mgcamb.f90'
                        stop
                    else
                        write(*,*) 'Please choose a model in params_MG.ini'
                        stop
                    end if

                end if

                ! Checking DE Model
                DE_model = Ini_Read_Int('DE_model', 0)

                write(*,*) 'DE_model:', DE_model

                if ( DE_model == 1 ) then
                    w0DE = Ini_Read_Double('w0DE', -1._dl)
                else if ( DE_model == 2 ) then
                    w0DE = Ini_Read_Double('w0DE', -1._dl)
                    waDE = Ini_Read_Double('waDE', 0._dl)
                else if ( DE_model == 3 ) then
                    write(*,*) 'This will contain the reconstruction of w_DE(a)'
                    write(*,*) 'Not implemented yet'
                    stop
                else if ( DE_model == 4 ) then
                    write(*,*) 'This will contain the reconstruction of rho_DE(a)'
                    write(*,*) 'Not implemented yet'
                    stop
                else if ( DE_model /= 0 ) then
                    write(*,*) 'Please choose a DE model'
                    stop
                end if



            else if ( MG_flag == 2 ) then
                alt_MG_flag = Ini_Read_Int('alt_MG_flag', 1)
                if ( alt_MG_flag == 1 ) then
                    write(*,*) '    MGCAMB: Linder Gamma'
                    Linder_gamma = Ini_Read_Double('Linder_gamma', 0._dl)
                else if ( alt_MG_flag == 2 ) then
                    write(*,*) 'Please write your alternative MG model in mgcamb.f90'
                    stop
                else
                    write(*,*) 'Please choose a model in params_MG.ini'
                    stop
                end if

                ! Checking DE Model
                DE_model = Ini_Read_Int('DE_model', 0)

                if ( DE_model /= 0 ) then
                    write(*,*) 'alternative MG models supported only with cosmological constant!'
                end if


            else if ( MG_flag == 3 ) then
                write(*,*) '    MGCAMB: quasi-static models'
                QSA_flag = Ini_Read_Int('QSA_flag', 1)
                if ( QSA_flag ==  1 ) then
                    write(*,*) '        QSA f(R)'
                    B1 = 4._dl/3._dl
                    lambda1_2= Ini_Read_Double('B0',0._dl) ! it is considered as the B0 parameter here
                    lambda1_2 = (lambda1_2*(299792458.d-3)**2)/(2._dl*mg_par_cache%H0**2)
                    B2 = 0.5d0
                    lambda2_2 = B1* lambda1_2
                    ss = 4._dl

                else if ( QSA_flag ==  2 ) then
                    write(*,*) '        QSA Symmetron'
                    beta_star = Ini_Read_Double('beta_star', 0._dl)
                    xi_star = Ini_Read_Double ('xi_star', 0._dl)
                    a_star = Ini_Read_Double('a_star', 0._dl)
                    GRtrans = a_star

                else if ( QSA_flag ==  3 ) then
                    write(*,*) '        QSA Dilaton'
                    ! GENERALIZED DILATON
                    beta0 = Ini_Read_Double('beta0', 0._dl)
                    xi0 = Ini_Read_Double('xi0', 0._dl)
                    DilR = Ini_Read_Double('DilR', 0._dl)
                    DilS = Ini_Read_Double('DilS', 0._dl)

                else if ( QSA_flag ==  4 ) then
                    write(*,*) '        QSA Hu-Sawicki f(R)'
                    F_R0 = Ini_Read_Double('F_R0', 0._dl)
                    FRn = Ini_Read_Double('FRn', 0._dl)
                    beta0 = 1._dl/sqrt(6._dl)
                else if ( QSA_flag ==  5 ) then
                    write(*,*) 'Please write your QSA model in mgcamb.f90'
                    stop

                end if

                ! Checking DE Model
                DE_model = Ini_Read_Int('DE_model', 0)

                if ( DE_model /= 0 ) then
                    write(*,*) 'QSA models supported only with cosmological constant!'
                end if

            else
                write(*,*) ' Please choose a model'
                stop
            end if



        end if


    end subroutine MGCAMB_read_model_params

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the MGCAMB header.
    subroutine print_MGCAMB_header

        implicit none

        ! print the header:
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') "     __  _________  ________   __  ______  "
        write(*,'(a)') "    /  \/  / ____/ / ___/ _ | /  |/  / _ ) "
        write(*,'(a)') "   / /\_/ / /_,-, / /__/ __ |/ /|_/ / _  | "
        write(*,'(a)') "  /_/  /_/_____/  \___/_/ |_/_/  /_/____/  "//" "//MGCAMB_version
        write(*,'(a)') "  "
        write(*,'(a)') "        Modified Growth with CAMB "
        write(*,'(a)') "  "
        write(*,'(a)') "***************************************************************"

    end subroutine print_MGCAMB_header

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the mgcamb_cache to zero
    subroutine MGCAMB_timestep_cache_nullify( mg_cache )
        use precision
        implicit none

        type(MGCAMB_timestep_cache),  intent(inout) :: mg_cache      !< cache containing the time-dependent quantities

        ! 1. Background quantities
        mg_cache%adotoa     = 0._dl
        mg_cache%Hdot       = 0._dl
        mg_cache%grho       = 0._dl
        mg_cache%gpres      = 0._dl
        mg_cache%grhob_t    = 0._dl
        mg_cache%grhoc_t    = 0._dl
        mg_cache%grhog_t    = 0._dl
        mg_cache%grhor_t    = 0._dl
        mg_cache%grhov_t    = 0._dl
        mg_cache%gpresv_t   = 0._dl
        mg_cache%grhonu_t   = 0._dl
        mg_cache%gpresnu_t  = 0._dl

        ! 2. Perturbation quantities
        mg_cache%k          = 0._dl
        mg_cache%k2         = 0._dl
        mg_cache%dgrho      = 0._dl
        mg_cache%dgq        = 0._dl
        mg_cache%pidot_sum  = 0._dl
        mg_cache%dgpi_w_sum = 0._dl
        mg_cache%dgpi       = 0._dl
        mg_cache%dgpi_diff  = 0._dl
        mg_cache%dgpidot    = 0._dl
        mg_cache%rhoDelta   = 0._dl
        mg_cache%rhoDeltadot= 0._dl

        ! 3. MG functions
        mg_cache%mu         = 0._dl
        mg_cache%mudot      = 0._dl
        mg_cache%gamma      = 0._dl
        mg_cache%gammadot   = 0._dl
        mg_cache%q          = 0._dl
        mg_cache%qdot       = 0._dl
        mg_cache%r          = 0._dl
        mg_cache%rdot       = 0._dl

        !> 4. Perturbations evolution variables
        mg_cache%z          = 0._dl
        mg_cache%sigma      = 0._dl
        mg_cache%sigmadot   = 0._dl
        mg_cache%etak       = 0._dl
        mg_cache%etadot     = 0._dl

        !> 5. ISW and lensing realted quantities
        mg_cache%MG_alpha   = 0._dl
        mg_cache%MG_alphadot= 0._dl
        mg_cache%MG_phi     = 0._dl
        mg_cache%MG_phidot  = 0._dl
        mg_cache%MG_psi     = 0._dl
        mg_cache%MG_psidot  = 0._dl
        mg_cache%MG_ISW     = 0._dl
        mg_cache%MG_lensing = 0._dl
        mg_cache%source1    = 0._dl
        mg_cache%source3    = 0._dl

    end subroutine

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that opens the MGCAMB cache files (for Debug)
    subroutine MGCAMB_open_cache_files
        use precision
        implicit none

        ! 1. Open sources file
        open(unit=111, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_sources.dat', status="new", &
            & action="write")
        write(111,*)  'k  ', 'a  ', 'MG_ISW  ', 'MG_Lensing  ', 'S_T  ', 'S_lensing'

        ! 2 Open MG functions file
        open(unit=222, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_MG_fncs.dat', status="new",&
            & action="write")
        write(222,*)  'k  ', 'a  ', 'mu  ', 'gamma ', 'Q ', 'R ', 'Phi ', 'Psi ', 'dPhi ', 'dPsi '

        ! 3. Open Einstein solutions file
        open(unit=333, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_EinsteinSol.dat', status="new",&
            & action="write")
        write(333,*) 'k  ', 'a  ', 'etak  ', 'z  ', 'sigma  ', 'etadot  ', 'sigmadot  '

        ! 4. Open Perturbation solution file
        open(unit=444, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_PerturbSol.dat', status="new",&
        & action="write")
        write(444,*)  'k  ', 'a  ', 'dgrho  ', 'dgq  ', 'rhoDelta  ', 'dgpi  ', 'pidot_sum  ', 'dgpi_w_sum  '

        ! 5. Open Background file
        open(unit=555, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_Background.dat', status="new",&
            & action="write")
        write(555,*)  'k  ', 'a  ', 'H  ', 'Hdot  ', 'grhov_t  ', 'gpresv_t  '

    end subroutine MGCAMB_open_cache_files


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that closes the MGCAMB cache files (for Debug)
    subroutine MGCAMB_close_cache_files
        use precision
        implicit none

        close(111);close(222); close(333);close(444);close(555)

    end subroutine MGCAMB_close_cache_files


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints the MGCAMB cache on a file
    subroutine MGCAMB_dump_cache( a, mg_cache )
        use precision
        implicit none

        real(dl), intent(in) :: a   !< scale factor
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache      !< cache containing the time-dependent quantities
        character(*), parameter :: cache_output_format = 'e18.8'


        ! 1. Write the sources
        write(111,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%MG_ISW, mg_cache%MG_Lensing,&
                                                    & mg_cache%source1, mg_cache%source3

        ! 2. Write the MG functions and the potentials
        write(222,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%mu, mg_cache%gamma, mg_cache%q, mg_cache%r, &
                                                & mg_cache%MG_phi, mg_cache%MG_psi, mg_cache%MG_phidot, mg_cache%MG_psidot

        ! 3. Write the Einstein equations solutions
        write(333,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%etak, mg_cache%z, mg_cache%sigma,&
                                                & mg_cache%etadot,mg_cache%sigmadot

        ! 4. Write the Perturbations Solutions
        write(444,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%dgrho, mg_cache%dgq, mg_cache%rhoDelta,&
                                                    & mg_cache%dgpi, mg_cache%pidot_sum, mg_cache%dgpi_w_sum

        !5. Write the background
        write(555,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%adotoa, mg_cache%Hdot, mg_cache%grhov_t,&
                                                    & mg_cache%gpresv_t



    end subroutine MGCAMB_dump_cache



end module MGCAMB

