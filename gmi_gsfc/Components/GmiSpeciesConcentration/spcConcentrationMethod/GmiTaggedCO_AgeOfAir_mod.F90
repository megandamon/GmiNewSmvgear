!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiTaggedCO_AgeOfAir_mod
!
! !INTERFACE:
!
  module GmiTaggedCO_AgeOfAir_mod
!
! !USES:
  use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
  use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
  use GmiSpcConcentrationMethod_mod, only : Get_concentration, Set_concentration
  use GmiSpcConcentrationMethod_mod, only : Get_tracer_opt   , Get_efol_time &
                                    , Get_SO3daily_infile_name, Get_SO3monthly_infile_name
  use GmiTimeControl_mod       , only : GmiSplitDateTime, t_GmiClock, Get_gmiTimeStep &
                                        , Get_gmiSeconds, Get_curGmiDate
  use GmiGrid_mod              , only : t_gmiGrid, Get_numSpecies
  use GmiGrid_mod              , only : Get_i1, Get_i2, Get_ju1, Get_j2, Get_k1, Get_k2
  use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_procID
  use GmiDiagnosticsMethod_mod  , only : t_Diagnostics, Get_pr_diag
  use GmiMetFieldsControl_mod   , only : t_metFields, Get_pbl, Get_gridBoxHeight, Get_press3c &
                                         , Get_kel, Get_humidity
  use GmiTracerMethod_mod, only : Update_strato3
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: calcTaggedCO_AgeOfAir
!
# include "GmiParameters.h"
# include "gmi_time_constants.h"
# include "gmi_phys_constants.h"
# include "setkin_par.h"
!
! !DESCRIPTION:
!   total_mass      : mass of the atmosphere (kg)
!   layer1_mass     : mass of the lowest model level of the atmosphere (kg)
!
! !AUTHOR:
!
! !REVISION HISTORY:
!
!EOP
!-------------------------------------------------------------------------
  CONTAINS
!-------------------------------------------------------------------------
!BOP
  subroutine calcTaggedCO_AgeOfAir (SpeciesConcentration, total_mass, layer1_mass, &
                gmiClock, gmiDomain, Diagnostics, metFields, mw, londeg, latdeg, coscen,  &
                numSpecies, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi,  &
                i1_gl, i2_gl, ju1_gl, j2_gl)
!
!
  implicit none
!
!  real*8           , intent(in) :: total_vol_layer1
  real*8           , intent(in) :: total_mass, layer1_mass
  type(t_gmiDomain), intent(in) :: gmiDomain
  type(t_GmiClock ), intent(in) :: gmiClock
  integer          , intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
  integer          , intent(in) :: numSpecies, i1_gl, i2_gl, ju1_gl, j2_gl
  real*8           , intent(in) :: londeg(i1_gl:), latdeg(ju1_gl:), coscen(ju1_gl:), mw(numSpecies)
  type(t_Diagnostics), intent(in) :: Diagnostics
  type (t_metFields ), intent(in) :: metFields
!  real*8           , intent(in) :: mcorGlob(i1_gl:i2_gl,ju1_gl:j2_gl)
!
!
! !INPUT/OUTPUT PARAMETERS:
  type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
  type (t_GmiArrayBundle), pointer :: concentration(:)
!
  real*8, parameter :: LAT_LIM  = 10.0d0           ! deg latitude
  real*8, parameter :: MINPRS   = 800.d0           ! hPa
  real*8, parameter :: TIME_LIM = 30.d0*86400.d0   ! sec (30 days)
!
  integer :: ij, il, ik, tracer_opt, procID
  integer :: ix, ic, id, nymd, idumyear, month, day
  integer :: ic_e90, ic_tm25, ic_clock, ic_strat_O3
!
  real*8, allocatable :: press3c(:,:,:)
  real*8, allocatable :: gridBoxHeight(:,:,:)
  real*8, allocatable :: pbl(:, :)
  real*8, allocatable :: kel(:,:,:), humidity(:,:,:)
  real*8  :: tdt, gmi_sec
  real*8  :: age_value
  real*8  :: efol_time
!
  logical :: pr_diag
!
  character (len=MAX_LENGTH_FILE_NAME) :: SO3daily_infile_name, SO3monthly_infile_name
!
!EOP
!-------------------------------------------------------------------------
!BOC
  call Get_procID    (gmiDomain, procID)
  call Get_pr_diag(Diagnostics, pr_diag)
!
  if (pr_diag) Write(6,*) 'calcTaggedCO_AgeOfAir called by ', procID
!
!
!... Obtain model time step
  call Get_gmiTimeStep (gmiClock, tdt)
  call Get_gmiSeconds  (gmiClock, gmi_sec)
  call Get_curGmiDate  (gmiClock, nymd          )
!... get date info
  call GmiSplitDateTime (nymd, idumyear, month, day)
!
  call Get_tracer_opt   (SpeciesConcentration, tracer_opt   )
  call Get_efol_time    (SpeciesConcentration, efol_time    )
  call Get_concentration(SpeciesConcentration, concentration)
!
  opt_tracer: select case (tracer_opt)
!... Add e-folding corrections for the tracers
!... This portion of the code is modified by Bigyani
    case (1) opt_tracer
      if (efol_time .ne. 0.0) then
         do ic=1,numSpecies
            concentration(ic)%pArray3D(:,:,:) = concentration(ic)%pArray3D(:,:,:) &
                                               * exp(-tdt/(86400.0*efol_time))
         enddo
      endif
!
!... Do Age of Air tracer experiment
    case (2) opt_tracer
!
      allocate(press3c(ilo:ihi, julo:jhi, k1:k2))
      call Get_press3c(metFields, press3c)
!... if past TIME_LIM (standard is 30 days) into run set age_value to 0.0
      if (gmi_sec .le. TIME_LIM) then
         age_value = 1.0
      else
         age_value = 0.0
      endif
!... set age tracer to "age_value" in lower tropospheic tropics
      do ik=k1,k2
         do ij=ju1,j2
            if (latdeg(ij) .ge. -LAT_LIM .and. latdeg(ij) .le. LAT_LIM) then
               do il=i1,i2
                  if (press3c(il,ij,ik) .ge. MINPRS) then
                     concentration(1)%pArray3D(il,ij,ik) = age_value
                  endif
               enddo
            endif
         enddo
      enddo
!
!... e90 tracer run (Prather)
      ic_e90 = Ie90
      call GMI_e90  &
         (concentration, ic_e90, tdt, total_mass, layer1_mass,  &
          pr_diag, procID, numSpecies, i1, i2, ju1, j2)
!
!... 25 day stratospheric tracer run (from polmip)
      ic_tm25 = Itm25
      call GMI_tm25  &
            (concentration, ic_tm25, coscen, press3c, tdt, pr_diag, procID, numSpecies,  &
             i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ju1_gl)
!
!... clock tracer for age of air (from GMAO add time step everywhere and total sink at the surface)
      ic_clock = Iclock
      concentration(ic_clock)%pArray3D(:,:,k1+1:k2) = concentration(ic_clock)%pArray3D(:,:,k1+1:k2)+tdt/86400.
      concentration(ic_clock)%pArray3D(:,:,k1) = 0.0
!
      deallocate(press3c)
!
!... CO tagged tracer run like HTAP TP1 specification
    case (3) opt_tracer
!
      call HTAP_TP1_Tagged_CO (concentration, tdt, mw, pr_diag, procID, numSpecies)
!
!... CH4 tagged tracer run as Yasuko Yoshida designed
    case (4) opt_tracer
!
      call GMI_Tagged_CH4 (concentration, tdt, pr_diag, procID, numSpecies)
!
!... CO tagged tracer run as Bryan Duncan designed
    case (5) opt_tracer
!
      call GMI_Tagged_CO(concentration, tdt, pr_diag, procID, numSpecies)
!
!... CO tagged tracer run as Bryan Duncan designed but with all sources doing 5, 10, 15 day efold decay
    case (6) opt_tracer
!
      call GMI_Tagged_CO_wDecay(concentration, tdt, pr_diag, procID, numSpecies)
!
!... CO tagged tracer run as Bryan Duncan designed but with all sources doing 5, 10, 15 day efold decay
    case (7) opt_tracer
!
      allocate(gridBoxHeight(i1:i2,ju1:j2,k1:k2))
      call Get_gridBoxHeight(metFields, gridBoxHeight)
      allocate(pbl(i1:i2,ju1:j2))
      call Get_pbl(metFields, pbl)
!      print *,  &
!       'call GMI_Waugh_tracers(concentration, tdt, pr_diag, procID, numSpecies)'
      call GMI_Waugh_Tracers  &
       (concentration, gridBoxHeight, pbl, latdeg, month, day, tdt, gmi_sec, pr_diag, procID,  &
        numSpecies, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ju1_gl)
!
      deallocate(gridBoxHeight)
      deallocate(pbl)
!
    case (8) opt_tracer
!
      allocate(press3c(ilo:ihi, julo:jhi, k1:k2))
      call Get_press3c(metFields, press3c)
!... e90 and Strat_O3 tracer run
      ic_e90 = 1
      ic_strat_O3 = 2
!
      call GMI_e90  &
         (concentration, ic_e90, tdt, total_mass, layer1_mass,  &
          pr_diag, procID, numSpecies, i1, i2, ju1, j2)
!
!... StratO3 tracer with emissions and parameterized chemistry
      call Get_SO3daily_infile_name   (SpeciesConcentration, SO3daily_infile_name )
      call Get_SO3monthly_infile_name (SpeciesConcentration, SO3monthly_infile_name )
      allocate (kel    (ilo:ihi, julo:jhi, k1:k2))
      call Get_kel  (metFields, kel  )
      allocate (humidity(i1:i2, ju1:j2,k1:k2))
      call Get_humidity(metFields, humidity)
!     =================
      call Update_stratO3  &
!     =================
        (SO3daily_infile_name, SO3monthly_infile_name, ic_strat_O3, ic_e90,  &
         concentration, humidity, kel, press3c, nymd, tdt, londeg, latdeg,  &
         pr_diag, procID, i1, i2, ju1, j2, k1, k2,  &
         ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl, numSpecies)
!
      deallocate(kel)
      deallocate(humidity)
      deallocate(press3c)
!
!... default - do nothing
    case default opt_tracer
  endselect opt_tracer
!
  call Set_concentration(SpeciesConcentration, concentration)
!
  return
!
  end subroutine calcTaggedCO_AgeOfAir
!
!-----------------------------------------------------------------------------
!
! ROUTINE
! GMI_e90
!
! DESCRIPTION
!   This routine calculates Prather e90 tracer in the age of air run
!     Characteristics:
!       - uniform surface emission
!       - global average = 100 ppb
!       - 90 day efold time
!       - tropopause is e90=90ppb (80.6+-0.5% of total mass)
!
! TO RUN USING THIS ROUTINE
!   This is run using tracer_opt=2, constituent index=2
!
!   In resource file set the following (in addition to other things)
!
!      numSpecies  = 2,
!      const_labels::
!      'Age'
!      'e90'
!      ::
!      const_init_val::
!      0.0d+0,
!      0.0d+0,
!      ::
!      mw::
!      1.0d0
!      1.0d0
!      ::
!      tracer_opt: 2
!
! SUBROUTINE ARGUMENTS
!   concentration   : species concentration, known at zone centers (mixing ratio)
!   ic_e90          : constituent number of the e90 tracer
!   tdt             : chemistry time step   (s)
!   total_mass      : mass of the atmosphere (kg)
!   layer1_mass     : mass of the lowest model level of the atmosphere (kg)
!   pr_diag         : diagnostic printing switch
!   procID          : local processor number
!   numSpecies      : dimension of concentration
!   i1, i2, ju1, j2 : horizontal dimensions of concentration
!
  subroutine GMI_e90  &
      (concentration, ic_e90, tdt, total_mass, layer1_mass,  &
       pr_diag, procID, numSpecies, i1, i2, ju1, j2)
!
!
  implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
  type (t_GmiArrayBundle), intent(inOut) :: concentration(numSpecies)
  integer, intent(in) :: ic_e90
  real*8,  intent(in) :: tdt
  real*8,  intent(in) :: total_mass
  real*8,  intent(in) :: layer1_mass
  logical, intent(in) :: pr_diag
  integer, intent(in) :: procID
  integer, intent(in) :: numSpecies
  integer, intent(in) :: i1, i2, ju1, j2
!
!
!     -----------------------
!     Parameter declarations.
!     -----------------------
!
!... efolding time in days for the constituent
  real*8, parameter :: tau = (90.)
!  real*8, parameter :: tau = (75.)
!  real*8, parameter :: tau = (60.)
  real*8, parameter :: lambda = (1.0/(tau * 86400.))
!
  real*8 :: source_e90(i1:i2,ju1:j2)
  real*8 :: tote90_loss
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
  if (pr_diag) then
    Write (6,*) 'GMI_e90 called by ', procID
  endif
!
!... efold loss
  concentration(ic_e90)%pArray3D(:,:,:) =  &
      concentration(ic_e90)%pArray3D(:,:,:) * exp(-tdt * lambda)
!
!... find total mass of atmosphere (kg) * 100 ppbv and efold rate to get total e90 mass lost (kg)
  tote90_loss = 100.0D-9 * total_mass * (1.0-exp(-tdt * lambda))
!
!... e90 lost converted to mixing ratio in bottom level
  source_e90(i1:i2,ju1:j2) = tote90_loss / layer1_mass
!
!... source of e90
  concentration(ic_e90)%pArray3D(i1:i2,ju1:j2,1) =  &
      concentration(ic_e90)%pArray3D(i1:i2,ju1:j2,1) + source_e90(i1:i2,ju1:j2)
!
  return
!
  end subroutine GMI_e90
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
! GMI_tm25
!
! DESCRIPTION
!   As decided in Paris we plan to simulate a stratospheric tracer to look at the strat-trop
!   exchange. As decided the concentration in the stratosphere should be fixed to 200ppbv, and have
!   a 25 day lifetime in the troposphere. As interface between troposphere and stratosphere we
!   decided to use an empirical pressure level as a function of the latitude where annual and zonal
!   mean ozone concentrations approximately exceed the 200ppbv level, based on the implementation
!   of TM5 chemistry in C-IFS. This function reads:
!
!    P = 230-148 (cos(lat))**4 (units: hPa)
!
!    (so P_tropopause = 82 hPa at the equator and 230 hPa at the poles.)
!
! TO RUN USING THIS ROUTINE
!   This is run using tracer_opt=2, constituent index=3
!
!   In resource file set the following (in addition to other things)
!
!      numSpecies  = 3,
!      const_labels::
!      'Age'
!      'e90'
!      'trop-e25'
!      ::
!      const_init_val::
!      0.0d+0,
!      100.0d-9,
!      0.0d+0,
!      ::
!      mw::
!      1.0d0
!      1.0d0
!      1.0d0
!      ::
!      tracer_opt: 2
!
! SUBROUTINE ARGUMENTS
!   concentration     : species concentration, known at zone centers (mixing ratio)
!   tdt       : chemistry time step   (s)
!   mw        : molecular weights of species in const
!   pr_diag   : diagnostic printing switch
!   procID  : local processor number
!   numSpecies   : dimensions of concentration
!
  subroutine GMI_tm25  &
    (concentration, ic_tm25, coscen, press3c, tdt, pr_diag, procID, numSpecies,  &
      i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ju1_gl)
!
!
  implicit none
!
# include "gmi_phys_constants.h"
!
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
  type (t_GmiArrayBundle), intent(inOut) :: concentration(numSpecies)
  integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
  integer, intent(in) :: numSpecies, ju1_gl
  real*8 , intent(in) :: coscen(ju1_gl:)
  logical, intent(in) :: pr_diag
  integer, intent(in) :: procID
  real*8 , intent(in) :: tdt
  real*8 , intent(in) :: press3c(ilo:ihi, julo:jhi, k1:k2)
!
!
!     -----------------------
!     Parameter declarations.
!     -----------------------
!
!... efolding time in days for the constituent
  real*8, parameter  :: tau = (25.)
  real*8, parameter  :: lambda = (1.0/(tau * 86400.))
!
  real*8  :: work(i1:i2, ju1:j2)
  real*8  :: p_limit(i1:i2, ju1:j2, k1:k2)
!
  integer :: ic, il, kk
  integer :: ic_tm25
!
  logical, save :: first = .true.
!
!     ----------------
!     Begin execution.
!     ----------------
!
  if (pr_diag) then
    Write (6,*) 'GMI_tm25 called by ', procID
  endif
!
!... setup array for applying 25 day efold below given pressure level ("tropopause")
  work(i1:i2,ju1:j2) = spread(230.d0 - 148.d0 * (coscen(ju1:j2))**4,1,(i2-i1+1))
  p_limit(i1:i2,ju1:j2,k1:k2) = spread(work,3,(k2-k1+1))
!
!... efold loss and fix stratosphere
   where(press3c(i1:i2,ju1:j2,k1:k2) <= p_limit(i1:i2,ju1:j2,k1:k2))
      concentration(ic_tm25)%pArray3D(i1:i2,ju1:j2,k1:k2) = 200.0D-09
   elsewhere
      concentration(ic_tm25)%pArray3D(i1:i2,ju1:j2,k1:k2) =  &
         concentration(ic_tm25)%pArray3D(i1:i2,ju1:j2,k1:k2) * exp(-tdt * lambda)
   endwhere
!
  return
!
  end subroutine GMI_tm25
!
!EOC
!-----------------------------------------------------------------------------
!
! ROUTINE
! HTAP_TP1_Tagged_CO
!
! DESCRIPTION
!   This routine calculates CO components (regional sources) with
!     transport and simplified chemistry (loss with 25 days and source
!     from fixed methane - as defined in HTAP-TP1
!
! TO RUN USING THIS ROUTINE
!   This is run using tracer_opt=3, start with no restart
!
!   In namelist set the following (in addition to other things)
!    &ACTM_CONTROL
!      numSpecies  = 10,
!    &ACTM_INPUT
!      const_labels(1)        = 'CO',
!      const_labels(2)        = 'AVOC',
!      const_labels(3)        = 'BVOC',
!      const_labels(4)        = 'CO_EA',
!      const_labels(5)        = 'CO_EU',
!      const_labels(6)        = 'CO_NA',
!      const_labels(7)        = 'CO_SA',
!      const_labels(8)        = 'CO_AVOC',
!      const_labels(9)        = 'CO_BVOC',
!      const_labels(10)       = 'CO_CH4',
!      const_init_val(1)      = 0.0d+0,
!      const_init_val(2)      = 0.0d+0,
!      const_init_val(3)      = 0.0d+0,
!      const_init_val(4)      = 0.0d+0,
!      const_init_val(5)      = 0.0d+0,
!      const_init_val(6)      = 0.0d+0,
!      const_init_val(7)      = 0.0d+0,
!      const_init_val(8)      = 0.0d+0,
!      const_init_val(9)      = 0.0d+0,
!      const_init_val(10)     = 0.0d+0,
!      mw(1)                  = 28.0d0,
!      mw(2)                  = 12.0d0,
!      mw(3)                  = 12.0d0,
!      mw(4)                  = 28.0d0,
!      mw(5)                  = 28.0d0,
!      mw(6)                  = 28.0d0,
!      mw(7)                  = 28.0d0,
!      mw(8)                  = 28.0d0,
!      mw(9)                  = 28.0d0,
!      mw(10)                 = 28.0d0
!    &ACTM_EMISS
!      emiss_opt    = 2,        ! Harvard emissions
!      emiss_timpyr = 12,
!      emiss_var_name =  'emiss',
!      emiss_infile_name = 'emiss_TP1_GMI.nc',
!      emiss_map(1)      =  1,  ! CO
!      emiss_map(2)      =  2,  ! AVOC
!      emiss_map(3)      =  3,  ! BVOC
!      emiss_map(4)      =  4,  ! CO_EA
!      emiss_map(5)      =  5,  ! CO_EU
!      emiss_map(6)      =  6,  ! CO_NA
!      emiss_map(7)      =  7,  ! CO_SA
!      emiss_conv_flag = 0,     ! 0:none 1:"emiss_conv_fac" 2-kg/km2-hr to kg/s
!      emiss_in_opt = 2,        ! 2:read in emissions
!    &ACTM_TRAC
!      tracer_opt       = 3
!
! SUBROUTINE ARGUMENTS
!   const     : species concentration, known at zone centers (mixing ratio)
!   tdt       : chemistry time step   (s)
!   mw        : molecular weights of species in const
!   pr_diag   : diagnostic printing switch
!   procID  : local processor number
!   numSpecies   : dimensions of const
!
      subroutine HTAP_TP1_Tagged_CO  &
     &  (concentration, tdt, mw, pr_diag, procID, numSpecies)
!
!
      implicit none
!
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: numSpecies
      real*8, intent(in)  :: mw (numSpecies)
      real*8, intent(in)  :: tdt
      type (t_GmiArrayBundle), intent(inOut) :: concentration(numSpecies)
!
!
!     -----------------------
!     Parameter declarations.
!     -----------------------
!... number of constituents
      integer, parameter :: NUM_HTAP_SP = 10
!
!... efolding time in days for the constituents
      real*8, parameter  :: tau(NUM_HTAP_SP)  &
     &  = (/ 25., 7., 1., 25., 25., 25., 25., 25., 25., 25. /)
!     &  = [25., 7., 1., 25., 25., 25., 25., 25., 25., 25.]
!     &  = [1.e20,1.e20,1.e20,1.e20,1.e20,1.e20,1.e20,1.e20,1.e20,1.e20]
!
      real*8, save  :: lambda(NUM_HTAP_SP)
!
!... decay rate of CH4 in days
      real*8, parameter  :: lambda_ch4 = 1.0/(8.5 * 365. * 86400.)
!... molecular weight of CH4
      real*8, parameter  :: mw_ch4 = 16.
!... fixed concentration of CH4
      real*8, parameter  :: ch4_mixrat = 1760.0e-9
!
      integer :: ic
!
      logical, save :: first = .true.
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'HTAP_TP1_Tagged_CO called by ', procID
      endif
!
!     ==========
      if (first) then
!     ==========
        first = .false.
        where (tau .ne. 0.0)
          lambda(:) = 1.0/(tau(:) * 86400.)
        elsewhere
          lambda(:) = 0.0
        endwhere
      endif
!
!     ---------------------------------------
!     do loss of constituents
!     ---------------------------------------
      do ic = 1, NUM_HTAP_SP
!... efold loss
        concentration(ic)%pArray3D(:,:,:) = concentration(ic)%pArray3D(:,:,:) &
                                             * exp(-tdt * lambda(ic))
      enddo
!
!     ---------------------------------------
!     do oxidation source of last 3 CO constituents
!     ---------------------------------------
!... AVOC source of CO
      concentration(8)%pArray3D(:,:,:) = concentration(8)%pArray3D(:,:,:)  &
     &                                 + concentration(2)%pArray3D(:,:,:)  &
     &                                 * (1.0 - exp(-tdt * lambda(2))) * 0.7 * mw(8)/mw(2)
!
!... BVOC source of CO
      concentration(9)%pArray3D(:,:,:) = concentration(9)%pArray3D(:,:,:)  &
     &                                 + concentration(3)%pArray3D(:,:,:)  &
     &                                 * (1.0 - exp(-tdt * lambda(3))) * 0.4 * mw(9)/mw(3)
!
!... CH4 source of CO
      concentration(10)%pArray3D(:,:,:) = concentration(10)%pArray3D(:,:,:)  &
     &                                  + ch4_mixrat * (1.0 - exp(-tdt * lambda_ch4)) &
     &                                  * 0.86 * mw(10)/mw_ch4
!
      return
!
      end subroutine HTAP_TP1_Tagged_CO
!
!-----------------------------------------------------------------------------
!
! ROUTINE
! GMI_Tagged_CH4
!
! DESCRIPTION
!   This routine calculates "tagged" CO in which different sources of CO
!     are carried as separate tracers.  The direct emissions are treated the
!     same as in a full chemistry simulation.  The oxidation of NMHC are
!     approximated as discussed in Duncan et al. (2007).  Monthly-averaged
!     methane, OH, and rate constants from a full chemistry run are used to
!     calculate the production of CO from methane oxidation and CO loss.
!
! TO RUN USING THIS ROUTINE
!   This is run using tracer_opt=4, I started with restart that had 0.0
!     for all the CO species and the first month's values for the model
!     CH4, OH, AD, kCO, kCH4. Note how these fixed species are not transported.
!
!   In the namelist, FF, BF, and BB stand for fossil fuels, biofuels, and biomass
!     burning, resp.  ROW = rest of world; SAs = Indonchina; NAm = N. America;
!     Eur = Europe; EAs = east Asia; Indo = Indonesia; NAf = northern Africa;
!     SAf = southern Africa; BO = boreal; CO_CH4_Ox = CO from methane oxidation;
!     CO_Biogen = CO from the oxidation of isoprene, methanol and monoterpenes;
!     AD = air density in molec/cm3; kCH4 and KCO are rate constants for reaction
!      w/OH. CO_70 and CO_FixedEm represent tracers of dynamics:  CO_70 uses a
!     boundary condition of 70 ppbv over the whole globe and CO_FixedEm uses a
!     fixed, uniform emission over the whole globe giving 2400 Tg CO/yr.  Both of
!     the dynamic tracers assume a 25 day lifetime.
!
!   The input file is created by the IDL code, tagCO_create_input.pro, which is
!     stored in the master input directory and can be obtained by contacting Bryan
!     Duncan. This code also accounts for the yield of CO per carbon atom oxidized
!     from each NMHC.
!
!   In resource file set the following (in addition to other things)
!     numSpecies  = 20,
!     const_init_val::
!     2.7d-7
!     9.6d-8
!     7.4d-8
!     8.0d-8
!     0.0d+0
!     1.4d-8
!     1.6d-7
!     0.0d+0
!     5.6d-8
!     1.7d-7
!     5.6d-8
!     1.1d-7
!     5.7d-7
!     1.9d-8
!     0.0d+0
!     2.0d-8
!     1.7d-6
!     0.0d+0
!     0.0d+0
!     0.0d+0
!     ::
!     fixed_const_infile_name: 'FIXED_SPECIES_FILENAME',
!     fixed_const_timpyr: 12,
!     fixedConcentrationSpeciesNames::
!     "K"
!     "OH"
!     "AD"
!     ::
!     const_labels::
!     "CH4ANIMLS"
!     "CH4COAL"
!     "CH4GASLEAK"
!     "CH4GASVENT"
!     "CH4HYDV"
!     "CH4HYDZ"
!     "CH4MSW"
!     "CH4SOILABS"
!     "CH4TRMITE"
!     "CH4BOGS"
!     "CH4BURN"
!     "CH4RICEC"
!     "CH4SWAMPS"
!     "CH4TUNDRA"
!     "CH4WETL"
!     "CH4BF"
!     "CH4TOT"
!     "K"
!     "OH"
!     "AD"
!     ::
!     mw::
!     16.0d0
!     16.0d0
!     16.0d0
!     16.0d0
!     16.0d0
!     16.0d0
!     16.0d0
!     16.0d0
!     16.0d0
!     16.0d0
!     16.0d0
!     16.0d0
!     16.0d0
!     16.0d0
!     16.0d0
!     16.0d0
!     16.0d0
!     1.0d0
!     1.0d0
!     1.0d0
!     ::
!     tracer_opt: 4
!     advec_flag::
!     1
!     1
!     1
!     1
!     1
!     1
!     1
!     1
!     1
!     1
!     1
!     1
!     1
!     1
!     1
!     1
!     1
!     0
!     0
!     0
!     ::
!     emiss_timpyr: 12
!     emiss_var_name:  'emiss'
!     emiss_infile_name: '/discover/nobackup/yyoshida/run_gmi/emiss_data/CH4_emiss_tag_2x25.nc'
!     emissionSpeciesNames::
!     "CH4ANIMLS"
!     "CH4COAL"
!     "CH4GASLEAK"
!     "CH4GASVENT"
!     "CH4HYDV"
!     "CH4HYDZ"
!     "CH4MSW"
!     "CH4SOILABS"
!     "CH4TRMITE"
!     "CH4BOGS"
!     "CH4BURN"
!     "CH4RICEC"
!     "CH4SWAMPS"
!     "CH4TUNDRA"
!     "CH4WETL"
!     "CH4BF"
!     "CH4TOT"
!     ::
!     emiss_conv_flag: 0           # 0:none 1:"emiss_conv_fac" 2-kg/km2-hr to kg/s
!     emiss_opt: 2                 # Harvard emissions?
!     emiss_in_opt: 2              # 2:read in emissions
!
! ARGUMENTS
!   tdt       : chemistry time step   (s)
!   const     : species concentration, known at zone centers (mixing ratio)
!   pr_diag   : diagnostic printing switch
!   procID  : local processor number
!   numSpecies   : dimensions of const
!
      subroutine GMI_Tagged_CH4  &
     &  (concentration, tdt, pr_diag, procID, numSpecies)
!
      implicit none
!
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: numSpecies
!
      real*8, intent(in)  :: tdt
!
      type (t_GmiArrayBundle), intent(inOut) :: concentration(numSpecies)
!
!
!     -----------------------
!     Local Parameter declarations.
!     -----------------------
!
!... number of CO constituents that use the k[CO][OH] loss
      integer, parameter :: NUM_BND_SP = 17
!
      integer, parameter :: iK           = 18  ! model CH4 loss reaction rate
      integer, parameter :: iOH          = 19  ! model OH
      integer, parameter :: iAD          = 20  ! model Air Density
!
!
      integer :: ic
!
!
!     ==========
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'GMI_Tagged_CH4 called by ', procID
      endif
!
!     ---------------------------------------
!     do loss of CO for the first 17 species in const using OH loss (const(
!     ---------------------------------------
!
!
!    CH4 loss=K*OH*AD*CH4 (yy 6/15/08)
      do ic = 1, NUM_BND_SP
!
        concentration(ic)%pArray3D(:,:,:) = concentration(ic  )%pArray3D(:,:,:)  &
     &                            - ( tdt * concentration(iK)%pArray3D(:,:,:)  &
     &                                    * concentration(ic  )%pArray3D(:,:,:)    &
     &                                    * concentration(iOH )%pArray3D(:,:,:)    &
     &                                    * concentration(iAD )%pArray3D(:,:,:)  )
!
      enddo
!
      return
!
      end subroutine GMI_Tagged_CH4
!-------------------------------------------------------------------------
!
! ROUTINE
! GMI_Tagged_CO
!
! DESCRIPTION
!   This routine calculates "tagged" CO in which different sources of CO
!     are carried as separate tracers.  The direct emissions are treated the
!     same as in a full chemistry simulation.  The oxidation of NMHC are
!     approximated as discussed in Duncan et al. (2007).  Monthly-averaged
!     methane, OH, and rate constants from a full chemistry run are used to
!     calculate the production of CO from methane oxidation and CO loss.
!
! TO RUN USING THIS ROUTINE
!   This is run using tracer_opt=5, I started with restart that had 0.0
!     for all the CO species and the first month's values for the model
!     CH4, OH, AD, kCO, kCH4. Note how these fixed species are not transported.
!
!   In the namelist, FF, BF, and BB stand for fossil fuels, biofuels, and biomass
!     burning, resp.  ROW = rest of world; SAs = Indonchina; NAm = N. America;
!     Eur = Europe; EAs = east Asia; Indo = Indonesia; NAf = northern Africa;
!     SAf = southern Africa; BO = boreal; CO_CH4_Ox = CO from methane oxidation;
!     CO_Biogen = CO from the oxidation of isoprene, methanol and monoterpenes;
!     AD = air density in molec/cm3; kCH4 and KCO are rate constants for reaction
!      w/OH. CO_70 and CO_FixedEm represent tracers of dynamics:  CO_70 uses a
!     boundary condition of 70 ppbv over the whole globe and CO_FixedEm uses a
!     fixed, uniform emission over the whole globe giving 2400 Tg CO/yr.  Both of
!     the dynamic tracers assume a 25 day lifetime.
!
!   The input file is created by the IDL code, tagCO_create_input.pro, which is
!     stored in the master input directory and can be obtained by contacting Bryan
!     Duncan. This code also accounts for the yield of CO per carbon atom oxidized
!     from each NMHC.
!
!   In namelist set the following (in addition to other things)
!    &nlGmiControl
!      numSpecies  = 24,
!    &nlGmiSpeciesConcentration
!      const_labels( 1)       = "CO_FF_ROW",
!      const_labels( 2)       = "CO_FF_SAs",
!      const_labels( 3)       = "CO_FF_NAm",
!      const_labels( 4)       = "CO_FF_Eur",
!      const_labels( 5)       = "CO_FF_EAs",
!      const_labels( 6)       = "CO_BF_ROW",
!      const_labels( 7)       = "CO_BF_SAs",
!      const_labels( 8)       = "CO_BF_EAs",
!      const_labels( 9)       = "CO_BB_ROW",
!      const_labels(10)       = "CO_BB_SAs",
!      const_labels(11)       = "CO_BB_Indo",
!      const_labels(12)       = "CO_BB_SAm",
!      const_labels(13)       = "CO_BB_SAf",
!      const_labels(14)       = "CO_BB_NAf",
!      const_labels(15)       = "CO_BB_BO",
!      const_labels(16)       = "CO_Biogen",
!      const_labels(17)       = "CO_CH4_Ox",
!      const_labels(18)       = "CH4",
!      const_labels(19)       = "OH",
!      const_labels(20)       = "AD",
!      const_labels(21)       = "kCH4",
!      const_labels(22)       = "kCO",
!      const_labels(23)       = "CO_70",
!      const_labels(24)       = "CO_FixedEm",
!      fixed_const_map( 1)     = 18,
!      fixed_const_map( 2)     = 19,
!      fixed_const_map( 3)     = 20,
!      fixed_const_map( 4)     = 21,
!      fixed_const_map( 5)     = 22,
!      fixed_const_infile_name = 'FIXED_SPECIES_FILENAME',
!      fixed_const_timpyr      = 12,
!      const_init_val( 1)     = 0.0d+0,
!      const_init_val( 2)     = 0.0d+0,
!      const_init_val( 3)     = 0.0d+0,
!      const_init_val( 4)     = 0.0d+0,
!      const_init_val( 5)     = 0.0d+0,
!      const_init_val( 6)     = 0.0d+0,
!      const_init_val( 7)     = 0.0d+0,
!      const_init_val( 8)     = 0.0d+0,
!      const_init_val( 9)     = 0.0d+0,
!      const_init_val(10)     = 0.0d+0,
!      const_init_val(11)     = 0.0d+0,
!      const_init_val(12)     = 0.0d+0,
!      const_init_val(13)     = 0.0d+0,
!      const_init_val(14)     = 0.0d+0,
!      const_init_val(15)     = 0.0d+0,
!      const_init_val(16)     = 0.0d+0,
!      const_init_val(17)     = 0.0d+0,
!      const_init_val(23)     = 0.0d+0,
!      const_init_val(24)     = 0.0d+0,
!      mw( 1)                 = 28.0d0,
!      mw( 2)                 = 28.0d0,
!      mw( 3)                 = 28.0d0,
!      mw( 4)                 = 28.0d0,
!      mw( 5)                 = 28.0d0,
!      mw( 6)                 = 28.0d0,
!      mw( 7)                 = 28.0d0,
!      mw( 8)                 = 28.0d0,
!      mw( 9)                 = 28.0d0,
!      mw(10)                 = 28.0d0,
!      mw(11)                 = 28.0d0,
!      mw(12)                 = 28.0d0,
!      mw(13)                 = 28.0d0,
!      mw(14)                 = 28.0d0,
!      mw(15)                 = 28.0d0,
!      mw(16)                 = 28.0d0,
!      mw(17)                 = 28.0d0,
!      mw(18)                 = 16.0d0,
!      mw(19)                 = 17.0d0,
!      mw(20)                 =  1.0d0,
!      mw(21)                 =  1.0d0,
!      mw(22)                 =  1.0d0
!      mw(23)                 = 28.0d0
!      mw(24)                 = 28.0d0
!    &nlGmiAdvection
!      advec_flag(1)  = 1,
!      advec_flag(2)  = 1,
!      advec_flag(3)  = 1,
!      advec_flag(4)  = 1,
!      advec_flag(5)  = 1,
!      advec_flag(6)  = 1,
!      advec_flag(7)  = 1,
!      advec_flag(8)  = 1,
!      advec_flag(9)  = 1,
!      advec_flag(10) = 1,
!      advec_flag(11) = 1,
!      advec_flag(12) = 1,
!      advec_flag(13) = 1,
!      advec_flag(14) = 1,
!      advec_flag(15) = 1,
!      advec_flag(16) = 1,
!      advec_flag(17) = 1,
!      advec_flag(18) = 0,
!      advec_flag(19) = 0,
!      advec_flag(20) = 0,
!      advec_flag(21) = 0,
!      advec_flag(22) = 0,
!      advec_flag(23) = 1,
!      advec_flag(24) = 1,
!    &nlGmiEmission
!      emiss_opt    = 2,        ! Harvard emissions?
!      emiss_timpyr = 12,
!      emiss_var_name =  'emiss',
!      emiss_infile_name = 'EMISS_FILENAME',
!      emiss_map(1)      =  1,  ! CO_FF_ROW"
!      emiss_map(2)      =  2,  ! CO_FF_SAs"
!      emiss_map(3)      =  3,  ! CO_FF_NAm"
!      emiss_map(4)      =  4,  ! CO_FF_Eur"
!      emiss_map(5)      =  5,  ! CO_FF_EAs"
!      emiss_map(6)      =  6,  ! CO_BF_ROW"
!      emiss_map(7)      =  7,  ! CO_BF_SAs"
!      emiss_map(8)      =  8,  ! CO_BF_EAs"
!      emiss_map(9)      =  9,  ! CO_BB_ROW"
!      emiss_map(10)     = 10,  ! CO_BB_SAs"
!      emiss_map(11)     = 11,  ! CO_BB_Indo"
!      emiss_map(12)     = 12,  ! CO_BB_SAm"
!      emiss_map(13)     = 13,  ! CO_BB_SAf"
!      emiss_map(14)     = 14,  ! CO_BB_NAf"
!      emiss_map(15)     = 15,  ! CO_BB_BO"
!      emiss_map(16)     = 16,  ! CO_Biogen"
!      emiss_map(18)     = 17,  ! CO_CH4_Ox"
!      emiss_map(17)     = 24,  ! FixedEm"
!      emiss_conv_flag = 0,     ! 0:none 1:"emiss_conv_fac" 2-kg/km2-hr to kg/s
!      emiss_in_opt = 2,        ! 2:read in emissions
!    &ACTM_TRAC
!      tracer_opt       = 5
!
! ARGUMENTS
!   tdt       : chemistry time step   (s)
!   const     : species concentration, known at zone centers (mixing ratio)
!   pr_diag   : diagnostic printing switch
!   procID  : local processor number
!   numSpecies   : dimensions of const
!
      subroutine GMI_Tagged_CO  &
     &  (concentration, tdt, pr_diag, procID, numSpecies)
!
      implicit none
!
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: numSpecies
!
      real*8, intent(in)  :: tdt
!
      type (t_GmiArrayBundle), intent(inOut) :: concentration(numSpecies)
!
!
!     -----------------------
!     Local Parameter declarations.
!     -----------------------
!
!... number of CO constituents that use the k[CO][OH] loss
      integer, parameter :: NUM_BND_SP = 17
!
      integer, parameter :: iCO_CH4_Ox  = 17
      integer, parameter :: iCH4        = 18  ! full chemistry model CH4
      integer, parameter :: iOH         = 19  ! model OH
      integer, parameter :: iAD         = 20  ! model Air Density
      integer, parameter :: ikCH4       = 21  ! model CH4 source of CO reaction rate
      integer, parameter :: ikCO        = 22  ! model CO loss reaction rate
      integer, parameter :: iCO_70      = 23  ! CO bottom layer fixed to sfc_CO_70
      integer, parameter :: iFixedEm    = 24  ! CO w/ fixed global emission of CO
!... decay rate of Co with 70 ppbv boundary condition - 25 days
      real*8, parameter  :: lambda_CO_70 = 1.0/(25. * 86400.)
      real*8, parameter  :: sfc_CO_70 = 70.0e-9
!
!
      integer :: ic
!
!
!     ==========
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'GMI_Tagged_CO called by ', procID
      endif
!
!     ---------------------------------------
!     do loss of CO for the first 17 species in const using OH loss (const(
!     ---------------------------------------
!
!... CO loss for tagged CO components
      do ic = 1, NUM_BND_SP
        concentration(ic)%pArray3D(:,:,:) = concentration(ic  )%pArray3D(:,:,:)  &
     &                            - ( tdt * concentration(ikCO)%pArray3D(:,:,:)  &
     &                                    * concentration(ic  )%pArray3D(:,:,:)    &
     &                                    * concentration(iOH )%pArray3D(:,:,:)    &
     &                                    * concentration(iAD )%pArray3D(:,:,:)  )
      enddo
!
!    ---------------------------------------------
!     do oxidation source of CO due to CH4 oxidation
!    ---------------------------------------------
!... CH4 source of CO
      concentration(iCO_CH4_Ox)%pArray3D(:,:,:) = concentration(iCO_CH4_Ox)%pArray3D(:,:,:)  &
     &                                  + ( tdt * concentration(ikCH4     )%pArray3D(:,:,:)  &
     &                                          * concentration(iCH4      )%pArray3D(:,:,:) &
     &                                          * concentration(iOH       )%pArray3D(:,:,:) &
     &                                          * concentration(iAD       )%pArray3D(:,:,:) )
!
!
!... Now do simple CO with bottom fixed to 70ppb - special case
!... loss for CO with fixed boundary condition
      concentration(iCO_70)%pArray3D(:,:,:) = concentration(iCO_70)%pArray3D(:,:,:) &
                                            * exp(-tdt * lambda_CO_70)
!... fix bottom level of iCO_70
      concentration(iCO_70)%pArray3D(:,:,1) = sfc_CO_70
!
!
!... Now do simple CO with fixed global emissions - special case
!... loss for CO with fixed emissions
      concentration(iFixedEm)%pArray3D(:,:,:) = concentration(iFixedEm)%pArray3D(:,:,:) &
                                              * exp(-tdt * lambda_CO_70)
!
!
      return
!
      end subroutine GMI_Tagged_CO
!
!-------------------------------------------------------------------------
!
! ROUTINE
! GMI_Tagged_CO_wDecay
!
! DESCRIPTION
!   This routine calculates "tagged" CO in which different sources of CO
!     are carried as separate tracers.  The direct emissions are treated the
!     same as in a full chemistry simulation.  The oxidation of NMHC are
!     approximated as discussed in Duncan et al. (2007).  Monthly-averaged
!     methane, OH, and rate constants from a full chemistry run are used to
!     calculate the production of CO from methane oxidation and CO loss.
!
! TO RUN USING THIS ROUTINE
!   This is run using tracer_opt=6, I started with restart that had 0.0
!     for all the CO species and the first month's values for the model
!     CH4, OH, AD, kCO, kCH4. Note how these fixed species are not transported.
!
!   In the namelist, FF, BF, and BB stand for fossil fuels, biofuels, and biomass
!     burning, resp.  ROW = rest of world; SAs = Indonchina; NAm = N. America;
!     Eur = Europe; EAs = east Asia; Indo = Indonesia; NAf = northern Africa;
!     SAf = southern Africa; BO = boreal; CO_CH4_Ox = CO from methane oxidation;
!     CO_Biogen = CO from the oxidation of isoprene, methanol and monoterpenes;
!     AD = air density in molec/cm3; kCH4 and KCO are rate constants for reaction
!      w/OH. CO_70 and CO_FixedEm represent tracers of dynamics:  CO_70 uses a
!     boundary condition of 70 ppbv over the whole globe and CO_FixedEm uses a
!     fixed, uniform emission over the whole globe giving 2400 Tg CO/yr.  Both of
!     the dynamic tracers assume a 25 day lifetime.
!
!   The input file is created by the IDL code, tagCO_create_input.pro, which is
!     stored in the master input directory and can be obtained by contacting Bryan
!     Duncan. This code also accounts for the yield of CO per carbon atom oxidized
!     from each NMHC.
!
!   In namelist set the following (in addition to other things)
!    &nlGmiControl
!      numSpecies  = 72,
!    &nlGmiSpeciesConcentration
!      const_labels( 1)       = "CO_FF_ROW",
!      const_labels( 2)       = "CO_FF_SAs",
!      const_labels( 3)       = "CO_FF_NAm",
!      const_labels( 4)       = "CO_FF_Eur",
!      const_labels( 5)       = "CO_FF_EAs",
!      const_labels( 6)       = "CO_BF_ROW",
!      const_labels( 7)       = "CO_BF_SAs",
!      const_labels( 8)       = "CO_BF_EAs",
!      const_labels( 9)       = "CO_BB_ROW",
!      const_labels(10)       = "CO_BB_SAs",
!      const_labels(11)       = "CO_BB_Indo",
!      const_labels(12)       = "CO_BB_SAm",
!      const_labels(13)       = "CO_BB_SAf",
!      const_labels(14)       = "CO_BB_NAf",
!      const_labels(15)       = "CO_BB_BO",
!      const_labels(16)       = "CO_Biogen",
!      const_labels(17)       = "CO_CH4_Ox",
!      const_labels(18)       = "CO_FF_ROW_5d",
!      const_labels(19)       = "CO_FF_SAs_5d",
!      const_labels(20)       = "CO_FF_NAm_5d",
!      const_labels(21)       = "CO_FF_Eur_5d",
!      const_labels(22)       = "CO_FF_EAs_5d",
!      const_labels(23)       = "CO_BF_ROW_5d",
!      const_labels(24)       = "CO_BF_SAs_5d",
!      const_labels(25)       = "CO_BF_EAs_5d",
!      const_labels(26)       = "CO_BB_ROW_5d",
!      const_labels(27)       = "CO_BB_SAs_5d",
!      const_labels(28)       = "CO_BB_Indo_5d",
!      const_labels(29)       = "CO_BB_SAm_5d",
!      const_labels(30)       = "CO_BB_SAf_5d",
!      const_labels(31)       = "CO_BB_NAf_5d",
!      const_labels(32)       = "CO_BB_BO_5d",
!      const_labels(33)       = "CO_Biogen_5d",
!      const_labels(34)       = "CO_FF_ROW_10d",
!      const_labels(35)       = "CO_FF_SAs_10d",
!      const_labels(36)       = "CO_FF_NAm_10d",
!      const_labels(37)       = "CO_FF_Eur_10d",
!      const_labels(38)       = "CO_FF_EAs_10d",
!      const_labels(39)       = "CO_BF_ROW_10d",
!      const_labels(40)       = "CO_BF_SAs_10d",
!      const_labels(41)       = "CO_BF_EAs_10d",
!      const_labels(42)       = "CO_BB_ROW_10d",
!      const_labels(43)       = "CO_BB_SAs_10d",
!      const_labels(44)       = "CO_BB_Indo_10d",
!      const_labels(45)       = "CO_BB_SAm_10d",
!      const_labels(46)       = "CO_BB_SAf_10d",
!      const_labels(47)       = "CO_BB_NAf_10d",
!      const_labels(48)       = "CO_BB_BO_10d",
!      const_labels(49)       = "CO_Biogen_10d",
!      const_labels(50)       = "CO_FF_ROW_15d",
!      const_labels(51)       = "CO_FF_SAs_15d",
!      const_labels(52)       = "CO_FF_NAm_15d",
!      const_labels(53)       = "CO_FF_Eur_15d",
!      const_labels(54)       = "CO_FF_EAs_15d",
!      const_labels(55)       = "CO_BF_ROW_15d",
!      const_labels(56)       = "CO_BF_SAs_15d",
!      const_labels(57)       = "CO_BF_EAs_15d",
!      const_labels(58)       = "CO_BB_ROW_15d",
!      const_labels(59)       = "CO_BB_SAs_15d",
!      const_labels(60)       = "CO_BB_Indo_15d",
!      const_labels(61)       = "CO_BB_SAm_15d",
!      const_labels(62)       = "CO_BB_SAf_15d",
!      const_labels(63)       = "CO_BB_NAf_15d",
!      const_labels(64)       = "CO_BB_BO_15d",
!      const_labels(65)       = "CO_Biogen_15d",
!      const_labels(66)       = "CH4",
!      const_labels(67)       = "OH",
!      const_labels(68)       = "AD",
!      const_labels(69)       = "kCH4",
!      const_labels(70)       = "kCO",
!      const_labels(71)       = "CO_70",
!      const_labels(72)       = "CO_FixedEm",
!      fixed_const_map(66)     = CH4",
!      fixed_const_map(67)     = OH",
!      fixed_const_map(68)     = AD",
!      fixed_const_map(69)     = kCH4,
!      fixed_const_map(70)     = kCO",
!      fixed_const_infile_name = 'FIXED_SPECIES_FILENAME',
!      fixed_const_timpyr      = 12,
!      const_init_val( 1)     = 0.0d+0,
!      const_init_val( 2)     = 0.0d+0,
!      const_init_val( 3)     = 0.0d+0,
!      const_init_val( 4)     = 0.0d+0,
!      const_init_val( 5)     = 0.0d+0,
!      const_init_val( 6)     = 0.0d+0,
!      const_init_val( 7)     = 0.0d+0,
!      const_init_val( 8)     = 0.0d+0,
!      const_init_val( 9)     = 0.0d+0,
!      const_init_val(10)     = 0.0d+0,
!      const_init_val(11)     = 0.0d+0,
!      const_init_val(12)     = 0.0d+0,
!      const_init_val(13)     = 0.0d+0,
!      const_init_val(14)     = 0.0d+0,
!      const_init_val(15)     = 0.0d+0,
!      const_init_val(16)     = 0.0d+0,
!      const_init_val(17)     = 0.0d+0,
!      ...
!      const_init_val(71)     = 0.0d+0,
!      const_init_val(72)     = 0.0d+0,
!      mw( 1)                 = 28.0d0,
!      mw( 2)                 = 28.0d0,
!      mw( 3)                 = 28.0d0,
!      mw( 4)                 = 28.0d0,
!      mw( 5)                 = 28.0d0,
!      mw( 6)                 = 28.0d0,
!      mw( 7)                 = 28.0d0,
!      mw( 8)                 = 28.0d0,
!      mw( 9)                 = 28.0d0,
!      mw(10)                 = 28.0d0,
!      mw(11)                 = 28.0d0,
!      mw(12)                 = 28.0d0,
!      mw(13)                 = 28.0d0,
!      mw(14)                 = 28.0d0,
!      mw(15)                 = 28.0d0,
!      mw(16)                 = 28.0d0,
!      mw(17)                 = 28.0d0,
!      ...
!      mw(65)                 = 28.0d0,
!      mw(66)                 = 16.0d0,
!      mw(67)                 = 17.0d0,
!      mw(68)                 =  1.0d0,
!      mw(69)                 =  1.0d0,
!      mw(70)                 =  1.0d0
!      mw(71)                 = 28.0d0
!      mw(72)                 = 28.0d0
!    &nlGmiAdvection
!      advec_flag(1)  = 1,
!      advec_flag(2)  = 1,
!      advec_flag(3)  = 1,
!      advec_flag(4)  = 1,
!      advec_flag(5)  = 1,
!      advec_flag(6)  = 1,
!      advec_flag(7)  = 1,
!      advec_flag(8)  = 1,
!      advec_flag(9)  = 1,
!      advec_flag(10) = 1,
!      advec_flag(11) = 1,
!      advec_flag(12) = 1,
!      advec_flag(13) = 1,
!      advec_flag(14) = 1,
!      advec_flag(15) = 1,
!      advec_flag(16) = 1,
!      advec_flag(17) = 1,
!      advec_flag(18) = 1,
!      advec_flag(19) = 1,
!      ...
!      advec_flag(65) = 1,
!      advec_flag(66) = 0,
!      advec_flag(67) = 0,
!      advec_flag(68) = 0,
!      advec_flag(69) = 0,
!      advec_flag(70) = 0,
!      advec_flag(71) = 1,
!      advec_flag(72) = 1,
!    &nlGmiEmission
!      emiss_opt    = 2,        ! Harvard emissions?
!      emiss_timpyr = 12,
!      emiss_var_name =  'emiss',
!      emiss_infile_name = 'EMISS_FILENAME',
!      emiss_map(1)      =  1,  ! CO_FF_ROW"
!      emiss_map(2)      =  2,  ! CO_FF_SAs"
!      emiss_map(3)      =  3,  ! CO_FF_NAm"
!      emiss_map(4)      =  4,  ! CO_FF_Eur"
!      emiss_map(5)      =  5,  ! CO_FF_EAs"
!      emiss_map(6)      =  6,  ! CO_BF_ROW"
!      emiss_map(7)      =  7,  ! CO_BF_SAs"
!      emiss_map(8)      =  8,  ! CO_BF_EAs"
!      emiss_map(9)      =  9,  ! CO_BB_ROW"
!      emiss_map(10)     = 10,  ! CO_BB_SAs"
!      emiss_map(11)     = 11,  ! CO_BB_Indo"
!      emiss_map(12)     = 12,  ! CO_BB_SAm"
!      emiss_map(13)     = 13,  ! CO_BB_SAf"
!      emiss_map(14)     = 14,  ! CO_BB_NAf"
!      emiss_map(15)     = 15,  ! CO_BB_BO"
!      emiss_map(16)     = 16,  ! CO_Biogen"
!      emiss_map(18)     = 17,  ! CO_CH4_Ox"
!      emiss_map(17)     = 72,  ! FixedEm"
!      emiss_conv_flag = 0,     ! 0:none 1:"emiss_conv_fac" 2-kg/km2-hr to kg/s
!      emiss_in_opt = 2,        ! 2:read in emissions
!    &ACTM_TRAC
!      tracer_opt       = 6
!
! ARGUMENTS
!   tdt       : chemistry time step   (s)
!   const     : species concentration, known at zone centers (mixing ratio)
!   pr_diag   : diagnostic printing switch
!   procID  : local processor number
!   numSpecies   : dimensions of const
!
      subroutine GMI_Tagged_CO_wDecay  &
     &  (concentration, tdt, pr_diag, procID, numSpecies)
!
      implicit none
!
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: numSpecies
!
      real*8, intent(in)  :: tdt
!
      type (t_GmiArrayBundle), intent(inOut) :: concentration(numSpecies)
!
!
!     -----------------------
!     Local Parameter declarations.
!     -----------------------
!
!... number of CO constituents that use the k[CO][OH] loss
      integer, parameter :: NUM_BND_SP = 17
      integer, parameter :: NUM_CO_SRC_SP = 16
      integer, parameter :: NUM_BEG_05D_SP = 18
      integer, parameter :: NUM_BEG_10D_SP = 34
      integer, parameter :: NUM_BEG_15D_SP = 50
!
      integer, parameter :: iCO_CH4_Ox  = 17
      integer, parameter :: iCH4        = 66  ! full chemistry model CH4
      integer, parameter :: iOH         = 67  ! model OH
      integer, parameter :: iAD         = 68  ! model Air Density
      integer, parameter :: ikCH4       = 69  ! model CH4 source of CO reaction rate
      integer, parameter :: ikCO        = 70  ! model CO loss reaction rate
      integer, parameter :: iCO_70      = 71  ! CO bottom layer fixed to sfc_CO_70
      integer, parameter :: iFixedEm    = 72  ! CO w/ fixed global emission of CO
!... decay rate of CO special cases
      real*8, parameter  :: lambda_CO_05 = 1.0/( 5. * 86400.)
      real*8, parameter  :: lambda_CO_10 = 1.0/(10. * 86400.)
      real*8, parameter  :: lambda_CO_15 = 1.0/(15. * 86400.)
!... decay rate of Co with 70 ppbv boundary condition - 25 days
      real*8, parameter  :: lambda_CO_70 = 1.0/(25. * 86400.)
      real*8, parameter  :: sfc_CO_70 = 70.0e-9
!
      integer :: ic
!
!
!     ==========
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'GMI_Tagged_CO_wDecay called by ', procID
      endif
!
!     ---------------------------------------
!     do loss of CO for the first 17 species in const using OH loss (const)
!     ---------------------------------------
!
!... CO loss for tagged source region CO components and from CH4 oxidation
      do ic = 1, NUM_CO_SRC_SP+1
        concentration(ic)%pArray3D(:,:,:) = concentration(ic  )%pArray3D(:,:,:)  &
     &                            - ( tdt * concentration(ikCO)%pArray3D(:,:,:)  &
     &                                    * concentration(ic  )%pArray3D(:,:,:)    &
     &                                    * concentration(iOH )%pArray3D(:,:,:)    &
     &                                    * concentration(iAD )%pArray3D(:,:,:)  )
      enddo
!
!    ---------------------------------------------
!     do oxidation source of CO due to CH4 oxidation
!    ---------------------------------------------
!... CH4 source of CO
      concentration(iCO_CH4_Ox)%pArray3D(:,:,:) = concentration(iCO_CH4_Ox)%pArray3D(:,:,:)  &
     &                                  + ( tdt * concentration(ikCH4     )%pArray3D(:,:,:)  &
     &                                          * concentration(iCH4      )%pArray3D(:,:,:) &
     &                                          * concentration(iOH       )%pArray3D(:,:,:) &
     &                                          * concentration(iAD       )%pArray3D(:,:,:) )
!
!
!... Now do simple CO with different sources and set lifetime
      do ic = 0, NUM_CO_SRC_SP-1
!... Now do simple CO with 5 day loss and different sources
        concentration(NUM_BEG_05D_SP+ic)%pArray3D(:,:,:) = concentration(NUM_BEG_05D_SP+ic)%pArray3D(:,:,:) &
                                            * exp(-tdt * lambda_CO_05)
!... Now do simple CO with 5 day loss and different sources
        concentration(NUM_BEG_10D_SP+ic)%pArray3D(:,:,:) = concentration(NUM_BEG_10D_SP+ic)%pArray3D(:,:,:) &
                                            * exp(-tdt * lambda_CO_10)
!... Now do simple CO with 5 day loss and different sources
        concentration(NUM_BEG_15D_SP+ic)%pArray3D(:,:,:) = concentration(NUM_BEG_15D_SP+ic)%pArray3D(:,:,:) &
                                            * exp(-tdt * lambda_CO_15)
      enddo
!
!
!... Now do simple CO with bottom fixed to 70ppb - special case
!... loss for CO with fixed boundary condition
      concentration(iCO_70)%pArray3D(:,:,:) = concentration(iCO_70)%pArray3D(:,:,:) &
                                            * exp(-tdt * lambda_CO_70)
!... fix bottom level of iCO_70
      concentration(iCO_70)%pArray3D(:,:,1) = sfc_CO_70
!
!
!... Now do simple CO with fixed global emissions - special case
!... loss for CO with fixed emissions
      concentration(iFixedEm)%pArray3D(:,:,:) = concentration(iFixedEm)%pArray3D(:,:,:) &
                                              * exp(-tdt * lambda_CO_70)
!
!
      return
!
      end subroutine GMI_Tagged_CO_wDecay
!
!-------------------------------------------------------------------------
!
! ROUTINE
! GMI_Waugh_Tracers
!
! DESCRIPTION
!   This routine calculates "tagged" CO in which different sources of CO
!     are carried as separate tracers.  The direct emissions are treated the
!     same as in a full chemistry simulation.  The oxidation of NMHC are
!     approximated as discussed in Duncan et al. (2007).  Monthly-averaged
!     methane, OH, and rate constants from a full chemistry run are used to
!     calculate the production of CO from methane oxidation and CO loss.
!
! TO RUN USING THIS ROUTINE
!   This is run using tracer_opt=6, I started with restart that had 0.0
!     for all the CO species and the first month's values for the model
!     CH4, OH, AD, kCO, kCH4. Note how these fixed species are not transported.
!
!   In the namelist, FF, BF, and BB stand for fossil fuels, biofuels, and biomass
!     burning, resp.  ROW = rest of world; SAs = Indonchina; NAm = N. America;
!     Eur = Europe; EAs = east Asia; Indo = Indonesia; NAf = northern Africa;
!     SAf = southern Africa; BO = boreal; CO_CH4_Ox = CO from methane oxidation;
!     CO_Biogen = CO from the oxidation of isoprene, methanol and monoterpenes;
!     AD = air density in molec/cm3; kCH4 and KCO are rate constants for reaction
!      w/OH. CO_70 and CO_FixedEm represent tracers of dynamics:  CO_70 uses a
!     boundary condition of 70 ppbv over the whole globe and CO_FixedEm uses a
!     fixed, uniform emission over the whole globe giving 2400 Tg CO/yr.  Both of
!     the dynamic tracers assume a 25 day lifetime.
!
!   The input file is created by the IDL code, tagCO_create_input.pro, which is
!     stored in the master input directory and can be obtained by contacting Bryan
!     Duncan. This code also accounts for the yield of CO per carbon atom oxidized
!     from each NMHC.
!
!   In resource file set the following (in addition to other things)
!  numSpecies: 11
!  const_opt: 1
!  const_init_val:: 
!  1.0d-30
!  1.0d-30
!  1.0d-30
!  1.0d-30
!  1.0d-30
!  1.0d-30
!  1.0d-30
!  1.0d-30
!  1.0d-30
!  1.0d-30
!  1.0d-30
!  ::
!  const_labels::
!  SVOC
!  MVOC
!  LVOC
!  SVOC_nh
!  MVOC_nh
!  LVOC_nh
!  clock_nh
!  Pulse_Jan_nh
!  Pulse_Apr_nh
!  Pulse_Jul_nh
!  Pulse_Oct_nh
!  ::
!  mw:: 
!  28.0d0
!  28.0d0
!  28.0d0
!  28.0d0
!  28.0d0
!  28.0d0
!  28.0d0
!  28.0d0
!  28.0d0
!  28.0d0
!  28.0d0
!  ::
!  
!  tracer_opt: 7
!
!  emiss_opt: 1
!  emiss_in_opt: 2,
!  emiss_timpyr: 12,
!  emiss_var_name:  'emiss_2d',
!  emiss_infile_name: 'EMISS_FILENAME'
!  emissionSpeciesNames::
!  SVOC
!  MVOC
!  LVOC
!  SVOC_nh
!  MVOC_nh
!  LVOC_nh
!  ::
!  veg_infile_name: /discover/nobackup/projects/gmi/gmidata2/input/emissions/common/vegtype_2x2.5_dao_corrected.asc 
!  lai_infile_name: /discover/nobackup/projects/gmi/gmidata2/input/emissions/common/lai_2x2.5_dao_corrected.asc 
!
! ARGUMENTS
!   concentration : species concentration, known at zone centers (mixing ratio)
!   gridBoxHeight : grid box height (m)
!   pbl           : planetary boundary layer (m)
!   latdeg        : grid latitudes (degrees north)
!   month         : simulation month number
!   day           : simulation day number
!   tdt           : chemistry time step (s)
!   gmi_sec       : elapsed time of simulation (s)
!   pr_diag       : diagnostic printing switch
!   procID        : local processor number
!   numSpecies    : species dimensions of concentration
!   i1, i2, ju1, j2, k1, k2 : sub-domain dimensions (without ghost zones)
!   ilo, ihi, julo, jhi     : sub-domain dimensions (with ghost zones)
!   ju1_gl                  : global grid latitude dimension
!
  subroutine GMI_Waugh_Tracers  &
    (concentration, gridBoxHeight, pbl, latdeg, month, day, tdt, gmi_sec, pr_diag, procID,  &
      numSpecies, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ju1_gl)
!
  implicit none
!
# include "gmi_time_constants.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
  type (t_GmiArrayBundle), intent(inOut) :: concentration(numSpecies)
!
  real*8, intent(in) :: gridBoxHeight(i1:i2, ju1:j2, k1:k2)
  real*8, intent(in) :: pbl(i1:i2, ju1:j2)
  real*8, intent(in) :: latdeg(ju1_gl:)
  real*8, intent(in) :: tdt
  real*8, intent(in) :: gmi_sec
!
  logical, intent(in) :: pr_diag
  integer, intent(in) :: month, day
  integer, intent(in) :: procID
  integer, intent(in) :: numSpecies
  integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ju1_gl
!
!
!     -----------------------
!     Local Parameter declarations.
!     -----------------------
!
!...
!
  real*8, parameter :: SLAT_LIM = 30.0d0           ! deg latitude southern limit for pulses and clock
  real*8, parameter :: NLAT_LIM = 50.0d0           ! deg latitude northern limit for pulses and clock
  real*8, parameter :: PULSE_DURATION = 1.d0       ! length of pulse duration (days)
!... decay rate of NMVOC cases short(S), medium(M) and long(L)
  real*8, parameter  :: lambda_5p6 = 1.0/(5.6 * 86400.)
  real*8, parameter  :: lambda_13  = 1.0/(13. * 86400.)
  real*8, parameter  :: lambda_64  = 1.0/(64. * 86400.)
!... local work variables
  integer :: ic, ij, il, ik
  integer, allocatable :: k_pbl(:,:,:)
  real*8 :: jan_pulse_value, apr_pulse_value, jul_pulse_value, oct_pulse_value
  real*8 :: xtmp(72)
!
!     ==========
!     ----------------
!     Begin execution.
!     ----------------
  if (pr_diag) then
    Write (6,*) 'GMI_Waugh_Tracers called by ', procID
  endif
!
!
!... Now do simple NMVOC with 5.6 day loss and different sources
  concentration( 1)%pArray3D(:,:,:) = concentration( 1)%pArray3D(:,:,:) &
                                      * exp(-tdt * lambda_5p6)
  concentration( 4)%pArray3D(:,:,:) = concentration( 4)%pArray3D(:,:,:) &
                                      * exp(-tdt * lambda_5p6)
!... Now do simple NMVOC with 13. day loss and different sources
  concentration( 2)%pArray3D(:,:,:) = concentration( 2)%pArray3D(:,:,:) &
                                      * exp(-tdt * lambda_13)
  concentration( 5)%pArray3D(:,:,:) = concentration( 5)%pArray3D(:,:,:) &
                                      * exp(-tdt * lambda_13)
!... Now do simple NMVOC with 64. day loss and different sources
  concentration( 3)%pArray3D(:,:,:) = concentration( 3)%pArray3D(:,:,:) &
                                      * exp(-tdt * lambda_64)
  concentration( 6)%pArray3D(:,:,:) = concentration( 6)%pArray3D(:,:,:) &
                                      * exp(-tdt * lambda_64)
!... clock_nh tracer 
  concentration( 7)%pArray3D(:,:,:) = concentration( 7)%pArray3D(:,:,:) + tdt/86400.
  do ij=ju1,j2
    if (latdeg(ij) .ge. SLAT_LIM .and. latdeg(ij) .le. NLAT_LIM) then
      concentration( 7)%pArray3D(:,ij,1) = 0.0d0
    endif
  enddo
!
!... 4 pulse tracers
!
!... find source/sink region
  allocate(k_pbl(i1:i2, ju1:j2, k1:k2))
  k_pbl(:,:,:) = 0
  do ij=ju1,j2
    do il=i1,i2
      if (latdeg(ij) .ge. SLAT_LIM .and. latdeg(ij) .le. NLAT_LIM) then
        k_pbl(il,ij,k1) = 1
        xtmp(k1) = gridBoxHeight(il,ij,k1)
        do ik=k1+1,k1+30
          xtmp(ik) = sum(gridBoxHeight(il,ij,k1:ik))
          if(xtmp(ik) .lt. pbl(il,ij)) k_pbl(il,ij,ik) = 1
        enddo
      endif
    enddo
  enddo
!
!... Do ( 8)-"Pulse_Jan_nh",
!... if past TIME_LIM (standard is 30 days) into run set pulse_value to 0.0
  if(gmi_sec .lt. SECPYR .and. month .eq.  1 .and. day .le. PULSE_DURATION) then
    jan_pulse_value = 1.0
  else
    jan_pulse_value = 0.0
  endif
!... set age tracer to "age_value" in lower troposphere
  where (k_pbl(:,:,:) .eq. 1) concentration( 8)%pArray3D(:,:,:) = jan_pulse_value
!
!... Do ( 9)-"Pulse_Apr_nh",
!... if past TIME_LIM (standard is 30 days) into run set pulse_value to 0.0
  if(gmi_sec .lt. SECPYR .and. month .eq.  4 .and. day .le. PULSE_DURATION) then
    apr_pulse_value = 1.0
  else
    apr_pulse_value = 0.0
  endif
!... set age tracer to "age_value" in lower troposphere
  where (k_pbl(:,:,:) .eq. 1) concentration( 9)%pArray3D(:,:,:) = apr_pulse_value
!
!... Do (10)-"Pulse_Jul_nh",
!... if past TIME_LIM (standard is 30 days) into run set pulse_value to 0.0
  if(gmi_sec .lt. SECPYR .and. month .eq.  7 .and. day .le. PULSE_DURATION) then
    jul_pulse_value = 1.0
  else
    jul_pulse_value = 0.0
  endif
!... set age tracer to "age_value" in lower tropospheic tropics
  where (k_pbl(:,:,:) .eq. 1) concentration(10)%pArray3D(:,:,:) = jul_pulse_value
!
!... Do (11)-"Pulse_Oct_nh",
!... if past TIME_LIM (standard is 30 days) into run set pulse_value to 0.0
  if(gmi_sec .lt. SECPYR .and. month .eq. 10 .and. day .le. PULSE_DURATION) then
    oct_pulse_value = 1.0
  else
    oct_pulse_value = 0.0
  endif
!... set age tracer to "age_value" in lower tropospheic tropics
  where (k_pbl(:,:,:) .eq. 1) concentration(11)%pArray3D(:,:,:) = oct_pulse_value
!
  deallocate(k_pbl)
!
  return
!
  end subroutine GMI_Waugh_Tracers
!
!-------------------------------------------------------------------------
!
  end module GmiTaggedCO_AgeOfAir_mod
!
