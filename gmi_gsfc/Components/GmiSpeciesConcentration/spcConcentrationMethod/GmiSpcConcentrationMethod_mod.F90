!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiSpcConcentrationMethod_mod
!
! !INTERFACE:
!
#include "GmiESMF_ErrLog.h"
!
  module GmiSpcConcentrationMethod_mod
!
! !USES:
      use ESMF_Mod
      use GmiESMF_ErrorChecking_mod
      use GmiESMFrcFileReading_mod, only : rcEsmfReadTable, rcEsmfReadLogical
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle, CleanArrayPointer
      use GmiArrayBundlePointer_mod, only : setArrayPointer
      use GmiGrid_mod, only : t_gmiGrid, Get_numSpecies, Get_i1, Get_i2,       &
     &       Get_ju1, Get_j2, Get_k1, Get_k2, Get_i2_gl, Get_j2_gl, Get_ilong, &
     &       Get_ilat, Get_ivert, Get_i1_gl, Get_ju1_gl
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_procID, Get_coscen
      use GmiDiagnosticsMethod_mod  , only : t_Diagnostics, Get_pr_diag
      use GmiPrintError_Mod, only : GmiPrintError
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: InitializeSpcConcentration, FinalizeSpcConcentration
      public  :: ResetFixedConcentration   , constructConcentrationBundle
      public  :: isfixedConcentration, readSpcConcentrationResourceFile
      public  :: Get_concentration      , Set_concentration
      public  :: Get_const_opt         
      public  :: Get_const_init_val     , Get_const_infile_name
      public  :: Get_const_var_name     
      public  :: Get_fixed_const_timpyr , Get_num_fixed_const
      public  :: Get_fixed_const_map    , Get_fixed_const_infile_name
      public  :: Get_fixed_const        , Get_io3_num           
      public  :: Get_tracer_opt         , Get_efol_time           
      public  :: Get_SO3daily_infile_name , Get_SO3monthly_infile_name
      public  :: Get_tr_source_land     , Get_tr_source_ocean        
      public  :: Get_concentrationSurf    , calcConcentrationSurf    
      public  :: Get_concentrationColTrop , calcConcentrationColTrop
      public  :: Get_concentrationColCombo, calcConcentrationColCombo
      public  :: Get_prod, Get_loss, Get_net_cum_tend
      public  :: Set_prod, Set_loss, Set_net_cum_tend
      public  :: Get_num_const_inrecs, Set_num_const_inrecs
!
! !PUBLIC MEMBER DATA:
!
      public  :: t_SpeciesConcentration
! !PRIVATE DATA MEMBERS:
      integer, pointer :: fixedConstMap(:)
      integer          :: numFixedConst
      integer          :: numMoleFrac
      integer          :: numChem

# include "GmiParameters.h"
!# include "gmi_time_constants.h"
# include "gmi_phys_constants.h"
# include "gmi_diag_constants_llnl.h"
#     include "setkin_par.h"

      real*8, parameter :: DUCONST = AVOGAD / (MWTAIR * 2.69d+16 * 10.0d0)

  type t_SpeciesConcentration
    private
    ! For tracer experiments
                                           ! tracer run option
    integer                                :: tracer_opt
                                           ! e-folding time of the tracer (in days)
    real*8                                 :: efol_time
                                           ! land  source of the tracer
    real*8                                 :: tr_source_land
                                           ! ocean source of the tracer
    real*8                                 :: tr_source_ocean
!
                                           ! spec. conc. option
    integer                                :: const_opt                  
                                           ! fixed spc. conc.  times per year
                                           ! (1 => yearly, 12 => monthly)
    integer                                :: fixed_const_timpyr         
            ! number of species to input
    integer                                :: num_const_inrecs
            ! number of fixed constituents
    integer                                :: num_fixed_const
                                           !  mapping of fixed spc. conc. number to spc. conc. number
    integer, pointer                       :: fixed_const_map(:)
                                           ! array of values to initialize each spc. conc. to mixing ratio.
    real*8                , pointer        :: const_init_val(:)
                                           ! Spc. conc. string labels
    character(len=MAX_LENGTH_VAR_NAME)     :: const_var_name    
                                           ! spc. conc. input file name
    character(len=MAX_LENGTH_FILE_NAME)    :: const_infile_name 
                                           ! fixed spc. conc. input file name
    character(len=MAX_LENGTH_FILE_NAME)    :: fixed_const_infile_name
                                           ! index of ozone
    integer                                :: io3_num  
                          ! value of production of a species
    real*8                , pointer        :: prod(:,:,:) => null()
                          ! value of loss       of a species
    real*8                , pointer        :: loss(:,:,:) => null()
    real*8                , pointer        :: fixed_const(:,:,:,:,:) => null()
                                           ! Species concentration of the surface
    real*8                , pointer        :: concentrationSurf    (:,:,:) => null()
                                           ! Column species concentration for troposphere
    real*8                , pointer        :: concentrationColTrop (:,:,:) => null()
                                           ! Column species concentration for strat/trop
    real*8                , pointer        :: concentrationColCombo(:,:,:) => null()
    type(t_GmiArrayBundle), pointer        :: concentration(:) => null()
    type(t_GmiArrayBundle), pointer        :: net_cum_tend(:,:) => null()
    character (len=MAX_LENGTH_FILE_NAME)   :: SO3daily_infile_name
    character (len=MAX_LENGTH_FILE_NAME)   :: SO3monthly_infile_name
  end type t_SpeciesConcentration

! !DESCRIPTION:
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
!
! !IROUTINE: calcConcentrationSurf
!
! !INTERFACE:
!
      subroutine calcConcentrationSurf (self, pr_diag, procID, numSpecies)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: numSpecies
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_SpeciesConcentration), intent(inOut) :: self
!
! !DESCRIPTION:
!  Calculates the surface species concentrations in Dobson units.
!
!  !LOCAL VARIABLES:
      integer :: ic
!EOP
!-------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'calcConcentrationSurf called by ', procID

      do ic = 1, numSpecies
         self%concentrationSurf(:,:,ic)  =  &
     &           self%concentration(ic)%pArray3D(:,:,1) * 1.0d9
      end do
      
      return

      end subroutine calcConcentrationSurf
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calcConcentrationColCombo
!
! !INTERFACE:
!
      subroutine calcConcentrationColCombo (self, mcor, mass, &
     &               pr_diag, procID, i1, i2, ju1, j2, k1, k2, numSpecies)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, numSpecies
      real*8 , intent(in) :: mcor(i1:i2,ju1:j2)
      real*8 , intent(in) :: mass(i1:i2,ju1:j2,k1:k2)
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_SpeciesConcentration), intent(inOut) :: self
!
! !DESCRIPTION:
!  Calculates the column species concentrations when the strat/trop
!  chemical mechanism is selected. Everything is done in Dobson units.
!
!  !LOCAL VARIABLES:
      integer :: ik, ic
!EOP
!-------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'calcConcentrationColCombo called by ', procID

      self%concentrationColCombo(:,:,:) = 0.0d0

      do ic = 1, numSpecies
         do ik = k1, k2
             self%concentrationColCombo(:,:,ic)  =  &
     &           self%concentrationColCombo(:,:,ic)  +  DUCONST*  &
     &           self%concentration(ic)%pArray3D(:,:,ik) * mass(:,:,ik)/mcor(:,:)
         end do
      end do
      
      return

      end subroutine calcConcentrationColCombo
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calcConcentrationColTrop
!
! !INTERFACE:
!
      subroutine calcConcentrationColTrop &
     &               (self, mcor, mass, synoz_threshold, isynoz_num, &
     &                pr_diag, procID, i1, i2, ju1, j2, k1, k2, numSpecies)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID, isynoz_num
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, numSpecies
      real*8 , intent(in) :: synoz_threshold
      real*8 , intent(in) :: mcor(i1:i2,ju1:j2)
      real*8 , intent(in) :: mass(i1:i2,ju1:j2,k1:k2)
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_SpeciesConcentration), intent(inOut) :: self
!
! !DESCRIPTION:
!  Calculates the column species concentrations when the troposphere or strat_trop
!  or strat_trop_aerosol chemical mechanism is selected. Everything is done in Dobson units.
!
!  !LOCAL VARIABLES:
      integer :: ik, ic
!EOP
!-------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'calcConcentrationColTrop called by ', procID

      self%concentrationColTrop(:,:,:) = 0.0d0

      do ic = 1, numSpecies
         do ik = k1, k2
             where(self%concentration(isynoz_num)%pArray3D(:,:,ik) <= synoz_threshold)
                self%concentrationColTrop(:,:,ic)  =  &
     &              self%concentrationColTrop(:,:,ic) + DUCONST*  &
     &              self%concentration(ic)%pArray3D(:,:,ik) * mass(:,:,ik)/mcor(:,:)
             end where
         end do
      end do

      return

      end subroutine calcConcentrationColTrop
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ReadSpcConcentrationResourceFile 
!
! !INTERFACE:
!
      subroutine ReadSpcConcentrationResourceFile (self, gmiGrid, gmiDomain,   &
     &               Diagnostics, config)
!
! !USES:
      use GmiCheckNamelistFile_mod, only : CheckNamelistOptionRange
      use GmiSpeciesRegistry_mod, only : getSpeciesIndex, UNKNOWN_SPECIES
      use GmiStringManipulation_mod, only : constructListNames
!
      implicit none
!
! !INPUT PARAMETERS:
      type (t_gmiGrid    ), intent(in) :: gmiGrid    
      type (t_gmiDomain  ), intent(in) :: gmiDomain  
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config)           , intent(inOut) :: config
      type (t_SpeciesConcentration), intent(inOut) :: self
!
! !DESCRIPTION:
! Reads in Advection related variables from the resource file.
!
! !LOCAL VARIABLES:
      integer :: numSpecies, procID, STATUS, RC, ic, neg_marker
      logical :: pr_diag
      character (len=MAX_LENGTH_SPECIES_NAME) :: tempWord
      character (len=MAX_LENGTH_SPECIES_NAME) :: tempListNames(NSP)
      character (len=MAX_STRING_LENGTH      ) :: fixedConcentrationSpeciesNames
      character(len=ESMF_MAXSTR) :: IAm, err_msg
!EOP
!-------------------------------------------------------------------------
!BOC
      IAm = "ReadSpcConcentrationResourceFile"

      call Get_procID (gmiDomain, procID)
      call Get_pr_diag (Diagnostics, pr_diag)

      if (pr_diag) Write(6,*) IAm, 'called by ', procID
     
      call Get_numSpecies(gmiGrid, numSpecies)

      !################################
      ! Begin reading the resource file
      !################################

      ! --------------------------------------------------------------
      ! const_opt
      !   1:  set const values to const_init_val
      !   2:  read in const values
      !   3:  solid body rotation
      !   4:  dummy test pattern with linear slope in each dimension
      !   5:  exponential in vertical (decays with height)
      !   6:  sin in latitude (largest at equator)
      !   7:  linear vertical gradient
      !   8:  sin in latitude (largest at equator) + vertical gradient
      ! --------------------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%const_opt, &
     &               label   = "const_opt:",&
     &               default = 2, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%const_infile_name, &
     &               label   = "const_infile_name:",&
     &               default = ' ', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%const_var_name, &
     &               label   = "const_var_name:",&
     &               default = 'const', rc=STATUS )
      VERIFY_(STATUS)

      allocate(self%const_init_val(numSpecies))
      self%const_init_val(:)  = 1.0d-30

      call rcEsmfReadTable(config, self%const_init_val, &
     &                     "const_init_val::", rc=STATUS)

      ! sets of fixed consts per year (1 => yearly, 12 => monthly)

      call ESMF_ConfigGetAttribute(config, self%fixed_const_timpyr, &
     &               label   = "fixed_const_timpyr:",&
     &               default = MONTHS_PER_YEAR, rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadTable(config, fixedConcentrationSpeciesNames, &
     &               "fixedConcentrationSpeciesNames::", rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%fixed_const_infile_name, &
     &                label   = "fixed_const_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%io3_num, &
     &                label   = "io3_num:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

!     -----------------------------------------------------
!     tracer_opt
!       0:  no tracer run
!       1:  tracer run with decay of efol_time
!       2:  Tracer suite (Age of Air and "e90" tracer run and 16 others)
!       3:  HTAP Tagged CO tracer run
!       4:  YY Tagged CH4 tracer run
!       5:  BND Tagged CO tracer run
!       6:  BND Tagged CO tracer with additional 5, 10, & 15 day decays also run
!       7:  Waugh tracers
!       8:  e90 and Strat_O3 only
!     -----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%tracer_opt, &
     &                label   = "tracer_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
      if(self%tracer_opt.eq.8) then
!
        call ESMF_ConfigGetAttribute(config, self%SO3daily_infile_name, &
     &                label   = "SO3daily_infile_name:", &
     &                default = 'gmic_MERRA_2004_daily_tracer_input.nc', rc=STATUS )
        VERIFY_(STATUS)
!
        call ESMF_ConfigGetAttribute(config, self%SO3monthly_infile_name, &
     &                label   = "SO3monthly_infile_name:", &
     &                default = 'gmic_MERRA_2004_monthly_tracer_input.nc', rc=STATUS )
        VERIFY_(STATUS)
!
      endif
!
      call ESMF_ConfigGetAttribute(config, self%efol_time, &
     &                label   = "efol_time:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%tr_source_land, &
     &                label   = "tr_source_land:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%tr_source_ocean, &
     &                label   = "tr_source_ocean:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      ! ---------------------------------------------------------------
      ! Check option ranges.  Note that as new options are added, these
      ! range checks will have to be modified.
      ! ---------------------------------------------------------------

      call CheckNamelistOptionRange ('const_opt', self%const_opt, 1, 8)
      call CheckNamelistOptionRange ('tracer_opt', self%tracer_opt, 0, 8)

      if (self%const_opt == 1) then
!       ---------------------------------------------------------
!       Check for any negative const_init_val() marker set in the
!       namelist input file.  If one is found, set all the
!       const_init_val's from the negative value on to the value
!       preceding the negative value.
!       ---------------------------------------------------------
        neg_marker = INT_DUM_VALUE
        do ic = 2, numSpecies
          if ((neg_marker == INT_DUM_VALUE) .and. (self%const_init_val(ic) < 0)) then
            neg_marker = ic
          end if
          if (neg_marker /= INT_DUM_VALUE) then
            self%const_init_val(ic) = self%const_init_val(neg_marker-1)
          end if
        end do
      end if

      self%num_fixed_const = 0

      ! Set the initial value of the list
      tempListNames(:) = ''

      ! Construct the list of names using the long string
      call constructListNames(tempListNames, fixedConcentrationSpeciesNames)

      self%num_fixed_const = count(tempListNames /= '')
      if (self%num_fixed_const > 0) then
         allocate(self%fixed_const_map(self%num_fixed_const))
         self%fixed_const_map(:) = 0

         do ic = 1, self%num_fixed_const
            self%fixed_const_map(ic) = getSpeciesIndex(tempListNames(ic))
         end do
      end if

      !###############
      ! Error Checking
      !###############

      if ((self%fixed_const_timpyr /= 1) .and.  &
     &    (self%fixed_const_timpyr /= MONTHS_PER_YEAR)) then
        err_msg = 'fixed_const_timpyr range problem in Check_Nlvalue.'
        call GmiPrintError  &
     &    (err_msg, .true., 1, self%fixed_const_timpyr, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine ReadSpcConcentrationResourceFile

!EOC
!-------------------------------------------------------------------------
!BOP

      subroutine InitializeSpcConcentration (self, gmiGrid, gmiDomain,         &
     &                     Diagnostics, do_full_chem,            & 
     &                     num_molefrac, num_chem, chem_mecha)
!
! !USES:
      use ReadConcentration_mod     , only : setSliceConcentration
      use ReadConcentration_mod     , only : setSliceConcentrationRst
      use ReadFixedConcentration_mod, only : readFixedConcentration
      use GmiDiagnosticsMethod_mod  , only : Get_rd_restart, Get_restart_inrec,    &
     &     Get_num_tend_outrecs, Get_restart_infile_name, Get_pr_ascii5,       &
     &     Get_pr_tend

!
  implicit none
!
      character (len=*) , intent(in) :: chem_mecha
      integer           , intent(in) :: num_molefrac, num_chem
      logical           , intent(in) :: do_full_chem
      type (t_gmiGrid  ), intent(in) :: gmiGrid 
      type (t_gmiDomain), intent(in) :: gmiDomain 
      type (t_Diagnostics), intent(in) :: Diagnostics 
  
      type (t_SpeciesConcentration)  , intent(inOut) :: self
!
! !LOCAL VARIABLES:
      integer :: ix, ic, id
      integer :: i1, i2, ju1, j2, k1, k2, ilong, ilat, ivert
      integer :: i1_gl, ju1_gl, i2_gl, j2_gl
      integer :: numSpecies, procID
      integer :: restart_inrec, num_tend_outrecs
      logical :: pr_diag, rd_restart, pr_ascii5, pr_tend
      character (len=MAX_LENGTH_FILE_NAME) :: restart_infile_name
      real*8 , allocatable :: coscen(:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID(gmiDomain, procID)
      call Get_pr_diag             (Diagnostics, pr_diag)

      if (pr_diag) Write(6,*) 'InitializeSpcConcentration called by ', procID

      call Get_pr_tend             (Diagnostics, pr_tend)
      call Get_num_tend_outrecs    (Diagnostics, num_tend_outrecs)
      call Get_pr_ascii5           (Diagnostics, pr_ascii5)
      call Get_rd_restart          (Diagnostics, rd_restart )
      call Get_restart_inrec       (Diagnostics, restart_inrec )
      call Get_restart_infile_name (Diagnostics, restart_infile_name )

      call Get_i1    (gmiGrid, i1)
      call Get_i2    (gmiGrid, i2)
      call Get_ju1   (gmiGrid, ju1)
      call Get_j2    (gmiGrid, j2)
      call Get_k1    (gmiGrid, k1)
      call Get_k2    (gmiGrid, k2)
      call Get_ilong (gmiGrid, ilong)
      call Get_ilat  (gmiGrid, ilat )
      call Get_ivert (gmiGrid, ivert)
      call Get_i1_gl (gmiGrid, i1_gl)
      call Get_i2_gl (gmiGrid, i2_gl)
      call Get_j2_gl (gmiGrid, j2_gl)
      call Get_ju1_gl(gmiGrid, ju1_gl)
      call Get_numSpecies(gmiGrid, numSpecies)

      allocate(coscen(ju1_gl:j2_gl))
      call Get_coscen(gmiDomain, coscen)

      call Allocate_concentration(self, i1, i2, ju1, j2, k1, k2, numSpecies)

      call Allocate_concentrationSurf(self, i1, i2, ju1, j2, numSpecies)

      if ((chem_mecha == 'strat_trop') .or. (chem_mecha == 'troposphere')  &
     &         .or. chem_mecha == 'strat_trop_aerosol' ) then
         call Allocate_concentrationColTrop(self, i1, i2, ju1, j2, numSpecies)
         if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then
            call Allocate_concentrationColCombo(self, i1, i2, ju1, j2, numSpecies)
         end if
      end if

      if (self%num_fixed_const > 0) then
         call Allocate_fixed_const    (self, i1, i2, ju1, j2, k1, k2)
      end if

      if (self%num_fixed_const /= 0) then
         call readFixedConcentration  &
               (self%fixed_const, self%fixed_const_infile_name, self%const_var_name, &
                self%num_fixed_const, self%fixed_const_timpyr, pr_diag, procID, &
                ilong, ilat, ivert, i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl)
      end if

      ! Initialization needed for the function isFixedConcentration

      allocate(fixedConstMap(1:self%num_fixed_const))
      fixedConstMap(:) = self%fixed_const_map(1:self%num_fixed_const)
      numFixedConst    = self%num_fixed_const
      numMoleFrac      = num_molefrac
      numChem          = num_chem

      if (pr_ascii5) then
         Allocate (self%loss(2, NUM_OPERATORS, 1:numSpecies))
         Allocate (self%prod(2, NUM_OPERATORS, 1:numSpecies))
         self%prod = 0.0d0
         self%loss = 0.0d0
      end if

      if (pr_tend) then
         call Allocate_net_cum_tend(self, i1, i2, ju1, j2, k1, k2, num_tend_outrecs, &
     &                 NUM_OPERATORS)
      end if

      if (rd_restart) then

          do ic = 1, numSpecies
             call setSliceConcentrationRst &
                    (ic, restart_infile_name, self%const_var_name, restart_inrec, &
                     self%concentration, self%tracer_opt, num_molefrac, &
                     pr_diag, procID, do_full_chem, &
                     i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl, &
                     numSpecies)
          end do

      else

          ! Set the concentration one species at the time
          do ic = 1, numSpecies
             call setSliceConcentration (ic, self%const_infile_name,               &
     &           self%const_var_name, self%const_init_val, self%concentration, &
     &           self%const_opt, self%num_const_inrecs, coscen, num_molefrac,  &
     &           pr_diag, procID, do_full_chem, i1, i2, ju1, j2, k1, k2, i1_gl,&
     &           i2_gl, ju1_gl, j2_gl, numSpecies)
          end do

      end if

      return

      end subroutine InitializeSpcConcentration
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FinalizeSpcConcentration
!
! !INTERFACE:
!
  subroutine FinalizeSpcConcentration (self)
!
  implicit none
!
! !INPUT/OUTPUT PARAMETERS:
  type (t_SpeciesConcentration)   , intent(inOut) :: self
!
! !DESCRIPTION:
! Deallocates species concentration related variables.
! 
!EOP
!-------------------------------------------------------------------------
!BOC
  deallocate(self%concentration)
  deallocate(self%concentrationSurf)
  if (self%num_fixed_const > 0) deallocate(self%fixed_const)

  return

  end subroutine FinalizeSpcConcentration
!EOC
!-------------------------------------------------------------------------
!BOP
  subroutine resetFixedConcentration (self, gmiClock, gmiGrid, numSpecies)
!
! !USES:
!
  use GmiTimeControl_mod, only : GmiSplitDateTime, Get_curGmiDate, t_gmiClock
      implicit none
! !INPUT PARAMETERS:
  integer            , intent(in) :: numSpecies
  type (t_gmiGrid )  , intent(in) :: gmiGrid 
  type (t_gmiClock)  , intent(in) :: gmiClock 
!
! !INPUT/OUTPUT PARAMETERS:
  type (t_SpeciesConcentration), intent(inOut) :: self
!
! !DESCRIPTION:
! This routines fixes values of certain species concentrations.
!
! !LOCAL VARIABLES:
  integer :: i1, i2, ju1, j2, k1, k2, nymd
  integer :: ic, im, icx, idumday, idumyear
!
! !AUTHOR:
! !HISTORY:
!BOP
!-------------------------------------------------------------------------
!BOC
  call Get_i1 (gmiGrid, i1 )
  call Get_i2 (gmiGrid, i2 )
  call Get_ju1(gmiGrid, ju1)
  call Get_j2 (gmiGrid, j2 )
  call Get_k1 (gmiGrid, k1 )
  call Get_k2 (gmiGrid, k2 )

  call Get_curGmiDate(gmiClock, nymd)

  if (self%num_fixed_const /= 0) then
     if (self%fixed_const_timpyr == MONTHS_PER_YEAR) then
        call GmiSplitDateTime(nymd, idumyear, im, idumday)
     else
        im = 1
     end if
     do ic = 1, self%num_fixed_const
        icx = self%fixed_const_map(ic)
        self%concentration(icx)%pArray3D(:,:,:) = self%fixed_const(:,:,:,ic,im)
     end do
  end if

  return

  end subroutine resetFixedConcentration

!EOC
!-------------------------------------------------------------------------
!BOP
  subroutine constructConcentrationBundle &
          (self, concArray, i1, i2, ju1, j2, k1, k2, numSpecies)
  implicit none
! !INPUT PARAMETERS:
  integer                      , intent(in) :: i1, i2, ju1, j2, k1, k2
  integer                      , intent(in) :: numSpecies
  real*8                       , intent(in) :: concArray(i1:i2,ju1:j2,k1:k2,numSpecies)
!
! !INPUT/OUTPUT PARAMETERS:
  type (t_SpeciesConcentration), intent(inOut) :: self
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
  integer :: i
!
!EOP
!-------------------------------------------------------------------------
!BOC
  do i = 1, numSpecies
     call setArrayPointer(self%concentration(i), concArray(i1,ju1,k1,i), &
                 i1 ,i2, ju1, j2, k1, k2)
  end do  
  return
  end subroutine constructConcentrationBundle
!EOC
!-------------------------------------------------------------------------
!BOP
  function isFixedConcentration (ic)
  implicit none
! !INPUT PARAMETERS:
  integer :: ic  ! specie concentration array index
  logical :: isFixedConcentration
! !DESCRIPTION:
! This function checks to see if a particular specie concentration
! has fixed values.
!
! !LOCAL VARIABLES:
  integer :: ix
!
!EOP
!-------------------------------------------------------------------------
!BOC
  isFixedConcentration = .false.
  ixloop: do ix = 1, numFixedConst
     if (fixedConstMap(ix) == ic) then
        isFixedConcentration = .true.
        exit ixloop
     end if
  end do ixloop
  if ((ic > numMoleFrac) .and. (ic <= numChem)) then
!------------------------------------------------------------------
!Take care of the special case of fixed species such as O2, N2, ad.
!------------------------------------------------------------------
     isFixedConcentration = .true.
   end if
   return
  end function isFixedConcentration
!EOC
!-------------------------------------------------------------------------
  subroutine Allocate_concentration (self, i1, i2, ju1, j2, k1, k2, numSpecies)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2, k1, k2, numSpecies
    type (t_SpeciesConcentration), intent(inOut) :: self
    integer :: ic
    Allocate(self%concentration(numSpecies))
    do ic = 1, numSpecies
       Allocate(self%concentration(ic)%pArray3D(i1:i2, ju1:j2, k1:k2))
    end do
    return
  end subroutine Allocate_concentration
!-------------------------------------------------------------------------
  subroutine Get_concentration (self, concentration)
    implicit none
    type (t_GmiArrayBundle), pointer :: concentration (:)
    type (t_SpeciesConcentration), intent(in)   :: self
    concentration => self%concentration
    return
  end subroutine Get_concentration
!-------------------------------------------------------------------------
  subroutine Set_concentration (self, concentration)
    implicit none
    type (t_GmiArrayBundle), pointer :: concentration (:)
    type (t_SpeciesConcentration), intent(inOut) :: self
    self%concentration => concentration
    return
  end subroutine Set_concentration

!-------------------------------------------------------------------------
  subroutine Allocate_net_cum_tend (self, i1, i2, ju1, j2, k1, k2, &
    num_tend_outrecs, NUM_OPERATORS)
    implicit none
    integer, intent(in) :: i1, i2, ju1, j2, k1, k2
    integer, intent(in) :: num_tend_outrecs, NUM_OPERATORS
    type (t_SpeciesConcentration), intent(inOut) :: self
    integer :: ic, it
    Allocate(self%net_cum_tend(num_tend_outrecs, NUM_OPERATORS))
    do it = 1, NUM_OPERATORS
       do ic = 1, num_tend_outrecs
          Allocate(self%net_cum_tend(ic,it)%pArray3D(i1:i2, ju1:j2, k1:k2))
       end do
    end do
    return
  end subroutine Allocate_net_cum_tend
!-------------------------------------------------------------------------
  subroutine Get_net_cum_tend (self, net_cum_tend)
    implicit none
    type (t_GmiArrayBundle), pointer :: net_cum_tend (:,:)
    type (t_SpeciesConcentration), intent(in)   :: self
    net_cum_tend => self%net_cum_tend
    return
  end subroutine Get_net_cum_tend
!-------------------------------------------------------------------------
  subroutine Set_net_cum_tend (self, net_cum_tend)
    implicit none
    type (t_GmiArrayBundle), pointer :: net_cum_tend (:,:)
    type (t_SpeciesConcentration), intent(inOut) :: self
    self%net_cum_tend => net_cum_tend
    return
  end subroutine Set_net_cum_tend
!-------------------------------------------------------------------------
  subroutine Allocate_concentrationSurf (self, i1, i2, ju1, j2, numSpecies)
    implicit none
    integer        , intent(in )  :: i1, i2, ju1, j2, numSpecies
    type (t_SpeciesConcentration), intent(inOut)   :: self
    Allocate(self%concentrationSurf(i1:i2, ju1:j2, numSpecies))
    self%concentrationSurf = 0.0d0
    return
  end subroutine Allocate_concentrationSurf
!-------------------------------------------------------------------------
  subroutine Set_concentrationSurf (self, concentrationSurf)
    implicit none
    real*8         , intent(in)  :: concentrationSurf (:,:,:)
    type (t_SpeciesConcentration), intent(inOut)   :: self
    self%concentrationSurf(:,:,:) = concentrationSurf(:,:,:)
    return
  end subroutine Set_concentrationSurf
!------------------------------------------------------------------------
  subroutine Get_concentrationSurf (self, concentrationSurf)
    implicit none
    real*8         , intent(out)  :: concentrationSurf (:,:,:)
    type (t_SpeciesConcentration), intent(in)   :: self
    concentrationSurf(:,:,:) = self%concentrationSurf(:,:,:)
    return
  end subroutine Get_concentrationSurf
!-------------------------------------------------------------------------
  subroutine Allocate_concentrationColTrop (self, i1, i2, ju1, j2, numSpecies)
    implicit none
    integer        , intent(in )  :: i1, i2, ju1, j2, numSpecies
    type (t_SpeciesConcentration), intent(inOut)   :: self
    Allocate(self%concentrationColTrop(i1:i2, ju1:j2, numSpecies))
    self%concentrationColTrop = 0.0d0
    return
  end subroutine Allocate_concentrationColTrop
!-------------------------------------------------------------------------
  subroutine Set_concentrationColTrop (self, concentrationColTrop)
    implicit none
    real*8         , intent(in)  :: concentrationColTrop (:,:,:)
    type (t_SpeciesConcentration), intent(inOut)   :: self
    self%concentrationColTrop(:,:,:) = concentrationColTrop(:,:,:)
    return
  end subroutine Set_concentrationColTrop
!------------------------------------------------------------------------
  subroutine Get_concentrationColTrop (self, concentrationColTrop)
    implicit none
    real*8         , intent(out)  :: concentrationColTrop (:,:,:)
    type (t_SpeciesConcentration), intent(in)   :: self
    concentrationColTrop(:,:,:) = self%concentrationColTrop(:,:,:)
    return
  end subroutine Get_concentrationColTrop
!-------------------------------------------------------------------------
  subroutine Allocate_concentrationColCombo (self, i1, i2, ju1, j2, numSpecies)
    implicit none
    integer        , intent(in )  :: i1, i2, ju1, j2, numSpecies
    type (t_SpeciesConcentration), intent(inOut)   :: self
    Allocate(self%concentrationColCombo(i1:i2, ju1:j2, numSpecies))
    self%concentrationColCombo = 0.0d0
    return
  end subroutine Allocate_concentrationColCombo
!-------------------------------------------------------------------------
  subroutine Set_concentrationColCombo (self, concentrationColCombo)
    implicit none
    real*8         , intent(in)  :: concentrationColCombo (:,:,:)
    type (t_SpeciesConcentration), intent(inOut)   :: self
    self%concentrationColCombo(:,:,:) = concentrationColCombo(:,:,:)
    return
  end subroutine Set_concentrationColCombo
!------------------------------------------------------------------------
  subroutine Get_concentrationColCombo (self, concentrationColCombo)
    implicit none
    real*8         , intent(out)  :: concentrationColCombo (:,:,:)
    type (t_SpeciesConcentration), intent(in)   :: self
    concentrationColCombo(:,:,:) = self%concentrationColCombo(:,:,:)
    return
  end subroutine Get_concentrationColCombo
!------------------------------------------------------------------------
  subroutine Allocate_fixed_const (self, i1, i2, ju1, j2, k1, k2)
    implicit none
    integer        , intent(in )  :: i1, i2, ju1, j2, k1, k2
    type (t_SpeciesConcentration), intent(inOut)   :: self
    Allocate(self%fixed_const(i1:i2, ju1:j2, k1:k2,  &
                        self%num_fixed_const, self%fixed_const_timpyr))
    self%fixed_const = 0.0d0
    return
  end subroutine Allocate_fixed_const
!-------------------------------------------------------------------------
  subroutine Get_fixed_const (self, fixed_const)
    implicit none
    real*8         , intent(out)  :: fixed_const (:,:,:,:,:)
    type (t_SpeciesConcentration), intent(in)   :: self
    fixed_const(:,:,:,:,:) = self%fixed_const(:,:,:,:,:)
    return
  end subroutine Get_fixed_const
!-------------------------------------------------------------------------
  subroutine Get_const_init_val (self, const_init_val, numSpecies)
    implicit none
    integer        , intent(in )  :: numSpecies
    real*8         , intent(out)  :: const_init_val (:)
    type (t_SpeciesConcentration), intent(in)   :: self
    const_init_val(1:numSpecies) = self%const_init_val(1:numSpecies)
    return
  end subroutine Get_const_init_val
!-------------------------------------------------------------------------
  subroutine Get_num_fixed_const (self, num_fixed_const)
    implicit none
    integer        , intent(out)  :: num_fixed_const 
    type (t_SpeciesConcentration), intent(in)   :: self
    num_fixed_const = self%num_fixed_const
    return
  end subroutine Get_num_fixed_const
!-------------------------------------------------------------------------
  subroutine Get_fixed_const_map (self, fixed_const_map, numSpecies)
    implicit none
    integer        , intent(in )  :: numSpecies
    integer        , intent(out)  :: fixed_const_map (:)
    type (t_SpeciesConcentration), intent(in)   :: self
    fixed_const_map(1:numSpecies) = self%fixed_const_map(1:numSpecies)
    return
  end subroutine Get_fixed_const_map
!-------------------------------------------------------------------------
  subroutine Get_fixed_const_timpyr (self, fixed_const_timpyr)
    implicit none
    integer        , intent(out)  :: fixed_const_timpyr 
    type (t_SpeciesConcentration), intent(in)   :: self
    fixed_const_timpyr = self%fixed_const_timpyr
    return
  end subroutine Get_fixed_const_timpyr
!-------------------------------------------------------------------------
  subroutine Get_const_opt (self, const_opt)
    implicit none
    integer        , intent(out)  :: const_opt 
    type (t_SpeciesConcentration), intent(in)   :: self
    const_opt = self%const_opt
    return
  end subroutine Get_const_opt
!-------------------------------------------------------------------------
  subroutine Get_io3_num (self, io3_num)
    implicit none
    integer        , intent(out)  :: io3_num
    type (t_SpeciesConcentration), intent(in)   :: self
    io3_num = self%io3_num
    return
  end subroutine Get_io3_num
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  subroutine Get_fixed_const_infile_name (self, fixed_const_infile_name)
    implicit none
    character (len=*), intent(out)  :: fixed_const_infile_name
    type (t_SpeciesConcentration)    , intent(in)   :: self
    fixed_const_infile_name = self%fixed_const_infile_name
    return
  end subroutine Get_fixed_const_infile_name
!-------------------------------------------------------------------------
  subroutine Get_const_infile_name (self, const_infile_name)
    implicit none
    character (len=*), intent(out)  :: const_infile_name
    type (t_SpeciesConcentration)    , intent(in)   :: self
    const_infile_name = self%const_infile_name
    return
  end subroutine Get_const_infile_name
!-------------------------------------------------------------------------
  subroutine Get_const_var_name (self, const_var_name)
    implicit none
    character (len=*), intent(out)  :: const_var_name
    type (t_SpeciesConcentration)    , intent(in)   :: self
    const_var_name = self%const_var_name
    return
  end subroutine Get_const_var_name
!-------------------------------------------------------------------------
  subroutine Get_tracer_opt (self, tracer_opt)
    implicit none
    integer        , intent(out)  :: tracer_opt
    type (t_SpeciesConcentration), intent(in)   :: self
    tracer_opt = self%tracer_opt
    return
  end subroutine Get_tracer_opt
!-------------------------------------------------------------------------
  subroutine Get_SO3daily_infile_name (self, SO3daily_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: SO3daily_infile_name
    type (t_SpeciesConcentration), intent(in)   :: self
    SO3daily_infile_name = self%SO3daily_infile_name
    return
  end subroutine Get_SO3daily_infile_name
!-------------------------------------------------------------------------
  subroutine Get_SO3monthly_infile_name (self, SO3monthly_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: SO3monthly_infile_name
    type (t_SpeciesConcentration), intent(in)   :: self
    SO3monthly_infile_name = self%SO3monthly_infile_name
    return
  end subroutine Get_SO3monthly_infile_name
!-------------------------------------------------------------------------
  subroutine Set_num_const_inrecs (self, num_const_inrecs)
    implicit none
    integer        , intent(in)  :: num_const_inrecs
    type (t_SpeciesConcentration), intent(inOut)   :: self
    self%num_const_inrecs = num_const_inrecs
    return
  end subroutine Set_num_const_inrecs
!-------------------------------------------------------------------------
  subroutine Get_num_const_inrecs (self, num_const_inrecs)
    implicit none
    integer        , intent(out)  :: num_const_inrecs
    type (t_SpeciesConcentration), intent(in)   :: self
    num_const_inrecs = self%num_const_inrecs
    return
  end subroutine Get_num_const_inrecs
!-------------------------------------------------------------------------
  subroutine Get_efol_time (self, efol_time)
    implicit none
    real*8         , intent(out)  :: efol_time
    type (t_SpeciesConcentration), intent(in)   :: self
    efol_time = self%efol_time
    return
  end subroutine Get_efol_time
!-------------------------------------------------------------------------
  subroutine Get_tr_source_land (self, tr_source_land)
    implicit none
    real*8         , intent(out)  :: tr_source_land
    type (t_SpeciesConcentration), intent(in)   :: self
    tr_source_land = self%tr_source_land
    return
  end subroutine Get_tr_source_land
!-------------------------------------------------------------------------
  subroutine Get_tr_source_ocean (self, tr_source_ocean)
    implicit none
    real*8         , intent(out)  :: tr_source_ocean
    type (t_SpeciesConcentration), intent(in)   :: self
    tr_source_ocean = self%tr_source_ocean
    return
  end subroutine Get_tr_source_ocean
!-------------------------------------------------------------------------
  subroutine Get_loss (self, loss)
    implicit none
    real*8         , intent(out)  :: loss(:,:,:)
    type (t_SpeciesConcentration), intent(in)   :: self
    loss(:,:,:) = self%loss(:,:,:)
    return
  end subroutine Get_loss
!-------------------------------------------------------------------------
  subroutine Set_loss (self, loss)
    implicit none
    real*8         , intent(in)  :: loss(:,:,:)
    type (t_SpeciesConcentration), intent(inOut)   :: self
    self%loss(:,:,:) = loss(:,:,:)
    return
  end subroutine Set_loss
!-------------------------------------------------------------------------
  subroutine Get_prod (self, prod)
    implicit none
    real*8         , intent(out)  :: prod(:,:,:)
    type (t_SpeciesConcentration), intent(in)   :: self
    prod(:,:,:) = self%prod(:,:,:)
    return
  end subroutine Get_prod
!-------------------------------------------------------------------------
  subroutine Set_prod (self, prod)
    implicit none
    real*8         , intent(in)  :: prod(:,:,:)
    type (t_SpeciesConcentration), intent(inOut)   :: self
    self%prod(:,:,:) = prod(:,:,:)
    return
  end subroutine Set_prod
!-------------------------------------------------------------------------
  end module GmiSpcConcentrationMethod_mod
