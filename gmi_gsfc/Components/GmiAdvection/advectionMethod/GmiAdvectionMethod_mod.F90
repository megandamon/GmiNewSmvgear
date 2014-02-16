!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiAdvectionMethod_mod
!
! !INTERFACE:
!
#include "GmiESMF_ErrLog.h"
!
      module GmiAdvectionMethod_mod
!
! USES:
      use ESMF_Mod
      use GmiESMF_ErrorChecking_mod
      use GmiESMFrcFileReading_mod, only : rcEsmfReadTable, rcEsmfReadLogical
      use GmiGrid_mod          , only : t_gmiGrid, Get_gmi_nborder
      use GmiGrid_mod          , only : Get_i1, Get_i2, Get_ju1, Get_j2, Get_k1, Get_k2
      use GmiGrid_mod          , only : Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl
      use GmiGrid_mod          , only : Get_ilo, Get_ihi, Get_julo, Get_jhi
      use GmiGrid_mod          , only : Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jhi_gl
      use GmiGrid_mod          , only : Get_ilong, Get_ivert
      use GmiGrid_mod          , only : Get_numSpecies
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_communicatorWorld
      use GmiDomainDecomposition_mod, only : Get_procID
      use GmiDomainDecomposition_mod, only : Get_communicatorSouthPole, Get_eastWestFlag
      use GmiDomainDecomposition_mod, only : Get_communicatorNorthPole, Get_northSouthPole
      use GmiDomainDecomposition_mod, only : Get_mapi_all
      use GmiDomainDecomposition_mod, only : Get_southBorder, Get_westBorder
      use GmiDomainDecomposition_mod, only : Get_northBorder, Get_eastBorder
      use GmiDomainDecomposition_mod, only : Get_southNeighbor, Get_westNeighbor
      use GmiDomainDecomposition_mod, only : Get_northNeighbor, Get_eastNeighbor
      use GmiDomainDecomposition_mod, only : Get_numDomains
      use GmiDomainDecomposition_mod, only : Get_geofac, Get_geofac_pc,        &
     &       Get_rel_area, Get_mcor, Get_cose
      use GmiDomainDecomposition_mod, only : Get_numLatDomains, Get_numLonDomains
      use GmiDiagnosticsMethod_mod  , only : t_Diagnostics, Get_pr_diag,       &
     &       Get_do_mean, Get_pr_flux,Get_pr_psf_flux, Get_pr_const_flux,      &
     &       Get_flux_species, Get_flux_species_num, Set_flux_species
      use GmiPrintError_mod, only : GmiPrintError
      use GmiMetFieldsControl_mod, only : t_metFields, Get_pt, Get_ai, Get_bi, &
     &       Get_dap, Get_dbk, Get_pctm1, Get_pctm2, Get_xmass, Get_ymass
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: initializeAdvection, RunAdvection, resetFlux
      public  :: Get_flux_x, Set_flux_x, Get_j1p, Get_j2p
      public  :: Get_flux_y, Set_flux_y
      public  :: Get_flux_z, Set_flux_z
      public  :: Get_psf_flux, Set_psf_flux
      public  :: Get_count_flux, Set_count_flux, Get_press_fix_opt
      public  :: Get_air_mass_flux, Set_air_mass_flux
      public  :: Get_numAdvectedSpecies, Get_advectedSpecies, Get_advec_opt
      public  :: Get_advectedSpeciesMap, Get_do_var_adv_tstp
      public  :: Get_num_adv_time_steps, Set_num_adv_time_steps
      public  :: Get_do_grav_set, Get_advec_consrv_opt, Get_pmet2_opt
!
! !PUBLIC MEMEMBER DATA:
      public  :: t_Advection
!
#     include "GmiParameters.h"
#     include "setkin_par.h"
!
      type t_Advection
         integer  :: advec_opt          ! advection          option
         integer  :: trans_opt
         integer  :: press_fix_opt      ! pressure fixer     option
         integer  :: pmet2_opt          ! pmet2              option
         integer  :: advec_consrv_opt   ! advection conserve option
         integer  :: j1p                ! determines size of the Polar cap
         integer  :: j2p                ! j2_gl - j1p + 1
         logical  :: do_grav_set        ! do gravitational settling of aerosols?
         logical  :: do_var_adv_tstp    ! take variable advection time steps?
         integer  :: num_adv_time_steps ! number of advection time steps
         integer  :: numAdvectedSpecies ! number of advected species
         integer  :: advec_flag(NSP)
         integer,                           pointer :: advectedSpeciesMap(:)
         character (len=MAX_STRING_LENGTH), pointer :: advectedSpecies(:)
!
         real*8, pointer :: air_mass_flux (:,:,:,:) => null()
         real*8, pointer :: flux_x        (:,:,:,:) => null()
         real*8, pointer :: flux_y        (:,:,:,:) => null()
         real*8, pointer :: flux_z        (:,:,:,:) => null()
         real*8, pointer :: psf_flux      (:,:)     => null()
         integer :: count_flux
      end type t_Advection
!
! !DESCRIPTION:
!  Interface routines for the advection component.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  12Jan2007 Initial code.
!
!EOP
!-------------------------------------------------------------------------
       CONTAINS
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initializeAdvection
!
! !INTERFACE:
!
      subroutine initializeAdvection(self, gmiGrid, gmiDomain, Diagnostics, &
     &                      config, speciesNames)
!
      implicit none
!
! !INPUT PARAMETERS:
      type (t_gmiGrid    ), intent(in) :: gmiGrid
      type (t_gmiDomain  ), intent(in) :: gmiDomain
      character (len=MAX_LENGTH_SPECIES_NAME), intent(in) :: speciesNames(:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Diagnostics), intent(inOut) :: Diagnostics
      type (ESMF_Config), intent(inOut) :: config
      type (t_Advection), intent(inOut) :: self
!
! !DESCRIPTION:
! Initialize the Advection component.
!
! !LOCAL VARIABLES:
      integer :: i1, i2, ju1, j2, k1, k2, flux_species_num
      logical :: pr_flux, pr_const_flux, pr_psf_flux
!EOP
!-------------------------------------------------------------------------
!BOC
      call readAdvectionResourceFile(self, gmiGrid, gmiDomain, Diagnostics, &
     &                      config, speciesNames)
!
      call Get_i1  (gmiGrid, i1)
      call Get_i2  (gmiGrid, i2)
      call Get_ju1 (gmiGrid, ju1)
      call Get_j2  (gmiGrid, j2)
      call Get_k1  (gmiGrid, k1)
      call Get_k2  (gmiGrid, k2)
!
      call Get_pr_flux (Diagnostics, pr_flux)
!
      if (pr_flux) then
!
         self%count_flux = 0
!
         Allocate (self%air_mass_flux(i1:i2, ju1:j2, k1:k2, 3))
         self%air_mass_flux = 0.0d0
!
         call Get_pr_const_flux (Diagnostics, pr_const_flux)
!
         if (pr_const_flux) then
            call Get_flux_species_num (Diagnostics, flux_species_num)
!
            Allocate (self%flux_x(i1:i2, ju1:j2, k1:k2, flux_species_num))
            self%flux_x = 0.0d0
!
            Allocate (self%flux_y(i1:i2, ju1:j2, k1:k2, flux_species_num))
            self%flux_y = 0.0d0
!
            Allocate (self%flux_z(i1:i2, ju1:j2, k1:k2, flux_species_num))
            self%flux_z = 0.0d0
         endif
!
         call Get_pr_psf_flux (Diagnostics, pr_psf_flux)
!
         if (pr_psf_flux) then
            Allocate (self%psf_flux (i1:i2, ju1:j2))
            self%psf_flux = 0.0d0
         endif
      endif
!
      return
!
      end subroutine initializeAdvection
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readAdvectionResourceFile
!
! !INTERFACE:
!
      subroutine readAdvectionResourceFile(self, gmiGrid, gmiDomain, Diagnostics, &
     &                      config, speciesNames)
!
! !USES:
      use GmiCheckNamelistFile_mod, only : CheckNamelistOptionRange
      use GmiSpeciesRegistry_mod, only : getSpeciesIndex
      use GmiStringManipulation_mod, only : constructListNames
!
      implicit none
!
! !INPUT PARAMETERS:
      type (t_gmiGrid    ), intent(in) :: gmiGrid
      type (t_gmiDomain  ), intent(in) :: gmiDomain
      character (len=MAX_LENGTH_SPECIES_NAME), intent(in) :: speciesNames(:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Diagnostics), intent(inOut) :: Diagnostics
      type (ESMF_Config), intent(inOut) :: config
      type (t_Advection), intent(inOut) :: self
!
! !DESCRIPTION:
! Reads in Advection related variables from the resource file.
!
! !LOCAL VARIABLES:
      integer :: numSpecies, procID, j2_gl, STATUS, RC, ic, ica, ju1_gl
      logical :: pr_diag, defFalse, endTable, pr_flux
      integer, allocatable :: flux_species(:)
      real*8  :: mdt, tdt
      character :: cValue
      character (len=MAX_LENGTH_SPECIES_NAME) :: tempWord
      character (len=MAX_LENGTH_SPECIES_NAME) :: tempListNames(NSP)
      character (len=MAX_STRING_LENGTH      ) :: advectedSpeciesNames
      character(len=ESMF_MAXSTR) :: IAm, err_msg
!EOP
!-------------------------------------------------------------------------
!BOC
      IAm = "readAdvectionResourceFile"
!
      call Get_procID (gmiDomain, procID)
      call Get_pr_diag (Diagnostics, pr_diag)
!
      if (pr_diag) Write(6,*) IAm, 'called by ', procID
!
      call Get_numSpecies(gmiGrid, numSpecies)
!
      call CheckNamelistOptionRange ('numSpecies       ', numSpecies,        0, NSP)
!
!
      !################################
      ! Begin reading the resource file
      !################################
!
      ! -----------------------------------------------------
      ! trans_opt
      !   1:  do LLNLTRANS transport
      ! -----------------------------------------------------
!
      call ESMF_ConfigGetAttribute(config, self%trans_opt, label="trans_opt:",&
     &               default=1, rc=STATUS )
      VERIFY_(STATUS)
!
      ! -----------------------
      ! advec_opt
      !   0:  no advection
      !   1:  do DAO2 advection
      !   2:  do advecCore advection
      ! -----------------------
!
      call ESMF_ConfigGetAttribute(config, self%advec_opt, label="advec_opt:",&
     &               default=1, rc=STATUS )
      VERIFY_(STATUS)
!
      ! ------------------------------
      ! press_fix_opt
      !   0:  no   pressure fixer used
      !   1:  LLNL pressure fixer used
      !   2:  UCI  pressure fixer used
      ! ------------------------------
!
      call ESMF_ConfigGetAttribute(config, self%press_fix_opt, &
     &               label="press_fix_opt:", default=1, rc=STATUS )
      VERIFY_(STATUS)
!
      ! ------------------------------------------------------------
      ! pmet2_opt
      !   0:  use pmet2
      !   1:  use (pmet2 - "global mean change in surface pressure")
      ! ------------------------------------------------------------
!
      call ESMF_ConfigGetAttribute(config, self%pmet2_opt, label="pmet2_opt:", &
     &               default=1, rc=STATUS )
      VERIFY_(STATUS)
!
      ! --------------------------------------------------
      ! advec_consrv_opt
      !   0:  conserve tracer conc.;             use pmet2
      !   1:  conserve tracer mass;              use pmet2
      !   2:  conserve both tracer conc. & mass; use pctm2
      ! --------------------------------------------------
!
      call ESMF_ConfigGetAttribute(config, self%advec_consrv_opt, &
     &                label="advec_consrv_opt:", default=2, rc=STATUS )
      VERIFY_(STATUS)
!
      !--------------------------------------------------
      ! Get the list of advected species as a long string
      !--------------------------------------------------
!
      call rcEsmfReadTable(config, advectedSpeciesNames, &
     &                     "advectedSpeciesNames::", rc=STATUS)
!
      !---------------------------------------
      ! do gravitational settling of aerosols?
      !---------------------------------------
!
      call rcEsmfReadLogical(config, self%do_grav_set, "do_grav_set:", &
     &                       default=.false., rc=STATUS)
!
      !------------------------------------
      ! Take variable advection time steps?
      !------------------------------------
!
      call rcEsmfReadLogical(config, self%do_var_adv_tstp, "do_var_adv_tstp:", &
     &                       default=.false., rc=STATUS)
!
      !##############################
      ! End reading the resource file
      !##############################
!
      !----------------------------------------------------------------
      ! Check option ranges.  Note that as new options are added, these
      ! range checks will have to be modified.
      !----------------------------------------------------------------
!
      call CheckNamelistOptionRange ('advec_opt       ', self%advec_opt,        0, 2)
      call CheckNamelistOptionRange ('advec_consrv_opt', self%advec_consrv_opt, 0, 2)
      call CheckNamelistOptionRange ('pmet2_opt       ', self%pmet2_opt,        0, 1)
      call CheckNamelistOptionRange ('press_fix_opt   ', self%press_fix_opt,    0, 2)
!
      ! ---------------------------------------------------------------
      ! advec_flag(): an array of flags that indicate whether or not to
      !               advect a given species; if not explicitly set for
      !               a particular species, advec_flag_default is used;
      !   0:  no advection
      !   1:  advection
      ! ---------------------------------------------------------------
!
      self%advec_flag(:)        = 0
!
      ! -------------------------------------
      ! j1p: determines size of the Polar cap
      ! -------------------------------------
!
      call Get_ju1_gl(gmiGrid, ju1_gl)
      call Get_j2_gl (gmiGrid, j2_gl)
!
!... read in polar cap size from resource file - default=3
!      self%j1p     = 3
      call ESMF_ConfigGetAttribute(config, self%j1p, &
     &                label="j1p:", default=3, rc=STATUS )
!
      self%j2p     = j2_gl - self%j1p + 1
!
      !-------------------------------------------------------------
      ! Construct the list of advected species using the long string
      !-------------------------------------------------------------
!
      ! Set the initial value of the list
      tempListNames(:) = ''
!
      call constructListNames(tempListNames, advectedSpeciesNames)
!
      self%numAdvectedSpecies = Count (tempListNames(:) /= '')
!
      ! For the case all the species are advected
!
      if (self%numAdvectedSpecies == 0) then
         self%numAdvectedSpecies = numSpecies
         tempListNames(1:self%numAdvectedSpecies) = &
     &                speciesNames(1:self%numAdvectedSpecies)
      end if
!
!
      ! Only consider advected species.
!
      allocate(self%advectedSpecies(self%numAdvectedSpecies))
      allocate(self%advectedSpeciesMap(self%numAdvectedSpecies))
!
      do ic = 1, self%numAdvectedSpecies
         ica                          = getSpeciesIndex(tempListNames(ic))
         self%advec_flag(ica)         = 1
         self%advectedSpecies(ic)     = tempListNames(ic)
         self%advectedSpeciesMap(ic ) = ica
      end do
!
      call Get_pr_flux(Diagnostics, pr_flux)
!
      if (pr_flux) then
         allocate(flux_species(numSpecies))
         call Get_flux_species(Diagnostics, flux_species)
!
         do ic = 1, numSpecies
            if ((flux_species(ic) == 1) .and. (self%advec_flag(ic) == 0)) then
               print *,'*******************************************'
               print *,'NOT saving advective flux diag for species ',ic &
     &                ,'because it is a non-advected species'
               print *,'*******************************************'
               flux_species(ic) = 0
           endif
         enddo
!
         call Set_flux_species(Diagnostics, flux_species)
      end if
!
      !###############
      ! Error Checking
      !###############
!
!
      if (pr_flux .and. (self%numAdvectedSpecies == 0)) then
        err_msg = 'pr_flux/numAdvectedSpecies problem in the rc File.'
        call GmiPrintError (err_msg, .true., 1, self%numAdvectedSpecies, &
     &          0, 0, 0.0d0, 0.0d0)
      end if
!
      call ESMF_ConfigGetAttribute(config, mdt, label="mdt:", &
     &               default=21600.0d0, rc=STATUS )
!
      call ESMF_ConfigGetAttribute(config, tdt, label="DT:", &
     &               default=180.0d0, rc=STATUS )
!
      if ((self%advec_opt /= 0) .and. (Mod (Nint (mdt), Nint(tdt)) /= 0)) then
        err_msg = 'mdt is not evenly divisible by tdt in Check_Nlvalue.'
        call GmiPrintError (err_msg, .true., 1, Nint(tdt), 0, 1, mdt, 0.0d0)
      end if
!
      if ((self%press_fix_opt == 0) .and. (self%advec_consrv_opt == 2)) then
        err_msg =  &
     &    'press_fix_opt/advec_consrv_opt problem in Check_Nlvalue.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, self%press_fix_opt, self%advec_consrv_opt,  &
     &     0, 0.0d0, 0.0d0)
      end if
!
!.sds      if (self%j1p == ju1_gl+1) then  ! Polar Cap NOT Enlarged.
!.sds!       --------------------------------------------------------
!.sds!       It is not known whether this option currently works or
!.sds!       not.  If you want to try it, just comment out this error
!.sds!       message, re-compile, and re-link.
!.sds!       --------------------------------------------------------
!.sds        err_msg = 'j1p problem in Check_Nlvalue.'
!.sds        call GmiPrintError (err_msg, .true., 1, self%j1p, 0, 0, 0.0d0, 0.0d0)
!.sds      end if
!
      return
!
      end subroutine readAdvectionResourceFile
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: RunAdvection
!
! !INTERFACE:
!
      subroutine RunAdvection  &
     &  (self, gmiGrid, gmiDomain, Diagnostics, metFields, &
     &   do_chem_grp, concentration, mw)
!
! !USES:
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_gmiDomain), intent(in   ) :: gmiDomain
      type(t_gmiGrid  ), intent(in   ) :: gmiGrid
      real*8 , intent(in) :: mw (1:)
      logical, intent(in) :: do_chem_grp
      type(t_metFields), intent(in) :: metFields
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_Advection), intent(inOut) :: self
      type(t_Diagnostics), intent(inOut) :: Diagnostics
      type (t_GmiArrayBundle), intent(inOut) :: concentration(:)
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
      integer :: ilong, ivert, ic, numSpecies
      integer :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer :: ilo, ihi, julo, jhi, ilo_gl, ihi_gl, julo_gl, jhi_gl
      integer :: i1, i2, ju1, j2, k1, k2
      integer :: numDomains, numLatDomains, numLonDomains
      integer :: communicatorWorld
      integer :: communicatorNorthPole, communicatorSouthPole
      integer nb, sb, eb, wb
      integer, allocatable :: ewflag(:)
      integer, allocatable :: nspoleu(:)
      integer, allocatable :: mapi_all(:,:)
      integer north_domain, south_domain, east_domain, west_domain
      integer :: procID
      integer :: gmi_nborder
      logical :: do_mean, pr_diag
      logical :: pr_flux, pr_const_flux, pr_psf_flux
      real*8,  allocatable :: const(:,:,:,:)
      real*8               :: geofac_pc
      real*8 , allocatable :: geofac(:)
      real*8 , allocatable :: cose(:)
      real*8 , allocatable :: mcor(:,:), rel_area(:,:)
      integer :: advec_consrv_opt
      real*8 , allocatable :: dap(:)
      real*8 , allocatable :: dbk(:)
      real*8               :: pt
      real*8 , allocatable :: ai(:)
      real*8 , allocatable :: bi(:)
      integer              :: flux_species_num
      integer, allocatable :: flux_species (:)
      real*8 , allocatable :: pctm1(:,:)
      real*8 , allocatable :: pctm2(:,:)
      real*8 , allocatable :: xmass(:,:,:)
      real*8 , allocatable :: ymass(:,:,:)
!
!EOP
!-------------------------------------------------------------------------------
!BOC
      call Get_procID       (gmiDomain, procID    )
      call Get_pr_diag(Diagnostics, pr_diag)
!
      if (pr_diag) Write(6,*) 'runAdvection called by ', procID
!
      ! Get the GMI grid information
      call Get_i1       (gmiGrid, i1   )
      call Get_i2       (gmiGrid, i2   )
      call Get_ju1      (gmiGrid, ju1  )
      call Get_j2       (gmiGrid, j2   )
      call Get_k1       (gmiGrid, k1   )
      call Get_k2       (gmiGrid, k2   )
      call Get_i1_gl    (gmiGrid, i1_gl)
      call Get_i2_gl    (gmiGrid, i2_gl)
      call Get_ju1_gl   (gmiGrid, ju1_gl)
      call Get_j2_gl    (gmiGrid, j2_gl )
      call Get_ilo      (gmiGrid, ilo  )
      call Get_ihi      (gmiGrid, ihi  )
      call Get_julo     (gmiGrid, julo )
      call Get_jhi      (gmiGrid, jhi  )
      call Get_ilo_gl   (gmiGrid, ilo_gl)
      call Get_ihi_gl   (gmiGrid, ihi_gl)
      call Get_julo_gl  (gmiGrid, julo_gl)
      call Get_jhi_gl   (gmiGrid, jhi_gl)
      call Get_numSpecies (gmiGrid, numSpecies)
      call Get_gmi_nborder (gmiGrid, gmi_nborder)
!
      ilong = i2 - i1 + 1
      ivert = k2 - k1 + 1
!
      call Get_do_mean(Diagnostics, do_mean)
      call Get_pr_flux(Diagnostics, pr_flux)
      call Get_pr_const_flux(Diagnostics, pr_const_flux)
      call Get_pr_psf_flux(Diagnostics, pr_psf_flux)
!
      call Get_flux_species_num(Diagnostics, flux_species_num)
      if (flux_species_num > 0) then
         allocate(flux_species(numSpecies))
         call Get_flux_species    (Diagnostics, flux_species    )
      end if
!
      call Get_numDomains   (gmiDomain, numDomains)
      call Get_numLonDomains(gmiDomain, numLonDomains)
      call Get_numLatDomains(gmiDomain, numLatDomains)
!
      call Get_pt (metFields, pt)
!
      allocate(ai(k1-1:k2))
      call Get_ai (metFields, ai)
!
      allocate(bi(k1-1:k2))
      call Get_bi (metFields, bi)
!
      allocate(dap(k1:k2))
      call Get_dap (metFields, dap)
!
      allocate(dbk(k1:k2))
      call Get_dbk (metFields, dbk)
!
      allocate(pctm1(ilo:ihi,julo:jhi))
      allocate(pctm2(ilo:ihi,julo:jhi))
      allocate(xmass(ilo:ihi,julo:jhi,k1:k2))
      allocate(ymass(ilo:ihi,julo:jhi,k1:k2))
!
      call Get_pctm1 (metFields, pctm1)
      call Get_pctm2 (metFields, pctm2)
      call Get_xmass (metFields, xmass)
      call Get_ymass (metFields, ymass)
!
      allocate(cose  (ju1_gl:j2_gl+1))
      allocate(mcor  (i1:i2,ju1:j2))
      allocate(rel_area(i1:i2,ju1:j2))
      allocate(geofac(ju1_gl:j2_gl))
!
      call Get_mcor     (gmiDomain, mcor     )
      call Get_cose     (gmiDomain, cose     )
      call Get_geofac   (gmiDomain, geofac   )
      call Get_rel_area (gmiDomain, rel_area )
      call Get_geofac_pc(gmiDomain, geofac_pc)
!
      ! Get the communicators
      call Get_communicatorWorld    (gmiDomain, communicatorWorld    )
      call Get_communicatorSouthPole(gmiDomain, communicatorSouthPole)
      call Get_communicatorNorthPole(gmiDomain, communicatorNorthPole)
!
      allocate(ewflag(ilo_gl:ihi_gl))
      call Get_eastWestFlag(gmiDomain, ewflag)
!
      allocate(nspoleu(julo_gl:jhi_gl))
      call Get_northSouthPole(gmiDomain, nspoleu)
!
      allocate(mapi_all(2,numDomains))
      call Get_mapi_all(gmiDomain, mapi_all)
!
      call Get_westBorder (gmiDomain, wb)
      call Get_eastBorder (gmiDomain, eb)
      call Get_northBorder(gmiDomain, nb)
      call Get_southBorder(gmiDomain, sb)
!
      call Get_westNeighbor (gmiDomain, west_domain)
      call Get_eastNeighbor (gmiDomain, east_domain)
      call Get_northNeighbor(gmiDomain, north_domain)
      call Get_southNeighbor(gmiDomain, south_domain)
!
      allocate(const(i1:i2, ju1:j2, k1:k2, numSpecies))
!
      do ic = 1, numSpecies
         const(:,:,:,ic) = concentration(ic)%pArray3D(:,:,:)
      end do
!
      call  Update_Advec  &
     &  (do_chem_grp, self%num_adv_time_steps, self%advec_consrv_opt, &
     &   self%advec_flag, pt, geofac_pc, geofac, cose, mcor, ai, bi, &
     &   rel_area, pctm1, pctm2, dap, dbk, xmass, ymass, const, &
     &   self%air_mass_flux, self%flux_x, self%flux_y, self%flux_z, self%psf_flux, &
     &   do_mean, pr_flux, pr_psf_flux, self%count_flux, &
     &   pr_diag, procID, numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   flux_species_num, flux_species, mw, pr_const_flux, &
     &   self%j1p, self%j2p, i1_gl, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
     &   i1, i2, ju1, j2, k1, k2, ilong, ivert, numSpecies, &
     &   communicatorNorthPole, communicatorSouthPole, &
     &   ewflag, nspoleu, nb, sb, eb, wb,  &
     &   north_domain, south_domain, east_domain, west_domain,  &
     &   communicatorWorld, mapi_all)
!
      do ic = 1, numSpecies
         concentration(ic)%pArray3D(:,:,:) = const(:,:,:,ic)
      end do
!
      deallocate(const)
      deallocate(ewflag)
      deallocate(nspoleu)
      deallocate(mapi_all)
      deallocate(ai)
      deallocate(bi)
      deallocate(dap)
      deallocate(dbk)
      deallocate(cose  )
      deallocate(mcor  )
      deallocate(rel_area)
      deallocate(geofac)
      deallocate(pctm2)
      deallocate(pctm1)
      deallocate(xmass)
      deallocate(ymass)
!
      if (flux_species_num > 0) then
         deallocate(flux_species)
      end if
!
      return
!
      end subroutine RunAdvection
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: resetFlux
!
! !INTERFACE:
!
      subroutine resetFlux (self, Diagnostics)
!
      implicit none
!
      type(t_Diagnostics), intent(in) :: Diagnostics
      type(t_Advection  ), intent(inOut) :: self
!
! !DESCRIPTION:
! Resets the flux diagnostics.
!
! !LOCAL VARIABLES:
      logical, save :: first = .true.
      logical :: pr_flux, pr_const_flux, pr_psf_flux
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_pr_flux      (Diagnostics, pr_flux)
      call Get_pr_psf_flux  (Diagnostics, pr_psf_flux )
      call Get_pr_const_flux(Diagnostics, pr_const_flux)
!
!... reset accumulated quantities
      if (pr_flux) then
!
!... zero air mass flux
         self%air_mass_flux = 0.0d0
!
!... reset counter
         self%count_flux = 0
!
!... constituent component fluxes
         if (pr_const_flux) then
            self%flux_x = 0.0d0
            self%flux_y = 0.0d0
            self%flux_z = 0.0d0
         endif
!
!... constituent component fluxes
         if (pr_psf_flux) then
            self%psf_flux = 0.0d0
         endif
!
       endif
!
      return
!
      end subroutine resetFlux
!EOC
!-----------------------------------------------------------------------------
      subroutine Get_advectedSpeciesMap(self, advectedSpeciesMap)
      implicit none
      integer, intent(out) :: advectedSpeciesMap(:)
      type(t_Advection), intent(in) :: self
      advectedSpeciesMap(:) = self%advectedSpeciesMap(:)
      return
      end subroutine Get_advectedSpeciesMap
!-------------------------------------------------------------------------------
      subroutine Get_numAdvectedSpecies(self, numAdvectedSpecies)
      implicit none
      integer, intent(out) :: numAdvectedSpecies
      type(t_Advection), intent(in) :: self
      numAdvectedSpecies = self%numAdvectedSpecies
      return
      end subroutine Get_numAdvectedSpecies
!-------------------------------------------------------------------------------
      subroutine Get_advectedSpecies(self, advectedSpecies)
      implicit none
      character (len=MAX_STRING_LENGTH), intent(out) :: advectedSpecies(:)
      type(t_Advection), intent(in) :: self
      advectedSpecies(:) = self%advectedSpecies(:)
      return
      end subroutine Get_advectedSpecies
!-------------------------------------------------------------------------------
      subroutine Get_advec_opt(self, advec_opt)
      implicit none
      integer, intent(out) :: advec_opt
      type(t_Advection), intent(in) :: self
      advec_opt = self%advec_opt
      return
      end subroutine Get_advec_opt
!-------------------------------------------------------------------------------
      subroutine Get_advec_consrv_opt(self, advec_consrv_opt)
      implicit none
      integer, intent(out) :: advec_consrv_opt
      type(t_Advection), intent(in) :: self
      advec_consrv_opt = self%advec_consrv_opt
      return
      end subroutine Get_advec_consrv_opt
!-------------------------------------------------------------------------------
      subroutine Get_pmet2_opt(self, pmet2_opt)
      implicit none
      integer, intent(out) :: pmet2_opt
      type(t_Advection), intent(in) :: self
      pmet2_opt = self%pmet2_opt
      return
      end subroutine Get_pmet2_opt
!-------------------------------------------------------------------------------
      subroutine Get_j1p(self, j1p)
      implicit none
      integer, intent(out) :: j1p
      type(t_Advection), intent(in) :: self
      j1p = self%j1p
      return
      end subroutine Get_j1p
!-------------------------------------------------------------------------------
      subroutine Get_j2p(self, j2p)
      implicit none
      integer, intent(out) :: j2p
      type(t_Advection), intent(in) :: self
      j2p = self%j2p
      return
      end subroutine Get_j2p
!-------------------------------------------------------------------------------
      subroutine Get_press_fix_opt(self, press_fix_opt)
      implicit none
      integer, intent(out) :: press_fix_opt
      type(t_Advection), intent(in) :: self
      press_fix_opt = self%press_fix_opt
      return
      end subroutine Get_press_fix_opt
!-------------------------------------------------------------------------------
      subroutine Get_num_adv_time_steps(self, num_adv_time_steps)
      implicit none
      integer, intent(out) :: num_adv_time_steps
      type(t_Advection), intent(in) :: self
      num_adv_time_steps = self%num_adv_time_steps
      return
      end subroutine Get_num_adv_time_steps
!-------------------------------------------------------------------------------
      subroutine Set_num_adv_time_steps(self, num_adv_time_steps)
      implicit none
      integer, intent(in) :: num_adv_time_steps
      type(t_Advection), intent(inOut) :: self
      self%num_adv_time_steps = num_adv_time_steps
      return
      end subroutine Set_num_adv_time_steps
!-------------------------------------------------------------------------------
      subroutine Get_do_grav_set(self, do_grav_set)
      implicit none
      logical, intent(out) :: do_grav_set
      type(t_Advection), intent(in) :: self
      do_grav_set = self%do_grav_set
      return
      end subroutine Get_do_grav_set
!-------------------------------------------------------------------------------
      subroutine Get_do_var_adv_tstp(self, do_var_adv_tstp)
      implicit none
      logical, intent(out) :: do_var_adv_tstp
      type(t_Advection), intent(in) :: self
      do_var_adv_tstp = self%do_var_adv_tstp
      return
      end subroutine Get_do_var_adv_tstp
!-------------------------------------------------------------------------------
      subroutine Get_count_flux(self, count_flux)
      implicit none
      integer, intent(out) :: count_flux
      type(t_Advection), intent(in) :: self
      count_flux = self%count_flux
      return
      end subroutine Get_count_flux
!-------------------------------------------------------------------------------
      subroutine Set_count_flux(self, count_flux)
      implicit none
      integer, intent(in) :: count_flux
      type(t_Advection), intent(inOut) :: self
      self%count_flux = count_flux
      return
      end subroutine Set_count_flux
!-------------------------------------------------------------------------------
      subroutine Get_flux_x(self, flux_x)
      implicit none
      real*8, intent(out) :: flux_x(:,:,:,:)
      type(t_Advection), intent(in) :: self
      flux_x(:,:,:,:) = self%flux_x(:,:,:,:)
      return
      end subroutine Get_flux_x
!-------------------------------------------------------------------------------
      subroutine Set_flux_x(self, flux_x)
      implicit none
      real*8, intent(in) :: flux_x(:,:,:,:)
      type(t_Advection), intent(inOut) :: self
      self%flux_x(:,:,:,:) = flux_x(:,:,:,:)
      return
      end subroutine Set_flux_x
!-------------------------------------------------------------------------------
      subroutine Get_flux_y(self, flux_y)
      implicit none
      real*8, intent(out) :: flux_y(:,:,:,:)
      type(t_Advection), intent(in) :: self
      flux_y(:,:,:,:) = self%flux_y(:,:,:,:)
      return
      end subroutine Get_flux_y
!-------------------------------------------------------------------------------
      subroutine Set_flux_y(self, flux_y)
      implicit none
      real*8, intent(in) :: flux_y(:,:,:,:)
      type(t_Advection), intent(inOut) :: self
      self%flux_y(:,:,:,:) = flux_y(:,:,:,:)
      return
      end subroutine Set_flux_y
!-------------------------------------------------------------------------------
      subroutine Get_flux_z(self, flux_z)
      implicit none
      real*8, intent(out) :: flux_z(:,:,:,:)
      type(t_Advection), intent(in) :: self
      flux_z(:,:,:,:) = self%flux_z(:,:,:,:)
      return
      end subroutine Get_flux_z
!-------------------------------------------------------------------------------
      subroutine Set_flux_z(self, flux_z)
      implicit none
      real*8, intent(in) :: flux_z(:,:,:,:)
      type(t_Advection), intent(inOut) :: self
      self%flux_z(:,:,:,:) = flux_z(:,:,:,:)
      return
      end subroutine Set_flux_z
!-------------------------------------------------------------------------------
      subroutine Get_air_mass_flux(self, air_mass_flux)
      implicit none
      real*8, intent(out) :: air_mass_flux(:,:,:,:)
      type(t_Advection), intent(in) :: self
      air_mass_flux(:,:,:,:) = self%air_mass_flux(:,:,:,:)
      return
      end subroutine Get_air_mass_flux
!-------------------------------------------------------------------------------
      subroutine Set_air_mass_flux(self, air_mass_flux)
      implicit none
      real*8, intent(in) :: air_mass_flux(:,:,:,:)
      type(t_Advection), intent(inOut) :: self
      self%air_mass_flux(:,:,:,:) = air_mass_flux(:,:,:,:)
      return
      end subroutine Set_air_mass_flux
!-------------------------------------------------------------------------------
      subroutine Get_psf_flux(self, psf_flux)
      implicit none
      real*8, intent(out) :: psf_flux(:,:)
      type(t_Advection), intent(in) :: self
      psf_flux(:,:) = self%psf_flux(:,:)
      return
      end subroutine Get_psf_flux
!-------------------------------------------------------------------------------
      subroutine Set_psf_flux(self, psf_flux)
      implicit none
      real*8, intent(in) :: psf_flux(:,:)
      type(t_Advection), intent(inOut) :: self
      self%psf_flux(:,:) = psf_flux(:,:)
      return
      end subroutine Set_psf_flux
!-------------------------------------------------------------------------------
      end module GmiAdvectionMethod_mod
!
