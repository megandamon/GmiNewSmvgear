!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlFreq2_mod
!
! !INTERFACE:
!
      module GmiControlFreq2_mod
!
! !USES:
      use  GmiSub2Glob_mod      , only : subDomain2Global
      use  GmiMessagePassing_mod, only : synchronizeGroup
      use GmiFileOperations_mod , only : makeOutfileName
      use m_netcdf_io_close     , only : Nccl, Nccl_Noerr
      use m_netcdf_io_write     , only : Ncwr_2d_Int, Ncwr_5d, Ncwr_3d, Ncwr_2d
      use m_netcdf_io_write     , only : Ncwr_1d, Ncwr_Scal, Ncwr_2d_Char, Ncwr_4d
      use m_netcdf_io_create    , only : Ncdo_Sync, Nccr_Wr
      use m_netcdf_io_define    , only : NcDef_glob_attributes, NcDef_dimension
      use m_netcdf_io_define    , only : NcDef_var_attributes, NcDef_variable
      use m_netcdf_io_define    , only : NcSetFill, NcEnd_def
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentrationSurf
      use GmiSpcConcentrationMethod_mod, only : Get_concentrationColTrop
      use GmiSpcConcentrationMethod_mod, only : Get_concentrationColCombo
      use GmiChemistryMethod_mod , only : t_Chemistry, Get_overheadO3col, &
     &       Get_const_labels
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      use GmiDomainDecomposition_mod, only : t_gmiDomain
      use GmiDomainDecomposition_mod, only : Get_rootProc, Get_procID, Get_iAmRootProc
      use GmiDomainDecomposition_mod, only : Get_numDomains
      use GmiDomainDecomposition_mod, only : Get_communicatorWorld
      use GmiDomainDecomposition_mod, only : Get_map1_u
!      use GmiDomainDecomposition_mod, only : Get_eastWestFlag, Get_northSouthPole
      use GmiDomainDecomposition_mod, only : Get_mcorGlob, Get_latdeg, Get_londeg

      use GmiGrid_mod, only : t_gmiGrid, Get_numSpecies, &
     &                        Get_i1, Get_i2, Get_ju1, Get_j2, Get_k1, Get_k2, &
     &                        Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl,     &
     &                        Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jhi_gl, &
     &                        Get_ilo, Get_ihi, Get_julo, Get_jhi

      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Get_pr_diag,         &
     &       Get_hdr_var_name, Get_rec_dim_name, Get_prs_dim_name,             &
     &       Get_problem_name, Get_freq2_name, Get_freq2_description,          &
     &       Get_pr_const_column_freq2, Get_pr_psf_freq2, Get_pr_kel_freq2,    &
     &       Get_pr_const_surface_freq2, Get_pr_freq2, Get_pr_const_freq2,     &
     &       Get_pr_tropopausePress_freq2, Get_pr_mass_freq2, Get_k1_freq2,    &
     &       Get_pr_grid_height_freq2, Get_pr_rel_hum_freq2, Get_k2_freq2,     &
     &       Get_pr_metwater_freq2, Get_do_mean_freq2, Get_do_day1_freq2,      &
     &       Get_do_last_tstep_freq2, Get_pr_overheadO3col_freq2,              &
     &       Get_freq2_species, Get_freq2_species_num, Get_pr_freq2_period,    &
     &       Get_pr_at_time_freq2, Get_iRange_freq2, Get_jRange_freq2,         &
     &       Get_pr_potentialVorticity_freq2, Get_hdf_dim_name,                &
     &       Get_lon_dim_name, Get_lat_dim_name

      use GmiTimeControl_mod, only : t_GmiClock, GmiSplitDateTime,            &
     &       Get_curGmiDate, Get_curGmiTime, Get_gmiSeconds, Get_gmiTimeStep, &
     &       Get_numTimeSteps, isOutFreqTime

      use GmiMetFieldsControl_mod, only : t_metFields, Get_pt, Get_mdt, Get_ai,&
     &       Get_bi, Get_am, Get_bm, Get_kel, Get_pctm2, Get_tropopausePress,  &
     &       Get_mass, Get_humidity, Get_potentialVorticity,            &
     &       Get_gridBoxHeight, Get_relativeHumidity

!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: controlOutputFreq2, initializeOutputFreq2
      public  :: finalizeOutputFreq2

#     include "gmi_phys_constants.h"
#     include "gmi_diag_constants_llnl.h"
#     include "GmiParameters.h"
!
! !DESCRIPTION:
! Contains the necessary routines for the freq2 outputs. Three routines
! are visible from the outside:
! %
! \begin{description}
! \item[initializeOutputFreq2:] called once in the initialization steps to
!      open the netCDF file, define variables in the file, write out header
!      information and allocate variables.
! \item[controlOutputFreq2:] called at the end of each model time step to
!      update variables and if needed for communications and writing out 
!      data in the file.
! \item[finalizeOutputFreq2:] called at the end of the integration to
!      deallocate all the variables and close the netCDF file.
! \end{description}
! %
! The above routines have as arguments derived types.
! %
! This module is self-contained. All the main routines are called by all
! the processors (master included) but each processor does not execute 
! all portions of the code. 
!
! There are two categories of arrays: (1) the ones manipulate by all the
! processors, and (2) those handle by worker processors only.
!
      ! Variables manipulate by all processors

      real*8 , pointer, save :: const_freq2_nc               (:,:,:) => null()
      real*8 , pointer, save :: constColTrop_freq2_nc        (:,:,:) => null()
      real*8 , pointer, save :: constColCombo_freq2_nc       (:,:,:) => null()
      real*8 , pointer, save :: constSurface_freq2_nc        (:,:,:) => null()
      real*8 , pointer, save :: psf_freq2_nc                   (:,:) => null()
      real*8 , pointer, save :: tropopausePress_freq2_nc       (:,:) => null()
      real*8 , pointer, save :: kel_freq2_nc                 (:,:,:) => null()
      real*8 , pointer, save :: mass_freq2_nc                (:,:,:) => null()
      real*8 , pointer, save :: grid_height_freq2_nc         (:,:,:) => null()
      real*8 , pointer, save :: rel_hum_freq2_nc             (:,:,:) => null()
      real*8 , pointer, save :: metwater_freq2_nc            (:,:,:) => null()
      real*8 , pointer, save :: overheadO3col_freq2_nc       (:,:,:) => null()
      real*8 , pointer, save :: potentialVorticity_freq2_nc  (:,:,:) => null()
!
      ! Variables manipulate by worker processors only
      real*8 , pointer, save :: psf_freq2_mean                 (:,:) => null()
      real*8 , pointer, save :: kel_freq2_mean               (:,:,:) => null()
      real*8 , pointer, save :: mass_freq2_mean              (:,:,:) => null()
      real*8 , pointer, save :: const_freq2_mean           (:,:,:,:) => null()
      real*8 , pointer, save :: rel_hum_freq2_mean           (:,:,:) => null()
      real*8 , pointer, save :: metwater_freq2_mean          (:,:,:) => null()
      real*8 , pointer, save :: grid_height_freq2_mean       (:,:,:) => null()
      real*8 , pointer, save :: constColTrop_freq2_mean      (:,:,:) => null()
      real*8 , pointer, save :: constColCombo_freq2_mean     (:,:,:) => null()
      real*8 , pointer, save :: constSurface_freq2_mean     (:,:,:) => null()
      real*8 , pointer, save :: overheadO3col_freq2_mean     (:,:,:) => null()
      real*8 , pointer, save :: tropopausePress_freq2_mean     (:,:) => null()
      real*8 , pointer, save :: potentialVorticity_freq2_mean(:,:,:) => null()

      ! netCDF file identifier
      integer,          save :: ncid_freq2      

      ! Counter for the number of records
      integer,          save :: rnum_out_freq2      

      ! Grid information of the variables to be written out
      integer,          save :: ix1, ix2     ! longitude
      integer,          save :: i1, i2
      integer,          save :: i1_gl, i2_gl
      integer,          save :: ilo_gl, ihi_gl, ilo, ihi
      integer,          save :: jx1, jx2   ! latitude
      integer,          save :: ju1, j2 
      integer,          save :: ju1_gl, j2_gl, julo, jhi
      integer,          save :: julo_gl, jhi_gl
      integer,          save :: kx1, kx2   ! vertical
      integer,          save :: k1, k2
      integer,          save :: numLon     ! number of longitudes
      integer,          save :: numLat     ! number of latitudes
      integer,          save :: numVert    ! number of vertical levels
      integer,          save :: numSpecies ! number of species

      logical,          save :: iAmRootProc
      integer,          save :: procID, rootProc
      integer,          save :: numDomains
      integer,          save :: commuWorld
      integer, pointer, save :: map1_u       (:,:,:) => null()
!      integer, pointer, save :: ewflag           (:) => null()
!      integer, pointer, save :: nspole           (:) => null()

      integer,            save :: nhdf

      logical,            save :: pr_diag
      character (len=GMI_MAXSTR), save :: chem_mecha
      character (len=MAX_LENGTH_VAR_NAME), save :: rec_dim_name, prs_dim_name, lon_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: hdr_var_name, hdf_dim_name, lat_dim_name
      character (len=80), save :: freq2_description
      logical           , save :: pr_const_column_freq2
      logical           , save :: pr_const_surface_freq2
      logical           , save :: pr_freq2
      logical           , save :: pr_const_freq2
      logical           , save :: pr_kel_freq2
      logical           , save :: pr_psf_freq2
      logical           , save :: pr_tropopausePress_freq2
      logical           , save :: pr_potentialVorticity_freq2
      logical           , save :: pr_mass_freq2
      logical           , save :: pr_grid_height_freq2
      logical           , save :: pr_rel_hum_freq2
      logical           , save :: pr_metwater_freq2
      logical           , save :: do_mean_freq2
      logical           , save :: do_day1_freq2
      logical           , save :: do_last_tstep_freq2
      logical           , save :: pr_overheadO3col_freq2
      integer           , save :: freq2_species    (MAX_NUM_CONST_GIO)
      integer           , save :: freq2_species_num
      real*8            , save :: pr_freq2_period
      real*8            , save :: pr_at_time_freq2
      real*8            , save :: mdt

!
! !AUTHOR:
! Jules Kouatchou, Jules.Kouatchou-1@nasa.gov
!
! !REVISION HISTORY:
! 12Nov2007: Original code.
!EOP
!-----------------------------------------------------------------------------
      contains
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initializeOutputFreq2
!
! !INTERFACE:
!
      subroutine initializeOutputFreq2(gmiGrid, gmiDomain, Diagnostics, &
     &                Chemistry, metFields, chemMecha)
!
! !USES:
!
      implicit none
!
! !INPUT PARAMETERS:
      character(len=*) , intent(in) :: chemMecha
      type(t_gmiGrid  ), intent(in) :: gmiGrid  
      type(t_gmiDomain), intent(in) :: gmiDomain
      type(t_Chemistry), intent(in) :: Chemistry
      type(t_metFields), intent(in) :: metFields
      type(t_Diagnostics), intent(in) :: Diagnostics
!
! !DESCRIPTION:
! The routines (1) opens a netCDF file, (2) defines variables in the file,
! (3) writes out header data in the file, and (4) allocates varaibles.
!
! !LOCAL VARIABLES:      
      character (len=75)  :: err_msg
      integer :: in1, k1_freq2, k2_freq2 
      integer :: iRange_freq2(2), jRange_freq2(2)
!      integer, allocatable :: ewflag_org(:)
!      integer, allocatable :: nspole_org(:)
      character (len=80                  ) :: freq2_name
      character (len=MAX_LENGTH_FILE_NAME) :: fname
      character (len=MAX_LENGTH_FILE_NAME) :: problem_name
      character (len=MAX_LENGTH_SPECIES_NAME), pointer :: const_labels(:)
!EOP
!-----------------------------------------------------------------------------
!BOC
      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write (6,*) 'initializeOutputFreq2 called by ', procID

      chem_mecha = trim(chemMecha)

      call Get_ilo (gmiGrid, ilo )
      call Get_ihi (gmiGrid, ihi )
      call Get_julo(gmiGrid, julo )
      call Get_jhi (gmiGrid, jhi )
      call Get_i1    (gmiGrid, i1   )
      call Get_i2    (gmiGrid, i2   )
      call Get_ju1   (gmiGrid, ju1  )
      call Get_j2    (gmiGrid, j2   )
      call Get_k1    (gmiGrid, k1   )
      call Get_k2    (gmiGrid, k2   )
      call Get_i1_gl (gmiGrid, i1_gl )
      call Get_i2_gl (gmiGrid, i2_gl )
      call Get_ju1_gl(gmiGrid, ju1_gl)
      call Get_j2_gl (gmiGrid, j2_gl )
      call Get_ilo_gl (gmiGrid, ilo_gl )
      call Get_ihi_gl (gmiGrid, ihi_gl )
      call Get_julo_gl(gmiGrid, julo_gl )
      call Get_jhi_gl (gmiGrid, jhi_gl )
      call Get_numSpecies (gmiGrid, numSpecies )

      call Get_rootProc         (gmiDomain, rootProc  )
      call Get_numDomains       (gmiDomain, numDomains )
      call Get_iAmRootProc      (gmiDomain, iAmRootProc )
      call Get_communicatorWorld(gmiDomain, commuWorld)

!      allocate(ewflag (i1:i2))
!      allocate(ewflag_org (ilo_gl:ihi_gl))
!      call Get_eastWestFlag(gmiDomain, ewflag_org)
!      ewflag(i1:i2) = ewflag_org(i1:i2)
!      deallocate(ewflag_org)

!      allocate(nspole(ju1:j2))
!      allocate(nspole_org(julo_gl:jhi_gl))
!      call Get_northSouthPole(gmiDomain, nspole_org)
!      nspole(ju1:j2) = nspole_org(ju1:j2)
!      deallocate(nspole_org)

      allocate(map1_u(2, 2, numDomains))
      call Get_map1_u (gmiDomain, map1_u)

      call Get_iRange_freq2(Diagnostics, iRange_freq2)
      call Get_jRange_freq2(Diagnostics, jRange_freq2)

      call Get_k1_freq2(Diagnostics, k1_freq2)
      call Get_k2_freq2(Diagnostics, k2_freq2)

      ix1     = iRange_freq2(1)
      ix2     = iRange_freq2(2)
      jx1     = jRange_freq2(1)
      jx2     = jRange_freq2(2)
      kx1     = k1_freq2
      kx2     = k2_freq2

      numLon  = ix2 - ix1 + 1
      numLat  = jx2 - jx1 + 1
      numVert = kx2 - kx1 + 1

      call Get_mdt(metFields, mdt)

      nhdf = NETCDF_HDF

      call Get_hdr_var_name(Diagnostics, hdr_var_name)
      call Get_rec_dim_name(Diagnostics, rec_dim_name)
      call Get_prs_dim_name(Diagnostics, prs_dim_name)
      call Get_hdf_dim_name (Diagnostics, hdf_dim_name)
      call Get_lat_dim_name (Diagnostics, lat_dim_name)
      call Get_lon_dim_name (Diagnostics, lon_dim_name)

      call Get_pr_freq2(Diagnostics, pr_freq2)
      call Get_pr_kel_freq2(Diagnostics, pr_kel_freq2)
      call Get_pr_psf_freq2(Diagnostics, pr_psf_freq2)
      call Get_pr_mass_freq2(Diagnostics, pr_mass_freq2)
      call Get_do_mean_freq2(Diagnostics, do_mean_freq2)
      call Get_do_day1_freq2(Diagnostics, do_day1_freq2)
      call Get_freq2_species(Diagnostics, freq2_species)
      call Get_pr_const_freq2(Diagnostics, pr_const_freq2)
      call Get_pr_freq2_period(Diagnostics, pr_freq2_period)
      call Get_pr_rel_hum_freq2(Diagnostics, pr_rel_hum_freq2)
      call Get_pr_at_time_freq2(Diagnostics, pr_at_time_freq2)
      call Get_freq2_description(Diagnostics, freq2_description)
      call Get_freq2_species_num(Diagnostics, freq2_species_num)
      call Get_pr_metwater_freq2(Diagnostics, pr_metwater_freq2)
      call Get_do_last_tstep_freq2(Diagnostics, do_last_tstep_freq2)
      call Get_pr_grid_height_freq2(Diagnostics, pr_grid_height_freq2)
      call Get_pr_const_column_freq2(Diagnostics, pr_const_column_freq2)
      call Get_pr_const_surface_freq2(Diagnostics, pr_const_surface_freq2)
      call Get_pr_overheadO3col_freq2(Diagnostics, pr_overheadO3col_freq2)
      call Get_pr_tropopausePress_freq2(Diagnostics, pr_tropopausePress_freq2)
      call Get_pr_potentialVorticity_freq2(Diagnostics, pr_potentialVorticity_freq2)

      allocate(const_labels(numSpecies))

      if (iAmRootProc) then

         !###########################
         ! Initialize the output file
         !###########################

         ! Determine the file name
         call Get_freq2_name  (Diagnostics, freq2_name  )
         call Get_problem_name(Diagnostics, problem_name)

         in1 = Len_Trim (freq2_name)
         call makeOutfileName (fname, '.' // freq2_name(1:in1) // '.nc', &
     &                         problem_name)

         ! Create the file and assign a file identifier
         call Nccr_Wr (ncid_freq2, fname)

         ! Define the variables in the file
         call Define_Netcdf_Out_Freq2()

         ! Write header data

         call Get_const_labels(Chemistry, const_labels)

         call Write_Netcdf_Hdr_Freq2 (gmiDomain, metFields, const_labels)

         call Ncdo_Sync (ncid_freq2)

         ! Initialize the counter for the number of records
         rnum_out_freq2 = 1
      end if

      !########################
      ! Allocation of variables
      !########################

      call allocateVariables_freq2()

      return

      end subroutine initializeOutputFreq2
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINES: allocateVariables_freq2
!
! !INTERFACE:
!
      subroutine allocateVariables_freq2()
!
      implicit none
!
! !DESCRIPTION:
! Allocates variables necessary to produce freq2 outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'allocateVariables_freq2 called by ', procID

      Allocate (const_freq2_nc(i1:i2, ju1:j2, kx1:kx2))
      const_freq2_nc = 0.0d0

      if (pr_const_column_freq2) then
         if ((chem_mecha == 'troposphere') .or. (chem_mecha == 'strat_trop')  &
     &         .or. chem_mecha == 'strat_trop_aerosol' ) then
            Allocate (constColTrop_freq2_nc(i1:i2, ju1:j2, 1:freq2_species_num))
            constColTrop_freq2_nc = 0.0d0

            if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol') then
               Allocate (constColCombo_freq2_nc(i1:i2, ju1:j2, freq2_species_num))
               constColCombo_freq2_nc = 0.0d0
            end if
         end if
      end if

      if (pr_const_surface_freq2) then
         Allocate (constSurface_freq2_nc(i1:i2, ju1:j2, 1:freq2_species_num))
         constSurface_freq2_nc = 0.0d0
      end if

      if (pr_psf_freq2) then
         Allocate (psf_freq2_nc (i1:i2, ju1:j2))
         psf_freq2_nc = 0.0d0
      end if

      if (pr_tropopausePress_freq2) then
         Allocate (tropopausePress_freq2_nc (i1:i2, ju1:j2))
         tropopausePress_freq2_nc = 0.0d0
      end if

      if (pr_potentialVorticity_freq2) then
         Allocate (potentialVorticity_freq2_nc (i1:i2, ju1:j2, kx1:kx2))
         potentialVorticity_freq2_nc = 0.0d0
      end if

      if (pr_kel_freq2) then
         Allocate (kel_freq2_nc (i1:i2, ju1:j2, kx1:kx2))
         kel_freq2_nc = 0.0d0
      end if

      if (pr_overheadO3col_freq2) then
         Allocate (overheadO3col_freq2_nc (i1:i2, ju1:j2, kx1:kx2))
         overheadO3col_freq2_nc = 0.0d0
      end if

      if (pr_mass_freq2) then
         Allocate (mass_freq2_nc(i1:i2, ju1:j2, kx1:kx2))
         mass_freq2_nc = 0.0d0
      end if

      if (pr_grid_height_freq2) then
         Allocate (grid_height_freq2_nc(i1:i2, ju1:j2, kx1:kx2))
         grid_height_freq2_nc = 0.0d0
      end if

      if (pr_rel_hum_freq2) then
         Allocate (rel_hum_freq2_nc(i1:i2, ju1:j2, kx1:kx2))
         rel_hum_freq2_nc = 0.0d0
      end if

      if (pr_metwater_freq2) then
         Allocate (metwater_freq2_nc (i1:i2, ju1:j2, kx1:kx2))
         metwater_freq2_nc = 0.0d0
      end if

      if (do_mean_freq2) then

         if (freq2_species_num > 0) then
            Allocate(const_freq2_mean(i1:i2,ju1:j2,kx1:kx2,freq2_species_num))
            const_freq2_mean = 0.0d0

            if (pr_const_column_freq2) then
               if ((chem_mecha == 'troposphere') .or. (chem_mecha == 'strat_trop') &
     &               .or. chem_mecha == 'strat_trop_aerosol' ) then
                   Allocate(constColTrop_freq2_mean(i1:i2, ju1:j2, freq2_species_num))
                   constColTrop_freq2_mean = 0.0d0
                   if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then
                      Allocate(constColCombo_freq2_mean(i1:i2, ju1:j2, freq2_species_num))
                      constColCombo_freq2_mean = 0.0d0
                   end if
                end if
            end if

            if (pr_const_surface_freq2) then
               Allocate(constSurface_freq2_mean(i1:i2, ju1:j2, freq2_species_num))
               constSurface_freq2_mean = 0.0d0
            end if
         endif

         if (pr_potentialVorticity_freq2) then
            Allocate (potentialVorticity_freq2_mean (i1:i2, ju1:j2, kx1:kx2))
            potentialVorticity_freq2_mean = 0.0d0
         end if

         if (pr_kel_freq2) then
            Allocate (kel_freq2_mean (i1:i2, ju1:j2, kx1:kx2))
            kel_freq2_mean = 0.0d0
         end if

         if (pr_overheadO3col_freq2) then
            Allocate (overheadO3col_freq2_mean (i1:i2, ju1:j2, kx1:kx2))
            overheadO3col_freq2_mean = 0.0d0
         end if

         if (pr_psf_freq2) then
            Allocate (psf_freq2_mean (i1:i2, ju1:j2))
            psf_freq2_mean = 0.0d0
         end if

         if (pr_tropopausePress_freq2) then
            Allocate (tropopausePress_freq2_mean (i1:i2, ju1:j2))
            tropopausePress_freq2_mean = 0.0d0
         end if

         if (pr_mass_freq2) then
            Allocate (mass_freq2_mean(i1:i2, ju1:j2, kx1:kx2))
            mass_freq2_mean = 0.0d0
         end if

         if (pr_grid_height_freq2) then
            Allocate (grid_height_freq2_mean(i1:i2, ju1:j2, kx1:kx2))
            grid_height_freq2_mean = 0.0d0
         end if

         if (pr_rel_hum_freq2) then
            Allocate (rel_hum_freq2_mean(i1:i2, ju1:j2, kx1:kx2))
            rel_hum_freq2_mean = 0.0d0
         end if

         if (pr_metwater_freq2) then
            Allocate (metwater_freq2_mean(i1:i2, ju1:j2, kx1:kx2))
            metwater_freq2_mean = 0.0d0
         end if

      end if

      return

      end subroutine allocateVariables_freq2
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINES: finalizeOutputFreq2
!
! !INTERFACE:
!
      subroutine finalizeOutputFreq2()
!
      implicit none
!
! !DESCRIPTION:
! Deallocates variables necessary to produce freq2 outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'finalizeOutputFreq2 called by ', procID

      deallocate (const_freq2_nc)

      if (pr_const_column_freq2) then
         if ((chem_mecha == 'troposphere') .or. (chem_mecha == 'strat_trop')  &
     &        .or. chem_mecha == 'strat_trop_aerosol' ) then
            deallocate (constColTrop_freq2_nc)

            if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then
               deallocate (constColCombo_freq2_nc)
            end if
         end if
      end if

      if (pr_const_surface_freq2) then
         deallocate (constSurface_freq2_nc)
      end if

      if (pr_psf_freq2) then
         deallocate (psf_freq2_nc )
      end if

      if (pr_tropopausePress_freq2) then
         deallocate (tropopausePress_freq2_nc )
      end if

      if (pr_potentialVorticity_freq2) then
         deallocate (potentialVorticity_freq2_nc )
      end if

      if (pr_kel_freq2) then
         deallocate (kel_freq2_nc )
      end if

      if (pr_overheadO3col_freq2) then
         deallocate (overheadO3col_freq2_nc )
      end if

      if (pr_mass_freq2) then
         deallocate (mass_freq2_nc)
      end if

      if (pr_grid_height_freq2) then
         deallocate (grid_height_freq2_nc)
      end if

      if (pr_rel_hum_freq2) then
         deallocate (rel_hum_freq2_nc)
      end if

      if (pr_metwater_freq2) then
         deallocate (metwater_freq2_nc )
      end if

      if (do_mean_freq2) then

         if (freq2_species_num > 0) then
            deallocate(const_freq2_mean)

            if (pr_const_column_freq2) then
               if ((chem_mecha == 'troposphere') .or. (chem_mecha == 'strat_trop')  &
     &               .or. chem_mecha == 'strat_trop_aerosol' ) then
                   deallocate(constColTrop_freq2_mean)
                   if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then
                      deallocate(constColCombo_freq2_mean)
                   end if
                end if
            end if

            if (pr_const_surface_freq2) then
               deallocate(constSurface_freq2_mean)
            end if
         endif

         if (pr_potentialVorticity_freq2) then
            deallocate (potentialVorticity_freq2_mean )
         end if

         if (pr_kel_freq2) then
            deallocate (kel_freq2_mean )
         end if

         if (pr_overheadO3col_freq2) then
            deallocate (overheadO3col_freq2_mean )
         end if

         if (pr_psf_freq2) then
            deallocate (psf_freq2_mean )
         end if

         if (pr_tropopausePress_freq2) then
            deallocate (tropopausePress_freq2_mean )
         end if

         if (pr_mass_freq2) then
            deallocate (mass_freq2_mean)
         end if

         if (pr_grid_height_freq2) then
            deallocate (grid_height_freq2_mean)
         end if

         if (pr_rel_hum_freq2) then
            deallocate (rel_hum_freq2_mean)
         end if

         if (pr_metwater_freq2) then
            deallocate (metwater_freq2_mean)
         end if

      end if

      if (iAmRootProc) call Nccl_Noerr (ncid_freq2)

      return

      end subroutine finalizeOutputFreq2
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: controlOutFreq2
!
! !INTERFACE:
!
      subroutine controlOutputFreq2 (last_tstp, Chemistry, SpeciesConcentration, &
     &                  gmiClock, metFields)

      use GmiTimeControl_mod, only : GmiSplitDateTime
      use GmiChemistryMethod_mod , only : t_Chemistry

      implicit none
!
! !INPUT PARAMETERS:
      logical                     , intent(in) :: last_tstp
      type(t_gmiClock )           , intent(in) :: gmiClock 
      type(t_Chemistry)           , intent(in) :: Chemistry
      type(t_metFields)           , intent(in) :: metFields
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !DESCRIPTION:
! This routine controls the freq2 file output. It is called at each time step.
!
! !LOCAL VARIABLES:
      logical :: time_for_freq2
      integer :: day, month, idumyear, ndt, num_time_steps, nhms, nymd
      real*8  :: gmi_sec, tstep
      integer, save :: month_save = -999
      logical, save :: printed_on_this_day = .false.
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'controlOutputFreq2 called by ', procID

      call Get_curGmiDate(gmiClock, nymd)
      call Get_curGmitime(gmiClock, nhms)
      call Get_numTimeSteps(gmiClock, num_time_steps)
      call Get_gmiSeconds(gmiClock, gmi_sec)
      call Get_gmiTimeStep(gmiClock, tstep)

      ndt = Nint(tstep)

      call GmiSplitDateTime (nymd, idumyear, month, day)

      if (month_save == -999) month_save = month

      time_for_freq2 = .false.

      call isOutFreqTime(time_for_freq2, printed_on_this_day, &
     &           month_save, month, day, nhms, ndt, &
     &           gmi_sec, pr_freq2_period, pr_at_time_freq2)

      month_save = month

      if (do_day1_freq2)then
         if ((day == 2) .and. (nhms == 0) ) then
            time_for_freq2 = .true.
         end if
      end if

      if (do_last_tstep_freq2) then
        if (last_tstp .and. (.not. time_for_freq2)) then
          if (Mod (num_time_steps, Nint (mdt) / ndt) == 0) then

!         ------------------------------------------------------------
!         Always update restart file after the last step if you are at
!         the end of a met record.
!         ------------------------------------------------------------

           time_for_freq2 = .true.
          endif
        end if
      end if


      if (time_for_freq2 .or. do_mean_freq2) then
         !====================
         call Prep_Netcdf_Freq2 (time_for_freq2, Chemistry, &
     &             SpeciesConcentration, metFields)
         !====================
      endif

!     ================
      if (time_for_freq2) then
!     ================

!        =====================
         call Write_Netcdf_Freq2 (nymd, nhms, gmi_sec)
!        =====================

!        ======================
         call bufferOutput_Freq2 (SpeciesConcentration)
!        ======================

        if (iAmRootProc) rnum_out_freq2 = rnum_out_freq2 + 1

!     ======
      end if
!     ======

      return

      end subroutine controlOutputFreq2
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutput_Freq2
!
! !INTERFACE:
!
      subroutine bufferOutput_Freq2 (SpeciesConcentration)
!
! !USES:
      use GmiNcOutputSlice_mod, only : prepSlice4d, prepSlice4d_Nomean
      use GmiNcOutputSlice_mod, only : writeSlice4d
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !DESCRIPTION:
! This routine buffers the freq2 data to the output file.  The Slaves send a
! slice of the data to the Master, the Master writes it out, they then send
! another slice and it is written out, and so on.
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_CONST = 3000
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_VAR_NAME), parameter :: CFREQ1_VNAM = 'const_freq2'
      integer :: ic, icx
      real*8, allocatable :: const_freq2 (:,:,:,:)
      real*8, allocatable :: arrayGlob3D (:,:,:)
      type (t_GmiArrayBundle), pointer :: concentration(:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'bufferOutput_Freq2 called by ', procID

      if ((freq2_species_num > 0) .and. (pr_const_freq2)) then
         call Get_concentration(SpeciesConcentration, concentration)
         Allocate (const_freq2(i1:i2, ju1:j2, kx1:kx2, freq2_species_num))

         if (iAmRootProc) allocate(arrayGlob3D(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))

         ic = 0
         do icx = 1, numSpecies
            if (freq2_species(icx) /= 0) then
               ic = ic + 1
               const_freq2(:,:,:,ic) = concentration(icx)%pArray3D(:,:,kx1:kx2)
            end if
         end do

         do ic = 1, freq2_species_num
            if (do_mean_freq2) then
               call prepSlice4d (do_mean_freq2, i1, i2, ju1, j2, &
     &                  kx1, kx2, ic, freq2_species_num, &
     &                  const_freq2, const_freq2_mean, const_freq2_nc)
            else
               call prepSlice4d_Nomean (i1, i2, ju1, j2, kx1, kx2, &
     &                  ic, freq2_species_num, const_freq2, const_freq2_nc)
            end if

            call subDomain2Global(arrayGlob3D, const_freq2_nc, &
     &              i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &              rootProc, procID, map1_u, numDomains,  &
     &              SG_CONST, commuWorld)

            if (iAmRootProc) then
               call writeSlice4d (ix1, ix2, jx1, jx2, kx1, kx2, &
     &                   arrayGlob3D, ncid_freq2, CFREQ1_VNAM, ic, rnum_out_freq2)
            end if

            call synchronizeGroup (commuWorld)
         enddo

         if (iAmRootProc) call Ncdo_Sync (ncid_freq2)
      endif

      return

      end subroutine bufferOutput_Freq2
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Define_Netcdf_Out_Freq2
!
! !INTERFACE:
!
      subroutine Define_Netcdf_Out_Freq2()
!
! !USES:
      use m_ncGeneralOpsOutput, only: Define_Netcdf_Out_Gen
!
      implicit none
!
#     include "netcdf.inc"
!
! !DESCRIPTION:
! This routine makes the necessary definitions for the freq2 species
! concentration NetCDF output file.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_VAR_NAME) :: prsp1_dim_name
      integer :: ierr, nchr1, nchr2, nstr, omode, pos1, varid
      integer :: chrd1(1), chrd2 (1), chrd3(1)
      integer :: hdfd (1)
      integer :: lond (1), latd  (1)
      integer :: freq2d(1)
      integer :: prsd (1), prsp1d(1)
      integer :: recd (1), spcd  (1)
      integer :: spcd_2(1), spcd_3(1), spcd_4(1)
      integer :: strd (1)
      integer :: var2(2), var3(3), var4(4), var5(5)
      character (len=200) :: AttrValue
!!EOP
!-----------------------------------------------------------------------------
!BOC

      if (pr_diag) Write (6,*) 'Define_Netcdf_Out_Freq2 called by ', procID

      nchr1 = MAX_LENGTH_SPECIES_NAME
      nchr2 = 50

      nstr  =  1

!     ==========================
      call Define_Netcdf_Out_Gen  &
!     ==========================
     &  (ncid_freq2, nhdf, numLon, numLat, numVert, hdfd, lond, latd, prsd,  &
     &   hdf_dim_name, lon_dim_name, lat_dim_name, prs_dim_name)

     ! ------------------
     ! Define dimensions.
     ! ------------------

      prsp1_dim_name = prs_dim_name
      pos1 = Len_Trim (prs_dim_name) + 1
      prsp1_dim_name(pos1:pos1+2) = 'p1'

      call NcDef_dimension(ncid_freq2, prsp1_dim_name, numVert+1, prsp1d(1))
                                      !--------------
      call NcDef_dimension(ncid_freq2, 'chr_dim1', nchr1, chrd1(1))
                                      !----------
      call NcDef_dimension(ncid_freq2, 'chr_dim2', nchr2, chrd2(1))
                                      !----------

      if (freq2_species_num > 0) then
         call NcDef_dimension  &
     &         (ncid_freq2, 'freq2_species', freq2_species_num, freq2d(1))
                           !---------------
      end if

      call NcDef_dimension(ncid_freq2, rec_dim_name, NF_UNLIMITED, recd(1))
                                      !------------

      ! -----------------------------------------
      ! Define variables and variable attributes.
      ! -----------------------------------------

      call NcDef_variable (ncid_freq2, 'am', NF_FLOAT, 1, prsd, varid)
                                      !----
      call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Hybrid pressure term')
      call NcDef_var_attributes (ncid_freq2, varid, 'units', 'unitless')

      call NcDef_variable (ncid_freq2, 'bm', NF_FLOAT, 1, prsd, varid)
                                      !----
      call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Hybrid sigma term ')
      call NcDef_var_attributes (ncid_freq2, varid, 'units', 'unitless')

      call NcDef_variable (ncid_freq2, 'ai', NF_FLOAT, 1, prsp1d, varid)
                                      !----
      call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Hybrid pressure edge term')
      call NcDef_var_attributes (ncid_freq2, varid, 'units', 'unitless')

      call NcDef_variable (ncid_freq2, 'bi', NF_FLOAT, 1, prsp1d, varid)
                                      !----
      call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Hybrid sigma edge term')
      call NcDef_var_attributes (ncid_freq2, varid, 'units', 'unitless')

      call NcDef_variable (ncid_freq2, 'pt', NF_FLOAT, 0, prsp1d, varid)
                                      !----
      var2(:) = (/ lond(1), latd(1) /)
      call NcDef_variable (ncid_freq2, 'mcor', NF_FLOAT, 2, var2, varid)
                                      !------
      call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Grid box area')
      call NcDef_var_attributes (ncid_freq2, varid, 'units', 'm^2')

!Freq2 Species
!
      if (freq2_species_num > 0) then
         call NcDef_variable (ncid_freq2, 'freq2_species', NF_FLOAT, 1, freq2d, varid)
                                         !---------------
         call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Freq2 Species index')
         call NcDef_var_attributes (ncid_freq2, varid, 'units', 'unitless')
         call NcDef_var_attributes (ncid_freq2, varid, 'coord_labels', 'freq2_labels')
         call NcDef_var_attributes (ncid_freq2, varid, 'selection_category', 'NULL')

         call NcDef_variable (ncid_freq2, 'freq2_description', NF_CHAR, 0, chrd2, varid)
                                         !-------------------
         call NcDef_var_attributes (ncid_freq2, varid, 'long_name', freq2_description)

         var2(:) = (/ chrd1(1), freq2d(1) /)
         call NcDef_variable (ncid_freq2, 'freq2_labels', NF_CHAR, 2, var2, varid)
                                         !--------------
         call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Freq2 Species name')
         call NcDef_var_attributes (ncid_freq2, varid, 'units', 'unitless')
         call NcDef_var_attributes (ncid_freq2, varid, 'selection_category', 'NULL')

         if (pr_const_freq2) then
            var5(:) = (/ lond(1), latd(1), prsd(1), freq2d(1), recd(1) /)
            call NcDef_variable (ncid_freq2, 'const_freq2', NF_FLOAT, 5, var5, varid)
                                            !-------------
            call NcDef_var_attributes (ncid_freq2, varid, 'long_name','Constituent at freq2')
            call NcDef_var_attributes (ncid_freq2, varid, 'units', 'volume mixing ratio')
         endif

         !  Put Const Column Values
         if (pr_const_column_freq2) then
            if ((chem_mecha == 'troposphere') .or. (chem_mecha == 'strat_trop')  &
     &            .or. chem_mecha == 'strat_trop_aerosol' ) then
               var4(:) = (/ lond(1), latd(1), freq2d(1), recd(1) /)
               call NcDef_variable (ncid_freq2, 'constColTrop_freq2', NF_FLOAT, 4, var4, varid)
                                               !--------------------
               call NcDef_var_attributes (ncid_freq2, varid, 'long_name',  &
     &                                  'Troposphere Constituent Column Freq2')
               call NcDef_var_attributes (ncid_freq2, varid, 'units', 'molec/cm^2')
       
               if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then
                  call NcDef_variable (ncid_freq2, 'constColCombo_freq2', NF_FLOAT, 4, var4, varid)
                                                  !---------------------
                  call NcDef_var_attributes (ncid_freq2, varid, 'long_name',  &
     &                                  'Combo Constituent Column Freq2')
                  call NcDef_var_attributes (ncid_freq2, varid, 'units', 'molec/cm^2')
               end if
            end if
         end if

         !  Put Const Surface Values
         if (pr_const_surface_freq2) then
            var4(:) = (/ lond(1), latd(1), freq2d(1), recd(1) /)
            call NcDef_variable (ncid_freq2, 'constSurface_freq2', NF_FLOAT, 4, var4, varid)
                                            !--------------------
            call NcDef_var_attributes (ncid_freq2, varid, 'long_name',  &
     &                'Constituent Surface Freq2')
            call NcDef_var_attributes (ncid_freq2, varid, 'units', 'ppbv')
         endif

      end if

      ! MetFields related variables
      if (pr_psf_freq2) then
         var3(:) = (/ lond(1), latd(1), recd(1) /)
         call NcDef_variable (ncid_freq2, 'psf_freq2', NF_FLOAT, 3, var3, varid)
                                         !-----------
         call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Surface pressure')
         call NcDef_var_attributes (ncid_freq2, varid, 'units', 'mb')
      end if

      if (pr_tropopausePress_freq2) then
         var3(:) = (/ lond(1), latd(1), recd(1) /)
         call NcDef_variable (ncid_freq2, 'tropopausePress_freq2', NF_FLOAT, 3, var3, varid)
                                         !-----------------------
         call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Tropopause pressure')
         call NcDef_var_attributes (ncid_freq2, varid, 'units', 'mb')
      end if

      if (pr_potentialVorticity_freq2) then
         var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
         call NcDef_variable (ncid_freq2, 'potentialVorticity_freq2', NF_FLOAT, 4, var4, varid)
!                                        !--------------------------
         call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Potential vorticity')
         call NcDef_var_attributes (ncid_freq2, varid, 'units', 'm2 s-1K kg-1')
      end if

      if (pr_kel_freq2) then
         var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
         call NcDef_variable (ncid_freq2, 'kel_freq2', NF_FLOAT, 4, var4, varid)
                                         !-----------
         call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Temperature')
         call NcDef_var_attributes (ncid_freq2, varid, 'units', 'K')
      end if

      if (pr_overheadO3col_freq2) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_freq2, 'overheadO3col_freq2', NF_FLOAT, 4, var4, varid)
!                                       !---------------------
        call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Photolysis overhead ozone column')
        call NcDef_var_attributes (ncid_freq2, varid, 'units', 'DU')
      end if

      if (pr_mass_freq2) then
         var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
         call NcDef_variable (ncid_freq2, 'mass_freq2', NF_FLOAT, 4, var4, varid)
!                                        !------------
         call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Mass')
         call NcDef_var_attributes (ncid_freq2, varid, 'units', 'kg')
      end if

      if (pr_grid_height_freq2) then
         var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
         call NcDef_variable (ncid_freq2, 'grid_height_freq2', NF_FLOAT, 4, var4, varid)
!                                        !-------------------
         call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Grid Box Height')
         call NcDef_var_attributes (ncid_freq2, varid, 'units', 'm')
      end if

      if (pr_rel_hum_freq2) then
         var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
         call NcDef_variable (ncid_freq2, 'rel_hum_freq2', NF_FLOAT, 4, var4, varid)
!                                        !---------------
         call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Relative Humidity')
         call NcDef_var_attributes (ncid_freq2, varid, 'units', '[0-1]')
      end if

      if (pr_metwater_freq2) then
         var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
         call NcDef_variable (ncid_freq2, 'metwater_freq2', NF_FLOAT, 4, var4, varid)
!                                        !----------------
         call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Meteorological Water')
         call NcDef_var_attributes (ncid_freq2, varid, 'units', 'volume mixing ratio')
      end if

! ==================

      var2(:) = (/ hdfd(1), recd(1) /)
      call NcDef_variable (ncid_freq2, hdr_var_name, NF_INT, 2, var2, varid)
!                                     !-------------
      call NcDef_var_attributes (ncid_freq2, varid, 'long_name', 'Header')
      call NcDef_var_attributes (ncid_freq2, varid, 'units', 'gmi_sec, nymd, nhms')
!
      ! -------------------------
      ! Define global attributes.
      ! -------------------------

      call NcDef_glob_attributes (ncid_freq2, 'title',  &
     &     'Gmimod freq2 species concentration file')

      AttrValue = 'DESCRIPTION OF WATER RELATED SPECIES AND VARIABLES'
      call NcDef_glob_attributes(ncid_freq2,'Note', Trim(AttrValue) )

      AttrValue = 'read in from meteorological fields used to drive &
                   transport - only used below tropopause.'

      call NcDef_glob_attributes(ncid_freq2,'metwater', Trim(AttrValue) )

      AttrValue = 'resets metwater to > 3 ppmv when less than 3 ppmv &
                   - not transported'

      call NcDef_glob_attributes(ncid_freq2,'H2O', Trim(AttrValue) )

      AttrValue = 'used in aircraft emissions experiments only &
                   - not transported simulations'

      call NcDef_glob_attributes(ncid_freq2,'H2OAIR', Trim(AttrValue) )

      AttrValue = 'accounts for removal/addition of water by PSC &
                   formation/dehydration and gravitational settling in &
                   stratosphere - modifies h2oclim in stratosphere - &
                   should always be transported'

      call NcDef_glob_attributes(ncid_freq2,'DEHYD', Trim(AttrValue) )

      AttrValue = 'stratospheric water climatology (UARS MLS) read in &
                   from file - not outputted to *const.nc'

      call NcDef_glob_attributes(ncid_freq2,'h2oclim', Trim(AttrValue) )

      call NcSetFill (ncid_freq2, NF_NOFILL, omode)

      call NcEnd_def (ncid_freq2)

      return

      end subroutine Define_Netcdf_Out_Freq2
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Prep_Netcdf_Freq2
!
      subroutine Prep_Netcdf_Freq2 (time_for_freq2, Chemistry, &
     &               SpeciesConcentration, metFields)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical :: time_for_freq2
      type(t_Chemistry)           , intent(in) :: Chemistry
      type(t_metFields)           , intent(in) :: metFields
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !DESCRIPTION:
! This routine prepares the freq2 species concentration NetCDF output.
!
! !LOCAL VARIABLES:
      logical, save :: first = .true.
      integer :: ic, icx, il
      integer, save :: counter = 0
      real*8  :: dcounter
      real*8  :: rhms
      real*8, allocatable :: mass(:,:,:)
      real*8, allocatable :: gridBoxHeight(:,:,:), kel(:,:,:)
      real*8, allocatable :: tropopausePress(:,:), pctm2(:,:)
      real*8, allocatable :: relativeHumidity(:,:,:), humidity(:,:,:)
      real*8, allocatable :: potentialVorticity(:,:,:)
      real*8, allocatable :: constColTrop(:,:,:)
      real*8, allocatable :: constColCombo(:,:,:)
      real*8, allocatable :: constSurface(:,:,:)
      real*8, allocatable :: overheadO3col(:,:,:)
      real*8, allocatable :: concentrationSurf(:,:,:)
      real*8, allocatable :: concentrationColTrop(:,:,:)
      real*8, allocatable :: concentrationColCombo(:,:,:)
      type (t_GmiArrayBundle), pointer :: concentration(:)
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Prep_Netcdf_Freq2 called by ', procID

      !--------------------------------------------
      ! Extract array values from the derived types
      !--------------------------------------------

      call Get_concentration(SpeciesConcentration, concentration)

      if (pr_const_surface_freq2) then
         allocate (constSurface(i1:i2, ju1:j2, freq2_species_num))
         allocate(concentrationSurf(i1:i2, ju1:j2, numSpecies))
         call Get_concentrationSurf(SpeciesConcentration, concentrationSurf)
      end if

      if (pr_const_column_freq2) then
         if ((chem_mecha == 'troposphere') .or. (chem_mecha == 'strat_trop')  &
     &         .or. chem_mecha == 'strat_trop_aerosol' ) then
            allocate (constColTrop(i1:i2, ju1:j2, freq2_species_num))
            allocate(concentrationColTrop(i1:i2, ju1:j2, numSpecies))
            call Get_concentrationColTrop(SpeciesConcentration, concentrationColTrop)

            if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then
               allocate (constColCombo(i1:i2, ju1:j2, freq2_species_num))
               allocate(concentrationColCombo(i1:i2, ju1:j2, numSpecies))
               call Get_concentrationColCombo(SpeciesConcentration, concentrationColCombo)
            end if
         end if
      end if

      if (pr_overheadO3col_freq2)  then
         allocate(overheadO3col(i1:i2, ju1:j2, k1:k2))
         call Get_overheadO3col(Chemistry, overheadO3col)
      end if

      if (pr_tropopausePress_freq2) then
         allocate(tropopausePress(i1:i2, ju1:j2))
         call Get_tropopausePress(metFields, tropopausePress)
      end if

      if (pr_psf_freq2) then
         allocate(pctm2(ilo:ihi, julo:jhi))
         call Get_pctm2(metFields, pctm2)
      end if

      if (pr_kel_freq2) then
         allocate(kel(ilo:ihi, julo:jhi, k1:k2))
         call Get_kel(metFields, kel)
      end if

      if (pr_metwater_freq2) then
         allocate(humidity(i1:i2, ju1:j2, k1:k2))
         call Get_humidity(metFields, humidity)
      end if

      if (pr_mass_freq2) then
         allocate(mass(i1:i2, ju1:j2, k1:k2))
         call Get_mass(metFields, mass)
      end if

      if (pr_potentialVorticity_freq2) then
         allocate(potentialVorticity(i1:i2, ju1:j2, k1:k2))
         call Get_potentialVorticity(metFields, potentialVorticity)
      end if

      if (pr_rel_hum_freq2) then
         allocate(relativeHumidity(i1:i2, ju1:j2, k1:k2))
         call Get_relativeHumidity(metFields, relativeHumidity)
      end if

      if (pr_grid_height_freq2) then
         allocate(gridBoxHeight(i1:i2, ju1:j2, k1:k2))
         call Get_gridBoxHeight(metFields, gridBoxHeight)
      end if

      if (freq2_species_num > 0) then
         ic = 0
         do icx = 1, numSpecies
           if (freq2_species(icx) /= 0) then
              ic = ic + 1
              
              if (pr_const_surface_freq2) then
                 constSurface(:,:,ic) = concentrationSurf(:,:,icx)
              end if

              if (pr_const_column_freq2) then
                 if ((chem_mecha == 'troposphere') .or. (chem_mecha == 'strat_trop')  &
     &                 .or. chem_mecha == 'strat_trop_aerosol' ) then
                    constColTrop(:,:,ic) = concentrationColTrop(:,:,icx)

                    if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then
                       constColCombo(:,:,ic) = concentrationColCombo(:,:,icx)
                    end if
                 end if
              end if

           end if
         end do
      end if


      if (do_mean_freq2) then
         if (counter == 0) then

            ! ----------------
            ! Reset mean sums.
            ! ----------------

            if (freq2_species_num > 0) then
              const_freq2_mean     (:,:,:,:) = 0.0d0

              if (pr_const_surface_freq2) then
                 constSurface_freq2_mean (:,:,:) = 0.0d0
              end if

              if (pr_const_column_freq2) then
                 if ((chem_mecha == 'troposphere') .or. (chem_mecha == 'strat_trop')  &
     &                 .or. chem_mecha == 'strat_trop_aerosol' ) then
                    constColTrop_freq2_mean (:,:,:) = 0.0d0

                    if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then
                       constColCombo_freq2_mean (:,:,:) = 0.0d0
                    end if
                 end if
              end if
            endif

            if (pr_tropopausePress_freq2) tropopausePress_freq2_mean (:,:) = 0.0d0
            if (pr_psf_freq2) psf_freq2_mean (:,:) = 0.0d0
            if (pr_potentialVorticity_freq2) &
     &           potentialVorticity_freq2_mean (:,:,:) = 0.0d0
            if (pr_kel_freq2) kel_freq2_mean (:,:,:) = 0.0d0
            if (pr_overheadO3col_freq2) overheadO3col_freq2_mean (:,:,:) = 0.0d0
            if (pr_mass_freq2) mass_freq2_mean (:,:,:) = 0.0d0
            if (pr_grid_height_freq2) grid_height_freq2_mean (:,:,:) = 0.0d0
            if (pr_rel_hum_freq2) rel_hum_freq2_mean (:,:,:) = 0.0d0
            if (pr_metwater_freq2) metwater_freq2_mean (:,:,:) = 0.0d0
         end if

         ! -----------------
         ! Update mean sums.
         ! -----------------

         if (freq2_species_num > 0) then
            ic = 0
            do icx = 1, numSpecies
               if (freq2_species(icx) /= 0) then
                  ic = ic + 1
                  const_freq2_mean(:,:,:,ic) =  &
     &                      const_freq2_mean(:,:,:,ic) +  &
     &                      concentration(icx)%pArray3D(:,:,kx1:kx2)
                  if (pr_const_surface_freq2) then
                     constSurface_freq2_mean(:,:,ic) =  &
     &                      constSurface_freq2_mean(:,:,ic) + &
     &                      concentrationSurf(:,:,icx)
                  end if

                  if (pr_const_column_freq2) then
                     if ((chem_mecha == 'troposphere') .or. (chem_mecha == 'strat_trop')  &
     &                     .or. chem_mecha == 'strat_trop_aerosol' ) then
                        constColTrop_freq2_mean(:,:,ic) =  &
     &                      constColTrop_freq2_mean(:,:,ic) +  &
     &                      concentrationColTrop(:,:,icx)

                        if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then
                            constColCombo_freq2_mean(:,:,ic) =  &
     &                          constColCombo_freq2_mean(:,:,ic) +  &
     &                          concentrationColCombo(:,:,icx)
                        end if
                     end if
                  end if
               end if
            end do
         end if

         if (pr_tropopausePress_freq2) then
            tropopausePress_freq2_mean(i1:i2,ju1:j2) =  &
     &                tropopausePress_freq2_mean (i1:i2,ju1:j2) + &
     &                tropopausePress(i1:i2,ju1:j2)
         endif

         if (pr_psf_freq2) psf_freq2_mean (i1:i2,ju1:j2)   =  &
     &       psf_freq2_mean (i1:i2,ju1:j2)   + pctm2(i1:i2,ju1:j2)

         if (pr_potentialVorticity_freq2) then
            potentialVorticity_freq2_mean (i1:i2,ju1:j2,:) =  &
     &               potentialVorticity_freq2_mean (i1:i2,ju1:j2,:) + &
     &               potentialVorticity  (i1:i2,ju1:j2,kx1:kx2)
         endif

         if (pr_kel_freq2) kel_freq2_mean (i1:i2,ju1:j2,:) =  &
     &      kel_freq2_mean(i1:i2,ju1:j2,:) + kel(i1:i2,ju1:j2,kx1:kx2)

         if (pr_overheadO3col_freq2) &
     &      overheadO3col_freq2_mean(:,:,:) = overheadO3col_freq2_mean(:,:,:) + &
     &                                    overheadO3col(:,:,kx1:kx2)

         if (pr_mass_freq2) then
            mass_freq2_mean (i1:i2,ju1:j2,:) =  &
     &            mass_freq2_mean(i1:i2,ju1:j2,:) + mass(i1:i2,ju1:j2,kx1:kx2)
         end if

         if (pr_grid_height_freq2) then
            grid_height_freq2_mean (i1:i2,ju1:j2,:) =  &
     &                   grid_height_freq2_mean (i1:i2,ju1:j2,:) + &
     &                   gridBoxHeight (i1:i2,ju1:j2,kx1:kx2)
         end if

         if (pr_rel_hum_freq2) then
            rel_hum_freq2_mean (i1:i2,ju1:j2,:) =  &
     &               rel_hum_freq2_mean (i1:i2,ju1:j2,:) + &
     &               relativeHumidity (i1:i2,ju1:j2,kx1:kx2)
         end if

         if (pr_metwater_freq2) metwater_freq2_mean (i1:i2,ju1:j2,:) =  &
     &      metwater_freq2_mean (i1:i2,ju1:j2,:) + &
     &      humidity (i1:i2,ju1:j2,kx1:kx2) * MWTAIR / (MWTH2O * GPKG)

         counter = counter + 1

         if (time_for_freq2) then

            ! ----------------
            ! Calculate means.
            ! ----------------

            dcounter = counter

            if (freq2_species_num > 0) then
               do ic = 1, freq2_species_num
                  const_freq2_mean(:,:,:,ic) = const_freq2_mean(:,:,:,ic) / &
     &                                         dcounter

                  if (pr_const_surface_freq2) then
                     constSurface_freq2_mean(:,:,ic) = &
     &                    constSurface_freq2_mean(:,:,ic)/dcounter
                  end if
            
                  if (pr_const_column_freq2) then
                     if ((chem_mecha == 'troposphere') .or. (chem_mecha == 'strat_trop')  &
     &                     .or. chem_mecha == 'strat_trop_aerosol' ) then
                        constColTrop_freq2_mean(:,:,ic) =  &
     &                       constColTrop_freq2_mean(:,:,ic)/dcounter

                        if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then
                           constColCombo_freq2_mean(:,:,ic) =  &
     &                          constColCombo_freq2_mean(:,:,ic)/dcounter
                        end if
                     end if
                  end if
               end do
            end if

            if (pr_tropopausePress_freq2) tropopausePress_freq2_mean (:,:) = &
               tropopausePress_freq2_mean(:,:) / dcounter

            if (pr_psf_freq2) psf_freq2_mean(:,:) = psf_freq2_mean(:,:) / dcounter

            if (pr_potentialVorticity_freq2) potentialVorticity_freq2_mean(:,:,:) = &
     &         potentialVorticity_freq2_mean(:,:,:) / dcounter

            if (pr_kel_freq2) kel_freq2_mean(:,:,:) = kel_freq2_mean(:,:,:) / dcounter

            if (pr_overheadO3col_freq2) then
               overheadO3col_freq2_mean(:,:,:) = &
     &                 overheadO3col_freq2_mean(:,:,:) / dcounter

               ! converting from mol/cm^2 to Dobson unit (DU)
               overheadO3col_freq2_mean (:,:,:) = &
                       overheadO3col_freq2_mean (:,:,:) / 2.690d+16
            end if

            if (pr_mass_freq2) mass_freq2_mean(:,:,:) = &
     &         mass_freq2_mean(:,:,:) / dcounter

            if (pr_grid_height_freq2) grid_height_freq2_mean (:,:,:) = &
               grid_height_freq2_mean (:,:,:) / dcounter

            if (pr_rel_hum_freq2) rel_hum_freq2_mean(:,:,:) = &
                rel_hum_freq2_mean (:,:,:) / dcounter

            if (pr_metwater_freq2) metwater_freq2_mean(:,:,:) = &
                metwater_freq2_mean (:,:,:) / dcounter

            if (pr_const_surface_freq2) constSurface_freq2_nc(:,:,:) = &
     &           constSurface_freq2_mean(:,:,:)

            if (pr_const_column_freq2) then
               if ((chem_mecha == 'troposphere') .or. (chem_mecha == 'strat_trop')  &
     &               .or. chem_mecha == 'strat_trop_aerosol' ) then
                  constColTrop_freq2_nc(:,:,:) = constColTrop_freq2_mean(:,:,:)

                  if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then
                     constColCombo_freq2_nc(:,:,:) = constColCombo_freq2_mean(:,:,:)
                  end if
               end if
            end if

            if (pr_tropopausePress_freq2) tropopausePress_freq2_nc(:,:) = &
     &          tropopausePress_freq2_mean (:,:)

            if (pr_psf_freq2) psf_freq2_nc(:,:) =   psf_freq2_mean (:,:)

            if (pr_kel_freq2) kel_freq2_nc(:,:,:) = kel_freq2_mean (:,:,:)

            if (pr_mass_freq2) mass_freq2_nc(:,:,:) = mass_freq2_mean (:,:,:)

            if (pr_overheadO3col_freq2) overheadO3col_freq2_nc(:,:,:) = &
               overheadO3col_freq2_mean (:,:,:)

            if (pr_potentialVorticity_freq2) &
     &         potentialVorticity_freq2_nc(:,:,:) = potentialVorticity_freq2_mean (:,:,:)

            if (pr_grid_height_freq2) grid_height_freq2_nc(:,:,:) = &
               grid_height_freq2_mean (:,:,:)

            if (pr_rel_hum_freq2) rel_hum_freq2_nc(:,:,:) = &
     &         rel_hum_freq2_mean(:,:,:)

            if (pr_metwater_freq2) metwater_freq2_nc(:,:,:) = &
               metwater_freq2_mean (:,:,:)

            ! Reset the counter
            counter = 0

         end if
!     ====
      else
!     ====
         ! --------------------------------------------------
         ! Fill output arrays with current "non-mean" values.
         ! --------------------------------------------------

         if (time_for_freq2) then
            if (pr_const_surface_freq2) constSurface_freq2_nc(:,:,:) = &
               constSurface(:,:,:)

            if (pr_const_column_freq2) then
               if ((chem_mecha == 'troposphere') .or. (chem_mecha == 'strat_trop')  &
     &               .or. chem_mecha == 'strat_trop_aerosol' ) then
                  constColTrop_freq2_nc(:,:,:) = constColTrop(:,:,:)

                  if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then
                     constColCombo_freq2_nc(:,:,:) = constColCombo(:,:,:)
                  end if
               end if
            end if

            if (pr_tropopausePress_freq2) &
     &         tropopausePress_freq2_nc(i1:i2,ju1:j2) = tropopausePress(i1:i2,ju1:j2)

            if (pr_psf_freq2) psf_freq2_nc(i1:i2,ju1:j2) = pctm2(i1:i2,ju1:j2)

            if (pr_potentialVorticity_freq2) &
               potentialVorticity_freq2_nc(i1:i2,ju1:j2,:) = &
                        potentialVorticity(i1:i2,ju1:j2,kx1:kx2)

            if (pr_kel_freq2) kel_freq2_nc(i1:i2,ju1:j2,:) = &
     &         kel(i1:i2,ju1:j2,kx1:kx2)

            if (pr_overheadO3col_freq2) overheadO3col_freq2_nc(:,:,:) = &
               overheadO3col(:,:,kx1:kx2)

            if (pr_mass_freq2) mass_freq2_nc(i1:i2,ju1:j2,:) = &
               mass(i1:i2,ju1:j2,kx1:kx2)

            if (pr_grid_height_freq2) grid_height_freq2_nc(i1:i2,ju1:j2,:) = &
     &         gridBoxHeight(i1:i2,ju1:j2,kx1:kx2)

            if (pr_rel_hum_freq2) rel_hum_freq2_nc(i1:i2,ju1:j2,:) = &
               relativeHumidity(i1:i2,ju1:j2,kx1:kx2)

            if (pr_metwater_freq2) metwater_freq2_nc(i1:i2,ju1:j2,:) = &
     &         humidity (i1:i2,ju1:j2,kx1:kx2) * MWTAIR / (MWTH2O * GPKG)

         end if
      end if

      if (pr_overheadO3col_freq2) deallocate(overheadO3col)
      if (pr_tropopausePress_freq2) deallocate(tropopausePress)
      if (pr_mass_freq2) deallocate(mass)
      if (pr_kel_freq2) deallocate(kel)
      if (pr_psf_freq2) deallocate(pctm2)
      if (pr_metwater_freq2) deallocate(humidity)
      if (pr_potentialVorticity_freq2) deallocate(potentialVorticity)
      if (pr_rel_hum_freq2) deallocate(relativeHumidity)
      if (pr_grid_height_freq2) deallocate(gridBoxHeight)

      return

      end subroutine Prep_Netcdf_Freq2
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Freq2
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Freq2 (nymd, nhms, gmi_sec)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: nymd, nhms
      real*8 , intent(in) :: gmi_sec
!
! !DESCRIPTION:
! This routine writes out variables that have at most 3 dimensions.
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_CONST_ColComb = 3001
      integer, parameter :: SG_CONST_ColTrop = 3002
      integer, parameter :: SG_CONST         = 3003
      integer, parameter :: SG_PSF           = 3004
      integer, parameter :: SG_CONST_SURF    = 3005
      integer, parameter :: SG_TROPP         = 3006
      integer, parameter :: SG_PV            = 3007
      integer, parameter :: SG_KEL           = 3008
      integer, parameter :: SG_OvO3col       = 3009
      integer, parameter :: SG_MASS          = 3010
      integer, parameter :: SG_GBH           = 3011
      integer, parameter :: SG_RelHum        = 3012
      integer, parameter :: SG_METWATER      = 3013
!
! !LOCAL VARIABLES:
      integer :: ic
      integer :: cnt2d (2), cnt3d (3), cnt4d (4)
      integer :: strt2d(2), strt3d(3), strt4d(4)
      integer :: hdr(NETCDF_HDF)
      real*8 , allocatable :: arrayGlob2D(:,:)
      real*8 , allocatable :: arrayGlob3D(:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Freq2 called by ', procID

      if (iAmRootProc) then
         strt2d(:) = (/ 1, rnum_out_freq2 /)
         strt3d(:) = (/ 1, 1, rnum_out_freq2 /)

         cnt2d (:) = (/ NETCDF_HDF, 1 /)
         cnt3d (:) = (/ numLon, numLat, 1 /)

         hdr(1) = Nint (gmi_sec)
         hdr(2) = nymd
         hdr(3) = nhms

         call Ncwr_2d_Int (hdr, ncid_freq2, hdr_var_name, strt2d, cnt2d)

         allocate(arrayGlob2D(numLon,numLat))
      end if

      if (pr_psf_freq2) then
         call subDomain2Global (arrayGlob2D, psf_freq2_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, &
     &           rootProc, procID, map1_u, numDomains, SG_PSF, commuWorld)

         if (iAmRootProc) call Ncwr_3d(arrayGlob2D(ix1:ix2, jx1:jx2), &
     &                            ncid_freq2, 'psf_freq2', strt3d, cnt3d)
      end if

      if (pr_tropopausePress_freq2) then
         call subDomain2Global (arrayGlob2D, tropopausePress_freq2_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, &
     &           rootProc, procID, map1_u, numDomains, SG_TROPP, commuWorld)

         if (iAmRootProc) call Ncwr_3d (arrayGlob2D(ix1:ix2, jx1:jx2), &
     &                   ncid_freq2, 'tropopausePress_freq2', strt3d, cnt3d)
      end if

      if (iAmRootProc) then
         allocate(arrayGlob3D(i1_gl:i2_gl,ju1_gl:j2_gl,freq2_species_num))

         cnt4d (:) = (/ numLon, numLat, 1,   1 /)
      end if

      if (pr_const_column_freq2) then

         if ((chem_mecha == 'troposphere') .or. (chem_mecha == 'strat_trop')  &
     &         .or. chem_mecha == 'strat_trop_aerosol' ) then

            call subDomain2Global (arrayGlob3D, constColTrop_freq2_nc, &
     &              i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, 1, &
     &              freq2_species_num, rootProc, procID, map1_u, numDomains, &
     &              SG_CONST_ColComb, commuWorld)

            if (iAmRootProc) then
               do ic = 1, freq2_species_num
                  strt4d(:) = (/ 1, 1, ic, rnum_out_freq2 /)
                  call Ncwr_4d (arrayGlob3D(ix1:ix2,jx1:jx2,ic), & 
     &                    ncid_freq2, 'constColTrop_freq2', strt4d, cnt4d)
               end do
            end if

            if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then
               call subDomain2Global (arrayGlob3D, constColCombo_freq2_nc, &
     &                 i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, 1, &
     &                 freq2_species_num, rootProc, procID, map1_u, numDomains, &
     &                 SG_CONST_ColTrop, commuWorld)

               if (iAmRootProc) then
                  do ic = 1, freq2_species_num
                     strt4d(:) = (/ 1, 1, ic, rnum_out_freq2 /)
                     call Ncwr_4d (arrayGlob3D(ix1:ix2,jx1:jx2,ic), & 
     &                    ncid_freq2, 'constColCombo_freq2', strt4d, cnt4d)
                  end do
               end if

            end if
         end if
      end if

      if (pr_const_surface_freq2) then
         call subDomain2Global (arrayGlob3D, constSurface_freq2_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, &
     &           1, freq2_species_num, rootProc, procID, map1_u, &
     &           numDomains, SG_CONST_SURF, commuWorld)

         if (iAmRootProc) then
            do ic = 1, freq2_species_num
               strt4d(:) = (/ 1, 1, ic, rnum_out_freq2 /)
               call Ncwr_4d (arrayGlob3D(ix1:ix2, jx1:jx2, ic), &
     &                   ncid_freq2,'constSurface_freq2', strt4d, cnt4d)
            end do
         end if
      end if

      if (iAmRootProc) then
         deallocate(arrayGlob3D)
         allocate(arrayGlob3D(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))

         strt4d(:) = (/ 1, 1, 1, rnum_out_freq2 /)
         cnt4d (:) = (/ numLon, numLat, numVert,  1 /)
      end if

      if (pr_potentialVorticity_freq2) then
         call subDomain2Global (arrayGlob3D, potentialVorticity_freq2_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_PV, commuWorld)

         if (iAmRootProc) call Ncwr_4d (arrayGlob3D(ix1:ix2, jx1:jx2,:),&
     &            ncid_freq2, 'potentialVorticity_freq2', strt4d, cnt4d)
      end if

      if (pr_kel_freq2) then
         call subDomain2Global (arrayGlob3D, kel_freq2_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_KEL, commuWorld)

         if (iAmRootProc) call Ncwr_4d (arrayGlob3D(ix1:ix2, jx1:jx2,:),&
     &            ncid_freq2, 'kel_freq2', strt4d, cnt4d)
      end if

      if (pr_overheadO3col_freq2) then
         call subDomain2Global (arrayGlob3D, overheadO3col_freq2_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_OvO3col, commuWorld)

         if (iAmRootProc) call Ncwr_4d (arrayGlob3D(ix1:ix2, jx1:jx2,:),&
     &            ncid_freq2, 'overheadO3col_freq2', strt4d, cnt4d)
      end if

      if (pr_mass_freq2) then
         call subDomain2Global (arrayGlob3D, mass_freq2_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_MASS, commuWorld)

         if (iAmRootProc) call Ncwr_4d (arrayGlob3D(ix1:ix2, jx1:jx2,:),&
     &            ncid_freq2, 'mass_freq2', strt4d, cnt4d)
      end if

      if (pr_grid_height_freq2) then
         call subDomain2Global (arrayGlob3D, grid_height_freq2_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_GBH, commuWorld)

         if (iAmRootProc) call Ncwr_4d (arrayGlob3D(ix1:ix2, jx1:jx2,:),&
     &            ncid_freq2, 'grid_height_freq2', strt4d, cnt4d)
      end if

      if (pr_rel_hum_freq2) then
         call subDomain2Global (arrayGlob3D, rel_hum_freq2_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_RelHum, commuWorld)

         if (iAmRootProc) call Ncwr_4d (arrayGlob3D(ix1:ix2, jx1:jx2,:),&
     &            ncid_freq2, 'rel_hum_freq2', strt4d, cnt4d)
      end if

      if (pr_metwater_freq2) then
         call subDomain2Global (arrayGlob3D, metwater_freq2_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_METWATER, commuWorld)

         if (iAmRootProc) call Ncwr_4d (arrayGlob3D(ix1:ix2, jx1:jx2,:),&
     &            ncid_freq2, 'metwater_freq2', strt4d, cnt4d)
      endif

      if (iAmRootProc) call Ncdo_Sync (ncid_freq2)

      return

      end subroutine Write_Netcdf_Freq2
!EOC
!-----------------------------------------------------------------------------
!BOP
! 
! !IROUTINE: Write_Netcdf_Hdr_Freq2
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Hdr_Freq2 (gmiDomain, metFields, const_labels)

      use m_ncGeneralOpsOutput, only : WriteNetcdfHdrGen

      implicit none
!
! !INPUT PARAMETERS:
      character (len=MAX_LENGTH_SPECIES_NAME), intent(in) :: const_labels(1:)
      type(t_metFields), intent(in) :: metFields
      type(t_gmiDomain), intent(in) :: gmiDomain
!
! !DESCRIPTION:
! This routine creates some header information for the netCDF output
! file and writes it out.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_SPECIES_NAME) :: spclab(freq2_species_num)
      character (len=50) :: metnam(1)
      integer :: ic, icx
      integer :: cnt1d (1), cnt2d (2)
      integer :: strt1d(1), strt2d(2)
      real*8  :: prsdat(k1:k2)
      real*8  :: spcdat(freq2_species_num)
      real*8 , allocatable :: locLatDeg(:)
      real*8 , allocatable :: locLonDeg(:)
      real*8 , allocatable :: locMCOR(:,:)
      real*8 , allocatable :: mcor(:,:)
      real*8 , allocatable :: londeg(:), latdeg(:)
      real*8 , allocatable :: ai(:), bi(:)
      real*8 , allocatable :: am(:), bm(:)
      real*8 :: pt
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Hdr_Freq2 called by ', procID

      allocate(latdeg(ju1_gl:j2_gl))
      call Get_latdeg(gmiDomain, latdeg)

      allocate(londeg(i1_gl:i2_gl))
      call Get_londeg(gmiDomain, londeg)

      allocate(mcor(i1_gl:i2_gl,ju1_gl:j2_gl))
      call Get_mcorGlob(gmiDomain, mcor)

      allocate(locLonDeg(1:numLon))
      allocate(locLatDeg(1:numLat))

      locLonDeg(1:numLon) = londeg(ix1:ix2)
      locLatDeg(1:numLat) = latdeg(jx1:jx2)

!     =========================
      call WriteNetcdfHdrGen  &
!     =========================
     &  (ncid_freq2, locLatDeg, locLonDeg, pr_diag, procID, numLon, numLat, &
     &   hdf_dim_name, lat_dim_name, lon_dim_name)

!     ---------
!     Pressure.
!     ---------

      allocate(ai(k1-1:k2))
      allocate(bi(k1-1:k2))
      allocate(am(k1:k2))
      allocate(bm(k1:k2))

      call Get_pt(metFields, pt)
      call Get_ai(metFields, ai)
      call Get_bi(metFields, bi)
      call Get_am(metFields, am)
      call Get_bm(metFields, bm)

      prsdat(:) = (am(:) * pt)  + (bm(:) * 1000.0d0)

      strt1d(1) = 1
      cnt1d (1) = numVert

      call Ncwr_1d (prsdat(kx1:kx2), ncid_freq2, prs_dim_name, strt1d, cnt1d)
      call Ncwr_1d (am(kx1:kx2), ncid_freq2, 'am', strt1d, cnt1d)
      call Ncwr_1d (bm(kx1:kx2), ncid_freq2, 'bm', strt1d, cnt1d)

      cnt1d(1) = numvert + 1

      call Ncwr_1d (ai(kx1-1:kx2), ncid_freq2, 'ai', strt1d, cnt1d)
      call Ncwr_1d (bi(kx1-1:kx2), ncid_freq2, 'bi', strt1d, cnt1d)

      call Ncwr_Scal (pt, ncid_freq2, 'pt')

!     ------------
!     Species map.
!     ------------

      icx = 0

      do ic = 1, numSpecies

        if (freq2_species(ic) /= 0) then
          icx = icx + 1
          spcdat(icx) = ic
        end if

      end do

      strt1d(1) = 1
      cnt1d (1) = freq2_species_num

      if (freq2_species_num > 0) then
        call Ncwr_1d (spcdat, ncid_freq2, 'freq2_species', strt1d, cnt1d)
      end if

!     ------------
!     Freq2 labels.
!     ------------

      icx = 0

      do ic = 1, numSpecies
        if (freq2_species(ic) /= 0) then
          icx = icx + 1
          spclab(icx) = const_labels(ic)
        end if
      end do

      strt2d(:) = (/  1, 1 /)
      cnt2d (:) = (/ MAX_LENGTH_SPECIES_NAME, freq2_species_num /)

      if (freq2_species_num > 0) then
         call Ncwr_2d_Char (spclab, ncid_freq2, 'freq2_labels', strt2d, cnt2d)
      end if

!     --------------
!     Grid box area.
!     --------------

      strt2d(:) = (/ 1, 1 /)
      cnt2d (:) = (/ numLon, numLat /)

      allocate(locMCOR(ix1:ix2, jx1:jx2))

      locMCOR(ix1:ix2,jx1:jx2) = mcor(ix1:ix2, jx1:jx2)
      call Ncwr_2d (locMCOR, ncid_freq2, 'mcor', strt2d, cnt2d)

      deallocate(locMCOR)
      deallocate(locLonDeg)
      deallocate(locLatDeg)

      return

      end subroutine Write_Netcdf_Hdr_Freq2
!EOC
!-----------------------------------------------------------------------------

      end module GmiControlFreq2_mod
