!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlConcentration_mod
!
! !INTERFACE:
!
module GmiControlConcentration_mod
!
! !USES:
      use  GmiSub2Glob_mod      , only : subDomain2Global
      use  GmiMessagePassing_mod, only : synchronizeGroup
      use GmiFileOperations_mod , only : makeOutfileName
      use m_netcdf_io_close , only : Nccl, Nccl_Noerr
      use m_netcdf_io_write , only : Ncwr_2d_Int, Ncwr_5d, Ncwr_3d, Ncwr_2d, &
     &       Ncwr_4d, Ncwr_3d_Int
      use m_netcdf_io_write , only : Ncwr_1d, Ncwr_Scal, Ncwr_2d_Char
      use m_netcdf_io_create, only : Ncdo_Sync, Nccr_Wr
      use m_netcdf_io_define, only : NcDef_glob_attributes, NcDef_dimension
      use m_netcdf_io_define, only : NcDef_var_attributes, NcDef_variable
      use m_netcdf_io_define, only : NcSetFill, NcEnd_def
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration,        &
     &       Get_concentration, Get_concentrationSurf, Get_const_var_name,     &
     &       Get_concentrationColTrop, Get_concentrationColCombo
      use GmiChemistryMethod_mod, only : t_Chemistry, Get_overheadO3col,       &
     &       Get_const_labels, Get_decay_3d_out, Set_decay_3d_out
      use GmiEmissionMethod_mod, only : t_Emission, Get_ndust, Get_naero, &
     &       Get_lightning_opt, Get_emiss_aero_opt, Get_emiss_dust_opt,   &
     &       Get_emiss_map_dust, Get_emiss_map_aero, Get_do_gcr
      use GmiDepositionMethod_mod      , only : t_Deposition
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_rootProc,        &
     &       Get_procID, Get_iAmRootProc,         &
     &       Get_numDomains, Get_communicatorWorld, Get_map1_u!,&
!     &       Get_eastWestFlag, Get_northSouthPole
      use GmiDomainDecomposition_mod, only : Get_mcorGlob, Get_latdeg, Get_londeg
      use GmiGrid_mod, only : t_gmiGrid, Get_numSpecies, Get_i1, Get_i2,       &
     &       Get_ju1, Get_j2, Get_k1, Get_k2, Get_i1_gl, Get_i2_gl, Get_ju1_gl,&
     &       Get_j2_gl, Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jhi_gl, &
     &       Get_ilo, Get_ihi, Get_julo, Get_jhi
      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Get_pr_diag,         &
     &       Get_num_drydep_outrecs, Get_num_const_outrecs, Get_pr_const,      &
     &       Get_num_wetdep_outrecs, Get_num_emiss_outrecs, Get_problem_name,  &
     &       Get_hdr_var_name, Get_hdf_dim_name, Get_rec_dim_name, Get_pr_psf, &
     &       Get_lat_dim_name, Get_lon_dim_name, Get_prs_dim_name, Get_pr_kel, &
     &       Get_spc_dim_name, Get_constOutputFrequency, Get_do_mean, Get_k1r_gl,      &
     &       Get_k2r_gl, Get_pr_sulf_src, Get_pr_surf_emiss, Get_pr_emiss_3d,  &
     &       Get_pr_dry_depos, Get_pr_wet_depos, Get_do_aerocom, Get_pr_mass,  &
     &       Get_pr_grid_height, Get_outmain_name, Get_pr_overheadO3col,       &
     &       Get_pr_tropopausePress, Get_pr_potentialVorticity, Get_pr_decay,  &
     &       Get_pr_metwater, Get_pr_relHumidity, Get_pr_const_surface,        &
     &       Get_pr_emiss_all, Get_pr_drydep_all, Get_pr_wetdep_all,           &
     &       Get_const_outrec_map, Get_emiss_outrec_map, Get_drydep_outrec_map,&
     &       Get_wetdep_outrec_map, Get_pr_const_column, Get_pr_level_all

      use GmiTimeControl_mod, only : t_GmiClock, GmiSplitDateTime, isOutTime, &
     &       Get_curGmiDate, Get_curGmiTime, Get_gmiSeconds, Get_gmiTimeStep, &
     &       Get_numTimeSteps

      use GmiMetFieldsControl_mod, only : t_metFields, Get_pt, Get_met_opt,    &
     &       Get_metdata_name, Get_mdt, Get_ai, Get_bi, Get_am, Get_bm,        &
     &       Get_mass, Get_tropopausePress, Get_potentialVorticity,            &
     &       Get_gridBoxHeight, Get_relativeHumidity, Get_kel, Get_pctm2,      &
     &       Get_humidity, Get_cmi_flags, Get_cmf
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: controlOutputConcentration, initializeOutputConcentration
      public  :: finalizeOutputConcentration

#     include "gmi_phys_constants.h"
#     include "gmi_diag_constants_llnl.h"
#     include "GmiParameters.h"
#     include "setkin_par.h"
!
! !DESCRIPTION:
! Contains the necessary routines for the species concentration outputs. 
! Three routines are visible from the outside:
! %
! \begin{description}
! \item[initializeOutputConcentration:] called once in the initialization steps to
!      open the netCDF file, define variables in the file, write out header
!      information and allocate variables.
! \item[controlOutputConcentration:] called at the end of each model time step to
!      update variables and if needed for communications and writing out
!      data in the file.
! \item[finalizeOutputConcentration:] called at the end of the integration to
!      deallocate all the variables and close the netCDF file.
! \end{description}
! %
! The above routines have as arguments derived types.
! %
! This module is self-contained. All the main routines are called by all
! the processors (root included) but each processor does not execute
! all portions of the code.
!
! There are two categories of arrays: (1) the ones manipulate by all the
! processors, and (2) those handle by worker processors only.
!
      ! Variables manipulate by all processors
      real*8 , pointer, save :: const_nc               (:,:,:) => null()
      real*8 , pointer, save :: psf_nc                   (:,:) => null()
      real*8 , pointer, save :: kel_nc                 (:,:,:) => null()
      real*8 , pointer, save :: cmf_nc                 (:,:,:) => null()
      real*8 , pointer, save :: ddep_nc                (:,:,:) => null()
      real*8 , pointer, save :: wdep_nc                (:,:,:) => null()
      real*8 , pointer, save :: mass_nc                (:,:,:) => null()
      real*8 , pointer, save :: dtrn_nc                (:,:,:) => null()
      real*8 , pointer, save :: cldmas0_nc               (:,:) => null()
      real*8 , pointer, save :: metwater_nc            (:,:,:) => null()
      real*8 , pointer, save :: cmi_flags_nc             (:,:) => null()
      real*8 , pointer, save :: cmi_flags1_nc            (:,:) => null()
      real*8 , pointer, save :: flashrate_nc             (:,:) => null()
      real*8 , pointer, save :: lightning_nc           (:,:,:) => null()
      real*8 , pointer, save :: GCR_NOx_nc             (:,:,:) => null()
      real*8 , pointer, save :: semiss_out_nc          (:,:,:) => null()
      real*8 , pointer, save :: semiss_out2_nc         (:,:,:) => null()
      real*8 , pointer, save :: relHumidity_nc         (:,:,:) => null()
      real*8 , pointer, save :: emiss_3d_out_nc        (:,:,:) => null()
      real*8 , pointer, save :: overheadO3col_nc       (:,:,:) => null()
      real*8 , pointer, save :: decay_3d_out_nc        (:,:,:) => null()
      real*8 , pointer, save :: aerosolEmiss3D_nc      (:,:,:) => null()
      real*8 , pointer, save :: gridbox_height_nc      (:,:,:) => null()
      real*8 , pointer, save :: tropopausePress_nc       (:,:) => null()
      real*8 , pointer, save :: aerosolSurfEmiss_nc    (:,:,:) => null()
      real*8 , pointer, save :: potentialVorticity_nc  (:,:,:) => null()
!      real*8 , pointer, save :: aerosolEmiss3D       (:,:,:,:) => null()

      real*8 , pointer, save :: dms_no3                (:,:,:) => null()
      real*8 , pointer, save :: dms_oh                 (:,:,:) => null()

      real*8 , pointer, save :: so2_h2o2               (:,:,:) => null()
      real*8 , pointer, save :: so2_o3                 (:,:,:) => null()
      real*8 , pointer, save :: so2_oh                 (:,:,:) => null()
!     
      ! Variables manipulate by worker processors only         
      real*8 , pointer, save :: psf_mean                 (:,:) => null()
      real*8 , pointer, save :: kel_mean               (:,:,:) => null()
      real*8 , pointer, save :: mass_mean              (:,:,:) => null()
      real*8 , pointer, save :: const_mean           (:,:,:,:) => null()
      real*8 , pointer, save :: relHumidity_mean       (:,:,:) => null()
      real*8 , pointer, save :: metwater_mean          (:,:,:) => null()
      real*8 , pointer, save :: grid_height_mean    (:,:,:) => null()
      real*8 , pointer, save :: overheadO3col_mean     (:,:,:) => null()
      real*8 , pointer, save :: tropopausePress_mean     (:,:) => null()
      real*8 , pointer, save :: potentialVorticity_mean(:,:,:) => null()
!
      ! netCDF file identifier
      integer,          save :: ncid_const

      ! Counter for the number of records
      integer,          save :: rnum_out_const

      ! Grid information of the variables to be written out
      integer,          save :: i1, i2     ! longitude
      integer,          save :: i1_gl, i2_gl
      integer,          save :: ilo_gl, ihi_gl, ilo, ihi
      integer,          save :: ju1        ! latitude
      integer,          save :: j2
      integer,          save :: ju1_gl, j2_gl
      integer,          save :: julo_gl, jhi_gl, julo, jhi
      integer,          save :: kx1, kx2   ! vertical
      integer,          save :: k1, k2
      integer,          save :: numLon     ! number of longitudes
      integer,          save :: numLat     ! number of latitudes
      integer,          save :: numVert    ! number of vertical levels
      integer,          save :: numSpecies ! number of species

      integer,          save :: numAero          ! # aerosol species
      integer,          save :: numDust          ! # dust    species
      integer,          save :: num_const_outrecs  ! # species to output
      integer,          save :: num_emiss_outrecs  ! # species selected for surface emiss output
      integer,          save :: num_drydep_outrecs ! # species selected for dry dep output
      integer,          save :: num_wetdep_outrecs ! # species selected for wet dep output
      
      logical,          save :: iAmRootProc
      integer,          save :: procID, rootProc
      integer,          save :: numDomains
      integer,          save :: commuWorld
      integer, pointer, save :: map1_u       (:,:,:) => null()
!      integer, pointer, save :: ewflag           (:) => null()
!      integer, pointer, save :: nspole           (:) => null()

      integer,          save :: nhdf
      logical,          save :: pr_diag
      logical,            save :: pr_const, pr_kel, pr_psf, do_mean
      logical,            save :: pr_overheadO3col      
      logical,            save :: pr_decay
      logical,            save :: pr_tropopausePress    
      logical,            save :: pr_potentialVorticity 
      logical,            save :: pr_mass               
      logical,            save :: pr_grid_height        
      logical,            save :: pr_relHumidity        
      logical,            save :: pr_metwater           
      logical,            save :: pr_dry_depos          
      logical,            save :: pr_wet_depos          
      logical,            save :: pr_surf_emiss         
      logical,            save :: pr_emiss_3d           
      logical,            save :: pr_emiss_all          
      logical,            save :: pr_drydep_all         
      logical,            save :: pr_wetdep_all         
      logical,            save :: pr_sulf_src           
      logical,            save :: do_aerocom
      logical,            save :: do_gcr
      integer,            save :: const_outrec_map(MAX_NUM_CONST_GIO)
      integer,            save :: emiss_outrec_map(MAX_NUM_CONST_GIO)
      integer,            save :: drydep_outrec_map(MAX_NUM_CONST_GIO)
      integer,            save :: wetdep_outrec_map(MAX_NUM_CONST_GIO)
      real*8 ,            save :: constOutputFrequency    
      logical,            save :: pr_const_column
      logical,            save :: pr_const_surface
      logical,            save :: pr_level_all
      character (len=MAX_LENGTH_VAR_NAME), save :: rec_dim_name, prs_dim_name, spc_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: hdr_var_name
      character (len=MAX_LENGTH_VAR_NAME), save :: hdf_dim_name, lon_dim_name, lat_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: const_var_name
      integer           , save :: met_opt, lightning_opt
      integer           , save :: emiss_aero_opt, emiss_dust_opt
      real*8            , save :: mdt
!
! !AUTHOR:
! Jules Kouatchou, Jules.Kouatchou-1@nasa.gov
!
! !REVISION HISTORY:
! 12Nov2007: Original code.
!EOP
!-----------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initializeOutputConcentration
!
! !INTERFACE:
!
      subroutine initializeOutputConcentration(gmiGrid, gmiDomain, Diagnostics, &
     &                    Chemistry, SpeciesConcentration, Emission, metFields)
!
! !USES:
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_gmiGrid  ), intent(in) :: gmiGrid
      type(t_gmiDomain), intent(in) :: gmiDomain
      type(t_metFields), intent(in) :: metFields
      type(t_Emission), intent(in) :: Emission
      type(t_Chemistry), intent(in) :: Chemistry
      type(t_Diagnostics), intent(in) :: Diagnostics
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !DESCRIPTION:
! The routines (1) opens a netCDF file, (2) defines variables in the file,
! (3) writes out header data in the file, and (4) allocates varaibles.
!
! !LOCAL VARIABLES:
      character (len=80), save :: outmain_name
      character (len=75)  :: err_msg
      integer :: in1, k1r_gl, k2r_gl
!      integer, allocatable :: ewflag_org(:)
!      integer, allocatable :: nspole_org(:)
      character (len=MAX_LENGTH_FILE_NAME) :: fname
      character (len=MAX_LENGTH_FILE_NAME) :: problem_name
      integer :: ndust, naero
      character (len=MAX_LENGTH_SPECIES_NAME), pointer :: const_labels(:)
!EOP
!-----------------------------------------------------------------------------
!BOC

      !#####################################################
      ! Initialization of variables used in all the routines
      !#####################################################

      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write (6,*) 'initializeOutputConcentration called by ', procID

      call Get_ndust(Emission, ndust)
      call Get_naero(Emission, naero)
      call Get_emiss_aero_opt(Emission, emiss_aero_opt)
      call Get_emiss_dust_opt(Emission, emiss_dust_opt)
      call Get_lightning_opt(Emission, lightning_opt)
      call Get_do_gcr(Emission, do_gcr)

      numAero          = naero
      numDust          = ndust

      call Get_pr_kel(Diagnostics, pr_kel)
      call Get_pr_psf(Diagnostics, pr_psf)
      call Get_k1r_gl(Diagnostics, k1r_gl)
      call Get_k2r_gl(Diagnostics, k2r_gl)
      call Get_do_mean(Diagnostics, do_mean)
      call Get_pr_mass(Diagnostics, pr_mass)
      call Get_pr_const(Diagnostics, pr_const)
      call Get_do_aerocom(Diagnostics, do_aerocom)
      call Get_pr_metwater(Diagnostics, pr_metwater)
      call Get_pr_emiss_3d(Diagnostics, pr_emiss_3d)
      call Get_pr_sulf_src(Diagnostics, pr_sulf_src)
      call Get_pr_emiss_all(Diagnostics, pr_emiss_all)
      call Get_pr_level_all(Diagnostics, pr_level_all)
      call Get_hdr_var_name(Diagnostics, hdr_var_name)
      call Get_spc_dim_name(Diagnostics, spc_dim_name)
      call Get_rec_dim_name(Diagnostics, rec_dim_name)
      call Get_prs_dim_name(Diagnostics, prs_dim_name)
      call Get_hdf_dim_name(Diagnostics, hdf_dim_name)
      call Get_lat_dim_name(Diagnostics, lat_dim_name)
      call Get_lon_dim_name(Diagnostics, lon_dim_name)
      call Get_constOutputFrequency(Diagnostics, constOutputFrequency)
      call Get_pr_dry_depos(Diagnostics, pr_dry_depos)
      call Get_pr_wet_depos(Diagnostics, pr_wet_depos)
      call Get_pr_surf_emiss(Diagnostics, pr_surf_emiss)
      call Get_pr_drydep_all(Diagnostics, pr_drydep_all)
      call Get_pr_wetdep_all(Diagnostics, pr_wetdep_all)
      call Get_pr_grid_height(Diagnostics, pr_grid_height)
      call Get_pr_relHumidity(Diagnostics, pr_relHumidity)
      call Get_pr_const_column(Diagnostics, pr_const_column)
      call Get_pr_overheadO3col(Diagnostics, pr_overheadO3col)
      call Get_pr_decay(Diagnostics, pr_decay)
      call Get_pr_const_surface(Diagnostics, pr_const_surface)
      call Get_const_outrec_map(Diagnostics, const_outrec_map)
      call Get_emiss_outrec_map(Diagnostics, emiss_outrec_map)
      call Get_drydep_outrec_map(Diagnostics, drydep_outrec_map)
      call Get_wetdep_outrec_map(Diagnostics, wetdep_outrec_map)
      call Get_num_const_outrecs(Diagnostics, num_const_outrecs)
      call Get_num_emiss_outrecs(Diagnostics, num_emiss_outrecs)
      call Get_num_drydep_outrecs(Diagnostics, num_drydep_outrecs)
      call Get_num_wetdep_outrecs(Diagnostics, num_wetdep_outrecs)
      call Get_pr_tropopausePress(Diagnostics, pr_tropopausePress)
      call Get_pr_potentialVorticity(Diagnostics, pr_potentialVorticity)

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

      call Get_const_var_name(SpeciesConcentration, const_var_name)

      call Get_met_opt(metFields, met_opt)
      call Get_mdt(metFields, mdt)

      allocate(map1_u(2, 2, numDomains))
      call Get_map1_u (gmiDomain, map1_u)

      kx1     = k1r_gl
      kx2     = k2r_gl

      numLon  = i2_gl - i1_gl + 1
      numLat  = j2_gl - ju1_gl + 1
      numVert = kx2 - kx1 + 1

      nhdf = NETCDF_HDF

      allocate(const_labels(numSpecies))

      if (iAmRootProc) then

         !###########################
         ! Initialize the output file
         !###########################

         ! Determine the file name
         call Get_outmain_name(Diagnostics, outmain_name)
         call Get_problem_name(Diagnostics, problem_name)

         in1 = Len_Trim (outmain_name)
         call makeOutfileName (fname, '.' // outmain_name(1:in1) // '.nc', &
     &                         problem_name)

         ! Create the file and assign a file identifier
         call Nccr_Wr (ncid_const, fname)

         ! Define the variables in the file
         call Define_Netcdf_Out_Const(const_var_name)

         ! Write header data
         call Get_const_labels(Chemistry, const_labels)

         call Write_Netcdf_Hdr_Const (gmiDomain, Emission, metFields, &
     &                                const_labels)

         call Ncdo_Sync (ncid_const)

         ! Initialize the counter for the number of records
         rnum_out_const = 1
      end if

      !########################
      ! Allocation of variables
      !########################
      call allocateVariablesConcentration()

      return

      end subroutine initializeOutputConcentration
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: controlOutConcentration
!
! !INTERFACE:
!
      subroutine controlOutputConcentration (last_tstp, Chemistry, Emission, &
     &                  Deposition, SpeciesConcentration, gmiClock, metFields)
!
! !USES:
      use GmiTimeControl_mod, only : GmiSplitDateTime

      implicit none
!
! !INPUT PARAMETERS:
      logical                     , intent(in) :: last_tstp
      type(t_gmiClock            ), intent(in) :: gmiClock            
      type(t_metFields           ), intent(in) :: metFields            
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_Emission)  , intent(inOut) :: Emission  
      type(t_Chemistry) , intent(inOut) :: Chemistry
      type(t_Deposition), intent(inOut) :: Deposition
!
! !DESCRIPTION:
! Controls the species concentration output file. It is called at each time step.
!
! !LOCAL VARIABLES:
      logical :: time_for_nc
      integer :: day, month, idumyear, ndt, num_time_steps, nhms, nymd
      real*8  :: gmi_sec, tstep
      integer, save :: month_save = -999
      logical, save :: printed_on_this_day = .false.
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'controlOutputConcentration called by ', procID

      call Get_curGmiDate(gmiClock, nymd)
      call Get_curGmitime(gmiClock, nhms)
      call Get_numTimeSteps(gmiClock, num_time_steps)
      call Get_gmiSeconds(gmiClock, gmi_sec) 
      call Get_gmiTimeStep(gmiClock, tstep)

      ndt = Nint(tstep)

      call GmiSplitDateTime (nymd, idumyear, month, day)

      if (month_save == -999) month_save = month

      time_for_nc = .false.

      call isOutTime(time_for_nc, printed_on_this_day, &
     &      month_save, month, day, nhms, ndt, gmi_sec, constOutputFrequency)

      month_save = month

      if (last_tstp .and. (.not. time_for_nc)) then
        if (Mod (num_time_steps, Nint (mdt) / ndt) == 0) then

          time_for_nc = .true.

        end if
      end if

      if (time_for_nc .or. do_mean) then
         !====================
         call Prep_Netcdf_Const &
         !====================
     &       (time_for_nc, Chemistry, Deposition, Emission, &
     &        SpeciesConcentration, metFields)
      endif

!     ================
      if (time_for_nc) then
!     ================

!        =====================
         call Write_Netcdf_Const (Chemistry, Emission, nymd, nhms, gmi_sec)
!        =====================

!       ======================
        call bufferOutput_Const (const_var_name, Emission, Chemistry, SpeciesConcentration)
!       ======================

        if (iAmRootProc) then
            rnum_out_const = rnum_out_const + 1
        end if

!     ======
      end if
!     ======

      return

      end subroutine controlOutputConcentration
!EOC
!-----------------------------------------------------------------------------
!BOP
!     
! !IROUTINES: finalizeOutputConcentration
!
! !INTERFACE:
!
      subroutine finalizeOutputConcentration()
!
      implicit none
! 
! !DESCRIPTION:
! Deallocates variables necessary to produce species concentration outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'finalizeOutputConcentration called by ', procID

      if (pr_const) then

         deallocate (const_nc)

         if (pr_psf) then
            deallocate (psf_nc )
         end if 

         if (pr_tropopausePress) then
            deallocate (tropopausePress_nc )
         end if

         if (pr_potentialVorticity) then
            deallocate (potentialVorticity_nc)
         end if

         if (pr_overheadO3col) then
            deallocate (overheadO3col_nc )
         end if

         if (pr_decay) then
            deallocate (decay_3d_out_nc)
         endif

         if (pr_kel) then
            deallocate (kel_nc )
         end if

         if (pr_mass) then
            deallocate (mass_nc)
         end if

         if (pr_relHumidity) then
            deallocate (relHumidity_nc)
         end if

         if (pr_grid_height) then
            deallocate (gridbox_height_nc)
         end if

         if (pr_emiss_3d) then
            deallocate (emiss_3d_out_nc)
            if ((emiss_dust_opt > 0) .or. (emiss_aero_opt > 0)) then
               deallocate (aerosolEmiss3D_nc)
            end if
         end if

         if (pr_metwater) then
            deallocate (metwater_nc)
         end if

         if (pr_dry_depos) then
            if (pr_drydep_all) then
               deallocate (ddep_nc)
            else
               deallocate (ddep_nc)
            end if
         end if
         if (pr_wet_depos) then
            if (pr_wetdep_all) then
               deallocate (wdep_nc)
            else
               deallocate (wdep_nc)
            endif
         end if

         if (pr_surf_emiss) then
            deallocate (semiss_out_nc)
            deallocate (semiss_out2_nc)

            if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
               deallocate (aerosolSurfEmiss_nc)
            end if
         end if
         if (do_gcr) then
            deallocate(GCR_NOx_nc)
         endif

         if (lightning_opt == 1) then
            deallocate(lightning_nc)
            deallocate(flashrate_nc)
            deallocate(cmf_nc)
            deallocate(cmi_flags_nc)
            deallocate(cmi_flags1_nc)
            deallocate(cldmas0_nc)
         endif

         if (pr_sulf_src) then
            deallocate (dms_no3)
            deallocate (dms_oh )
            deallocate (so2_h2o2)
            deallocate (so2_o3  )
            deallocate (so2_oh  )
         end if

         if (do_mean) then
            deallocate (const_mean)

            if (pr_psf) then
               deallocate (psf_mean )
            end if
  
            if (pr_tropopausePress) then
               deallocate (tropopausePress_mean )
            end if
  
            if (pr_potentialVorticity) then
               deallocate (potentialVorticity_mean )
            end if
  
            if (pr_overheadO3col) then
               deallocate (overheadO3col_mean )
            end if
  
            if (pr_kel) then
               deallocate (kel_mean )
            end if
  
            if (pr_mass) then
               deallocate (mass_mean)
            end if
  
            if (pr_relHumidity) then
               deallocate (relHumidity_mean)
            end if
  
            if (pr_grid_height) then
               deallocate (grid_height_mean)
            end if
  
            if (pr_metwater) then
               deallocate (metwater_mean)
            end if
         end if

         if (iAmRootProc) call Nccl_Noerr (ncid_const)
      end if

      return

      end subroutine finalizeOutputConcentration
!EOC
!-----------------------------------------------------------------------------
!BOP
!     
! !IROUTINES: allocateVariablesConcentration
!
! !INTERFACE:
!
      subroutine allocateVariablesConcentration()
!
      implicit none
! 
! !DESCRIPTION:
! Allocates variables necessary to produce species concentration outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'allocateVariablesConcentration called by ', procID

      if (pr_const) then

         Allocate (const_nc(i1:i2, ju1:j2, kx1:kx2))
         const_nc = 0.0d0

         if (pr_psf) then
            Allocate (psf_nc (i1:i2, ju1:j2))
            psf_nc = 0.0d0
         end if 

         if (pr_tropopausePress) then
            Allocate (tropopausePress_nc (i1:i2, ju1:j2))
            tropopausePress_nc = 0.0d0
         end if

         if (pr_potentialVorticity) then
            Allocate (potentialVorticity_nc(i1:i2, ju1:j2, kx1:kx2))
            potentialVorticity_nc = 0.0d0 
         end if

         if (pr_overheadO3col) then
            Allocate (overheadO3col_nc (i1:i2, ju1:j2, kx1:kx2))
            overheadO3col_nc = 0.0d0
         end if

         if (pr_decay) then
           Allocate (decay_3d_out_nc(i1:i2, ju1:j2, kx1:kx2))
           decay_3d_out_nc = 0.0d0
         endif

         if (pr_kel) then
            Allocate (kel_nc (i1:i2, ju1:j2, kx1:kx2))
            kel_nc = 0.0d0
         end if

         if (pr_mass) then
            Allocate (mass_nc(i1:i2, ju1:j2, kx1:kx2))
            mass_nc = 0.0d0
         end if

         if (pr_relHumidity) then
            Allocate (relHumidity_nc(i1:i2, ju1:j2, kx1:kx2))
            relHumidity_nc = 0.0d0
         end if

         if (pr_grid_height) then
            Allocate (gridbox_height_nc(i1:i2, ju1:j2, kx1:kx2))
            gridbox_height_nc = 0.0d0
         end if

         if (pr_emiss_3d) then
            !Allocate (emiss_3d_out(i1:i2, ju1:j2, kx1:kx2, numSpecies))
            !emiss_3d_out = 0.0d0

            Allocate (emiss_3d_out_nc(i1:i2, ju1:j2, kx1:kx2))
            emiss_3d_out_nc = 0.0d0
            if ((emiss_dust_opt > 0) .or. (emiss_aero_opt > 0)) then
               Allocate (aerosolEmiss3D_nc(i1:i2, ju1:j2, kx1:kx2))
               aerosolEmiss3D_nc = 0.0d0

               !Allocate(aerosolEmiss3D(i1:i2,ju1:j2,kx1:kx2,5))
               !aerosolEmiss3D = 0.0d0
            end if
         end if

         if (pr_metwater) then
            Allocate (metwater_nc(i1:i2, ju1:j2, kx1:kx2))
            metwater_nc = 0.0d0
         end if

         if (pr_dry_depos) then
            if (pr_drydep_all) then
               Allocate (ddep_nc(i1:i2, ju1:j2, 1:numSpecies))
            else
               Allocate (ddep_nc(i1:i2, ju1:j2, 1:num_drydep_outrecs))
            end if
            ddep_nc = 0.0d0
         end if
         if (pr_wet_depos) then
            if (pr_wetdep_all) then
               Allocate (wdep_nc(i1:i2, ju1:j2, 1:numSpecies))
            else
               Allocate (wdep_nc(i1:i2, ju1:j2, 1:num_wetdep_outrecs))
            endif
            wdep_nc = 0.0d0
         end if

         if (pr_surf_emiss) then
            if (pr_emiss_all) then
               Allocate (semiss_out_nc(i1:i2, ju1:j2, 1:numSpecies))
            else
               Allocate (semiss_out_nc(i1:i2, ju1:j2, 1:num_emiss_outrecs))
            endif
            semiss_out_nc = 0.0d0
            Allocate (semiss_out2_nc(i1:i2, ju1:j2, 1:6))
            semiss_out2_nc= 0.0d0

            if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
               Allocate (aerosolSurfEmiss_nc(i1:i2, ju1:j2, numDust+numAero+5))
            end if
         end if

         if (do_gcr) then
            allocate(GCR_NOx_nc(i1:i2, ju1:j2, kx1:kx2))
            GCR_NOx_nc = 0.0d0
         endif

         if (lightning_opt == 1) then
            allocate(lightning_nc(i1:i2, ju1:j2, kx1:kx2))
            lightning_nc = 0.0d0

            allocate(flashrate_nc(i1:i2, ju1:j2))
            flashrate_nc = 0.0d0

            allocate(cmf_nc(i1:i2, ju1:j2, kx1:kx2))
            cmf_nc = 0.0d0

            allocate(cmi_flags_nc(i1:i2, ju1:j2))
            cmi_flags_nc = 0.0d0

            allocate(cmi_flags1_nc(i1:i2, ju1:j2))
            cmi_flags1_nc = 0.0d0

            allocate(cldmas0_nc(i1:i2, ju1:j2))
            cldmas0_nc = 0.0d0
         endif

         if (pr_sulf_src) then
            Allocate (dms_no3(i1:i2, ju1:j2, kx1:kx2))
            Allocate (dms_oh (i1:i2, ju1:j2, kx1:kx2))
            dms_no3 = 0.0d0
            dms_oh = 0.0d0

            Allocate (so2_h2o2(i1:i2, ju1:j2, kx1:kx2))
            Allocate (so2_o3  (i1:i2, ju1:j2, kx1:kx2))
            Allocate (so2_oh  (i1:i2, ju1:j2, kx1:kx2))
            so2_h2o2 = 0.0d0
            so2_o3 = 0.0d0
            so2_oh = 0.0d0
         end if

         if (do_mean) then
            Allocate (const_mean(i1:i2, ju1:j2, kx1:kx2, num_const_outrecs))
            const_mean = 0.0d0

            if (pr_psf) then
               Allocate (psf_mean (i1:i2, ju1:j2))
               psf_mean = 0.0d0
            end if
  
            if (pr_tropopausePress) then
               Allocate (tropopausePress_mean (i1:i2, ju1:j2))
               tropopausePress_mean = 0.0d0
            end if
  
            if (pr_potentialVorticity) then
               Allocate (potentialVorticity_mean (i1:i2, ju1:j2, kx1:kx2))
               potentialVorticity_mean = 0.0d0
            end if
  
            if (pr_overheadO3col) then
               Allocate (overheadO3col_mean (i1:i2, ju1:j2, kx1:kx2))
               overheadO3col_mean = 0.0d0
            end if
 
            if (pr_kel) then
               Allocate (kel_mean (i1:i2, ju1:j2, kx1:kx2))
               kel_mean = 0.0d0
            end if
  
            if (pr_mass) then
               Allocate (mass_mean(i1:i2, ju1:j2, kx1:kx2))
               mass_mean = 0.0d0
            end if
  
            if (pr_relHumidity) then
               Allocate (relHumidity_mean(i1:i2, ju1:j2, kx1:kx2))
               relHumidity_mean = 0.0d0
            end if
  
            if (pr_grid_height) then
               Allocate (grid_height_mean(i1:i2, ju1:j2, kx1:kx2))
               grid_height_mean = 0.0d0
            end if
  
            if (pr_metwater) then
               Allocate (metwater_mean(i1:i2, ju1:j2, kx1:kx2))
               metwater_mean = 0.0d0
            end if
         end if

      end if

      return

      end subroutine allocateVariablesConcentration
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Define_Netcdf_Out_Const
!
! !INTERFACE:
!
      subroutine Define_Netcdf_Out_Const (const_var_name)
!
! !USES:
      use m_ncGeneralOpsOutput, only : Define_Netcdf_Out_Gen

      implicit none

#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
      character (len=*), intent(in) :: const_var_name
!
! !DESCRIPTION:
!  Makes the necessary definitions for the const NetCDF output file.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_VAR_NAME) :: prsp1_dim_name
      integer :: ierr
      integer :: nchr1, nchr2
      integer :: nstr
      integer :: omode
      integer :: pos1
      integer :: varid
      integer :: chrd1(1), chrd2 (1)
      integer :: hdfd (1), lond (1), latd  (1)
      integer :: prsd (1), prsp1d(1), recd (1), spcd  (1)
      integer :: spcd_2(1), spcd_3(1), spcd_4(1), spcd_6(1), spcd_7(1), spcd_8(1)
      integer :: strd (1)
      integer :: new_num_emiss

      integer :: var2(2), var3(3), var4(4), var5(5)
      character(len=200) :: AttrValue
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) then
        Write (6,*) 'Define_Netcdf_Out_Const called by ', procID
      end if

      nchr1 = MAX_LENGTH_SPECIES_NAME
      nchr2 = 50

      nstr  =  1

!     ==========================
      call Define_Netcdf_Out_Gen  &
!     ==========================
     &  (ncid_const, nhdf, numLon, numLat, numVert, hdfd, lond, latd, prsd,  &
     &   hdf_dim_name, lon_dim_name, lat_dim_name, prs_dim_name)

!     ------------------
!     Define dimensions.
!     ------------------

      prsp1_dim_name = prs_dim_name
      pos1 = Len_Trim (prs_dim_name) + 1
      prsp1_dim_name(pos1:pos1+2) = 'p1'

      call NcDef_dimension(ncid_const, prsp1_dim_name, numVert+1, prsp1d(1))
!                                     --------------
      call NcDef_dimension(ncid_const, 'chr_dim1', nchr1, chrd1(1))
!                                     ----------
      if (met_opt /= 1) then
         call NcDef_dimension(ncid_const, 'chr_dim2', nchr2, chrd2(1))
!                                        ----------

         call NcDef_dimension(ncid_const, 'str_dim', nstr, strd(1))
!                                       ---------
      end if

      call NcDef_dimension(ncid_const, spc_dim_name, num_const_outrecs, spcd(1))
!                                     ------------

      IF (pr_surf_emiss) THEN
         IF (.NOT. pr_emiss_all) THEN
            call NcDef_dimension  &
     &         (ncid_const,'emiss_spc_dim', num_emiss_outrecs, spcd_2(1))
!                        ---------------
         END IF

         new_num_emiss = 6
         call NcDef_dimension  &
     &       (ncid_const,'emiss_spc_dim2', new_num_emiss, spcd_6(1))
!                      ----------------
         if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
            call NcDef_dimension  &
     &       (ncid_const,'aerosolSurfEmissDim', numDust+numAero+5, spcd_7(1))
         end if
      END IF

      if (pr_emiss_3d) then
         if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
            call NcDef_dimension (ncid_const,'aerosolEmiss3dDim', 5, spcd_8(1))
         end if
      end if

      IF (pr_dry_depos .AND. .NOT. pr_drydep_all) THEN
            call NcDef_dimension  &
     &         (ncid_const,'drydep_spc_dim', num_drydep_outrecs, spcd_3(1))
!                        ----------------
      END IF

      IF (pr_wet_depos .AND. .NOT. pr_wetdep_all) THEN
            call NcDef_dimension  &
     &         (ncid_const,'wetdep_spc_dim', num_wetdep_outrecs, spcd_4(1))
!                        ----------------
      END IF


      call NcDef_dimension(ncid_const, rec_dim_name, NF_UNLIMITED, recd(1))
!                                     ------------

!     -----------------------------------------
!     Define variables and variable attributes.
!     -----------------------------------------

      call NcDef_variable (ncid_const, spc_dim_name, NF_FLOAT, 1, spcd, varid)
!                                     ------------
      call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Species index')
      call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_const, varid, 'coord_labels', 'const_labels')
      call NcDef_var_attributes (ncid_const, varid, 'selection_category', 'NULL')

      call NcDef_variable (ncid_const, 'am', NF_FLOAT, 1, prsd, varid)
!                                     ----
      call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Hybrid pressure term')
      call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')

      call NcDef_variable (ncid_const, 'bm', NF_FLOAT, 1, prsd, varid)
!                                     ----
      call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Hybrid sigma term')
      call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')

      call NcDef_variable (ncid_const, 'ai', NF_FLOAT, 1, prsp1d, varid)
!                                     ----
      call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Hybrid pressure edge term')
      call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')

      call NcDef_variable (ncid_const, 'bi', NF_FLOAT, 1, prsp1d, varid)
!                                     ----
      call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Hybrid sigma edge term')
      call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')

      call NcDef_variable (ncid_const, 'pt', NF_FLOAT, 0, prsp1d, varid)
!                                     ----
      if (met_opt /= 1) then
        var2(:) = (/ chrd2(1), strd(1) /)
        call NcDef_variable (ncid_const, 'metdata_name', NF_CHAR, 2, var2, varid)
!                                       --------------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Metdata name')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
      end if


      var2(:) = (/ chrd1(1), spcd(1) /)
      call NcDef_variable (ncid_const, 'const_labels', NF_CHAR, 2, var2, varid)
!                                     --------------
      call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Species name')
      call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_const, varid, 'selection_category', 'NULL')

      var2(:) = (/ lond(1), latd(1) /)
      call NcDef_variable (ncid_const, 'mcor', NF_FLOAT, 2, var2, varid)
!                                     ------
      call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Grid box area')
      call NcDef_var_attributes (ncid_const, varid, 'units', 'm^2')

      var2(:) = (/ hdfd(1), recd(1) /)
      call NcDef_variable (ncid_const, hdr_var_name, NF_INT, 2, var2, varid)
!                                     ------------
      call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Header')
      call NcDef_var_attributes (ncid_const, varid, 'units', 'gmi_sec, nymd, nhms')

      var5(:) = (/ lond(1), latd(1), prsd(1), spcd(1), recd(1) /)
      call NcDef_variable (ncid_const, const_var_name, NF_FLOAT, 5, var5, varid)
!                                     --------------
      call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Constituent')
      call NcDef_var_attributes (ncid_const, varid, 'units', 'volume mixing ratio')

!    -----------------
      if (pr_psf) then
        var3(:) = (/ lond(1), latd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'psf', NF_FLOAT, 3, var3, varid)
!                                       -----
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Surface pressure')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'mb')
      end if

      if (pr_tropopausePress) then
         var3(:) = (/ lond(1), latd(1), recd(1) /)
         call NcDef_variable (ncid_const, 'tropopausePress', NF_FLOAT, 3, var3, varid)
!                                     -----
         call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Tropopause pressure')
         call NcDef_var_attributes (ncid_const, varid, 'units', 'mb')
      end if

      if (pr_potentialVorticity) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'potentialVorticity', NF_FLOAT, 4, var4, varid)
!                                       -----
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Potential vorticity')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'm2 s-1K kg-1')
      end if

      if (pr_kel) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'kel', NF_FLOAT, 4, var4, varid)
!                                       -----
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Temperature')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'K')
      end if

      if (pr_overheadO3col) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'overheadO3col', NF_FLOAT, 4, var4, varid)
!                                       -----
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Photolysis overhead ozone column')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'DU')
      end if

!... decay diagnostic for tracer mechanism. (jules, hyl, 1/8/11)
      if (pr_decay) then
        var5(:) = (/ lond(1), latd(1), prsd(1), spcd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'decay_3d_out', NF_FLOAT, 5, var5, varid)
!                                       -----
        call NcDef_var_attributes (ncid_const, varid, 'long_name', '3D Decay')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/m^2/constOutputFrequency')
      end if

      if (pr_mass) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'mass', NF_FLOAT, 4, var4, varid)
!                                       ------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Mass')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'kg')
      end if

      if (pr_grid_height) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'gridbox_height', NF_FLOAT, 4, var4, varid)
!                                       ------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Grid Box Height')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'm')
      end if 

      if (pr_relHumidity) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'relHumidity', NF_FLOAT, 4, var4, varid)
!                                       ------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Relative Humidity')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
      end if 

      if (do_gcr) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'GCR_NOx_nc', NF_FLOAT, 4,  &
     &                                                 var4, varid)
        call NcDef_var_attributes (ncid_const, varid, 'long_name',  &
     &                                              'Galactic_Cosmic_Ray_NOx')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/output interval')
!                                    -----
      endif

      if (lightning_opt == 1) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'lightning_nc', NF_FLOAT, 4,  &
     &                                                 var4, varid)
        call NcDef_var_attributes (ncid_const, varid, 'long_name',  &
     &                                              'lightning_no')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/sec')
!                                    -----

        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'cmf_nc', NF_FLOAT, 4,  &
     &                                                 var4, varid)
        call NcDef_var_attributes (ncid_const, varid, 'long_name',  &
     &                                              'Convective mass flux')
        call NcDef_var_attributes (ncid_const, varid, 'units',  &
     &                                               'kg/m^2 sec.')
!                                    -----

        var3(:) = (/ lond(1), latd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'flashrate_nc', NF_FLOAT, 3,  &
     &                                                 var3, varid)
        call NcDef_var_attributes (ncid_const, varid, 'long_name',  &
     &                                                   'flashrate')
        call NcDef_var_attributes (ncid_const, varid, 'units',  &
     &                                    'flashes/sec')
!                                    -----
        var3(:) = (/ lond(1), latd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'cldmas0_nc', NF_FLOAT, 3,  &
     &                                                 var3, varid)
        call NcDef_var_attributes (ncid_const, varid, 'long_name',  &
     &                                                   'cldmas')
        call NcDef_var_attributes (ncid_const, varid, 'units',  &
     &                                    '????')

        var3(:) = (/ lond(1), latd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'cmi_flags1_nc', NF_INT, 3,  &
     &                                                 var3, varid)
        call NcDef_var_attributes (ncid_const, varid, 'long_name',  &
     &                                                   'cmi_flags1')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')

        var3(:) = (/ lond(1), latd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'cmi_flags_nc', NF_INT, 3,  &
     &                                                 var3, varid)
        call NcDef_var_attributes (ncid_const, varid, 'long_name',  &
     &                                                   'Continental/Marine/Ice flags')
        call NcDef_var_attributes (ncid_const, varid, 'units',  &
     &                                    'unitless')
      end if


      if (pr_metwater) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_const, 'metwater', NF_FLOAT, 4, var4, varid)
!                                       ----------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Meteorological Water')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'volume mixing ratio')
      end if

      if (pr_dry_depos) then
         if (.not. pr_drydep_all) then
            call NcDef_variable (ncid_const, 'drydep_spc_dim', NF_FLOAT, 1, spcd_3, varid)
!                                           ---------------
            call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Species index for dry_dep')
            call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
            call NcDef_var_attributes (ncid_const, varid, 'coord_labels', 'drydep_spc_labels')
            call NcDef_var_attributes (ncid_const, varid, 'selection_category', 'NULL')

            var2(:) = (/ chrd1(1), spcd_3(1) /)
            call NcDef_variable (ncid_const, 'drydep_spc_labels', NF_CHAR, 2, var2, varid)
!                                            -----------------
            call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Species name for Dry Deposition')
            call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
            call NcDef_var_attributes (ncid_const, varid, 'selection_category', 'NULL')
         end if
         if (pr_drydep_all) then
            var4(:) = (/ lond(1), latd(1), spcd(1), recd(1) /)
         else
            var4(:) = (/ lond(1), latd(1), spcd_3(1), recd(1) /)
         end if
         call NcDef_variable (ncid_const, 'dry_depos', NF_FLOAT, 4, var4, varid)
!                                        -----------
         call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Dry deposition')
         call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/m^2/constOutputFrequency')
      end if

      if (pr_wet_depos) then
         if (.not. pr_wetdep_all) then
            call NcDef_variable (ncid_const, 'wetdep_spc_dim', NF_FLOAT, 1, spcd_4, varid)
!                                           ---------------
            call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Species index for wet_dep')
            call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
            call NcDef_var_attributes (ncid_const, varid, 'coord_labels', 'wetdep_spc_labels')
            call NcDef_var_attributes (ncid_const, varid, 'selection_category', 'NULL')

            var2(:) = (/ chrd1(1), spcd_4(1) /)
            call NcDef_variable (ncid_const, 'wetdep_spc_labels', NF_CHAR, 2, var2, varid)
!                                            -----------------

            call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Species name for Wet Deposition')
            call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
            call NcDef_var_attributes (ncid_const, varid, 'selection_category', 'NULL')
         end if

         if (pr_wetdep_all) then
            var4(:) = (/ lond(1), latd(1), spcd(1), recd(1) /)
         else
            var4(:) = (/ lond(1), latd(1), spcd_4(1), recd(1) /)
         end if

         call NcDef_variable (ncid_const, 'wet_depos', NF_FLOAT, 4, var4, varid)
!                                        -----------
         call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Wet deposition')
         call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/m^2/constOutputFrequency')
      end if

      if (pr_surf_emiss) then

        IF (.NOT. pr_emiss_all) THEN
           if (num_emiss_outrecs > 0) then

              call NcDef_variable (ncid_const, 'emiss_spc_dim', NF_FLOAT, 1, spcd_2, varid)
!                                             ---------------
              call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Species index for surf_emiss')
              call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
              call NcDef_var_attributes (ncid_const, varid, 'coord_labels', 'emiss_spc_labels')
              call NcDef_var_attributes (ncid_const, varid, 'selection_category', 'NULL')

              var2(:) = (/ chrd1(1), spcd_2(1) /)
              call NcDef_variable (ncid_const, 'emiss_spc_labels', NF_CHAR, 2, var2, varid)
!                                              -----------------
              call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Species name for Surface Emission')
              call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
              call NcDef_var_attributes (ncid_const, varid, 'selection_category', 'NULL')
           end if

        END IF

        if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
           call NcDef_variable (ncid_const, 'aerosolSurfEmissDim', NF_FLOAT, 1, spcd_7, varid)
           call NcDef_var_attributes (ncid_const, varid, 'long_name', &
                                       'Species index for aerosol surface emission')
           call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
           call NcDef_var_attributes (ncid_const, varid, 'coord_labels', 'aerosolSurfEmissLabels')
           call NcDef_var_attributes (ncid_const, varid, 'selection_category', 'NULL')

           var2(:) = (/ chrd1(1), spcd_7(1) /)
           call NcDef_variable (ncid_const, 'aerosolSurfEmissLabels', NF_CHAR, 2, var2, varid)
           call NcDef_var_attributes (ncid_const, varid, 'long_name', &
                                       'Species name for Aerosol Surface Emission')
           call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
           call NcDef_var_attributes (ncid_const, varid, 'selection_category', 'NULL')

           var4(:) = (/ lond(1), latd(1), spcd_7(1), recd(1) /)
           call NcDef_variable (ncid_const, 'aerosolSurfEmiss', NF_FLOAT, 4, var4, varid)
           call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Surface emission for aerosols')
           call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/m^2/constOutputFrequency')
        end if

        IF (pr_emiss_all) THEN
           var4(:) = (/ lond(1), latd(1), spcd(1), recd(1) /)
        ELSE
           var4(:) = (/ lond(1), latd(1), spcd_2(1), recd(1) /)
        ENDIF

        call NcDef_variable (ncid_const, 'surf_emiss', NF_FLOAT, 4, var4, varid)
!                                        -----------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Surface emission')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/m^2/constOutputFrequency')
!!!!!
!!!!!
        call NcDef_variable  &
     &      (ncid_const, 'emiss_spc_dim2', NF_FLOAT, 1, spcd_6, varid)
!                      ---------------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Species index for surf_emiss')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
        call NcDef_var_attributes (ncid_const, varid, 'coord_labels', 'emiss_spc_labels2')
        call NcDef_var_attributes (ncid_const, varid, 'selection_category', 'NULL')

        var2(:) = (/ chrd1(1), spcd_6(1) /)
        call NcDef_variable (ncid_const, 'emiss_spc_labels2', NF_CHAR, 2, var2, varid)
!                                        -----------------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Species name for Surface Emission')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
        call NcDef_var_attributes (ncid_const, varid, 'selection_category', 'NULL')

        var4(:) = (/ lond(1), latd(1), spcd_6(1), recd(1) /)
        call NcDef_variable (ncid_const, 'surf_emiss2', NF_FLOAT, 4, var4, varid)
!                                        -----------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'Surface emission')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/m^2/constOutputFrequency')
      end if

!     ------- 3D Emission ---------------
      if (pr_emiss_3d) then

        IF (pr_emiss_all) THEN
            var5(:) = (/ lond(1), latd(1), prsd(1), spcd(1), recd(1) /)
        ELSE
            var5(:) = (/ lond(1), latd(1), prsd(1), spcd_2(1), recd(1) /)
        ENDIF

        call NcDef_variable (ncid_const, 'emiss_3d_out', NF_FLOAT, 5, var5, varid)
!                                       --------------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', '3D Emission')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/m^2/constOutputFrequency')

        if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
           call NcDef_variable (ncid_const, 'aerosolEmiss3dDim', NF_FLOAT, 1, spcd_8, varid)
           call NcDef_var_attributes (ncid_const, varid, 'long_name', &
                                       'Species index for aerosol 3d emission')
           call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
           call NcDef_var_attributes (ncid_const, varid, 'coord_labels', 'aerosolEmiss3dLabels')
           call NcDef_var_attributes (ncid_const, varid, 'selection_category', 'NULL')

           var2(:) = (/ chrd1(1), spcd_8(1) /)
           call NcDef_variable (ncid_const, 'aerosolEmiss3dLabels', NF_CHAR, 2, var2, varid)
           call NcDef_var_attributes (ncid_const, varid, 'long_name', &
                                       'Species name for Aerosol 3D Emission')
           call NcDef_var_attributes (ncid_const, varid, 'units', 'unitless')
           call NcDef_var_attributes (ncid_const, varid, 'selection_category', 'NULL')

           var5(:) = (/ lond(1), latd(1), prsd(1), spcd_8(1), recd(1) /)
           call NcDef_variable (ncid_const, 'aerosolEmiss3d', NF_FLOAT, 5, var5, varid)
           call NcDef_var_attributes (ncid_const, varid, 'long_name', &
                                       '3D emission for aerosols')
           call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/m^2/constOutputFrequency')
        end if
      endif
!     -------------------------------------

      if (pr_sulf_src) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)

        call NcDef_variable (ncid_const, 'dms_oh', NF_FLOAT, 4, var4, varid)
!                                       --------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'SO2 from DMS+OH')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/constOutputFrequency')

        call NcDef_variable (ncid_const, 'dms_no3', NF_FLOAT, 4, var4, varid)
!                                       ---------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'SO2 from DMS+NO3')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/constOutputFrequency')

        call NcDef_variable (ncid_const, 'so2_oh', NF_FLOAT, 4, var4, varid)
!                                       --------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'SO4 from SO2+OH')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/constOutputFrequency')

        call NcDef_variable (ncid_const, 'so2_h2o2', NF_FLOAT, 4, var4, varid)
!                                       ----------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'SO4 from SO2+H2O2')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/constOutputFrequency')

        call NcDef_variable (ncid_const, 'so2_o3', NF_FLOAT, 4, var4, varid)
!                                       --------
        call NcDef_var_attributes (ncid_const, varid, 'long_name', 'SO4 from SO2+O3')
        call NcDef_var_attributes (ncid_const, varid, 'units', 'kg/constOutputFrequency')
      end if

!     -------------------------
!     Define global attributes.
!     -------------------------

      call NcDef_glob_attributes(ncid_const,'title','Gmimod constituent file')
!
      AttrValue = 'DESCRIPTION OF WATER RELATED SPECIES AND VARIABLES'
      call NcDef_glob_attributes(ncid_const,'Note', Trim(AttrValue) )

      AttrValue = 'read in from meteorological fields used to drive &
                   transport - only used below tropopause.'

      call NcDef_glob_attributes(ncid_const,'metwater', Trim(AttrValue) )

      AttrValue = 'resets metwater to > 3 ppmv when less than 3 ppmv - &
                   not transported'

      call NcDef_glob_attributes(ncid_const,'H2O', Trim(AttrValue) )

      AttrValue = 'used in aircraft emissions experiments only - &
                   not transported simulations'

      call NcDef_glob_attributes(ncid_const,'H2OAIR', Trim(AttrValue) )

      AttrValue = 'accounts for removal/addition of water by PSC &
                   formation/dehydration and gravitational settling &
                   in stratosphere - modifies h2oclim in stratosphere &
                   - should always be transported'

      call NcDef_glob_attributes(ncid_const,'DEHYD', Trim(AttrValue) )

      AttrValue = 'stratospheric water climatology (UARS MLS) read &
                   in from file - not outputted to *const.nc'

      call NcDef_glob_attributes(ncid_const,'h2oclim', Trim(AttrValue) )
!
      call NcSetFill (ncid_const, NF_NOFILL, omode)

      call NcEnd_def (ncid_const)

      return

      end subroutine Define_Netcdf_Out_Const
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Hdr_Const
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Hdr_Const(gmiDomain, Emission, metFields, &
     &                 const_labels)
!
! !USES:
      use m_ncGeneralOpsOutput, only : Write_Netcdf_Hdr_Gen
!
      implicit none
!
! !INPUT PARAMETERS:
      character (len=MAX_LENGTH_SPECIES_NAME), intent(in) :: const_labels(1:)
      type(t_Emission ), intent(in) :: Emission 
      type(t_gmiDomain), intent(in) :: gmiDomain
      type(t_metFields), intent(in) :: metFields
!
! !DESCRIPTION:
!  This routine creates some header information for the const NetCDF output
!  file and writes it out.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_SPECIES_NAME) :: spclab(num_const_outrecs)
      character (len=MAX_LENGTH_SPECIES_NAME) :: spclab_2(num_emiss_outrecs)
      character (len=MAX_LENGTH_SPECIES_NAME) :: spclab_3(num_drydep_outrecs)
      character (len=MAX_LENGTH_SPECIES_NAME) :: spclab_4(num_wetdep_outrecs)
      character (len=50) :: metnam(1)

      integer :: ic, icx
      integer, allocatable :: aerosolSurfEmissMap(:)
      integer, allocatable :: emiss_map_dust(:)
      integer, allocatable :: emiss_map_aero(:)

      integer :: cnt1d (1), cnt2d (2)
      integer :: strt1d(1), strt2d(2)

      real*8  :: prsdat(k1:k2)

      real*8  :: spcdat(num_const_outrecs)
      real*8  :: spcdat_2(num_emiss_outrecs)
      real*8  :: spcdat_3(num_drydep_outrecs)
      real*8  :: spcdat_4(num_wetdep_outrecs)

      character (len=MAX_LENGTH_SPECIES_NAME) :: spclab_6(6), aerosol_spclab(numDust+numAero+5)
      real*8             :: spcdat_6(6), aerosol_spcdat(numDust+numAero+5)
      character (len=MAX_LENGTH_SPECIES_NAME) :: aerosol3d_spclab(5)
      real*8             :: aerosol3d_spcdat(5)
      integer            :: new_num_emiss
      real*8 , allocatable :: mcor(:,:)
      real*8 , allocatable :: londeg(:), latdeg(:)
      real*8 , allocatable :: ai(:), bi(:)
      real*8 , allocatable :: am(:), bm(:)
      real*8 :: pt
      character (len=50) :: metdata_name
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Hdr_Const called by ', procID

      allocate(latdeg(ju1_gl:j2_gl))
      call Get_latdeg(gmiDomain, latdeg)

      allocate(londeg(i1_gl:i2_gl))
      call Get_londeg(gmiDomain, londeg)

      allocate(mcor(i1_gl:i2_gl,ju1_gl:j2_gl))
      call Get_mcorGlob(gmiDomain, mcor)

!     =========================
      call Write_Netcdf_Hdr_Gen  &
!     =========================
     &  (ncid_const, latdeg, londeg, pr_diag, procID, i1_gl, i2_gl, ju1_gl, j2_gl,&
     &   hdf_dim_name, lat_dim_name, lon_dim_name)

!     ---------
!     Pressure.
!     ---------

      allocate(ai(k1-1:k2))
      allocate(bi(k1-1:k2))
      allocate(am(k1:k2))
      allocate(bm(k1:k2))

      call Get_metdata_name(metFields, metdata_name)
      call Get_pt(metFields, pt)
      call Get_ai(metFields, ai)
      call Get_bi(metFields, bi)
      call Get_am(metFields, am)
      call Get_bm(metFields, bm)

      prsdat(:) = (am(:) * pt)  + (bm(:) * 1000.0d0)

      strt1d(1) = 1
      cnt1d (1) = numVert

      call Ncwr_1d (prsdat(kx1:kx2), ncid_const, prs_dim_name, strt1d, cnt1d)
      call Ncwr_1d (am(kx1:kx2), ncid_const, 'am', strt1d, cnt1d)
      call Ncwr_1d (bm(kx1:kx2), ncid_const, 'bm', strt1d, cnt1d)

      cnt1d(1) = numvert + 1

      call Ncwr_1d (ai(kx1-1:kx2), ncid_const, 'ai', strt1d, cnt1d)
      call Ncwr_1d (bi(kx1-1:kx2), ncid_const, 'bi', strt1d, cnt1d)

      call Ncwr_Scal (pt, ncid_const, 'pt')

!     -------------
!     Metdata name.
!     -------------

      if (met_opt /= 1) then
        metnam(1) = metdata_name

        strt2d(:) = (/  1, 1 /)
        cnt2d (:) = (/ 50, 1 /)

        call Ncwr_2d_Char (metnam, ncid_const, 'metdata_name', strt2d, cnt2d)
      end if

!     ------------
!     Species map.
!     ------------

      if (pr_diag) then
        Write (6,*) 'Const NetCDF output info:'
        Write (6,*) '  num_const_outrecs:         ', num_const_outrecs

        do ic = 1, num_const_outrecs
          Write (6,900) ic, const_outrec_map(ic), trim(const_labels(const_outrec_map(ic)))
        end do

 900    format ('   const_outrec to species mapping:  ', 2i4,a16)
      end if

      do ic = 1, num_const_outrecs
        spcdat(ic) = const_outrec_map(ic)
      end do

      strt1d(1) = 1
      cnt1d (1) = num_const_outrecs

      call Ncwr_1d (spcdat, ncid_const, spc_dim_name, strt1d, cnt1d)

!     ---------------
!     Species labels.
!     ---------------

      do ic = 1, num_const_outrecs
        spclab(ic) = const_labels(const_outrec_map(ic))
      end do

      strt2d(:) = (/  1, 1 /)
      cnt2d (:) = (/ MAX_LENGTH_SPECIES_NAME, num_const_outrecs /)

      call Ncwr_2d_Char (spclab, ncid_const, 'const_labels', strt2d, cnt2d)

!     --------------
!     Grid box area.
!     --------------

      strt2d(:) = (/ 1, 1 /)
      cnt2d (:) = (/ numLon, numLat /)

      call Ncwr_2d (mcor, ncid_const, 'mcor', strt2d, cnt2d)

!     ------------------------------------------------------------------
!     Surface Emission, dry deposition and wet deposition Species labels.
!     ------------------------------------------------------------------
!     This portion of the code is only used if the user has
!     selected in the namelist file the species for which
!     he/she wants surface emission, dry deposition, and
!     wet deposition diagnostic information to be saved.
!
      if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
         allocate(aerosolSurfEmissMap(numDust+numAero+5))
         allocate(emiss_map_dust(numDust))
         allocate(emiss_map_aero(numAero))

         call Get_emiss_map_aero(Emission, emiss_map_aero, numAero)
         call Get_emiss_map_dust(Emission, emiss_map_dust, numDust)

         aerosolSurfEmissMap(1)                     = IFSO2
         aerosolSurfEmissMap(2)                     = INSO2
         aerosolSurfEmissMap(3)                     = INDMS
         aerosolSurfEmissMap(4)                     = INSO4A
         aerosolSurfEmissMap(5)                     = IFSO4A
         aerosolSurfEmissMap(6:numAero+5)             = emiss_map_aero(1:numAero)
         aerosolSurfEmissMap(numAero+6:numAero+numDust+5) = emiss_map_dust(1:numDust)

      end if

      IF (pr_surf_emiss) THEN
         IF (.NOT. pr_emiss_all) THEN

            do ic = 1, num_emiss_outrecs
              spcdat_2(ic) = emiss_outrec_map(ic)
            end do

            strt1d(1) = 1
            cnt1d (1) = num_emiss_outrecs

            call Ncwr_1d(spcdat_2,ncid_const,'emiss_spc_dim',strt1d,cnt1d)

            do ic = 1, num_emiss_outrecs
              spclab_2(ic) = const_labels(emiss_outrec_map(ic))
            end do

            strt2d(:) = (/  1, 1 /)
            cnt2d (:) = (/ MAX_LENGTH_SPECIES_NAME, num_emiss_outrecs /)

            call Ncwr_2d_Char  &
     &          (spclab_2, ncid_const, 'emiss_spc_labels', strt2d, cnt2d)
         END IF

         new_num_emiss = 6

         do ic = 1, new_num_emiss
            spcdat_6(ic) = ic
         end do

         strt1d(1) = 1
         cnt1d (1) = new_num_emiss

         call Ncwr_1d(spcdat_6,ncid_const,'emiss_spc_dim2',strt1d,cnt1d)

         spclab_6(1) = 'CO_methanol'
         spclab_6(2) = 'CO_monoterpene'
         spclab_6(3) = 'biogenic_propene'
         spclab_6(4) = 'soil_NOx'
         spclab_6(5) = 'HNO3'
         spclab_6(6) = 'O3'

         strt2d(:) = (/  1, 1 /)
         cnt2d (:) = (/ MAX_LENGTH_SPECIES_NAME, new_num_emiss /)

         call Ncwr_2d_Char  &
     &       (spclab_6, ncid_const, 'emiss_spc_labels2', strt2d, cnt2d)

         if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
            do ic = 1, numDust+numAero+5
               aerosol_spcdat(ic) = ic
            end do

            strt1d(1) = 1
            cnt1d (1) = numDust+numAero+5

            call Ncwr_1d(aerosol_spcdat,ncid_const,'aerosolSurfEmissDim',strt1d,cnt1d)

            do ic = 1, numDust+numAero+5
               aerosol_spclab(ic) = const_labels(aerosolSurfEmissMap(ic))
            end do

            strt2d(:) = (/  1, 1 /)
            cnt2d (:) = (/ MAX_LENGTH_SPECIES_NAME, numDust+numAero+5 /)

            call Ncwr_2d_Char  &
     &         (aerosol_spclab, ncid_const, 'aerosolSurfEmissLabels', strt2d, cnt2d)

         end if

      END IF

      if (pr_emiss_3d) then
         if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
            do ic = 1,5
               aerosol3d_spcdat(ic) = ic
            end do

            strt1d(1) = 1
            cnt1d (1) = 5

            call Ncwr_1d(aerosol3d_spcdat,ncid_const,'aerosolEmiss3dDim',strt1d,cnt1d)

            do ic = 1, 5
               aerosol3d_spclab(ic) = const_labels(aerosolSurfEmissMap(ic))
            end do

            strt2d(:) = (/  1, 1 /)
            cnt2d (:) = (/ MAX_LENGTH_SPECIES_NAME, 5 /)

            call Ncwr_2d_Char  &
     &         (aerosol3d_spclab, ncid_const, 'aerosolEmiss3dLabels', strt2d, cnt2d)

         end if

      end if

      IF (pr_dry_depos .AND. (.NOT. pr_drydep_all)) THEN

         do ic = 1, num_drydep_outrecs
           spcdat_3(ic) = drydep_outrec_map(ic)
         end do

         strt1d(1) = 1
         cnt1d (1) = num_drydep_outrecs

         call Ncwr_1d(spcdat_3,ncid_const,'drydep_spc_dim',strt1d,cnt1d)

         do ic = 1, num_drydep_outrecs
           spclab_3(ic) = const_labels(drydep_outrec_map(ic))
         end do

         strt2d(:) = (/  1, 1 /)
         cnt2d (:) = (/ MAX_LENGTH_SPECIES_NAME, num_drydep_outrecs /)

         call Ncwr_2d_Char(spclab_3, ncid_const, 'drydep_spc_labels', strt2d, cnt2d)
      END IF

      IF (pr_wet_depos .AND. (.NOT. pr_wetdep_all)) THEN

         do ic = 1, num_wetdep_outrecs
           spcdat_4(ic) = wetdep_outrec_map(ic)
         end do

         strt1d(1) = 1
         cnt1d (1) = num_wetdep_outrecs

         call Ncwr_1d(spcdat_4,ncid_const,'wetdep_spc_dim',strt1d,cnt1d)

         do ic = 1, num_wetdep_outrecs
           spclab_4(ic) = const_labels(wetdep_outrec_map(ic))
         end do

         strt2d(:) = (/  1, 1 /)
         cnt2d (:) = (/ MAX_LENGTH_SPECIES_NAME, num_wetdep_outrecs /)

         call Ncwr_2d_Char(spclab_4, ncid_const, 'wetdep_spc_labels', strt2d, cnt2d)

      END IF

      return

      end subroutine Write_Netcdf_Hdr_Const
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Prep_Netcdf_Const
!
! !INTERFACE:
!
      subroutine Prep_Netcdf_Const  &
     &  (time_for_nc, Chemistry, Deposition, Emission, SpeciesConcentration, &
     &   metFields)
!
! !USES:
      use GmiDepositionMethod_mod, only : Get_wet_depos, Get_dry_depos
      use GmiDepositionMethod_mod, only : Set_wet_depos, Set_dry_depos
      use GmiChemistryMethod_mod , only : Get_overheadO3col
      use GmiEmissionMethod_mod  , only : Get_aerosolSurfEmiss
      use GmiEmissionMethod_mod  , only : Get_aerosolEmiss3D, Get_emiss_3d_out
      use GmiEmissionMethod_mod  , only : Get_surf_emiss_out, Get_surf_emiss_out2
      use GmiEmissionMethod_mod  , only : Get_lightning_no, Get_flashrate
      use GmiEmissionMethod_mod  , only : Get_cldmas0, Get_cmi_flags1, Get_GCR_NOx

      implicit none
!
! !INPUT PARAMETERS:
      logical                     , intent(in) :: time_for_nc
      type(t_metFields)           , intent(in) :: metFields
      type(t_Chemistry)           , intent(in) :: Chemistry
      type(t_Emission )           , intent(in) :: Emission 
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_Deposition), intent(inOut) :: Deposition
!
! !DESCRIPTION:
! Prepares the species concentration NetCDF output.
!
! !LOCAL VARIABLES:
      integer :: ic, icx, il
      integer, save :: counter = 0
      real*8  :: rcounter
      integer, allocatable :: cmi_flags1 (:,:), cmi_flags(:,:)
      real*8, allocatable :: cldmas0 (:,:)
      real*8, allocatable :: GCR_NOx(:,:,:)
      real*8, allocatable :: lightning_no(:,:,:)
      real*8, allocatable :: flashrate   (:,:)
      real*8, allocatable :: wet_depos(:,:,:)
      real*8, allocatable :: dry_depos(:,:,:)
      real*8, allocatable :: overheadO3col(:,:,:)
      real*8, allocatable :: mass(:,:,:), kel(:,:,:), cmf(:,:,:)
      real*8, allocatable :: gridBoxHeight(:,:,:)
      real*8, allocatable :: tropopausePress(:,:), pctm2(:,:)
      real*8, allocatable :: relativeHumidity(:,:,:), humidity(:,:,:)
      real*8, allocatable :: potentialVorticity(:,:,:)
      real*8, allocatable :: aerosolSurfEmiss(:,:,:)
      real*8, allocatable :: surf_emiss_out  (:,:,:)
      real*8, allocatable :: surf_emiss_out2 (:,:,:)
      real*8, allocatable :: emiss_3d_out  (:,:,:,:)
      type (t_GmiArrayBundle), pointer :: concentration(:)
      
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Prep_Netcdf_Const called by ', procID

      call Get_concentration(SpeciesConcentration, concentration)

      if (pr_overheadO3col) then
         allocate(overheadO3col(i1:i2, ju1:j2, k1:k2))
         call Get_overheadO3col(Chemistry, overheadO3col)
      end if

      if (pr_potentialVorticity) then
         allocate(potentialVorticity(i1:i2, ju1:j2, k1:k2))
         call Get_potentialVorticity(metFields, potentialVorticity)
      end if

      if (pr_tropopausePress) then
         allocate(tropopausePress(i1:i2, ju1:j2))
         call Get_tropopausePress(metFields, tropopausePress)
      end if

      if (pr_psf) then
         allocate(pctm2(ilo:ihi, julo:jhi))
         call Get_pctm2(metFields, pctm2)
      end if

      if (pr_kel) then
         allocate(kel(ilo:ihi, julo:jhi, k1:k2))
         call Get_kel(metFields, kel)
      end if

      if (pr_mass) then
         allocate(mass(i1:i2, ju1:j2, k1:k2))
         call Get_mass(metFields, mass)
      end if

      if (pr_relHumidity) then
         allocate(relativeHumidity(i1:i2, ju1:j2, k1:k2))
         call Get_relativeHumidity(metFields, relativeHumidity)
      end if

      if (pr_grid_height) then
         allocate(gridBoxHeight(i1:i2, ju1:j2, k1:k2))
         call Get_gridBoxHeight(metFields, gridBoxHeight)
      end if

      if (pr_metwater) then
         allocate(humidity(i1:i2, ju1:j2, k1:k2))
         call Get_humidity(metFields, humidity)
      end if

      if (lightning_opt == 1) then
         allocate(cmi_flags(i1:i2, ju1:j2))
         call Get_cmi_flags(metFields, cmi_flags)

         allocate(cmf(i1:i2, ju1:j2, k1:k2))
         call Get_cmf(metFields, cmf)
      end if 

!     ============
      if (do_mean) then
!     ============

        if (counter == 0) then

!         ----------------
!         Reset mean sums.
!         ----------------

          const_mean     (:,:,:,:) = 0.0d0

          if (pr_tropopausePress)    tropopausePress_mean     (:,:) = 0.0d0
          if (pr_overheadO3col)      overheadO3col_mean     (:,:,:) = 0.0d0
          if (pr_potentialVorticity) potentialVorticity_mean(:,:,:) = 0.0d0
          if (pr_psf)      psf_mean     (:,:)   = 0.0d0
          if (pr_kel)      kel_mean     (:,:,:) = 0.0d0
          if (pr_mass)     mass_mean    (:,:,:) = 0.0d0
          if (pr_relHumidity)     relHumidity_mean    (:,:,:) = 0.0d0
          if (pr_grid_height)     grid_height_mean    (:,:,:) = 0.0d0
          if (pr_metwater)        metwater_mean(:,:,:) = 0.0d0

          if(do_gcr)       GCR_NOx_nc(:,:,:) = 0.0d0

          if(lightning_opt == 1) then
             lightning_nc(:,:,:) = 0.0d0
             flashrate_nc(:,:  ) = 0.0d0
                   cmf_nc(:,:,:) = 0.0d0
             cmi_flags_nc(:,:  ) = 0.0d0
            cmi_flags1_nc(:,:  ) = 0.0d0
               cldmas0_nc(:,:  ) = 0.0d0
          endif

        end if

!       -----------------
!       Update mean sums.
!       -----------------

        do icx = 1, num_const_outrecs
           ic = const_outrec_map(icx)
           const_mean(:,:,:,icx) =  const_mean(:,:,:,icx) + &
     &                              concentration(ic)%pArray3D(:,:,kx1:kx2)
        end do

        if(do_gcr) then
             allocate(GCR_NOx(i1:i2,ju1:j2,k1:k2))
             call Get_GCR_NOx(Emission, GCR_NOx)

             GCR_NOx_nc(:,:,kx1:kx2) = GCR_NOx_nc(:,:,kx1:kx2) + GCR_NOx(:,:,k1:k2)

             deallocate(GCR_NOx)
        endif

        if(lightning_opt == 1) then
             allocate(lightning_no(i1:i2,ju1:j2,k1:k2))
             call Get_lightning_no(Emission, lightning_no)

             lightning_nc(:,:,kx1:kx2) = lightning_nc(:,:,kx1:kx2) +  &
     &                                   lightning_no(:,:,k1:k2)

             deallocate(lightning_no)

             allocate(flashrate(i1:i2,ju1:j2))
             call Get_flashrate(Emission, flashrate)

             flashrate_nc(i1:i2,ju1:j2) =  &
     &           flashrate_nc(i1:i2,ju1:j2) + flashrate(i1:i2,ju1:j2)

             deallocate(flashrate)

             cmi_flags_nc(i1:i2,ju1:j2) = cmi_flags(i1:i2,ju1:j2)

             cmf_nc(:,:,kx1:kx2) =  cmf_nc(:,:,kx1:kx2) + cmf(:,:,kx1:kx2)

             allocate(   cldmas0(i1:i2,ju1:j2))
             cldmas0 = 0.d0
             call Get_cldmas0(Emission, cldmas0)
             cldmas0_nc(i1:i2,ju1:j2) = cldmas0(i1:i2,ju1:j2)
             deallocate (cldmas0)

             allocate(cmi_flags1(i1:i2,ju1:j2))
             cmi_flags1 = 0
             call Get_cmi_flags1(Emission, cmi_flags1)
             cmi_flags1_nc(i1:i2,ju1:j2) =  real (cmi_flags1(i1:i2,ju1:j2))
             deallocate(cmi_flags1)

        endif

!       -----------------
!       Update sums.
!       -----------------

        if (pr_tropopausePress) then
           tropopausePress_mean(:,:) = tropopausePress_mean(:,:) + &
     &                                 tropopausePress(:,:)
        end if

        if (pr_potentialVorticity) then
           potentialVorticity_mean (i1:i2,ju1:j2,:) = &
     &       potentialVorticity_mean (i1:i2,ju1:j2,:) +&
     &       potentialVorticity (i1:i2,ju1:j2,kx1:kx2)
        end if

        if (pr_overheadO3col) then
           overheadO3col_mean(:,:,:) = overheadO3col_mean(:,:,:) + overheadO3col(:,:,kx1:kx2)
        end if

        if (pr_psf) &
           psf_mean (i1:i2,ju1:j2) =  psf_mean (i1:i2,ju1:j2)   + pctm2(i1:i2,ju1:j2)

        if (pr_kel) &
           kel_mean (i1:i2,ju1:j2,:) = kel_mean (i1:i2,ju1:j2,:) + kel (i1:i2,ju1:j2,kx1:kx2)

        if (pr_mass) then
           mass_mean(i1:i2,ju1:j2,:) =  &
     &      mass_mean(i1:i2,ju1:j2,:) + mass (i1:i2,ju1:j2,kx1:kx2)
        end if

        if (pr_relHumidity) then
           relHumidity_mean(:,:,:) = relHumidity_mean(:,:,:) + &
     &              relativeHumidity (:,:,kx1:kx2)
        end if

        if (pr_grid_height) then
          grid_height_mean(:,:,:) = grid_height_mean(:,:,:) + &
     &                              gridBoxHeight (:,:,kx1:kx2)
        end if

        if (pr_metwater) then
           metwater_mean(i1:i2,ju1:j2,:) = metwater_mean(i1:i2,ju1:j2,:) +  &
     &        humidity (i1:i2,ju1:j2,kx1:kx2) * MWTAIR / (MWTH2O * GPKG)

        end if

        counter = counter + 1

        if (time_for_nc) then

!         ----------------
!         Calculate means.
!         ----------------

          rcounter = counter

          do ic = 1, num_const_outrecs
            const_mean(:,:,:,ic) = const_mean(:,:,:,ic) / rcounter
          end do

          if (pr_tropopausePress) &
             tropopausePress_mean(:,:) = tropopausePress_mean(:,:) / rcounter
          if (pr_overheadO3col) then
             overheadO3col_mean (:,:,:) = overheadO3col_mean (:,:,:) / rcounter
             ! converting from mol/cm^2 to Dobson unit (DU)
             overheadO3col_mean (:,:,:) = overheadO3col_mean (:,:,:) / 2.690d+16
          end if
          if (pr_potentialVorticity) &
     &       potentialVorticity_mean (:,:,:) = potentialVorticity_mean (:,:,:) / rcounter
          if (pr_psf)  psf_mean (:,:)   = psf_mean (:,:)   / rcounter
          if (pr_kel)  kel_mean (:,:,:) = kel_mean (:,:,:) / rcounter
          if (pr_mass) mass_mean(:,:,:) = mass_mean(:,:,:) / rcounter
          if (pr_relHumidity) relHumidity_mean(:,:,:) = relHumidity_mean(:,:,:) / rcounter
          if (pr_grid_height) grid_height_mean(:,:,:) = grid_height_mean(:,:,:) / rcounter

          if(do_gcr) then
             GCR_NOx_nc(:,:,:) = GCR_NOx_nc(:,:,:) / rcounter
          endif

          if(lightning_opt == 1) then
             flashrate_nc  (:,:) = flashrate_nc  (:,:) / rcounter
             cmf_nc      (:,:,:) = cmf_nc      (:,:,:) / rcounter
             cldmas0_nc    (:,:) = cldmas0_nc    (:,:) / rcounter
             cmi_flags1_nc (:,:) = cmi_flags1_nc (:,:) / rcounter
             lightning_nc(:,:,:) = lightning_nc(:,:,:) / rcounter
          endif

          if (pr_metwater) then
             metwater_mean(:,:,:) = metwater_mean(:,:,:) / rcounter
          end if

          counter = 0

!         -----------------------------------------------
!         Fill output arrays with calculated mean values.
!         -----------------------------------------------

          if (pr_tropopausePress) tropopausePress_nc(:,:) = tropopausePress_mean(:,:)
          if (pr_overheadO3col)  overheadO3col_nc (:,:,:) = overheadO3col_mean (:,:,:)
          if (pr_potentialVorticity) potentialVorticity_nc(:,:,:) = &
     &           potentialVorticity_mean (:,:,:)
          if (pr_psf)  psf_nc (:,:)   = psf_mean (:,:)
          if (pr_kel)  kel_nc (:,:,:) = kel_mean (:,:,:)
          if (pr_mass) mass_nc(:,:,:) = mass_mean(:,:,:)
          if (pr_relHumidity) relHumidity_nc(:,:,:) = relHumidity_mean(:,:,:)
          if (pr_grid_height) gridbox_height_nc(:,:,:) = grid_height_mean(:,:,:)

          if (pr_metwater) metwater_nc(:,:,:) = metwater_mean(:,:,:)

        end if

!     ====
      else
!     ====

        if (time_for_nc) then

!         --------------------------------------------------
!         Fill output arrays with current "non-mean" values.
!         --------------------------------------------------

          if (pr_tropopausePress) tropopausePress_nc(:,:) = tropopausePress(:,:)

          if (pr_psf) psf_nc(i1:i2,ju1:j2) = pctm2(i1:i2,ju1:j2)

          if (pr_potentialVorticity) potentialVorticity_nc(:,:,:) =&
     &          potentialVorticity(i1:i2,ju1:j2,kx1:kx2)
          if (pr_overheadO3col) then
             overheadO3col_nc(:,:,:) = overheadO3col(:,:,kx1:kx2)
          end if
          if (pr_kel) kel_nc(:,:,:) = kel(i1:i2,ju1:j2,kx1:kx2)

          if (pr_mass) mass_nc(:,:,:) = mass(:,:,kx1:kx2)
          if (pr_relHumidity) relHumidity_nc(:,:,:) = relativeHumidity(:,:,kx1:kx2)
          if (pr_grid_height) gridbox_height_nc(:,:,:) = gridBoxHeight(:,:,kx1:kx2)

          if (pr_metwater) then
            metwater_nc(:,:,:) =  &
     &        humidity(:,:,kx1:kx2) * MWTAIR / (MWTH2O * GPKG)
          end if

        end if

      end if

      if (time_for_nc) then

        if (pr_dry_depos) then

           allocate(dry_depos(i1:i2, ju1:j2, numSpecies))
           call Get_dry_depos(Deposition, dry_depos)

           if (pr_drydep_all) then
              ddep_nc(:,:,:) = dry_depos(:,:,:)
           else
              do ic = 1, num_drydep_outrecs
                 ddep_nc(:,:,ic)  &
     &             = dry_depos(:,:,drydep_outrec_map(ic))
              end do
           end if

           dry_depos(:,:,:) = 0.0d0
           call Set_dry_depos(Deposition, dry_depos)
        end if

        if (pr_wet_depos) then

           allocate(wet_depos(i1:i2, ju1:j2, numSpecies))
           call Get_wet_depos(Deposition, wet_depos)

           if (pr_wetdep_all) then
              wdep_nc(:,:,:) = wet_depos(:,:,:)
           else
              do ic = 1, num_wetdep_outrecs
                 wdep_nc(:,:,ic)  &
     &             = wet_depos(:,:,wetdep_outrec_map(ic))
              end do
           end if

           wet_depos(:,:,:) = 0.0d0
           call Set_wet_depos(Deposition, wet_depos)
        end if

        if (pr_surf_emiss) then
           allocate(surf_emiss_out(i1:i2, ju1:j2, numSpecies))
           call Get_surf_emiss_out(Emission, surf_emiss_out)

           if (pr_emiss_all) then
              semiss_out_nc(:,:,:) = surf_emiss_out(:,:,:)
           else
              do ic = 1, num_emiss_outrecs
                 semiss_out_nc(:,:,ic)  &
     &             = surf_emiss_out(:,:,emiss_outrec_map(ic))
              end do
           end if

           allocate(surf_emiss_out2(i1:i2, ju1:j2, 6))
           call Get_surf_emiss_out2(Emission, surf_emiss_out2)
           semiss_out2_nc(:,:,:) = surf_emiss_out2(:,:,:)
          
           if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
              allocate(aerosolSurfEmiss(i1:i2, ju1:j2, numDust+numAero+5))
              call Get_aerosolSurfEmiss(Emission, aerosolSurfEmiss)
              aerosolSurfEmiss_nc(:,:,:) = aerosolSurfEmiss(:,:,:)
           end if

        end if
      end if

      if (pr_kel) deallocate(kel)
      if (pr_psf) deallocate(pctm2)
      if (pr_mass) deallocate(mass)
      if (pr_metwater) deallocate(humidity)
      if (pr_relHumidity) deallocate(relativeHumidity)
      if (pr_grid_height) deallocate(gridBoxHeight)
      if (pr_overheadO3col) deallocate(overheadO3col)
      if (pr_tropopausePress) deallocate(tropopausePress)
      if (pr_potentialVorticity) deallocate(potentialVorticity)

      return

      end subroutine Prep_Netcdf_Const
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Const
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Const (Chemistry, Emission, nymd, nhms, gmi_sec)
!
! USES:
      use GmiChemistryMethod_mod, only : Get_dms_oh  , Set_dms_oh
      use GmiChemistryMethod_mod, only : Get_dms_no3 , Set_dms_no3
      use GmiChemistryMethod_mod, only : Get_so2_h2o2, Set_so2_h2o2
      use GmiChemistryMethod_mod, only : Get_so2_oh  , Set_so2_oh
      use GmiChemistryMethod_mod, only : Get_so2_o3  , Set_so2_o3
      use GmiEmissionMethod_mod , only : Set_surf_emiss_out, Set_surf_emiss_out2
      use GmiEmissionMethod_mod , only : Set_aerosolSurfEmiss

!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: nymd, nhms
      real*8 , intent(in) :: gmi_sec
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Emission ), intent(inOut) :: Emission 
      type (t_Chemistry), intent(inOut) :: Chemistry
!
! !DESCRIPTION:
! This routine writes the const NetCDF output.
!
! !DEFINED PARAMETERS:
!      integer, parameter :: SG_CONST     = 4000
!      integer, parameter :: SG_EMISS3D   = 4001
      integer, parameter :: SG_PSF       = 4002
      integer, parameter :: SG_TROPP     = 4003
      integer, parameter :: SG_OvO3col   = 4004
      integer, parameter :: SG_PV        = 4005
      integer, parameter :: SG_KEL       = 4006
      integer, parameter :: SG_MASS      = 4007
      integer, parameter :: SG_GBH       = 4008
      integer, parameter :: SG_RH        = 4009
      integer, parameter :: SG_DRY       = 4010
      integer, parameter :: SG_WET       = 4011
      integer, parameter :: SG_METWATER  = 4012
      integer, parameter :: SG_EMISS     = 4013
      integer, parameter :: SG_EMISS2    = 4014
      integer, parameter :: SG_DMS_OH    = 4015
      integer, parameter :: SG_DMS_NO3   = 4016
      integer, parameter :: SG_SO2_H2O2  = 4017
      integer, parameter :: SG_SO2_O3    = 4018
      integer, parameter :: SG_LIGHTNING = 4019
      integer, parameter :: SG_FLASHRATE = 4020
      integer, parameter :: SG_CLDMASS   = 4021
      integer, parameter :: SG_CMF       = 4022
      integer, parameter :: SG_CMI1      = 4023
      integer, parameter :: SG_CMI2      = 4024
      integer, parameter :: SG_SO2_OH    = 4025
      integer, parameter :: SG_GCR       = 4026
!
! !LOCAL VARIABLES:
      integer :: ic
      integer :: cnt2d (2), cnt3d (3), cnt4d (4)
      integer :: strt2d(2), strt3d(3), strt4d(4)
      integer :: hdr(NETCDF_HDF)
      INTEGER :: Numb_Spec
      real*8 , allocatable :: zero3d(:,:,:)
      real*8 , allocatable :: loc_dms_oh  (:,:,:)
      real*8 , allocatable :: loc_dms_no3 (:,:,:)
      real*8 , allocatable :: loc_so2_oh  (:,:,:)
      real*8 , allocatable :: loc_so2_o3  (:,:,:)
      real*8 , allocatable :: loc_so2_h2o2(:,:,:)
      real*8 , allocatable :: arrayGlob3D(:,:,:)
      real*8 , allocatable :: arrayGlob2D(:,:)
      integer, allocatable :: arrayGlob2Dint(:,:)
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Const called by ', procID

      if (iAmRootProc) then

          strt2d(:) = (/ 1, rnum_out_const /)

          cnt2d (:) = (/ NETCDF_HDF, 1 /)

          hdr(1) = Nint (gmi_sec)
          hdr(2) = nymd
          hdr(3) = nhms

          call Ncwr_2d_Int (hdr, ncid_const, hdr_var_name, strt2d, cnt2d)

          strt3d(:) = (/ 1, 1, rnum_out_const /)
          cnt3d (:) = (/ numLon, numLat, 1 /)


          strt4d(:) = (/ 1, 1, 1, rnum_out_const /)
          cnt4d (:) = (/ numLon, numLat, numVert, 1 /)

          allocate(arrayGlob2D(i1_gl:i2_gl,ju1_gl:j2_gl))
          allocate(arrayGlob2Dint(i1_gl:i2_gl,ju1_gl:j2_gl))
          allocate(arrayGlob3D(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))
      end if

      if (pr_psf) then
          call subDomain2Global (arrayGlob2D, psf_nc, &
     &            i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, &
     &            rootProc, procID, map1_u, numDomains, SG_PSF, commuWorld)

         if (iAmRootProc)  &
     &      call Ncwr_3d (arrayGlob2D, ncid_const, 'psf', strt3d, cnt3d)
      end if

      if (pr_tropopausePress) then
         call subDomain2Global (arrayGlob2D, tropopausePress_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, &
     &           rootProc, procID, map1_u, numDomains, SG_TROPP, commuWorld)

         if (iAmRootProc)  &
     &      call Ncwr_3d (arrayGlob2D, ncid_const, 'tropopausePress', strt3d, cnt3d)
      end if

      strt4d(:) = (/ 1, 1, 1, rnum_out_const /)
      cnt4d (:) = (/ numLon, numLat, numVert,   1 /)

!      if (pr_emiss_3d) then
!         call subDomain2Global (arrayGlob3D, emiss_3d_out_nc, &
!     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
!     &           rootProc, procID, map1_u, numDomains, SG_EMISS3D, commuWorld)
!
!         if (iAmRootProc)  &
!     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'emiss_3d_out', strt4d, cnt4d)
!
!         if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
!            call subDomain2Global (arrayGlob3D, aerosolEmiss3D_nc, &
!     &              i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
!     &              rootProc, procID, map1_u, numDomains, SG_EMISS3D, commuWorld)
!
!            if (iAmRootProc)  &
!     &         call Ncwr_4d (arrayGlob3D, ncid_const, 'aerosolEmiss3d', strt4d, cnt4d)
!         end if
!      endif

      if (pr_overheadO3col) then
         call subDomain2Global (arrayGlob3D, overheadO3col_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u,numDomains, SG_OvO3col, commuWorld)

         if (iAmRootProc)  &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'overheadO3col', strt4d, cnt4d)
      end if

      if (pr_potentialVorticity) then
         call subDomain2Global (arrayGlob3D, potentialVorticity_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_PV, commuWorld)

         if (iAmRootProc)  &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'potentialVorticity', strt4d, cnt4d)
      end if

      if (pr_kel) then
         call subDomain2Global (arrayGlob3D, kel_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_KEL, commuWorld)

         if (iAmRootProc)  &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'kel', strt4d, cnt4d)
      end if

      if (pr_mass) then
         call subDomain2Global (arrayGlob3D, mass_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_MASS, commuWorld)

         if (iAmRootProc)  &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'mass', strt4d, cnt4d)
      end if

      if (pr_grid_height) then
         call subDomain2Global (arrayGlob3D, gridbox_height_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_GBH, commuWorld)

         if (iAmRootProc)  &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'gridbox_height', strt4d, cnt4d)
      end if

      if (pr_relHumidity) then
         call subDomain2Global (arrayGlob3D, relHumidity_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_RH, commuWorld)

         if (iAmRootProc)  &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'relHumidity', strt4d, cnt4d)
      end if

      if (pr_metwater) then
         call subDomain2Global (arrayGlob3D, metwater_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_METWATER, commuWorld)

         if (iAmRootProc)  &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'metwater', strt4d, cnt4d)
      end if

      if (pr_sulf_src) then

         allocate(loc_dms_oh(i1:i2, ju1:j2, k1:k2))
         allocate(loc_dms_no3(i1:i2, ju1:j2, k1:k2))
         allocate(loc_so2_oh(i1:i2, ju1:j2, k1:k2))
         allocate(loc_so2_h2o2(i1:i2, ju1:j2, k1:k2))
         allocate(loc_so2_o3(i1:i2, ju1:j2, k1:k2))

         call Get_so2_o3(Chemistry, loc_so2_o3)
         call Get_dms_oh(Chemistry, loc_dms_oh)
         call Get_so2_h2o2(Chemistry, loc_so2_h2o2)
         call Get_so2_oh(Chemistry, loc_so2_oh)
         call Get_dms_no3(Chemistry, loc_dms_no3)

         dms_oh  (:,:,:) = loc_dms_oh  (:,:,kx1:kx2)
         dms_no3 (:,:,:) = loc_dms_no3 (:,:,kx1:kx2)
         so2_oh  (:,:,:) = loc_so2_oh  (:,:,kx1:kx2)
         so2_o3  (:,:,:) = loc_so2_o3  (:,:,kx1:kx2)
         so2_h2o2(:,:,:) = loc_so2_h2o2(:,:,kx1:kx2)

         call subDomain2Global (arrayGlob3D, dms_oh, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_DMS_OH, commuWorld)

         if (iAmRootProc)  &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'dms_oh', strt4d, cnt4d)

         call subDomain2Global (arrayGlob3D, dms_no3, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_DMS_NO3, commuWorld)

         if (iAmRootProc)  &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'dms_no3', strt4d, cnt4d)

        call subDomain2Global (arrayGlob3D, so2_oh, &
     &          i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &          rootProc, procID, map1_u, numDomains, SG_SO2_OH, commuWorld)

         if (iAmRootProc)  &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'so2_oh', strt4d, cnt4d)

        call subDomain2Global (arrayGlob3D, so2_h2o2, &
     &          i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &          rootProc, procID, map1_u, numDomains, SG_SO2_H2O2, commuWorld)

         if (iAmRootProc)  &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'so2_h2o2', strt4d, cnt4d)

        call subDomain2Global (arrayGlob3D, so2_o3, &
     &          i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &     rootProc, procID, map1_u, numDomains, SG_SO2_O3, commuWorld)

         if (iAmRootProc)  &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'so2_o3', strt4d, cnt4d)

        deallocate(loc_dms_no3)
        deallocate(loc_dms_oh)
        deallocate(loc_so2_h2o2)
        deallocate(loc_so2_oh)
        deallocate(loc_so2_o3)

        allocate(zero3d(i1:i2, ju1:j2, k1:k2))
        zero3d = 0.0d0
        call Set_dms_oh  (Chemistry, zero3d)
        call Set_dms_no3 (Chemistry, zero3d)
        call Set_so2_oh  (Chemistry, zero3d)
        call Set_so2_h2o2(Chemistry, zero3d)
        call Set_so2_o3  (Chemistry, zero3d)
        deallocate(zero3d)
      end if

      if (do_gcr) then

         call subDomain2Global (arrayGlob3D, GCR_NOx_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_GCR, commuWorld)

         if (iAmRootProc) &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'GCR_NOx_nc', strt4d, cnt4d)

      endif

      if (lightning_opt == 1) then

         call subDomain2Global (arrayGlob3D, lightning_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_LIGHTNING, commuWorld)

         if (iAmRootProc) &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'lightning_nc', strt4d, cnt4d)

         call subDomain2Global (arrayGlob3D, cmf_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_CMF, commuWorld)

         if (iAmRootProc) &
     &      call Ncwr_4d (arrayGlob3D, ncid_const, 'cmf_nc', strt4d, cnt4d)

         call subDomain2Global (arrayGlob2D, flashrate_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, &
     &           rootProc, procID, map1_u, numDomains, SG_FLASHRATE, commuWorld)

         if (iAmRootProc) &
     &      call Ncwr_3d (arrayGlob2D, ncid_const, 'flashrate_nc', strt3d, cnt3d)

         call subDomain2Global (arrayGlob2D, cldmas0_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, &
     &           rootProc, procID, map1_u, numDomains, SG_CLDMASS, commuWorld)

         if (iAmRootProc) &
     &      call Ncwr_3d (arrayGlob2D, ncid_const, 'cldmas0_nc', strt3d, cnt3d)

         call subDomain2Global (arrayGlob2D, cmi_flags_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, &
     &           rootProc, procID, map1_u, numDomains, SG_CMI1, commuWorld)

         if (iAmRootProc) then
            arrayGlob2Dint(:,:) = arrayGlob2D(:,:) + 0.1
            call Ncwr_3d_Int (arrayGlob2Dint, ncid_const, 'cmi_flags_nc', strt3d, cnt3d)
         end if

         call subDomain2Global (arrayGlob2D, cmi_flags1_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2,  &
     &           rootProc, procID, map1_u, numDomains, SG_CMI2, commuWorld)

         if (iAmRootProc) then
            arrayGlob2Dint(:,:) = arrayGlob2D(:,:) + 0.1
            call Ncwr_3d_Int (arrayGlob2Dint, ncid_const, 'cmi_flags1_nc', strt3d, cnt3d)
         end if

      end if

      cnt4d (:) = (/ numLon, numLat, 1, 1 /)

      if (pr_dry_depos) then

         if (pr_drydep_all) then
            Numb_Spec = numSpecies
         else
            Numb_Spec = num_drydep_outrecs
         endif

         if (iAmRootProc) then
            deallocate(arrayGlob3D)
            allocate(arrayGlob3D(numLon, numLat, Numb_Spec))
         end if

         call subDomain2Global (arrayGlob3D, ddep_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, 1, Numb_Spec, &
     &           rootProc, procID, map1_u, numDomains, SG_DRY, commuWorld)

         if (iAmRootProc) then
            if (pr_drydep_all) then
               do ic = 1, num_const_outrecs
                  strt4d(:) = (/ 1, 1, ic, rnum_out_const /)

                  call Ncwr_4d (arrayGlob3D(:,:,const_outrec_map(ic)), ncid_const,  &
     &                    'dry_depos', strt4d, cnt4d)
               end do
            else
               do ic = 1, num_drydep_outrecs
                  strt4d(:) = (/ 1, 1, ic, rnum_out_const /)

                  call Ncwr_4d (arrayGlob3D(:,:,ic), ncid_const, 'dry_depos', &
     &                          strt4d, cnt4d)
               end do
            end if
         end if
      end if

      if (pr_wet_depos) then
         if (pr_wetdep_all) then
            Numb_Spec = numSpecies
         else
            Numb_Spec = num_wetdep_outrecs
         endif

         if (iAmRootProc) then
            deallocate(arrayGlob3D)
            allocate(arrayGlob3D(numLon, numLat, Numb_Spec))
         end if

         call subDomain2Global (arrayGlob3D, wdep_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, 1, Numb_Spec, &
     &           rootProc, procID, map1_u, numDomains, SG_WET, commuWorld)

         if (iAmRootProc) then
            if (pr_wetdep_all) then
               do ic = 1, num_const_outrecs
                  strt4d(:) = (/ 1, 1, ic, rnum_out_const /)

                  call Ncwr_4d (arrayGlob3D(:,:,const_outrec_map(ic)), ncid_const,  &
     &                    'wet_depos', strt4d, cnt4d)
               end do
            else
               do ic = 1, num_wetdep_outrecs
                  strt4d(:) = (/ 1, 1, ic, rnum_out_const /)

                  call Ncwr_4d (arrayGlob3D(:,:,ic), ncid_const, 'wet_depos', &
     &                          strt4d, cnt4d)
               end do
            end if
         end if

      end if

      if (pr_surf_emiss) then

         if (pr_emiss_all) then
            Numb_Spec = numSpecies
         else
            Numb_Spec = num_emiss_outrecs
         endif

         if (iAmRootProc) then
            deallocate(arrayGlob3D)
            allocate(arrayGlob3D(numLon, numLat, Numb_Spec))
         end if

         call subDomain2Global (arrayGlob3D, semiss_out_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, 1, Numb_Spec, &
     &           rootProc, procID, map1_u, numDomains, SG_EMISS, commuWorld)

         if (iAmRootProc) then
            if (pr_emiss_all) then
               do ic = 1, num_const_outrecs
                  strt4d(:) = (/ 1, 1, ic, rnum_out_const /)

                  call Ncwr_4d (arrayGlob3D(:,:,const_outrec_map(ic)), ncid_const,  &
     &                          'surf_emiss', strt4d, cnt4d)
               end do
            else
               do ic = 1, num_emiss_outrecs
                  strt4d(:) = (/ 1, 1, ic, rnum_out_const /)

                  call Ncwr_4d (arrayGlob3D(:,:,ic), ncid_const,  &
     &                          'surf_emiss', strt4d, cnt4d)
               end do
            end if
         end if

         allocate(zero3d(i1:i2,ju1:j2,numSpecies))
         zero3d = 0.d0
         call Set_surf_emiss_out(Emission, zero3d)
         deallocate(zero3d)

         Numb_Spec = 6

         if (iAmRootProc) then
            deallocate(arrayGlob3D)
            allocate(arrayGlob3D(numLon, numLat, Numb_Spec))
         end if

         call subDomain2Global (arrayGlob3D, semiss_out2_nc, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, 1, Numb_Spec, &
     &           rootProc, procID, map1_u, numDomains, SG_EMISS2, commuWorld)

         if (iAmRootProc) then
            do ic = 1, Numb_Spec
               strt4d(:) = (/ 1, 1, ic, rnum_out_const /)
               call Ncwr_4d (arrayGlob3D(:,:,ic), ncid_const,  &
     &           'surf_emiss2', strt4d, cnt4d)
            end do
         end if

         allocate(zero3d(i1:i2,ju1:j2,Numb_Spec))
         zero3d = 0.d0
         call Set_surf_emiss_out2(Emission, zero3d)
         deallocate(zero3d)

         if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
            Numb_Spec = numDust+numAero+5

            if (iAmRootProc) then
               deallocate(arrayGlob3D)
               allocate(arrayGlob3D(numLon, numLat, Numb_Spec))
            end if

            call subDomain2Global (arrayGlob3D, aerosolSurfEmiss_nc, &
     &               i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, &
     &               1, Numb_Spec, rootProc, procID, map1_u, &
     &               numDomains, SG_EMISS2, commuWorld)

            if (iAmRootProc) then
               do ic = 1, Numb_Spec
                  strt4d(:) = (/ 1, 1, ic, rnum_out_const /)
                  call Ncwr_4d (arrayGlob3D(:,:,ic), ncid_const,  &
     &              'aerosolSurfEmiss', strt4d, cnt4d)
               end do
            end if

             allocate(zero3d(i1:i2,ju1:j2,Numb_Spec))
             zero3d = 0.d0
             call Set_aerosolSurfEmiss(Emission, zero3d)
             deallocate(zero3d)
          end if

      end if

      if (iAmRootProc) then
         call Ncdo_Sync (ncid_const)
      end if

      return

      end subroutine Write_Netcdf_Const
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutput_Const
!
! !INTERFACE:
!
      subroutine bufferOutput_Const (const_var_name, Emission, Chemistry, SpeciesConcentration)
!
! !USES:
      use GmiEmissionMethod_mod , only : Get_emiss_3d_out, Get_aerosolEmiss3D
      use GmiEmissionMethod_mod , only : Set_emiss_3d_out, Set_aerosolEmiss3D
      use GmiBuffer_mod         , only : bufferOutput4d_Map, bufferOutput4d_Emiss
      use GmiBuffer_mod         , only : bufferOutput4d_Nomean
      use m_netcdf_io_create    , only : Ncdo_Sync
      use GmiNcOutputSlice_mod  , only : writeSlice4d
!
      implicit none
!
! !INPUT PARAMETERS:
      character (len=*), intent(in) :: const_var_name
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_Emission), intent(inOut) :: Emission
      type(t_Chemistry), intent(inOut) :: Chemistry
!
! !DESCRIPTION:
!   This routine buffers the const data to the output file.  The Slaves send
!   a slice of the data to the Master, the Master writes it out, they then
!   send another slice and it is written out, and so on.
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_CONST  = 4050
      integer, parameter :: SG_EMISS1 = 4051
      integer, parameter :: SG_EMISS2 = 4052
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_VAR_NAME), parameter :: emiss3d_var_name = 'emiss_3d_out'
      real*8, allocatable :: loc_emiss_3d_out   (:,:,:,:)
      real*8, allocatable :: emiss_3d_out   (:,:,:,:)
      real*8, allocatable :: loc_aerosolEmiss3D (:,:,:,:)
      real*8, allocatable :: aerosolEmiss3D (:,:,:,:)
      real*8, allocatable :: decay_3d_out (:,:,:,:)
      real*8, allocatable :: arrayGlob3D (:,:,:)
      integer :: ic, icx
      type (t_GmiArrayBundle), pointer :: concentration(:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (num_const_outrecs > 0) then
         if (iAmRootProc) allocate(arrayGlob3D(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))

         call Get_concentration(SpeciesConcentration, concentration)

         do ic = 1, num_const_outrecs
            if (do_mean) then
               const_nc(:,:,:) = const_mean(:,:,:,ic)
            else
               icx = const_outrec_map(ic)
               const_nc(:,:,:) = concentration(icx)%pArray3D(:,:,kx1:kx2)
            end if

            call subDomain2Global (arrayGlob3D, const_nc, &
     &               i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, &
     &               kx1, kx2, rootProc, procID, map1_u, &
     &               numDomains, SG_CONST, commuWorld)

            if (iAmRootProc) then
               call writeSlice4d (i1_gl, i2_gl, ju1_gl, j2_gl, kx1, kx2, &
     &                   arrayGlob3D, ncid_const, const_var_name, ic, rnum_out_const)
            end if

            call synchronizeGroup (commuWorld)
         end do

         if (iAmRootProc) then
            deallocate(arrayGlob3D)

            call Ncdo_Sync (ncid_const)
         end if
      end if

!... write out radioacive decay diagnostic for tracer package
      if (pr_decay) then
         if (iAmRootProc) allocate(arrayGlob3D(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))

         allocate(decay_3d_out(i1:i2, ju1:j2, k1:k2, 1:numSpecies))
         call Get_decay_3d_out(Chemistry, decay_3d_out)

         do ic = 1, numSpecies
            decay_3d_out_nc(:,:,:) = decay_3d_out(:,:,kx1:kx2,ic)

            call subDomain2Global (arrayGlob3D, decay_3d_out_nc, &
     &                  i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, &
     &                  kx1, kx2, rootProc, procID, map1_u, &
     &                  numDomains, SG_CONST, commuWorld)

            if (iAmRootProc) then
               call writeSlice4d (i1_gl, i2_gl, ju1_gl, j2_gl, kx1, kx2, &
     &                   arrayGlob3D, ncid_const, 'decay_3d_out', ic, rnum_out_const)
            end if
            call synchronizeGroup (commuWorld)
         end do

         if (iAmRootProc) then
            deallocate(arrayGlob3D)

            call Ncdo_Sync (ncid_const)
         end if

         decay_3d_out = 0.0d0
         call Set_decay_3d_out(Chemistry, decay_3d_out)
         deallocate(decay_3d_out)
      endif

!!!! End Jules

      if (pr_emiss_3d) then
         allocate(emiss_3d_out(i1:i2, ju1:j2, kx1:kx2, numSpecies))
         emiss_3d_out = 0.0d0
         allocate(loc_emiss_3d_out(i1:i2, ju1:j2, k1:k2, numSpecies))
         call Get_emiss_3d_out(Emission, loc_emiss_3d_out)
         emiss_3d_out(:,:,:,:) = loc_emiss_3d_out(:,:,kx1:kx2,:)

         loc_emiss_3d_out = 0.0d0
         call Set_emiss_3d_out(Emission, loc_emiss_3d_out)
         deallocate(loc_emiss_3d_out)

         if (num_emiss_outrecs > 0) then
            call bufferOutput4d_Emiss (emiss3d_var_name, kx1, kx2, commuWorld, &
     &                 SG_EMISS1, ncid_const, numSpecies, rnum_out_const,  &
     &                 num_emiss_outrecs, emiss_outrec_map, emiss_3d_out,  &
     &                 emiss_3d_out_nc, map1_u, numDomains, rootProc, procID, &
     &                 i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)
         endif

         if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
            allocate(aerosolEmiss3D(i1:i2, ju1:j2, kx1:kx2, 5))
            aerosolEmiss3D = 0.0d0

            allocate(loc_aerosolEmiss3D(i1:i2, ju1:j2, k1:k2, 5))
            call Get_aerosolEmiss3D(Emission, loc_aerosolEmiss3D)

            aerosolEmiss3D(:,:,:, 5) = loc_aerosolEmiss3D(:,:, kx1:kx2, 5)

            loc_aerosolEmiss3D = 0.0d0
            call Set_aerosolEmiss3D(Emission, loc_aerosolEmiss3D)
            deallocate(loc_aerosolEmiss3D)

            call bufferOutput4d_Nomean ('aerosolEmiss3d', kx1, kx2, commuWorld, &
     &                 SG_EMISS2, ncid_const, 5, rnum_out_const, &
     &                 aerosolEmiss3D, aerosolEmiss3D_nc, map1_u, numDomains, &
     &                 rootProc, procID, &
     &                 i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)
         endif
      endif

      return

      end subroutine bufferOutput_Const
!EOC
!------------------------------------------------------------------------------

end module GmiControlConcentration_mod
