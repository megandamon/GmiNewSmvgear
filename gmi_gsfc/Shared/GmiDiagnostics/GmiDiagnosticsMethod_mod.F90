!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiDiagnosticsMethod_mod
!
! !INTERFACE:
!
#include "GmiESMF_ErrLog.h"
!
      module GmiDiagnosticsMethod_mod
!
! !USES:
      use ESMF_Mod, only : ESMF_Config, ESMF_MAXSTR, ESMF_ConfigGetAttribute
      use GmiESMF_ErrorChecking_mod
      use GmiESMFrcFileReading_mod, only : rcEsmfReadTable, rcEsmfReadLogical, &
     &       reconstructPhrase
      use GmiGrid_mod, only : t_gmiGrid, Get_i1, Get_i2, Get_ju1, Get_jv1,     &
     &                        Get_j2, Get_k1, Get_k2, Get_i1_gl, Get_i2_gl,    &
     &                        Get_ju1_gl, Get_j2_gl, Get_k1_gl,    &
     &                        Get_k2_gl, Get_ilo, Get_ihi, Get_julo, Get_jvlo, &
     &                        Get_jhi, Get_ilo_gl, Get_ihi_gl, Get_julo_gl,    &
     &                        Get_jvlo_gl, Get_jhi_gl, Get_numSpecies
      use GmiReadList_mod, only : Read_List
      use GmiFileOperations_mod, only : finishFileName, makeOutfileName
      use GmiMessagePassing_mod, only : stopCode
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_procID
      use m_netcdf_io_open , only : Ncop_Rd
      use m_netcdf_io_close , only : Nccl
      use m_netcdf_io_get_dimlen , only : Ncget_Unlim_Dimlen

      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: initializeDiagnostics
      public  :: Get_pr_netcdf
      public  :: Get_do_ftiming, Get_problem_name, Get_pr_diag
      public  :: Get_hdr_var_name, Get_hdf_dim_name, Get_rec_dim_name
      public  :: Get_lat_dim_name, Get_lon_dim_name, Get_prs_dim_name
      public  :: Get_spc_dim_name, Get_tim_dim_name
      public  :: Get_constOutputFrequency, Get_aerdustOutputFrequency
      public  :: Get_qjOutputFrequency, Get_qkOutputFrequency, Get_qqjkOutputFrequency
      public  :: Get_sadOutputFrequency, Get_tendOutputFrequency
      public  :: Get_cloudOutputFrequency, Get_fluxOutputFrequency
      public  :: Get_do_mean, Get_k1r_gl, Get_k2r_gl
      public  :: Get_pr_const, Get_pr_sulf_src, Get_pr_surf_emiss
      public  :: Get_pr_emiss_3d, Get_pr_dry_depos, Get_pr_wet_depos
      public  :: Get_num_tend_outrecs, Get_do_aerocom, Get_do_dust_emiss
      public  :: Get_num_emiss_outrecs, Get_num_const_outrecs
      public  :: Get_num_wetdep_outrecs, Get_num_drydep_outrecs
      public  :: Get_pr_mass, Get_pr_kel, Get_pr_psf, Get_pr_grid_height

      public  :: Get_outmain_name, Get_pr_tend, Get_tend_outrec_map
      public  :: Get_pr_overheadO3col, Get_pr_tropopausePress
      public  :: Get_pr_decay
      public  :: Get_pr_potentialVorticity, Get_pr_metwater
      public  :: Get_pr_relHumidity, Get_pr_const_surface, Get_pr_emiss_all
      public  :: Get_pr_drydep_all, Get_pr_wetdep_all, Get_const_outrec_map
      public  :: Get_emiss_outrec_map, Get_drydep_outrec_map
      public  :: Get_wetdep_outrec_map, Get_pr_const_column, Get_pr_level_all

      public  :: Get_sad_dim_name, Get_sad_var_name
      public  :: Get_qj_dim_name, Get_qj_var_name
      public  :: Get_qk_dim_name, Get_qk_var_name
      public  :: Get_qqj_dim_name, Get_qqj_var_name
      public  :: Get_qqk_dim_name, Get_qqk_var_name
      public  :: Get_pr_sad, Get_pr_smv2
      public  :: Get_pr_qj, Get_pr_qk, Get_pr_qqjk, Get_pr_cloud
      public  :: Get_pr_qj_o3_o1d, Get_pr_qj_opt_depth, Get_do_qqjk_inchem
      public  :: Set_do_qqjk_reset, Get_do_qqjk_reset

      public  :: Get_outaerdust_name, Get_AerDust_var_name, Get_pr_AerDust

      public  :: Get_pr_flux, Get_pr_psf_flux, Get_pr_const_flux
      public  :: Get_do_flux_reset, Set_do_flux_reset, Get_flux_name
      public  :: Get_do_day1_flux, Get_do_mean_flux, Get_flux_species
      public  :: Get_flux_species_num, Get_flux_var_name, Set_flux_species

      public  :: Get_pr_overpass1, Get_pr_kel_overpass1
      public  :: Get_pr_psf_overpass1, Get_pr_qj_overpass1
      public  :: Get_pr_qqjk_overpass1, Get_pr_const_overpass1
      public  :: Get_pr_metwater_overpass1, Get_pr_totalMass_overpass1
      public  :: Get_pr_relHumidity_overpass1, Get_pr_gridBoxHeight_overpass1
      public  :: Get_pr_cloudOptDepth_overpass1, Get_pr_overheadO3col_overpass1
      public  :: Get_pr_cloudFraction_overpass1, Get_pr_tropopausePress_overpass1
      public  :: Get_pr_lightningNO_overpass1
      public  :: Get_numSpecies_overpass1, Get_species_overpass1
      public  :: Get_begTime_overpass1, Get_endTime_overpass1
      public  :: Get_pr_overpass1_period

      public  :: Get_pr_overpass2, Get_pr_kel_overpass2
      public  :: Get_pr_psf_overpass2, Get_pr_qj_overpass2
      public  :: Get_pr_qqjk_overpass2, Get_pr_const_overpass2
      public  :: Get_pr_metwater_overpass2, Get_pr_totalMass_overpass2
      public  :: Get_pr_relHumidity_overpass2, Get_pr_gridBoxHeight_overpass2
      public  :: Get_pr_cloudOptDepth_overpass2, Get_pr_overheadO3col_overpass2
      public  :: Get_pr_cloudFraction_overpass2, Get_pr_tropopausePress_overpass2
      public  :: Get_pr_lightningNO_overpass2
      public  :: Get_numSpecies_overpass2, Get_species_overpass2
      public  :: Get_begTime_overpass2, Get_endTime_overpass2
      public  :: Get_pr_overpass2_period

      public  :: Get_pr_overpass3, Get_pr_kel_overpass3
      public  :: Get_pr_psf_overpass3, Get_pr_qj_overpass3
      public  :: Get_pr_qqjk_overpass3, Get_pr_const_overpass3
      public  :: Get_pr_metwater_overpass3, Get_pr_totalMass_overpass3
      public  :: Get_pr_relHumidity_overpass3, Get_pr_gridBoxHeight_overpass3
      public  :: Get_pr_cloudOptDepth_overpass3, Get_pr_overheadO3col_overpass3
      public  :: Get_pr_cloudFraction_overpass3, Get_pr_tropopausePress_overpass3
      public  :: Get_pr_lightningNO_overpass3
      public  :: Get_numSpecies_overpass3, Get_species_overpass3
      public  :: Get_begTime_overpass3, Get_endTime_overpass3
      public  :: Get_pr_overpass3_period

      public  :: Get_pr_overpass4, Get_pr_kel_overpass4
      public  :: Get_pr_psf_overpass4, Get_pr_qj_overpass4
      public  :: Get_pr_qqjk_overpass4, Get_pr_const_overpass4
      public  :: Get_pr_metwater_overpass4, Get_pr_totalMass_overpass4
      public  :: Get_pr_relHumidity_overpass4, Get_pr_gridBoxHeight_overpass4
      public  :: Get_pr_cloudOptDepth_overpass4, Get_pr_overheadO3col_overpass4
      public  :: Get_pr_cloudFraction_overpass4, Get_pr_tropopausePress_overpass4
      public  :: Get_pr_lightningNO_overpass4
      public  :: Get_numSpecies_overpass4, Get_species_overpass4
      public  :: Get_begTime_overpass4, Get_endTime_overpass4
      public  :: Get_pr_overpass4_period

      public  :: Get_rd_restart, Get_pr_restart, Get_do_overwrt_rst
      public  :: Get_restart_inrec, Get_pr_rst_period
      public  :: Get_restart_infile_name

      public  :: Get_freq1_name, Get_freq1_description, Get_pr_const_column_freq1
      public  :: Get_pr_const_surface_freq1, Get_pr_freq1, Get_pr_const_freq1
      public  :: Get_pr_kel_freq1, Get_pr_psf_freq1, Get_pr_tropopausePress_freq1
      public  :: Get_pr_mass_freq1, Get_pr_grid_height_freq1, Get_pr_rel_hum_freq1
      public  :: Get_pr_metwater_freq1, Get_do_mean_freq1, Get_do_day1_freq1
      public  :: Get_do_last_tstep_freq1, Get_pr_overheadO3col_freq1
      public  :: Get_k1_freq1, Get_k2_freq1, Get_freq1_species
      public  :: Get_freq1_species_num, Get_pr_freq1_period
      public  :: Get_pr_at_time_freq1, Get_iRange_freq1, Get_jRange_freq1
      public  :: Get_pr_potentialVorticity_freq1

      public  :: Get_freq2_name, Get_freq2_description, Get_pr_const_column_freq2
      public  :: Get_pr_const_surface_freq2, Get_pr_freq2, Get_pr_const_freq2
      public  :: Get_pr_kel_freq2, Get_pr_psf_freq2, Get_pr_tropopausePress_freq2
      public  :: Get_pr_mass_freq2, Get_pr_grid_height_freq2, Get_pr_rel_hum_freq2
      public  :: Get_pr_metwater_freq2, Get_do_mean_freq2, Get_do_day1_freq2
      public  :: Get_do_last_tstep_freq2, Get_pr_overheadO3col_freq2
      public  :: Get_k1_freq2, Get_k2_freq2, Get_freq2_species
      public  :: Get_freq2_species_num, Get_pr_freq2_period
      public  :: Get_pr_at_time_freq2, Get_iRange_freq2, Get_jRange_freq2
      public  :: Get_pr_potentialVorticity_freq2

      public  :: Get_freq3_name, Get_freq3_description, Get_pr_const_column_freq3
      public  :: Get_pr_const_surface_freq3, Get_pr_freq3, Get_pr_const_freq3
      public  :: Get_pr_kel_freq3, Get_pr_psf_freq3, Get_pr_tropopausePress_freq3
      public  :: Get_pr_mass_freq3, Get_pr_grid_height_freq3, Get_pr_rel_hum_freq3
      public  :: Get_pr_metwater_freq3, Get_do_mean_freq3, Get_do_day1_freq3
      public  :: Get_do_last_tstep_freq3, Get_pr_overheadO3col_freq3
      public  :: Get_k1_freq3, Get_k2_freq3, Get_freq3_species
      public  :: Get_freq3_species_num, Get_pr_freq3_period
      public  :: Get_pr_at_time_freq3, Get_iRange_freq3, Get_jRange_freq3
      public  :: Get_pr_potentialVorticity_freq3

      public  :: Get_freq4_name, Get_freq4_description, Get_pr_const_column_freq4
      public  :: Get_pr_const_surface_freq4, Get_pr_freq4, Get_pr_const_freq4
      public  :: Get_pr_kel_freq4, Get_pr_psf_freq4, Get_pr_tropopausePress_freq4
      public  :: Get_pr_mass_freq4, Get_pr_grid_height_freq4, Get_pr_rel_hum_freq4
      public  :: Get_pr_metwater_freq4, Get_do_mean_freq4, Get_do_day1_freq4
      public  :: Get_do_last_tstep_freq4, Get_pr_overheadO3col_freq4
      public  :: Get_k1_freq4, Get_k2_freq4, Get_freq4_species
      public  :: Get_freq4_species_num, Get_pr_freq4_period
      public  :: Get_pr_at_time_freq4, Get_iRange_freq4, Get_jRange_freq4
      public  :: Get_pr_potentialVorticity_freq4

      public  :: Get_pr_ascii, Get_pr_ascii1, Get_pr_ascii2, Get_pr_ascii3
      public  :: Get_pr_ascii4, Get_pr_ascii5, Get_ascii_out_i, Get_ascii_out_n
      public  :: Get_asclun, Get_pr_ascii_step_interval, Get_pr_time

      public  :: Get_col_diag_site, Get_pr_col_diag, Get_col_diag_num
      public  :: Get_col_diag_pres_num, Get_col_diag_species_num
      public  :: Get_col_diag_species, Get_ncid_col, Get_col_diag_period
      public  :: Get_col_diag_pres, Get_col_diag_lat_lon, Set_col_diag_lat_lon
      public  :: Set_ncid_col

      public  :: t_Diagnostics

#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
#     include "GmiParameters.h"
#     include "setkin_par.h"

      type t_Diagnostics
         character (len=MAX_LENGTH_FILE_NAME) :: problem_name ! problem name

         logical :: do_ftiming          ! do fine timing?

         logical :: pr_diag
         logical :: pr_time

         !------------------
         ! ASCII output file
         !------------------
         logical :: pr_ascii  ! should ASCII output file be written at all?
         logical :: pr_ascii1 ! should first  section of ASCII output file be
                              ! written (species mass data)?
         logical :: pr_ascii2 ! should second section of ASCII output file be
                              ! written (species concentration data)?
         logical :: pr_ascii3 ! should third  section of ASCII output file be
                              ! written (species concentration min/maxs)?
         logical :: pr_ascii4 ! should fourth section of ASCII output file be
                              ! written (total mass of each species)?
         logical :: pr_ascii5 ! should fifth  section of ASCII output file be
                              ! written (total production and loss of each 
                              ! species)
         integer :: ascii_out_i ! longitude index to use in ASCII output file
         integer :: ascii_out_n ! species   index to use in ASCII output file
         integer :: asclun
         integer :: pr_ascii_step_interval  ! interval for ASCII output
         
         logical :: pr_smv2

         character (len=MAX_LENGTH_VAR_NAME) :: hdr_var_name 
         character (len=MAX_LENGTH_VAR_NAME) :: hdf_dim_name ! netCDF header    dimension  name
         character (len=MAX_LENGTH_VAR_NAME) :: lat_dim_name ! netCDF latitude  dimension  name
         character (len=MAX_LENGTH_VAR_NAME) :: lon_dim_name ! netCDF longitude dimension  name
         character (len=MAX_LENGTH_VAR_NAME) :: prs_dim_name ! netCDF pressure  dimension  name
         character (len=MAX_LENGTH_VAR_NAME) :: rec_dim_name ! netCDF record    dimension  name
         character (len=MAX_LENGTH_VAR_NAME) :: spc_dim_name ! netCDF species   dimension  name
         character (len=MAX_LENGTH_VAR_NAME) :: tim_dim_name ! netCDF time      dimension  name
         logical :: pr_netcdf ! should any periodic output files be written at all?
         logical :: do_mean   ! should means or current values be put in periodic output file?

         !-----------------
         ! Main output file
         !-----------------
         character (len=80)  ::  outmain_name
         logical :: pr_const              ! should the species concentrations file be written?
         logical :: pr_kel                ! should temperatures be written?
         logical :: pr_psf                ! Should surface pressures also be written?
         logical :: pr_overheadO3col      ! Should the overhead ozone column be written?
         logical :: pr_decay              ! Should the radioactive decay be written?
         logical :: pr_tropopausePress    ! Should tropopause pressure be written?
         logical :: pr_potentialVorticity ! Should potential vorticity be written?
         logical :: pr_mass               ! Should mass be written?
         logical :: pr_grid_height        ! Should grid box height be written?
         logical :: pr_relHumidity        ! Should relative humidity be written?
         logical :: pr_metwater           ! Should met water be written?
         logical :: pr_dry_depos          ! Should dry depositions be written?
         logical :: pr_wet_depos          ! Should wet depositions be written?
         logical :: pr_surf_emiss         ! should surface emissions be written?
         logical :: pr_emiss_3d           ! should 3D emissions also be written?
         logical :: pr_emiss_all          ! should surface emissions of all the species be written?
         logical :: pr_drydep_all         ! should dry deposition of all the species be written?
         logical :: pr_wetdep_all         ! should wet deposition of all the species be written?
         logical :: pr_sulf_src           ! should sulfur budgets be written?
         logical :: do_aerocom
         logical :: do_dust_emiss
         integer :: num_const_outrecs    ! number of species to output
         integer :: num_emiss_outrecs    ! number of species selected for surface emission output
         integer :: num_drydep_outrecs   ! number of species selected for dry deposition output
         integer :: num_wetdep_outrecs   ! number of species selected for wet deposition output
         integer :: k1r_gl   ! lowest  level for ouput
         integer :: k2r_gl   ! highest level for output
         integer :: const_outrec_map(MAX_NUM_CONST_GIO)
         integer :: emiss_outrec_map(MAX_NUM_CONST_GIO)
         integer :: drydep_outrec_map(MAX_NUM_CONST_GIO)
         integer :: wetdep_outrec_map(MAX_NUM_CONST_GIO)
         real*8  :: constOutputFrequency   ! frequency for the const   output file
         logical :: pr_const_column 
         logical :: pr_const_surface 
         logical :: pr_level_all 

         !---------------------
         ! Station output files
         !---------------------
         character (len=ESMF_MAXSTR)  ::  col_diag_site(MAX_COL_DIAG_SITES) ! name of station locations for column diagnostics
         logical :: pr_col_diag     ! should the periodic column diagnostics be written?
         integer :: col_diag_num         ! number of different locations for column diag.
         integer :: col_diag_pres_num    ! number of pressure levels     for column diag.
         integer :: col_diag_species_num ! number of species             for column diag.
         integer :: col_diag_species(MAX_NUM_CONST_GIO)
         integer :: ncid_col        (MAX_COL_DIAG_SITES)
         real*8  :: col_diag_period ! column diagnostics output period (s)
         real*8  :: col_diag_pres   (MAX_COL_DIAG_PRES) ! pressure levels for column diagnostics
         real*8  :: col_diag_lat_lon(2, MAX_COL_DIAG_SITES) ! lat/lon location of each column diagnostics site


         !-------------------------
         ! Aerosol/Dust output file
         !-------------------------
         character (len=80)  ::  outaerdust_name
         character (len=MAX_LENGTH_VAR_NAME)  ::  AerDust_var_name
         logical :: pr_AerDust ! should the periodic aerosol  diagnostics be written?
         real*8  :: aerdustOutputFrequency ! frequency for the aerdust output file

         !-----------------------
         ! Tendencies output file
         !-----------------------
         logical :: pr_tend        ! should the tendency diagnostics be written?
         logical :: pr_tend_all    ! should the tendency diagnostics be written for all the species?
         real*8  :: tendOutputFrequency    ! frequency for the tend    output file
         integer :: tend_outrec_map (MAX_NUM_CONST_GIO)
         integer :: num_tend_outrecs ! number of species selected for tendencies to output

         !-----------------
         ! Flux output file
         !-----------------
         character (len=80)  ::  flux_name
         logical :: pr_flux        ! should the periodic flux diagnostics be written?
         logical :: do_flux_reset
         logical :: pr_const_flux
         logical :: pr_psf_flux
         logical :: do_day1_flux
         logical :: do_mean_flux
         integer , pointer :: flux_species(:) => null()
         integer :: flux_species_num
         real*8  :: fluxOutputFrequency    ! frequency for the flux    output file
         character(len=MAX_LENGTH_VAR_NAME) :: flux_var_name

         !-------------------------------
         ! Georgia Tech Cloud output file
         !-------------------------------
         logical :: pr_cloud             ! Should cloud output data written?
         real*8  :: cloudOutputFrequency ! frequency for the cloud   output file

         !----------------------
         ! Overpass1 output file
         !----------------------
         logical :: pr_overpass1
         logical :: pr_kel_overpass1
         logical :: pr_psf_overpass1
         logical :: pr_qj_overpass1
         logical :: pr_qqjk_overpass1
         logical :: pr_const_overpass1
         logical :: pr_metwater_overpass1
         logical :: pr_totalMass_overpass1
         logical :: pr_relHumidity_overpass1
         logical :: pr_gridBoxHeight_overpass1
         logical :: pr_cloudOptDepth_overpass1
         logical :: pr_overheadO3col_overpass1
         logical :: pr_cloudFraction_overpass1
         logical :: pr_tropopausePress_overpass1
         logical :: pr_lightningNO_overpass1
         integer :: numSpecies_overpass1
         integer :: species_overpass1 (MAX_NUM_CONST_GIO)
         real*8  :: begTime_overpass1, endTime_overpass1
         real*8  :: pr_overpass1_period

         !----------------------
         ! Overpass2 output file
         !----------------------
         logical :: pr_overpass2
         logical :: pr_kel_overpass2
         logical :: pr_psf_overpass2
         logical :: pr_qj_overpass2
         logical :: pr_qqjk_overpass2
         logical :: pr_const_overpass2
         logical :: pr_metwater_overpass2
         logical :: pr_totalMass_overpass2
         logical :: pr_relHumidity_overpass2
         logical :: pr_gridBoxHeight_overpass2
         logical :: pr_cloudOptDepth_overpass2
         logical :: pr_overheadO3col_overpass2
         logical :: pr_cloudFraction_overpass2
         logical :: pr_tropopausePress_overpass2
         logical :: pr_lightningNO_overpass2
         integer :: numSpecies_overpass2
         integer :: species_overpass2 (MAX_NUM_CONST_GIO)
         real*8  :: begTime_overpass2, endTime_overpass2
         real*8  :: pr_overpass2_period

         !----------------------
         ! Overpass3 output file
         !----------------------
         logical :: pr_overpass3
         logical :: pr_kel_overpass3
         logical :: pr_psf_overpass3
         logical :: pr_qj_overpass3
         logical :: pr_qqjk_overpass3
         logical :: pr_const_overpass3
         logical :: pr_metwater_overpass3
         logical :: pr_totalMass_overpass3
         logical :: pr_relHumidity_overpass3
         logical :: pr_gridBoxHeight_overpass3
         logical :: pr_cloudOptDepth_overpass3
         logical :: pr_overheadO3col_overpass3
         logical :: pr_cloudFraction_overpass3
         logical :: pr_tropopausePress_overpass3
         logical :: pr_lightningNO_overpass3
         integer :: numSpecies_overpass3
         integer :: species_overpass3 (MAX_NUM_CONST_GIO)
         real*8  :: begTime_overpass3, endTime_overpass3
         real*8  :: pr_overpass3_period

         !----------------------
         ! Overpass4 output file
         !----------------------
         logical :: pr_overpass4
         logical :: pr_kel_overpass4
         logical :: pr_psf_overpass4
         logical :: pr_qj_overpass4
         logical :: pr_qqjk_overpass4
         logical :: pr_const_overpass4
         logical :: pr_metwater_overpass4
         logical :: pr_totalMass_overpass4
         logical :: pr_relHumidity_overpass4
         logical :: pr_gridBoxHeight_overpass4
         logical :: pr_cloudOptDepth_overpass4
         logical :: pr_overheadO3col_overpass4
         logical :: pr_cloudFraction_overpass4
         logical :: pr_tropopausePress_overpass4
         logical :: pr_lightningNO_overpass4
         integer :: numSpecies_overpass4
         integer :: species_overpass4 (MAX_NUM_CONST_GIO)
         real*8  :: begTime_overpass4, endTime_overpass4
         real*8  :: pr_overpass4_period

         !---------------
         ! Qj output file
         !---------------
         logical :: pr_qj ! should the periodic qj   output file be written?
         logical :: pr_qj_o3_o1d ! should the o1d reaction be added to the end of the rate constants?
         logical ::  pr_qj_opt_depth ! should the optical depths be added to the end of the rate constants?
         character (len=MAX_LENGTH_VAR_NAME) :: qj_var_name
         character (len=MAX_LENGTH_VAR_NAME) :: qj_dim_name
         real*8  :: qjOutputFrequency      ! frequency for the qj      output file

         !---------------
         ! Qk output file
         !---------------
         logical ::  pr_qk    ! should the periodic qk   output file be written?
         character (len=MAX_LENGTH_VAR_NAME) :: qk_var_name
         character (len=MAX_LENGTH_VAR_NAME) :: qk_dim_name
         real*8  :: qkOutputFrequency      ! frequency for the qk      output file

         !-----------------
         ! Qqjk output file
         !-----------------
         logical ::  pr_qqjk  ! should the periodic qqjk output file be written?
         logical :: do_qqjk_inchem 
         logical :: do_qqjk_reset 
         character (len=MAX_LENGTH_VAR_NAME) :: qqk_var_name
         character (len=MAX_LENGTH_VAR_NAME) :: qqk_dim_name
         character (len=MAX_LENGTH_VAR_NAME) :: qqj_var_name
         character (len=MAX_LENGTH_VAR_NAME) :: qqj_dim_name
         real*8  :: qqjkOutputFrequency    ! frequency for the qqjk    output file

         !-----------------------------------
         ! Surface Area Densities output file
         !-----------------------------------
         logical ::  pr_sad   ! should the periodic sad  output file be written?
         character (len=MAX_LENGTH_VAR_NAME) :: sad_var_name
         character (len=MAX_LENGTH_VAR_NAME) :: sad_dim_name
         real*8  :: sadOutputFrequency     ! frequency for the sad     output file

         !------------------
         ! Freq1 output file
         !------------------
         character (len=80) :: freq1_name
         character (len=80) :: freq1_description
         logical :: pr_const_column_freq1 
         logical :: pr_const_surface_freq1
         logical :: pr_freq1
         logical :: pr_const_freq1
         logical :: pr_kel_freq1 
         logical :: pr_psf_freq1 
         logical :: pr_tropopausePress_freq1
         logical :: pr_potentialVorticity_freq1
         logical :: pr_mass_freq1 
         logical :: pr_grid_height_freq1
         logical :: pr_rel_hum_freq1
         logical :: pr_metwater_freq1
         logical :: do_mean_freq1 
         logical :: do_day1_freq1  
         logical :: do_last_tstep_freq1 
         logical :: pr_overheadO3col_freq1 
         integer :: k1_freq1      ! lowest  level for the freq1 output (>= k1)
         integer :: k2_freq1      ! highest level for the freq1 output (<= k2)
         integer :: freq1_species    (MAX_NUM_CONST_GIO)
         integer :: freq1_species_num
         real*8  :: pr_freq1_period   
         real*8  :: pr_at_time_freq1 
         ! longitude and latitude ranges for freq# output
         ! iRange_freq#(1) = west  longitude grid point
         ! iRange_freq#(2) = east  longitude grid point
         ! jRange_freq#(1) = south latitude  grid point
         ! jRange_freq#(2) = north latitude  grid point
         integer :: iRange_freq1(2), jRange_freq1(2)
!
         !------------------
         ! Freq2 output file
         !------------------
         character (len=80)  ::  freq2_name
         character (len=80)  ::  freq2_description
         logical :: pr_const_column_freq2   
         logical :: pr_const_surface_freq2   
         logical :: pr_freq2            
         logical :: pr_const_freq2     
         logical :: pr_kel_freq2  
         logical :: pr_psf_freq2  
         logical :: pr_tropopausePress_freq2 
         logical :: pr_potentialVorticity_freq2 
         logical :: pr_mass_freq2 
         logical :: pr_grid_height_freq2 
         logical :: pr_rel_hum_freq2 
         logical :: pr_metwater_freq2 
         logical :: do_mean_freq2  
         logical :: do_day1_freq2  
         logical :: do_last_tstep_freq2 
         logical :: pr_overheadO3col_freq2 
         integer :: k1_freq2      ! lowest  level for the freq2 output (>= k1)
         integer :: k2_freq2      ! highest level for the freq2 output (<= k2)
         integer :: freq2_species_num
         integer :: freq2_species    (MAX_NUM_CONST_GIO)
         real*8  :: pr_freq2_period   
         real*8  :: pr_at_time_freq2 
         ! longitude and latitude ranges for freq# output
         ! iRange_freq#(1) = west  longitude grid point
         ! iRange_freq#(2) = east  longitude grid point
         ! jRange_freq#(1) = south latitude  grid point
         ! jRange_freq#(2) = north latitude  grid point
         integer :: iRange_freq2(2), jRange_freq2(2)
!
         !------------------
         ! Freq3 output file
         !------------------
         character (len=80)  ::  freq3_name
         character (len=80)  ::  freq3_description
         logical :: pr_const_column_freq3   
         logical :: pr_const_surface_freq3   
         logical :: pr_freq3            
         logical :: pr_const_freq3     
         logical :: pr_kel_freq3  
         logical :: pr_psf_freq3  
         logical :: pr_tropopausePress_freq3 
         logical :: pr_potentialVorticity_freq3 
         logical :: pr_mass_freq3  
         logical :: pr_grid_height_freq3 
         logical :: pr_rel_hum_freq3 
         logical :: pr_metwater_freq3 
         logical :: do_mean_freq3  
         logical :: do_day1_freq3  
         logical :: do_last_tstep_freq3 
         logical :: pr_overheadO3col_freq3 
         integer :: k1_freq3      ! lowest  level for the freq3 output (>= k1)
         integer :: k2_freq3      ! highest level for the freq3 output (<= k2)
         integer :: freq3_species_num
         integer :: freq3_species    (MAX_NUM_CONST_GIO)
         real*8  :: pr_freq3_period   
         real*8  :: pr_at_time_freq3 
         ! longitude and latitude ranges for freq# output
         ! iRange_freq#(1) = west  longitude grid point
         ! iRange_freq#(2) = east  longitude grid point
         ! jRange_freq#(1) = south latitude  grid point
         ! jRange_freq#(2) = north latitude  grid point
         integer :: iRange_freq3(2), jRange_freq3(2)
!
         !------------------
         ! Freq4 output file
         !------------------
         character (len=80)  ::  freq4_name
         character (len=80)  ::  freq4_description
         logical :: pr_const_column_freq4   
         logical :: pr_const_surface_freq4   
         logical :: pr_freq4            
         logical :: pr_const_freq4     
         logical :: pr_kel_freq4  
         logical :: pr_psf_freq4  
         logical :: pr_tropopausePress_freq4 
         logical :: pr_potentialVorticity_freq4 
         logical :: pr_mass_freq4  
         logical :: pr_grid_height_freq4 
         logical :: pr_rel_hum_freq4 
         logical :: pr_metwater_freq4 
         logical :: do_mean_freq4  
         logical :: do_day1_freq4  
         logical :: do_last_tstep_freq4 
         logical :: pr_overheadO3col_freq4
         integer :: freq4_species_num
         integer :: freq4_species    (MAX_NUM_CONST_GIO)
         integer :: k1_freq4      ! lowest  level for the freq4 output (>= k1)
         integer :: k2_freq4      ! highest level for the freq4 output (<= k2)
         real*8  :: pr_freq4_period   
         real*8  :: pr_at_time_freq4 
         ! longitude and latitude ranges for freq# output
         ! iRange_freq#(1) = west  longitude grid point
         ! iRange_freq#(2) = east  longitude grid point
         ! jRange_freq#(1) = south latitude  grid point
         ! jRange_freq#(2) = north latitude  grid point
         integer :: iRange_freq4(2), jRange_freq4(2)

         !-------------
         ! Restart file
         !-------------
         logical :: rd_restart     ! should a restart file be read?
         integer :: restart_inrec  ! record number in restart input file to read from
         logical :: pr_restart     ! should a restart file be written?
         logical :: do_overwrt_rst ! should the restart file be over-written?
         real*8  :: pr_rst_period  ! restart output period
         character (len=MAX_LENGTH_FILE_NAME) :: restart_infile_name ! restart input file name
      end type t_Diagnostics

!EOP
!=============================================================================
      contains
!=============================================================================
!BOP
!
! !IROUTINE: initializeDiagnostics
!
! !INTERFACE:
!
      subroutine initializeDiagnostics (self, config, gmiGrid, gmiDomain)
!
! !USES:
      use GmiFileUnit_mod, only : GetFileUnitNumber
      use GmiSpeciesRegistry_mod, only : getSpeciesIndex, setNumberSpecies
      use GmiSpeciesRegistry_mod, only : UNKNOWN_SPECIES
      use GmiStringManipulation_mod, only : constructListNames
      use GmiStationFinder_mod     , only : getStationPosition
      use GmiPrintError_mod, only : GmiPrintError

      implicit none

!
! !INPUT PARAMETERS:
      type (t_gmiGrid  ), intent(in) :: gmiGrid
      type (t_gmiDomain), intent(in) :: gmiDomain
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Diagnostics), intent(inOut) :: self
      type (ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
!  Initializes the diagnostics derived type.
!
! !DEFINED PARAMETERS:
      logical, parameter :: GET_NL_PR_DIAG     = .true.
!
! !LOCAL VARIABLES:
      integer :: numSpecies, procID, STATUS, RC, ic
      integer :: i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, k1, k2
      logical :: pr_diag
      character (len=MAX_LENGTH_SPECIES_NAME), allocatable :: tempListNames(:)
!      character (len=MAX_STRING_LENGTH      ) :: emissionSpeciesNames
!      character (len=MAX_STRING_LENGTH      ) :: emissionDustSpeciesNames
!      character (len=MAX_STRING_LENGTH      ) :: emissionAeroSpeciesNames
      character(len=ESMF_MAXSTR) :: IAm, err_msg
      character(len=ESMF_MAXSTR) :: stationsInputFileName
!
      logical :: pr_const_all
      logical :: Tpr_overpass1, Tpr_overpass2
      logical :: Tpr_overpass3, Tpr_overpass4

      integer :: ica
      integer :: ierr, ios
      integer :: max_list_items
      integer :: neg_marker
      integer :: nllun

      integer :: cnt1d (1)
      integer :: strt1d(1)

      integer :: pr_const_rec_flag(MAX_NUM_CONST_GIO)

      real*8  :: pr_nc_period_days
      real*8  :: pr_aerdust_period_days
      real*8  :: pr_qj_period_days
      real*8  :: pr_qk_period_days
      real*8  :: pr_qqjk_period_days
      real*8  :: pr_sad_period_days
      real*8  :: pr_tend_period_days
      real*8  :: pr_flux_period_days
      real*8  :: pr_cloud_period_days
      real*8  :: pr_overpass1_period_days, pr_overpass2_period_days
      real*8  :: pr_overpass3_period_days, pr_overpass4_period_days
      real*8  :: tmet1

      integer :: pr_emiss_rec_flag(MAX_NUM_CONST_GIO)
      integer :: pr_drydep_rec_flag(MAX_NUM_CONST_GIO)
      integer :: pr_wetdep_rec_flag(MAX_NUM_CONST_GIO)
      integer :: pr_tend_rec_flag (MAX_NUM_CONST_GIO)


      character (len=MAX_STRING_LENGTH      ) :: colDiagSpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: dryDepSpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: wetDepSpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: tendSpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: surfEmissSpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: concentrationSpeciesNames

      character (len=MAX_STRING_LENGTH      ) :: fluxSpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: freq1SpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: freq2SpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: freq3SpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: freq4SpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: overpass1SpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: overpass2SpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: overpass3SpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: overpass4SpeciesNames

      character (len=MAX_STRING_LENGTH      ) :: colDiagStationsNames

      logical :: Tpr_freq1
      real*8  :: pr_nc_freq1
      real*8  :: at_time_freq1

      logical :: Tpr_freq2
      real*8  :: pr_nc_freq2
      real*8  :: at_time_freq2

      logical :: Tpr_freq3
      real*8  :: pr_nc_freq3
      real*8  :: at_time_freq3

      logical :: Tpr_freq4
      real*8  :: pr_nc_freq4
      real*8  :: at_time_freq4

      ! longitude and latitude ranges for freq# output
      ! lonRange_freq#(1) = west  longitude
      ! lonRange_freq#(2) = east  longitude
      ! latRange_freq#(1) = south latitude
      ! latRange_freq#(2) = north latitude

      real*8  :: lonRange_freq1(2), latRange_freq1(2)
      real*8  :: lonRange_freq2(2), latRange_freq2(2)
      real*8  :: lonRange_freq3(2), latRange_freq3(2)
      real*8  :: lonRange_freq4(2), latRange_freq4(2)

      real*8  :: ri2_gl, rj2_gl
      real*8  :: latRes, lonRes

!     fffffffffffffffffffffffff

      integer :: ncid_rst
      real*8  :: pr_rst_period_days
      real*8  :: tdt, mdt
      integer :: ndt
!
!EOP
!------------------------------------------------------------------------------
!BOC
      IAm = "initializeDiagnostics"

      call Get_procID(gmiDomain, procID)

      if (procID == 0) then
         if (GET_NL_PR_DIAG) Write (6,*) 'Calling ',IAm
      end if

      call Get_i1_gl (gmiGrid, i1_gl)
      call Get_i2_gl (gmiGrid, i2_gl)
      call Get_ju1_gl (gmiGrid, ju1_gl)
      call Get_j2_gl (gmiGrid, j2_gl)
      call Get_k1 (gmiGrid, k1)
      call Get_k2 (gmiGrid, k2)
      call Get_numSpecies(gmiGrid, numSpecies)

      allocate(tempListNames(numSpecies))

      !################################
      ! Begin reading the resource file
      !################################

      call ESMF_ConfigGetAttribute(config, self%problem_name, &
     &                label   = "problem_name:",&
     &                default = 'gmiExperiment', rc=STATUS )
      VERIFY_(STATUS)

      ! do fine timing?

      call rcEsmfReadLogical(config, self%do_ftiming, "do_ftiming:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%hdr_var_name, &
     &                label   = "hdr_var_name:",&
     &                default = 'hdr', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%hdf_dim_name, &
     &                label   = "hdf_dim_name:",&
     &                default = 'hdf_dim', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%lat_dim_name, &
     &                label   = "lat_dim_name:",&
     &                default = 'latitude_dim', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%lon_dim_name, &
     &                label   = "lon_dim_name:",&
     &                default = 'longitude_dim', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%prs_dim_name, &
     &                label   = "prs_dim_name:",&
     &                default = 'pressure_dim', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%spc_dim_name, &
     &                label   = "spc_dim_name:",&
     &                default = 'species_dim', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%rec_dim_name, &
     &                label   = "rec_dim_name:",&
     &                default = 'rec_dim', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%tim_dim_name, &
     &                label   = "tim_dim_name:",&
     &                default = 'time_dim', rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_diag, "pr_diag:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_time, "pr_time:", &
     &           default=.true., rc=STATUS )
      VERIFY_(STATUS)

      ! ASCII output

      call rcEsmfReadLogical(config, self%pr_ascii, "pr_ascii:", &
     &           default=.true., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_ascii1, "pr_ascii1:", &
     &           default=.true., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_ascii2, "pr_ascii2:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_ascii3, "pr_ascii3:", &
     &           default=.true., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_ascii4, "pr_ascii4:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_ascii5, "pr_ascii5:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%ascii_out_n, &
     &                label   = "ascii_out_n:",&
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%ascii_out_i, &
     &                label   = "ascii_out_i:",&
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

!     ----------------------------------------------
!     pr_ascii_step_interval
!       >0:  ascii output at specified step interval
!       -1:  ascii output at monthly intervals
!     ----------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%pr_ascii_step_interval, &
     &                label   = "pr_ascii_step_interval:",&
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

!     ----------------------
!     General NetCDF output:
!     ----------------------

      ! Output all vertical levels?

      call rcEsmfReadLogical(config, self%pr_level_all, "pr_level_all:", &
     &           default=.true., rc=STATUS )
      VERIFY_(STATUS)

      ! do_mean: should means or current values be put in the
      !          periodic output files?

      call rcEsmfReadLogical(config, self%do_mean, &
     &          "do_mean:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%k1r_gl, &
     &                label   = "k1r_gl:",&
     &                default = k1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%k2r_gl, &
     &                label   = "k2r_gl:",&
     &                default = k2, rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_netcdf, "pr_netcdf:", &
     &           default=.true., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_const, "pr_const:", &
     &           default=.true., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_const_column, "pr_const_column:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_const_surface, "pr_const_surface:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_tropopausePress, &
     &           "pr_tropopausePress:",  default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_potentialVorticity, &
     &           "pr_potentialVorticity:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_psf, "pr_psf:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_kel, "pr_kel:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_overheadO3col, &
     &           "pr_overheadO3col:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_decay, &
     &           "pr_decay:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_mass, "pr_mass:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_relHumidity, "pr_relHumidity:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_grid_height, "pr_grid_height:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_metwater, "pr_metwater:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_dry_depos, "pr_dry_depos:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_wet_depos, "pr_wet_depos:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_surf_emiss, "pr_surf_emiss:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_emiss_3d, "pr_emiss_3d:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_sulf_src, "pr_sulf_src:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%do_aerocom, "do_aerocom:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%do_dust_emiss, "do_dust_emiss:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      !----------------
      ! Overpass Output
      !----------------

      call rcEsmfReadLogical(config, self%pr_overpass1, &
     &          "pr_overpass1:", default=.true., rc=STATUS )
      VERIFY_(STATUS)

      Tpr_overpass1           = .false.

      call rcEsmfReadLogical(config, self%pr_qj_overpass1, &
     &          "qj_overpass1:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_qqjk_overpass1, &
     &          "qqjk_overpass1:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_psf_overpass1, &
     &          "psf_overpass1:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_kel_overpass1, &
     &          "kel_overpass1:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_const_overpass1, &
     &          "const_overpass1:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_metwater_overpass1, &
     &          "metwater_overpass1:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_totalMass_overpass1, &
     &          "totalMass_overpass1:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_relHumidity_overpass1, &
     &          "relHumidity_overpass1:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_gridBoxHeight_overpass1, &
     &          "gridBoxHeight_overpass1:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_cloudOptDepth_overpass1, &
     &          "cloudOptDepth_overpass1:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_cloudFraction_overpass1, &
     &          "cloudFraction_overpass1:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_overheadO3col_overpass1, &
     &          "overheadO3col_overpass1:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_tropopausePress_overpass1, &
     &          "tropopausePress_overpass1:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_lightningNO_overpass1, &
     &          "lightningNO_overpass1:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_overpass2, &
     &          "pr_overpass2:", default=.true., rc=STATUS )
      VERIFY_(STATUS)

      Tpr_overpass2           = .false.

      call rcEsmfReadLogical(config, self%pr_qj_overpass2, &
     &          "qj_overpass2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_qqjk_overpass2, &
     &          "qqjk_overpass2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_psf_overpass2, &
     &          "psf_overpass2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_kel_overpass2, &
     &          "kel_overpass2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_const_overpass2, &
     &          "const_overpass2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_metwater_overpass2, &
     &          "metwater_overpass2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_totalMass_overpass2, &
     &          "totalMass_overpass2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_relHumidity_overpass2, &
     &          "relHumidity_overpass2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_gridBoxHeight_overpass2, &
     &          "gridBoxHeight_overpass2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_cloudOptDepth_overpass2, &
     &          "cloudOptDepth_overpass2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_cloudFraction_overpass2, &
     &          "cloudFraction_overpass2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_overheadO3col_overpass2, &
     &          "overheadO3col_overpass2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_tropopausePress_overpass2, &
     &          "tropopausePress_overpass2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_lightningNO_overpass2, &
     &          "lightningNO_overpass2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_overpass3, &
     &          "pr_overpass3:", default=.true., rc=STATUS )
      VERIFY_(STATUS)

      Tpr_overpass3           = .false.

      call rcEsmfReadLogical(config, self%pr_qj_overpass3, &
     &          "qj_overpass3:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_qqjk_overpass3, &
     &          "qqjk_overpass3:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_psf_overpass3, &
     &          "psf_overpass3:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_kel_overpass3, &
     &          "kel_overpass3:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_const_overpass3, &
     &          "const_overpass3:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_metwater_overpass3, &
     &          "metwater_overpass3:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_totalMass_overpass3, &
     &          "totalMass_overpass3:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_relHumidity_overpass3, &
     &          "relHumidity_overpass3:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_gridBoxHeight_overpass3, &
     &          "gridBoxHeight_overpass3:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_cloudOptDepth_overpass3, &
     &          "cloudOptDepth_overpass3:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_cloudFraction_overpass3, &
     &          "cloudFraction_overpass3:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_overheadO3col_overpass3, &
     &          "overheadO3col_overpass3:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_tropopausePress_overpass3, &
     &          "tropopausePress_overpass3:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_lightningNO_overpass3, &
     &          "lightningNO_overpass3:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_overpass4, &
     &          "pr_overpass4:", default=.true., rc=STATUS )
      VERIFY_(STATUS)

      Tpr_overpass4           = .false.

      call rcEsmfReadLogical(config, self%pr_qj_overpass4, &
     &          "qj_overpass4:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_qqjk_overpass4, &
     &          "qqjk_overpass4:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_psf_overpass4, &
     &          "psf_overpass4:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_kel_overpass4, &
     &          "kel_overpass4:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_const_overpass4, &
     &          "const_overpass4:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_metwater_overpass4, &
     &          "metwater_overpass4:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_totalMass_overpass4, &
     &          "totalMass_overpass4:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_relHumidity_overpass4, &
     &          "relHumidity_overpass4:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_gridBoxHeight_overpass4, &
     &          "gridBoxHeight_overpass4:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_cloudOptDepth_overpass4, &
     &          "cloudOptDepth_overpass4:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_cloudFraction_overpass4, &
     &          "cloudFraction_overpass4:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_overheadO3col_overpass4, &
     &          "overheadO3col_overpass4:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_tropopausePress_overpass4, &
     &          "tropopausePress_overpass4:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_lightningNO_overpass4, &
     &          "lightningNO_overpass4:", default=.false., rc=STATUS )
      VERIFY_(STATUS)


      !------------
      ! Flux Output
      !------------

      call rcEsmfReadLogical(config, self%pr_flux, &
     &          "pr_flux:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_const_flux, &
     &          "pr_const_flux:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_psf_flux, &
     &          "pr_psf_flux:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadTable(config, fluxSpeciesNames, &
     &                     "fluxSpeciesNames::", rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%flux_name, &
     &                label   = "flux_name:", &
     &                default = 'flux', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%flux_var_name, &
     &                label   = "flux_var_name:", &
     &                default = 'mf', rc=STATUS )
      VERIFY_(STATUS)

     !---------------------
     ! qj, qk, qqjk outputs
     !---------------------

      call rcEsmfReadLogical(config, self%pr_qj, &
     &          "pr_qj:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_qk, &
     &          "pr_qk:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_qqjk, &
     &          "pr_qqjk:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_qj_o3_o1d, &
     &          "pr_qj_o3_o1d:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_qj_opt_depth, &
     &          "pr_qj_opt_depth:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%qj_var_name, &
     &                label   = "qj_var_name:", &
     &                default = 'qj', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%qqj_var_name, &
     &                label   = "qqj_var_name:", &
     &                default = 'qqj', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%qj_dim_name, &
     &                label   = "qj_dim_name:", &
     &                default = 'qj_dim', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%qqj_dim_name, &
     &                label   = "qqj_dim_name:", &
     &                default = 'qqj_dim', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%qk_var_name, &
     &                label   = "qk_var_name:", &
     &                default = 'qk', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%qqk_var_name, &
     &                label   = "qqk_var_name:", &
     &                default = 'qqk', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%qk_dim_name, &
     &                label   = "qk_dim_name:", &
     &                default = 'qk_dim', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%qqk_dim_name, &
     &                label   = "qqk_dim_name:", &
     &                default = 'qqk_dim', rc=STATUS )
      VERIFY_(STATUS)

      ! do_qqjk_inchem : if pr_qqjk is on, should qqj's & qqk's be
      !                  determined inside the chemistry solver, or outside?

      call rcEsmfReadLogical(config, self%do_qqjk_inchem, &
     &          "do_qqjk_inchem:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      !-----------
      ! SAD output
      !-----------

      call rcEsmfReadLogical(config, self%pr_sad, &
     &          "pr_sad:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%sad_var_name, &
     &                label   = "sad_var_name:", &
     &                default = 'sad', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%sad_dim_name, &
     &                label   = "sad_dim_name:", &
     &                default = 'sad_dim', rc=STATUS )
      VERIFY_(STATUS)

      !---------------
      ! AerDust output
      !---------------

      call rcEsmfReadLogical(config, self%pr_AerDust, &
     &          "pr_AerDust:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%outaerdust_name, &
     &                label   = "outaerdust_name:", &
     &                default = 'aerdust', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%AerDust_var_name, &
     &                label   = "AerDust_var_name:", &
     &                default = 'OptDepth', rc=STATUS )
      VERIFY_(STATUS)

      !------------------
      ! Tendencies output
      !------------------

      call rcEsmfReadLogical(config, self%pr_tend, &
     &          "pr_tend:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_tend_all, &
     &          "pr_tend_all:", default=.true., rc=STATUS )
      VERIFY_(STATUS)

      pr_tend_rec_flag(1)  = 1
      pr_tend_rec_flag(2:) = INT_DUM_VALUE

      call rcEsmfReadTable(config, tendSpeciesNames, &
     &                     "tendSpeciesNames::", rc=STATUS)

      !-------------
      ! Cloud output
      !-------------

      call rcEsmfReadLogical(config, self%pr_cloud, &
     &          "pr_cloud:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, pr_const_all, &
     &          "pr_const_all:", default=.true., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_emiss_all, &
     &          "pr_emiss_all:", default=.true., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_drydep_all, &
     &          "pr_drydep_all:", default=.true., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_wetdep_all, &
     &          "pr_wetdep_all:", default=.true., rc=STATUS )
      VERIFY_(STATUS)


      self%species_overpass1(:) = 0

      call rcEsmfReadTable(config, overpass1SpeciesNames, &
     &                     "overpass1SpeciesNames::", rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%begTime_overpass1, &
     &                label   = "begTime_overpass1:", &
     &                default = 11.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%endTime_overpass1, &
     &                label   = "endTime_overpass1:", &
     &                default = 13.0d0, rc=STATUS )
!      VERIFY_(STATUS)


      self%species_overpass2(:) = 0

      call rcEsmfReadTable(config, overpass2SpeciesNames, &
     &                     "overpass2SpeciesNames::", rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%begTime_overpass2, &
     &                label   = "begTime_overpass2:", &
     &                default = 11.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%endTime_overpass2, &
     &                label   = "endTime_overpass2:", &
     &                default = 13.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      self%species_overpass3(:) = 0

      call rcEsmfReadTable(config, overpass3SpeciesNames, &
     &                     "overpass3SpeciesNames::", rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%begTime_overpass3, &
     &                label   = "begTime_overpass3:", &
     &                default = 11.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%endTime_overpass3, &
     &                label   = "endTime_overpass3:", &
     &                default = 13.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      self%species_overpass4(:) = 0

      call rcEsmfReadTable(config, overpass4SpeciesNames, &
     &                     "overpass4SpeciesNames::", rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%begTime_overpass4, &
     &                label   = "begTime_overpass4:", &
     &                default = 11.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%endTime_overpass4, &
     &                label   = "endTime_overpass4:", &
     &                default = 13.0d0, rc=STATUS )
!      VERIFY_(STATUS)

!     -----------------------------------------------------
!     pr_const_rec_flag only used if pr_const_all is false.
!     -----------------------------------------------------

      pr_const_rec_flag(1)  = 1
      pr_const_rec_flag(2:) = INT_DUM_VALUE

      call rcEsmfReadTable(config, concentrationSpeciesNames, &
     &                     "concentrationSpeciesNames::", rc=STATUS)

!     -----------------------------------------------------
!     pr_nc_period_days
!       >0.0d0:  NetCDF output at specified interval (days)
!       -1.0d0:  NetCDF output at monthly intervals
!       -2.0d0:  NetCDF output at 1st & 15th of each month
!     -----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, pr_nc_period_days, &
     &                label   = "pr_nc_period_days:", &
     &                default = 1.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, pr_flux_period_days, &
     &                label   = "pr_flux_period_days:", &
     &                default = 1.0d0, rc=STATUS )
!      VERIFY_(STATUS)

!     -----------------------------------------------------
!     period_days
!       >0.0d0:  netCDF output at specified interval (days)
!       -1.0d0:  netCDF output at monthly intervals
!       -2.0d0:  netCDF output at 1st & 15th of each month
!     By default, they will be set to the same frequency as
!     pr_nc_period_days.
!     -----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, pr_qj_period_days, &
     &                label   = "pr_qj_period_days:", &
     &                default = -5.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, pr_qk_period_days, &
     &                label   = "pr_qk_period_days:", &
     &                default = -5.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, pr_sad_period_days, &
     &                label   = "pr_sad_period_days:", &
     &                default = -5.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, pr_qqjk_period_days, &
     &                label   = "pr_qqjk_period_days:", &
     &                default = -5.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, pr_tend_period_days, &
     &                label   = "pr_tend_period_days:", &
     &                default = -5.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, pr_cloud_period_days, &
     &                label   = "pr_cloud_period_days:", &
     &                default = -5.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, pr_aerdust_period_days, &
     &                label   = "pr_aerdust_period_days:", &
     &                default = -5.0d0, rc=STATUS )
!      VERIFY_(STATUS)

!     -----------------------------------------------------
!     pr_overpass_period_days (for overpass variables)
!       >0.0d0:  NetCDF output at specified interval (days)
!       -1.0d0:  NetCDF output at monthly intervals
!       -2.0d0:  NetCDF output at 1st & 15th of each month
!     -----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, pr_overpass1_period_days, &
     &                label   = "pr_overpass1_period_days:", &
     &                default = 1.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, pr_overpass2_period_days, &
     &                label   = "pr_overpass2_period_days:", &
     &                default = 1.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, pr_overpass3_period_days, &
     &                label   = "pr_overpass3_period_days:", &
     &                default = 1.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, pr_overpass4_period_days, &
     &                label   = "pr_overpass4_period_days:", &
     &                default = 1.0d0, rc=STATUS )
!      VERIFY_(STATUS)

      !----------------------------------------------------------------------
      !column diagnostics
      !  colDiagStationNames   : List of the selected sites as a long string
      !  stationsInputFileName : File name containing all the possible sites
      !                          and locations.
      !  colDiagSpeciesNames   : List of the selected species as a long string
      !  col_diag_num          : Number of column diagnostics sites.
      !                          Determine at run time.
      !  col_diag_period       : column diagnostics output period (s)
      !  col_diag_site         : Names of site locations for column diagnostics
      !                          This is construct at run time using the long
      !                          string colDiagStationNames.
      !  col_diag_species      : Should species be included in column diagnostics?
      !                          This is construct at run time using the long
      !                          string colDiagSpeciesNames.
      !  col_diag_pres         : pressure levels for column diagnostics (mb)
      !  col_diag_lat_lon      : lat/lon location of each column diagnostics site
      !                          This is construct at run time using information
      !                          from the input file stationsInputFileName.
      !----------------------------------------------------------------------

      self%ncid_col(:)      = 0
      self%col_diag_num     = 0

      call ESMF_ConfigGetAttribute(config, self%col_diag_period, &
     &                label   = "col_diag_period:", &
     &                default = 3600.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, stationsInputFileName, &
     &                label   = "stationsInputFileName:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

!      call rcEsmfReadTable(config, colDiagStationsNames, &
!     &                     "colDiagStationsNames::", rc=STATUS)

      self%col_diag_site(:)      = ' '
      call rcEsmfReadTable(config, self%col_diag_site, &
     &                     "colDiagStationsNames::", rc=STATUS)

      call rcEsmfReadTable(config, colDiagSpeciesNames, &
     &                     "colDiagSpeciesNames::", rc=STATUS)

      self%col_diag_species(1)   = 0
      self%col_diag_species(2:MAX_NUM_CONST_GIO) = 0

      self%col_diag_pres(1:10)   =  &
     &  (/ 1000.0d0, 900.0d0, 800.0d0, 700.0d0, 600.0d0,  &
     &      500.0d0, 400.0d0, 300.0d0, 200.0d0, 100.0d0 /)

      call rcEsmfReadTable(config, self%col_diag_pres, &
     &                     "col_diag_pres::", rc=STATUS)

      self%col_diag_lat_lon(:,:) = 0.0d0

!     -----------------------------------------------------
!     pr_emiss_rec_flag only used if pr_emiss_all is false.
!     -----------------------------------------------------

      pr_emiss_rec_flag(1)  = 1
      pr_emiss_rec_flag(2:) = INT_DUM_VALUE

      call rcEsmfReadTable(config, surfEmissSpeciesNames, &
     &                     "surfEmissSpeciesNames::", rc=STATUS)

!     -----------------------------------------------------
!     pr_drydep_rec_flag only used if pr_drydep_all is false.
!     -----------------------------------------------------

      pr_drydep_rec_flag(1)  = 1
      pr_drydep_rec_flag(2:) = INT_DUM_VALUE

      call rcEsmfReadTable(config, dryDepSpeciesNames, &
     &                     "dryDepSpeciesNames::", rc=STATUS)

!     -----------------------------------------------------
!     pr_wetdep_rec_flag only used if pr_wetdep_all is false.
!     -----------------------------------------------------

      pr_wetdep_rec_flag(1)  = 1
      pr_wetdep_rec_flag(2:) = INT_DUM_VALUE

      call rcEsmfReadTable(config, wetDepSpeciesNames, &
     &                     "wetDepSpeciesNames::", rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%outmain_name, &
     &                label   = "outmain_name:", &
     &                default = 'const', rc=STATUS )
!      VERIFY_(STATUS)

      !------------------
      ! Frequency Outputs
      !------------------

      ! Set default ranges for longitude/latitude.

      lonRange_freq1 = (/ 0.0d0, 360.0d0 /)
      latRange_freq1 = (/ -90.0d0, 90.0d0 /)

      lonRange_freq2 = (/ 0.0d0, 360.0d0 /)
      latRange_freq2 = (/ -90.0d0, 90.0d0 /)

      lonRange_freq3 = (/ 0.0d0, 360.0d0 /)
      latRange_freq3 = (/ -90.0d0, 90.0d0 /)

      lonRange_freq4 = (/ 0.0d0, 360.0d0 /)
      latRange_freq4 = (/ -90.0d0, 90.0d0 /)

      ! Set the default range of vertical levels for output

      call ESMF_ConfigGetAttribute(config, self%k1_freq1, &
     &                label   = "k1_freq1:", &
     &                default = k1, rc=STATUS )
!      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%k2_freq1, &
     &                label   = "k2_freq1:", &
     &                default = k2, rc=STATUS )
!      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_const_freq1, &
     &           "pr_const_freq1:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_const_column_freq1, &
     &           "pr_const_column_freq1:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_const_surface_freq1, &
     &           "pr_const_surface_freq1:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_tropopausePress_freq1, &
     &           "pr_tropopausePress_freq1:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_potentialVorticity_freq1, &
     &           "pr_potentialVorticity_freq1:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_psf_freq1, &
     &           "pr_psf_freq1:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_kel_freq1, &
     &           "pr_kel_freq1:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_mass_freq1, &
     &           "pr_mass_freq1:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_grid_height_freq1, &
     &           "pr_grid_height_freq1:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_rel_hum_freq1, &
     &           "pr_rel_hum_freq1:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_metwater_freq1, &
     &           "pr_metwater_freq1:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_overheadO3col_freq1, &
     &           "pr_overheadO3col_freq1:", default=.false., rc=STATUS)

      Tpr_freq1       = .false.

      call rcEsmfReadLogical(config, self%do_mean_freq1, &
     &           "do_mean_freq1:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%do_day1_freq1, &
     &           "do_day1_freq1:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_freq1, &
     &           "pr_freq1:", default=.true., rc=STATUS)

      call rcEsmfReadLogical(config, self%do_last_tstep_freq1, &
     &           "do_last_tstep_freq1:", default=.true., rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%freq1_description, &
     &                label   = "freq1_description:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call reconstructPhrase (self%freq1_description)

      call ESMF_ConfigGetAttribute(config, self%freq1_name, &
     &                label   = "freq1_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, pr_nc_freq1, &
     &                label   = "pr_nc_freq1:", &
     &                default = 1.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, at_time_freq1, &
     &                label   = "at_time_freq1:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      self%freq1_species(:) = 0

      call rcEsmfReadTable(config, freq1SpeciesNames, &
     &                     "freq1SpeciesNames::", rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%k1_freq2, &
     &                label   = "k1_freq2:", &
     &                default = k1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%k2_freq2, &
     &                label   = "k2_freq2:", &
     &                default = k2, rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_const_freq2, &
     &           "pr_const_freq2:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_const_column_freq2, &
     &           "pr_const_column_freq2:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_const_surface_freq2, &
     &           "pr_const_surface_freq2:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_tropopausePress_freq2, &
     &           "pr_tropopausePress_freq2:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_potentialVorticity_freq2, &
     &           "pr_potentialVorticity_freq2:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_psf_freq2, &
     &           "pr_psf_freq2:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_kel_freq2, &
     &           "pr_kel_freq2:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_mass_freq2, &
     &           "pr_mass_freq2:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_grid_height_freq2, &
     &           "pr_grid_height_freq2:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_rel_hum_freq2, &
     &           "pr_rel_hum_freq2:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_metwater_freq2, &
     &           "pr_metwater_freq2:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_overheadO3col_freq2, &
     &           "pr_overheadO3col_freq2:", default=.false., rc=STATUS)

      Tpr_freq2       = .false.

      call rcEsmfReadLogical(config, self%do_mean_freq2, &
     &           "do_mean_freq2:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%do_day1_freq2, &
     &           "do_day1_freq2:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_freq2, &
     &           "pr_freq2:", default=.true., rc=STATUS)

      call rcEsmfReadLogical(config, self%do_last_tstep_freq2, &
     &           "do_last_tstep_freq2:", default=.true., rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%freq2_description, &
     &                label   = "freq2_description:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call reconstructPhrase (self%freq2_description)

      call ESMF_ConfigGetAttribute(config, self%freq2_name, &
     &                label   = "freq2_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, pr_nc_freq2, &
     &                label   = "pr_nc_freq2:", &
     &                default = 1.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, at_time_freq2, &
     &                label   = "at_time_freq2:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      self%freq2_species(:) = 0

      call rcEsmfReadTable(config, freq2SpeciesNames, &
     &                     "freq2SpeciesNames::", rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%k1_freq3, &
     &                label   = "k1_freq3:", &
     &                default = k1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%k2_freq3, &
     &                label   = "k2_freq3:", &
     &                default = k2, rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_const_freq3, &
     &           "pr_const_freq3:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_const_column_freq3, &
     &           "pr_const_column_freq3:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_const_surface_freq3, &
     &           "pr_const_surface_freq3:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_tropopausePress_freq3, &
     &           "pr_tropopausePress_freq3:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_potentialVorticity_freq3, &
     &           "pr_potentialVorticity_freq3:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_psf_freq3, &
     &           "pr_psf_freq3:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_kel_freq3, &
     &           "pr_kel_freq3:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_mass_freq3, &
     &           "pr_mass_freq3:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_grid_height_freq3, &
     &           "pr_grid_height_freq3:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_rel_hum_freq3, &
     &           "pr_rel_hum_freq3:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_metwater_freq3, &
     &           "pr_metwater_freq3:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_overheadO3col_freq3, &
     &           "pr_overheadO3col_freq3:", default=.false., rc=STATUS)

      Tpr_freq3       = .false.

      call rcEsmfReadLogical(config, self%do_mean_freq3, &
     &           "do_mean_freq3:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%do_day1_freq3, &
     &           "do_day1_freq3:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_freq3, &
     &           "pr_freq3:", default=.true., rc=STATUS)

      call rcEsmfReadLogical(config, self%do_last_tstep_freq3, &
     &           "do_last_tstep_freq3:", default=.true., rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%freq3_description, &
     &                label   = "freq3_description:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call reconstructPhrase (self%freq3_description)

      call ESMF_ConfigGetAttribute(config, self%freq3_name, &
     &                label   = "freq3_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, pr_nc_freq3, &
     &                label   = "pr_nc_freq3:", &
     &                default = 1.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, at_time_freq3, &
     &                label   = "at_time_freq3:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      self%freq3_species(:) = 0

      call rcEsmfReadTable(config, freq3SpeciesNames, &
     &                     "freq3SpeciesNames::", rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%k1_freq4, &
     &                label   = "k1_freq4:", &
     &                default = k1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%k2_freq4, &
     &                label   = "k2_freq4:", &
     &                default = k2, rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%pr_const_freq4, &
     &           "pr_const_freq4:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_const_column_freq4, &
     &           "pr_const_column_freq4:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_const_surface_freq4, &
     &           "pr_const_surface_freq4:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_tropopausePress_freq4, &
     &           "pr_tropopausePress_freq4:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_potentialVorticity_freq4, &
     &           "pr_potentialVorticity_freq4:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_psf_freq4, &
     &           "pr_psf_freq4:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_kel_freq4, &
     &           "pr_kel_freq4:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_mass_freq4, &
     &           "pr_mass_freq4:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_grid_height_freq4, &
     &           "pr_grid_height_freq4:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_rel_hum_freq4, &
     &           "pr_rel_hum_freq4:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_metwater_freq4, &
     &           "pr_metwater_freq4:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_overheadO3col_freq4, &
     &           "pr_overheadO3col_freq4:", default=.false., rc=STATUS)

      Tpr_freq4       = .false.

      call rcEsmfReadLogical(config, self%do_mean_freq4, &
     &           "do_mean_freq4:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%do_day1_freq4, &
     &           "do_day1_freq4:", default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%pr_freq4, &
     &           "pr_freq4:", default=.true., rc=STATUS)

      call rcEsmfReadLogical(config, self%do_last_tstep_freq4, &
     &           "do_last_tstep_freq4:", default=.true., rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%freq4_description, &
     &                label   = "freq4_description:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call reconstructPhrase (self%freq4_description)

      call ESMF_ConfigGetAttribute(config, self%freq4_name, &
     &                label   = "freq4_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, pr_nc_freq4, &
     &                label   = "pr_nc_freq4:", &
     &                default = 1.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, at_time_freq4, &
     &                label   = "at_time_freq4:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      self%freq4_species(:) = 0

      call rcEsmfReadTable(config, freq4SpeciesNames, &
     &                     "freq4SpeciesNames::", rc=STATUS)

      !================
      ! Restart Section
      !================

      call rcEsmfReadLogical(config, self%pr_restart, "pr_restart:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%rd_restart, "rd_restart:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%do_overwrt_rst, "do_overwrt_rst:", &
     &           default=.true., rc=STATUS )
      VERIFY_(STATUS)

      !------------------------------------------------------
      !pr_rst_period_days
      !  >0.0d0:  restart output at specified interval (days)
      !  -1.0d0:  restart output at monthly intervals
      !  -2.0d0:  restart output at 1st & 15th of each month
      !------------------------------------------------------

      call ESMF_ConfigGetAttribute(config, pr_rst_period_days, &
     &                label   = "pr_rst_period_days:", &
     &                default = 7.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%restart_inrec, &
     &                label   = "restart_inrec:", &
     &                default = INT_DUM_VALUE, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%restart_infile_name, &
     &                label   = "restart_infile_name:", &
     &                default = 'gmi.rst.nc', rc=STATUS )
      VERIFY_(STATUS)

      !##############################
      ! End reading the resource file
      !##############################

      if (.not. self%pr_ascii) then
        self%pr_ascii1    = .false.
        self%pr_ascii2    = .false.
        self%pr_ascii3    = .false.
        self%pr_ascii4    = .false.
        self%pr_ascii5    = .false.
      end if

      if (.not. self%pr_netcdf) then
        self%pr_const      = .false.
        self%pr_flux       = .false.
        self%pr_tend       = .false.
        self%pr_qj         = .false.
        self%pr_qk         = .false.
        self%pr_qqjk       = .false.
        self%pr_sad        = .false.
        self%pr_AerDust    = .false.
      end if

      if (.not. self%pr_const) then
        self%pr_tropopausePress = .false.
        self%pr_potentialVorticity = .false.
        self%pr_psf        = .false.
        self%pr_kel        = .false.
        self%pr_overheadO3col = .false.
        self%pr_decay      = .false.
        self%pr_mass       = .false.
        self%pr_relHumidity = .false.
        self%pr_grid_height = .false.
        self%pr_metwater   = .false.
        self%pr_dry_depos  = .false.
        self%pr_wet_depos  = .false.
        self%pr_surf_emiss = .false.
        self%pr_emiss_3d   = .false.
        self%pr_sulf_src   = .false.
        self%do_aerocom    = .false.
        self%do_dust_emiss = .false.
      end if

      if ((.not. self%pr_const) .and.  &
     &    (.not. self%pr_flux)  .and.  &
     &    (.not. self%pr_tend)) then
        pr_const_all = .false.
        self%pr_emiss_all = .false.
        self%pr_drydep_all = .false.
        self%pr_wetdep_all = .false.
        self%pr_tend_all  = .false.
      end if

!     ------------------------------
!     const/flux/tend Netcdf output.
!     ------------------------------

      self%num_const_outrecs   = 0
      self%const_outrec_map(:) = INT_DUM_VALUE

      if (self%pr_const .or. self%pr_flux .or. self%pr_tend) then
         if (pr_const_all) then
            self%num_const_outrecs = numSpecies
            do ic = 1, numSpecies
               self%const_outrec_map(ic) = ic
            end do
         else
            ! Set the initial value of the list
            tempListNames(:) = ''

            ! Construct the list of names using the long string
            call constructListNames(tempListNames, concentrationSpeciesNames)

            self%num_const_outrecs = count(tempListNames /= '')
            do ic = 1, self%num_const_outrecs
               self%const_outrec_map(ic) = getSpeciesIndex(tempListNames(ic))
            end do
         end if
      end if

      self%num_emiss_outrecs   = 0
      self%emiss_outrec_map(:) = INT_DUM_VALUE

      if (self%pr_surf_emiss) then
         if (self%pr_emiss_all) then
            self%num_emiss_outrecs = count(tempListNames /= '')
            do ic = 1, self%num_emiss_outrecs
               self%emiss_outrec_map(ic) = ic
            end do         
         else 
         ! Set the initial value of the list
            tempListNames(:) = ''

         ! Construct the list of names using the long string
            call constructListNames(tempListNames, surfEmissSpeciesNames)

            self%num_emiss_outrecs = count(tempListNames /= '')
            do ic = 1, self%num_emiss_outrecs
               self%emiss_outrec_map(ic) = getSpeciesIndex(tempListNames(ic))
            end do
         endif
      end if

      self%num_drydep_outrecs   = 0
      self%drydep_outrec_map(:) = INT_DUM_VALUE

      if (self%pr_dry_depos .and. (.not.self%pr_drydep_all)) then
         ! Set the initial value of the list
         tempListNames(:) = ''

         ! Construct the list of names using the long string
         call constructListNames(tempListNames, dryDepSpeciesNames)

         self%num_drydep_outrecs = count(tempListNames /= '')
         do ic = 1, self%num_drydep_outrecs
            self%drydep_outrec_map(ic) = getSpeciesIndex(tempListNames(ic))
         end do
      end if

      self%num_wetdep_outrecs   = 0
      self%wetdep_outrec_map(:) = INT_DUM_VALUE

      if (self%pr_wet_depos .and. (.not.self%pr_wetdep_all)) then
         ! Set the initial value of the list
         tempListNames(:) = ''

         ! Construct the list of names using the long string
         call constructListNames(tempListNames, wetDepSpeciesNames)

         self%num_wetdep_outrecs = count(tempListNames /= '')
         do ic = 1, self%num_wetdep_outrecs
            self%wetdep_outrec_map(ic) = getSpeciesIndex(tempListNames(ic))
         end do
      end if

      self%num_tend_outrecs   = 0
      self%tend_outrec_map(:) = INT_DUM_VALUE

      if (self%pr_tend) then
         if (self%pr_tend_all) then
            self%num_tend_outrecs = self%num_const_outrecs
            self%tend_outrec_map(:) = self%const_outrec_map(:)
         else
            ! Set the initial value of the list
            tempListNames(:) = ''

            ! Construct the list of names using the long string
            call constructListNames(tempListNames, tendSpeciesNames)

            self%num_tend_outrecs = count(tempListNames /= '')
            do ic = 1, self%num_tend_outrecs
               self%tend_outrec_map(ic) = getSpeciesIndex(tempListNames(ic))
            end do
         end if
      end if

!     -------------------------------------------------------
!     Find out if a overpass variable netCDF file will be created
!     -------------------------------------------------------

      ! Set the initial value of the list
      tempListNames(:) = ''

      ! Construct the list of names using the long string
      call constructListNames(tempListNames, overpass1SpeciesNames)

      self%numSpecies_overpass1 = Count (tempListNames(:) /= '')
      if (self%numSpecies_overpass1 > 0) then
         Tpr_overpass1 = .true.
         do ic = 1, self%numSpecies_overpass1 
            self%species_overpass1(getSpeciesIndex(tempListNames(ic))) = 1
         end do
      end if

      if (self%pr_overpass1 .and. Tpr_overpass1 .and. self%do_mean) then
         self%pr_overpass1       = .true.
      else
         self%pr_overpass1       = .false.
         self%pr_const_overpass1 = .false.
      end if

      if (self%pr_overpass1) then
        self%pr_psf_overpass1   = .true.
        self%pr_kel_overpass1   = .true.
        self%pr_tropopausePress_overpass1 = .true.
        self%pr_const_overpass1 = .true.
        self%pr_totalMass_overpass1 = .true.
        self%pr_relHumidity_overpass1 = .true.
        self%pr_gridBoxHeight_overpass1 = .true.
        self%pr_cloudOptDepth_overpass1 = .true.
        self%pr_cloudFraction_overpass1 = .true.
        self%pr_overheadO3col_overpass1 = .true.
        self%pr_tropopausePress_overpass1 = .true.
        if (self%pr_lightningNO_overpass1 /= .false.) then
           self%pr_lightningNO_overpass1 = .true.
        endif
      end if

      ! Set the initial value of the list
      tempListNames(:) = ''

      ! Construct the list of names using the long string
      call constructListNames(tempListNames, overpass2SpeciesNames)

      self%numSpecies_overpass2 = Count (tempListNames(:) /= '')
      if (self%numSpecies_overpass2 > 0) then
         Tpr_overpass2 = .true.
         do ic = 1, self%numSpecies_overpass2
            self%species_overpass2(getSpeciesIndex(tempListNames(ic))) = 1
         end do
      end if

      if (self%pr_overpass2 .and. Tpr_overpass2 .and. self%do_mean) then
         self%pr_overpass2       = .true.
      else
         self%pr_overpass2       = .false.
         self%pr_const_overpass2 = .false.
      end if

      if (self%pr_overpass2) then
        self%pr_psf_overpass2   = .true.
        self%pr_kel_overpass2   = .true.
        self%pr_tropopausePress_overpass2 = .true.
        self%pr_const_overpass2 = .true.
        self%pr_totalMass_overpass2 = .true.
        self%pr_relHumidity_overpass2 = .true.
        self%pr_gridBoxHeight_overpass2 = .true.
        self%pr_cloudOptDepth_overpass2 = .true.
        self%pr_cloudFraction_overpass2 = .true.
        self%pr_overheadO3col_overpass2 = .true.
        self%pr_tropopausePress_overpass2 = .true.
        if (self%pr_lightningNO_overpass2 /= .false.) then
           self%pr_lightningNO_overpass2 = .true.
        endif
      end if

      ! Set the initial value of the list
      tempListNames(:) = ''
      
      ! Construct the list of names using the long string
      call constructListNames(tempListNames, overpass3SpeciesNames)

      self%numSpecies_overpass3 = Count (tempListNames(:) /= '')
      if (self%numSpecies_overpass3 > 0) then
         Tpr_overpass3 = .true.
         do ic = 1, self%numSpecies_overpass3
            self%species_overpass3(getSpeciesIndex(tempListNames(ic))) = 1
         end do
      end if

      if (self%pr_overpass3 .and. Tpr_overpass3 .and. self%do_mean) then
         self%pr_overpass3       = .true.
      else
         self%pr_overpass3       = .false.
         self%pr_const_overpass3 = .false.
      end if

      if (self%pr_overpass3) then
        self%pr_psf_overpass3   = .true.
        self%pr_kel_overpass3   = .true.
        self%pr_tropopausePress_overpass3 = .true.
        self%pr_const_overpass3 = .true.
        self%pr_totalMass_overpass3 = .true.
        self%pr_relHumidity_overpass3 = .true.
        self%pr_gridBoxHeight_overpass3 = .true.
        self%pr_cloudOptDepth_overpass3 = .true.
        self%pr_cloudFraction_overpass3 = .true.
        self%pr_overheadO3col_overpass3 = .true.
        self%pr_tropopausePress_overpass3 = .true.
        self%pr_lightningNO_overpass3 = .true.
      end if

      ! Set the initial value of the list
      tempListNames(:) = ''

      ! Construct the list of names using the long string
      call constructListNames(tempListNames, overpass4SpeciesNames)

      self%numSpecies_overpass4 = Count (tempListNames(:) /= '')
      if (self%numSpecies_overpass4 > 0) then
         Tpr_overpass4 = .true.
         do ic = 1, self%numSpecies_overpass4
            self%species_overpass4(getSpeciesIndex(tempListNames(ic))) = 1
         end do
      end if

      if (self%pr_overpass4 .and. Tpr_overpass4 .and. self%do_mean) then
         self%pr_overpass4       = .true.
      else
         self%pr_overpass4       = .false.
         self%pr_const_overpass4 = .false.
      end if

      if (self%pr_overpass4) then
        self%pr_psf_overpass4   = .true.
        self%pr_kel_overpass4   = .true.
        self%pr_tropopausePress_overpass4 = .true.
        self%pr_const_overpass4 = .true.
        self%pr_totalMass_overpass4 = .true.
        self%pr_relHumidity_overpass4 = .true.
        self%pr_gridBoxHeight_overpass4 = .true.
        self%pr_cloudOptDepth_overpass4 = .true.
        self%pr_cloudFraction_overpass4 = .true.
        self%pr_overheadO3col_overpass4 = .true.
        self%pr_tropopausePress_overpass4 = .true.
       self%pr_lightningNO_overpass4 = .true.
      end if

!     -------------------------------------------------------
!     Find out if a hourly variable netCDF file will be created
!     -------------------------------------------------------
!     Column diagnostic output.
!     -------------------------

      self%pr_col_diag = .false.

      ! Construct the list of station using the long string
      !call constructListNames(self%col_diag_site, colDiagStationsNames)

      self%col_diag_num = Count (self%col_diag_site(:) /= '')

      if (self%col_diag_num /= 0) then

         do ic = 1, self%col_diag_num
            call getStationPosition(self%col_diag_site(ic), &
                                    self%col_diag_lat_lon(1,ic), &
                                    self%col_diag_lat_lon(2,ic), &
                                    trim(stationsInputFileName))
         end do

        self%pr_col_diag = .true.

        self%col_diag_pres_num    = count(self%col_diag_pres(:) /= 0.0d0)
      
        ! Set the initial value of the list
        tempListNames(:) = ''

        ! Construct the list of names using the long string
        call constructListNames(tempListNames, colDiagSpeciesNames)

        self%col_diag_species_num = Count (tempListNames(:) /= '')

        if (self%col_diag_species_num > 0) then
           do ic = 1, self%col_diag_species_num
              self%col_diag_species(getSpeciesIndex(tempListNames(ic))) = 1
           end do
        end if

!       ----------------------------------------------------------
!       If longitude was input as negative (West) shift it to East
!       longitude (0-360).
!       ----------------------------------------------------------
!
        where(self%col_diag_lat_lon(2,1:self%col_diag_num) < 0.0d0)  &
     &    self%col_diag_lat_lon(2,1:self%col_diag_num) =  &
     &    360.0d0 + self%col_diag_lat_lon(2,1:self%col_diag_num)

      end if

!     ffffffffffffffffffffffffffffffffffffffffffffffff

      !-------------------------------------------------
      ! Determine the horizontal domain for freq# output
      !-------------------------------------------------

      ! If longitude was input as negative (West) shift it to East
      ! longitude (0-360)

      where (lonRange_freq1(1:2) < 0.0d0) &
             lonRange_freq1(1:2) = 360.0d0 + lonRange_freq1(1:2)
      where (lonRange_freq2(1:2) < 0.0d0) &
             lonRange_freq2(1:2) = 360.0d0 + lonRange_freq2(1:2)
      where (lonRange_freq3(1:2) < 0.0d0) &
             lonRange_freq3(1:2) = 360.0d0 + lonRange_freq3(1:2)
      where (lonRange_freq4(1:2) < 0.0d0) &
             lonRange_freq4(1:2) = 360.0d0 + lonRange_freq4(1:2)

      ! Get the model horizontal resolution
      ri2_gl = i2_gl
      lonRes = 360.0d0 / ri2_gl

      rj2_gl = j2_gl
      latRes = 180.0d0 / (rj2_gl - 1.0d0)

      self%jRange_freq1(1:2) = (90.0d0 + latRange_freq1(1:2)) / latRes + 1
      self%iRange_freq1(1:2) = lonRange_freq1(1:2) / lonRes + 1
      if (lonRange_freq1(2) == 360.0d0) self%iRange_freq1(2) = self%iRange_freq1(2) - 1

      self%jRange_freq2(1:2) = (90.0d0 + latRange_freq2(1:2)) / latRes + 1
      self%iRange_freq2(1:2) = lonRange_freq2(1:2) / lonRes + 1
      if (lonRange_freq2(2) == 360.0d0) self%iRange_freq2(2) = self%iRange_freq2(2) - 1

      self%jRange_freq3(1:2) = (90.0d0 + latRange_freq3(1:2)) / latRes + 1
      self%iRange_freq3(1:2) = lonRange_freq3(1:2) / lonRes + 1
      if (lonRange_freq3(2) == 360.0d0) self%iRange_freq3(2) = self%iRange_freq3(2) - 1

      self%jRange_freq4(1:2) = (90.0d0 + latRange_freq4(1:2)) / latRes + 1
      self%iRange_freq4(1:2) = lonRange_freq4(1:2) / lonRes + 1
      if (lonRange_freq4(2) == 360.0d0) self%iRange_freq4(2) = self%iRange_freq4(2) - 1

      ! Set the initial value of the list
      tempListNames(:) = ''

      ! Construct the list of names using the long string
      call constructListNames(tempListNames, freq1SpeciesNames)

      self%freq1_species_num = Count (tempListNames(:) /= '')
      if (self%freq1_species_num > 0) then
         Tpr_freq1 = .true.
         do ic = 1, self%freq1_species_num 
            self%freq1_species(getSpeciesIndex(tempListNames(ic))) = 1
         end do
      end if

      if ((self%pr_freq1) .and. (Tpr_freq1)) then
         self%pr_freq1 = .true.
      else
         self%pr_freq1 = .false.
      end if

      self%pr_freq1_period = Nint (pr_nc_freq1 * SECPDY)
      self%pr_at_time_freq1 = Nint (at_time_freq1 * SECPDY)
!

      ! Set the initial value of the list
      tempListNames(:) = ''

      ! Construct the list of names using the long string
      call constructListNames(tempListNames, freq2SpeciesNames)

      self%freq2_species_num = Count (tempListNames(:) /= '')
      if (self%freq2_species_num > 0) then
         Tpr_freq2 = .true.
         do ic = 1, self%freq2_species_num
            self%freq2_species(getSpeciesIndex(tempListNames(ic))) = 1
         end do
      end if

      if ((self%pr_freq2) .and. (Tpr_freq2)) then
         self%pr_freq2 = .true.
      else
         self%pr_freq2 = .false.
      end if

      self%pr_freq2_period = Nint (pr_nc_freq2 * SECPDY)
      self%pr_at_time_freq2 = Nint (at_time_freq2 * SECPDY)
!
      ! Set the initial value of the list
      tempListNames(:) = ''

      ! Construct the list of names using the long string
      call constructListNames(tempListNames, freq3SpeciesNames)

      self%freq3_species_num = Count (tempListNames(:) /= '')
      if (self%freq3_species_num > 0) then
         Tpr_freq3 = .true.
         do ic = 1, self%freq3_species_num
            self%freq3_species(getSpeciesIndex(tempListNames(ic))) = 1
         end do
      end if

      if ((self%pr_freq3) .and. (Tpr_freq3)) then
         self%pr_freq3 = .true.
      else
         self%pr_freq3 = .false.
      end if

      self%pr_freq3_period = Nint (pr_nc_freq3 * SECPDY)
      self%pr_at_time_freq3 = Nint (at_time_freq3 * SECPDY)

!
      ! Set the initial value of the list
      tempListNames(:) = ''

      ! Construct the list of names using the long string
      call constructListNames(tempListNames, freq4SpeciesNames)

      self%freq4_species_num = Count (tempListNames(:) /= '')
      if (self%freq4_species_num > 0) then
         Tpr_freq4 = .true.
         do ic = 1, self%freq4_species_num
            self%freq4_species(getSpeciesIndex(tempListNames(ic))) = 1
         end do
      end if

      if ((self%pr_freq4) .and. (Tpr_freq4)) then
         self%pr_freq4 = .true.
      else
         self%pr_freq4 = .false.
      end if

      self%pr_freq4_period = Nint (pr_nc_freq4 * SECPDY)
      self%pr_at_time_freq4 = Nint (at_time_freq4 * SECPDY)

!     ----------------
!     Flux Diagnostics
!     ----------------

      ! Set the initial value of the list
      tempListNames(:) = ''

      ! Construct the list of names using the long string
      call constructListNames(tempListNames, fluxSpeciesNames)

      self%flux_species_num = Count (tempListNames(:) /= '')
      if (self%flux_species_num > 0) then
         allocate(self%flux_species(numSpecies))
         self%flux_species(:) = 0

         do ic = 1, self%flux_species_num
            self%flux_species(getSpeciesIndex(tempListNames(ic))) = 1
         end do
      end if

!      do ic = 1, numSpecies
!        if(flux_species(ic).eq.1.and.advec_flag(ic).eq.0) then
!          print *,'*******************************************'
!          print *,'NOT saving advective flux diag for species ',ic &
!     &     ,'because it is a non-advected species'
!          print *,'*******************************************'
!          flux_species(ic) = 0
!        endif
!      enddo

      if (self%pr_flux .and. (self%flux_species_num > 0)) then
         self%pr_flux = .true.
      else
         self%pr_flux = .false.
         self%pr_psf_flux = .false.
         self%pr_const_flux = .false.
      end if

!     ------------------------
!     Set misc. control stuff.
!     ------------------------

      if (pr_qj_period_days      == -5.0d0) &
     &       pr_qj_period_days      = pr_nc_period_days
      if (pr_qk_period_days      == -5.0d0) &
     &       pr_qk_period_days      = pr_nc_period_days
      if (pr_sad_period_days     == -5.0d0) &
     &        pr_sad_period_days     = pr_nc_period_days
      if (pr_tend_period_days    == -5.0d0) &
     &        pr_tend_period_days    = pr_nc_period_days
      if (pr_qqjk_period_days    == -5.0d0) &
     &        pr_qqjk_period_days    = pr_nc_period_days
      if (pr_cloud_period_days   == -5.0d0) &
     &        pr_cloud_period_days   = pr_nc_period_days
      if (pr_aerdust_period_days == -5.0d0) &
     &        pr_aerdust_period_days = pr_nc_period_days

      self%constOutputFrequency   = Nint (pr_nc_period_days      * SECPDY)
      self%qjOutputFrequency      = Nint (pr_qj_period_days      * SECPDY)
      self%qkOutputFrequency      = Nint (pr_qk_period_days      * SECPDY)
      self%qqjkOutputFrequency    = Nint (pr_qqjk_period_days    * SECPDY)
      self%sadOutputFrequency     = Nint (pr_sad_period_days     * SECPDY)
      self%fluxOutputFrequency    = Nint (pr_flux_period_days    * SECPDY)
      self%tendOutputFrequency    = Nint (pr_tend_period_days    * SECPDY)
      self%cloudOutputFrequency   = Nint (pr_cloud_period_days   * SECPDY)
      self%aerdustOutputFrequency = Nint (pr_aerdust_period_days * SECPDY)

      self%pr_freq1_period = Nint (pr_nc_freq1 * SECPDY)
      self%pr_freq2_period = Nint (pr_nc_freq2 * SECPDY)
      self%pr_freq3_period = Nint (pr_nc_freq3 * SECPDY)
      self%pr_freq4_period = Nint (pr_nc_freq4 * SECPDY)

      self%pr_at_time_freq1    = Nint (at_time_freq1 * SECPDY)
      self%pr_at_time_freq2    = Nint (at_time_freq2 * SECPDY)
      self%pr_at_time_freq3    = Nint (at_time_freq3 * SECPDY)
      self%pr_at_time_freq4    = Nint (at_time_freq4 * SECPDY)

      self%pr_overpass1_period = Nint (pr_overpass1_period_days * SECPDY)
      self%pr_overpass2_period = Nint (pr_overpass2_period_days * SECPDY)
      self%pr_overpass3_period = Nint (pr_overpass3_period_days * SECPDY)
      self%pr_overpass4_period = Nint (pr_overpass4_period_days * SECPDY)

      self%pr_rst_period  = Nint (pr_rst_period_days  * SECPDY)

      if (self%rd_restart) then
         if (self%restart_inrec == INT_DUM_VALUE) then
            call Ncop_Rd (ncid_rst, self%restart_infile_name)
            call Ncget_Unlim_Dimlen(ncid_rst, self%restart_inrec)
            call Nccl (ncid_rst)
         end if

      end if

      !###############
      ! Error Checking
      !###############

      if ((numSpecies < self%ascii_out_n) .and. self%pr_ascii) then
         err_msg = 'numSpecies/ascii_out_n problem in rc File.'
         call GmiPrintError (err_msg, .true., 2, numSpecies, &
     &           self%ascii_out_n,  0, 0.0d0, 0.0d0)
      end if

      !----------------------------------------------------
      ! Check that some things are evenly divisible by ndt.
      !----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, mdt, &
     &                label   = "mdt:", &
     &                default = 3600.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, tdt, &
     &                label   = "RUN_DT:", &
     &                default = 3600.0d0, rc=STATUS )
      VERIFY_(STATUS)

      ndt = Nint(tdt)

      if (Mod (Nint (self%constOutputFrequency), ndt) /= 0) then
        err_msg =  &
     &    'constOutputFrequency is not evenly div. by ndt in the rc File.'
        call GmiPrintError  &
     &    (err_msg, .true., 1, ndt, 0, 1, self%constOutputFrequency, 0.0d0)
      end if

      if (Mod (Nint (self%pr_rst_period), ndt) /= 0) then
        err_msg =  &
     &    'pr_rst_period is not evenly div. by ndt in the rc File.'
        call GmiPrintError  &
     &    (err_msg, .true., 1, ndt, 0, 1, self%pr_rst_period, 0.0d0)
      end if

!      if ((self%pr_rst_period > 0.0d0) .and.  &
!     &    ((self%pr_rst_period < mdt) .or.  &
!     &     (Mod (Nint (self%pr_rst_period), Nint (mdt)) /= 0))) then
!        err_msg = 'pr_rst_period/mdt problem in the rc File.'
!        call GmiPrintError  &
!     &    (err_msg, .true., 0, 0, 0, 2, self%pr_rst_period, mdt)
!      end if

      if (.not. self%pr_level_all) then
         if (self%k1r_gl .lt. k1 .or. (self%k2r_gl .gt. k2)) then
            err_msg = 'k1r_gl/k2r_gl problem.'
            call GmiPrintError  &
     &    (err_msg, .true., 2, self%k1r_gl, self%k2r_gl, 0, 0.0d0, 0.0d0)
         end if
      end if

      if (self%pr_const .and. (self%num_const_outrecs == 0)) then
        err_msg = 'pr_const/num_const_outrecs problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 1, self%num_const_outrecs, 0, 0, 0.0d0, 0.0d0)
      end if

      if (self%pr_flux .and. (self%num_const_outrecs == 0)) then
        err_msg = 'pr_flux/num_const_outrecs problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 1, self%num_const_outrecs, 0, 0, 0.0d0, 0.0d0)
      end if

      if (self%pr_tend .and. (self%num_const_outrecs == 0)) then
        err_msg = 'pr_tend/num_const_outrecs problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 1, self%num_const_outrecs, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((.not. self%pr_const) .and.  &
     &    (self%pr_psf        .or. self%pr_kel .or. self%pr_relHumidity .or. self%pr_mass .or.  &
     &     self%pr_metwater   .or. self%pr_dry_depos .or. self%pr_wet_depos .or.  &
     &     self%pr_surf_emiss .or. self%pr_grid_height .or. self%pr_sulf_src)) then
        err_msg = 'pr_const/pr_xxx problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((.not. self%pr_netcdf) .and. self%do_mean) then
        err_msg = 'pr_netcdf/do_mean problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((.not. self%do_mean) .and. (self%numSpecies_overpass1 > 0)) then
        err_msg = 'numSpecies_overpass1/do_mean problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      ! Overpass outputs

      if ((self%pr_overpass1) .and. (self%numSpecies_overpass1 == 0)) then
        err_msg = 'numSpecies_overpass1/pr_overpass1 problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%begTime_overpass1 .lt. 00.d0) .or.  &
     &    (self%endTime_overpass1 .gt. 24.d0) .or.  &
     &    (self%begTime_overpass1 .ge. self%endTime_overpass1) ) then
        err_msg = 'begTime_overpass1/endTime_overpass1  &
     &             problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((.not. self%do_mean) .and. (self%numSpecies_overpass2 > 0)) then
        err_msg = 'numSpecies_overpass2/do_mean problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%pr_overpass2) .and. (self%numSpecies_overpass2 == 0)) then
        err_msg = 'numSpecies_overpass2/pr_overpass2 problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%begTime_overpass2 .lt. 00.d0) .or.  &
     &    (self%endTime_overpass2 .gt. 24.d0) .or.  &
     &    (self%begTime_overpass2 .ge. self%endTime_overpass2) ) then
        err_msg = 'begTime_overpass2/endTime_overpass2  &
     &             problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((.not. self%do_mean) .and. (self%numSpecies_overpass3 > 0)) then
        err_msg = 'numSpecies_overpass3/do_mean problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%pr_overpass3) .and. (self%numSpecies_overpass3 == 0)) then
        err_msg = 'numSpecies_overpass3/pr_overpass3 problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%begTime_overpass3 .lt. 00.d0) .or.  &
     &    (self%endTime_overpass3 .gt. 24.d0) .or.  &
     &    (self%begTime_overpass3 .ge. self%endTime_overpass3) ) then
        err_msg = 'begTime_overpass3/endTime_overpass3  &
     &             problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((.not. self%do_mean) .and. (self%numSpecies_overpass4 > 0)) then
        err_msg = 'numSpecies_overpass4/do_mean problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%pr_overpass4) .and. (self%numSpecies_overpass4 == 0)) then
        err_msg = 'numSpecies_overpass4/pr_overpass4 problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%begTime_overpass4 .lt. 00.d0) .or.  &
     &    (self%endTime_overpass4 .gt. 24.d0) .or.  &
     &    (self%begTime_overpass4 .ge. self%endTime_overpass4) ) then
        err_msg = 'begTime_overpass4/endTime_overpass4  &
     &             problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      ! Freq outputs

      if ((self%iRange_freq1(1) > self%iRange_freq1(2)) .or. &
          (self%iRange_freq1(1) < i1_gl) .or. (self%iRange_freq1(2) > i2_gl)) then
        err_msg = 'lonRange_freq1  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%jRange_freq1(1) > self%jRange_freq1(2)) .or. &
          (self%jRange_freq1(1) < ju1_gl) .or. (self%jRange_freq1(2) > j2_gl)) then
        err_msg = 'latRange_freq1  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%iRange_freq2(1) > self%iRange_freq2(2)) .or. &
          (self%iRange_freq2(1) < i1_gl) .or. (self%iRange_freq2(2) > i2_gl)) then
        err_msg = 'lonRange_freq2  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%jRange_freq2(1) > self%jRange_freq2(2)) .or. &
          (self%jRange_freq2(1) < ju1_gl) .or. (self%jRange_freq2(2) > j2_gl)) then
        err_msg = 'latRange_freq2  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%iRange_freq3(1) > self%iRange_freq3(2)) .or. &
          (self%iRange_freq3(1) < i1_gl) .or. (self%iRange_freq3(2) > i2_gl)) then
        err_msg = 'lonRange_freq3  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%jRange_freq3(1) > self%jRange_freq3(2)) .or. &
          (self%jRange_freq3(1) < ju1_gl) .or. (self%jRange_freq3(2) > j2_gl)) then
        err_msg = 'latRange_freq3  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%iRange_freq4(1) > self%iRange_freq4(2)) .or. &
          (self%iRange_freq4(1) < i1_gl) .or. (self%iRange_freq4(2) > i2_gl)) then
        err_msg = 'lonRange_freq4  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%jRange_freq4(1) > self%jRange_freq4(2)) .or. &
          (self%jRange_freq4(1) < ju1_gl) .or. (self%jRange_freq4(2) > j2_gl)) then
        err_msg = 'latRange_freq4  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%k1_freq1 > self%k2_freq1) .or. (self%k1_freq1 < k1) .or. (self%k2_freq1 > k2)) then
        err_msg = 'k1_freq1/k2_freq1  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%k1_freq2 > self%k2_freq2) .or. (self%k1_freq2 < k1) .or. (self%k2_freq2 > k2)) then
        err_msg = 'k1_freq2/k2_freq2  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%k1_freq2 > self%k2_freq3) .or. (self%k1_freq3 < k1) .or. (self%k2_freq3 > k2)) then
        err_msg = 'k1_freq3/k2_freq3  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%k1_freq4 > self%k2_freq4) .or. (self%k1_freq4 < k1) .or. (self%k2_freq4 > k2)) then
        err_msg = 'k1_freq4/k2_freq4  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%pr_freq1) .and. (self%freq1_species_num == 0)) then
        err_msg = 'freq1_species_num/pr_freq1  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%do_mean_freq1) .and. (self%do_day1_freq1)) then
        err_msg = 'do_mean_freq1/do_day1_freq1  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%do_last_tstep_freq1) .and. (self%pr_at_time_freq1 > 0.0d0)) then
        err_msg = 'do_last_tstep_freq1/at_time_freq1  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if
!
      if ((self%pr_freq2) .and. (self%freq2_species_num == 0)) then
        err_msg = 'freq2_species_num/pr_freq2  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if
      if ((self%do_mean_freq2) .and. (self%do_day1_freq2)) then
        err_msg = 'do_mean_freq2/do_day1_freq2  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%do_last_tstep_freq2) .and. (self%pr_at_time_freq2 > 0.0d0)) then
        err_msg = 'do_last_tstep_freq2/at_time_freq2  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

!
      if ((self%pr_freq3) .and. (self%freq3_species_num == 0)) then
        err_msg = 'freq3_species_num/pr_freq3  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if
      if ((self%do_mean_freq3) .and. (self%do_day1_freq3)) then
        err_msg = 'do_mean_freq3/do_day1_freq3  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%do_last_tstep_freq3) .and. (self%pr_at_time_freq3 > 0.0d0)) then
        err_msg = 'do_last_tstep_freq3/at_time_freq3  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

!
      if ((self%pr_freq4) .and. (self%freq4_species_num == 0)) then
        err_msg = 'freq4_species_num/pr_freq4  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if
      if ((self%do_mean_freq4) .and. (self%do_day1_freq4)) then
        err_msg = 'do_mean_freq4/do_day1_freq4  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if ((self%do_last_tstep_freq4) .and. (self%pr_at_time_freq4 > 0.0d0)) then
        err_msg = 'do_last_tstep_freq4/at_time_freq4  problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

     if (self%pr_surf_emiss .and. (.not. self%pr_emiss_all) .and.  &
     &   self%num_emiss_outrecs == 0) then
        err_msg =  &
     &  'pr_surf_emiss/pr_emiss_all/pr_emiss_rec_flag problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if (self%pr_dry_depos .and. (.not. self%pr_drydep_all) .and.  &
     &   self%num_drydep_outrecs == 0) then
        err_msg =  &
     &  'pr_dry_depos/pr_drydep_all/pr_drydep_rec_flag problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if (self%pr_wet_depos .and. (.not. self%pr_wetdep_all) .and.  &
     &   self%num_wetdep_outrecs == 0) then
        err_msg =  &
     &  'pr_wet_depos/pr_wetdep_all/pr_wetdep_rec_flag problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if (self%pr_tend .and. (.not. self%pr_tend_all) .and.  &
     &   self%num_tend_outrecs == 0) then
        err_msg =  &
     &  'pr_tend/pr_tend_all/pr_tend_rec_flag problem.'
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if


      return

      end subroutine initializeDiagnostics
!EOC
!------------------------------------------------------------------------------
!##############################################################################
! Routines for obtaining derived type member variables
!##############################################################################
!------------------------------------------------------------------------------
      subroutine Get_problem_name(self, problem_name)
        implicit none
        character(len=*), intent(out) :: problem_name
        type(t_Diagnostics), intent(in) :: self
        problem_name = self%problem_name
        return
      end subroutine Get_problem_name
!------------------------------------------------------------------------------
      subroutine Get_pr_netcdf(self, pr_netcdf)
        implicit none
        logical, intent(out) :: pr_netcdf
        type(t_Diagnostics), intent(in) :: self
        pr_netcdf = self%pr_netcdf
        return
      end subroutine Get_pr_netcdf
!------------------------------------------------------------------------------
      subroutine Set_do_flux_reset(self, do_flux_reset)
        implicit none
        logical, intent(in) :: do_flux_reset
        type(t_Diagnostics), intent(inOut) :: self
        self%do_flux_reset = do_flux_reset
        return
      end subroutine Set_do_flux_reset
!------------------------------------------------------------------------------
      subroutine Get_do_flux_reset(self, do_flux_reset)
        implicit none
        logical, intent(out) :: do_flux_reset
        type(t_Diagnostics), intent(in) :: self
        do_flux_reset = self%do_flux_reset
        return
      end subroutine Get_do_flux_reset
!------------------------------------------------------------------------------
      subroutine Get_flux_name(self, flux_name)
        implicit none
        character(len=*), intent(out) :: flux_name
        type(t_Diagnostics), intent(in) :: self
        flux_name = self%flux_name
        return
      end subroutine Get_flux_name
!------------------------------------------------------------------------------
      subroutine Get_flux_var_name(self, flux_var_name)
        implicit none
        character(len=*), intent(out) :: flux_var_name
        type(t_Diagnostics), intent(in) :: self
        flux_var_name = self%flux_var_name
        return
      end subroutine Get_flux_var_name
!------------------------------------------------------------------------------
      subroutine Get_flux_species_num(self, flux_species_num)
        implicit none
        integer, intent(out) :: flux_species_num
        type(t_Diagnostics), intent(in) :: self
        flux_species_num = self%flux_species_num
        return
      end subroutine Get_flux_species_num
!------------------------------------------------------------------------------
      subroutine Set_flux_species(self, flux_species)
        implicit none
        integer, intent(in) :: flux_species(:)
        type(t_Diagnostics), intent(inOut) :: self
        self%flux_species(:) = flux_species(:)
        return
      end subroutine Set_flux_species
!------------------------------------------------------------------------------
      subroutine Get_flux_species(self, flux_species)
        implicit none
        integer, intent(out) :: flux_species(:)
        type(t_Diagnostics), intent(in) :: self
        flux_species(:) = self%flux_species(:)
        return
      end subroutine Get_flux_species
!------------------------------------------------------------------------------
      subroutine Get_do_mean_flux(self, do_mean_flux)
        implicit none
        logical, intent(out) :: do_mean_flux
        type(t_Diagnostics), intent(in) :: self
        do_mean_flux = self%do_mean_flux
        return
      end subroutine Get_do_mean_flux
!------------------------------------------------------------------------------
      subroutine Get_do_day1_flux(self, do_day1_flux)
        implicit none
        logical, intent(out) :: do_day1_flux
        type(t_Diagnostics), intent(in) :: self
        do_day1_flux = self%do_day1_flux
        return
      end subroutine Get_do_day1_flux
!------------------------------------------------------------------------------
      subroutine Get_pr_flux(self, pr_flux)
        implicit none
        logical, intent(out) :: pr_flux
        type(t_Diagnostics), intent(in) :: self
        pr_flux = self%pr_flux
        return
      end subroutine Get_pr_flux
!------------------------------------------------------------------------------
      subroutine Get_pr_psf_flux(self, pr_psf_flux)
        implicit none
        logical, intent(out) :: pr_psf_flux
        type(t_Diagnostics), intent(in) :: self
        pr_psf_flux = self%pr_psf_flux
        return
      end subroutine Get_pr_psf_flux
!------------------------------------------------------------------------------
      subroutine Get_pr_const_flux(self, pr_const_flux)
        implicit none
        logical, intent(out) :: pr_const_flux
        type(t_Diagnostics), intent(in) :: self
        pr_const_flux = self%pr_const_flux
        return
      end subroutine Get_pr_const_flux
!------------------------------------------------------------------------------
      subroutine Get_do_ftiming(self, do_ftiming)
        implicit none
        logical, intent(out) :: do_ftiming
        type(t_Diagnostics), intent(in) :: self
        do_ftiming = self%do_ftiming
        return
      end subroutine Get_do_ftiming
!------------------------------------------------------------------------------
      subroutine Get_pr_diag(self, pr_diag)
        implicit none
        logical, intent(out) :: pr_diag
        type(t_Diagnostics), intent(in) :: self
        pr_diag = self%pr_diag
        return
      end subroutine Get_pr_diag
!------------------------------------------------------------------------------
      subroutine Get_pr_smv2(self, pr_smv2)
        implicit none
        logical, intent(out) :: pr_smv2
        type(t_Diagnostics), intent(in) :: self
        pr_smv2 = self%pr_smv2
        return
      end subroutine Get_pr_smv2
!------------------------------------------------------------------------------
      subroutine Get_do_mean(self, do_mean)
        implicit none
        logical, intent(out) :: do_mean
        type(t_Diagnostics), intent(in) :: self
        do_mean = self%do_mean
        return
      end subroutine Get_do_mean
!------------------------------------------------------------------------------
      subroutine Get_do_aerocom(self, do_aerocom)
        implicit none
        logical, intent(out) :: do_aerocom
        type(t_Diagnostics), intent(in) :: self
        do_aerocom = self%do_aerocom
        return
      end subroutine Get_do_aerocom
!------------------------------------------------------------------------------
      subroutine Get_do_dust_emiss(self, do_dust_emiss)
        implicit none
        logical, intent(out) :: do_dust_emiss
        type(t_Diagnostics), intent(in) :: self
        do_dust_emiss = self%do_dust_emiss
        return
      end subroutine Get_do_dust_emiss
!------------------------------------------------------------------------------
      subroutine Get_hdr_var_name(self, hdr_var_name)
        implicit none
        character(len=*), intent(out) :: hdr_var_name
        type(t_Diagnostics), intent(in) :: self
        hdr_var_name = self%hdr_var_name
        return
      end subroutine Get_hdr_var_name
!------------------------------------------------------------------------------
      subroutine Get_hdf_dim_name(self, hdf_dim_name)
        implicit none
        character(len=*), intent(out) :: hdf_dim_name
        type(t_Diagnostics), intent(in) :: self
        hdf_dim_name = self%hdf_dim_name
        return
      end subroutine Get_hdf_dim_name
!------------------------------------------------------------------------------
      subroutine Get_lat_dim_name(self, lat_dim_name)
        implicit none
        character(len=*), intent(out) :: lat_dim_name
        type(t_Diagnostics), intent(in) :: self
        lat_dim_name = self%lat_dim_name
        return
      end subroutine Get_lat_dim_name
!------------------------------------------------------------------------------
      subroutine Get_lon_dim_name(self, lon_dim_name)
        implicit none
        character(len=*), intent(out) :: lon_dim_name
        type(t_Diagnostics), intent(in) :: self
        lon_dim_name = self%lon_dim_name
        return
      end subroutine Get_lon_dim_name
!------------------------------------------------------------------------------
      subroutine Get_prs_dim_name(self, prs_dim_name)
        implicit none
        character(len=*), intent(out) :: prs_dim_name
        type(t_Diagnostics), intent(in) :: self
        prs_dim_name = self%prs_dim_name
        return
      end subroutine Get_prs_dim_name
!------------------------------------------------------------------------------
      subroutine Get_spc_dim_name(self, spc_dim_name)
        implicit none
        character(len=*), intent(out) :: spc_dim_name
        type(t_Diagnostics), intent(in) :: self
        spc_dim_name = self%spc_dim_name
        return
      end subroutine Get_spc_dim_name
!------------------------------------------------------------------------------
      subroutine Get_rec_dim_name(self, rec_dim_name)
        implicit none
        character(len=*), intent(out) :: rec_dim_name
        type(t_Diagnostics), intent(in) :: self
        rec_dim_name = self%rec_dim_name
        return
      end subroutine Get_rec_dim_name
!------------------------------------------------------------------------------
      subroutine Get_tim_dim_name(self, tim_dim_name)
        implicit none
        character(len=*), intent(out) :: tim_dim_name
        type(t_Diagnostics), intent(in) :: self
        tim_dim_name = self%tim_dim_name
        return
      end subroutine Get_tim_dim_name
!------------------------------------------------------------------------------
      subroutine Get_sad_dim_name(self, sad_dim_name)
        implicit none
        character(len=*), intent(out) :: sad_dim_name
        type(t_Diagnostics), intent(in) :: self
        sad_dim_name = self%sad_dim_name
        return
      end subroutine Get_sad_dim_name
!------------------------------------------------------------------------------
      subroutine Get_sad_var_name(self, sad_var_name)
        implicit none
        character(len=*), intent(out) :: sad_var_name
        type(t_Diagnostics), intent(in) :: self
        sad_var_name = self%sad_var_name
        return
      end subroutine Get_sad_var_name
!------------------------------------------------------------------------------
      subroutine Get_qj_dim_name(self, qj_dim_name)
        implicit none
        character(len=*), intent(out) :: qj_dim_name
        type(t_Diagnostics), intent(in) :: self
        qj_dim_name = self%qj_dim_name
        return
      end subroutine Get_qj_dim_name
!------------------------------------------------------------------------------
      subroutine Get_qj_var_name(self, qj_var_name)
        implicit none
        character(len=*), intent(out) :: qj_var_name
        type(t_Diagnostics), intent(in) :: self
        qj_var_name = self%qj_var_name
        return
      end subroutine Get_qj_var_name
!------------------------------------------------------------------------------
      subroutine Get_qk_dim_name(self, qk_dim_name)
        implicit none
        character(len=*), intent(out) :: qk_dim_name
        type(t_Diagnostics), intent(in) :: self
        qk_dim_name = self%qk_dim_name
        return
      end subroutine Get_qk_dim_name
!------------------------------------------------------------------------------
      subroutine Get_qk_var_name(self, qk_var_name)
        implicit none
        character(len=*), intent(out) :: qk_var_name
        type(t_Diagnostics), intent(in) :: self
        qk_var_name = self%qk_var_name
        return
      end subroutine Get_qk_var_name
!------------------------------------------------------------------------------
      subroutine Get_qqj_dim_name(self, qqj_dim_name)
        implicit none
        character(len=*), intent(out) :: qqj_dim_name
        type(t_Diagnostics), intent(in) :: self
        qqj_dim_name = self%qqj_dim_name
        return
      end subroutine Get_qqj_dim_name
!------------------------------------------------------------------------------
      subroutine Get_qqj_var_name(self, qqj_var_name)
        implicit none
        character(len=*), intent(out) :: qqj_var_name
        type(t_Diagnostics), intent(in) :: self
        qqj_var_name = self%qqj_var_name
        return
      end subroutine Get_qqj_var_name
!------------------------------------------------------------------------------
      subroutine Get_qqk_dim_name(self, qqk_dim_name)
        implicit none
        character(len=*), intent(out) :: qqk_dim_name
        type(t_Diagnostics), intent(in) :: self
        qqk_dim_name = self%qqk_dim_name
        return
      end subroutine Get_qqk_dim_name
!------------------------------------------------------------------------------
      subroutine Get_qqk_var_name(self, qqk_var_name)
        implicit none
        character(len=*), intent(out) :: qqk_var_name
        type(t_Diagnostics), intent(in) :: self
        qqk_var_name = self%qqk_var_name
        return
      end subroutine Get_qqk_var_name
!------------------------------------------------------------------------------
      subroutine Get_k1r_gl(self, k1r_gl)
        implicit none
        integer, intent(out) :: k1r_gl
        type(t_Diagnostics), intent(in) :: self
        k1r_gl = self%k1r_gl
        return
      end subroutine Get_k1r_gl
!------------------------------------------------------------------------------
      subroutine Get_k2r_gl(self, k2r_gl)
        implicit none
        integer, intent(out) :: k2r_gl
        type(t_Diagnostics), intent(in) :: self
        k2r_gl = self%k2r_gl
        return
      end subroutine Get_k2r_gl
!------------------------------------------------------------------------------
      subroutine Get_cloudOutputFrequency(self, cloudOutputFrequency)
        implicit none
        real*8, intent(out) :: cloudOutputFrequency
        type(t_Diagnostics), intent(in) :: self
        cloudOutputFrequency = self%cloudOutputFrequency
        return
      end subroutine Get_cloudOutputFrequency
!------------------------------------------------------------------------------
      subroutine Get_fluxOutputFrequency(self, fluxOutputFrequency)
        implicit none
        real*8, intent(out) :: fluxOutputFrequency
        type(t_Diagnostics), intent(in) :: self
        fluxOutputFrequency = self%fluxOutputFrequency
        return
      end subroutine Get_fluxOutputFrequency
!------------------------------------------------------------------------------
      subroutine Get_sadOutputFrequency(self, sadOutputFrequency)
        implicit none
        real*8, intent(out) :: sadOutputFrequency
        type(t_Diagnostics), intent(in) :: self
        sadOutputFrequency = self%sadOutputFrequency
        return
      end subroutine Get_sadOutputFrequency
!------------------------------------------------------------------------------
      subroutine Get_tendOutputFrequency(self, tendOutputFrequency)
        implicit none
        real*8, intent(out) :: tendOutputFrequency
        type(t_Diagnostics), intent(in) :: self
        tendOutputFrequency = self%tendOutputFrequency
        return
      end subroutine Get_tendOutputFrequency
!------------------------------------------------------------------------------
      subroutine Get_qkOutputFrequency(self, qkOutputFrequency)
        implicit none
        real*8, intent(out) :: qkOutputFrequency
        type(t_Diagnostics), intent(in) :: self
        qkOutputFrequency = self%qkOutputFrequency
        return
      end subroutine Get_qkOutputFrequency
!------------------------------------------------------------------------------
      subroutine Get_qjOutputFrequency(self, qjOutputFrequency)
        implicit none
        real*8, intent(out) :: qjOutputFrequency
        type(t_Diagnostics), intent(in) :: self
        qjOutputFrequency = self%qjOutputFrequency
        return
      end subroutine Get_qjOutputFrequency
!------------------------------------------------------------------------------
      subroutine Get_qqjkOutputFrequency(self, qqjkOutputFrequency)
        implicit none
        real*8, intent(out) :: qqjkOutputFrequency
        type(t_Diagnostics), intent(in) :: self
        qqjkOutputFrequency = self%qqjkOutputFrequency
        return
      end subroutine Get_qqjkOutputFrequency
!------------------------------------------------------------------------------
      subroutine Get_aerdustOutputFrequency(self, aerdustOutputFrequency)
        implicit none
        real*8, intent(out) :: aerdustOutputFrequency
        type(t_Diagnostics), intent(in) :: self
        aerdustOutputFrequency = self%aerdustOutputFrequency
        return
      end subroutine Get_aerdustOutputFrequency
!------------------------------------------------------------------------------
      subroutine Get_constOutputFrequency(self, constOutputFrequency)
        implicit none
        real*8, intent(out) :: constOutputFrequency
        type(t_Diagnostics), intent(in) :: self
        constOutputFrequency = self%constOutputFrequency
        return
      end subroutine Get_constOutputFrequency
!------------------------------------------------------------------------------
      subroutine Get_pr_cloud(self, pr_cloud)
        implicit none
        logical, intent(out) :: pr_cloud
        type(t_Diagnostics), intent(in) :: self
        pr_cloud = self%pr_cloud
        return
      end subroutine Get_pr_cloud
!------------------------------------------------------------------------------
      subroutine Get_outmain_name(self, outmain_name)
        implicit none
        character (len=*), intent(out) :: outmain_name
        type(t_Diagnostics), intent(in) :: self
        outmain_name = self%outmain_name
        return
      end subroutine Get_outmain_name
!------------------------------------------------------------------------------
      subroutine Get_pr_tend(self, pr_tend)
        implicit none
        logical, intent(out) :: pr_tend
        type(t_Diagnostics), intent(in) :: self
        pr_tend = self%pr_tend
        return
      end subroutine Get_pr_tend
!------------------------------------------------------------------------------
      subroutine Get_num_tend_outrecs(self, num_tend_outrecs)
        implicit none
        integer, intent(out) :: num_tend_outrecs
        type(t_Diagnostics), intent(in) :: self
        num_tend_outrecs = self%num_tend_outrecs
        return
      end subroutine Get_num_tend_outrecs
!------------------------------------------------------------------------------
      subroutine Get_tend_outrec_map(self, tend_outrec_map)
        implicit none
        integer, intent(out) :: tend_outrec_map(:)
        type(t_Diagnostics), intent(in) :: self
        tend_outrec_map(:) = self%tend_outrec_map(:)
        return
      end subroutine Get_tend_outrec_map
!------------------------------------------------------------------------------
      subroutine Get_const_outrec_map(self, const_outrec_map)
        implicit none
        integer, intent(out) :: const_outrec_map(:)
        type(t_Diagnostics), intent(in) :: self
        const_outrec_map(:) = self%const_outrec_map(:)
        return
      end subroutine Get_const_outrec_map
!------------------------------------------------------------------------------
      subroutine Get_emiss_outrec_map(self, emiss_outrec_map)
        implicit none
        integer, intent(out) :: emiss_outrec_map(:)
        type(t_Diagnostics), intent(in) :: self
        emiss_outrec_map(:) = self%emiss_outrec_map(:)
        return
      end subroutine Get_emiss_outrec_map
!------------------------------------------------------------------------------
      subroutine Get_drydep_outrec_map(self, drydep_outrec_map)
        implicit none
        integer, intent(out) :: drydep_outrec_map(:)
        type(t_Diagnostics), intent(in) :: self
        drydep_outrec_map(:) = self%drydep_outrec_map(:)
        return
      end subroutine Get_drydep_outrec_map
!------------------------------------------------------------------------------
      subroutine Get_wetdep_outrec_map(self, wetdep_outrec_map)
        implicit none
        integer, intent(out) :: wetdep_outrec_map(:)
        type(t_Diagnostics), intent(in) :: self
        wetdep_outrec_map(:) = self%wetdep_outrec_map(:)
        return
      end subroutine Get_wetdep_outrec_map
!------------------------------------------------------------------------------
      subroutine Get_pr_overheadO3col(self, pr_overheadO3col)
        implicit none
        logical, intent(out) :: pr_overheadO3col
        type(t_Diagnostics), intent(in) :: self
        pr_overheadO3col = self%pr_overheadO3col
        return
      end subroutine Get_pr_overheadO3col
!------------------------------------------------------------------------------
      subroutine Get_pr_decay(self, pr_decay)
        implicit none
        logical, intent(out) :: pr_decay
        type(t_Diagnostics), intent(in) :: self
        pr_decay = self%pr_decay
        return
      end subroutine Get_pr_decay
!------------------------------------------------------------------------------
      subroutine Get_pr_tropopausePress(self, pr_tropopausePress)
        implicit none
        logical, intent(out) :: pr_tropopausePress
        type(t_Diagnostics), intent(in) :: self
        pr_tropopausePress = self%pr_tropopausePress
        return
      end subroutine Get_pr_tropopausePress
!------------------------------------------------------------------------------
      subroutine Get_pr_potentialVorticity(self, pr_potentialVorticity)
        implicit none
        logical, intent(out) :: pr_potentialVorticity
        type(t_Diagnostics), intent(in) :: self
        pr_potentialVorticity = self%pr_potentialVorticity
        return
      end subroutine Get_pr_potentialVorticity
!------------------------------------------------------------------------------
      subroutine Get_pr_metwater(self, pr_metwater)
        implicit none
        logical, intent(out) :: pr_metwater
        type(t_Diagnostics), intent(in) :: self
        pr_metwater = self%pr_metwater
        return
      end subroutine Get_pr_metwater
!------------------------------------------------------------------------------
      subroutine Get_pr_relHumidity(self, pr_relHumidity)
        implicit none
        logical, intent(out) :: pr_relHumidity
        type(t_Diagnostics), intent(in) :: self
        pr_relHumidity = self%pr_relHumidity
        return
      end subroutine Get_pr_relHumidity
!------------------------------------------------------------------------------
      subroutine Get_pr_const_surface(self, pr_const_surface)
        implicit none
        logical, intent(out) :: pr_const_surface
        type(t_Diagnostics), intent(in) :: self
        pr_const_surface = self%pr_const_surface
        return
      end subroutine Get_pr_const_surface
!------------------------------------------------------------------------------
      subroutine Get_pr_emiss_all(self, pr_emiss_all)
        implicit none
        logical, intent(out) :: pr_emiss_all
        type(t_Diagnostics), intent(in) :: self
        pr_emiss_all = self%pr_emiss_all
        return
      end subroutine Get_pr_emiss_all
!------------------------------------------------------------------------------
      subroutine Get_pr_const_column(self, pr_const_column)
        implicit none
        logical, intent(out) :: pr_const_column
        type(t_Diagnostics), intent(in) :: self
        pr_const_column = self%pr_const_column
        return
      end subroutine Get_pr_const_column
!------------------------------------------------------------------------------
      subroutine Get_pr_drydep_all(self, pr_drydep_all)
        implicit none
        logical, intent(out) :: pr_drydep_all
        type(t_Diagnostics), intent(in) :: self
        pr_drydep_all = self%pr_drydep_all
        return
      end subroutine Get_pr_drydep_all
!------------------------------------------------------------------------------
      subroutine Get_pr_wetdep_all(self, pr_wetdep_all)
        implicit none
        logical, intent(out) :: pr_wetdep_all
        type(t_Diagnostics), intent(in) :: self
        pr_wetdep_all = self%pr_wetdep_all
        return
      end subroutine Get_pr_wetdep_all
!------------------------------------------------------------------------------
      subroutine Get_pr_level_all(self, pr_level_all)
        implicit none
        logical, intent(out) :: pr_level_all
        type(t_Diagnostics), intent(in) :: self
        pr_level_all = self%pr_level_all
        return
      end subroutine Get_pr_level_all
!------------------------------------------------------------------------------
      subroutine Get_pr_const(self, pr_const)
        implicit none
        logical, intent(out) :: pr_const
        type(t_Diagnostics), intent(in) :: self
        pr_const = self%pr_const
        return
      end subroutine Get_pr_const
!------------------------------------------------------------------------------
      subroutine Get_pr_emiss_3d(self, pr_emiss_3d)
        implicit none
        logical, intent(out) :: pr_emiss_3d
        type(t_Diagnostics), intent(in) :: self
        pr_emiss_3d = self%pr_emiss_3d
        return
      end subroutine Get_pr_emiss_3d
!------------------------------------------------------------------------------
      subroutine Get_pr_surf_emiss(self, pr_surf_emiss)
        implicit none
        logical, intent(out) :: pr_surf_emiss
        type(t_Diagnostics), intent(in) :: self
        pr_surf_emiss = self%pr_surf_emiss
        return
      end subroutine Get_pr_surf_emiss
!------------------------------------------------------------------------------
      subroutine Get_pr_sulf_src(self, pr_sulf_src)
        implicit none
        logical, intent(out) :: pr_sulf_src
        type(t_Diagnostics), intent(in) :: self
        pr_sulf_src = self%pr_sulf_src
        return
      end subroutine Get_pr_sulf_src
!------------------------------------------------------------------------------
      subroutine Get_pr_wet_depos(self, pr_wet_depos)
        implicit none
        logical, intent(out) :: pr_wet_depos
        type(t_Diagnostics), intent(in) :: self
        pr_wet_depos = self%pr_wet_depos
        return
      end subroutine Get_pr_wet_depos
!------------------------------------------------------------------------------
      subroutine Get_pr_dry_depos(self, pr_dry_depos)
        implicit none
        logical, intent(out) :: pr_dry_depos
        type(t_Diagnostics), intent(in) :: self
        pr_dry_depos = self%pr_dry_depos
        return
      end subroutine Get_pr_dry_depos
!------------------------------------------------------------------------------
      subroutine Get_pr_psf(self, pr_psf)
        implicit none
        logical, intent(out) :: pr_psf
        type(t_Diagnostics), intent(in) :: self
        pr_psf = self%pr_psf
        return
      end subroutine Get_pr_psf
!------------------------------------------------------------------------------
      subroutine Get_pr_mass(self, pr_mass)
        implicit none
        logical, intent(out) :: pr_mass
        type(t_Diagnostics), intent(in) :: self
        pr_mass = self%pr_mass
        return
      end subroutine Get_pr_mass
!------------------------------------------------------------------------------
      subroutine Get_pr_grid_height(self, pr_grid_height)
        implicit none
        logical, intent(out) :: pr_grid_height
        type(t_Diagnostics), intent(in) :: self
        pr_grid_height = self%pr_grid_height
        return
      end subroutine Get_pr_grid_height
!------------------------------------------------------------------------------
      subroutine Get_pr_kel(self, pr_kel)
        implicit none
        logical, intent(out) :: pr_kel
        type(t_Diagnostics), intent(in) :: self
        pr_kel = self%pr_kel
        return
      end subroutine Get_pr_kel
!------------------------------------------------------------------------------
      subroutine Get_num_emiss_outrecs(self, num_emiss_outrecs)
        implicit none
        integer, intent(out) :: num_emiss_outrecs
        type(t_Diagnostics), intent(in) :: self
        num_emiss_outrecs = self%num_emiss_outrecs
        return
      end subroutine Get_num_emiss_outrecs
!------------------------------------------------------------------------------
      subroutine Get_num_const_outrecs(self, num_const_outrecs)
        implicit none
        integer, intent(out) :: num_const_outrecs
        type(t_Diagnostics), intent(in) :: self
        num_const_outrecs = self%num_const_outrecs
        return
      end subroutine Get_num_const_outrecs
!------------------------------------------------------------------------------
      subroutine Get_num_drydep_outrecs(self, num_drydep_outrecs)
        implicit none
        integer, intent(out) :: num_drydep_outrecs
        type(t_Diagnostics), intent(in) :: self
        num_drydep_outrecs = self%num_drydep_outrecs
        return
      end subroutine Get_num_drydep_outrecs
!------------------------------------------------------------------------------
      subroutine Get_num_wetdep_outrecs(self, num_wetdep_outrecs)
        implicit none
        integer, intent(out) :: num_wetdep_outrecs
        type(t_Diagnostics), intent(in) :: self
        num_wetdep_outrecs = self%num_wetdep_outrecs
        return
      end subroutine Get_num_wetdep_outrecs
!------------------------------------------------------------------------------
      subroutine Get_pr_sad(self, pr_sad)
        implicit none
        logical, intent(out) :: pr_sad
        type(t_Diagnostics), intent(in) :: self
        pr_sad = self%pr_sad
        return
      end subroutine Get_pr_sad
!------------------------------------------------------------------------------
      subroutine Get_pr_qj(self, pr_qj)
        implicit none
        logical, intent(out) :: pr_qj
        type(t_Diagnostics), intent(in) :: self
        pr_qj = self%pr_qj
        return
      end subroutine Get_pr_qj
!------------------------------------------------------------------------------
      subroutine Get_pr_qk(self, pr_qk)
        implicit none
        logical, intent(out) :: pr_qk
        type(t_Diagnostics), intent(in) :: self
        pr_qk = self%pr_qk
        return
      end subroutine Get_pr_qk
!------------------------------------------------------------------------------
      subroutine Get_pr_qqjk(self, pr_qqjk)
        implicit none
        logical, intent(out) :: pr_qqjk
        type(t_Diagnostics), intent(in) :: self
        pr_qqjk = self%pr_qqjk
        return
      end subroutine Get_pr_qqjk
!------------------------------------------------------------------------------
      subroutine Get_do_qqjk_inchem(self, do_qqjk_inchem)
        implicit none
        logical, intent(out) :: do_qqjk_inchem
        type(t_Diagnostics), intent(in) :: self
        do_qqjk_inchem = self%do_qqjk_inchem
        return
      end subroutine Get_do_qqjk_inchem
!------------------------------------------------------------------------------
      subroutine Get_pr_qj_opt_depth(self, pr_qj_opt_depth)
        implicit none
        logical, intent(out) :: pr_qj_opt_depth
        type(t_Diagnostics), intent(in) :: self
        pr_qj_opt_depth = self%pr_qj_opt_depth
        return
      end subroutine Get_pr_qj_opt_depth
!------------------------------------------------------------------------------
      subroutine Get_pr_qj_o3_o1d(self, pr_qj_o3_o1d)
        implicit none
        logical, intent(out) :: pr_qj_o3_o1d
        type(t_Diagnostics), intent(in) :: self
        pr_qj_o3_o1d = self%pr_qj_o3_o1d
        return
      end subroutine Get_pr_qj_o3_o1d
!------------------------------------------------------------------------------
      subroutine Get_pr_AerDust(self, pr_AerDust)
        implicit none
        logical, intent(out) :: pr_AerDust
        type(t_Diagnostics), intent(in) :: self
        pr_AerDust = self%pr_AerDust
        return
      end subroutine Get_pr_AerDust
!------------------------------------------------------------------------------
      subroutine Get_AerDust_var_name(self, AerDust_var_name)
        implicit none
        character(len=*), intent(out) :: AerDust_var_name
        type(t_Diagnostics), intent(in) :: self
        AerDust_var_name = self%AerDust_var_name
        return
      end subroutine Get_AerDust_var_name
!------------------------------------------------------------------------------
      subroutine Get_outaerdust_name(self, outaerdust_name)
        implicit none
        character(len=*), intent(out) :: outaerdust_name
        type(t_Diagnostics), intent(in) :: self
        outaerdust_name = self%outaerdust_name
        return
      end subroutine Get_outaerdust_name
!------------------------------------------------------------------------------
      subroutine Get_do_qqjk_reset(self, do_qqjk_reset)
        implicit none
        logical, intent(out) :: do_qqjk_reset
        type(t_Diagnostics), intent(in) :: self
        do_qqjk_reset = self%do_qqjk_reset
        return
      end subroutine Get_do_qqjk_reset
!------------------------------------------------------------------------------
      subroutine Set_do_qqjk_reset(self, do_qqjk_reset)
        implicit none
        logical, intent(in) :: do_qqjk_reset
        type(t_Diagnostics), intent(inOut) :: self
        self%do_qqjk_reset = do_qqjk_reset
        return
      end subroutine Set_do_qqjk_reset
!------------------------------------------------------------------------------
      subroutine Get_pr_overpass1(self, pr_overpass1)
        implicit none
        logical, intent(out) :: pr_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_overpass1 = self%pr_overpass1
        return
      end subroutine Get_pr_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_kel_overpass1(self, pr_kel_overpass1)
        implicit none
        logical, intent(out) :: pr_kel_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_kel_overpass1 = self%pr_kel_overpass1
        return
      end subroutine Get_pr_kel_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_psf_overpass1(self, pr_psf_overpass1)
        implicit none
        logical, intent(out) :: pr_psf_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_psf_overpass1 = self%pr_psf_overpass1
        return
      end subroutine Get_pr_psf_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_qj_overpass1(self, pr_qj_overpass1)
        implicit none
        logical, intent(out) :: pr_qj_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_qj_overpass1 = self%pr_qj_overpass1
        return
      end subroutine Get_pr_qj_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_qqjk_overpass1(self, pr_qqjk_overpass1)
        implicit none
        logical, intent(out) :: pr_qqjk_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_qqjk_overpass1 = self%pr_qqjk_overpass1
        return
      end subroutine Get_pr_qqjk_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_const_overpass1(self, pr_const_overpass1)
        implicit none
        logical, intent(out) :: pr_const_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_const_overpass1 = self%pr_const_overpass1
        return
      end subroutine Get_pr_const_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_tropopausePress_overpass1(self, pr_tropopausePress_overpass1)
        implicit none
        logical, intent(out) :: pr_tropopausePress_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_tropopausePress_overpass1 = self%pr_tropopausePress_overpass1
        return
      end subroutine Get_pr_tropopausePress_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_cloudFraction_overpass1(self, pr_cloudFraction_overpass1)
        implicit none
        logical, intent(out) :: pr_cloudFraction_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_cloudFraction_overpass1 = self%pr_cloudFraction_overpass1
        return
      end subroutine Get_pr_cloudFraction_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_overheadO3col_overpass1(self, pr_overheadO3col_overpass1)
        implicit none
        logical, intent(out) :: pr_overheadO3col_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_overheadO3col_overpass1 = self%pr_overheadO3col_overpass1
        return
      end subroutine Get_pr_overheadO3col_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_cloudOptDepth_overpass1(self, pr_cloudOptDepth_overpass1)
        implicit none
        logical, intent(out) :: pr_cloudOptDepth_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_cloudOptDepth_overpass1 = self%pr_cloudOptDepth_overpass1
        return
      end subroutine Get_pr_cloudOptDepth_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_lightningNO_overpass1(self, pr_lightningNO_overpass1)
        implicit none
        logical, intent(out) :: pr_lightningNO_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_lightningNO_overpass1 = self%pr_lightningNO_overpass1
        return
      end subroutine Get_pr_lightningNO_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_gridBoxHeight_overpass1(self, pr_gridBoxHeight_overpass1)
        implicit none
        logical, intent(out) :: pr_gridBoxHeight_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_gridBoxHeight_overpass1 = self%pr_gridBoxHeight_overpass1
        return
      end subroutine Get_pr_gridBoxHeight_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_relHumidity_overpass1(self, pr_relHumidity_overpass1)
        implicit none
        logical, intent(out) :: pr_relHumidity_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_relHumidity_overpass1 = self%pr_relHumidity_overpass1
        return
      end subroutine Get_pr_relHumidity_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_totalMass_overpass1(self, pr_totalMass_overpass1)
        implicit none
        logical, intent(out) :: pr_totalMass_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_totalMass_overpass1 = self%pr_totalMass_overpass1
        return
      end subroutine Get_pr_totalMass_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_metwater_overpass1(self, pr_metwater_overpass1)
        implicit none
        logical, intent(out) :: pr_metwater_overpass1
        type(t_Diagnostics), intent(in) :: self
        pr_metwater_overpass1 = self%pr_metwater_overpass1
        return
      end subroutine Get_pr_metwater_overpass1
!------------------------------------------------------------------------------
      subroutine Get_numSpecies_overpass1(self, numSpecies_overpass1)
        implicit none
        integer, intent(out) :: numSpecies_overpass1
        type(t_Diagnostics), intent(in) :: self
        numSpecies_overpass1 = self%numSpecies_overpass1
        return
      end subroutine Get_numSpecies_overpass1
!------------------------------------------------------------------------------
      subroutine Get_species_overpass1(self, species_overpass1)
        implicit none
        integer, intent(out) :: species_overpass1(:)
        type(t_Diagnostics), intent(in) :: self
        species_overpass1 = self%species_overpass1
        return
      end subroutine Get_species_overpass1
!------------------------------------------------------------------------------
      subroutine Get_begTime_overpass1(self, begTime_overpass1)
        implicit none
        real*8, intent(out) :: begTime_overpass1
        type(t_Diagnostics), intent(in) :: self
        begTime_overpass1 = self%begTime_overpass1
        return
      end subroutine Get_begTime_overpass1
!------------------------------------------------------------------------------
      subroutine Get_endTime_overpass1(self, endTime_overpass1)
        implicit none
        real*8, intent(out) :: endTime_overpass1
        type(t_Diagnostics), intent(in) :: self
        endTime_overpass1 = self%endTime_overpass1
        return
      end subroutine Get_endTime_overpass1
!------------------------------------------------------------------------------
      subroutine Get_pr_overpass1_period(self, pr_overpass1_period)
        implicit none
        real*8, intent(out) :: pr_overpass1_period
        type(t_Diagnostics), intent(in) :: self
        pr_overpass1_period = self%pr_overpass1_period
        return
      end subroutine Get_pr_overpass1_period
!------------------------------------------------------------------------------

      subroutine Get_pr_overpass2(self, pr_overpass2)
        implicit none
        logical, intent(out) :: pr_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_overpass2 = self%pr_overpass2
        return
      end subroutine Get_pr_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_kel_overpass2(self, pr_kel_overpass2)
        implicit none
        logical, intent(out) :: pr_kel_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_kel_overpass2 = self%pr_kel_overpass2
        return
      end subroutine Get_pr_kel_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_psf_overpass2(self, pr_psf_overpass2)
        implicit none
        logical, intent(out) :: pr_psf_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_psf_overpass2 = self%pr_psf_overpass2
        return
      end subroutine Get_pr_psf_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_qj_overpass2(self, pr_qj_overpass2)
        implicit none
        logical, intent(out) :: pr_qj_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_qj_overpass2 = self%pr_qj_overpass2
        return
      end subroutine Get_pr_qj_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_qqjk_overpass2(self, pr_qqjk_overpass2)
        implicit none
        logical, intent(out) :: pr_qqjk_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_qqjk_overpass2 = self%pr_qqjk_overpass2
        return
      end subroutine Get_pr_qqjk_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_const_overpass2(self, pr_const_overpass2)
        implicit none
        logical, intent(out) :: pr_const_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_const_overpass2 = self%pr_const_overpass2
        return
      end subroutine Get_pr_const_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_tropopausePress_overpass2(self, pr_tropopausePress_overpass2)
        implicit none
        logical, intent(out) :: pr_tropopausePress_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_tropopausePress_overpass2 = self%pr_tropopausePress_overpass2
        return
      end subroutine Get_pr_tropopausePress_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_cloudFraction_overpass2(self, pr_cloudFraction_overpass2)
        implicit none
        logical, intent(out) :: pr_cloudFraction_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_cloudFraction_overpass2 = self%pr_cloudFraction_overpass2
        return
      end subroutine Get_pr_cloudFraction_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_overheadO3col_overpass2(self, pr_overheadO3col_overpass2)
        implicit none
        logical, intent(out) :: pr_overheadO3col_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_overheadO3col_overpass2 = self%pr_overheadO3col_overpass2
        return
      end subroutine Get_pr_overheadO3col_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_cloudOptDepth_overpass2(self, pr_cloudOptDepth_overpass2)
        implicit none
        logical, intent(out) :: pr_cloudOptDepth_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_cloudOptDepth_overpass2 = self%pr_cloudOptDepth_overpass2
        return
      end subroutine Get_pr_cloudOptDepth_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_lightningNO_overpass2(self, pr_lightningNO_overpass2)
        implicit none
        logical, intent(out) :: pr_lightningNO_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_lightningNO_overpass2 = self%pr_lightningNO_overpass2
        return
      end subroutine Get_pr_lightningNO_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_gridBoxHeight_overpass2(self, pr_gridBoxHeight_overpass2)
        implicit none
        logical, intent(out) :: pr_gridBoxHeight_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_gridBoxHeight_overpass2 = self%pr_gridBoxHeight_overpass2
        return
      end subroutine Get_pr_gridBoxHeight_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_relHumidity_overpass2(self, pr_relHumidity_overpass2)
        implicit none
        logical, intent(out) :: pr_relHumidity_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_relHumidity_overpass2 = self%pr_relHumidity_overpass2
        return
      end subroutine Get_pr_relHumidity_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_totalMass_overpass2(self, pr_totalMass_overpass2)
        implicit none
        logical, intent(out) :: pr_totalMass_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_totalMass_overpass2 = self%pr_totalMass_overpass2
        return
      end subroutine Get_pr_totalMass_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_metwater_overpass2(self, pr_metwater_overpass2)
        implicit none
        logical, intent(out) :: pr_metwater_overpass2
        type(t_Diagnostics), intent(in) :: self
        pr_metwater_overpass2 = self%pr_metwater_overpass2
        return
      end subroutine Get_pr_metwater_overpass2
!------------------------------------------------------------------------------
      subroutine Get_numSpecies_overpass2(self, numSpecies_overpass2)
        implicit none
        integer, intent(out) :: numSpecies_overpass2
        type(t_Diagnostics), intent(in) :: self
        numSpecies_overpass2 = self%numSpecies_overpass2
        return
      end subroutine Get_numSpecies_overpass2
!------------------------------------------------------------------------------
      subroutine Get_species_overpass2(self, species_overpass2)
        implicit none
        integer, intent(out) :: species_overpass2(:)
        type(t_Diagnostics), intent(in) :: self
        species_overpass2 = self%species_overpass2
        return
      end subroutine Get_species_overpass2
!------------------------------------------------------------------------------
      subroutine Get_begTime_overpass2(self, begTime_overpass2)
        implicit none
        real*8, intent(out) :: begTime_overpass2
        type(t_Diagnostics), intent(in) :: self
        begTime_overpass2 = self%begTime_overpass2
        return
      end subroutine Get_begTime_overpass2
!------------------------------------------------------------------------------
      subroutine Get_endTime_overpass2(self, endTime_overpass2)
        implicit none
        real*8, intent(out) :: endTime_overpass2
        type(t_Diagnostics), intent(in) :: self
        endTime_overpass2 = self%endTime_overpass2
        return
      end subroutine Get_endTime_overpass2
!------------------------------------------------------------------------------
      subroutine Get_pr_overpass2_period(self, pr_overpass2_period)
        implicit none
        real*8, intent(out) :: pr_overpass2_period
        type(t_Diagnostics), intent(in) :: self
        pr_overpass2_period = self%pr_overpass2_period
        return
      end subroutine Get_pr_overpass2_period
!------------------------------------------------------------------------------

      subroutine Get_pr_overpass3(self, pr_overpass3)
        implicit none
        logical, intent(out) :: pr_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_overpass3 = self%pr_overpass3
        return
      end subroutine Get_pr_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_kel_overpass3(self, pr_kel_overpass3)
        implicit none
        logical, intent(out) :: pr_kel_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_kel_overpass3 = self%pr_kel_overpass3
        return
      end subroutine Get_pr_kel_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_psf_overpass3(self, pr_psf_overpass3)
        implicit none
        logical, intent(out) :: pr_psf_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_psf_overpass3 = self%pr_psf_overpass3
        return
      end subroutine Get_pr_psf_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_qj_overpass3(self, pr_qj_overpass3)
        implicit none
        logical, intent(out) :: pr_qj_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_qj_overpass3 = self%pr_qj_overpass3
        return
      end subroutine Get_pr_qj_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_qqjk_overpass3(self, pr_qqjk_overpass3)
        implicit none
        logical, intent(out) :: pr_qqjk_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_qqjk_overpass3 = self%pr_qqjk_overpass3
        return
      end subroutine Get_pr_qqjk_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_const_overpass3(self, pr_const_overpass3)
        implicit none
        logical, intent(out) :: pr_const_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_const_overpass3 = self%pr_const_overpass3
        return
      end subroutine Get_pr_const_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_tropopausePress_overpass3(self, pr_tropopausePress_overpass3)
        implicit none
        logical, intent(out) :: pr_tropopausePress_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_tropopausePress_overpass3 = self%pr_tropopausePress_overpass3
        return
      end subroutine Get_pr_tropopausePress_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_cloudFraction_overpass3(self, pr_cloudFraction_overpass3)
        implicit none
        logical, intent(out) :: pr_cloudFraction_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_cloudFraction_overpass3 = self%pr_cloudFraction_overpass3
        return
      end subroutine Get_pr_cloudFraction_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_overheadO3col_overpass3(self, pr_overheadO3col_overpass3)
        implicit none
        logical, intent(out) :: pr_overheadO3col_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_overheadO3col_overpass3 = self%pr_overheadO3col_overpass3
        return
      end subroutine Get_pr_overheadO3col_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_cloudOptDepth_overpass3(self, pr_cloudOptDepth_overpass3)
        implicit none
        logical, intent(out) :: pr_cloudOptDepth_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_cloudOptDepth_overpass3 = self%pr_cloudOptDepth_overpass3
        return
      end subroutine Get_pr_cloudOptDepth_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_lightningNO_overpass3(self, pr_lightningNO_overpass3)
        implicit none
        logical, intent(out) :: pr_lightningNO_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_lightningNO_overpass3 = self%pr_lightningNO_overpass3
        return
      end subroutine Get_pr_lightningNO_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_gridBoxHeight_overpass3(self, pr_gridBoxHeight_overpass3)
        implicit none
        logical, intent(out) :: pr_gridBoxHeight_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_gridBoxHeight_overpass3 = self%pr_gridBoxHeight_overpass3
        return
      end subroutine Get_pr_gridBoxHeight_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_relHumidity_overpass3(self, pr_relHumidity_overpass3)
        implicit none
        logical, intent(out) :: pr_relHumidity_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_relHumidity_overpass3 = self%pr_relHumidity_overpass3
        return
      end subroutine Get_pr_relHumidity_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_totalMass_overpass3(self, pr_totalMass_overpass3)
        implicit none
        logical, intent(out) :: pr_totalMass_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_totalMass_overpass3 = self%pr_totalMass_overpass3
        return
      end subroutine Get_pr_totalMass_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_metwater_overpass3(self, pr_metwater_overpass3)
        implicit none
        logical, intent(out) :: pr_metwater_overpass3
        type(t_Diagnostics), intent(in) :: self
        pr_metwater_overpass3 = self%pr_metwater_overpass3
        return
      end subroutine Get_pr_metwater_overpass3
!------------------------------------------------------------------------------
      subroutine Get_numSpecies_overpass3(self, numSpecies_overpass3)
        implicit none
        integer, intent(out) :: numSpecies_overpass3
        type(t_Diagnostics), intent(in) :: self
        numSpecies_overpass3 = self%numSpecies_overpass3
        return
      end subroutine Get_numSpecies_overpass3
!------------------------------------------------------------------------------
      subroutine Get_species_overpass3(self, species_overpass3)
        implicit none
        integer, intent(out) :: species_overpass3(:)
        type(t_Diagnostics), intent(in) :: self
        species_overpass3 = self%species_overpass3
        return
      end subroutine Get_species_overpass3
!------------------------------------------------------------------------------
      subroutine Get_begTime_overpass3(self, begTime_overpass3)
        implicit none
        real*8, intent(out) :: begTime_overpass3
        type(t_Diagnostics), intent(in) :: self
        begTime_overpass3 = self%begTime_overpass3
        return
      end subroutine Get_begTime_overpass3
!------------------------------------------------------------------------------
      subroutine Get_endTime_overpass3(self, endTime_overpass3)
        implicit none
        real*8, intent(out) :: endTime_overpass3
        type(t_Diagnostics), intent(in) :: self
        endTime_overpass3 = self%endTime_overpass3
        return
      end subroutine Get_endTime_overpass3
!------------------------------------------------------------------------------
      subroutine Get_pr_overpass3_period(self, pr_overpass3_period)
        implicit none
        real*8, intent(out) :: pr_overpass3_period
        type(t_Diagnostics), intent(in) :: self
        pr_overpass3_period = self%pr_overpass3_period
        return
      end subroutine Get_pr_overpass3_period
!------------------------------------------------------------------------------
      subroutine Get_pr_overpass4(self, pr_overpass4)
        implicit none
        logical, intent(out) :: pr_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_overpass4 = self%pr_overpass4
        return
      end subroutine Get_pr_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_kel_overpass4(self, pr_kel_overpass4)
        implicit none
        logical, intent(out) :: pr_kel_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_kel_overpass4 = self%pr_kel_overpass4
        return
      end subroutine Get_pr_kel_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_psf_overpass4(self, pr_psf_overpass4)
        implicit none
        logical, intent(out) :: pr_psf_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_psf_overpass4 = self%pr_psf_overpass4
        return
      end subroutine Get_pr_psf_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_qj_overpass4(self, pr_qj_overpass4)
        implicit none
        logical, intent(out) :: pr_qj_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_qj_overpass4 = self%pr_qj_overpass4
        return
      end subroutine Get_pr_qj_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_qqjk_overpass4(self, pr_qqjk_overpass4)
        implicit none
        logical, intent(out) :: pr_qqjk_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_qqjk_overpass4 = self%pr_qqjk_overpass4
        return
      end subroutine Get_pr_qqjk_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_const_overpass4(self, pr_const_overpass4)
        implicit none
        logical, intent(out) :: pr_const_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_const_overpass4 = self%pr_const_overpass4
        return
      end subroutine Get_pr_const_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_tropopausePress_overpass4(self, pr_tropopausePress_overpass4)
        implicit none
        logical, intent(out) :: pr_tropopausePress_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_tropopausePress_overpass4 = self%pr_tropopausePress_overpass4
        return
      end subroutine Get_pr_tropopausePress_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_cloudFraction_overpass4(self, pr_cloudFraction_overpass4)
        implicit none
        logical, intent(out) :: pr_cloudFraction_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_cloudFraction_overpass4 = self%pr_cloudFraction_overpass4
        return
      end subroutine Get_pr_cloudFraction_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_overheadO3col_overpass4(self, pr_overheadO3col_overpass4)
        implicit none
        logical, intent(out) :: pr_overheadO3col_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_overheadO3col_overpass4 = self%pr_overheadO3col_overpass4
        return
      end subroutine Get_pr_overheadO3col_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_cloudOptDepth_overpass4(self, pr_cloudOptDepth_overpass4)
        implicit none
        logical, intent(out) :: pr_cloudOptDepth_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_cloudOptDepth_overpass4 = self%pr_cloudOptDepth_overpass4
        return
      end subroutine Get_pr_cloudOptDepth_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_lightningNO_overpass4(self, pr_lightningNO_overpass4)
        implicit none
        logical, intent(out) :: pr_lightningNO_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_lightningNO_overpass4 = self%pr_lightningNO_overpass4
        return
      end subroutine Get_pr_lightningNO_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_gridBoxHeight_overpass4(self, pr_gridBoxHeight_overpass4)
        implicit none
        logical, intent(out) :: pr_gridBoxHeight_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_gridBoxHeight_overpass4 = self%pr_gridBoxHeight_overpass4
        return
      end subroutine Get_pr_gridBoxHeight_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_relHumidity_overpass4(self, pr_relHumidity_overpass4)
        implicit none
        logical, intent(out) :: pr_relHumidity_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_relHumidity_overpass4 = self%pr_relHumidity_overpass4
        return
      end subroutine Get_pr_relHumidity_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_totalMass_overpass4(self, pr_totalMass_overpass4)
        implicit none
        logical, intent(out) :: pr_totalMass_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_totalMass_overpass4 = self%pr_totalMass_overpass4
        return
      end subroutine Get_pr_totalMass_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_metwater_overpass4(self, pr_metwater_overpass4)
        implicit none
        logical, intent(out) :: pr_metwater_overpass4
        type(t_Diagnostics), intent(in) :: self
        pr_metwater_overpass4 = self%pr_metwater_overpass4
        return
      end subroutine Get_pr_metwater_overpass4
!------------------------------------------------------------------------------
      subroutine Get_numSpecies_overpass4(self, numSpecies_overpass4)
        implicit none
        integer, intent(out) :: numSpecies_overpass4
        type(t_Diagnostics), intent(in) :: self
        numSpecies_overpass4 = self%numSpecies_overpass4
        return
      end subroutine Get_numSpecies_overpass4
!------------------------------------------------------------------------------
      subroutine Get_species_overpass4(self, species_overpass4)
        implicit none
        integer, intent(out) :: species_overpass4(:)
        type(t_Diagnostics), intent(in) :: self
        species_overpass4 = self%species_overpass4
        return
      end subroutine Get_species_overpass4
!------------------------------------------------------------------------------
      subroutine Get_begTime_overpass4(self, begTime_overpass4)
        implicit none
        real*8, intent(out) :: begTime_overpass4
        type(t_Diagnostics), intent(in) :: self
        begTime_overpass4 = self%begTime_overpass4
        return
      end subroutine Get_begTime_overpass4
!------------------------------------------------------------------------------
      subroutine Get_endTime_overpass4(self, endTime_overpass4)
        implicit none
        real*8, intent(out) :: endTime_overpass4
        type(t_Diagnostics), intent(in) :: self
        endTime_overpass4 = self%endTime_overpass4
        return
      end subroutine Get_endTime_overpass4
!------------------------------------------------------------------------------
      subroutine Get_pr_overpass4_period(self, pr_overpass4_period)
        implicit none
        real*8, intent(out) :: pr_overpass4_period
        type(t_Diagnostics), intent(in) :: self
        pr_overpass4_period = self%pr_overpass4_period
        return
      end subroutine Get_pr_overpass4_period
!------------------------------------------------------------------------------
      subroutine Get_do_overwrt_rst(self, do_overwrt_rst)
        implicit none
        logical, intent(out) :: do_overwrt_rst
        type(t_Diagnostics), intent(in) :: self
        do_overwrt_rst = self%do_overwrt_rst
        return
      end subroutine Get_do_overwrt_rst
!------------------------------------------------------------------------------
      subroutine Get_pr_restart(self, pr_restart)
        implicit none
        logical, intent(out) :: pr_restart
        type(t_Diagnostics), intent(in) :: self
        pr_restart = self%pr_restart
        return
      end subroutine Get_pr_restart
!------------------------------------------------------------------------------
      subroutine Get_rd_restart(self, rd_restart)
        implicit none
        logical, intent(out) :: rd_restart
        type(t_Diagnostics), intent(in) :: self
        rd_restart = self%rd_restart
        return
      end subroutine Get_rd_restart
!------------------------------------------------------------------------------
      subroutine Get_restart_inrec(self, restart_inrec)
        implicit none
        integer, intent(out) :: restart_inrec
        type(t_Diagnostics), intent(in) :: self
        restart_inrec = self%restart_inrec
        return
      end subroutine Get_restart_inrec
!------------------------------------------------------------------------------
      subroutine Get_pr_rst_period(self, pr_rst_period)
        implicit none
        real*8 , intent(out) :: pr_rst_period
        type(t_Diagnostics), intent(in) :: self
        pr_rst_period = self%pr_rst_period
        return
      end subroutine Get_pr_rst_period
!------------------------------------------------------------------------------
      subroutine Get_restart_infile_name(self, restart_infile_name)
        implicit none
        character(len=*), intent(out) :: restart_infile_name
        type(t_Diagnostics), intent(in) :: self
        restart_infile_name = self%restart_infile_name
        return
      end subroutine Get_restart_infile_name
!------------------------------------------------------------------------------
      subroutine Get_k1_freq1(self, k1_freq1)
        implicit none
        integer, intent(out) :: k1_freq1
        type(t_Diagnostics), intent(in) :: self
        k1_freq1 = self%k1_freq1
        return
      end subroutine Get_k1_freq1
!------------------------------------------------------------------------------
      subroutine Get_k2_freq1(self, k2_freq1)
        implicit none
        integer, intent(out) :: k2_freq1
        type(t_Diagnostics), intent(in) :: self
        k2_freq1 = self%k2_freq1
        return
      end subroutine Get_k2_freq1
!------------------------------------------------------------------------------
      subroutine Get_do_day1_freq1(self, do_day1_freq1)
        implicit none
        logical, intent(out) :: do_day1_freq1
        type(t_Diagnostics), intent(in) :: self
        do_day1_freq1 = self%do_day1_freq1
        return
      end subroutine Get_do_day1_freq1
!------------------------------------------------------------------------------
      subroutine Get_pr_potentialVorticity_freq1(self, pr_potentialVorticity_freq1)
        implicit none
        logical, intent(out) :: pr_potentialVorticity_freq1
        type(t_Diagnostics), intent(in) :: self
        pr_potentialVorticity_freq1 = self%pr_potentialVorticity_freq1
        return
      end subroutine Get_pr_potentialVorticity_freq1
!------------------------------------------------------------------------------
      subroutine Get_pr_freq1_period(self, pr_freq1_period)
        implicit none
        real*8, intent(out) :: pr_freq1_period
        type(t_Diagnostics), intent(in) :: self
        pr_freq1_period = self%pr_freq1_period
        return
      end subroutine Get_pr_freq1_period
!------------------------------------------------------------------------------
      subroutine Get_pr_at_time_freq1(self, pr_at_time_freq1)
        implicit none
        real*8, intent(out) :: pr_at_time_freq1
        type(t_Diagnostics), intent(in) :: self
        pr_at_time_freq1 = self%pr_at_time_freq1
        return
      end subroutine Get_pr_at_time_freq1
!------------------------------------------------------------------------------
      subroutine Get_freq1_description(self, freq1_description)
        implicit none
        character(len=*), intent(out) :: freq1_description
        type(t_Diagnostics), intent(in) :: self
        freq1_description(:) = self%freq1_description(:)
        return
      end subroutine Get_freq1_description
!------------------------------------------------------------------------------
      subroutine Get_freq1_name(self, freq1_name)
        implicit none
        character(len=*), intent(out) :: freq1_name
        type(t_Diagnostics), intent(in) :: self
        freq1_name(:) = self%freq1_name(:)
        return
      end subroutine Get_freq1_name
!------------------------------------------------------------------------------
      subroutine Get_iRange_freq1(self, iRange_freq1)
        implicit none
        integer, intent(out) :: iRange_freq1(:)
        type(t_Diagnostics), intent(in) :: self
        iRange_freq1(:) = self%iRange_freq1(:)
        return
      end subroutine Get_iRange_freq1
!------------------------------------------------------------------------------
      subroutine Get_jRange_freq1(self, jRange_freq1)
        implicit none
        integer, intent(out) :: jRange_freq1(:)
        type(t_Diagnostics), intent(in) :: self
        jRange_freq1(:) = self%jRange_freq1(:)
        return
      end subroutine Get_jRange_freq1
!------------------------------------------------------------------------------
      subroutine Get_freq1_species_num(self, freq1_species_num)
        implicit none
        integer, intent(out) :: freq1_species_num
        type(t_Diagnostics), intent(in) :: self
        freq1_species_num = self%freq1_species_num
        return
      end subroutine Get_freq1_species_num
!------------------------------------------------------------------------------
      subroutine Get_freq1_species(self, freq1_species)
        implicit none
        integer, intent(out) :: freq1_species(:)
        type(t_Diagnostics), intent(in) :: self
        freq1_species(:) = self%freq1_species(:)
        return
      end subroutine Get_freq1_species
!------------------------------------------------------------------------------
      subroutine Get_pr_freq1(self, pr_freq1)
        implicit none
        logical, intent(out) :: pr_freq1
        type(t_Diagnostics), intent(in) :: self
        pr_freq1 = self%pr_freq1
        return
      end subroutine Get_pr_freq1
!------------------------------------------------------------------------------
      subroutine Get_do_mean_freq1(self, do_mean_freq1)
        implicit none
        logical, intent(out) :: do_mean_freq1
        type(t_Diagnostics), intent(in) :: self
        do_mean_freq1 = self%do_mean_freq1
        return
      end subroutine Get_do_mean_freq1
!------------------------------------------------------------------------------
      subroutine Get_pr_mass_freq1(self, pr_mass_freq1)
        implicit none
        logical, intent(out) :: pr_mass_freq1
        type(t_Diagnostics), intent(in) :: self
        pr_mass_freq1 = self%pr_mass_freq1
        return
      end subroutine Get_pr_mass_freq1
!------------------------------------------------------------------------------
      subroutine Get_pr_psf_freq1(self, pr_psf_freq1)
        implicit none
        logical, intent(out) :: pr_psf_freq1
        type(t_Diagnostics), intent(in) :: self
        pr_psf_freq1 = self%pr_psf_freq1
        return
      end subroutine Get_pr_psf_freq1
!------------------------------------------------------------------------------
      subroutine Get_pr_kel_freq1(self, pr_kel_freq1)
        implicit none
        logical, intent(out) :: pr_kel_freq1
        type(t_Diagnostics), intent(in) :: self
        pr_kel_freq1 = self%pr_kel_freq1
        return
      end subroutine Get_pr_kel_freq1
!------------------------------------------------------------------------------
      subroutine Get_pr_const_freq1(self, pr_const_freq1)
        implicit none
        logical, intent(out) :: pr_const_freq1
        type(t_Diagnostics), intent(in) :: self
        pr_const_freq1 = self%pr_const_freq1
        return
      end subroutine Get_pr_const_freq1
!------------------------------------------------------------------------------
      subroutine Get_pr_rel_hum_freq1(self, pr_rel_hum_freq1)
        implicit none
        logical, intent(out) :: pr_rel_hum_freq1
        type(t_Diagnostics), intent(in) :: self
        pr_rel_hum_freq1 = self%pr_rel_hum_freq1
        return
      end subroutine Get_pr_rel_hum_freq1
!------------------------------------------------------------------------------
      subroutine Get_pr_metwater_freq1(self, pr_metwater_freq1)
        implicit none
        logical, intent(out) :: pr_metwater_freq1
        type(t_Diagnostics), intent(in) :: self
        pr_metwater_freq1 = self%pr_metwater_freq1
        return
      end subroutine Get_pr_metwater_freq1
!------------------------------------------------------------------------------
      subroutine Get_do_last_tstep_freq1(self, do_last_tstep_freq1)
        implicit none
        logical, intent(out) :: do_last_tstep_freq1
        type(t_Diagnostics), intent(in) :: self
        do_last_tstep_freq1 = self%do_last_tstep_freq1
        return
      end subroutine Get_do_last_tstep_freq1
!------------------------------------------------------------------------------
      subroutine Get_pr_grid_height_freq1(self, pr_grid_height_freq1)
        implicit none
        logical, intent(out) :: pr_grid_height_freq1
        type(t_Diagnostics), intent(in) :: self
        pr_grid_height_freq1 = self%pr_grid_height_freq1
        return
      end subroutine Get_pr_grid_height_freq1
!------------------------------------------------------------------------------
      subroutine Get_pr_const_column_freq1(self, pr_const_column_freq1)
        implicit none
        logical, intent(out) :: pr_const_column_freq1
        type(t_Diagnostics), intent(in) :: self
        pr_const_column_freq1 = self%pr_const_column_freq1
        return
      end subroutine Get_pr_const_column_freq1
!------------------------------------------------------------------------------
      subroutine Get_pr_const_surface_freq1(self, pr_const_surface_freq1)
        implicit none
        logical, intent(out) :: pr_const_surface_freq1
        type(t_Diagnostics), intent(in) :: self
        pr_const_surface_freq1 = self%pr_const_surface_freq1
        return
      end subroutine Get_pr_const_surface_freq1
!------------------------------------------------------------------------------
      subroutine Get_pr_overheadO3col_freq1(self, pr_overheadO3col_freq1)
        implicit none
        logical, intent(out) :: pr_overheadO3col_freq1
        type(t_Diagnostics), intent(in) :: self
        pr_overheadO3col_freq1 = self%pr_overheadO3col_freq1
        return
      end subroutine Get_pr_overheadO3col_freq1
!------------------------------------------------------------------------------
      subroutine Get_pr_tropopausePress_freq1(self, pr_tropopausePress_freq1)
        implicit none
        logical, intent(out) :: pr_tropopausePress_freq1
        type(t_Diagnostics), intent(in) :: self
        pr_tropopausePress_freq1 = self%pr_tropopausePress_freq1
        return
      end subroutine Get_pr_tropopausePress_freq1
!------------------------------------------------------------------------------
      subroutine Get_k1_freq2(self, k1_freq2)
        implicit none
        integer, intent(out) :: k1_freq2
        type(t_Diagnostics), intent(in) :: self
        k1_freq2 = self%k1_freq2
        return
      end subroutine Get_k1_freq2
!------------------------------------------------------------------------------
      subroutine Get_k2_freq2(self, k2_freq2)
        implicit none
        integer, intent(out) :: k2_freq2
        type(t_Diagnostics), intent(in) :: self
        k2_freq2 = self%k2_freq2
        return
      end subroutine Get_k2_freq2
!------------------------------------------------------------------------------
      subroutine Get_do_day1_freq2(self, do_day1_freq2)
        implicit none
        logical, intent(out) :: do_day1_freq2
        type(t_Diagnostics), intent(in) :: self
        do_day1_freq2 = self%do_day1_freq2
        return
      end subroutine Get_do_day1_freq2
!------------------------------------------------------------------------------
      subroutine Get_pr_potentialVorticity_freq2(self, pr_potentialVorticity_freq2)
        implicit none
        logical, intent(out) :: pr_potentialVorticity_freq2
        type(t_Diagnostics), intent(in) :: self
        pr_potentialVorticity_freq2 = self%pr_potentialVorticity_freq2
        return
      end subroutine Get_pr_potentialVorticity_freq2
!------------------------------------------------------------------------------
      subroutine Get_pr_freq2_period(self, pr_freq2_period)
        implicit none
        real*8, intent(out) :: pr_freq2_period
        type(t_Diagnostics), intent(in) :: self
        pr_freq2_period = self%pr_freq2_period
        return
      end subroutine Get_pr_freq2_period
!------------------------------------------------------------------------------
      subroutine Get_pr_at_time_freq2(self, pr_at_time_freq2)
        implicit none
        real*8, intent(out) :: pr_at_time_freq2
        type(t_Diagnostics), intent(in) :: self
        pr_at_time_freq2 = self%pr_at_time_freq2
        return
      end subroutine Get_pr_at_time_freq2
!------------------------------------------------------------------------------
      subroutine Get_freq2_description(self, freq2_description)
        implicit none
        character(len=*), intent(out) :: freq2_description
        type(t_Diagnostics), intent(in) :: self
        freq2_description(:) = self%freq2_description(:)
        return
      end subroutine Get_freq2_description
!------------------------------------------------------------------------------
      subroutine Get_freq2_name(self, freq2_name)
        implicit none
        character(len=*), intent(out) :: freq2_name
        type(t_Diagnostics), intent(in) :: self
        freq2_name(:) = self%freq2_name(:)
        return
      end subroutine Get_freq2_name
!------------------------------------------------------------------------------
      subroutine Get_iRange_freq2(self, iRange_freq2)
        implicit none
        integer, intent(out) :: iRange_freq2(:)
        type(t_Diagnostics), intent(in) :: self
        iRange_freq2(:) = self%iRange_freq2(:)
        return
      end subroutine Get_iRange_freq2
!------------------------------------------------------------------------------
      subroutine Get_jRange_freq2(self, jRange_freq2)
        implicit none
        integer, intent(out) :: jRange_freq2(:)
        type(t_Diagnostics), intent(in) :: self
        jRange_freq2(:) = self%jRange_freq2(:)
        return
      end subroutine Get_jRange_freq2
!------------------------------------------------------------------------------
      subroutine Get_freq2_species_num(self, freq2_species_num)
        implicit none
        integer, intent(out) :: freq2_species_num
        type(t_Diagnostics), intent(in) :: self
        freq2_species_num = self%freq2_species_num
        return
      end subroutine Get_freq2_species_num
!------------------------------------------------------------------------------
      subroutine Get_freq2_species(self, freq2_species)
        implicit none
        integer, intent(out) :: freq2_species(:)
        type(t_Diagnostics), intent(in) :: self
        freq2_species(:) = self%freq2_species(:)
        return
      end subroutine Get_freq2_species
!------------------------------------------------------------------------------
      subroutine Get_pr_freq2(self, pr_freq2)
        implicit none
        logical, intent(out) :: pr_freq2
        type(t_Diagnostics), intent(in) :: self
        pr_freq2 = self%pr_freq2
        return
      end subroutine Get_pr_freq2
!------------------------------------------------------------------------------
      subroutine Get_do_mean_freq2(self, do_mean_freq2)
        implicit none
        logical, intent(out) :: do_mean_freq2
        type(t_Diagnostics), intent(in) :: self
        do_mean_freq2 = self%do_mean_freq2
        return
      end subroutine Get_do_mean_freq2
!------------------------------------------------------------------------------
      subroutine Get_pr_mass_freq2(self, pr_mass_freq2)
        implicit none
        logical, intent(out) :: pr_mass_freq2
        type(t_Diagnostics), intent(in) :: self
        pr_mass_freq2 = self%pr_mass_freq2
        return
      end subroutine Get_pr_mass_freq2
!------------------------------------------------------------------------------
      subroutine Get_pr_psf_freq2(self, pr_psf_freq2)
        implicit none
        logical, intent(out) :: pr_psf_freq2
        type(t_Diagnostics), intent(in) :: self
        pr_psf_freq2 = self%pr_psf_freq2
        return
      end subroutine Get_pr_psf_freq2
!------------------------------------------------------------------------------
      subroutine Get_pr_kel_freq2(self, pr_kel_freq2)
        implicit none
        logical, intent(out) :: pr_kel_freq2
        type(t_Diagnostics), intent(in) :: self
        pr_kel_freq2 = self%pr_kel_freq2
        return
      end subroutine Get_pr_kel_freq2
!------------------------------------------------------------------------------
      subroutine Get_pr_const_freq2(self, pr_const_freq2)
        implicit none
        logical, intent(out) :: pr_const_freq2
        type(t_Diagnostics), intent(in) :: self
        pr_const_freq2 = self%pr_const_freq2
        return
      end subroutine Get_pr_const_freq2
!------------------------------------------------------------------------------
      subroutine Get_pr_rel_hum_freq2(self, pr_rel_hum_freq2)
        implicit none
        logical, intent(out) :: pr_rel_hum_freq2
        type(t_Diagnostics), intent(in) :: self
        pr_rel_hum_freq2 = self%pr_rel_hum_freq2
        return
      end subroutine Get_pr_rel_hum_freq2
!------------------------------------------------------------------------------
      subroutine Get_pr_metwater_freq2(self, pr_metwater_freq2)
        implicit none
        logical, intent(out) :: pr_metwater_freq2
        type(t_Diagnostics), intent(in) :: self
        pr_metwater_freq2 = self%pr_metwater_freq2
        return
      end subroutine Get_pr_metwater_freq2
!------------------------------------------------------------------------------
      subroutine Get_do_last_tstep_freq2(self, do_last_tstep_freq2)
        implicit none
        logical, intent(out) :: do_last_tstep_freq2
        type(t_Diagnostics), intent(in) :: self
        do_last_tstep_freq2 = self%do_last_tstep_freq2
        return
      end subroutine Get_do_last_tstep_freq2
!------------------------------------------------------------------------------
      subroutine Get_pr_grid_height_freq2(self, pr_grid_height_freq2)
        implicit none
        logical, intent(out) :: pr_grid_height_freq2
        type(t_Diagnostics), intent(in) :: self
        pr_grid_height_freq2 = self%pr_grid_height_freq2
        return
      end subroutine Get_pr_grid_height_freq2
!------------------------------------------------------------------------------
      subroutine Get_pr_const_column_freq2(self, pr_const_column_freq2)
        implicit none
        logical, intent(out) :: pr_const_column_freq2
        type(t_Diagnostics), intent(in) :: self
        pr_const_column_freq2 = self%pr_const_column_freq2
        return
      end subroutine Get_pr_const_column_freq2
!------------------------------------------------------------------------------
      subroutine Get_pr_const_surface_freq2(self, pr_const_surface_freq2)
        implicit none
        logical, intent(out) :: pr_const_surface_freq2
        type(t_Diagnostics), intent(in) :: self
        pr_const_surface_freq2 = self%pr_const_surface_freq2
        return
      end subroutine Get_pr_const_surface_freq2
!------------------------------------------------------------------------------
      subroutine Get_pr_overheadO3col_freq2(self, pr_overheadO3col_freq2)
        implicit none
        logical, intent(out) :: pr_overheadO3col_freq2
        type(t_Diagnostics), intent(in) :: self
        pr_overheadO3col_freq2 = self%pr_overheadO3col_freq2
        return
      end subroutine Get_pr_overheadO3col_freq2
!------------------------------------------------------------------------------
      subroutine Get_pr_tropopausePress_freq2(self, pr_tropopausePress_freq2)
        implicit none
        logical, intent(out) :: pr_tropopausePress_freq2
        type(t_Diagnostics), intent(in) :: self
        pr_tropopausePress_freq2 = self%pr_tropopausePress_freq2
        return
      end subroutine Get_pr_tropopausePress_freq2
!------------------------------------------------------------------------------
      subroutine Get_k1_freq3(self, k1_freq3)
        implicit none
        integer, intent(out) :: k1_freq3
        type(t_Diagnostics), intent(in) :: self
        k1_freq3 = self%k1_freq3
        return
      end subroutine Get_k1_freq3
!------------------------------------------------------------------------------
      subroutine Get_k2_freq3(self, k2_freq3)
        implicit none
        integer, intent(out) :: k2_freq3
        type(t_Diagnostics), intent(in) :: self
        k2_freq3 = self%k2_freq3
        return
      end subroutine Get_k2_freq3
!------------------------------------------------------------------------------
      subroutine Get_do_day1_freq3(self, do_day1_freq3)
        implicit none
        logical, intent(out) :: do_day1_freq3
        type(t_Diagnostics), intent(in) :: self
        do_day1_freq3 = self%do_day1_freq3
        return
      end subroutine Get_do_day1_freq3
!------------------------------------------------------------------------------
      subroutine Get_pr_potentialVorticity_freq3(self, pr_potentialVorticity_freq3)
        implicit none
        logical, intent(out) :: pr_potentialVorticity_freq3
        type(t_Diagnostics), intent(in) :: self
        pr_potentialVorticity_freq3 = self%pr_potentialVorticity_freq3
        return
      end subroutine Get_pr_potentialVorticity_freq3
!------------------------------------------------------------------------------
      subroutine Get_pr_freq3_period(self, pr_freq3_period)
        implicit none
        real*8, intent(out) :: pr_freq3_period
        type(t_Diagnostics), intent(in) :: self
        pr_freq3_period = self%pr_freq3_period
        return
      end subroutine Get_pr_freq3_period
!------------------------------------------------------------------------------
      subroutine Get_pr_at_time_freq3(self, pr_at_time_freq3)
        implicit none
        real*8, intent(out) :: pr_at_time_freq3
        type(t_Diagnostics), intent(in) :: self
        pr_at_time_freq3 = self%pr_at_time_freq3
        return
      end subroutine Get_pr_at_time_freq3
!------------------------------------------------------------------------------
      subroutine Get_freq3_description(self, freq3_description)
        implicit none
        character(len=*), intent(out) :: freq3_description
        type(t_Diagnostics), intent(in) :: self
        freq3_description(:) = self%freq3_description(:)
        return
      end subroutine Get_freq3_description
!------------------------------------------------------------------------------
      subroutine Get_freq3_name(self, freq3_name)
        implicit none
        character(len=*), intent(out) :: freq3_name
        type(t_Diagnostics), intent(in) :: self
        freq3_name(:) = self%freq3_name(:)
        return
      end subroutine Get_freq3_name
!------------------------------------------------------------------------------
      subroutine Get_iRange_freq3(self, iRange_freq3)
        implicit none
        integer, intent(out) :: iRange_freq3(:)
        type(t_Diagnostics), intent(in) :: self
        iRange_freq3(:) = self%iRange_freq3(:)
        return
      end subroutine Get_iRange_freq3
!------------------------------------------------------------------------------
      subroutine Get_jRange_freq3(self, jRange_freq3)
        implicit none
        integer, intent(out) :: jRange_freq3(:)
        type(t_Diagnostics), intent(in) :: self
        jRange_freq3(:) = self%jRange_freq3(:)
        return
      end subroutine Get_jRange_freq3
!------------------------------------------------------------------------------
      subroutine Get_freq3_species_num(self, freq3_species_num)
        implicit none
        integer, intent(out) :: freq3_species_num
        type(t_Diagnostics), intent(in) :: self
        freq3_species_num = self%freq3_species_num
        return
      end subroutine Get_freq3_species_num
!------------------------------------------------------------------------------
      subroutine Get_freq3_species(self, freq3_species)
        implicit none
        integer, intent(out) :: freq3_species(:)
        type(t_Diagnostics), intent(in) :: self
        freq3_species(:) = self%freq3_species(:)
        return
      end subroutine Get_freq3_species
!------------------------------------------------------------------------------
      subroutine Get_pr_freq3(self, pr_freq3)
        implicit none
        logical, intent(out) :: pr_freq3
        type(t_Diagnostics), intent(in) :: self
        pr_freq3 = self%pr_freq3
        return
      end subroutine Get_pr_freq3
!------------------------------------------------------------------------------
      subroutine Get_do_mean_freq3(self, do_mean_freq3)
        implicit none
        logical, intent(out) :: do_mean_freq3
        type(t_Diagnostics), intent(in) :: self
        do_mean_freq3 = self%do_mean_freq3
        return
      end subroutine Get_do_mean_freq3
!------------------------------------------------------------------------------
      subroutine Get_pr_mass_freq3(self, pr_mass_freq3)
        implicit none
        logical, intent(out) :: pr_mass_freq3
        type(t_Diagnostics), intent(in) :: self
        pr_mass_freq3 = self%pr_mass_freq3
        return
      end subroutine Get_pr_mass_freq3
!------------------------------------------------------------------------------
      subroutine Get_pr_psf_freq3(self, pr_psf_freq3)
        implicit none
        logical, intent(out) :: pr_psf_freq3
        type(t_Diagnostics), intent(in) :: self
        pr_psf_freq3 = self%pr_psf_freq3
        return
      end subroutine Get_pr_psf_freq3
!------------------------------------------------------------------------------
      subroutine Get_pr_kel_freq3(self, pr_kel_freq3)
        implicit none
        logical, intent(out) :: pr_kel_freq3
        type(t_Diagnostics), intent(in) :: self
        pr_kel_freq3 = self%pr_kel_freq3
        return
      end subroutine Get_pr_kel_freq3
!------------------------------------------------------------------------------
      subroutine Get_pr_const_freq3(self, pr_const_freq3)
        implicit none
        logical, intent(out) :: pr_const_freq3
        type(t_Diagnostics), intent(in) :: self
        pr_const_freq3 = self%pr_const_freq3
        return
      end subroutine Get_pr_const_freq3
!------------------------------------------------------------------------------
      subroutine Get_pr_rel_hum_freq3(self, pr_rel_hum_freq3)
        implicit none
        logical, intent(out) :: pr_rel_hum_freq3
        type(t_Diagnostics), intent(in) :: self
        pr_rel_hum_freq3 = self%pr_rel_hum_freq3
        return
      end subroutine Get_pr_rel_hum_freq3
!------------------------------------------------------------------------------
      subroutine Get_pr_metwater_freq3(self, pr_metwater_freq3)
        implicit none
        logical, intent(out) :: pr_metwater_freq3
        type(t_Diagnostics), intent(in) :: self
        pr_metwater_freq3 = self%pr_metwater_freq3
        return
      end subroutine Get_pr_metwater_freq3
!------------------------------------------------------------------------------
      subroutine Get_do_last_tstep_freq3(self, do_last_tstep_freq3)
        implicit none
        logical, intent(out) :: do_last_tstep_freq3
        type(t_Diagnostics), intent(in) :: self
        do_last_tstep_freq3 = self%do_last_tstep_freq3
        return
      end subroutine Get_do_last_tstep_freq3
!------------------------------------------------------------------------------
      subroutine Get_pr_grid_height_freq3(self, pr_grid_height_freq3)
        implicit none
        logical, intent(out) :: pr_grid_height_freq3
        type(t_Diagnostics), intent(in) :: self
        pr_grid_height_freq3 = self%pr_grid_height_freq3
        return
      end subroutine Get_pr_grid_height_freq3
!------------------------------------------------------------------------------
      subroutine Get_pr_const_column_freq3(self, pr_const_column_freq3)
        implicit none
        logical, intent(out) :: pr_const_column_freq3
        type(t_Diagnostics), intent(in) :: self
        pr_const_column_freq3 = self%pr_const_column_freq3
        return
      end subroutine Get_pr_const_column_freq3
!------------------------------------------------------------------------------
      subroutine Get_pr_const_surface_freq3(self, pr_const_surface_freq3)
        implicit none
        logical, intent(out) :: pr_const_surface_freq3
        type(t_Diagnostics), intent(in) :: self
        pr_const_surface_freq3 = self%pr_const_surface_freq3
        return
      end subroutine Get_pr_const_surface_freq3
!------------------------------------------------------------------------------
      subroutine Get_pr_overheadO3col_freq3(self, pr_overheadO3col_freq3)
        implicit none
        logical, intent(out) :: pr_overheadO3col_freq3
        type(t_Diagnostics), intent(in) :: self
        pr_overheadO3col_freq3 = self%pr_overheadO3col_freq3
        return
      end subroutine Get_pr_overheadO3col_freq3
!------------------------------------------------------------------------------
      subroutine Get_pr_tropopausePress_freq3(self, pr_tropopausePress_freq3)
        implicit none
        logical, intent(out) :: pr_tropopausePress_freq3
        type(t_Diagnostics), intent(in) :: self
        pr_tropopausePress_freq3 = self%pr_tropopausePress_freq3
        return
      end subroutine Get_pr_tropopausePress_freq3
!------------------------------------------------------------------------------
      subroutine Get_k1_freq4(self, k1_freq4)
        implicit none
        integer, intent(out) :: k1_freq4
        type(t_Diagnostics), intent(in) :: self
        k1_freq4 = self%k1_freq4
        return
      end subroutine Get_k1_freq4
!------------------------------------------------------------------------------
      subroutine Get_k2_freq4(self, k2_freq4)
        implicit none
        integer, intent(out) :: k2_freq4
        type(t_Diagnostics), intent(in) :: self
        k2_freq4 = self%k2_freq4
        return
      end subroutine Get_k2_freq4
!------------------------------------------------------------------------------
      subroutine Get_do_day1_freq4(self, do_day1_freq4)
        implicit none
        logical, intent(out) :: do_day1_freq4
        type(t_Diagnostics), intent(in) :: self
        do_day1_freq4 = self%do_day1_freq4
        return
      end subroutine Get_do_day1_freq4
!------------------------------------------------------------------------------
      subroutine Get_pr_potentialVorticity_freq4(self, pr_potentialVorticity_freq4)
        implicit none
        logical, intent(out) :: pr_potentialVorticity_freq4
        type(t_Diagnostics), intent(in) :: self
        pr_potentialVorticity_freq4 = self%pr_potentialVorticity_freq4
        return
      end subroutine Get_pr_potentialVorticity_freq4
!------------------------------------------------------------------------------
      subroutine Get_pr_freq4_period(self, pr_freq4_period)
        implicit none
        real*8, intent(out) :: pr_freq4_period
        type(t_Diagnostics), intent(in) :: self
        pr_freq4_period = self%pr_freq4_period
        return
      end subroutine Get_pr_freq4_period
!------------------------------------------------------------------------------
      subroutine Get_pr_at_time_freq4(self, pr_at_time_freq4)
        implicit none
        real*8, intent(out) :: pr_at_time_freq4
        type(t_Diagnostics), intent(in) :: self
        pr_at_time_freq4 = self%pr_at_time_freq4
        return
      end subroutine Get_pr_at_time_freq4
!------------------------------------------------------------------------------
      subroutine Get_freq4_description(self, freq4_description)
        implicit none
        character(len=*), intent(out) :: freq4_description
        type(t_Diagnostics), intent(in) :: self
        freq4_description(:) = self%freq4_description(:)
        return
      end subroutine Get_freq4_description
!------------------------------------------------------------------------------
      subroutine Get_freq4_name(self, freq4_name)
        implicit none
        character(len=*), intent(out) :: freq4_name
        type(t_Diagnostics), intent(in) :: self
        freq4_name(:) = self%freq4_name(:)
        return
      end subroutine Get_freq4_name
!------------------------------------------------------------------------------
      subroutine Get_iRange_freq4(self, iRange_freq4)
        implicit none
        integer, intent(out) :: iRange_freq4(:)
        type(t_Diagnostics), intent(in) :: self
        iRange_freq4(:) = self%iRange_freq4(:)
        return
      end subroutine Get_iRange_freq4
!------------------------------------------------------------------------------
      subroutine Get_jRange_freq4(self, jRange_freq4)
        implicit none
        integer, intent(out) :: jRange_freq4(:)
        type(t_Diagnostics), intent(in) :: self
        jRange_freq4(:) = self%jRange_freq4(:)
        return
      end subroutine Get_jRange_freq4
!------------------------------------------------------------------------------
      subroutine Get_freq4_species_num(self, freq4_species_num)
        implicit none
        integer, intent(out) :: freq4_species_num
        type(t_Diagnostics), intent(in) :: self
        freq4_species_num = self%freq4_species_num
        return
      end subroutine Get_freq4_species_num
!------------------------------------------------------------------------------
      subroutine Get_freq4_species(self, freq4_species)
        implicit none
        integer, intent(out) :: freq4_species(:)
        type(t_Diagnostics), intent(in) :: self
        freq4_species(:) = self%freq4_species(:)
        return
      end subroutine Get_freq4_species
!------------------------------------------------------------------------------
      subroutine Get_pr_freq4(self, pr_freq4)
        implicit none
        logical, intent(out) :: pr_freq4
        type(t_Diagnostics), intent(in) :: self
        pr_freq4 = self%pr_freq4
        return
      end subroutine Get_pr_freq4
!------------------------------------------------------------------------------
      subroutine Get_do_mean_freq4(self, do_mean_freq4)
        implicit none
        logical, intent(out) :: do_mean_freq4
        type(t_Diagnostics), intent(in) :: self
        do_mean_freq4 = self%do_mean_freq4
        return
      end subroutine Get_do_mean_freq4
!------------------------------------------------------------------------------
      subroutine Get_pr_mass_freq4(self, pr_mass_freq4)
        implicit none
        logical, intent(out) :: pr_mass_freq4
        type(t_Diagnostics), intent(in) :: self
        pr_mass_freq4 = self%pr_mass_freq4
        return
      end subroutine Get_pr_mass_freq4
!------------------------------------------------------------------------------
      subroutine Get_pr_psf_freq4(self, pr_psf_freq4)
        implicit none
        logical, intent(out) :: pr_psf_freq4
        type(t_Diagnostics), intent(in) :: self
        pr_psf_freq4 = self%pr_psf_freq4
        return
      end subroutine Get_pr_psf_freq4
!------------------------------------------------------------------------------
      subroutine Get_pr_kel_freq4(self, pr_kel_freq4)
        implicit none
        logical, intent(out) :: pr_kel_freq4
        type(t_Diagnostics), intent(in) :: self
        pr_kel_freq4 = self%pr_kel_freq4
        return
      end subroutine Get_pr_kel_freq4
!------------------------------------------------------------------------------
      subroutine Get_pr_const_freq4(self, pr_const_freq4)
        implicit none
        logical, intent(out) :: pr_const_freq4
        type(t_Diagnostics), intent(in) :: self
        pr_const_freq4 = self%pr_const_freq4
        return
      end subroutine Get_pr_const_freq4
!------------------------------------------------------------------------------
      subroutine Get_pr_rel_hum_freq4(self, pr_rel_hum_freq4)
        implicit none
        logical, intent(out) :: pr_rel_hum_freq4
        type(t_Diagnostics), intent(in) :: self
        pr_rel_hum_freq4 = self%pr_rel_hum_freq4
        return
      end subroutine Get_pr_rel_hum_freq4
!------------------------------------------------------------------------------
      subroutine Get_pr_metwater_freq4(self, pr_metwater_freq4)
        implicit none
        logical, intent(out) :: pr_metwater_freq4
        type(t_Diagnostics), intent(in) :: self
        pr_metwater_freq4 = self%pr_metwater_freq4
        return
      end subroutine Get_pr_metwater_freq4
!------------------------------------------------------------------------------
      subroutine Get_do_last_tstep_freq4(self, do_last_tstep_freq4)
        implicit none
        logical, intent(out) :: do_last_tstep_freq4
        type(t_Diagnostics), intent(in) :: self
        do_last_tstep_freq4 = self%do_last_tstep_freq4
        return
      end subroutine Get_do_last_tstep_freq4
!------------------------------------------------------------------------------
      subroutine Get_pr_grid_height_freq4(self, pr_grid_height_freq4)
        implicit none
        logical, intent(out) :: pr_grid_height_freq4
        type(t_Diagnostics), intent(in) :: self
        pr_grid_height_freq4 = self%pr_grid_height_freq4
        return
      end subroutine Get_pr_grid_height_freq4
!------------------------------------------------------------------------------
      subroutine Get_pr_const_column_freq4(self, pr_const_column_freq4)
        implicit none
        logical, intent(out) :: pr_const_column_freq4
        type(t_Diagnostics), intent(in) :: self
        pr_const_column_freq4 = self%pr_const_column_freq4
        return
      end subroutine Get_pr_const_column_freq4
!------------------------------------------------------------------------------
      subroutine Get_pr_const_surface_freq4(self, pr_const_surface_freq4)
        implicit none
        logical, intent(out) :: pr_const_surface_freq4
        type(t_Diagnostics), intent(in) :: self
        pr_const_surface_freq4 = self%pr_const_surface_freq4
        return
      end subroutine Get_pr_const_surface_freq4
!------------------------------------------------------------------------------
      subroutine Get_pr_overheadO3col_freq4(self, pr_overheadO3col_freq4)
        implicit none
        logical, intent(out) :: pr_overheadO3col_freq4
        type(t_Diagnostics), intent(in) :: self
        pr_overheadO3col_freq4 = self%pr_overheadO3col_freq4
        return
      end subroutine Get_pr_overheadO3col_freq4
!------------------------------------------------------------------------------
      subroutine Get_pr_tropopausePress_freq4(self, pr_tropopausePress_freq4)
        implicit none
        logical, intent(out) :: pr_tropopausePress_freq4
        type(t_Diagnostics), intent(in) :: self
        pr_tropopausePress_freq4 = self%pr_tropopausePress_freq4
        return
      end subroutine Get_pr_tropopausePress_freq4
!------------------------------------------------------------------------------
      subroutine Get_pr_time(self, pr_time)
        implicit none
        logical, intent(out) :: pr_time
        type(t_Diagnostics), intent(in) :: self
        pr_time = self%pr_time
        return
      end subroutine Get_pr_time
!------------------------------------------------------------------------------
      subroutine Get_asclun(self, asclun)
        implicit none
        integer, intent(out) :: asclun
        type(t_Diagnostics), intent(in) :: self
        asclun = self%asclun
        return
      end subroutine Get_asclun
!------------------------------------------------------------------------------
      subroutine Get_ascii_out_i(self, ascii_out_i)
        implicit none
        integer, intent(out) :: ascii_out_i
        type(t_Diagnostics), intent(in) :: self
        ascii_out_i = self%ascii_out_i
        return
      end subroutine Get_ascii_out_i
!------------------------------------------------------------------------------
      subroutine Get_ascii_out_n(self, ascii_out_n)
        implicit none
        integer, intent(out) :: ascii_out_n
        type(t_Diagnostics), intent(in) :: self
        ascii_out_n = self%ascii_out_n
        return
      end subroutine Get_ascii_out_n
!------------------------------------------------------------------------------
      subroutine Get_pr_ascii_step_interval(self, pr_ascii_step_interval)
        implicit none
        integer, intent(out) :: pr_ascii_step_interval
        type(t_Diagnostics), intent(in) :: self
        pr_ascii_step_interval = self%pr_ascii_step_interval
        return
      end subroutine Get_pr_ascii_step_interval
!------------------------------------------------------------------------------
      subroutine Get_pr_ascii(self, pr_ascii)
        implicit none
        logical, intent(out) :: pr_ascii
        type(t_Diagnostics), intent(in) :: self
        pr_ascii = self%pr_ascii
        return
      end subroutine Get_pr_ascii
!------------------------------------------------------------------------------
      subroutine Get_pr_ascii1(self, pr_ascii1)
        implicit none
        logical, intent(out) :: pr_ascii1
        type(t_Diagnostics), intent(in) :: self
        pr_ascii1 = self%pr_ascii1
        return
      end subroutine Get_pr_ascii1
!------------------------------------------------------------------------------
      subroutine Get_pr_ascii2(self, pr_ascii2)
        implicit none
        logical, intent(out) :: pr_ascii2
        type(t_Diagnostics), intent(in) :: self
        pr_ascii2 = self%pr_ascii2
        return
      end subroutine Get_pr_ascii2
!------------------------------------------------------------------------------
      subroutine Get_pr_ascii3(self, pr_ascii3)
        implicit none
        logical, intent(out) :: pr_ascii3
        type(t_Diagnostics), intent(in) :: self
        pr_ascii3 = self%pr_ascii3
        return
      end subroutine Get_pr_ascii3
!------------------------------------------------------------------------------
      subroutine Get_pr_ascii4(self, pr_ascii4)
        implicit none
        logical, intent(out) :: pr_ascii4
        type(t_Diagnostics), intent(in) :: self
        pr_ascii4 = self%pr_ascii4
        return
      end subroutine Get_pr_ascii4
!------------------------------------------------------------------------------
      subroutine Get_pr_ascii5(self, pr_ascii5)
        implicit none
        logical, intent(out) :: pr_ascii5
        type(t_Diagnostics), intent(in) :: self
        pr_ascii5 = self%pr_ascii5
        return
      end subroutine Get_pr_ascii5
!------------------------------------------------------------------------------
      subroutine Get_ncid_col(self, ncid, isite)
      implicit none
      integer, intent(in) :: isite
      type(t_Diagnostics), intent(in) :: self
      integer, intent(out) :: ncid
      ncid = self%ncid_col(isite)
      return
      end subroutine Get_ncid_col
!------------------------------------------------------------------------------
      subroutine Set_ncid_col(self, ncid, isite)
      implicit none
      integer, intent(in) :: isite
      integer, intent(in) :: ncid
      type(t_Diagnostics), intent(inOut) :: self
      self%ncid_col(isite) = ncid
      return
      end subroutine Set_ncid_col
!------------------------------------------------------------------------------
      subroutine Get_pr_col_diag(self, pr_col_diag)
      implicit none
      type(t_Diagnostics), intent(in) :: self
      logical, intent(out) :: pr_col_diag
      pr_col_diag = self%pr_col_diag
      return
      end subroutine Get_pr_col_diag
!------------------------------------------------------------------------------
      subroutine Get_col_diag_num(self, col_diag_num)
      implicit none
      type(t_Diagnostics), intent(in) :: self
      integer, intent(out) :: col_diag_num
      col_diag_num = self%col_diag_num
      return
      end subroutine Get_col_diag_num
!------------------------------------------------------------------------------
      subroutine Get_col_diag_pres_num(self, col_diag_pres_num)
      implicit none
      type(t_Diagnostics), intent(in) :: self
      integer, intent(out) :: col_diag_pres_num
      col_diag_pres_num = self%col_diag_pres_num
      return
      end subroutine Get_col_diag_pres_num
!------------------------------------------------------------------------------
      subroutine Get_col_diag_species(self, col_diag_species)
      implicit none
      type(t_Diagnostics), intent(in) :: self
      integer, intent(out) :: col_diag_species(:)
      col_diag_species(:) = self%col_diag_species(:)
      return
      end subroutine Get_col_diag_species
!------------------------------------------------------------------------------
      subroutine Get_col_diag_species_num(self, col_diag_species_num)
      implicit none
      type(t_Diagnostics), intent(in) :: self
      integer, intent(out) :: col_diag_species_num
      col_diag_species_num = self%col_diag_species_num
      return
      end subroutine Get_col_diag_species_num
!------------------------------------------------------------------------------
      subroutine Get_col_diag_period(self, col_diag_period)
      implicit none
      type(t_Diagnostics), intent(in) :: self
      real*8, intent(out) :: col_diag_period
      col_diag_period = self%col_diag_period
      return
      end subroutine Get_col_diag_period
!------------------------------------------------------------------------------
      subroutine Get_col_diag_pres(self, col_diag_pres)
      implicit none
      type(t_Diagnostics), intent(in) :: self
      real*8, intent(out) :: col_diag_pres(:)
      col_diag_pres(:) = self%col_diag_pres(:)
      return
      end subroutine Get_col_diag_pres
!------------------------------------------------------------------------------
      subroutine Set_col_diag_lat_lon(self, col_diag_lat_lon)
      implicit none
      type(t_Diagnostics), intent(inOut) :: self
      real*8, intent(in) :: col_diag_lat_lon(:,:)
      self%col_diag_lat_lon(:,:) = col_diag_lat_lon(:,:)
      return
      end subroutine Set_col_diag_lat_lon
!------------------------------------------------------------------------------
      subroutine Get_col_diag_lat_lon(self, col_diag_lat_lon)
      implicit none
      type(t_Diagnostics), intent(in) :: self
      real*8, intent(out) :: col_diag_lat_lon(:,:)
      col_diag_lat_lon(:,:) = self%col_diag_lat_lon(:,:)
      return
      end subroutine Get_col_diag_lat_lon
!------------------------------------------------------------------------------
      subroutine Get_col_diag_site(self, col_diag_site)
      implicit none
      type(t_Diagnostics), intent(in) :: self
      character(len=24), intent(out) :: col_diag_site(:)
      col_diag_site(:) = self%col_diag_site(:)
      return
      end subroutine Get_col_diag_site
!------------------------------------------------------------------------------
      end module GmiDiagnosticsMethod_mod
