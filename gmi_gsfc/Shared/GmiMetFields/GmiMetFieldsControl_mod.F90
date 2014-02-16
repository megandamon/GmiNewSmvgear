!-----------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!-----------------------------------------------------------------------------
!BOP
!
! !IMODULE: GmiMetFieldsControl_mod
!
! !INTERFACE:
!
#include "GmiESMF_ErrLog.h"
!
      module GmiMetFieldsControl_mod
!
! !USES:
      use Esmf_Mod
      use GmiESMF_ErrorChecking_mod
      use GmiTotalMass_mod       , only : calcTotalMass
      use GmiPrintError_mod, only : GmiPrintError
      use GmiReadList_mod  , only : Read_List
      use GmiWrapMaster_mod, only : wrapMaster_2d, wrapMaster_3du, wrapMaster_3dv
      use m_netcdf_io_read , only : Ncrd_3d, Ncrd_4d
      use m_netcdf_io_get_dimlen, only : Ncget_Unlim_Dimlen
      use m_netcdf_io_open      , only : Ncop_Rd
      use m_netcdf_io_close     , only : Nccl
      use m_netcdf_io_checks    , only : Ncdoes_Var_Exist
      use GmiUtilsMetFields_mod , only : Set_Met1_Simple, Read_Met1, &
     &            Read_Met2, Check_Met1_Range, Get_Next_Met_File, &
     &            Calc_Var_Adv_Tstp
      use GmiTimeControl_mod, only : t_GmiClock, Get_gmiTimeStep, Get_curGmiDate
      use GmiTimeControl_mod, only : Get_totNumDays, Get_begGmiDate
      use GmiGrid_mod, only : t_gmiGrid, Get_numSpecies, Get_i1, Get_i2,       &
     &       Get_ju1, Get_jv1, Get_j2, Get_k1, Get_k2, Get_ilo, Get_ihi, Get_julo, &
     &       Get_jvlo, Get_jhi, Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl,   &
     &       Get_jv1_gl, Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jvlo_gl,     &
     &       Get_jhi_gl, Get_gmi_nborder
      use GmiDomainDecomposition_mod, only : t_gmiDomain
      use GmiDomainDecomposition_mod, only : Get_numLonDomains
      use GmiDomainDecomposition_mod, only : Get_procID, Get_communicatorNorthPole, &
     &       Get_communicatorSouthPole
      use GmiDomainDecomposition_mod, only : Get_mcor, Get_latdeg, Get_londeg, &
     &       Get_dlatr, Get_cosp, Get_coscen, Get_cose, Get_rel_areaGlob, &
     &       Get_geofac, Get_geofac_pc
      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Get_pr_diag
      use GmiESMFrcFileReading_mod, only : rcEsmfReadLogical, rcEsmfReadTable
      use GmiModelData_mod, only : Set_Vert_Dao, Calc_Layer_Dsigma
      use GmiLandWaterIce_mod, only : setLWIflags, setCMIflags
      USE GmiGenerateMetFileList_mod, ONLY : createMetFileList
!
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: t_MetFields, Control_Met1_File_Input, Control_Met2_File_Input
      public  :: initializeMetFields, Set_Met1, Update_Met1, met1Glob2Sub
      public  :: Set_Met2, Update_Met2, finishMet
      public  :: Get_ai, Get_bi, Get_am, Get_bm, Get_dap, Get_dbk
      public  :: Get_pt, Get_ptop, Get_met_opt, Get_mdt, Get_do_cycle_met
      public  :: Get_metdata_name, Get_metdata_name_org, Get_metdata_name_dims
      public  :: Get_metdata_name_model, Get_met_grid_type
      public  :: Get_met_num_infiles, Get_metNumMEGAN, Get_met_infile_num
      public  :: Get_mrnum_in, Get_gwet_opt, Get_do_wind_pole
      public  :: Get_do_timinterp_winds, Get_do_timinterp_met
      public  :: Set_cmi_flags, Get_cmi_flags, Get_lwi_flags, Get_met_infile_names
      public  :: Get_surfTemp15DayAvg, Set_surfTemp15DayAvg
      public  :: Get_m1num_recs, Set_m1num_recs
      public  :: Get_m1rnum_in, Set_m1rnum_in
      public  :: Get_met1_infile_num, Set_met1_infile_num
      public  :: Get_ncid_met1, Set_ncid_met1
      public  :: Get_tmet1_1, Set_tmet1_1
      public  :: Get_tmet1_2, Set_tmet1_2
      public  :: Get_m2num_recs, Set_m2num_recs
      public  :: Get_m2rnum_in, Set_m2rnum_in
      public  :: Get_met2_infile_num, Set_met2_infile_num
      public  :: Get_ncid_met2, Set_ncid_met2
      public  :: Get_tmet2_1, Set_tmet2_1
      public  :: Get_tmet2_2, Set_tmet2_2
      public  :: Get_mass, Set_mass
      public  :: Get_gridBoxHeight, Set_gridBoxHeight
      public  :: Get_tropopausePress, Set_tropopausePress
      public  :: Get_potentialVorticity, Set_potentialVorticity
      public  :: Get_potentialTemp, Set_potentialTemp
      public  :: Get_relativeHumidity, Set_relativeHumidity
      public  :: Get_press3c, Set_press3c, Get_press3e, Set_press3e
      public  :: Get_kel, Get_pctm1, Get_pctm2, Get_uux, Get_vvx
      public  :: Get_xmass, Get_ymass, Get_zmass
!
      public  :: Get_humidity, Get_tau_cloud, Get_max_cloud, Set_pctm2Glob, Set_pctm1Glob
      public  :: Get_totalCloudFraction, Set_totalCloudFraction, Get_pctm2Glob, Get_pctm1Glob
      public  :: Get_kelGlob
      public  :: Get_fracCloudCover, Set_fracCloudCover, Set_tau_cloud
      public  :: Get_radswg, Get_surf_air_temp, Get_surf_rough
      public  :: Get_con_precip, Get_tot_precip, Get_ustar, Get_pbl, Get_cmf
      public  :: Get_dtrn, Get_u10m, Get_v10m, Get_gwet, Get_pardif, Get_pardir
      public  :: Get_kzz, Get_eu, Get_ed, Get_md, Get_zmdu, Get_zmeu, Get_zmed
      public  :: Get_zmmd, Get_zmmu, Get_hkdu, Get_hkeu, Get_hkmu, Set_kzz
      public  :: Set_moistq, Get_moistq, Get_clwc, Get_surf_alb_uv
      public  :: Get_rain, Get_rain_zm, Get_rain_hk, Get_rain_ls, Get_ran_cloud
      public  :: Set_surf_alb_uv, Get_taucli, Get_tauclw
      public  :: Get_sasdir, Get_sasdif, Get_saldir, Get_saldif
      public  :: Set_sasdir, Set_sasdif, Set_saldir, Set_saldif
!
#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"
!
      type t_MetFields
         private
!
         logical :: do_read_met_list   !
         logical :: do_cycle_met       ! cycle metFields data?
         logical :: do_timinterp_met   ! time interpolate met  fields,
                                       ! except winds?
         logical :: do_timinterp_winds ! time interpolate wind fields?
         logical :: do_wind_pole       ! make winds blow over Pole?
         integer :: met_opt            ! metFields input option
         integer :: gwet_opt           ! gwet option
         integer :: m1rnum_in          ! next NetCDF met1 record to read
         integer :: m2rnum_in          ! next NetCDF met2 record to read
         integer :: met_num_infiles    ! total number of metFields input files
         integer :: met1_num_infiles   ! index of current met1 input file
         integer :: met2_num_infiles   ! index of current met2 input file
         real*8  :: mdt                ! time increment for reading new
                                       ! metFields data (s)
         real*8  :: tmet1              ! Time tag of the mrnum_in (s).
         integer :: met_infile_num
         integer :: met1_infile_num    ! index of current met1 input file
         integer :: met2_infile_num    ! index of current met2 input file
         integer :: metNumMEGAN        ! points to the metFields file number
                                       ! to be  read in for the MEGAN emissions
         integer :: mrnum_in
         integer :: m1num_recs         ! number of netCDF met1 records in file
         integer :: m2num_recs         ! number of netCDF met2 records in file
!
         character (len=1):: met_grid_type  ! met grid type, A or C
         character (len=50) :: metdata_name ! met data netcdf file attribute
                                            ! "Met_Dat_Name", e.g.,
                                            ! "NCAR_MATCH_4x5x52"
         character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org ! first  part of metdata_name,
                                                ! e.g., "NCAR"
         character (len=MAX_LENGTH_MET_NAME) :: metdata_name_model ! second part of metdata_name,
                                                  ! e.g., "MATCH"
         character (len=MAX_LENGTH_MET_NAME) :: metdata_name_dims ! third  part of metdata_name,
                                                 ! e.g., "4x5x52"
         character (len=MAX_LENGTH_FILE_NAME) :: met_infile_names(MAX_INFILE_NAMES)
                                      ! metFields input file names
         character (len=MAX_LENGTH_FILE_NAME) :: met_filnam_list ! Name of file to get
                                                        ! names of met input
                                                        ! files from.
         INTEGER :: metStartDate ! starting date (YYYYMMDD) for the metFields
         INTEGER :: metRefTime   ! reference time (HHMMSS)
         LOGICAL :: do_manualReadMetList
         CHARACTER(len=MAX_LENGTH_FILE_NAME) :: metFileTemplate ! file template to be used to create
                                            ! the list of needed metFields files.
!
         integer :: ncid_met1, ncid_met2
         real*8  :: tmet1_1          ! time tag of first  met1 record (s)
         real*8  :: tmet1_2          ! time tag of second met1 record (s)
         real*8  :: tmet2_1          ! time tag of first  met2 record (s)
         real*8  :: tmet2_2          ! time tag of second met2 record (s)
!
         integer :: ktr              ! number of vertical tropospheric  levels
         real*8  :: pt               ! pressure = (am * pt) + (bm * psx) (mb)
         real*8  :: ptop             ! pressure at top edge of top zone  (mb)
!
         real*8, pointer :: ai(:)  => null() ! pressure = (ai * pt) + (bi * psx)
                                             ! ai at zone interface
         real*8, pointer :: bi(:)  => null() ! pressure = (ai * pt) + (bi * psx)
                                             ! bi at zone interface
         real*8, pointer :: am(:)  => null() ! pressure = (am * pt) + (bm * psx)
                                             ! am at zone midpoint
         real*8, pointer :: bm(:)  => null() ! pressure = (am * pt) + (bm * psx)
                                             ! bm at zone midpoint
         real*8, pointer :: dap(:) => null() ! pressure difference across layer
                                             ! from log pressure term
         real*8, pointer :: dbk(:) => null() ! pressure difference across layer
                                             ! from sigma term
!
         !------
         ! Set 1
         !------
         real*8 , pointer :: pctm1  (:,:) => null() ! CTM surface pressure at t1     (mb)
         real*8 , pointer :: pctm2  (:,:) => null() ! CTM surface pressure at t1+tdt (mb)
         real*8 , pointer :: psx    (:,:) => null() ! surface pressure at beginning
                                                    ! of met time period (mb)
         real*8 , pointer :: del_psx(:,:) => null() ! surface pressure at end
                                                    ! of met time period (mb)
         real*8 , pointer :: uux    (:,:,:) => null() ! current horizontal velocity
                                                      ! in zonal (longitude) direction;
                                                      ! known at edges in longitude
                                                      ! direction and centers in latitude
                                                      ! direction (m/s)
         real*8 , pointer :: del_uux(:,:,:) => null() ! horizontal velocity in zonal
                                                      ! (longitude) direction at end of
                                                      ! met time period; known at edges
                                                      ! in longitude direction and
                                                      ! centers in latitude direction (m/s)
         real*8 , pointer :: vvx    (:,:,:) => null() ! current horizontal velocity
                                                      ! in meridional (latitude) direction;
                                                      ! known at edges in latitude
                                                      ! direction and centers in longitude
                                                      ! direction (m/s)
         real*8 , pointer :: del_vvx(:,:,:) => null() ! horizontal velocity in meridional
                                                      ! (latitude) direction at end of
                                                      ! met time period; known at edges
                                                      ! in latitude direction and
                                                      ! centers in longitude direction (m/s)
         real*8 , pointer :: kel    (:,:,:) => null() ! current temperature (degK)
         real*8 , pointer :: del_kel(:,:,:) => null() ! temperature at end of met time
                                                      ! period (degK)
         real*8 , pointer :: xmass  (:,:,:) => null() ! horizontal mass flux in E-W
                                                      ! direction (mb)
         real*8 , pointer :: ymass  (:,:,:) => null() ! horizontal mass flux in N-S
                                                      ! direction (mb)
!
         real*8 , pointer :: zmass  (:,:,:) => null()
!
         real*8 , pointer :: pctm1Glob  (:,:) => null()
         real*8 , pointer :: pctm2Glob  (:,:) => null()
         real*8 , pointer :: del_psxGlob  (:,:) => null()
         real*8 , pointer :: psxGlob  (:,:) => null()
         real*8 , pointer :: uuxGlob  (:,:,:) => null()
         real*8 , pointer :: del_uuxGlob  (:,:,:) => null()
         real*8 , pointer :: vvxGlob  (:,:,:) => null()
         real*8 , pointer :: del_vvxGlob  (:,:,:) => null()
         real*8 , pointer :: kelGlob  (:,:,:) => null()
         real*8 , pointer :: del_kelGlob  (:,:,:) => null()
         real*8 , pointer :: xmassGlob  (:,:,:) => null()
         real*8 , pointer :: ymassGlob  (:,:,:) => null()
         real*8 , pointer :: zmassGlob  (:,:,:) => null()
!
         !------
         ! Set 2
         !------
                          ! Land, water, ice flags.
         integer, pointer :: lwi_flags           (:,:) => null()
         integer, pointer :: cmi_flags           (:,:) => null()
!
         ! NCAR convection =>
!
                          ! entrainment into convective downdraft (s^-1)
         real*8 , pointer :: ed                (:,:,:) => null()
         real*8 , pointer :: del_ed            (:,:,:) => null()
                          ! entrainment into convective updraft   (s^-1)
         real*8 , pointer :: eu                (:,:,:) => null()
         real*8 , pointer :: del_eu            (:,:,:) => null()
                          ! convective mass flux in     downdraft (kg/m^2*s)
         real*8 , pointer :: md                (:,:,:) => null()
         real*8 , pointer :: del_md            (:,:,:) => null()
!
         ! GEOS4 convection =>
                          ! Z-M detrainment into convective updraft   (Pa/s)
         real*8 , pointer :: zmdu              (:,:,:) => null()
         real*8 , pointer :: del_zmdu          (:,:,:) => null()
                          ! Z-M entrainment into convective updraft   (Pa/s)
         real*8 , pointer :: zmeu              (:,:,:) => null()
         real*8 , pointer :: del_zmeu          (:,:,:) => null()
                          ! Z-M entrainment into convective downdraft (Pa/s)
         real*8 , pointer :: zmed              (:,:,:) => null()
         real*8 , pointer :: del_zmed          (:,:,:) => null()
                          ! Z-M convective mass flux in     downdraft (Pa/s)
         real*8 , pointer :: zmmd              (:,:,:) => null()
         real*8 , pointer :: del_zmmd          (:,:,:) => null()
                          ! Z-M convective mass flux in     updraft   (Pa/s)
         real*8 , pointer :: zmmu              (:,:,:) => null()
         real*8 , pointer :: del_zmmu          (:,:,:) => null()
!
                          ! Hack detrainment into convective updraft (Pa/s)
         real*8 , pointer :: hkdu              (:,:,:) => null()
         real*8 , pointer :: del_hkdu          (:,:,:) => null()
                          ! Hack entrainment into convective updraft (Pa/s)
         real*8 , pointer :: hkeu              (:,:,:) => null()
         real*8 , pointer :: del_hkeu          (:,:,:) => null()
                          ! Hack convective mass flux in     updraft (Pa/s)
         real*8 , pointer :: hkmu              (:,:,:) => null()
         real*8 , pointer :: del_hkmu          (:,:,:) => null()
!
                          ! 10 meter Zonal Wind      (m/s)
         real*8 , pointer :: u10m                (:,:) => null()
         real*8 , pointer :: del_u10m            (:,:) => null()
                          ! 10 meter Meridional Wind (m/s)
         real*8 , pointer :: v10m                (:,:) => null()
         real*8 , pointer :: del_v10m            (:,:) => null()
                          ! Root zone soil wetness   (fraction)
         real*8 , pointer :: gwet                (:,:) => null()
         real*8 , pointer :: del_gwet            (:,:) => null()
                          ! Diffuse photosynthetically active radiation (0.35-0.70 um)
         real*8 , pointer :: pardif              (:,:) => null()
         real*8 , pointer :: del_pardif          (:,:) => null()
                          ! Direct  photosynthetically active radiation (0.35-0.70 um)
         real*8 , pointer :: pardir              (:,:) => null()
         real*8 , pointer :: del_pardir          (:,:) => null()
!
                          ! convective precipitation (mm/day)
         real*8 , pointer :: con_precip          (:,:) => null()
         real*8 , pointer :: del_con_precip      (:,:) => null()
                          ! total precipitation      (mm/day)
         real*8 , pointer :: tot_precip          (:,:) => null()
         real*8 , pointer :: del_tot_precip      (:,:) => null()
!
                          ! ground temperature       (degK)
         real*8 , pointer :: grnd_temp           (:,:) => null()
         real*8 , pointer :: del_grnd_temp       (:,:) => null()
                          ! planetary boundary layer depth (mb)
         real*8 , pointer :: pbl                 (:,:) => null()
         real*8 , pointer :: del_pbl             (:,:) => null()
                          ! net downward shortwave radiation at ground (W/m^2)
         real*8 , pointer :: radswg              (:,:) => null()
         real*8 , pointer :: del_radswg          (:,:) => null()
!
                          ! surface albedo for diffuse light (nr IR)  (0-1)
         real*8 , pointer :: saldif              (:,:) => null()
         real*8 , pointer :: del_saldif          (:,:) => null()
                          ! surface albedo for direct  light (nr IR)  (0-1)
         real*8 , pointer :: saldir              (:,:) => null()
         real*8 , pointer :: del_saldir          (:,:) => null()
                          ! surface albedo for diffuse light (uv/vis) (0-1)
         real*8 , pointer :: sasdif              (:,:) => null()
         real*8 , pointer :: del_sasdif          (:,:) => null()
                          ! surface albedo for direct  light (uv/vis) (0-1)
         real*8 , pointer :: sasdir              (:,:) => null()
         real*8 , pointer :: del_sasdir          (:,:) => null()
!
                          ! surface air temperature  (degK)
         real*8 , pointer :: surf_air_temp       (:,:) => null()
         real*8 , pointer :: del_surf_air_temp   (:,:) => null()
                          ! bulk surface albedo      (fraction 0-1)
         real*8 , pointer :: surf_alb_uv         (:,:) => null()
         real*8 , pointer :: del_surf_alb_uv     (:,:) => null()
                          ! surface roughness        (m)
         real*8 , pointer :: surf_rough          (:,:) => null()
         real*8 , pointer :: del_surf_rough      (:,:) => null()
                          ! ustar (m/s)
         real*8 , pointer :: ustar               (:,:) => null()
         real*8 , pointer :: del_ustar           (:,:) => null()
!
                          ! convective cloud fraction from  CCM3
         real*8 , pointer :: aconv             (:,:,:) => null()
         real*8 , pointer :: del_aconv         (:,:,:) => null()
                          ! stratiform cloud fraction from  CCM3
         real*8 , pointer :: astrat            (:,:,:) => null()
         real*8 , pointer :: del_astrat        (:,:,:) => null()
                          ! cloud liquid water content from CCM3 (gm/m^3)
         real*8 , pointer :: clwc              (:,:,:) => null()
         real*8 , pointer :: del_clwc          (:,:,:) => null()
                          ! convective mass flux (kg/m^2*s)
         real*8 , pointer :: cmf               (:,:,:) => null()
         real*8 , pointer :: del_cmf           (:,:,:) => null()
                          ! detrainment rate  (DAO:kg/m^2*s, NCAR:s^-1))
         real*8 , pointer :: dtrn              (:,:,:) => null()
         real*8 , pointer :: del_dtrn          (:,:,:) => null()
                          ! specific humidity (g/kg)
         real*8 , pointer :: humidity          (:,:,:) => null()
         real*8 , pointer :: del_humidity      (:,:,:) => null()
                          ! array of vertical diffusion coefficients (m^2/s)
         real*8 , pointer :: kzz               (:,:,:) => null()
         real*8 , pointer :: del_kzz           (:,:,:) => null()
!
                          ! convective mass flux in updraft (Pa/s) (from CCM3)
         real*8 , pointer :: cldmas            (:,:,:) => null()
         real*8 , pointer :: del_cldmas        (:,:,:) => null()
                          ! vertical pressure velocity      (Pa/s) (from CCM3)
         real*8 , pointer :: omega             (:,:,:) => null()
         real*8 , pointer :: del_omega         (:,:,:) => null()
!
                          ! maximum overlap cloud fraction for LW
         real*8 , pointer :: max_cloud         (:,:,:) => null()
         real*8 , pointer :: del_max_cloud     (:,:,:) => null()
                          ! random  overlap cloud fraction for LW
         real*8 , pointer :: ran_cloud         (:,:,:) => null()
         real*8 , pointer :: del_ran_cloud     (:,:,:) => null()
                          ! cloud optical thickness in visible
         real*8 , pointer :: tau_cloud         (:,:,:) => null()
         real*8 , pointer :: del_tau_cloud     (:,:,:) => null()
                          ! fractional cloud cover
         real*8 , pointer :: fracCloudCover      (:,:) => null()
         real*8 , pointer :: totalCloudFraction(:,:,:) => null()
!
                          ! moisture changes due to wet processes (g/kg/day)
         real*8 , pointer :: moistq            (:,:,:) => null()
         real*8 , pointer :: del_moistq        (:,:,:) => null()
                          ! rainfall across cell edges (mm/day)
         real*8 , pointer :: rain              (:,:,:) => null()
         real*8 , pointer :: del_rain          (:,:,:) => null()
                          ! rain production due to lasrge-scale processes (kg/kg/sec)
         real*8 , pointer :: rain_ls           (:,:,:) => null()
         real*8 , pointer :: del_rain_ls       (:,:,:) => null()
                          ! rain production due to deep conv processes (kg/kg/sec)
         real*8 , pointer :: rain_zm           (:,:,:) => null()
         real*8 , pointer :: del_rain_zm       (:,:,:) => null()
                          ! rain production due to shal conv processes (kg/kg/sec)
         real*8 , pointer :: rain_hk           (:,:,:) => null()
         real*8 , pointer :: del_rain_hk       (:,:,:) => null()
!
                          ! optical_thickness_for_ice_clouds
         real*8 , pointer :: taucli            (:,:,:) => null()
         real*8 , pointer :: del_taucli        (:,:,:) => null()
                          ! optical_thickness_for_liquid_clouds
         real*8 , pointer :: tauclw            (:,:,:) => null()
         real*8 , pointer :: del_tauclw        (:,:,:) => null()
                          ! relative humidity from CCM3 (fraction)
         real*8 , pointer :: relhum            (:,:,:) => null()
         real*8 , pointer :: del_relhum        (:,:,:) => null()
                          ! relative humidity of clear region from CCM3 (fraction)
         real*8 , pointer :: rhclear           (:,:,:) => null()
         real*8 , pointer :: del_rhclear       (:,:,:) => null()
!
        !-------------------
        ! Derived quantities
        !-------------------
         real*8 , pointer :: surfTemp15DayAvg  (:,:)   => null()
         real*8 , pointer :: tropopausePress   (:,:)   => null()
         real*8 , pointer :: mass              (:,:,:) => null()
         real*8 , pointer :: gridBoxHeight     (:,:,:) => null()
         real*8 , pointer :: potentialTemp     (:,:,:) => null()
         real*8 , pointer :: relativeHumidity  (:,:,:) => null()
         real*8 , pointer :: potentialVorticity(:,:,:) => null()
!
                           ! atmospheric pressure at the center of each grid box (mb)
         real*8 , pointer :: press3c           (:,:,:) => null()
                           ! atmospheric pressure at the edge   of each grid box (mb)
         real*8 , pointer :: press3e           (:,:,:) => null()
!
      end type t_MetFields
!
      integer, save :: ilo, ihi, julo, jvlo, jhi, k1, k2, i1, i2, ju1, jv1, j2
      integer, save :: ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl
      integer, save :: i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl
      integer, save :: procID, numLonDomains, gmi_nborder
      logical, save :: pr_diag, do_wetdep, do_drydep
      integer, save :: chem_opt, sfalbedo_opt, uvalbedo_opt
      integer, save :: convec_opt, emiss_in_opt, lightning_opt
      integer, save :: commu_npole, commu_spole
      logical       :: doMEGANemission
      real*8 , save :: tdt
!EOP
!-----------------------------------------------------------------------------
      contains
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
!
      subroutine initializeMetFields(self, gmiGrid, gmiDomain, gmiClock, &
     &                     Diagnostics, config)
!
      implicit none
!
! !INPUT PARAMETERS:
      type (t_gmiClock), intent(in) :: gmiClock
      type (t_gmiGrid), intent(in) :: gmiGrid
      type (t_gmiDomain), intent(in) :: gmiDomain
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT
      type (Esmf_Config), intent(inOut) :: config
      type (t_metFields), intent(inOut) :: self
!
! !LOCAL VARIABLES:
      integer :: RC, STATUS
      character(len=ESMF_MAXSTR) :: IAm , err_msg
!
!EOP
!-----------------------------------------------------------------------------
!BOC
      Iam = "initializeMetFields"
!
      call Get_procID (gmiDomain, procID)
      call Get_pr_diag (Diagnostics, pr_diag)
!
      if (pr_diag) Write(6,*) IAm, 'called by ', procID
!
      call Get_i1     (gmiGrid, i1   )
      call Get_i2     (gmiGrid, i2   )
      call Get_ju1    (gmiGrid, ju1  )
      call Get_jv1    (gmiGrid, jv1  )
      call Get_j2     (gmiGrid, j2   )
      call Get_k1     (gmiGrid, k1   )
      call Get_k2     (gmiGrid, k2   )
      call Get_ilo    (gmiGrid, ilo   )
      call Get_ihi    (gmiGrid, ihi   )
      call Get_julo   (gmiGrid, julo  )
      call Get_jvlo   (gmiGrid, jvlo  )
      call Get_jhi    (gmiGrid, jhi   )
      call Get_i1_gl  (gmiGrid, i1_gl )
      call Get_i2_gl  (gmiGrid, i2_gl )
      call Get_ju1_gl (gmiGrid, ju1_gl )
      call Get_jv1_gl (gmiGrid, jv1_gl )
      call Get_j2_gl  (gmiGrid, j2_gl )
      call Get_ilo_gl (gmiGrid, ilo_gl )
      call Get_ihi_gl (gmiGrid, ihi_gl )
      call Get_julo_gl(gmiGrid, julo_gl )
      call Get_jvlo_gl(gmiGrid, jvlo_gl )
      call Get_jhi_gl (gmiGrid, jhi_gl)
      call Get_gmi_nborder(gmiGrid, gmi_nborder)
!
      call Get_numLonDomains        (gmiDomain, numLonDomains)
      call Get_communicatorNorthPole(gmiDomain, commu_npole)
      call Get_communicatorSouthPole(gmiDomain, commu_spole)
!
      call Get_gmiTimeStep(gmiClock, tdt)
!
      !###################################################
      ! Get non metFields variables from the resource file
      !###################################################
!
      call rcEsmfReadLogical(config, doMEGANemission, "doMEGANemission:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)
!
      call rcEsmfReadLogical(config, do_drydep, "do_drydep:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)
!
      call rcEsmfReadLogical(config, do_wetdep, "do_wetdep:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, sfalbedo_opt, &
     &                label   = "sfalbedo_opt:",&
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, uvalbedo_opt, &
     &                label   = "uvalbedo_opt:",&
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, chem_opt, &
     &                label   = "chem_opt:",&
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, convec_opt, &
     &                label  = "convec_opt:",&
     &               default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, emiss_in_opt, &
     &                label  = "emiss_in_opt:",&
     &               default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, lightning_opt, &
     &                label  = "lightning_opt:",&
     &               default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
      !#############################
      ! Read resource file variables
      !#############################
!
      call readMetFieldsResourceFile (self, gmiDomain, Diagnostics, &
                                      gmiClock, config)
!
      !###################
      ! Allocate variables
      !###################
!
      call allocateMetFields(self)
!
      !#####################
      ! Initialize variables
      !#####################
!
      call Set_Vert_Dao (self%metdata_name, self%ktr, self%ptop, self%pt, &
     &                   self%ai, self%bi, self%am, self%bm, k1, k2)
!
      call Calc_Layer_Dsigma (self%pt, self%ai, self%bi, self%dap, self%dbk, &
     &                   k1, k2)
!
      return
!
      end subroutine initializeMetFields
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
!
      subroutine allocateMetFields(self)
!
      implicit none
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_metFields), intent(inOut) :: self
!
! !LOCAL VARIABLES:
      integer :: RC, STATUS
      character(len=ESMF_MAXSTR) :: IAm
!EOP
!-----------------------------------------------------------------------------
!BOC
      Iam = "allocateMetFields"
!
      if (pr_diag) Write(6,*) IAm, 'called by ', procID
!
      !###################
      ! Allocate Variables
      !###################
!
      Allocate (self%ai (k1-1:k2))
      Allocate (self%bi (k1-1:k2))
      Allocate (self%am (k1:k2))
      Allocate (self%bm (k1:k2))
      Allocate (self%dap(k1:k2))
      Allocate (self%dbk(k1:k2))
      self%ai  = 0.0d0
      self%bi  = 0.0d0
      self%am  = 0.0d0
      self%bm  = 0.0d0
      self%dap = 0.0d0
      self%dbk = 0.0d0
!
      Allocate (self%pctm1(ilo:ihi, julo:jhi))
      Allocate (self%pctm2(ilo:ihi, julo:jhi))
      self%pctm1 = 0.0d0
      self%pctm2 = 0.0d0
!
      Allocate (self%xmass(ilo:ihi, julo:jhi, k1:k2))
      Allocate (self%ymass(ilo:ihi, julo:jhi, k1:k2))
      Allocate (self%zmass(i1:i2, ju1:j2, k1:k2))
      self%xmass = 0.0d0
      self%ymass = 0.0d0
      self%zmass = 0.0d0
!
      Allocate (self%psx    (ilo:ihi, julo:jhi))
      Allocate (self%del_psx(ilo:ihi, julo:jhi))
      self%psx = 0.0d0
      self%del_psx = 0.0d0
!
      Allocate (self%kel    (ilo:ihi, julo:jhi, k1:k2))
      Allocate (self%del_kel(ilo:ihi, julo:jhi, k1:k2))
      Allocate (self%uux    (ilo:ihi, julo:jhi, k1:k2))
      Allocate (self%del_uux(ilo:ihi, julo:jhi, k1:k2))
      self%kel = 0.0d0
      self%del_kel = 0.0d0
      self%uux = 0.0d0
      self%del_uux = 0.0d0
!
      Allocate (self%vvx    (ilo:ihi, jvlo:jhi, k1:k2))
      Allocate (self%del_vvx(ilo:ihi, jvlo:jhi, k1:k2))
      self%vvx = 0.0d0
      self%del_vvx = 0.0d0
!
      Allocate (self%pctm1Glob(ilo_gl:ihi_gl, julo_gl:jhi_gl))
      Allocate (self%pctm2Glob(ilo_gl:ihi_gl, julo_gl:jhi_gl))
      self%pctm1Glob = 0.0d0
      self%pctm2Glob = 0.0d0
!
      Allocate (self%xmassGlob(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2))
      Allocate (self%ymassGlob(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2))
      Allocate (self%zmassGlob(i1_gl:i2_gl, ju1_gl:j2_gl, k1:k2))
      self%xmassGlob = 0.0d0
      self%ymassGlob = 0.0d0
      self%zmassGlob = 0.0d0
!
      Allocate (self%psxGlob    (ilo_gl:ihi_gl, julo_gl:jhi_gl))
      Allocate (self%del_psxGlob(ilo_gl:ihi_gl, julo_gl:jhi_gl))
      self%psxGlob = 0.0d0
      self%del_psxGlob = 0.0d0
!
      Allocate (self%kelGlob    (ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2))
      Allocate (self%del_kelGlob(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2))
      Allocate (self%uuxGlob    (ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2))
      Allocate (self%del_uuxGlob(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2))
      self%kelGlob = 0.0d0
      self%del_kelGlob = 0.0d0
      self%uuxGlob = 0.0d0
      self%del_uuxGlob = 0.0d0
!
      Allocate (self%vvxGlob    (ilo_gl:ihi_gl, jvlo_gl:jhi_gl, k1:k2))
      Allocate (self%del_vvxGlob(ilo_gl:ihi_gl, jvlo_gl:jhi_gl, k1:k2))
      self%vvxGlob = 0.0d0
      self%del_vvxGlob = 0.0d0
!
      !----------------
      ! Set 2 variables
      !----------------
!
      if ((chem_opt     /= 0) .or.  &
     &     (convec_opt   /= 0) .or. (emiss_in_opt /= 0)) then
          Allocate (self%lwi_flags(i1:i2, ju1:j2))
          self%lwi_flags = 0
      end if
!
      if(lightning_opt == 1) then
         Allocate (self%cmi_flags(i1:i2, ju1:j2))
         self%cmi_flags = 0
      endif
!
      if (self%met_opt /= 1) then
         Allocate (self%surf_alb_uv      (i1:i2, ju1:j2))
         Allocate (self%del_surf_alb_uv  (i1:i2, ju1:j2))
         self%surf_alb_uv   = 0.0d0
         self%del_surf_alb_uv   = 0.0d0
      end if
!
!     =================
      if (self%met_opt /= 3) then
!     =================
!
        if (chem_opt == 2) then
          Allocate (self%humidity(i1:i2, ju1:j2, k1:k2))
          self%humidity = 0.0d0
        end if
!
!     ====
      else
!     ====
!
        if (doMEGANemission) then
           Allocate (self%surfTemp15DayAvg (i1:i2, ju1:j2))
           self%surfTemp15DayAvg = 0.0d0
        end if
!
        Allocate (self%con_precip       (i1:i2, ju1:j2))
        Allocate (self%del_con_precip   (i1:i2, ju1:j2))
        Allocate (self%grnd_temp        (i1:i2, ju1:j2))
        Allocate (self%del_grnd_temp    (i1:i2, ju1:j2))
        self%con_precip = 0.0d0
        self%del_con_precip = 0.0d0
        self%grnd_temp  = 0.0d0
        self%del_grnd_temp  = 0.0d0
!
        Allocate (self%pbl              (i1:i2, ju1:j2))
        Allocate (self%del_pbl          (i1:i2, ju1:j2))
        Allocate (self%radswg           (i1:i2, ju1:j2))
        Allocate (self%del_radswg       (i1:i2, ju1:j2))
        self%pbl    = 0.0d0
        self%del_pbl    = 0.0d0
        self%radswg = 0.0d0
        self%del_radswg = 0.0d0
!
        Allocate (self%saldif           (i1:i2, ju1:j2))
        Allocate (self%del_saldif       (i1:i2, ju1:j2))
        Allocate (self%saldir           (i1:i2, ju1:j2))
        Allocate (self%del_saldir       (i1:i2, ju1:j2))
        self%saldif = 0.0d0
        self%del_saldif = 0.0d0
        self%saldir = 0.0d0
        self%del_saldir = 0.0d0
!
        Allocate (self%sasdif           (i1:i2, ju1:j2))
        Allocate (self%del_sasdif       (i1:i2, ju1:j2))
        Allocate (self%sasdir           (i1:i2, ju1:j2))
        Allocate (self%del_sasdir       (i1:i2, ju1:j2))
        self%sasdif = 0.0d0
        self%del_sasdif = 0.0d0
        self%sasdir = 0.0d0
        self%del_sasdir = 0.0d0
!
        Allocate (self%surf_air_temp    (i1:i2, ju1:j2))
        Allocate (self%del_surf_air_temp(i1:i2, ju1:j2))
        self%surf_air_temp = 0.0d0
        self%del_surf_air_temp = 0.0d0
!
        Allocate (self%surf_rough       (i1:i2, ju1:j2))
        Allocate (self%del_surf_rough   (i1:i2, ju1:j2))
        Allocate (self%tot_precip       (i1:i2, ju1:j2))
        Allocate (self%del_tot_precip   (i1:i2, ju1:j2))
        self%surf_rough = 0.0d0
        self%del_surf_rough = 0.0d0
        self%tot_precip = 0.0d0
        self%del_tot_precip = 0.0d0
!
        Allocate (self%ustar            (i1:i2, ju1:j2))
        Allocate (self%del_ustar        (i1:i2, ju1:j2))
        self%ustar = 0.0d0
        self%del_ustar = 0.0d0
!
!
        Allocate (self%humidity     (i1:i2, ju1:j2, k1:k2))
        Allocate (self%del_humidity (i1:i2, ju1:j2, k1:k2))
        Allocate (self%kzz          (i1:i2, ju1:j2, k1:k2))
        Allocate (self%del_kzz      (i1:i2, ju1:j2, k1:k2))
        self%humidity = 0.0d0
        self%del_humidity = 0.0d0
        self%kzz      = 0.0d0
        self%del_kzz      = 0.0d0
!
        Allocate (self%max_cloud    (i1:i2, ju1:j2, k1:k2))
        Allocate (self%del_max_cloud(i1:i2, ju1:j2, k1:k2))
        Allocate (self%ran_cloud    (i1:i2, ju1:j2, k1:k2))
        Allocate (self%del_ran_cloud(i1:i2, ju1:j2, k1:k2))
        self%max_cloud = 0.0d0
        self%del_max_cloud = 0.0d0
        self%ran_cloud = 0.0d0
        self%del_ran_cloud = 0.0d0
!
        Allocate (self%tau_cloud    (i1:i2, ju1:j2, k1:k2))
        Allocate (self%del_tau_cloud(i1:i2, ju1:j2, k1:k2))
        self%tau_cloud = 0.0d0
        self%del_tau_cloud = 0.0d0
        if ((chem_opt ==2) .or. (chem_opt ==7) .or. (chem_opt == 8)) then
!.sds        if ((chem_opt ==7) .or. (chem_opt == 8)) then
           allocate(self%totalCloudFraction(i1:i2, ju1:j2, k1:k2))
           self%totalCloudFraction = 0.0d0
        end if
!
        Allocate (self%fracCloudCover (i1:i2, ju1:j2))
        self%fracCloudCover = 0.0d0
!
        Allocate (self%moistq       (i1:i2, ju1:j2, k1:k2))
        self%moistq = 0.0d0
!
        if ((convec_opt /= 0) .or. do_wetdep) then
          Allocate (self%cmf    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_cmf(i1:i2, ju1:j2, k1:k2))
          self%cmf = 0.0d0
          self%del_cmf = 0.0d0
        endif
!
        if (convec_opt /= 0 .or. do_drydep) then
          Allocate (self%dtrn    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_dtrn(i1:i2, ju1:j2, k1:k2))
          self%dtrn = 0.0d0
          self%del_dtrn = 0.0d0
        endif
!
!... read in convective fluxes
        if (convec_opt == 2) then
          Allocate (self%ed    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_ed(i1:i2, ju1:j2, k1:k2))
          Allocate (self%eu    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_eu(i1:i2, ju1:j2, k1:k2))
          Allocate (self%md    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_md(i1:i2, ju1:j2, k1:k2))
          self%ed = 0.0d0
          self%del_ed = 0.0d0
          self%eu = 0.0d0
          self%del_eu = 0.0d0
          self%md = 0.0d0
          self%del_md = 0.0d0
        endif
!
        if (convec_opt == 3) then
          Allocate (self%zmdu    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_zmdu(i1:i2, ju1:j2, k1:k2))
          Allocate (self%zmeu    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_zmeu(i1:i2, ju1:j2, k1:k2))
          Allocate (self%zmed    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_zmed(i1:i2, ju1:j2, k1:k2))
          Allocate (self%zmmd    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_zmmd(i1:i2, ju1:j2, k1:k2))
          Allocate (self%zmmu    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_zmmu(i1:i2, ju1:j2, k1:k2))
          Allocate (self%hkdu    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_hkdu(i1:i2, ju1:j2, k1:k2))
          Allocate (self%hkeu    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_hkeu(i1:i2, ju1:j2, k1:k2))
          Allocate (self%hkmu    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_hkmu(i1:i2, ju1:j2, k1:k2))
          self%zmdu = 0.0d0
          self%del_zmdu = 0.0d0
          self%zmeu = 0.0d0
          self%del_zmeu = 0.0d0
          self%zmed = 0.0d0
          self%del_zmed = 0.0d0
          self%zmmd = 0.0d0
          self%del_zmmd = 0.0d0
          self%zmmu = 0.0d0
          self%del_zmmu = 0.0d0
          self%hkdu = 0.0d0
          self%del_hkdu = 0.0d0
          self%hkeu = 0.0d0
          self%del_hkeu = 0.0d0
          self%hkmu = 0.0d0
          self%del_hkmu = 0.0d0
        endif
!
        if (self%metdata_name_model(1:5) == 'GEOS4' .or. &
     &      self%metdata_name_model(1:5) == 'GEOS5') then
          Allocate (self%u10m    (i1:i2, ju1:j2))
          Allocate (self%del_u10m(i1:i2, ju1:j2))
          Allocate (self%v10m    (i1:i2, ju1:j2))
          Allocate (self%del_v10m(i1:i2, ju1:j2))
          Allocate (self%gwet    (i1:i2, ju1:j2))
          Allocate (self%del_gwet(i1:i2, ju1:j2))
          Allocate (self%pardif    (i1:i2, ju1:j2))
          Allocate (self%del_pardif(i1:i2, ju1:j2))
          Allocate (self%pardir    (i1:i2, ju1:j2))
          Allocate (self%del_pardir(i1:i2, ju1:j2))
          self%u10m     = 0.0d0
          self%v10m     = 0.0d0
          self%gwet     = 0.0d0
          self%pardif     = 0.0d0
          self%pardir     = 0.0d0
          self%del_u10m = 0.0d0
          self%del_v10m = 0.0d0
          self%del_gwet = 0.0d0
          self%del_pardif = 0.0d0
          self%del_pardir = 0.0d0
        end if
!
        if ((self%metdata_name_org(1:4) == 'GISS') .or.  &
     &      ((self%metdata_name_org  (1:4) == 'NCAR') .and.  &
     &       (self%metdata_name_model(1:4) == 'CCM3'))) then
          Allocate (self%clwc    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_clwc(i1:i2, ju1:j2, k1:k2))
          self%clwc = 0.0d0
          self%del_clwc = 0.0d0
        endif
!
        if (self%metdata_name_org(1:2) == 'EC') then
          Allocate (self%u10m    (i1:i2, ju1:j2))
          self%u10m     = 0.0d0
          Allocate (self%del_u10m(i1:i2, ju1:j2))
          self%del_u10m = 0.0d0
          Allocate (self%v10m    (i1:i2, ju1:j2))
          self%v10m     = 0.0d0
          Allocate (self%del_v10m(i1:i2, ju1:j2))
          self%del_v10m = 0.0d0
          Allocate (self%rain      (i1:i2, ju1:j2, k1:k2))
          self%rain = 0.0d0
          Allocate (self%del_rain  (i1:i2, ju1:j2, k1:k2))
          self%del_rain = 0.0d0
          Allocate (self%clwc    (i1:i2, ju1:j2, k1:k2))
          self%clwc = 0.0d0
          Allocate (self%del_clwc(i1:i2, ju1:j2, k1:k2))
          self%del_clwc = 0.0d0
          Allocate (self%del_moistq(i1:i2, ju1:j2, k1:k2))
          self%del_moistq = 0.0d0
        end if
!
        if ((self%metdata_name_org(1:4) == 'GISS') .or.  &
     &      ((self%metdata_name_org  (1:4) == 'NCAR' ) .and.  &
     &       (self%metdata_name_model(1:5) == 'MATCH'))) then
          Allocate (self%rain      (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_rain  (i1:i2, ju1:j2, k1:k2))
          self%rain = 0.0d0
          self%del_rain = 0.0d0
!
        elseif (self%metdata_name_org(1:4) == 'GMAO' .and.  &
     &         (self%metdata_name_model(1:5) == 'GEOS4' .or. &
     &          self%metdata_name_model(1:5) == 'GEOS5')) then
!
          !---------------
          ! For FastJX 6.5
          !---------------
! Needed for GEOS-4 too? Megan Damon July 20th 2012
!          if (self%metdata_name_model(1:5) == 'GEOS5') then
             Allocate (self%taucli   (i1:i2, ju1:j2, k1:k2))
             Allocate (self%del_taucli (i1:i2, ju1:j2, k1:k2))
             self%taucli = 0.0d0
             self%del_taucli = 0.0d0
!
             Allocate (self%tauclw   (i1:i2, ju1:j2, k1:k2))
             Allocate (self%del_tauclw (i1:i2, ju1:j2, k1:k2))
             self%tauclw = 0.0d0
             self%del_tauclw = 0.0d0
!          end if
          !---------------
!
          Allocate (self%rain_ls   (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_rain_ls (i1:i2, ju1:j2, k1:k2))
          self%rain_ls = 0.0d0
          self%del_rain_ls = 0.0d0
          Allocate (self%rain_zm   (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_rain_zm (i1:i2, ju1:j2, k1:k2))
          self%rain_zm = 0.0d0
          self%del_rain_zm = 0.0d0
          Allocate (self%rain_hk   (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_rain_hk (i1:i2, ju1:j2, k1:k2))
          self%rain_hk = 0.0d0
          self%del_rain_hk = 0.0d0
          Allocate (self%rain      (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_rain    (i1:i2, ju1:j2, k1:k2))
          self%rain = 0.0d0
          self%del_rain = 0.0d0
!... need to store sum of rain prod into moistq in read_met2
          Allocate (self%del_moistq(i1:i2, ju1:j2, k1:k2))
          self%del_moistq = 0.0d0
!
        elseif (self%metdata_name_org(1:2) == 'EC') then
          Allocate (self%rain_ls   (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_rain_ls (i1:i2, ju1:j2, k1:k2))
          self%rain_ls = 0.0d0
          self%del_rain_ls = 0.0d0
          Allocate (self%rain_zm   (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_rain_zm (i1:i2, ju1:j2, k1:k2))
          self%rain_zm = 0.0d0
          self%del_rain_zm = 0.0d0
!
        else
          Allocate (self%del_moistq(i1:i2, ju1:j2, k1:k2))
          self%del_moistq = 0.0d0
        endif
!
!
        if ((self%metdata_name_org  (1:4) == 'NCAR') .and.  &
     &      (self%metdata_name_model(1:4) == 'CCM3')) then
          Allocate (self%aconv      (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_aconv  (i1:i2, ju1:j2, k1:k2))
          Allocate (self%astrat     (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_astrat (i1:i2, ju1:j2, k1:k2))
          self%aconv  = 0.0d0
          self%del_aconv  = 0.0d0
          self%astrat = 0.0d0
          self%del_astrat = 0.0d0
!
          Allocate (self%cldmas     (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_cldmas (i1:i2, ju1:j2, k1:k2))
          Allocate (self%omega      (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_omega  (i1:i2, ju1:j2, k1:k2))
          self%cldmas = 0.0d0
          self%del_cldmas = 0.0d0
          self%omega  = 0.0d0
          self%del_omega  = 0.0d0
!
          Allocate (self%relhum     (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_relhum (i1:i2, ju1:j2, k1:k2))
          Allocate (self%rhclear    (i1:i2, ju1:j2, k1:k2))
          Allocate (self%del_rhclear(i1:i2, ju1:j2, k1:k2))
          self%relhum  = 0.0d0
          self%del_relhum  = 0.0d0
          self%rhclear = 0.0d0
          self%del_rhclear = 0.0d0
        end if
!
!     ======
      end if
!     ======
!
      !----------------------------
      ! Allocate derived quantities
      !----------------------------
!
      Allocate (self%tropopausePress(i1:i2, ju1:j2))
      self%tropopausePress = 0.0d0
!
      Allocate (self%potentialVorticity(i1:i2, ju1:j2, k1:k2))
      self%potentialVorticity = 0.0d0
!
      Allocate (self%potentialTemp(i1:i2, ju1:j2, k1:k2))
      self%potentialTemp = 0.0d0
!
      Allocate (self%mass(i1:i2, ju1:j2, k1:k2))
      self%mass = 0.0d0
!
      Allocate (self%gridBoxHeight(i1:i2, ju1:j2, k1:k2))
      self%gridBoxHeight = 0.0d0
!
      Allocate (self%relativeHumidity(i1:i2, ju1:j2, k1:k2))
      self%relativeHumidity = 0.0d0
!
      Allocate (self%press3c(ilo:ihi, julo:jhi, k1:k2))
      Allocate (self%press3e(ilo:ihi, julo:jhi, k1-1:k2))
      self%press3c = 0.0d0
      self%press3e = 0.0d0
!
      return
!
      end subroutine allocateMetFields
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readMetFieldsResourceFile
!
! !INTERFACE:
!
      subroutine readMetFieldsResourceFile (self, gmiDomain, Diagnostics,&
                    gmiClock, config)
!
! !USES:
      use GmiCheckNamelistFile_mod, only : CheckNamelistOptionRange
      use GmiMetFieldsAttribute_mod, only : getMetFieldsAttribute
!
      implicit none
!
! !INPUT PARAMETERS:
      type (t_gmiDomain), intent(in) :: gmiDomain
      type (t_gmiClock), intent(in) :: gmiClock
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config), intent(inOut) :: config
      type (t_MetFields), intent(inOut) :: self
!
! !DESCRIPTION:
! Reads the MetFields section of the resource file.
!
! !LOCAL VARIABLES:
      integer :: ncid_met, STATUS, RC, ic
      integer :: max_list_items, totNumberDays, begGmiDate
      logical :: doMEGANemission
      real*8  :: tmet1
      character(len=ESMF_MAXSTR) :: err_msg, IAm
!EOP
!-----------------------------------------------------------------------------
!BOC
      Iam = "readMetFieldsResourceFile"
!
      if (pr_diag) Write(6,*) IAm, 'called by ', procID
!
      !################################
      ! Begin reading the resource file
      !################################
!
!     ----------------------------------------------------------
!     met_opt
!       1:  use values fixed in code for        u, v, ps, & kel;
!           no other met data set
!       2:  read in a minimal set of met data:  u, v, ps, & kel;
!           no other met data set
!       3:  read in a full set of met data
!     ----------------------------------------------------------
!
      call ESMF_ConfigGetAttribute(config, self%met_opt, &
     &                label   = "met_opt:",&
     &                default = 3, rc=STATUS )
      VERIFY_(STATUS)
!
      ! 'A'  ! 'A' => A grid, 'C' => C grid
!
      call ESMF_ConfigGetAttribute(config, self%met_grid_type, &
     &                label   = "met_grid_type:",&
     &                default = 'A', rc=STATUS )
      VERIFY_(STATUS)
!
      ! time increment for reading new met data must be a multiple of tdt (s)
!
      call ESMF_ConfigGetAttribute(config, self%mdt, &
     &                label   = "mdt:",&
     &                default = 21600.0d0, rc=STATUS )
      VERIFY_(STATUS)
!
      call rcEsmfReadLogical(config, self%do_cycle_met, "do_cycle_met:", &
     &              default=.false., rc=STATUS )
      VERIFY_(STATUS)
!
!     ----------------------------------------------------------------
!     do_timinterp_met   : time interpolate met  fields, except winds?
!     do_timinterp_winds : time interpolate wind fields?
!     Note that pressure fields are always interpolated.
!     ----------------------------------------------------------------
!
      call rcEsmfReadLogical(config, self%do_timinterp_met, &
                "do_timinterp_met:", default=.true., rc=STATUS )
      VERIFY_(STATUS)
!
      call rcEsmfReadLogical(config, self%do_timinterp_winds, &
     &               "do_timinterp_winds:", default=.true., rc=STATUS )
      VERIFY_(STATUS)
!
      call rcEsmfReadLogical(config, self%do_wind_pole, "do_wind_pole:", &
     &              default=.false., rc=STATUS )
      VERIFY_(STATUS)
!
      ! index of NetCDF file name to start reading met input data from
!
      call ESMF_ConfigGetAttribute(config, self%met_infile_num, &
     &                label   = "met_infile_num:",&
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)
!
      ! NetCDF file record to start reading met input data from
!
      call ESMF_ConfigGetAttribute(config, self%mrnum_in, &
     &                label   = "mrnum_in:",&
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)
!
      ! time tag of mrnum_in (s)
!
      call ESMF_ConfigGetAttribute(config, self%tmet1, &
     &                label   = "tmet1:",&
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)
!
      call rcEsmfReadLogical(config, self%do_read_met_list, &
                 "do_read_met_list:", default=.false., rc=STATUS )
      VERIFY_(STATUS)
!
      if (self%do_read_met_list) then
         call rcEsmfReadLogical(config, self%do_manualReadMetList, &
                     "do_manualReadMetList:", default=.true., rc=STATUS )
         VERIFY_(STATUS)
!
         if (self%do_manualReadMetList) then
            call ESMF_ConfigGetAttribute(config, self%met_filnam_list, &
     &                label   = "met_filnam_list:",&
     &                default = 'met_filnam_list.in', rc=STATUS )
            VERIFY_(STATUS)
         else
            call Get_begGmiDate(gmiClock, begGmiDate)
!
            call rcEsmfReadLogical(config, doMEGANemission, &
                        "doMEGANemission:", default=.false., rc=STATUS )
            VERIFY_(STATUS)
!
            call ESMF_ConfigGetAttribute(config, self%metStartDate, &
                            label  = "metStartDate:",&
                           default = begGmiDate, rc=STATUS )
            VERIFY_(STATUS)
!
            call ESMF_ConfigGetAttribute(config, self%metRefTime, &
                            label  = "metRefTime:",&
                           default = 120000, rc=STATUS )
            VERIFY_(STATUS)
!
            call ESMF_ConfigGetAttribute(config, self%metFileTemplate, &
                            label  = "metFileTemplate:",&
                           default = " ", rc=STATUS )
            VERIFY_(STATUS)
         end if
      else
         self%met_infile_names(:) = CHAR_DUM_VALUE
         call rcEsmfReadTable (config, self%met_infile_names, &
     &              "met_infile_names::", rc=STATUS)
         VERIFY_(STATUS)
      end if
!
      ! -----------------------------------------------------
      ! gwet_opt = 0:
      !            1:
      ! -----------------------------------------------------
!
      call ESMF_ConfigGetAttribute(config, self%gwet_opt, &
     &                label   = "gwet_opt:",&
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
      !##############################
      ! End reading the resource file
      !##############################
!
      ! ---------------------------------------------------------------
      ! Check option ranges.  Note that as new options are added, these
      ! range checks will have to be modified.
      ! ---------------------------------------------------------------
!
      call CheckNamelistOptionRange ('met_opt', self%met_opt, 1, 3)
      call CheckNamelistOptionRange ('gwet_opt', self%gwet_opt, 0, 1)
!
      self%met1_infile_num = self%met_infile_num
      self%met2_infile_num = self%met_infile_num
!
      self%metNumMEGAN     = self%met_infile_num - 1 ! for MEGAN emission
!
      self%m1rnum_in = self%mrnum_in
      self%m2rnum_in = self%mrnum_in
!
      self%tmet1_1 = self%tmet1
      self%tmet2_1 = self%tmet1
!
      self%tmet1_2 = 0.0d0
      self%tmet2_2 = 0.0d0
!
      self%met_num_infiles = 0
      max_list_items  = MAX_INFILE_NAMES
!
      if (self%do_read_met_list) then
         IF (self%do_manualReadMetList) THEN
            call Read_List (self%met_infile_names, self%met_num_infiles, &
                         self%met_filnam_list,  max_list_items)
         ELSE
            !###############################################
            ! Automatically generate the metFields file list
            !###############################################
            call Get_totNumDays(GmiClock, totNumberDays)
!
            call createMetFileList(self%met_infile_names, &
                       self%metFileTemplate, self%metStartDate, &
                       self%metRefTime, totNumberDays, doMEGANemission, &
                       self%met_num_infiles)
         ENDIF
      else
         aloop: do ic = 1, MAX_INFILE_NAMES
            if (self%met_infile_names(ic)(1:3) /= CHAR_DUM_VALUE) then
               self%met_num_infiles = ic
            else
               exit aloop
            end if
         end do aloop
      end if
!
      call Ncop_Rd (ncid_met, self%met_infile_names(self%met_infile_num))
!
      call getMetFieldsAttribute &
     &    (ncid_met, self%metdata_name, self%metdata_name_org,  &
     &     self%metdata_name_model, self%metdata_name_dims, pr_diag, procID)
!
      call Nccl (ncid_met)
!
     !------------------------
     ! Checking for mismatches
     !------------------------
!
      if (self%mdt <= 0.0d0) then
         err_msg = 'mdt <= 0.0 in the resource file.'
         call GmiPrintError (err_msg, .true., 0, 0, 0, 1, self%mdt, 0.0d0)
      end if
!
      if ((self%met_grid_type /= 'A') .and. (self%met_grid_type /= 'C')) then
        err_msg =  &
     &    'met_grid_type problem in Check_Nlvalue:  ' // self%met_grid_type
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if
!
      if (do_drydep .and. (self%met_opt /= 3)) then
        err_msg = 'do_drydep/met_opt problem in Check_Nlvalue.'
        call GmiPrintError  &
     &    (err_msg, .true., 1, self%met_opt, 0, 0, 0.0d0, 0.0d0)
      end if
!
      if (do_wetdep .and. (self%met_opt /= 3)) then
        err_msg = 'do_wetdep/met_opt problem in Check_Nlvalue.'
        call GmiPrintError  &
     &    (err_msg, .true., 1, self%met_opt, 0, 0, 0.0d0, 0.0d0)
      end if
!
      return
!
      end subroutine readMetFieldsResourceFile
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Control_Met1_File_Input
!
! !INTERFACE:
!
      subroutine Control_Met1_File_Input (self, out_of_data)
!
      implicit none
!
! !OUTPUT PARAMETERS:
      logical, intent(out) :: out_of_data ! out of data?
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_MetFields), intent(inOut) :: self
!
! !DESCRIPTION:
! Controls the met1 (ps, u, v, & kel) file input.
!EOP
!-----------------------------------------------------------------------------
!BOC
!
      out_of_data = .false.
!
      if (self%m1rnum_in > self%m1num_recs) then
!
        out_of_data = .false.
!
!       ======================
        call Get_Next_Met_File  &
!       ======================
     &    (self%do_cycle_met, out_of_data, self%met1_infile_num,  &
     &     self%met_num_infiles, self%m1num_recs, self%m1rnum_in, &
     &     self%ncid_met1, self%met_infile_names, pr_diag, procID)
!
!       =======================
        if (out_of_data) return  ! Early return.
!       =======================
!
      end if
!
      self%tmet1_1 = self%tmet1_2
!
!     ==============
      call Save_Met1 (self)
!     ==============
!
      call Read_Met1 (self%metdata_name_model, &
    &           self%ncid_met1, self%m1rnum_in, self%del_psxGlob, &
    &           self%del_uuxGlob, self%del_vvxGlob, self%del_kelGlob, &
    &           pr_diag, procID, i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, k1, k2, &
    &           ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl)
!
      call Check_Met1_Range (self%metdata_name_org, self%del_psxGlob, &
     &          self%del_uuxGlob, self%del_vvxGlob, self%del_kelGlob, &
     &          procID, i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, ilo_gl,  &
     &          ihi_gl, julo_gl, jvlo_gl, jhi_gl, k1, k2)
!
      self%m1rnum_in = self%m1rnum_in + 1
!
      self%tmet1_2   = Nint (self%tmet1_1 + self%mdt)
!
!     ------------------------------------------------
!     Always done on Master, so only "wrap" necessary.
!     ------------------------------------------------
!
      call wrapMaster_2d  (self%psxGlob, i1_gl, i2_gl, ju1_gl, j2_gl, ilo_gl, &
     &         ihi_gl, julo_gl, jhi_gl, gmi_nborder)
      call wrapMaster_2d  (self%del_psxGlob, i1_gl, i2_gl, ju1_gl, j2_gl, ilo_gl, &
     &         ihi_gl, julo_gl, jhi_gl, gmi_nborder)
      call wrapMaster_3du (self%kelGlob, i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2,     &
     &         ilo_gl, ihi_gl, julo_gl, jhi_gl, gmi_nborder)
      call wrapMaster_3du (self%del_kelGlob, i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2, &
     &         ilo_gl, ihi_gl, julo_gl, jhi_gl, gmi_nborder)
      call wrapMaster_3du (self%uuxGlob, i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2,     &
     &         ilo_gl, ihi_gl, julo_gl, jhi_gl, gmi_nborder)
      call wrapMaster_3du (self%del_uuxGlob, i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2, &
     &         ilo_gl, ihi_gl, julo_gl, jhi_gl, gmi_nborder)
      call wrapMaster_3dv (self%vvxGlob, i1_gl, i2_gl, jv1_gl, j2_gl, k1, k2, &
     &         ilo_gl, ihi_gl, jvlo_gl, jhi_gl, gmi_nborder)
      call wrapMaster_3dv (self%del_vvxGlob, i1_gl, i2_gl, jv1_gl, j2_gl, k1, k2, &
     &         ilo_gl, ihi_gl, jvlo_gl, jhi_gl, gmi_nborder)
!
      return
!
      end subroutine Control_Met1_File_Input
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Met1
!
! !INTERFACE:
!
      subroutine Set_Met1 (self, gmiDomain)
!
      implicit none
!
! !INPUT PARAMETERS:
      type (t_gmiDomain), intent(in) :: gmiDomain
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_metFields), intent(inOut) :: self
!
! !DESCRIPTION:
! Sets the met1 data (ps, u, v, & kel).
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      logical :: out_of_data
      real*8 , allocatable :: dlatr (:)
      real*8 , allocatable :: coscen(:)
!EOC
!-----------------------------------------------------------------------------
!BOP
      if (pr_diag) Write (6,*) 'Set_Met1 called by ', procID
!
!     =================
      if (self%met_opt == 1) then
!     =================
!
        allocate(dlatr (ju1_gl:j2_gl))
        call Get_dlatr(gmiDomain, dlatr)
!
        allocate(coscen(ju1_gl:j2_gl))
        call Get_coscen(gmiDomain, coscen)
!
        call Set_Met1_Simple (self%psxGlob, self%del_psxGlob, self%uuxGlob, &
     &           self%del_uuxGlob, self%vvxGlob, self%del_vvxGlob, &
     &           self%kelGlob, self%del_kelGlob, coscen, dlatr, &
     &           self%do_wind_pole, pr_diag, procID, i1_gl, i2_gl, ju1_gl, &
     &           jv1_gl, j2_gl, ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl, k1, k2)
!
        self%tmet1_2 = self%tmet1_1 + 1.0d99
!
!     ====
      else
!     ====
!
        call Ncop_Rd (self%ncid_met1, self%met_infile_names(self%met1_infile_num))
!
        call Ncget_Unlim_Dimlen (self%ncid_met1, self%m1num_recs)
!
        if (self%m1rnum_in > self%m1num_recs) then
          err_msg = 'm1rnum_in/m1num_recs problem in Set_Met1.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%m1rnum_in, self%m1num_recs, 0, 0.0d0, 0.0d0)
        end if
!
        call Read_Met1 (self%metdata_name_model, &
     &           self%ncid_met1, self%m1rnum_in, self%psxGlob, &
     &           self%uuxGlob, self%vvxGlob, self%kelGlob, pr_diag, procID, &
     &           i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, k1, k2, ilo_gl, ihi_gl, &
     &           julo_gl, jvlo_gl, jhi_gl)
!
        call Check_Met1_Range (self%metdata_name_org, self%psxGlob, &
     &           self%uuxGlob, self%vvxGlob, self%kelGlob, procID, &
     &           i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, &
     &           ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl, k1, k2)
!
        self%m1rnum_in = self%m1rnum_in + 1
!
        if (self%m1rnum_in > self%m1num_recs) then
!
          out_of_data = .false.
!
!         ======================
          call Get_Next_Met_File  &
!         ======================
     &      (self%do_cycle_met, out_of_data, self%met1_infile_num,  &
     &       self%met_num_infiles, self%m1num_recs, self%m1rnum_in, &
     &       self%ncid_met1, self%met_infile_names, pr_diag, procID)
!
          if (out_of_data) then
            err_msg = 'out_of_data problem in Set_Met1.'
            call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
          end if
!
        end if
!
        call Read_Met1 (self%metdata_name_model, &
     &           self%ncid_met1, self%m1rnum_in, self%del_psxGlob, &
     &           self%del_uuxGlob, self%del_vvxGlob, self%del_kelGlob, &
     &           pr_diag, procID, i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, &
     &           k1, k2, ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl)
!
        call Check_Met1_Range (self%metdata_name_org, self%del_psxGlob, &
     &             self%del_uuxGlob, self%del_vvxGlob, self%del_kelGlob, &
     &             procID, i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, ilo_gl, &
     &             ihi_gl, julo_gl, jvlo_gl, jhi_gl, k1, k2)
!
        self%m1rnum_in = self%m1rnum_in + 1
        self%tmet1_2   = Nint (self%tmet1_1 + self%mdt)
!
      end if
!
!     ------------------------------------------------
!     Always done on Master, so only "wrap" necessary.
!     ------------------------------------------------
!
      call wrapMaster_2d  (self%psxGlob, i1_gl, i2_gl, ju1_gl, j2_gl, ilo_gl, &
     &         ihi_gl, julo_gl, jhi_gl, gmi_nborder)
      call wrapMaster_2d  (self%del_psxGlob, i1_gl, i2_gl, ju1_gl, j2_gl, ilo_gl, &
     &         ihi_gl, julo_gl, jhi_gl, gmi_nborder)
      call wrapMaster_3du (self%kelGlob, i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2,     &
     &         ilo_gl, ihi_gl, julo_gl, jhi_gl, gmi_nborder)
      call wrapMaster_3du (self%del_kelGlob, i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2, &
     &         ilo_gl, ihi_gl, julo_gl, jhi_gl, gmi_nborder)
      call wrapMaster_3du (self%uuxGlob, i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2,     &
     &         ilo_gl, ihi_gl, julo_gl, jhi_gl, gmi_nborder)
      call wrapMaster_3du (self%del_uuxGlob, i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2, &
     &         ilo_gl, ihi_gl, julo_gl, jhi_gl, gmi_nborder)
      call wrapMaster_3dv (self%vvxGlob, i1_gl, i2_gl, jv1_gl, j2_gl, k1, k2, &
     &         ilo_gl, ihi_gl, jvlo_gl, jhi_gl, gmi_nborder)
      call wrapMaster_3dv (self%del_vvxGlob, i1_gl, i2_gl, jv1_gl, j2_gl, k1, k2, &
     &         ilo_gl, ihi_gl, jvlo_gl, jhi_gl, gmi_nborder)
!
      return
!
      end subroutine Set_Met1
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Save_Met1
!
! !INTERFACE:
!
      subroutine Save_Met1 (self)
!
      implicit none
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_metFields), intent(inOut) :: self
!
! !DESCRIPTION:
! Saves the met1 data (ps, u, v, & kel) for the current step
! before a new set of met1 data is read in.
!EOP
!-----------------------------------------------------------------------------
!BOC
      self%psxGlob(:,:)   = self%del_psxGlob(:,:)
      self%kelGlob(:,:,:) = self%del_kelGlob(:,:,:)
      self%uuxGlob(:,:,:) = self%del_uuxGlob(:,:,:)
      self%vvxGlob(:,:,:) = self%del_vvxGlob(:,:,:)
!
      return
!
      end subroutine Save_Met1
!EOC
!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update_Met1
!
! !INTERFACE:
!
      subroutine Update_Met1 (self, gmiDomain, new_met_rec, &
     &                j1p, j2p, rd_restart, num_adv_time_steps, do_var_adv_tstp, &
     &                advec_consrv_opt, pmet2_opt, press_fix_opt)
!
! !USES:
      use GmiPressureFixer_mod , only : adjustPressFixer
      use GmiUtilsMetFields_mod, only : Check_Vert_Courant
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: rd_restart, do_var_adv_tstp
      integer, intent(in) :: j1p, j2p
      integer, intent(in) :: advec_consrv_opt, pmet2_opt, press_fix_opt
      logical, intent(in) :: new_met_rec
      type (t_gmiDomain), intent(in) :: gmiDomain
!
! !OUTPUT PARAMETERS:
      integer, intent(out) :: num_adv_time_steps
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_metFields), intent(inOut) :: self
!
! !DESCRIPTION:
!  Updates the met1 data (ps, u, v, & kel) each time step.
!
! !LOCAL VARIABLES:
      logical, save :: first = .true.
      integer       :: ilat, ilong
      real*8  :: alpha
      real*8  :: weight1, weight2
      real*8,  save :: rstep          = 0.0d0
      real*8,  save :: tot_met1_steps = 0.0d0
!
!     ---------------------------------------
!     pmet1 : surface pressure at t1     (mb)
!     pmet2 : surface pressure at t1+tdt (mb)
!     ---------------------------------------
!
      real*8, allocatable, save :: pmet1(:,:)
      real*8, allocatable, save :: pmet2(:,:)
!
!      real*8 , allocatable :: dlatr (:)
!      real*8 , allocatable :: coscen(:)
      real*8 , allocatable :: cose  (:)
      real*8 , allocatable :: cosp  (:)
      real*8 , allocatable :: rel_area(:,:)
      real*8 , allocatable :: geofac(:)
      real*8               :: geofac_pc
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6, *) 'Update_Met1 called by ', procID
!
!      allocate(dlatr (ju1_gl:j2_gl))
!      call Get_dlatr(gmiDomain, dlatr)
!
!      allocate(coscen(ju1_gl:j2_gl))
!      call Get_coscen(gmiDomain, coscen)
!
      allocate(cose  (ju1_gl:j2_gl+1))
      call Get_cose(gmiDomain, cose)
!
      allocate(cosp  (ju1_gl:j2_gl))
      call Get_cosp(gmiDomain, cosp)
!
      allocate(rel_area(i1_gl:i2_gl,ju1_gl:j2_gl))
      call Get_rel_areaGlob(gmiDomain, rel_area)
!
      allocate(geofac  (ju1_gl:j2_gl))
      call Get_geofac   (gmiDomain, geofac)
      call Get_geofac_pc(gmiDomain, geofac_pc)
!
      ilong = i2_gl - i1_gl + 1
      ilat  = j2_gl - ju1_gl + 1
!
      if (first) then
        Allocate (pmet1(ilo_gl:ihi_gl, julo_gl:jhi_gl))
        Allocate (pmet2(ilo_gl:ihi_gl, julo_gl:jhi_gl))
        pmet1 = 0.0d0; pmet2 = 0.0d0
      end if
!
!
!     -------------
!     Update pmet1.
!     -------------
!
      if (new_met_rec) then
         pmet1(:,:) = self%psxGlob(:,:)
      else
         pmet1(:,:) = pmet2(:,:)
      end if
!
!     -------------
!     Update pmet2.
!     -------------
!
      if (self%met_opt == 1) then
!
        pmet2(:,:) = self%psxGlob(:,:)
!
      else
!
        if (new_met_rec) then
          rstep = 1.0d0
        else
          rstep = rstep + 1.0d0
        end if
!
        alpha = rstep * (tdt / self%mdt)
!
        pmet2(:,:) = self%psxGlob(:,:) +  &
     &               alpha * (self%del_psxGlob(:,:) - self%psxGlob(:,:))
!
      end if
!
!     -------------
!     Update pctm1.
!     -------------
!
      if ((.not. first) .or. (rd_restart .and. (self%pctm2Glob(1,1) > 0.0d0))) then
         self%pctm1Glob(:,:) = self%pctm2Glob(:,:)
      else
         self%pctm1Glob(:,:) = pmet1(:,:)
      end if
!
      if (first) first = .false.
!
!
!     =============================================
      if (self%do_timinterp_met .or. self%do_timinterp_winds) then
!     =============================================
!
!       ------------------------
!       Update kel &/or uux/vvx.
!       ------------------------
!
        if (new_met_rec) then
          tot_met1_steps = self%mdt / tdt
!
          weight1 = (tot_met1_steps - 0.5d0) / tot_met1_steps
          weight2 =                   0.5d0  / tot_met1_steps
        else
          weight1 = (tot_met1_steps - 0.5d0) / (tot_met1_steps + 0.5d0)
          weight2 =                   1.0d0  / (tot_met1_steps + 0.5d0)
        end if
!
        if (self%do_timinterp_met) then
          self%kelGlob(:,:,:) = weight1 * self%kelGlob    (:,:,:) +  &
     &                          weight2 * self%del_kelGlob(:,:,:)
        end if
!
        if (self%do_timinterp_winds) then
          self%uuxGlob(:,:,:) = weight1 * self%uuxGlob(:,:,:) +  &
     &                          weight2 * self%del_uuxGlob(:,:,:)
          self%vvxGlob(:,:,:) = weight1 * self%vvxGlob(:,:,:) +  &
     &                          weight2 * self%del_vvxGlob(:,:,:)
        end if
!
        tot_met1_steps = tot_met1_steps - 1.0d0
!
!     ======
      end if
!     ======
!
!     =================
      call adjustPressFixer  &
!     =================
     &  (self%metdata_name_org, self%met_grid_type, self%do_timinterp_winds,  &
     &   new_met_rec, advec_consrv_opt, pmet2_opt, press_fix_opt,  &
     &   tdt, geofac_pc, geofac, cose, cosp, rel_area, self%dap, self%dbk,  &
     &   self%pctm1Glob, self%pctm2Glob, pmet2, self%uuxGlob, self%vvxGlob, &
     &   self%xmassGlob, self%ymassGlob, self%zmassGlob, &
     &   pr_diag, procID, numLonDomains, ilong, ilat, j1p, j2p, &
     &   i1_gl, i2_gl, ju1_gl, j2_gl, ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl, i1_gl, &
     &   i2_gl, ju1_gl, j2_gl, k1, k2, commu_npole, commu_spole, gmi_nborder)
!
!
      if (do_var_adv_tstp) then
!
!       ======================
        call Calc_Var_Adv_Tstp  &
!       ======================
     &    (num_adv_time_steps, cose, self%dap, self%dbk, geofac_pc, geofac,  &
     &     self%pctm1Glob, self%pctm2Glob, self%xmassGlob, self%ymassGlob, &
     &     self%zmassGlob,  &
     &     commu_npole, commu_spole, pr_diag, procID, numLonDomains, &
     &   gmi_nborder, j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, k1, k2, ilo_gl, ihi_gl, julo_gl, jhi_gl)
!
        if (pr_diag) then
          Write (6,*) "num_adv_time_steps =", num_adv_time_steps
        end if
!
!
!
        if (num_adv_time_steps > 10) then
          Write (6,*) "WARNING:  num_adv_time_steps is large."
          Write (6,*) "          This may indicate a problem."
          Write (6,*) "          num_adv_time_steps = ",  &
     &                num_adv_time_steps
!c?       call Stopcode ("STOPPING:  num_adv_time_steps > 10.")
        end if
!
!       =======================
        call Check_Vert_Courant  &
!       =======================
     &    (.false., num_adv_time_steps, geofac_pc, geofac,  &
     &     self%dap, self%dbk, self%pctm1Glob, self%pctm2Glob, self%xmassGlob, &
     &     self%ymassGlob, &
     &   pr_diag, procID, numLonDomains, j1p, j2p, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, ilo_gl, ihi_gl, julo_gl, jhi_gl, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, k1, k2, commu_npole, commu_spole)
!
      else
!
        num_adv_time_steps = 1
!
      end if
!
      return
!
      end subroutine Update_Met1
!EOC
!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: met1Glob2Sub
!
! !INTERFACE:
!
     subroutine met1Glob2Sub (self)
!
     implicit none
!
! !INPUT/OUTPUT PARAMETERS:
     type(t_metFields), intent(inOut) :: self
!EOP
!----------------------------------------------------------------------------
!BOC
!
     !---------------------------------
     ! From global domain to subdomains
     !---------------------------------
!
      self%pctm1  (ilo:ihi,julo:jhi)   = self%pctm1Glob  (ilo:ihi,julo:jhi)
      self%pctm2  (ilo:ihi,julo:jhi)   = self%pctm2Glob  (ilo:ihi,julo:jhi)
      self%psx    (ilo:ihi,julo:jhi)   = self%psxGlob    (ilo:ihi,julo:jhi)
      self%del_psx(ilo:ihi,julo:jhi)   = self%del_psxGlob(ilo:ihi,julo:jhi)
      self%uux    (ilo:ihi,julo:jhi,:) = self%uuxGlob    (ilo:ihi,julo:jhi,:)
      self%del_uux(ilo:ihi,julo:jhi,:) = self%del_uuxGlob(ilo:ihi,julo:jhi,:)
      self%vvx    (ilo:ihi,jvlo:jhi,:) = self%vvxGlob    (ilo:ihi,jvlo:jhi,:)
      self%del_vvx(ilo:ihi,jvlo:jhi,:) = self%del_vvxGlob(ilo:ihi,jvlo:jhi,:)
      self%kel    (ilo:ihi,julo:jhi,:) = self%kelGlob    (ilo:ihi,julo:jhi,:)
      self%del_kel(ilo:ihi,julo:jhi,:) = self%del_kelGlob(ilo:ihi,julo:jhi,:)
      self%xmass  (ilo:ihi,julo:jhi,:) = self%xmassGlob  (ilo:ihi,julo:jhi,:)
      self%ymass  (ilo:ihi,julo:jhi,:) = self%ymassGlob  (ilo:ihi,julo:jhi,:)
      self%zmass  (i1:i2  ,ju1:j2  ,:) = self%zmassGlob  (i1:i2  ,ju1:j2  ,:)
!
      return
!
      end subroutine met1Glob2Sub
!EOC
!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Met2
!
! !INTERFACE:
!
      subroutine Set_Met2 (self, gmiDomain)
!
      implicit none
!
! !INPUT PARAMETERS:
      type (t_gmiDomain), intent(in) :: gmiDomain
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_MetFields), intent(inOut) :: self
!
!
! !DESCRIPTION:
!  Sets the met2 data (everything except ps, u, v, & kel).
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      logical :: out_of_data
      real*8, allocatable :: mcor(:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
      if (pr_diag) Write (6,*) 'Set_Met2 called by ', procID
!
      if (self%met_opt == 3) then
!
        call Ncop_Rd (self%ncid_met2, self%met_infile_names(self%met2_infile_num))
!
        call Ncget_Unlim_Dimlen (self%ncid_met2, self%m2num_recs)
!
        if (self%m2rnum_in > self%m2num_recs) then
          err_msg = 'm2rnum_in/m2num_recs problem in Set_Met2.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%m2rnum_in, self%m2num_recs, 0, 0.0d0, 0.0d0)
        end if
!
        allocate(mcor(i1:i2,ju1:j2))
        call Get_mcor(gmiDomain, mcor)
!
!
!       ==============
        call Read_Met2  &
!       ==============
     &    (mcor, self%mass, &
     &     self%metdata_name_org, self%metdata_name_model, do_wetdep, &
     &     self%ncid_met2, self%m2rnum_in, chem_opt, convec_opt, emiss_in_opt, &
     &     self%gwet_opt, sfalbedo_opt, uvalbedo_opt, self%lwi_flags, &
     &     self%con_precip, self%grnd_temp, self%pbl, self%radswg, self%saldif,&
     &     self%saldir, self%sasdif, self%sasdir, self%surf_air_temp, &
     &     self%surf_alb_uv, self%surf_rough, self%tot_precip, self%ustar, &
     &     self%md, self%cmf, self%dtrn, self%ed, self%eu, self%humidity, &
     &     self%kzz, self%max_cloud, self%moistq, self%rain, self%ran_cloud,  &
     &     self%tau_cloud, self%aconv, self%astrat, self%clwc, self%relhum, &
     &     self%rhclear, self%cldmas, self%omega, self%rain_zm, self%rain_hk, &
     &     self%rain_ls, self%taucli, self%tauclw, &
     &     self%zmdu, self%zmeu, self%zmed, self%zmmd, self%zmmu,&
     &     self%hkdu, self%hkeu, self%hkmu, self%u10m, self%v10m, self%gwet, &
     &     self%pardif, self%pardir, pr_diag, procID, i1, i2, ju1, j2, k1, k2, &
     &     i1_gl, ju1_gl, j2_gl)
!
        self%m2rnum_in = self%m2rnum_in + 1
!
        if (self%m2rnum_in > self%m2num_recs) then
!
          out_of_data = .false.
!
!         ======================
          call Get_Next_Met_File  &
!         ======================
     &      (self%do_cycle_met, out_of_data, self%met2_infile_num,  &
     &       self%met_num_infiles, self%m2num_recs, self%m2rnum_in, self%ncid_met2,  &
     &       self%met_infile_names, pr_diag, procID)
!
          if (out_of_data) then
            err_msg = 'out_of_data problem in Set_Met2.'
            call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
          end if
!
        end if
!
!
!       ==============
        call Read_Met2  &
!       ==============
     &    (mcor, self%mass, &
     &     self%metdata_name_org, self%metdata_name_model, do_wetdep, &
     &     self%ncid_met2, self%m2rnum_in, chem_opt, convec_opt, emiss_in_opt, &
     &     self%gwet_opt, sfalbedo_opt, uvalbedo_opt, self%lwi_flags, &
     &     self%del_con_precip, self%del_grnd_temp, self%del_pbl, &
     &     self%del_radswg, self%del_saldif, self%del_saldir, self%del_sasdif, &
     &     self%del_sasdir, self%del_surf_air_temp, self%del_surf_alb_uv,  &
     &     self%del_surf_rough, self%del_tot_precip, self%del_ustar, &
     &     self%del_md, self%del_cmf, self%del_dtrn, self%del_ed, self%del_eu, &
     &     self%del_humidity, self%del_kzz, self%del_max_cloud, &
     &     self%del_moistq, self%del_rain, self%del_ran_cloud,&
     &     self%del_tau_cloud, self%del_aconv, self%del_astrat, self%del_clwc, &
     &     self%del_relhum, self%del_rhclear, self%del_cldmas, self%del_omega, &
     &     self%del_rain_zm, self%del_rain_hk, self%del_rain_ls, &
     &     self%del_taucli, self%del_tauclw, self%del_zmdu,&
     &     self%del_zmeu, self%del_zmed, self%del_zmmd, self%del_zmmu, &
     &     self%del_hkdu, self%del_hkeu, self%del_hkmu, self%del_u10m, &
     &     self%del_v10m, self%del_gwet, self%del_pardif, self%del_pardir, &
     &     pr_diag, procID, i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl, j2_gl)
!
!
        self%m2rnum_in = self%m2rnum_in + 1
!
        self%tmet2_2   = Nint (self%tmet2_1 + self%mdt)
!
      end if
!
      !-----------------------------
      ! Setting of LWI and CMI flags
      !-----------------------------
!
      if (((chem_opt    /= 0) .or.  &
     &    (convec_opt  /= 0) .or. (emiss_in_opt/= 0))) then
         if (.not.((self%metdata_name_org(1:2)   == 'EC') .or.  &
     &             (self%metdata_name_org(1:4)   == 'GMAO') .or.  &
     &             (self%metdata_name_org(1:3)   == 'DAO' .and.  &
     &              self%metdata_name_model(1:5) == 'GEOS3'))) then
            call setLWIflags (self%lwi_flags, i1, i2, ju1, j2, &
     &                          i1_gl, i2_gl, ju1_gl, j2_gl)
         end if
      end if
!
      ! Initialize CMI flags variables
!
      if (lightning_opt == 1) then
         call setCMIflags (self%cmi_flags, i1, i2, ju1, j2, &
     &                       i1_gl, i2_gl, ju1_gl, j2_gl)
      end if
!
      return
!
      end subroutine Set_Met2
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Save_Met2
!
! !INTERFACE:
!
      subroutine Save_Met2  (self)
!
      implicit none
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_MetFields), intent(inOut) :: self
!
! !DESCRIPTION:
! Saves the met2 data (everything except ps, u, v, & kel) for
! the current step before a new set of met2 data is read in.
!
!
!EOP
!-----------------------------------------------------------------------------
!BOC
!
      self%grnd_temp    (:,:) = self%del_grnd_temp    (:,:)
      self%pbl          (:,:) = self%del_pbl          (:,:)
      self%radswg       (:,:) = self%del_radswg       (:,:)
      self%surf_rough   (:,:) = self%del_surf_rough   (:,:)
      self%ustar        (:,:) = self%del_ustar        (:,:)
!
      self%con_precip   (:,:) = self%del_con_precip   (:,:)
      self%tot_precip   (:,:) = self%del_tot_precip   (:,:)
!
      self%surf_air_temp(:,:) = self%del_surf_air_temp(:,:)
!
      if ( self%metdata_name_model(1:5) == 'GEOS4' .or.  &
     &     self%metdata_name_model(1:5) == 'GEOS5') then
         self%gwet(:,:) = self%del_gwet
         self%u10m(:,:) = self%del_u10m
         self%v10m(:,:) = self%del_v10m
         self%pardif(:,:) = self%del_pardif
         self%pardir(:,:) = self%del_pardir
      end if
      if (sfalbedo_opt == 3) then
        self%sasdir     (:,:) = self%del_sasdir       (:,:)
        self%sasdif     (:,:) = self%del_sasdif       (:,:)
        self%saldir     (:,:) = self%del_saldir       (:,:)
        self%saldif     (:,:) = self%del_saldif       (:,:)
      end if
!
      if (uvalbedo_opt == 3) then
        self%surf_alb_uv(:,:) = self%del_surf_alb_uv  (:,:)
      end if
!
      if (convec_opt /= 0 .or. do_wetdep) then
        self%cmf(:,:,:) = self%del_cmf(:,:,:)
        if (convec_opt /= 0 ) then
          self%dtrn(:,:,:) = self%del_dtrn(:,:,:)
        end if
      end if
!
!.sds.. need these arrays for convec_opt=2
      if (convec_opt == 2) then
        self%ed(:,:,:) = self%del_ed(:,:,:)
        self%eu(:,:,:) = self%del_eu(:,:,:)
        self%md(:,:,:) = self%del_md(:,:,:)
      end if
!
!.sds.. need these arrays for convec_opt=3
      if (convec_opt == 3) then
        self%zmdu(:,:,:) = self%del_zmdu(:,:,:)
        self%zmeu(:,:,:) = self%del_zmeu(:,:,:)
        self%zmed(:,:,:) = self%del_zmed(:,:,:)
        self%zmmd(:,:,:) = self%del_zmmd(:,:,:)
        self%zmmu(:,:,:) = self%del_zmmu(:,:,:)
        self%hkdu(:,:,:) = self%del_hkdu(:,:,:)
        self%hkeu(:,:,:) = self%del_hkeu(:,:,:)
        self%hkmu(:,:,:) = self%del_hkmu(:,:,:)
      end if
!
      if ((self%metdata_name_org  (1:4) == 'NCAR'  .and.  &
     &     self%metdata_name_model(1:5) == 'MATCH') .or.  &
     &     self%metdata_name_org(1:4) == 'GISS') then
!
        self%rain   (:,:,:) = self%del_rain  (:,:,:)
!
      elseif (self%metdata_name_model(1:5) == 'GEOS4') then
!... GEOS4 3D rain fields
        self%rain_ls (:,:,:) = self%del_rain_ls (:,:,:)
        self%rain_zm (:,:,:) = self%del_rain_zm (:,:,:)
        self%rain_hk (:,:,:) = self%del_rain_hk (:,:,:)
!... put total rain production into moistq
        self%moistq(:,:,:) = self%rain_hk(:,:,:) + self%rain_ls(:,:,:) + &
     &                       self%rain_zm(:,:,:)
        self%moistq(:,:,:) = -self%moistq(:,:,:)
!
      elseif (self%metdata_name_model(1:5) == 'GEOS5') then
!... GEOS5 3D rain fields
        self%rain_ls (:,:,:) = self%del_rain_ls (:,:,:)
        self%rain_zm (:,:,:) = self%del_rain_zm (:,:,:)
        self%rain_hk (:,:,:) = 0.0d0
        self%moistq (:,:,:) = self%del_moistq(:,:,:)
!
        ! For FastJX 6.5
        !---------------
        self%taucli (:,:,:) = self%del_taucli (:,:,:)
        self%tauclw (:,:,:) = self%del_tauclw (:,:,:)
        !---------------
      else
!
        self%moistq (:,:,:) = self%del_moistq(:,:,:)
!
      end if
!
      if ((self%metdata_name_org  (1:4) == 'NCAR') .and.  &
     &    (self%metdata_name_model(1:4) == 'CCM3')) then
!
        self%astrat (:,:,:) = self%del_astrat (:,:,:)
        self%clwc   (:,:,:) = self%del_clwc   (:,:,:)
        self%aconv  (:,:,:) = self%del_aconv  (:,:,:)
        self%relhum (:,:,:) = self%del_relhum (:,:,:)
        self%rhclear(:,:,:) = self%del_rhclear(:,:,:)
        self%cldmas (:,:,:) = self%del_cldmas (:,:,:)
        self%omega  (:,:,:) = self%del_omega  (:,:,:)
!
      end if
!
      self%max_cloud  (:,:,:) = self%del_max_cloud  (:,:,:)
      self%ran_cloud  (:,:,:) = self%del_ran_cloud  (:,:,:)
      self%humidity   (:,:,:) = self%del_humidity   (:,:,:)
      self%tau_cloud  (:,:,:) = self%del_tau_cloud  (:,:,:)
      self%kzz        (:,:,:) = self%del_kzz        (:,:,:)
!
      return
!
      end subroutine Save_Met2
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update_Met2
!
! !INTERFACE:
!
      subroutine Update_Met2 (self, new_met_rec)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: new_met_rec
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_MetFields), intent(inOut) :: self
!
! !DESCRIPTION:
! Updates the met2 data (everything except ps, u, v, & kel) each time step.
!
! !LOCAL VARIABLES:
      real*8  :: weight1, weight2
      real*8, save :: tot_met2_steps
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6, *) 'Update_Met2 called by ', procID
!
!     =====================
      if (self%do_timinterp_met) then
!     =====================
!
        if (new_met_rec) then
!
          tot_met2_steps = self%mdt / tdt
!
          weight1 = (tot_met2_steps - 0.5d0) / tot_met2_steps
          weight2 =                   0.5d0  / tot_met2_steps
!
        else
!
          weight1 = (tot_met2_steps - 0.5d0) / (tot_met2_steps + 0.5d0)
          weight2 =                   1.0d0  / (tot_met2_steps + 0.5d0)
!
        end if
!
!
        self%con_precip   (:,:) = weight1 *     self%con_precip   (:,:) +  &
     &                            weight2 * self%del_con_precip   (:,:)
        self%grnd_temp    (:,:) = weight1 *     self%grnd_temp    (:,:) +  &
     &                            weight2 * self%del_grnd_temp    (:,:)
        self%pbl          (:,:) = weight1 *     self%pbl          (:,:) +  &
     &                            weight2 * self%del_pbl          (:,:)
        self%radswg       (:,:) = weight1 *     self%radswg       (:,:) +  &
     &                            weight2 * self%del_radswg       (:,:)
        self%surf_air_temp(:,:) = weight1 *     self%surf_air_temp(:,:) +  &
     &                            weight2 * self%del_surf_air_temp(:,:)
        self%surf_rough   (:,:) = weight1 *     self%surf_rough   (:,:) +  &
     &                            weight2 * self%del_surf_rough   (:,:)
        self%tot_precip   (:,:) = weight1 *     self%tot_precip   (:,:) +  &
     &                            weight2 * self%del_tot_precip   (:,:)
        self%ustar        (:,:) = weight1 *     self%ustar        (:,:) +  &
     &                            weight2 * self%del_ustar        (:,:)
!
!
       if (sfalbedo_opt == 3) then
          self%saldif     (:,:) = weight1 *     self%saldif       (:,:) +  &
     &                            weight2 * self%del_saldif       (:,:)
          self%saldir     (:,:) = weight1 *     self%saldir       (:,:) +  &
     &                            weight2 * self%del_saldir       (:,:)
          self%sasdif     (:,:) = weight1 *     self%sasdif       (:,:) +  &
     &                            weight2 * self%del_sasdif       (:,:)
          self%sasdir     (:,:) = weight1 *     self%sasdir       (:,:) +  &
     &                            weight2 * self%del_sasdir       (:,:)
        end if
!
!
        if (uvalbedo_opt == 3) then
          self%surf_alb_uv(:,:) = weight1 *     self%surf_alb_uv  (:,:) +  &
     &                            weight2 * self%del_surf_alb_uv  (:,:)
        end if
!
!
        self%humidity   (:,:,:) = weight1 *     self%humidity   (:,:,:) +  &
     &                            weight2 * self%del_humidity   (:,:,:)
        self%kzz        (:,:,:) = weight1 *     self%kzz        (:,:,:) +  &
     &                            weight2 * self%del_kzz        (:,:,:)
        self%max_cloud  (:,:,:) = weight1 *     self%max_cloud  (:,:,:) +  &
     &                            weight2 * self%del_max_cloud  (:,:,:)
        self%ran_cloud  (:,:,:) = weight1 *     self%ran_cloud  (:,:,:) +  &
     &                            weight2 * self%del_ran_cloud  (:,:,:)
        self%tau_cloud  (:,:,:) = weight1 *     self%tau_cloud  (:,:,:) +  &
     &                            weight2 * self%del_tau_cloud  (:,:,:)
!
!
        if (convec_opt /= 0 .or. do_wetdep) then
!
          self%cmf(:,:,:) = weight1 * self%cmf    (:,:,:) +  &
     &                      weight2 * self%del_cmf(:,:,:)
!
        end if
!
        if (convec_opt /= 0 .or. do_drydep) then
!
          self%dtrn(:,:,:) = weight1 * self%dtrn    (:,:,:) +  &
     &                       weight2 * self%del_dtrn(:,:,:)
!
        end if
!
!
        if (convec_opt == 2) then
!
          self%ed(:,:,:) = weight1 * self%ed    (:,:,:) +  &
     &                     weight2 * self%del_ed(:,:,:)
          self%eu(:,:,:) = weight1 * self%eu    (:,:,:) +  &
     &                     weight2 * self%del_eu(:,:,:)
          self%md(:,:,:) = weight1 * self%md    (:,:,:) +  &
     &                     weight2 * self%del_md(:,:,:)
!
        end if
!
        if (self%metdata_name_model(1:5) == 'GEOS4' .or. &
     &      self%metdata_name_model(1:5) == 'GEOS5') then
           self%u10m(:,:) = weight1 * self%u10m    (:,:) +  &
     &                      weight2 * self%del_u10m(:,:)
           self%v10m(:,:) = weight1 * self%v10m    (:,:) +  &
     &                      weight2 * self%del_v10m(:,:)
           self%gwet(:,:) = weight1 * self%gwet    (:,:) +  &
     &                      weight2 * self%del_gwet(:,:)
           self%pardif(:,:) = weight1 * self%pardif    (:,:) +  &
     &                      weight2 * self%del_pardif(:,:)
           self%pardir(:,:) = weight1 * self%pardir    (:,:) +  &
     &                      weight2 * self%del_pardir(:,:)
        end if
!
        if (convec_opt == 3) then
!
          self%zmdu(:,:,:) = weight1 * self%zmdu    (:,:,:) +  &
     &                       weight2 * self%del_zmdu(:,:,:)
          self%zmeu(:,:,:) = weight1 * self%zmeu    (:,:,:) +  &
     &                       weight2 * self%del_zmeu(:,:,:)
          self%zmed(:,:,:) = weight1 * self%zmed    (:,:,:) +  &
     &                       weight2 * self%del_zmed(:,:,:)
          self%zmmd(:,:,:) = weight1 * self%zmmd    (:,:,:) +  &
     &                       weight2 * self%del_zmmd(:,:,:)
          self%zmmu(:,:,:) = weight1 * self%zmmu    (:,:,:) +  &
     &                       weight2 * self%del_zmmu(:,:,:)
          self%hkdu(:,:,:) = weight1 * self%hkdu    (:,:,:) +  &
     &                       weight2 * self%del_hkdu(:,:,:)
          self%hkeu(:,:,:) = weight1 * self%hkeu    (:,:,:) +  &
     &                       weight2 * self%del_hkeu(:,:,:)
          self%hkmu(:,:,:) = weight1 * self%hkmu    (:,:,:) +  &
     &                       weight2 * self%del_hkmu(:,:,:)
!
        end if
!
!
!
        if (((self%metdata_name_org  (1:4) == 'NCAR' ) .and.  &
     &       (self%metdata_name_model(1:5) == 'MATCH')) .or.  &
     &      (self%metdata_name_org(1:4) == 'GISS')) then
!
          self%rain  (:,:,:) = weight1 *     self%rain  (:,:,:) +  &
     &                         weight2 * self%del_rain  (:,:,:)
!
        elseif ((self%metdata_name_org(1:4) /= 'GMAO') .and.  &
     &      (self%metdata_name_model(1:5) == 'GEOS4')) then
!
          self%rain_ls  (:,:,:) = weight1 *  self%rain_ls  (:,:,:) +  &
     &                            weight2 * self%del_rain_ls  (:,:,:)
          self%rain_zm  (:,:,:) = weight1 *  self%rain_zm  (:,:,:) +  &
     &                            weight2 * self%del_rain_zm  (:,:,:)
          self%rain_hk  (:,:,:) = weight1 *  self%rain_hk  (:,:,:) +  &
     &                            weight2 * self%del_rain_hk  (:,:,:)
!.... this is not correct for geos5
          self%moistq = self%rain_ls + self%rain_zm + self%rain_hk
          self%moistq = -self%moistq
        elseif ((self%metdata_name_org(1:4) /= 'GMAO') .and.  &
     &      (self%metdata_name_model(1:5) == 'GEOS5')) then
!
          self%rain_ls  (:,:,:) = weight1 *  self%rain_ls  (:,:,:) +  &
     &                            weight2 * self%del_rain_ls  (:,:,:)
          self%rain_zm  (:,:,:) = weight1 *  self%rain_zm  (:,:,:) +  &
     &                            weight2 * self%del_rain_zm  (:,:,:)
!
          self%moistq  (:,:,:) = weight1 *     self%moistq  (:,:,:) +  &
     &                           weight2 * self%del_moistq  (:,:,:)
!
        ! For fastJX 6.5
        !---------------
          self%taucli  (:,:,:) = weight1 *  self%taucli  (:,:,:) +  &
     &                            weight2 * self%del_taucli  (:,:,:)
          self%tauclw  (:,:,:) = weight1 *  self%tauclw  (:,:,:) +  &
     &                            weight2 * self%del_tauclw  (:,:,:)
        !---------------
        else
!
          self%moistq(:,:,:) = weight1 *     self%moistq(:,:,:) +  &
     &                         weight2 * self%del_moistq(:,:,:)
!
        end if
!
!
        if ((self%metdata_name_org  (1:4) == 'NCAR') .and.  &
     &      (self%metdata_name_model(1:4) == 'CCM3')) then
!
          self%aconv  (:,:,:) = weight1 *     self%aconv  (:,:,:) +  &
     &                          weight2 * self%del_aconv  (:,:,:)
          self%astrat (:,:,:) = weight1 *     self%astrat (:,:,:) +  &
     &                          weight2 * self%del_astrat (:,:,:)
          self%clwc   (:,:,:) = weight1 *     self%clwc   (:,:,:) +  &
     &                          weight2 * self%del_clwc   (:,:,:)
          self%relhum (:,:,:) = weight1 *     self%relhum (:,:,:) +  &
     &                          weight2 * self%del_relhum (:,:,:)
          self%rhclear(:,:,:) = weight1 *     self%rhclear(:,:,:) +  &
     &                          weight2 * self%del_rhclear(:,:,:)
!
          self%cldmas (:,:,:) = weight1 *     self%cldmas (:,:,:) +  &
     &                          weight2 * self%del_cldmas (:,:,:)
          self%md     (:,:,:) = weight1 *     self%md     (:,:,:) +  &
     &                          weight2 * self%del_md     (:,:,:)
          self%omega  (:,:,:) = weight1 *     self%omega  (:,:,:) +  &
     &                          weight2 * self%del_omega  (:,:,:)
!
        end if
!
        tot_met2_steps = tot_met2_steps - 1.0d0
!
!     ======
      end if
!     ======
!
      return
!
      end subroutine Update_Met2
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Control_Met2_File_Input
!
! !INTERFACE:
!
      subroutine Control_Met2_File_Input (self, gmiDomain, new_met_rec, out_of_data)
!
      implicit none
!
! !INPUT PARAMETERS:
      type (t_gmiDomain), intent(in) :: gmiDomain
      logical :: new_met_rec   ! new metFields record?
!
! !OUTPUT PARAMETERS:
      logical :: out_of_data   ! out of data?
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_MetFields), intent(inOut) :: self
!
! !DESCRIPTION:
!   This routine controls the met2 (everything except ps, u, v, & kel) file
!   input.
!
! !LOCAL VARIABLES:
    real*8 , allocatable :: mcor(:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
      out_of_data = .false.
!
      if (self%m2rnum_in > self%m2num_recs) then
!
        out_of_data = .false.
!
!       ======================
        call Get_Next_Met_File  &
!       ======================
     &    (self%do_cycle_met, out_of_data, self%met2_infile_num,  &
     &     self%met_num_infiles, self%m2num_recs, self%m2rnum_in, self%ncid_met2,  &
     &     self%met_infile_names, pr_diag, procID)
!
!       =======================
        if (out_of_data) return  ! Early return.
!       =======================
!
      end if
!
      if (new_met_rec) then
!
        self%tmet2_1 = self%tmet2_2
!
!       ==============
        call Save_Met2  (self)
!       ==============
!
        if (self%metdata_name_model(1:5) == 'GEOS5'  &
     &       .or. self%metdata_name_org(1:2) == 'EC') then
           allocate(mcor(i1:i2,ju1:j2))
           call Get_mcor(gmiDomain, mcor)
!
           call calcTotalMass (self%pt, self%ai, self%bi, self%pctm1, mcor, &
     &          self%mass, pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, &
     &          ihi, julo, jhi)
        end if
!
!       ==============
        call Read_Met2  &
!       ==============
     &    (mcor, self%mass, &
     &     self%metdata_name_org, self%metdata_name_model, do_wetdep, &
     &     self%ncid_met2, self%m2rnum_in, chem_opt, convec_opt, emiss_in_opt, &
     &     self%gwet_opt, sfalbedo_opt, uvalbedo_opt, self%lwi_flags, &
     &     self%del_con_precip, self%del_grnd_temp, self%del_pbl, &
     &     self%del_radswg, self%del_saldif, self%del_saldir, self%del_sasdif,  &
     &     self%del_sasdir, self%del_surf_air_temp, self%del_surf_alb_uv,  &
     &     self%del_surf_rough, self%del_tot_precip, self%del_ustar, &
     &     self%del_md, self%del_cmf, self%del_dtrn, self%del_ed, self%del_eu, &
     &     self%del_humidity, self%del_kzz, self%del_max_cloud, &
     &     self%del_moistq, self%del_rain, self%del_ran_cloud,  &
     &     self%del_tau_cloud, self%del_aconv, self%del_astrat, self%del_clwc, &
     &     self%del_relhum, self%del_rhclear, self%del_cldmas, self%del_omega,  &
     &     self%del_rain_zm, self%del_rain_hk, self%del_rain_ls,  &
     &     self%del_taucli, self%del_tauclw, &
     &     self%del_zmdu, self%del_zmeu, self%del_zmed, self%del_zmmd, &
     &     self%del_zmmu, self%del_hkdu, self%del_hkeu, self%del_hkmu, &
     &     self%del_u10m, self%del_v10m, self%del_gwet, self%del_pardif, &
     &     self%del_pardir, pr_diag, procID, i1, i2, ju1, j2, k1, k2, i1_gl, &
     &     ju1_gl, j2_gl)
!
        self%m2rnum_in = self%m2rnum_in + 1
!
        self%tmet2_2   = Nint (self%tmet2_1 + self%mdt)
!
      end if
!
      return
!
      end subroutine Control_Met2_File_Input
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finishMet
!
! !INTERFACE:
!
      subroutine finishMet (self, new_met_rec, t_cloud_ice)
!
! !USES:
      use GmiUtilsMetFields_mod, only : Calc_Humidity, Convert_Pbl
      implicit none
!
! !INPUT PARAMETERS:
!      logical, intent(in) :: pr_diag
!      integer, intent(in) :: procID
!      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      logical, intent(in) :: new_met_rec
      real*8 , intent(in) :: t_cloud_ice
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_metFields), intent(inOut) :: self
!
! !DESCRIPTION:
! Does some final calculations/conversions of humidity and pbl, if necessary.
!
! !LOCAL VARIABLES:
      logical, save :: first_pbl = .true.
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'finishMet called by ', procID
!
     if (new_met_rec) then
!
        if ((self%metdata_name_org  (1:4) == 'NCAR') .and.  &
     &      (self%metdata_name_model(1:4) == 'CCM3')) then
!
!         ==================
          call Calc_Humidity  &
!         ==================
     &      (t_cloud_ice, self%kel, self%press3e, self%relhum, self%humidity, &
     &       pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
!         ==================
          call Calc_Humidity  &
!         ==================
     &      (t_cloud_ice, self%del_kel, self%press3e, self%del_relhum, &
     &       self%del_humidity, &
     &       pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
        end if
!
        if (self%metdata_name_org(1:4) /= 'GISS' .and.  &
     &      self%metdata_name_org(1:2) /= 'EC' .and.  &
     &      (self%metdata_name_org  (1:4) /= 'GMAO' .and.  &
     &       (self%metdata_name_model(1:5) /= 'GEOS4' .or. &
     &        self%metdata_name_model(1:5) == 'GEOS5'))) then
!
          if (first_pbl) then
!
            first_pbl = .false.
!
!           ================
            call Convert_Pbl  &
!           ================
     &        (self%psx, self%surf_air_temp, self%humidity, self%pbl, &
     &         pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
          end if
!
!         ================
          call Convert_Pbl  &
!         ================
     &      (self%del_psx, self%del_surf_air_temp, self%del_humidity, &
     &       self%del_pbl, &
     &       pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
        end if
!
      end if
!
      return
!
      end subroutine finishMet
!EOC
!------------------------------------------------------------------------------
! Set/Get routines
!------------------------------------------------------------------------------
      subroutine Get_met_infile_names (self, met_infile_names)
        implicit none
        character (len=MAX_LENGTH_FILE_NAME), intent(out) :: met_infile_names(:)
        type(t_metFields), intent(in) :: self
        met_infile_names(1:self%met_num_infiles) = &
     &               self%met_infile_names(1:self%met_num_infiles)
        return
      end subroutine Get_met_infile_names
!------------------------------------------------------------------------------
      subroutine Get_metdata_name (self, metdata_name)
        implicit none
        character(len=*), intent(out) :: metdata_name
        type(t_metFields), intent(in) :: self
        metdata_name = self%metdata_name
        return
      end subroutine Get_metdata_name
!------------------------------------------------------------------------------
      subroutine Get_metdata_name_org (self, metdata_name_org)
        implicit none
        character(len=*), intent(out) :: metdata_name_org
        type(t_metFields), intent(in) :: self
        metdata_name_org = self%metdata_name_org
        return
      end subroutine Get_metdata_name_org
!------------------------------------------------------------------------------
      subroutine Get_metdata_name_model (self, metdata_name_model)
        implicit none
        character(len=*), intent(out) :: metdata_name_model
        type(t_metFields), intent(in) :: self
        metdata_name_model = self%metdata_name_model
        return
      end subroutine Get_metdata_name_model
!------------------------------------------------------------------------------
      subroutine Get_metdata_name_dims (self, metdata_name_dims)
        implicit none
        character(len=*), intent(out) :: metdata_name_dims
        type(t_metFields), intent(in) :: self
        metdata_name_dims = self%metdata_name_dims
        return
      end subroutine Get_metdata_name_dims
!------------------------------------------------------------------------------
      subroutine Get_met_grid_type (self, met_grid_type)
        implicit none
        character(len=*), intent(out) :: met_grid_type
        type(t_metFields), intent(in) :: self
        met_grid_type = self%met_grid_type
        return
      end subroutine Get_met_grid_type
!------------------------------------------------------------------------------
      subroutine Get_do_wind_pole (self, do_wind_pole)
        implicit none
        logical, intent(out) :: do_wind_pole
        type(t_metFields), intent(in) :: self
        do_wind_pole = self%do_wind_pole
        return
      end subroutine Get_do_wind_pole
!------------------------------------------------------------------------------
      subroutine Get_do_timinterp_met (self, do_timinterp_met)
        implicit none
        logical, intent(out) :: do_timinterp_met
        type(t_metFields), intent(in) :: self
        do_timinterp_met = self%do_timinterp_met
        return
      end subroutine Get_do_timinterp_met
!------------------------------------------------------------------------------
      subroutine Get_do_timinterp_winds (self, do_timinterp_winds)
        implicit none
        logical, intent(out) :: do_timinterp_winds
        type(t_metFields), intent(in) :: self
        do_timinterp_winds = self%do_timinterp_winds
        return
      end subroutine Get_do_timinterp_winds
!------------------------------------------------------------------------------
      subroutine Get_do_cycle_met (self, do_cycle_met)
        implicit none
        logical, intent(out) :: do_cycle_met
        type(t_metFields), intent(in) :: self
        do_cycle_met = self%do_cycle_met
        return
      end subroutine Get_do_cycle_met
!------------------------------------------------------------------------------
      subroutine Get_gwet_opt (self, gwet_opt)
        implicit none
        integer, intent(out) :: gwet_opt
        type(t_metFields), intent(in) :: self
        gwet_opt = self%gwet_opt
        return
      end subroutine Get_gwet_opt
!------------------------------------------------------------------------------
      subroutine Set_m2rnum_in (self, m2rnum_in)
        implicit none
        integer, intent(in) :: m2rnum_in
        type(t_metFields), intent(inOut) :: self
        self%m2rnum_in = m2rnum_in
        return
      end subroutine Set_m2rnum_in
!------------------------------------------------------------------------------
      subroutine Get_m2rnum_in (self, m2rnum_in)
        implicit none
        integer, intent(out) :: m2rnum_in
        type(t_metFields), intent(in) :: self
        m2rnum_in = self%m2rnum_in
        return
      end subroutine Get_m2rnum_in
!------------------------------------------------------------------------------
      subroutine Set_met2_infile_num (self, met2_infile_num)
        implicit none
        integer, intent(in) :: met2_infile_num
        type(t_metFields), intent(inOut) :: self
        self%met2_infile_num = met2_infile_num
        return
      end subroutine Set_met2_infile_num
!------------------------------------------------------------------------------
      subroutine Get_met2_infile_num (self, met2_infile_num)
        implicit none
        integer, intent(out) :: met2_infile_num
        type(t_metFields), intent(in) :: self
        met2_infile_num = self%met2_infile_num
        return
      end subroutine Get_met2_infile_num
!------------------------------------------------------------------------------
      subroutine Set_ncid_met2 (self, ncid_met2)
        implicit none
        integer, intent(in) :: ncid_met2
        type(t_metFields), intent(inOut) :: self
        self%ncid_met2 = ncid_met2
        return
      end subroutine Set_ncid_met2
!------------------------------------------------------------------------------
      subroutine Get_ncid_met2 (self, ncid_met2)
        implicit none
        integer, intent(out) :: ncid_met2
        type(t_metFields), intent(in) :: self
        ncid_met2 = self%ncid_met2
        return
      end subroutine Get_ncid_met2
!------------------------------------------------------------------------------
      subroutine Set_tmet2_1 (self, tmet2_1)
        implicit none
        real*8, intent(in) :: tmet2_1
        type(t_metFields), intent(inOut) :: self
        self%tmet2_1 = tmet2_1
        return
      end subroutine Set_tmet2_1
!------------------------------------------------------------------------------
      subroutine Get_tmet2_1 (self, tmet2_1)
        implicit none
        real*8, intent(out) :: tmet2_1
        type(t_metFields), intent(in) :: self
        tmet2_1 = self%tmet2_1
        return
      end subroutine Get_tmet2_1
!------------------------------------------------------------------------------
      subroutine Set_tmet2_2 (self, tmet2_2)
        implicit none
        real*8, intent(in) :: tmet2_2
        type(t_metFields), intent(inOut) :: self
        self%tmet2_2 = tmet2_2
        return
      end subroutine Set_tmet2_2
!------------------------------------------------------------------------------
      subroutine Get_tmet2_2 (self, tmet2_2)
        implicit none
        real*8, intent(out) :: tmet2_2
        type(t_metFields), intent(in) :: self
        tmet2_2 = self%tmet2_2
        return
      end subroutine Get_tmet2_2
!------------------------------------------------------------------------------
      subroutine Set_m2num_recs (self, m2num_recs)
        implicit none
        integer, intent(in) :: m2num_recs
        type(t_metFields), intent(inOut) :: self
        self%m2num_recs = m2num_recs
        return
      end subroutine Set_m2num_recs
!------------------------------------------------------------------------------
      subroutine Get_m2num_recs (self, m2num_recs)
        implicit none
        integer, intent(out) :: m2num_recs
        type(t_metFields), intent(in) :: self
        m2num_recs = self%m2num_recs
        return
      end subroutine Get_m2num_recs
!------------------------------------------------------------------------------
      subroutine Set_m1rnum_in (self, m1rnum_in)
        implicit none
        integer, intent(in) :: m1rnum_in
        type(t_metFields), intent(inOut) :: self
        self%m1rnum_in = m1rnum_in
        return
      end subroutine Set_m1rnum_in
!------------------------------------------------------------------------------
      subroutine Get_m1rnum_in (self, m1rnum_in)
        implicit none
        integer, intent(out) :: m1rnum_in
        type(t_metFields), intent(in) :: self
        m1rnum_in = self%m1rnum_in
        return
      end subroutine Get_m1rnum_in
!------------------------------------------------------------------------------
      subroutine Set_met1_infile_num (self, met1_infile_num)
        implicit none
        integer, intent(in) :: met1_infile_num
        type(t_metFields), intent(inOut) :: self
        self%met1_infile_num = met1_infile_num
        return
      end subroutine Set_met1_infile_num
!------------------------------------------------------------------------------
      subroutine Get_met1_infile_num (self, met1_infile_num)
        implicit none
        integer, intent(out) :: met1_infile_num
        type(t_metFields), intent(in) :: self
        met1_infile_num = self%met1_infile_num
        return
      end subroutine Get_met1_infile_num
!------------------------------------------------------------------------------
      subroutine Set_ncid_met1 (self, ncid_met1)
        implicit none
        integer, intent(in) :: ncid_met1
        type(t_metFields), intent(inOut) :: self
        self%ncid_met1 = ncid_met1
        return
      end subroutine Set_ncid_met1
!------------------------------------------------------------------------------
      subroutine Get_ncid_met1 (self, ncid_met1)
        implicit none
        integer, intent(out) :: ncid_met1
        type(t_metFields), intent(in) :: self
        ncid_met1 = self%ncid_met1
        return
      end subroutine Get_ncid_met1
!------------------------------------------------------------------------------
      subroutine Set_tmet1_1 (self, tmet1_1)
        implicit none
        real*8, intent(in) :: tmet1_1
        type(t_metFields), intent(inOut) :: self
        self%tmet1_1 = tmet1_1
        return
      end subroutine Set_tmet1_1
!------------------------------------------------------------------------------
      subroutine Get_tmet1_1 (self, tmet1_1)
        implicit none
        real*8, intent(out) :: tmet1_1
        type(t_metFields), intent(in) :: self
        tmet1_1 = self%tmet1_1
        return
      end subroutine Get_tmet1_1
!------------------------------------------------------------------------------
      subroutine Set_tmet1_2 (self, tmet1_2)
        implicit none
        real*8, intent(in) :: tmet1_2
        type(t_metFields), intent(inOut) :: self
        self%tmet1_2 = tmet1_2
        return
      end subroutine Set_tmet1_2
!------------------------------------------------------------------------------
      subroutine Get_tmet1_2 (self, tmet1_2)
        implicit none
        real*8, intent(out) :: tmet1_2
        type(t_metFields), intent(in) :: self
        tmet1_2 = self%tmet1_2
        return
      end subroutine Get_tmet1_2
!------------------------------------------------------------------------------
      subroutine Set_m1num_recs (self, m1num_recs)
        implicit none
        integer, intent(in) :: m1num_recs
        type(t_metFields), intent(inOut) :: self
        self%m1num_recs = m1num_recs
        return
      end subroutine Set_m1num_recs
!------------------------------------------------------------------------------
      subroutine Get_m1num_recs (self, m1num_recs)
        implicit none
        integer, intent(out) :: m1num_recs
        type(t_metFields), intent(in) :: self
        m1num_recs = self%m1num_recs
        return
      end subroutine Get_m1num_recs
!------------------------------------------------------------------------------
      subroutine Get_mrnum_in (self, mrnum_in)
        implicit none
        integer, intent(out) :: mrnum_in
        type(t_metFields), intent(in) :: self
        mrnum_in = self%mrnum_in
        return
      end subroutine Get_mrnum_in
!------------------------------------------------------------------------------
      subroutine Get_met_infile_num (self, met_infile_num)
        implicit none
        integer, intent(out) :: met_infile_num
        type(t_metFields), intent(in) :: self
        met_infile_num = self%met_infile_num
        return
      end subroutine Get_met_infile_num
!------------------------------------------------------------------------------
      subroutine Get_met_num_infiles (self, met_num_infiles)
        implicit none
        integer, intent(out) :: met_num_infiles
        type(t_metFields), intent(in) :: self
        met_num_infiles = self%met_num_infiles
        return
      end subroutine Get_met_num_infiles
!------------------------------------------------------------------------------
      subroutine Get_metNumMEGAN (self, metNumMEGAN)
        implicit none
        integer, intent(out) :: metNumMEGAN
        type(t_metFields), intent(in) :: self
        metNumMEGAN = self%metNumMEGAN
        return
      end subroutine Get_metNumMEGAN
!------------------------------------------------------------------------------
      subroutine Get_met_opt (self, met_opt)
        implicit none
        integer, intent(out) :: met_opt
        type(t_metFields), intent(in) :: self
        met_opt = self%met_opt
        return
      end subroutine Get_met_opt
!------------------------------------------------------------------------------
      subroutine Get_mdt (self, mdt)
        implicit none
        real*8, intent(out) :: mdt
        type(t_metFields), intent(in) :: self
        mdt = self%mdt
        return
      end subroutine Get_mdt
!------------------------------------------------------------------------------
      subroutine Get_pt (self, pt)
        implicit none
        real*8, intent(out) :: pt
        type(t_metFields), intent(in) :: self
        pt = self%pt
        return
      end subroutine Get_pt
!------------------------------------------------------------------------------
      subroutine Get_ptop (self, ptop)
        implicit none
        real*8, intent(out) :: ptop
        type(t_metFields), intent(in) :: self
        ptop = self%ptop
        return
      end subroutine Get_ptop
!------------------------------------------------------------------------------
      subroutine Get_dap (self, dap)
        implicit none
        real*8, intent(out) :: dap(:)
        type(t_metFields), intent(in) :: self
        dap(:) = self%dap(:)
        return
      end subroutine Get_dap
!------------------------------------------------------------------------------
      subroutine Get_dbk (self, dbk)
        implicit none
        real*8, intent(out) :: dbk(:)
        type(t_metFields), intent(in) :: self
        dbk(:) = self%dbk(:)
        return
      end subroutine Get_dbk
!------------------------------------------------------------------------------
      subroutine Get_ai (self, ai)
        implicit none
        real*8, intent(out) :: ai(:)
        type(t_metFields), intent(in) :: self
        ai(:) = self%ai(:)
        return
      end subroutine Get_ai
!------------------------------------------------------------------------------
      subroutine Get_bi (self, bi)
        implicit none
        real*8, intent(out) :: bi(:)
        type(t_metFields), intent(in) :: self
        bi(:) = self%bi(:)
        return
      end subroutine Get_bi
!------------------------------------------------------------------------------
      subroutine Get_am (self, am)
        implicit none
        real*8, intent(out) :: am(:)
        type(t_metFields), intent(in) :: self
        am(:) = self%am(:)
        return
      end subroutine Get_am
!------------------------------------------------------------------------------
      subroutine Get_bm (self, bm)
        implicit none
        real*8, intent(out) :: bm(:)
        type(t_metFields), intent(in) :: self
        bm(:) = self%bm(:)
        return
      end subroutine Get_bm
!------------------------------------------------------------------------------
      subroutine Get_lwi_flags (self, lwi_flags)
        implicit none
        integer, intent(out) :: lwi_flags(:,:)
        type(t_metFields), intent(in) :: self
        lwi_flags(:,:) = self%lwi_flags(:,:)
        return
      end subroutine Get_lwi_flags
!------------------------------------------------------------------------------
      subroutine Set_cmi_flags (self, cmi_flags)
        implicit none
        integer, intent(in) :: cmi_flags(:,:)
        type(t_metFields), intent(inOut) :: self
        self%cmi_flags(:,:) = cmi_flags(:,:)
        return
      end subroutine Set_cmi_flags
!------------------------------------------------------------------------------
      subroutine Get_cmi_flags (self, cmi_flags)
        implicit none
        integer, intent(out) :: cmi_flags(:,:)
        type(t_metFields), intent(in) :: self
        cmi_flags(:,:) = self%cmi_flags(:,:)
        return
      end subroutine Get_cmi_flags
!------------------------------------------------------------------------------
      subroutine Set_mass (self, mass)
        implicit none
        real*8 , intent(in) :: mass(:,:,:)
        type(t_metFields), intent(inOut) :: self
        self%mass(:,:,:) = mass(:,:,:)
        return
      end subroutine Set_mass
!------------------------------------------------------------------------------
      subroutine Get_mass (self, mass)
        implicit none
        real*8 , intent(out) :: mass(:,:,:)
        type(t_metFields), intent(in) :: self
        mass(:,:,:) = self%mass(:,:,:)
        return
      end subroutine Get_mass
!------------------------------------------------------------------------------
      subroutine Set_tropopausePress (self, tropopausePress)
        implicit none
        real*8 , intent(in) :: tropopausePress(:,:)
        type(t_metFields), intent(inOut) :: self
        self%tropopausePress(:,:) = tropopausePress(:,:)
        return
      end subroutine Set_tropopausePress
!------------------------------------------------------------------------------
      subroutine Get_pctm2Glob (self, pctm2Glob)
        implicit none
        real*8 , intent(out) :: pctm2Glob(:,:)
        type(t_metFields), intent(in) :: self
        pctm2Glob(:,:) = self%pctm2Glob(:,:)
        return
      end subroutine Get_pctm2Glob
!------------------------------------------------------------------------------
      subroutine Set_pctm2Glob (self, pctm2Glob)
        implicit none
        real*8 , intent(in) :: pctm2Glob(:,:)
        type(t_metFields), intent(inOut) :: self
        self%pctm2Glob(:,:) = pctm2Glob(:,:)
        return
      end subroutine Set_pctm2Glob
!------------------------------------------------------------------------------
      subroutine Get_pctm1Glob (self, pctm1Glob)
        implicit none
        real*8 , intent(out) :: pctm1Glob(:,:)
        type(t_metFields), intent(in) :: self
        pctm1Glob(:,:) = self%pctm1Glob(:,:)
        return
      end subroutine Get_pctm1Glob
!------------------------------------------------------------------------------
      subroutine Get_kelGlob (self, kelGlob)
        implicit none
        real*8 , intent(out) :: kelGlob(:,:,:)
        type(t_metFields), intent(in) :: self
        kelGlob(:,:,:) = self%kelGlob(:,:,:)
        return
      end subroutine Get_kelGlob
!------------------------------------------------------------------------------
      subroutine Set_pctm1Glob (self, pctm1Glob)
        implicit none
        real*8 , intent(in) :: pctm1Glob(:,:)
        type(t_metFields), intent(inOut) :: self
        self%pctm1Glob(:,:) = pctm1Glob(:,:)
        return
      end subroutine Set_pctm1Glob
!------------------------------------------------------------------------------
      subroutine Get_tropopausePress (self, tropopausePress)
        implicit none
        real*8 , intent(out) :: tropopausePress(:,:)
        type(t_metFields), intent(in) :: self
        tropopausePress(:,:) = self%tropopausePress(:,:)
        return
      end subroutine Get_tropopausePress
!------------------------------------------------------------------------------
      subroutine Set_potentialVorticity (self, potentialVorticity)
        implicit none
        real*8 , intent(in) :: potentialVorticity(:,:,:)
        type(t_metFields), intent(inOut) :: self
        self%potentialVorticity(:,:,:) = potentialVorticity(:,:,:)
        return
      end subroutine Set_potentialVorticity
!------------------------------------------------------------------------------
      subroutine Get_potentialVorticity (self, potentialVorticity)
        implicit none
        real*8 , intent(out) :: potentialVorticity(:,:,:)
        type(t_metFields), intent(in) :: self
        potentialVorticity(:,:,:) = self%potentialVorticity(:,:,:)
        return
      end subroutine Get_potentialVorticity
!------------------------------------------------------------------------------
      subroutine Set_potentialTemp (self, potentialTemp)
        implicit none
        real*8 , intent(in) :: potentialTemp(:,:,:)
        type(t_metFields), intent(inOut) :: self
        self%potentialTemp(:,:,:) = potentialTemp(:,:,:)
        return
      end subroutine Set_potentialTemp
!------------------------------------------------------------------------------
      subroutine Get_potentialTemp (self, potentialTemp)
        implicit none
        real*8 , intent(out) :: potentialTemp(:,:,:)
        type(t_metFields), intent(in) :: self
        potentialTemp(:,:,:) = self%potentialTemp(:,:,:)
        return
      end subroutine Get_potentialTemp
!------------------------------------------------------------------------------
      subroutine Set_relativeHumidity (self, relativeHumidity)
        implicit none
        real*8 , intent(in) :: relativeHumidity(:,:,:)
        type(t_metFields), intent(inOut) :: self
        self%relativeHumidity(:,:,:) = relativeHumidity(:,:,:)
        return
      end subroutine Set_relativeHumidity
!------------------------------------------------------------------------------
      subroutine Get_relativeHumidity (self, relativeHumidity)
        implicit none
        real*8 , intent(out) :: relativeHumidity(:,:,:)
        type(t_metFields), intent(in) :: self
        relativeHumidity(:,:,:) = self%relativeHumidity(:,:,:)
        return
      end subroutine Get_relativeHumidity
!------------------------------------------------------------------------------
      subroutine Set_gridBoxHeight (self, gridBoxHeight)
        implicit none
        real*8 , intent(in) :: gridBoxHeight(:,:,:)
        type(t_metFields), intent(inOut) :: self
        self%gridBoxHeight(:,:,:) = gridBoxHeight(:,:,:)
        return
      end subroutine Set_gridBoxHeight
!------------------------------------------------------------------------------
      subroutine Get_gridBoxHeight (self, gridBoxHeight)
        implicit none
        real*8 , intent(out) :: gridBoxHeight(:,:,:)
        type(t_metFields), intent(in) :: self
        gridBoxHeight(:,:,:) = self%gridBoxHeight(:,:,:)
        return
      end subroutine Get_gridBoxHeight
!------------------------------------------------------------------------------
      subroutine Set_press3e (self, press3e)
        implicit none
        real*8 , intent(in) :: press3e(:,:,:)
        type(t_metFields), intent(inOut) :: self
        self%press3e(:,:,:) = press3e(:,:,:)
        return
      end subroutine Set_press3e
!------------------------------------------------------------------------------
      subroutine Get_press3e (self, press3e)
        implicit none
        real*8 , intent(out) :: press3e(:,:,:)
        type(t_metFields), intent(in) :: self
        press3e(:,:,:) = self%press3e(:,:,:)
        return
      end subroutine Get_press3e
!------------------------------------------------------------------------------
      subroutine Set_press3c (self, press3c)
        implicit none
        real*8 , intent(in) :: press3c(:,:,:)
        type(t_metFields), intent(inOut) :: self
        self%press3c(:,:,:) = press3c(:,:,:)
        return
      end subroutine Set_press3c
!------------------------------------------------------------------------------
      subroutine Get_press3c (self, press3c)
        implicit none
        real*8 , intent(out) :: press3c(:,:,:)
        type(t_metFields), intent(in) :: self
        press3c(:,:,:) = self%press3c(:,:,:)
        return
      end subroutine Get_press3c
!------------------------------------------------------------------------------
      subroutine Set_surfTemp15DayAvg (self, surfTemp15DayAvg)
        implicit none
        real*8 , intent(in) :: surfTemp15DayAvg(:,:)
        type(t_metFields), intent(inOut) :: self
        self%surfTemp15DayAvg(:,:) = surfTemp15DayAvg(:,:)
        return
      end subroutine Set_surfTemp15DayAvg
!------------------------------------------------------------------------------
      subroutine Get_surfTemp15DayAvg (self, surfTemp15DayAvg)
        implicit none
        real*8 , intent(out) :: surfTemp15DayAvg(:,:)
        type(t_metFields), intent(in) :: self
        surfTemp15DayAvg(:,:) = self%surfTemp15DayAvg(:,:)
        return
      end subroutine Get_surfTemp15DayAvg
!------------------------------------------------------------------------------
      subroutine Get_pctm2 (self, pctm2)
        implicit none
        real*8 , intent(out) :: pctm2(:,:)
        type(t_metFields), intent(in) :: self
        pctm2(:,:) = self%pctm2(:,:)
        return
      end subroutine Get_pctm2
!------------------------------------------------------------------------------
      subroutine Get_pctm1 (self, pctm1)
        implicit none
        real*8 , intent(out) :: pctm1(:,:)
        type(t_metFields), intent(in) :: self
        pctm1(:,:) = self%pctm1(:,:)
        return
      end subroutine Get_pctm1
!------------------------------------------------------------------------------
      subroutine Get_kel (self, kel)
        implicit none
        real*8 , intent(out) :: kel(:,:,:)
        type(t_metFields), intent(in) :: self
        kel(:,:,:) = self%kel(:,:,:)
        return
      end subroutine Get_kel
!------------------------------------------------------------------------------
      subroutine Get_uux (self, uux)
        implicit none
        real*8 , intent(out) :: uux(:,:,:)
        type(t_metFields), intent(in) :: self
        uux(:,:,:) = self%uux(:,:,:)
        return
      end subroutine Get_uux
!------------------------------------------------------------------------------
      subroutine Get_vvx (self, vvx)
        implicit none
        real*8 , intent(out) :: vvx(:,:,:)
        type(t_metFields), intent(in) :: self
        vvx(:,:,:) = self%vvx(:,:,:)
        return
      end subroutine Get_vvx
!------------------------------------------------------------------------------
      subroutine Get_xmass (self, xmass)
        implicit none
        real*8 , intent(out) :: xmass(:,:,:)
        type(t_metFields), intent(in) :: self
        xmass(:,:,:) = self%xmass(:,:,:)
        return
      end subroutine Get_xmass
!------------------------------------------------------------------------------
      subroutine Get_ymass (self, ymass)
        implicit none
        real*8 , intent(out) :: ymass(:,:,:)
        type(t_metFields), intent(in) :: self
        ymass(:,:,:) = self%ymass(:,:,:)
        return
      end subroutine Get_ymass
!------------------------------------------------------------------------------
      subroutine Get_zmass (self, zmass)
        implicit none
        real*8 , intent(out) :: zmass(:,:,:)
        type(t_metFields), intent(in) :: self
        zmass(:,:,:) = self%zmass(:,:,:)
        return
      end subroutine Get_zmass
!------------------------------------------------------------------------------
      subroutine Get_con_precip (self, con_precip)
        implicit none
        real*8 , intent(out) :: con_precip(:,:)
        type(t_metFields), intent(in) :: self
        con_precip(:,:) = self%con_precip(:,:)
        return
      end subroutine Get_con_precip
!------------------------------------------------------------------------------
      subroutine Get_grnd_temp (self, grnd_temp)
        implicit none
        real*8 , intent(out) :: grnd_temp(:,:)
        type(t_metFields), intent(in) :: self
        grnd_temp(:,:) = self%grnd_temp(:,:)
        return
      end subroutine Get_grnd_temp
!------------------------------------------------------------------------------
      subroutine Get_pbl (self, pbl)
        implicit none
        real*8 , intent(out) :: pbl(:,:)
        type(t_metFields), intent(in) :: self
        pbl(:,:) = self%pbl(:,:)
        return
      end subroutine Get_pbl
!------------------------------------------------------------------------------
      subroutine Get_radswg (self, radswg)
        implicit none
        real*8 , intent(out) :: radswg(:,:)
        type(t_metFields), intent(in) :: self
        radswg(:,:) = self%radswg(:,:)
        return
      end subroutine Get_radswg
!------------------------------------------------------------------------------
      subroutine Get_surf_air_temp (self, surf_air_temp)
        implicit none
        real*8 , intent(out) :: surf_air_temp(:,:)
        type(t_metFields), intent(in) :: self
        surf_air_temp(:,:) = self%surf_air_temp(:,:)
        return
      end subroutine Get_surf_air_temp
!------------------------------------------------------------------------------
      subroutine Get_surf_rough (self, surf_rough)
        implicit none
        real*8 , intent(out) :: surf_rough(:,:)
        type(t_metFields), intent(in) :: self
        surf_rough(:,:) = self%surf_rough(:,:)
        return
      end subroutine Get_surf_rough
!------------------------------------------------------------------------------
      subroutine Get_tot_precip (self, tot_precip)
        implicit none
        real*8 , intent(out) :: tot_precip(:,:)
        type(t_metFields), intent(in) :: self
        tot_precip(:,:) = self%tot_precip(:,:)
        return
      end subroutine Get_tot_precip
!------------------------------------------------------------------------------
      subroutine Get_ustar (self, ustar)
        implicit none
        real*8 , intent(out) :: ustar(:,:)
        type(t_metFields), intent(in) :: self
        ustar(:,:) = self%ustar(:,:)
        return
      end subroutine Get_ustar
!------------------------------------------------------------------------------
      subroutine Set_saldif (self, saldif)
        implicit none
        real*8 , intent(in) :: saldif(:,:)
        type(t_metFields), intent(inOut) :: self
        self%saldif(:,:) = saldif(:,:)
        return
      end subroutine Set_saldif
!------------------------------------------------------------------------------
      subroutine Get_saldif (self, saldif)
        implicit none
        real*8 , intent(out) :: saldif(:,:)
        type(t_metFields), intent(in) :: self
        saldif(:,:) = self%saldif(:,:)
        return
      end subroutine Get_saldif
!------------------------------------------------------------------------------
      subroutine Set_saldir (self, saldir)
        implicit none
        real*8 , intent(in) :: saldir(:,:)
        type(t_metFields), intent(inOut) :: self
        self%saldir(:,:) = saldir(:,:)
        return
      end subroutine Set_saldir
!------------------------------------------------------------------------------
      subroutine Get_saldir (self, saldir)
        implicit none
        real*8 , intent(out) :: saldir(:,:)
        type(t_metFields), intent(in) :: self
        saldir(:,:) = self%saldir(:,:)
        return
      end subroutine Get_saldir
!------------------------------------------------------------------------------
      subroutine Set_sasdif (self, sasdif)
        implicit none
        real*8 , intent(in) :: sasdif(:,:)
        type(t_metFields), intent(inOut) :: self
        self%sasdif(:,:) = sasdif(:,:)
        return
      end subroutine Set_sasdif
!------------------------------------------------------------------------------
      subroutine Get_sasdif (self, sasdif)
        implicit none
        real*8 , intent(out) :: sasdif(:,:)
        type(t_metFields), intent(in) :: self
        sasdif(:,:) = self%sasdif(:,:)
        return
      end subroutine Get_sasdif
!------------------------------------------------------------------------------
      subroutine Set_sasdir (self, sasdir)
        implicit none
        real*8 , intent(in) :: sasdir(:,:)
        type(t_metFields), intent(inOut) :: self
        self%sasdir(:,:) = sasdir(:,:)
        return
      end subroutine Set_sasdir
!------------------------------------------------------------------------------
      subroutine Get_sasdir (self, sasdir)
        implicit none
        real*8 , intent(out) :: sasdir(:,:)
        type(t_metFields), intent(in) :: self
        sasdir(:,:) = self%sasdir(:,:)
        return
      end subroutine Get_sasdir
!------------------------------------------------------------------------------
      subroutine Set_surf_alb_uv (self, surf_alb_uv)
        implicit none
        real*8 , intent(in) :: surf_alb_uv(:,:)
        type(t_metFields), intent(inOut) :: self
        self%surf_alb_uv(:,:) = surf_alb_uv(:,:)
        return
      end subroutine Set_surf_alb_uv
!------------------------------------------------------------------------------
      subroutine Get_surf_alb_uv (self, surf_alb_uv)
        implicit none
        real*8 , intent(out) :: surf_alb_uv(:,:)
        type(t_metFields), intent(in) :: self
        surf_alb_uv(:,:) = self%surf_alb_uv(:,:)
        return
      end subroutine Get_surf_alb_uv
!------------------------------------------------------------------------------
      subroutine Get_fracCloudCover (self, fracCloudCover)
        implicit none
        real*8 , intent(out) :: fracCloudCover(:,:)
        type(t_metFields), intent(in) :: self
        fracCloudCover(:,:) = self%fracCloudCover(:,:)
        return
      end subroutine Get_fracCloudCover
!------------------------------------------------------------------------------
      subroutine Set_fracCloudCover (self, fracCloudCover)
        implicit none
        real*8 , intent(in) :: fracCloudCover(:,:)
        type(t_metFields), intent(inOut) :: self
        self%fracCloudCover(:,:) = fracCloudCover(:,:)
        return
      end subroutine Set_fracCloudCover
!------------------------------------------------------------------------------
      subroutine Get_totalCloudFraction (self, totalCloudFraction)
        implicit none
        real*8 , intent(out) :: totalCloudFraction(:,:,:)
        type(t_metFields), intent(in) :: self
        totalCloudFraction(:,:,:) = self%totalCloudFraction(:,:,:)
        return
      end subroutine Get_totalCloudFraction
!------------------------------------------------------------------------------
      subroutine Set_totalCloudFraction (self, totalCloudFraction)
        implicit none
        real*8 , intent(in) :: totalCloudFraction(:,:,:)
        type(t_metFields), intent(inOut) :: self
        self%totalCloudFraction(:,:,:) = totalCloudFraction(:,:,:)
        return
      end subroutine Set_totalCloudFraction
!------------------------------------------------------------------------------
      subroutine Set_kzz (self, kzz)
        implicit none
        real*8 , intent(in) :: kzz(:,:,:)
        type(t_metFields), intent(inOut) :: self
        self%kzz(:,:,:) = kzz(:,:,:)
        return
      end subroutine Set_kzz
!------------------------------------------------------------------------------
      subroutine Get_kzz (self, kzz)
        implicit none
        real*8 , intent(out) :: kzz(:,:,:)
        type(t_metFields), intent(in) :: self
        kzz(:,:,:) = self%kzz(:,:,:)
        return
      end subroutine Get_kzz
!------------------------------------------------------------------------------
      subroutine Get_max_cloud (self, max_cloud)
        implicit none
        real*8 , intent(out) :: max_cloud(:,:,:)
        type(t_metFields), intent(in) :: self
        max_cloud(:,:,:) = self%max_cloud(:,:,:)
        return
      end subroutine Get_max_cloud
!------------------------------------------------------------------------------
      subroutine Get_ran_cloud (self, ran_cloud)
        implicit none
        real*8 , intent(out) :: ran_cloud(:,:,:)
        type(t_metFields), intent(in) :: self
        ran_cloud(:,:,:) = self%ran_cloud(:,:,:)
        return
      end subroutine Get_ran_cloud
!------------------------------------------------------------------------------
      subroutine Set_tau_cloud (self, tau_cloud)
        implicit none
        real*8 , intent(in) :: tau_cloud(:,:,:)
        type(t_metFields), intent(inOut) :: self
        self%tau_cloud(:,:,:) = tau_cloud(:,:,:)
        return
      end subroutine Set_tau_cloud
!------------------------------------------------------------------------------
      subroutine Get_tau_cloud (self, tau_cloud)
        implicit none
        real*8 , intent(out) :: tau_cloud(:,:,:)
        type(t_metFields), intent(in) :: self
        tau_cloud(:,:,:) = self%tau_cloud(:,:,:)
        return
      end subroutine Get_tau_cloud
!------------------------------------------------------------------------------
      subroutine Get_cmf (self, cmf)
        implicit none
        real*8 , intent(out) :: cmf(:,:,:)
        type(t_metFields), intent(in) :: self
        cmf(:,:,:) = self%cmf(:,:,:)
        return
      end subroutine Get_cmf
!------------------------------------------------------------------------------
      subroutine Get_dtrn (self, dtrn)
        implicit none
        real*8 , intent(out) :: dtrn(:,:,:)
        type(t_metFields), intent(in) :: self
        dtrn(:,:,:) = self%dtrn(:,:,:)
        return
      end subroutine Get_dtrn
!------------------------------------------------------------------------------
      subroutine Get_ed (self, ed)
        implicit none
        real*8 , intent(out) :: ed(:,:,:)
        type(t_metFields), intent(in) :: self
        ed(:,:,:) = self%ed(:,:,:)
        return
      end subroutine Get_ed
!------------------------------------------------------------------------------
      subroutine Get_eu (self, eu)
        implicit none
        real*8 , intent(out) :: eu(:,:,:)
        type(t_metFields), intent(in) :: self
        eu(:,:,:) = self%eu(:,:,:)
        return
      end subroutine Get_eu
!------------------------------------------------------------------------------
      subroutine Get_md (self, md)
        implicit none
        real*8 , intent(out) :: md(:,:,:)
        type(t_metFields), intent(in) :: self
        md(:,:,:) = self%md(:,:,:)
        return
      end subroutine Get_md
!------------------------------------------------------------------------------
      subroutine Get_u10m (self, u10m)
        implicit none
        real*8 , intent(out) :: u10m(:,:)
        type(t_metFields), intent(in) :: self
        u10m(:,:) = self%u10m(:,:)
        return
      end subroutine Get_u10m
!------------------------------------------------------------------------------
      subroutine Get_v10m (self, v10m)
        implicit none
        real*8 , intent(out) :: v10m(:,:)
        type(t_metFields), intent(in) :: self
        v10m(:,:) = self%v10m(:,:)
        return
      end subroutine Get_v10m
!------------------------------------------------------------------------------
      subroutine Get_gwet (self, gwet)
        implicit none
        real*8 , intent(out) :: gwet(:,:)
        type(t_metFields), intent(in) :: self
        gwet(:,:) = self%gwet(:,:)
        return
      end subroutine Get_gwet
!------------------------------------------------------------------------------
      subroutine Get_pardif (self, pardif)
        implicit none
        real*8 , intent(out) :: pardif(:,:)
        type(t_metFields), intent(in) :: self
        pardif(:,:) = self%pardif(:,:)
        return
      end subroutine Get_pardif
!------------------------------------------------------------------------------
      subroutine Get_pardir (self, pardir)
        implicit none
        real*8 , intent(out) :: pardir(:,:)
        type(t_metFields), intent(in) :: self
        pardir(:,:) = self%pardir(:,:)
        return
      end subroutine Get_pardir
!------------------------------------------------------------------------------
      subroutine Get_zmdu (self, zmdu)
        implicit none
        real*8 , intent(out) :: zmdu(:,:,:)
        type(t_metFields), intent(in) :: self
        zmdu(:,:,:) = self%zmdu(:,:,:)
        return
      end subroutine Get_zmdu
!------------------------------------------------------------------------------
      subroutine Get_zmeu (self, zmeu)
        implicit none
        real*8 , intent(out) :: zmeu(:,:,:)
        type(t_metFields), intent(in) :: self
        zmeu(:,:,:) = self%zmeu(:,:,:)
        return
      end subroutine Get_zmeu
!------------------------------------------------------------------------------
      subroutine Get_zmmd (self, zmmd)
        implicit none
        real*8 , intent(out) :: zmmd(:,:,:)
        type(t_metFields), intent(in) :: self
        zmmd(:,:,:) = self%zmmd(:,:,:)
        return
      end subroutine Get_zmmd
!------------------------------------------------------------------------------
      subroutine Get_zmmu (self, zmmu)
        implicit none
        real*8 , intent(out) :: zmmu(:,:,:)
        type(t_metFields), intent(in) :: self
        zmmu(:,:,:) = self%zmmu(:,:,:)
        return
      end subroutine Get_zmmu
!------------------------------------------------------------------------------
      subroutine Get_hkdu (self, hkdu)
        implicit none
        real*8 , intent(out) :: hkdu(:,:,:)
        type(t_metFields), intent(in) :: self
        hkdu(:,:,:) = self%hkdu(:,:,:)
        return
      end subroutine Get_hkdu
!------------------------------------------------------------------------------
      subroutine Get_hkeu (self, hkeu)
        implicit none
        real*8 , intent(out) :: hkeu(:,:,:)
        type(t_metFields), intent(in) :: self
        hkeu(:,:,:) = self%hkeu(:,:,:)
        return
      end subroutine Get_hkeu
!------------------------------------------------------------------------------
      subroutine Get_hkmu (self, hkmu)
        implicit none
        real*8 , intent(out) :: hkmu(:,:,:)
        type(t_metFields), intent(in) :: self
        hkmu(:,:,:) = self%hkmu(:,:,:)
        return
      end subroutine Get_hkmu
!------------------------------------------------------------------------------
      subroutine Get_clwc (self, clwc)
        implicit none
        real*8 , intent(out) :: clwc(:,:,:)
        type(t_metFields), intent(in) :: self
        clwc(:,:,:) = self%clwc(:,:,:)
        return
      end subroutine Get_clwc
!------------------------------------------------------------------------------
      subroutine Get_rain (self, rain)
        implicit none
        real*8 , intent(out) :: rain(:,:,:)
        type(t_metFields), intent(in) :: self
        rain(:,:,:) = self%rain(:,:,:)
        return
      end subroutine Get_rain
!------------------------------------------------------------------------------
      subroutine Get_taucli (self, taucli)
        implicit none
        real*8 , intent(out) :: taucli(:,:,:)
        type(t_metFields), intent(in) :: self
        taucli(:,:,:) = self%taucli(:,:,:)
        return
      end subroutine Get_taucli
!------------------------------------------------------------------------------
      subroutine Get_tauclw (self, tauclw)
        implicit none
        real*8 , intent(out) :: tauclw(:,:,:)
        type(t_metFields), intent(in) :: self
        tauclw(:,:,:) = self%tauclw(:,:,:)
        return
      end subroutine Get_tauclw
!------------------------------------------------------------------------------
      subroutine Get_rain_ls (self, rain_ls)
        implicit none
        real*8 , intent(out) :: rain_ls(:,:,:)
        type(t_metFields), intent(in) :: self
        rain_ls(:,:,:) = self%rain_ls(:,:,:)
        return
      end subroutine Get_rain_ls
!------------------------------------------------------------------------------
      subroutine Get_rain_zm (self, rain_zm)
        implicit none
        real*8 , intent(out) :: rain_zm(:,:,:)
        type(t_metFields), intent(in) :: self
        rain_zm(:,:,:) = self%rain_zm(:,:,:)
        return
      end subroutine Get_rain_zm
!------------------------------------------------------------------------------
      subroutine Get_rain_hk (self, rain_hk)
        implicit none
        real*8 , intent(out) :: rain_hk(:,:,:)
        type(t_metFields), intent(in) :: self
        rain_hk(:,:,:) = self%rain_hk(:,:,:)
        return
      end subroutine Get_rain_hk
!------------------------------------------------------------------------------
      subroutine Get_aconv (self, aconv)
        implicit none
        real*8 , intent(out) :: aconv(:,:,:)
        type(t_metFields), intent(in) :: self
        aconv(:,:,:) = self%aconv(:,:,:)
        return
      end subroutine Get_aconv
!------------------------------------------------------------------------------
      subroutine Get_astrat (self, astrat)
        implicit none
        real*8 , intent(out) :: astrat(:,:,:)
        type(t_metFields), intent(in) :: self
        astrat(:,:,:) = self%astrat(:,:,:)
        return
      end subroutine Get_astrat
!------------------------------------------------------------------------------
      subroutine Get_cldmas (self, cldmas)
        implicit none
        real*8 , intent(out) :: cldmas(:,:,:)
        type(t_metFields), intent(in) :: self
        cldmas(:,:,:) = self%cldmas(:,:,:)
        return
      end subroutine Get_cldmas
!------------------------------------------------------------------------------
      subroutine Get_omega (self, omega)
        implicit none
        real*8 , intent(out) :: omega(:,:,:)
        type(t_metFields), intent(in) :: self
        omega(:,:,:) = self%omega(:,:,:)
        return
      end subroutine Get_omega
!------------------------------------------------------------------------------
      subroutine Get_relhum (self, relhum)
        implicit none
        real*8 , intent(out) :: relhum(:,:,:)
        type(t_metFields), intent(in) :: self
        relhum(:,:,:) = self%relhum(:,:,:)
        return
      end subroutine Get_relhum
!------------------------------------------------------------------------------
      subroutine Get_rhclear (self, rhclear)
        implicit none
        real*8 , intent(out) :: rhclear(:,:,:)
        type(t_metFields), intent(in) :: self
        rhclear(:,:,:) = self%rhclear(:,:,:)
        return
      end subroutine Get_rhclear
!------------------------------------------------------------------------------
      subroutine Get_zmed (self, zmed)
        implicit none
        real*8 , intent(out) :: zmed(:,:,:)
        type(t_metFields), intent(in) :: self
        zmed(:,:,:) = self%zmed(:,:,:)
        return
      end subroutine Get_zmed
!------------------------------------------------------------------------------
      subroutine Get_humidity (self, humidity)
        implicit none
        real*8 , intent(out) :: humidity(:,:,:)
        type(t_metFields), intent(in) :: self
        humidity(:,:,:) = self%humidity(:,:,:)
        return
      end subroutine Get_humidity
!------------------------------------------------------------------------------
      subroutine Set_moistq (self, moistq)
        implicit none
        real*8 , intent(in) :: moistq(:,:,:)
        type(t_metFields), intent(inOut) :: self
        self%moistq(:,:,:) = moistq(:,:,:)
        return
      end subroutine Set_moistq
!------------------------------------------------------------------------------
      subroutine Get_moistq (self, moistq)
        implicit none
        real*8 , intent(out) :: moistq(:,:,:)
        type(t_metFields), intent(in) :: self
        moistq(:,:,:) = self%moistq(:,:,:)
        return
      end subroutine Get_moistq
!------------------------------------------------------------------------------
      end module GmiMetFieldsControl_mod
!
