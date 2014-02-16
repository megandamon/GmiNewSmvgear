!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlColumn_mod
!
! !INTERFACE:
!
module GmiControlColumn_mod
!
      use GmiInterpolation_mod, only : Interp
      use GmiFileOperations_mod, only : makeOutfileName
      use m_netcdf_io_close     , only : Nccl, Nccl_Noerr
      use m_netcdf_io_write     , only : Ncwr_2d_Int, Ncwr_3d, Ncwr_2d
      use m_netcdf_io_write     , only : Ncwr_1d, Ncwr_Scal, Ncwr_2d_Char
      use m_netcdf_io_create    , only : Ncdo_Sync, Nccr_Wr
      use m_netcdf_io_define    , only : NcDef_glob_attributes, NcDef_dimension
      use m_netcdf_io_define    , only : NcDef_var_attributes, NcDef_variable
      use m_netcdf_io_define    , only : NcSetFill, NcEnd_def
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration
!
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiChemistryMethod_mod, only : t_Chemistry, Get_const_labels
!
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_rootProc,        &
     &       Get_procID
      use GmiDomainDecomposition_mod, only : Get_latdeg, Get_londeg
      use GmiGrid_mod, only : t_gmiGrid, Get_i1, Get_i2, Get_ju1, Get_j2,      &
     &       Get_i2_gl, Get_j2_gl, Get_k1, Get_k2, Get_numSpecies, Get_i1_gl,  &
     &       Get_ju1_gl, Get_ilo, Get_ihi, Get_julo, Get_jhi
      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Get_problem_name,    &
     &       Get_pr_diag, Get_hdf_dim_name, Get_hdr_var_name, Get_spc_dim_name,&
     &       Get_col_diag_site, Get_pr_col_diag, Get_col_diag_num,             &
     &       Get_col_diag_pres_num, Get_col_diag_species_num,                  &
     &       Get_col_diag_species, Get_col_diag_period, Get_col_diag_pres,     &
     &       Get_col_diag_lat_lon, Get_rec_dim_name, Set_col_diag_lat_lon
!
      use GmiTimeControl_mod, only : t_GmiClock, GmiSplitDateTime,            &
     &       Get_curGmiDate, Get_curGmiTime, Get_gmiSeconds, Get_gmiTimeStep, &
     &       Get_numTimeSteps
!
      use GmiMetFieldsControl_mod, only : t_metFields, Get_pt, Get_met_opt,    &
     &       Get_metdata_name, Get_ai, Get_bi, Get_am, Get_bm, Get_humidity,   &
     &       Get_pctm2, Get_kel
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: initializeOutputColumn, controlOutputColumn, finalizeOutputColumn
!
#     include "GmiParameters.h"
!
!
! !DESCRIPTION:
!  Provides routines to manipulate the file <problem\_name>.column.nc
!
      ! Grid information of the variables to be written out
      integer,          save :: i1, i2     ! longitude
      integer,          save :: i1_gl, i2_gl
      integer,          save :: ilo_gl, ihi_gl, ilo, ihi
      integer,          save :: ju1        ! latitude
      integer,          save :: j2
      integer,          save :: ju1_gl, j2_gl, julo, jhi
      integer,          save :: k1, k2
      integer,          save :: numVert
!
      logical,          save :: pr_diag
      integer,          save :: procID
      integer,          save :: numSpecies
      real*8 ,          save :: ri2_gl, rj2_gl
      real*8 ,          save :: half_lat, half_lon
!
      integer,            save :: rnum_out_col
!
      integer, pointer,   save :: ncid_col(:) => null()
      character (len=24), save :: col_diag_site(MAX_COL_DIAG_SITES)
      logical,            save :: pr_col_diag
      integer,            save :: col_diag_num
      integer,            save :: col_diag_pres_num
      integer,            save :: col_diag_species_num
      integer,            save :: col_diag_species(MAX_NUM_CONST_GIO)
      real*8 ,            save :: col_diag_period
      real*8 ,            save :: col_diag_pres(MAX_COL_DIAG_PRES)
      real*8 ,            save :: col_diag_lat_lon(2, MAX_COL_DIAG_SITES)
!      real*8 , pointer,   save :: latdeg(:) => null()
!      real*8 , pointer,   save :: londeg(:) => null()
!
      character (len=MAX_LENGTH_VAR_NAME), save :: hdr_var_name, hdf_dim_name, spc_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: rec_dim_name
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!  Jan 20 2012: Megan Damon - changes for station output on all model levels.
!
!EOP
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initializeOutputColumn
!
! !INTERFACE:
!
      subroutine initializeOutputColumn(gmiGrid, gmiDomain, Diagnostics, &
     &                                  Chemistry, metfields)
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_gmiGrid    ), intent(in) :: gmiGrid
      type(t_gmiDomain  ), intent(in) :: gmiDomain
      type(t_Chemistry  ), intent(in) :: Chemistry
      type(t_Diagnostics), intent(in) :: Diagnostics
      type(t_metFields)  , intent(in) :: metFields
!
! !LOCAL VARIABLES:
      integer :: in1, isite
      character (len=MAX_LENGTH_SPECIES_NAME), pointer :: const_labels(:)
      real*8 , allocatable :: latdeg(:), londeg(:)
      character (len=MAX_LENGTH_FILE_NAME) :: fname
      character (len=MAX_LENGTH_FILE_NAME) :: problem_name
!
!EOP
!-------------------------------------------------------------------------------
!BOC
      call Get_procID (gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)
!
      if (pr_diag) Write(6,*) 'initializeOutputColumn called by ', procID
!
      call Get_pr_col_diag         (Diagnostics, pr_col_diag)
!
      call Get_problem_name        (Diagnostics, problem_name)
      call Get_hdf_dim_name        (Diagnostics, hdf_dim_name)
      call Get_spc_dim_name        (Diagnostics, spc_dim_name)
      call Get_rec_dim_name        (Diagnostics, rec_dim_name)
      call Get_hdr_var_name        (Diagnostics, hdr_var_name)
      call Get_col_diag_num        (Diagnostics, col_diag_num)
      call Get_col_diag_site       (Diagnostics, col_diag_site)
!
!     MRD: station levels now on all model levels
!      call Get_col_diag_pres       (Diagnostics, col_diag_pres)
      call Get_col_diag_period     (Diagnostics, col_diag_period)
      call Get_col_diag_lat_lon    (Diagnostics, col_diag_lat_lon)
      call Get_col_diag_species    (Diagnostics, col_diag_species)
!     MRD: station levels now on all model levels
!      call Get_col_diag_pres_num   (Diagnostics, col_diag_pres_num)
      call Get_col_diag_species_num(Diagnostics, col_diag_species_num)
!
      allocate(ncid_col(col_diag_num))
      ncid_col = 0
!
      call Get_ilo   (gmiGrid, ilo  )
      call Get_ihi   (gmiGrid, ihi  )
      call Get_julo  (gmiGrid, julo )
      call Get_jhi   (gmiGrid, jhi  )
      call Get_i1    (gmiGrid, i1   )
      call Get_i2    (gmiGrid, i2   )
      call Get_ju1   (gmiGrid, ju1  )
      call Get_j2    (gmiGrid, j2   )
      call Get_k1    (gmiGrid, k1   )
      call Get_k2    (gmiGrid, k2   )
      call Get_i1_gl (gmiGrid, i1_gl )
      call Get_i2_gl (gmiGrid, i2_gl )
      call Get_j2_gl (gmiGrid, j2_gl )
      call Get_ju1_gl(gmiGrid, ju1_gl)
      call Get_numSpecies (gmiGrid, numSpecies)
!
      numvert = k2 - k1 + 1
!
      ri2_gl = i2_gl
      rj2_gl = j2_gl
!
      half_lat = 180.0d0 / (rj2_gl - 1.0d0) * 0.5d0
      half_lon = 360.0d0 /  ri2_gl          * 0.5d0
!
      allocate(latdeg(ju1_gl:j2_gl))
      call Get_latdeg(gmiDomain, latdeg)
!
      allocate(londeg(i1_gl:i2_gl))
      call Get_londeg(gmiDomain, londeg)
!
      allocate(const_labels(numSpecies))
      call Get_const_labels(Chemistry, const_labels)
!
      do isite = 1, col_diag_num
!
         !-------------------------------------------------------------
         ! Longitude is greater than maximum longitude plus half a grid
         ! box, so make the longitude a small negative longitude
         ! (west longitude).
         !-------------------------------------------------------------
         if (col_diag_lat_lon(2,isite) >= (londeg(i2_gl) + half_lon)) then
             col_diag_lat_lon(2,isite) = col_diag_lat_lon(2,isite) - 360.0d0
         end if
!
         if ((col_diag_lat_lon(1,isite) >= (latdeg(ju1) - half_lat)) .and.  &
     &       (col_diag_lat_lon(1,isite) <  (latdeg(j2)  + half_lat)) .and.  &
     &       (col_diag_lat_lon(2,isite) >= (londeg(i1)  - half_lon)) .and.  &
     &       (col_diag_lat_lon(2,isite) <  (londeg(i2)  + half_lon))) then
!
            in1 = Len_Trim (col_diag_site(isite))
!
            call makeOutfileName (fname, '_' // &
     &           col_diag_site(isite)(1:in1) // '.profile.nc', problem_name)
!
            ! Create the file and assign a file identifier
            call Nccr_Wr (ncid_col(isite), fname)
!
            ! Define the variables in the file
            call Define_Netcdf_Out_Col (ncid_col(isite), isite)
!
            ! Write header data
            call Write_Netcdf_Hdr_Col(ncid_col(isite), isite, col_diag_lat_lon, const_labels, &
                 metfields, gmiDomain)
!
            call Ncdo_Sync (ncid_col(isite))
         end if
      end do
!
      ! Initialize the counter for the number of records
      rnum_out_col = 1
!
      return
!
      end subroutine initializeOutputColumn
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: controlOutputColumn
!
! !INTERFACE:
!
      subroutine controlOutputColumn (gmiDomain, SpeciesConcentration, &
     &                  metFields, gmiClock)
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_gmiClock )           , intent(in) :: gmiClock
      type(t_gmiDomain)           , intent(in) :: gmiDomain
      type(t_metFields)           , intent(in) :: metFields
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !LOCAL VARIABLES:
      logical :: time_for_col
      integer :: ndt, num_time_steps, nhms, nymd
      real*8  :: gmi_sec, tstep
      integer :: isite
      real*8 , allocatable :: latdeg(:), londeg(:)
      type(t_GmiArrayBundle), pointer :: concentration(:)
!
!EOP
!-------------------------------------------------------------------------------
!BOC
     if (pr_diag) Write (6,*) 'controlOutputColumn called by ', procID
!
      call Get_curGmiDate(gmiClock, nymd)
      call Get_curGmitime(gmiClock, nhms)
      call Get_numTimeSteps(gmiClock, num_time_steps)
      call Get_gmiSeconds(gmiClock, gmi_sec)
      call Get_gmiTimeStep(gmiClock, tstep)
!
      ndt = Nint(tstep)
!
     if (Mod (Nint (gmi_sec), Nint (col_diag_period)) < ndt) then
       time_for_col = .true.
     else
       time_for_col = .false.
     end if
!
!
!    =================
     if (time_for_col) then
!    =================
!
       allocate(latdeg(ju1_gl:j2_gl))
       call Get_latdeg(gmiDomain, latdeg)
!
       allocate(londeg(i1_gl:i2_gl))
       call Get_londeg(gmiDomain, londeg)
!
       do isite = 1, col_diag_num
!
         if ((col_diag_lat_lon(1,isite) >= (latdeg(ju1) - half_lat)) .and.  &
     &       (col_diag_lat_lon(1,isite) <  (latdeg(j2)  + half_lat)) .and.  &
     &       (col_diag_lat_lon(2,isite) >= (londeg(i1) - half_lon)) .and.  &
     &       (col_diag_lat_lon(2,isite) <  (londeg(i2) + half_lon))) then
!
           call Get_concentration(SpeciesConcentration, concentration)
!
           call Write_Netcdf_Col (isite, metFields, concentration, &
     &                 nymd, nhms, gmi_sec, londeg, latdeg)
!
         end if
!
       end do
!
       rnum_out_col = rnum_out_col + 1
!
!     ======
      end if
!     ======
!
      return
!
      end subroutine controlOutputColumn
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finalizeOutputColumn
!
! !INTERFACE:
!
      subroutine finalizeOutputColumn()
!
      implicit none
!
! !DESCRIPTION:
! Closes all the netCDF files created for the column diagnostics.
!
! !LOCAL VARIABLES:
      integer :: isite
!
!EOP
!-------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'finalizeOutputColumn called by ', procID
!
      do isite = 1, col_diag_num
         if (ncid_col(isite) > 0) call Nccl_Noerr (ncid_col(isite))
      end do
!
      return
!
      end subroutine finalizeOutputColumn
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Define_Netcdf_Out_Col
!
! !INTERFACE:
!
      subroutine Define_Netcdf_Out_Col (ncid, isite)
!
! !USES:
       use m_ncGeneralOpsOutput, only : Define_Netcdf_Out_Gen
!
      implicit none
!
#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
      integer, intent(in) :: ncid  ! netCDF file id
      integer, intent(in) :: isite ! column site index
!
! !DESCRIPTION:
!   Makes the necessary definitions for a netCDF column diagnostic output file.
!
! !LOCAL VARIABLES:
      integer :: collatid, collonid
      integer :: ierr
      integer :: nchr, nhdf
      integer :: omode
      integer :: varid
      integer :: chrd(1), hdfd(1)
      integer :: prsd(1), prsd2(1), scal(1)
      integer :: spcd(1), recd(1)
      integer :: var1(1), var2(2), var3(3)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
!     ------------------
!     Define dimensions.
!     ------------------
!
      nhdf = NETCDF_HDF
      nchr = MAX_LENGTH_SPECIES_NAME
!
      call NcDef_dimension(ncid, 'scalar', 1, scal(1))
!                                 ---------
      call NcDef_dimension(ncid, 'chr_dim', nchr, chrd(1))
!                                 ---------
      call NcDef_dimension(ncid, 'pressure', k2, prsd(1))
!
      call NcDef_dimension(ncid, 'pressure_edges', k2+1, prsd2(1))
!                                 ----------
      call NcDef_dimension (ncid, spc_dim_name, col_diag_species_num, spcd(1))
!                                 ------------
!
      call NcDef_dimension(ncid, hdf_dim_name, nhdf, hdfd(1))
!                                 ------------
!
      call NcDef_dimension(ncid, rec_dim_name, NF_UNLIMITED, recd(1))
!                                 ------------
!     -----------------------------------------
!     Define variables and variable attributes.
!     -----------------------------------------
!
      call NcDef_variable (ncid, hdf_dim_name, NF_INT, 1, hdfd, varid)
!                                 ------------
      call NcDef_var_attributes (ncid, varid, 'long_name', 'Header index')
      call NcDef_var_attributes (ncid, varid, 'units', 'unitless')
!
      var2(:) = (/ hdfd(1), recd(1) /)
      call NcDef_variable (ncid, hdr_var_name, NF_INT, 2, var2, varid)
!                               !-------------
      call NcDef_var_attributes (ncid, varid, 'long_name', 'Header')
      call NcDef_var_attributes (ncid, varid, 'units', 'gmi_sec, nymd, nhms')
!
      var1(1) = scal(1)
      call NcDef_variable (ncid, 'col_latitude', NF_FLOAT, 1, var1, collatid)
!                                 --------------
      call NcDef_variable (ncid, 'col_longitude', NF_FLOAT, 1, var1,collonid)
!                                 ---------------
      call NcDef_variable (ncid, spc_dim_name, NF_FLOAT, 1, spcd, varid)
!                                 ------------
      call NcDef_var_attributes (ncid, varid, 'long_name', 'Species index')
      call NcDef_var_attributes (ncid, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid, varid, 'coord_labels', 'const_labels')
      call NcDef_var_attributes (ncid, varid, 'selection_category', 'NULL')
!
      var2(:) = (/ chrd(1), spcd(1) /)
      call NcDef_variable (ncid, 'const_labels', NF_CHAR, 2, var2, varid)
!                                 --------------
      call NcDef_var_attributes (ncid, varid, 'long_name', 'Species name')
      call NcDef_var_attributes (ncid, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid, varid, 'selection_category', 'NULL')
!
      call NcDef_variable (ncid, 'pressure', NF_FLOAT, 1, prsd, varid)
!                                 ----------
      call NcDef_var_attributes (ncid, varid, 'long_name', 'pressure')
      call NcDef_var_attributes (ncid, varid, 'units', 'mb')
!
      call NcDef_variable (ncid, 'am', NF_FLOAT, 1, prsd, varid)
!                                 --
      call NcDef_var_attributes (ncid, varid, 'long_name', 'Hybrid pressure term')
      call NcDef_var_attributes (ncid, varid, 'units', 'unitless')
!
      call NcDef_variable (ncid, 'bm', NF_FLOAT, 1, prsd, varid)
!                                 --
      call NcDef_var_attributes (ncid, varid, 'long_name', 'Hybrid sigma term')
      call NcDef_var_attributes (ncid, varid, 'units', 'unitless')
!
      call NcDef_variable (ncid, 'ai', NF_FLOAT, 1, prsd2, varid)
!                                 --
      call NcDef_var_attributes (ncid, varid, 'long_name', 'Hybrid pressure edge term')
      call NcDef_var_attributes (ncid, varid, 'units', 'unitless')
!
      call NcDef_variable (ncid, 'bi', NF_FLOAT, 1, prsd2, varid)
!                                 --
      call NcDef_var_attributes (ncid, varid, 'long_name', 'Hybrid sigma edge term')
      call NcDef_var_attributes (ncid, varid, 'units', 'unitless')
!
      call NcDef_variable (ncid, 'surf_pres', NF_FLOAT, 1, recd, varid)
!                                 -----------
      call NcDef_var_attributes (ncid, varid, 'long_name', 'surf_pres')
      call NcDef_var_attributes (ncid, varid, 'units', 'mb')
!
      var2(:) = (/ prsd(1), recd(1) /)
      call NcDef_variable (ncid, 'pot_temp', NF_FLOAT, 2, var2, varid)
!                                 ----------
!
      call NcDef_var_attributes (ncid, varid, 'long_name', 'pot_temp')
      call NcDef_var_attributes (ncid, varid, 'units', 'deg K')
!
      var2(:) = (/ prsd(1), recd(1) /)
      call NcDef_variable (ncid, 'humidity', NF_FLOAT, 2, var2, varid)
!                                 ----------
      call NcDef_var_attributes (ncid, varid, 'long_name', 'humidity')
      call NcDef_var_attributes (ncid, varid, 'units', 'g/kg')
!
      var3(:) = (/ prsd(1), spcd(1), recd(1) /)
      call NcDef_variable (ncid, 'const', NF_FLOAT, 3, var3, varid)
!                                 -------
      call NcDef_var_attributes (ncid, varid, 'long_name', 'Constituent')
      call NcDef_var_attributes (ncid, varid, 'units', 'volume mixing ratio')
!
      var2(:) = (/ spcd(1), recd(1) /)
      call NcDef_variable (ncid, 'const_surf', NF_FLOAT, 2, var2, varid)
!                                 ----------
      call NcDef_var_attributes (ncid, varid, 'long_name', 'Constituent at the surface')
      call NcDef_var_attributes (ncid, varid, 'units', 'volume mixing ratio')
!
!     -------------------------
!     Define global attributes.
!     -------------------------
!
      call NcDef_glob_attributes (ncid, 'title', 'Gmimod column diagnostic file')
!
      call NcSetFill (ncid, NF_NOFILL, omode)
!
      call NcEnd_def (ncid)
!
      return
!
      end subroutine Define_Netcdf_Out_Col
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Hdr_Col
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Hdr_Col (ncid, isite, col_diag_lat_lon, const_labels, &
           metfields, gmiDomain)
!
! !USES:
      use m_ncGeneralOpsOutput, only : Write_Netcdf_Hdr_Gen
!
      implicit none
!
! !INPUT PARAMETERS:
      integer,            intent(in) :: ncid
      integer,            intent(in) :: isite
      character (len=MAX_LENGTH_SPECIES_NAME), intent(in) :: const_labels(1:)
      real*8 ,            intent(in) :: col_diag_lat_lon(2, MAX_COL_DIAG_SITES)
      type (t_metFields), intent(in) :: metFields
      type(t_gmiDomain),  intent(in) :: gmiDomain
!
! !DESCRIPTION:
!  Creates some header information for a column diagnostics netCDF output file
!  and writes it out.
!
! !LOCAL VARIABLES:
      integer :: ic, icx, il, ij, numVert
      integer :: cnt1d (1), cnt2d (2)
      integer :: strt1d(1), strt2d(2)
      real*8  :: spcdat(numSpecies)
      character (len=MAX_LENGTH_SPECIES_NAME) :: spclab(numSpecies)
      real*8 , allocatable :: am(:), bm(:), ai(:), bi(:)
      real*8 , allocatable :: pctm2(:,:)
      real*8 , allocatable :: latdeg(:), londeg(:)
      real*8 :: pt
      real*8  :: pres      (k1:k2)
!
!
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Hdr_Col called by ', procID
!
      call Ncwr_Scal (col_diag_lat_lon(1,isite), ncid, 'col_latitude')
      call Ncwr_Scal (col_diag_lat_lon(2,isite), ncid, 'col_longitude')
!
      allocate(am(k1:k2))
      allocate(bm(k1:k2))
      allocate(ai(k1:k2+1))
      allocate(bi(k1:k2+1))
!
      call Get_am(metFields, am)
      call Get_bm(metFields, bm)
      call Get_ai(metFields, ai)
      call Get_bi(metFields, bi)
!
      allocate(pctm2(ilo:ihi,julo:jhi))
      call Get_pctm2(metFields, pctm2)
      call Get_pt(metFields, pt)
!
      allocate(latdeg(ju1_gl:j2_gl))
      call Get_latdeg(gmiDomain, latdeg)
!
      allocate(londeg(i1_gl:i2_gl))
      call Get_londeg(gmiDomain, londeg)
!
!
      icx = 0
!
      do ic = 1, numSpecies
!
        if (col_diag_species(ic) /= 0) then
          icx = icx + 1
          spcdat(icx) = ic
          spclab(icx) = const_labels(ic)
        end if
!
      end do
!
      pres(:) = (am(:) * pt)  + (bm(:) * 1000.0d0)
      strt1d(1) = 1
      cnt1d (1) = k2 - k1 + 1
!
      call Ncwr_1d (pres, ncid, 'pressure', strt1d, cnt1d)
      call Ncwr_1d (am, ncid, 'am', strt1d, cnt1d)
      call Ncwr_1d (bm, ncid, 'bm', strt1d, cnt1d)
      call Ncwr_1d (ai, ncid, 'ai', strt1d, cnt1d + 1)
      call Ncwr_1d (bi, ncid, 'bi', strt1d, cnt1d + 1)
!
      strt1d(1) = 1
      cnt1d (1) = col_diag_species_num
!
      call Ncwr_1d (spcdat, ncid, spc_dim_name, strt1d, cnt1d)
!
!     ---------------
!     Species labels.
!     ---------------
!
      strt2d(:) = (/ 1, 1 /)
      cnt2d (:) = (/ MAX_LENGTH_SPECIES_NAME, col_diag_species_num /)
!
      call Ncwr_2d_Char(spclab, ncid, 'const_labels', strt2d, cnt2d)
!
      return
!
      end subroutine Write_Netcdf_Hdr_Col
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Col
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Col (isite, metFields, concentration, &
     &                 nymd, nhms, gmi_sec, londeg, latdeg)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: nymd, nhms
      real*8 , intent(in) :: gmi_sec
      integer,                 intent(in) :: isite
      real*8 ,                 intent(in) :: londeg(i1_gl:), latdeg(ju1_gl:)
      type (t_metFields     ), intent(in) :: metFields
      type (t_GmiArrayBundle), intent(in) :: concentration(numSpecies)
!
! !DESCRIPTION:
!  Writes the col. diag. NetCDF output.
!
! !LOCAL VARIABLES:
      integer :: il, ij, ic
      integer :: in1
      integer :: nc
      integer :: ncid
      integer :: hdr(NETCDF_HDF)
      integer :: cnt1d (1), cnt2d (2), cnt3d (3)
      integer :: strt1d(1), strt2d(2), strt3d(3)
      integer :: count2d(2)
      integer :: start2d(2)
      real*8  :: rtmp(1)
      real*8  :: const_surf(1)
      real*8  :: col_diag_data(MAX_COL_DIAG_PRES)
!      real*8  :: theta_data   (MAX_COL_DIAG_PRES)
      real*8  :: humidity_data(MAX_COL_DIAG_PRES)
      real*8  :: temp_col_pres(MAX_COL_DIAG_PRES)
      real*8  :: pres      (k1:k2)
      real*8  :: theta_grid(k1:k2)
      real*8 , allocatable :: am(:), bm(:), ai(:), bi(:)
      real*8 , allocatable :: pctm2(:,:)
      real*8 , allocatable :: humidity(:,:,:), kel(:,:,:)
      real*8 :: pt
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Col called by ', procID
!
      ncid = ncid_col(isite)
!
      strt2d(:) = (/ 1, rnum_out_col /)
      cnt2d (:) = (/ NETCDF_HDF, 1 /)
!
      hdr(1) = Nint (gmi_sec)
      hdr(2) = nymd
      hdr(3) = nhms
!
      call Ncwr_2d_Int (hdr, ncid, hdr_var_name, strt2d, cnt2d)
!
      strt1d(:) = (/ rnum_out_col /)
      strt2d(:) = (/ 1, rnum_out_col /)
!
      cnt1d (:) = (/ 1 /)
      cnt2d (:) = (/ k2, 1 /)
      cnt3d (:) = (/ k2, 1, 1 /)
!
      start2d(:) = (/ 1, rnum_out_col /)
      count2d(:) = (/ 1, 1 /)
!
      allocate(ai(k1-1:k2))
      allocate(bi(k1-1:k2))
      allocate(am(k1:k2))
      allocate(bm(k1:k2))
!
      call Get_pt(metFields, pt)
      call Get_am(metFields, am)
      call Get_bm(metFields, bm)
      call Get_ai(metFields, ai)
      call Get_bi(metFields, bi)
!
      allocate(pctm2(ilo:ihi,julo:jhi))
      call Get_pctm2(metFields, pctm2)
!
      allocate(kel(ilo:ihi,julo:jhi,k1:k2))
      call Get_kel(metFields, kel)
!
      allocate(humidity(i1:i2,ju1:j2,k1:k2))
      call Get_humidity(metFields, humidity)
!
!     ================
      do ij = ju1, j2
        do il = i1, i2
!     ================
!
          if ((col_diag_lat_lon(1,isite) >= (latdeg(ij) - half_lat)) .and.  &
     &        (col_diag_lat_lon(1,isite) <  (latdeg(ij) + half_lat)) .and.  &
     &        (col_diag_lat_lon(2,isite) >= (londeg(il) - half_lon)) .and.  &
     &        (col_diag_lat_lon(2,isite) <  (londeg(il) + half_lon))) then
!
!           -------------------------------------------------
!           Grid box found for this site.
!
!           Note that the calls to Interp below do a vertical
!           interpolation to pressure surface.
!           -------------------------------------------------
!
            pres(:) = (am(k1:) * pt) + (bm(k1:) * pctm2(il,ij))
!
            rtmp(1) = pctm2(il,ij)
            call Ncwr_1d (rtmp, ncid, 'surf_pres', strt1d, cnt1d)
!
            theta_grid = kel(il,ij,:) * (1013.25d0 / pres(:))**0.286d0
!
!
!     MRD: We used to interpolate station data to user-specified pressure
!          We now output station data on all model levels.
!            temp_col_pres(1:col_diag_pres_num:1) = &
!     &                             col_diag_pres(col_diag_pres_num:1:-1)
!
!           ===========
!            call Interp  &
!           ===========
!     &        (pres(k2:k1:-1), theta_grid(k2:k1:-1), numVert,  &
!     &         col_diag_pres(col_diag_pres_num:1:-1), theta_data,  &
!     &         col_diag_pres_num)
!
!          WHERE (PCTM2(IL,IJ) < temp_col_pres(:)) THETA_DATA(:) = -999.d0
!
!           call Ncwr_2d (theta_data(col_diag_pres_num:1:-1), ncid, 'pot_temp',  &
!     &         strt2d, cnt2d)
!
            call Ncwr_2d (theta_grid(k1:k2), ncid, 'pot_temp',  &
     &         strt2d, cnt2d)
!
!     MRD: We used to interpolate station data to user-specified pressure
!          We now output station data on all model levels.
!
!           ===========
!            call Interp  &
!           ===========
!     &        (pres(k2:k1:-1), humidity(il,ij,k2:k1:-1), numVert,  &
!     &         col_diag_pres(col_diag_pres_num:1:-1), humidity_data,  &
!     &         col_diag_pres_num)
!
!         WHERE (PCTM2(IL,IJ) < temp_col_pres(:)) HUMIDITY_DATA(:) = -999.d0
!
!            call Ncwr_2d  &
!     &        (humidity(il,ij,k2:k1:-1), ncid, 'humidity',  &
!     &         strt2d, cnt2d)
!
            call Ncwr_2d  &
     &        (humidity(il,ij,k1:k2), ncid, 'humidity',  &
     &         strt2d, cnt2d)
!
            nc = 0
!
            do ic = 1, numSpecies
              if (col_diag_species(ic) /= 0) then
!
                 nc = nc + 1
!
                 ! constituents at the surface
                 start2d(:) = (/ nc, rnum_out_col /)
                 const_surf(1) = concentration(ic)%pArray3D(il,ij,k1)
                 call Ncwr_2d (const_surf, ncid, 'const_surf', start2d, count2d)
!
!     MRD: We used to interpolate station data to user-specified pressure
!          We now output station data on all model levels.
!               ===========
!                call Interp  &
!               ===========
!     &            (pres(k2:k1:-1), concentration(ic)%pArray3D(il,ij,k2:k1:-1), numVert,  &
!     &             col_diag_pres(col_diag_pres_num:1:-1),  &
!     &             col_diag_data, col_diag_pres_num)
!
 !         WHERE (PCTM2(IL,IJ) < temp_col_pres(:)) COL_DIAG_DATA(:) = -999.d0
!
                strt3d(:) = (/ 1, nc, rnum_out_col /)
!
!                call Ncwr_3d (concentration(ic)%pArray3D(il,ij,k2:k1:-1), ncid,  &
!     &             'const', strt3d, cnt3d)
!
                call Ncwr_3d (concentration(ic)%pArray3D(il,ij,k1:k2), ncid,  &
     &             'const', strt3d, cnt3d)
!
              end if
            end do
!
          end if
!
!     ========
        end do
      end do
!     ========
!
!
!     ==============
      call Ncdo_Sync (ncid)
!     ==============
!
      return
!
      end subroutine Write_Netcdf_Col
!EOC
!-------------------------------------------------------------------------------
end module GmiControlColumn_mod
