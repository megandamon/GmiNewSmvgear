!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlTendencies_mod
!
! !INTERFACE:
!
module GmiControlTendencies_mod
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use  GmiSub2Glob_mod, only : subDomain2Global
      use GmiFileOperations_mod , only : makeOutfileName
      use m_netcdf_io_close , only : Nccl, Nccl_Noerr
      use m_netcdf_io_write , only : Ncwr_2d_Int, Ncwr_5d, Ncwr_3d, Ncwr_2d
      use m_netcdf_io_write , only : Ncwr_1d, Ncwr_Scal, Ncwr_2d_Char
      use m_netcdf_io_create, only : Ncdo_Sync, Nccr_Wr
      use m_netcdf_io_define, only : NcDef_glob_attributes, NcDef_dimension
      use m_netcdf_io_define, only : NcDef_var_attributes, NcDef_variable
      use m_netcdf_io_define, only : NcSetFill, NcEnd_def
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration,        &
     &       Get_net_cum_tend, Set_net_cum_tend
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_rootProc,        &
     &       Get_procID, Get_iAmRootProc, Get_numDomains, &
     &       Get_communicatorWorld, Get_map1_u
      use GmiDomainDecomposition_mod, only : Get_mcorGlob, Get_latdeg, Get_londeg

      use GmiGrid_mod, only : t_gmiGrid, Get_numSpecies, Get_i1, Get_i2,       &
     &       Get_ju1, Get_j2, Get_k1, Get_k2, Get_i1_gl, Get_i2_gl, Get_ju1_gl,&
     &       Get_j2_gl, Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jhi_gl
      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Get_pr_diag,         &
     &       Get_num_tend_outrecs, Get_tend_outrec_map,                        &
     &       Get_num_const_outrecs, Get_const_outrec_map,                      &
     &       Get_problem_name, Get_hdr_var_name, Get_hdf_dim_name,             &
     &       Get_rec_dim_name, Get_lat_dim_name, Get_lon_dim_name, Get_k1r_gl, &
     &       Get_prs_dim_name, Get_spc_dim_name, Get_tendOutputFrequency, Get_k2r_gl

      use GmiChemistryMethod_mod, only : t_Chemistry, Get_const_labels, Get_mw

      use GmiMetFieldsControl_mod, only : t_metFields, Get_pt, Get_mdt, &
     &       Get_met_opt, Get_metdata_name, Get_ai, Get_bi, Get_am, Get_bm

      use GmiTimeControl_mod, only : t_GmiClock, GmiSplitDateTime, isOutTime, &
     &       Get_curGmiDate, Get_curGmiTime, Get_gmiSeconds, Get_gmiTimeStep, &
     &       Get_numTimeSteps
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: controlOutputTendencies, initializeOutputTendencies
      public  :: finalizeOutputTendencies

#     include "gmi_diag_constants_llnl.h"
#     include "gmi_phys_constants.h"
#     include "GmiParameters.h"
!
! !DESCRIPTION:
! Contains the necessary routines for the tendencies outputs.
! Three routines are visible from the outside:
! %
! \begin{description}
! \item[initializeOutputTendencies:] called once in the initialization steps to
!      open the netCDF file, define variables in the file, write out header
!      information and allocate variables.
! \item[controlOutputTendencies:] called at the end of each model time step to
!      update variables and if needed for communications and writing out
!      data in the file.
! \item[finalizeOutputTendencies:] called at the end of the integration to
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

      ! netCDF file identifier
      integer,          save :: ncid_tnd

      ! Counter for the number of records
      integer,          save :: rnum_out_tnd

      ! Grid information of the variables to be written out
      integer,          save :: i1, i2     ! longitude
      integer,          save :: i1_gl, i2_gl
      integer,          save :: ilo_gl, ihi_gl
      integer,          save :: ju1        ! latitude
      integer,          save :: j2
      integer,          save :: ju1_gl, j2_gl
      integer,          save :: julo_gl, jhi_gl
      integer,          save :: kx1, kx2   ! vertical
      integer,          save :: k1, k2
      integer,          save :: numLon     ! number of longitudes
      integer,          save :: numLat     ! number of latitudes
      integer,          save :: numVert    ! number of vertical levels
      integer,          save :: numSpecies ! number of species
      integer,          save :: num_const_outrecs ! # species to output
      integer,          save :: num_tend_outrecs  ! # species selected for output

      logical,          save :: iAmRootProc
      integer,          save :: procID, rootProc
      integer,          save :: numDomains
      integer,          save :: commuWorld
      integer, pointer, save :: map1_u       (:,:,:) => null()

      integer,          save :: nhdf
      logical,          save :: pr_diag
      real*8 ,            save :: tendOutputFrequency
      integer,            save :: const_outrec_map(MAX_NUM_CONST_GIO)
      integer,            save :: tend_outrec_map(MAX_NUM_CONST_GIO)
      character (len=MAX_LENGTH_VAR_NAME), save :: rec_dim_name, prs_dim_name, spc_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: hdr_var_name
      character (len=MAX_LENGTH_VAR_NAME), save :: hdf_dim_name, lon_dim_name, lat_dim_name
      real*8 ,            save :: mdt
      integer           , save :: met_opt
      real*8 , pointer, save :: mw       (:) => null()
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initializeOutputTendencies
!
! !INTERFACE:
!
      subroutine initializeOutputTendencies(gmiGrid, gmiDomain, Diagnostics, &
     &                     Chemistry, metFields)
!
! !USES:
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_gmiGrid    ), intent(in) :: gmiGrid
      type(t_gmiDomain  ), intent(in) :: gmiDomain
      type(t_metFields  ), intent(in) :: metFields
      type(t_Chemistry  ), intent(in) :: Chemistry  
      type(t_Diagnostics), intent(in) :: Diagnostics
!
! !DESCRIPTION:
! The routines (1) opens a netCDF file, (2) defines variables in the file,
! (3) writes out header data in the file, and (4) allocates varaibles.
!
! !LOCAL VARIABLES:
      character (len=80), save :: outmain_name
      integer :: in1, k1r_gl, k2r_gl
      integer, allocatable :: ewflag_org(:)
      integer, allocatable :: nspole_org(:)
      character (len=MAX_LENGTH_FILE_NAME) :: fname
      character (len=MAX_LENGTH_FILE_NAME) :: problem_name
      character (len=MAX_LENGTH_SPECIES_NAME), pointer :: const_labels(:)
!EOP
!-----------------------------------------------------------------------------
!BOC

      !#####################################################
      ! Initialization of variables used in all the routines
      !#####################################################

      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write (6,*) 'initializeOutputTendencies called by ', procID

      call Get_k1r_gl(Diagnostics, k1r_gl)
      call Get_k2r_gl(Diagnostics, k2r_gl)
      call Get_hdr_var_name(Diagnostics, hdr_var_name)
      call Get_spc_dim_name(Diagnostics, spc_dim_name)
      call Get_rec_dim_name(Diagnostics, rec_dim_name)
      call Get_prs_dim_name(Diagnostics, prs_dim_name)
      call Get_hdf_dim_name(Diagnostics, hdf_dim_name)
      call Get_lat_dim_name(Diagnostics, lat_dim_name)
      call Get_lon_dim_name(Diagnostics, lon_dim_name)
      call Get_tendOutputFrequency(Diagnostics, tendOutputFrequency)
      call Get_tend_outrec_map(Diagnostics, tend_outrec_map)
      call Get_num_tend_outrecs(Diagnostics, num_tend_outrecs)
      call Get_const_outrec_map(Diagnostics, const_outrec_map)
      call Get_num_const_outrecs(Diagnostics, num_const_outrecs)

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

      allocate(map1_u(2, 2, numDomains))
      call Get_map1_u (gmiDomain, map1_u)

      kx1     = k1r_gl
      kx2     = k2r_gl

      numLon  = i2_gl - i1_gl + 1
      numLat  = j2_gl - ju1_gl + 1
      numVert = kx2 - kx1 + 1

      nhdf = NETCDF_HDF

      call Get_mdt(metFields, mdt)
      call Get_met_opt(metFields, met_opt)

      allocate(mw(numSpecies))
      call Get_mw(Chemistry, mw)

      allocate(const_labels(numSpecies))

      if (iAmRootProc) then

         !###########################
         ! Initialize the output file
         !###########################

         ! Determine the file name
         call Get_problem_name(Diagnostics, problem_name)

         call makeOutfileName (fname, '.tend.nc', problem_name)

         ! Create the file and assign a file identifier
         call Nccr_Wr (ncid_tnd, fname)

         ! Define the variables in the file
         call Define_Netcdf_Out_Tend()

         ! Write header data
         call Get_const_labels(Chemistry, const_labels)

         call Write_Netcdf_Hdr_Tend (gmiDomain, const_labels, metFields)

         call Ncdo_Sync (ncid_tnd)

         ! Initialize the counter for the number of records
         rnum_out_tnd = 1
      end if

      !########################
      ! Allocation of variables
      !########################

      return

      end subroutine initializeOutputTendencies
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: controlOutTendencies
!
! !INTERFACE:
!
      subroutine controlOutputTendencies (last_tstp, &
     &                  SpeciesConcentration, gmiClock)
!
! !USES:
      use GmiTimeControl_mod, only : GmiSplitDateTime

      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: last_tstp
      type(t_gmiClock), intent(in) :: gmiClock
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !INPUT/OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! Controls the tendencies output file. It is called at each time step.
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
      if (pr_diag) Write (6,*) 'controlOutputTendencies called by ', procID

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
     &      month_save, month, day, nhms, ndt, gmi_sec, tendOutputFrequency)

      month_save = month

      if (last_tstp .and. (.not. time_for_nc)) then
        if (Mod (num_time_steps, Nint (mdt) / ndt) == 0) then

          time_for_nc = .true.

        end if
      end if

      if (time_for_nc) then
         !====================
         call Prep_Netcdf_Tend (time_for_nc, SpeciesConcentration)
         !====================
      endif

!     ================
      if (time_for_nc) then
!     ================

         if (iAmRootProc) call Write_Netcdf_Tend (nymd, nhms, gmi_sec)

!       ======================
        call bufferOutput_Tend (SpeciesConcentration)
!       ======================

        if (iAmRootProc) rnum_out_tnd = rnum_out_tnd + 1

!     ======
      end if
!     ======

      return

      end subroutine controlOutputTendencies
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINES: finalizeOutputTendencies
!
! !INTERFACE:
!
      subroutine finalizeOutputTendencies()
!
      implicit none
!
! !DESCRIPTION:
! Deallocates variables necessary to produce tendencies outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'finalizeOutputTendencies called by ', procID

      if (iAmRootProc) then
         call Nccl_Noerr (ncid_tnd)
      end if

      end subroutine finalizeOutputTendencies
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Define_Netcdf_Out_Tend
!
! !INTERFACE:
!

      subroutine Define_Netcdf_Out_Tend ()
!
! !USES:
      use m_ncGeneralOpsOutput, only : Define_Netcdf_Out_Gen

      implicit none

#     include "netcdf.inc"
!
! !DESCRIPTION:
!   Makes the necessary definitions for the tendency
!   diagnostic NetCDF output file.
!
! !DECLARED PARAMETERS:
      character (len=MAX_LENGTH_VAR_NAME), parameter :: OPR_DNAM = 'operator_dim'
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_VAR_NAME) :: prsp1_dim_name

      integer :: ierr
      integer :: nchr1, nchr2
      integer :: nstr
      integer :: omode
      integer :: pos1
      integer :: varid

      integer :: chrd1(1), chrd2(1)
      integer :: hdfd (1)
      integer :: lond (1), latd (1)
      integer :: oprd (1)
      integer :: prsd (1), prsp1d(1)
      integer :: spcd (1), recd (1), spcd_2(1)
      integer :: strd (1)

      integer :: var1 (1)
      integer :: var2 (2)
      integer :: var2c(2), var2l(2), var2m(2)
      integer :: var5 (5)
      integer :: var6 (6)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Define_Netcdf_Out_Tend called by ', procID

      nchr1 = MAX_LENGTH_SPECIES_NAME
      nchr2 = 50

      nstr =   1

!     ==========================
      call Define_Netcdf_Out_Gen  &
!     ==========================
     &  (ncid_tnd, nhdf, numLon, numLat, numVert, hdfd, lond, latd, prsd, &
     &   hdf_dim_name, lon_dim_name, lat_dim_name, prs_dim_name)

!     ------------------
!     Define dimensions.
!     ------------------

      prsp1_dim_name = prs_dim_name
      pos1 = Len_Trim (prs_dim_name) + 1
      prsp1_dim_name(pos1:pos1+2) = 'p1'

      call NcDef_dimension(ncid_tnd, prsp1_dim_name, numVert+1, prsp1d(1))
!                                     --------------
      call NcDef_dimension(ncid_tnd, 'chr_dim1', nchr1, chrd1(1))
!                                     ----------
      if (met_opt /= 1) then
        call NcDef_dimension(ncid_tnd, 'chr_dim2', nchr2, chrd2(1))
!                                       ----------
        call NcDef_dimension(ncid_tnd, 'str_dim', nstr, strd(1))
!                                       ---------
      end if

      call NcDef_dimension(ncid_tnd, spc_dim_name, num_const_outrecs, spcd(1))
!                                     ------------
      call NcDef_dimension (ncid_tnd,'tend_spc_dim', num_tend_outrecs,spcd_2(1))
!                                    -------------

      call NcDef_dimension(ncid_tnd, OPR_DNAM, NUM_OPERATORS, oprd(1))
!                                     --------
      call NcDef_dimension(ncid_tnd, rec_dim_name, NF_UNLIMITED, recd(1))
!                                     ------------
!     -----------------------------------------
!     Define variables and variable attributes.
!     -----------------------------------------

      call NcDef_variable  &
     &         (ncid_tnd, prsp1_dim_name, NF_FLOAT, 1, prsp1d, varid)
!                         --------------
      call NcDef_var_attributes (ncid_tnd, varid, 'long_name', 'Pressure + 1')
      call NcDef_var_attributes (ncid_tnd, varid, 'units', 'mb')

      call NcDef_variable  &
     &         (ncid_tnd, spc_dim_name, NF_FLOAT, 1, spcd, varid)
!                         ------------
      call NcDef_var_attributes (ncid_tnd, varid, 'long_name', 'Species index')
      call NcDef_var_attributes (ncid_tnd, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_tnd, varid, 'coord_labels', 'const_labels')
      call NcDef_var_attributes (ncid_tnd, varid, 'selection_category', 'NULL')

      call NcDef_variable  &
     &         (ncid_tnd, 'operator_dim', NF_FLOAT, 1, oprd, varid)
!                         --------------
      call NcDef_var_attributes (ncid_tnd, varid, 'long_name', 'Operator index')
      call NcDef_var_attributes (ncid_tnd, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_tnd, varid, 'coord_labels', 'operator_labels')
      call NcDef_var_attributes (ncid_tnd, varid, 'selection_category', 'NULL')

      if (met_opt /= 1) then
        var2m(:) = (/ chrd2(1), strd(1) /)
        call NcDef_variable  &
     &           (ncid_tnd, 'metdata_name', NF_CHAR, 2, var2m, varid)
!                           --------------
        call NcDef_var_attributes (ncid_tnd, varid, 'long_name', 'Metdata name')
        call NcDef_var_attributes (ncid_tnd, varid, 'units', 'unitless')
      end if

      var2c(:) = (/ chrd1(1), spcd(1) /)
      call NcDef_variable  &
     &         (ncid_tnd, 'const_labels', NF_CHAR, 2, var2c, varid)
!                         --------------
      call NcDef_var_attributes (ncid_tnd, varid, 'long_name', 'Species name')
      call NcDef_var_attributes (ncid_tnd, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_tnd, varid, 'selection_category', 'NULL')

      call NcDef_variable (ncid_tnd, 'tend_spc_dim', NF_FLOAT, 1, spcd_2, varid)
!                                    --------------
      call NcDef_var_attributes (ncid_tnd, varid, 'long_name', 'Species index for tendencies')
      call NcDef_var_attributes (ncid_tnd, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_tnd, varid,'coord_labels', 'tend_spc_labels')
      call NcDef_var_attributes (ncid_tnd, varid, 'selection_category', 'NULL')

      var2c(:) = (/ chrd1(1), spcd_2(1) /)
      call NcDef_variable (ncid_tnd, 'tend_spc_labels', NF_CHAR, 2, var2c, varid)
!                                    -----------------
      call NcDef_var_attributes (ncid_tnd,varid,'long_name','Species names for tendencies')
      call NcDef_var_attributes (ncid_tnd, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_tnd, varid, 'selection_category', 'NULL')

      var2c(:) = (/ chrd1(1), oprd(1) /)
      call NcDef_variable (ncid_tnd, 'operator_labels', NF_CHAR, 2, var2c, varid)
!                                    -----------------
      call NcDef_var_attributes (ncid_tnd, varid, 'long_name', 'Operator name')
      call NcDef_var_attributes (ncid_tnd, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_tnd, varid, 'selection_category', 'NULL')

      var1 (:) = (/ spcd(1) /)
      call NcDef_variable (ncid_tnd, 'mw', NF_FLOAT, 1, var1, varid)
!                                     ----
      call NcDef_var_attributes (ncid_tnd, varid, 'long_name', 'Molecular weight')
      call NcDef_var_attributes (ncid_tnd, varid, 'units', 'g/mol')

      var2l(:) = (/ lond(1), latd(1) /)
      call NcDef_variable (ncid_tnd, 'mcor', NF_FLOAT, 2, var2l, varid)
!                                     ------
      call NcDef_var_attributes (ncid_tnd, varid, 'long_name', 'Area of grid box')
      call NcDef_var_attributes (ncid_tnd, varid, 'units', 'm^2')

      var2 (:) = (/ hdfd(1), recd(1) /)
      call NcDef_variable  &
     &         (ncid_tnd, hdr_var_name, NF_INT, 2, var2, varid)
!                         ------------
      call NcDef_var_attributes (ncid_tnd, varid, 'long_name', 'Header')
      call NcDef_var_attributes (ncid_tnd, varid, 'units', 'gmi_sec, nymd, nhms')

      var6(:) = (/ lond(1), latd(1), prsd(1), spcd_2(1), oprd(1), recd(1) /)

      call NcDef_variable (ncid_tnd, 'net_cum_tnd', NF_FLOAT, 6, var6, varid)
!                                    -------------
      call NcDef_var_attributes (ncid_tnd, varid, 'long_name', 'Net cumulative tendencies')
      call NcDef_var_attributes (ncid_tnd, varid, 'units', 'kg')

!     -------------------------
!     Define global attributes.
!     -------------------------

      call NcDef_glob_attributes (ncid_tnd, 'title',  &
     &          'Gmimod net cumulative tendencies diagnostic file')

      call NcSetFill (ncid_tnd, NF_NOFILL, omode)

      call NcEnd_def (ncid_tnd)

      return

      end subroutine Define_Netcdf_Out_Tend
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Hdr_Tend
!
! !INTEFACE:
!
      subroutine Write_Netcdf_Hdr_Tend  (gmiDomain, const_labels, metFields)
!
! !USES:
      use m_ncGeneralOpsOutput, only : Write_Netcdf_Hdr_Gen
!
      implicit none
!
! !INPUT PARAMETERS:
      character (len=MAX_LENGTH_SPECIES_NAME), intent(in) :: const_labels(:)
      type(t_gmiDomain), intent(in) :: gmiDomain
      type(t_metFields), intent(in) :: metFields
!
! !DESCRIPTION:
!  Creates some header information for the net cumulative
!   tendencies NetCDF output file and writes it out.
!
!LOCAL VARIABLES:
      character (len=50) :: metdata_name
      character (len=MAX_LENGTH_VAR_NAME) :: prsp1_dim_name

      character (len=MAX_LENGTH_SPECIES_NAME) :: spclab  (num_const_outrecs)
      character (len=MAX_LENGTH_SPECIES_NAME) :: spclab_2(num_tend_outrecs )
      character (len=16) :: oprlab(NUM_OPERATORS)
      character (len=50) :: metnam(1)

      integer :: ic
      integer :: pos1

      integer :: cnt1d (1), cnt2d (2)
      integer :: strt1d(1), strt2d(2)

      real*8  :: prsdat  (1:numVert)
      real*8  :: prsdatp1(1:numVert+1)

      real*8, allocatable :: prsdat_r  (:)
      real*8, allocatable :: prsdatp1_r(:)

      real*8  :: spcdat  (num_const_outrecs)
      real*8  :: spcdat_2(num_tend_outrecs )
      real*8  :: oprdat(NUM_OPERATORS)
      real*8 , allocatable :: mcor(:,:)
      real*8 , allocatable :: londeg(:), latdeg(:)
      real*8 , allocatable :: am(:), bm(:)
      real*8 , allocatable :: ai(:), bi(:)
      real*8 :: pt
      real*8 :: loc_mw(1)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Hdr_Tend called by ', procID

      allocate(latdeg(ju1_gl:j2_gl))
      call Get_latdeg(gmiDomain, latdeg)

      allocate(londeg(i1_gl:i2_gl))
      call Get_londeg(gmiDomain, londeg)

      allocate(mcor(i1_gl:i2_gl,ju1_gl:j2_gl))
      call Get_mcorGlob(gmiDomain, mcor)

!     =========================
      call Write_Netcdf_Hdr_Gen  &
!     =========================
     &  (ncid_tnd, latdeg, londeg, pr_diag, procID, i1_gl, i2_gl, ju1_gl, j2_gl, &
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

      call Ncwr_1d (prsdat(kx1:kx2), ncid_tnd, prs_dim_name, strt1d, cnt1d)

      cnt1d(1) = numvert + 1

      prsp1_dim_name = prs_dim_name
      pos1 = Len_Trim (prs_dim_name) + 1
      prsp1_dim_name(pos1:pos1+2) = 'p1'

      prsdatp1(1:numVert+1) = (ai(kx1-1:kx2) * pt)  + (bi(kx1-1:kx2) * 1000.0d0)

      call Ncwr_1d (prsdatp1, ncid_tnd, prsp1_dim_name, strt1d, cnt1d)

!     -------------
!     Metdata name.
!     -------------

      if (met_opt /= 1) then
         call Get_metdata_name (metFields, metdata_name)

        metnam(1) = metdata_name

        strt2d(:) = (/  1, 1 /)
        cnt2d (:) = (/ 50, 1 /)

        call Ncwr_2d_Char (metnam, ncid_tnd, 'metdata_name', strt2d, cnt2d)
      end if

!     ------------
!     Species map.
!     ------------

      if (pr_diag) then
        Write (6,*) 'Tendency NetCDF output info:'
        Write (6,*) '  num_const_outrecs:         ', num_const_outrecs

        do ic = 1, num_const_outrecs
          Write (6,900) ic, const_outrec_map(ic)
        end do

 900    format ('   const_outrec to species mapping:  ', i4, ' => ', i4)
      end if

      do ic = 1, num_const_outrecs
        spcdat(ic) = const_outrec_map(ic)
      end do

      strt1d(1) = 1
      cnt1d (1) = num_const_outrecs

      call Ncwr_1d (spcdat, ncid_tnd, spc_dim_name, strt1d, cnt1d)

!     ---------------
!     Species labels.
!     ---------------

      do ic = 1, num_const_outrecs
        spclab(ic) = const_labels(const_outrec_map(ic))
      end do

      strt2d(:) = (/  1, 1 /)
      cnt2d (:) = (/ MAX_LENGTH_SPECIES_NAME, num_const_outrecs /)

      call Ncwr_2d_Char (spclab, ncid_tnd, 'const_labels', strt2d, cnt2d)

!     -----------------
!     Operator indices.
!     -----------------

      if (pr_diag) then
        Write (6,*) 'Tendency NetCDF output info:'
        Write (6,*) '  NUM_OPERATORS:         ', NUM_OPERATORS

        do ic = 1, NUM_OPERATORS
          Write (6,910) ic, OPERATOR_NAME(ic)
        end do

 910    format ('   operators:  ', i4, ' => ', a16)
      end if

      do ic = 1, NUM_OPERATORS
        oprdat(ic) = ic
      end do

      strt1d(1) = 1
      cnt1d (1) = NUM_OPERATORS

      call Ncwr_1d (oprdat, ncid_tnd, 'operator_dim', strt1d, cnt1d)

!     ----------------
!     Operator labels.
!     ----------------

      do ic = 1, NUM_OPERATORS
        oprlab(ic) = OPERATOR_NAME(ic)
      end do

      strt2d(:) = (/  1, 1 /)
      cnt2d (:) = (/ 16, NUM_OPERATORS /)

      call Ncwr_2d_Char (oprlab, ncid_tnd, 'operator_labels', strt2d, cnt2d)

!     -----------------
!     Molecular weight.
!     -----------------

      cnt1d (1) = 1

      do ic = 1, num_const_outrecs

        strt1d(1) = ic
        loc_mw(1) = mw(const_outrec_map(ic))
        call Ncwr_1d (loc_mw, ncid_tnd, 'mw', strt1d, cnt1d)

      end do

!     --------------
!     Grid box area.
!     --------------

      strt2d(:) = (/ 1, 1 /)
      cnt2d (:) = (/ numLon, numLat /)

      call Ncwr_2d (mcor, ncid_tnd, 'mcor', strt2d, cnt2d)

!     ---------------------------------------------
!     Cum. Tendency Species map and Species labels.
!     ---------------------------------------------
!     This portion of the code is only used if the user has
!     selected in the namelist file the species for which
!     he/she wants tendency diagnostic to be saved.
!

      do ic = 1, num_tend_outrecs
         spcdat_2(ic) = tend_outrec_map(ic)
      end do

      strt1d(1) = 1
      cnt1d (1) = num_tend_outrecs

      call Ncwr_1d (spcdat_2,ncid_tnd,'tend_spc_dim',strt1d,cnt1d)

      do ic = 1, num_tend_outrecs
         spclab_2(ic) = const_labels(tend_outrec_map(ic))
      end do

      strt2d(:) = (/  1, 1 /)
      cnt2d (:) = (/ MAX_LENGTH_SPECIES_NAME, num_tend_outrecs /)

      call Ncwr_2d_Char (spclab_2, ncid_tnd, 'tend_spc_labels', strt2d, cnt2d)

      return

      end subroutine Write_Netcdf_Hdr_Tend
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Prep_Netcdf_Tend
!
! !INTERFACE:
!
      subroutine Prep_Netcdf_Tend (time_for_nc, SpeciesConcentration)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: time_for_nc ! time for NetCDF output?
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !DESCRIPTION:
! Prepares the tend (net cumulative tendencies) NetCDF output.
!
! !LOCAL VARIABLES:
      integer :: ic, io, icx
      real*8  :: mw_fac
      type (t_GmiArrayBundle), pointer :: net_cum_tend(:,:)
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Prep_Netcdf_Tend called by ', procID

      if (time_for_nc) then

         call Get_net_cum_tend(SpeciesConcentration, net_cum_tend)

!       ------------------------
!       Scale values for output.
!       ------------------------

         do icx = 1, NUM_OPERATORS
            do ic = 1, num_tend_outrecs
               mw_fac = mw(tend_outrec_map(ic)) / MWTAIR
               net_cum_tend(ic,icx)%pArray3D(:,:,:) = net_cum_tend(ic,icx)%pArray3D(:,:,:) * mw_fac
            end do
         end do

         call Set_net_cum_tend(SpeciesConcentration, net_cum_tend)

      end if

      return

      end subroutine Prep_Netcdf_Tend
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Tend
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Tend (nymd, nhms, gmi_sec)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: nymd, nhms
      real*8 , intent(in) :: gmi_sec
!
! !DESCRIPTION:
! Writes the tend (net cumulative tendencies) NetCDF output.
!
! !LOCAL VARIABLES:
      integer :: cnt2d (2)
      integer :: strt2d(2)
      integer :: hdr(NETCDF_HDF)
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Tend called by ', procID

      strt2d(:) = (/ 1, rnum_out_tnd /)

      cnt2d (:) = (/ NETCDF_HDF, 1 /)

      hdr(1) = Nint (gmi_sec)
      hdr(2) = nymd
      hdr(3) = nhms

      call Ncwr_2d_Int (hdr, ncid_tnd, hdr_var_name, strt2d, cnt2d)

!     ==============
      call Ncdo_Sync (ncid_tnd)
!     ==============

      return

      end subroutine Write_Netcdf_Tend
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutput_Tend
!
! !INTERFACE:
!
      subroutine bufferOutput_Tend (SpeciesConcentration)
!
! !USES:
      use GmiBuffer_mod, only : bufferOutput5d_Tend

      implicit none
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !DESCRIPTION:
!  Buffers the tendencies data to the output file. The worker PEs send a
!  slice of the data to the root PE, the root PE writes it out, they then send
!  another slice and it is written out, and so on.
!
! !DEFINED PARAMETERS:
      character (len=MAX_LENGTH_VAR_NAME), parameter :: TEND_VNAM = 'net_cum_tnd'
      integer           , parameter :: SG_NCUMT2 = 6100
      integer                       :: ic, icx
      type (t_GmiArrayBundle), pointer :: net_cum_tend(:,:)
      real*8, allocatable :: net_cum_tend_nc(:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'bufferOutput_Tend called by ', procID

      if ((NUM_OPERATORS > 0) .and. (num_tend_outrecs > 0)) then

         Allocate (net_cum_tend_nc(i1:i2, ju1:j2, kx1:kx2))
         net_cum_tend_nc = 0.0d0
         
         call Get_net_cum_tend(SpeciesConcentration, net_cum_tend)

!        =========================
         call bufferOutput5d_Tend  &
!        =========================
     &     (TEND_VNAM, kx1, kx2, commuWorld, SG_NCUMT2,  &
     &     ncid_tnd, num_tend_outrecs, NUM_OPERATORS, rnum_out_tnd,  &
     &     net_cum_tend, net_cum_tend_nc, &
     &    map1_u, numDomains, rootProc, procID, &
     &   i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)

         do ic = 1, NUM_OPERATORS
            do icx = 1, num_tend_outrecs
               net_cum_tend(icx,ic)%pArray3D(:,:,:) = 0.0d0
            end do
         end do

         call Set_net_cum_tend(SpeciesConcentration, net_cum_tend)

      end if

      return

      end subroutine bufferOutput_Tend
!EOC
!-----------------------------------------------------------------------------
end module GmiControlTendencies_mod
