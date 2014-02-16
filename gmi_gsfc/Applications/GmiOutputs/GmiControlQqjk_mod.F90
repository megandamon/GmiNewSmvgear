!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlQqjk_mod
!
! !INTERFACE:
!
module GmiControlQqjk_mod
!
! !USES:
      use  GmiSub2Glob_mod      , only : subDomain2Global
      use  GmiMessagePassing_mod, only : synchronizeGroup
      use GmiFileOperations_mod , only : makeOutfileName
      use m_netcdf_io_close     , only : Nccl, Nccl_Noerr
      use m_netcdf_io_write     , only : Ncwr_2d_Int, Ncwr_5d, Ncwr_3d, Ncwr_2d
      use m_netcdf_io_write     , only : Ncwr_1d, Ncwr_Scal, Ncwr_2d_Char
      use m_netcdf_io_create    , only : Ncdo_Sync, Nccr_Wr
      use m_netcdf_io_define    , only : NcDef_glob_attributes, NcDef_dimension
      use m_netcdf_io_define    , only : NcDef_var_attributes, NcDef_variable
      use m_netcdf_io_define    , only : NcSetFill, NcEnd_def
      use GmiChemistryMethod_mod, only : t_Chemistry, Get_qjgmi, Get_num_qjs,  &
     &       Get_num_qjo, Get_num_qks, Get_num_active, Get_qj_labels,          &
     &       Get_qk_labels
      use GmiChemistryMethod_mod, only : Get_yda, Get_qqjda, Get_qqkda
      use GmiChemistryMethod_mod, only : Get_qqjgmi, Get_qqkgmi
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      use GmiDomainDecomposition_mod, only : t_gmiDomain
      use GmiDomainDecomposition_mod, only : Get_rootProc, Get_procID, Get_iAmRootProc
      use GmiDomainDecomposition_mod, only : Get_numDomains
      use GmiDomainDecomposition_mod, only : Get_communicatorWorld
      use GmiDomainDecomposition_mod, only : Get_map1_u
      use GmiDomainDecomposition_mod, only : Get_latdeg, Get_londeg

      use GmiGrid_mod, only : t_gmiGrid, Get_i1, Get_i2, Get_ju1, Get_j2,      &
     &       Get_k1, Get_k2, Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl,      &
     &       Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jhi_gl

      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Set_do_qqjk_reset,   &
     &       Get_qqjkOutputFrequency, Get_hdr_var_name, Get_rec_dim_name,             &
     &       Get_prs_dim_name, Get_do_qqjk_inchem, Get_qqj_var_name,           &
     &       Get_qqj_dim_name, Get_qqk_var_name, Get_qqk_dim_name,             &
     &       Get_hdf_dim_name, Get_lon_dim_name, Get_lat_dim_name, Get_do_mean,&
     &       Get_k1r_gl, Get_k2r_gl, Get_pr_qqjk, Get_pr_diag, Get_problem_name

      use GmiTimeControl_mod, only : t_GmiClock, GmiSplitDateTime, isOutTime, &
     &       Get_curGmiDate, Get_curGmiTime, Get_gmiSeconds, Get_gmiTimeStep, &
     &       Get_numTimeSteps

      use GmiMetFieldsControl_mod, only : t_metFields, Get_pt,  &
     &       Get_metdata_name, Get_mdt, Get_ai, Get_bi, Get_am, Get_bm
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: controlOutputQqjk, initializeOutputQqjk
      public  :: finalizeOutputQqjk

#     include "gmi_phys_constants.h"
#     include "GmiParameters.h"
!
! !DESCRIPTION:
! Contains the necessary routines for the qqjk outputs. 
! Three routines are visible from the outside:
! %
! \begin{description}
! \item[initializeOutputQqjk:] called once in the initialization steps to
!      open the netCDF file, define variables in the file, write out header
!      information and allocate variables.
! \item[controlOutputQqjk:] called at the end of each model time step to
!      update variables and if needed for communications and writing out
!      data in the file.
! \item[finalizeOutputQqjk:] called at the end of the integration to
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
      real*8 , pointer, save :: qqj_nc   (:,:,:) => null() ! rates of photolytic processes (molecules/cm^3*s)
      real*8 , pointer, save :: qqk_nc   (:,:,:) => null() ! rates of thermal    processes (molecules/cm^3*s)
      real*8 , pointer, save :: yda_nc   (:,:,:) => null() ! corresponding species concentration to qqj_nc/qqk_nc
!     
      ! Variables manipulate by worker processors only         
      real*8 , pointer, save :: qqj_mean (:,:,:,:) => null()
      real*8 , pointer, save :: qqk_mean (:,:,:,:) => null()
!
      ! netCDF file identifier
      integer,          save :: ncid_qqjk

      ! Counter for the number of records
      integer,          save :: rnum_out_qqjk

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

      integer,          save :: numQjo
      integer,          save :: numQjs
      integer,          save :: numQks
      integer,          save :: numActiveChem
      
      logical,          save :: iAmRootProc
      logical,          save :: pr_diag
      integer,          save :: procID, rootProc
      integer,          save :: numDomains
      integer,          save :: commuWorld
      integer, pointer, save :: map1_u       (:,:,:) => null()

      integer,          save :: nhdf
      character (len=MAX_LENGTH_VAR_NAME), save :: qqj_dim_name, qqj_var_name
      character (len=MAX_LENGTH_VAR_NAME), save :: qqk_dim_name, qqk_var_name
      character (len=MAX_LENGTH_VAR_NAME), save :: rec_dim_name, prs_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: hdf_dim_name, lon_dim_name, lat_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: hdr_var_name
      logical           , save :: pr_qqjk, do_mean, do_qqjk_inchem
      real*8            , save :: qqjkOutputFrequency
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
! !IROUTINE: initializeOutputQqjk
!
! !INTERFACE:
!
      subroutine initializeOutputQqjk(gmiGrid, gmiDomain, Diagnostics, &
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
      character (len=75)  :: err_msg
      integer :: in1, k1r_gl, k2r_gl
      character (len=MAX_LENGTH_FILE_NAME) :: fname
      character (len=MAX_LENGTH_FILE_NAME) :: problem_name
      character (len=MAX_LENGTH_LABELS) :: qj_labels(MAX_NUM_QJ)
      character (len=MAX_LENGTH_LABELS) :: qk_labels(MAX_NUM_QK)
!EOP
!-----------------------------------------------------------------------------
!BOC

      !#####################################################
      ! Initialization of variables used in all the routines
      !#####################################################

      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write (6,*) 'initializeOutputQqjk called by ', procID

      call Get_num_qjo(Chemistry, numQjo)
      call Get_num_qjs(Chemistry, numQjs)
      call Get_num_qks(Chemistry, numQks)
      call Get_num_active(Chemistry, numActiveChem)

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

      call Get_rootProc         (gmiDomain, rootProc  )
      call Get_numDomains       (gmiDomain, numDomains )
      call Get_iAmRootProc      (gmiDomain, iAmRootProc )
      call Get_communicatorWorld(gmiDomain, commuWorld)

      allocate(map1_u(2, 2, numDomains))
      call Get_map1_u (gmiDomain, map1_u)

      call Get_k1r_gl      (Diagnostics, k1r_gl)
      call Get_k2r_gl      (Diagnostics, k2r_gl)
      
      kx1     = k1r_gl
      kx2     = k2r_gl
      
      numLon  = i2_gl - i1_gl + 1
      numLat  = j2_gl - ju1_gl + 1
      numVert = kx2 - kx1 + 1
      
      nhdf = NETCDF_HDF
      
      call Get_do_mean     (Diagnostics, do_mean)
      call Get_qqjkOutputFrequency(Diagnostics, qqjkOutputFrequency)
      call Get_pr_qqjk     (Diagnostics, pr_qqjk)
      call Get_hdr_var_name(Diagnostics, hdr_var_name)
      call Get_rec_dim_name(Diagnostics, rec_dim_name)
      call Get_prs_dim_name(Diagnostics, prs_dim_name)
      call Get_qqj_var_name(Diagnostics, qqj_var_name)
      call Get_qqj_dim_name(Diagnostics, qqj_dim_name)
      call Get_qqk_var_name(Diagnostics, qqk_var_name)
      call Get_qqk_dim_name(Diagnostics, qqk_dim_name)
      call Get_hdf_dim_name (Diagnostics, hdf_dim_name)
      call Get_lat_dim_name (Diagnostics, lat_dim_name)
      call Get_lon_dim_name (Diagnostics, lon_dim_name)
      call Get_do_qqjk_inchem(Diagnostics, do_qqjk_inchem)

      call Get_mdt (metFields, mdt)

      if (iAmRootProc) then

         !###########################
         ! Initialize the output file
         !###########################

         ! Determine the file name
         call Get_problem_name(Diagnostics, problem_name)
         call makeOutfileName (fname, '.qqjk.nc', problem_name)

         ! Create the file and assign a file identifier
         call Nccr_Wr (ncid_qqjk, fname)

         ! Define the variables in the file
         call Define_Netcdf_Out_Qqjk( )

         ! Write header data
         call Get_qj_labels(Chemistry, qj_labels)

         call Get_qk_labels(Chemistry, qk_labels)

         call Write_Netcdf_Hdr_Qqjk (gmiDomain, metFields, qj_labels, qk_Labels)

         call Ncdo_Sync (ncid_qqjk)

         ! Initialize the counter for the number of records
         rnum_out_qqjk = 1
      end if

      !########################
      ! Allocation of variables
      !########################
      call allocateVariablesQqjk()

      return

      end subroutine initializeOutputQqjk
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: controlOutQqjk
!
! !INTERFACE:
!
      subroutine controlOutputQqjk (last_tstp, Diagnostics, Chemistry, gmiClock)
!
! !USES:
      use GmiTimeControl_mod, only : GmiSplitDateTime

      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: last_tstp
      type(t_gmiClock) , intent(in) :: gmiClock
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_Diagnostics) , intent(inOut) :: Diagnostics
      type(t_Chemistry) , intent(inOut) :: Chemistry
!
! !DESCRIPTION:
! Controls the qj output file. It is called at each time step.
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
      if (pr_diag) Write (6,*) 'controlOutputQqjk called by ', procID

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
     &         month_save, month, day, nhms, ndt, gmi_sec, qqjkOutputFrequency)

      month_save = month

      if (last_tstp .and. (.not. time_for_nc)) then
        if (Mod (num_time_steps, Nint (mdt) / ndt) == 0) then

          time_for_nc = .true.

        end if
      end if

      if (time_for_nc .or. do_mean) then
         !====================
         call Prep_Netcdf_Qqjk  (time_for_nc, Chemistry)
         !====================
      endif

!     ================
      if (time_for_nc) then
!     ================

!        =====================
         call Write_Netcdf_Qqjk (nymd, nhms, gmi_sec)
!        =====================

!        ======================
         call bufferOutput_Qqjk (Chemistry)
!        ======================

         if (iAmRootProc) then
            rnum_out_qqjk = rnum_out_qqjk + 1
         end if

         call Set_do_qqjk_reset(Diagnostics, .true.)

!     ======
      end if
!     ======

      return

      end subroutine controlOutputQqjk
!EOC
!-----------------------------------------------------------------------------
!BOP
!     
! !IROUTINES: finalizeOutputQqjk
!
! !INTERFACE:
!
      subroutine finalizeOutputQqjk()
!
      implicit none
! 
! !DESCRIPTION:
! Deallocates variables necessary to produce qj outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'finalizeOutputQqjk called by ', procID

      deallocate (qqj_nc)
      deallocate (qqk_nc)
      if (do_qqjk_inchem) then
         deallocate (yda_nc)
      end if

      if (do_mean) then
         deallocate (qqj_mean)
         deallocate (qqk_mean)
      end if

      if (iAmRootProc) call Nccl_Noerr (ncid_qqjk)

      return

      end subroutine finalizeOutputQqjk
!EOC
!-----------------------------------------------------------------------------
!BOP
!     
! !IROUTINES: allocateVariablesQqjk
!
! !INTERFACE:
!
      subroutine allocateVariablesQqjk()
!
      implicit none
! 
! !DESCRIPTION:
! Allocates variables necessary to produce qj outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'allocateVariablesQqjk called by ', procID

      Allocate (qqj_nc(i1:i2, ju1:j2, kx1:kx2))
      Allocate (qqk_nc(i1:i2, ju1:j2, kx1:kx2))
      qqj_nc = 0.0d0; qqk_nc = 0.0d0

      if (do_qqjk_inchem) then
         Allocate (yda_nc(i1:i2, ju1:j2, kx1:kx2))
         yda_nc = 0.0d0
      end if

      if (do_mean) then
         Allocate (qqj_mean(i1:i2, ju1:j2, kx1:kx2, numQjo))
         Allocate (qqk_mean(i1:i2, ju1:j2, kx1:kx2, numQks))
         qqj_mean = 0.0d0
         qqk_mean = 0.0d0
      end if

      return

      end subroutine allocateVariablesQqjk
!EOC
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Define_Netcdf_Out_Qqjk
!
! DESCRIPTION
!   This routine makes the necessary definitions for the qj NetCDF output
!   file.

      subroutine Define_Netcdf_Out_Qqjk ()

      use m_ncGeneralOpsOutput, only : Define_Netcdf_Out_Gen

      implicit none

#     include "netcdf.inc"


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ierr
      integer :: omode
      integer :: varid

      integer :: chrd(1)
      integer :: hdfd(1)
      integer :: lond(1), latd(1), prsd(1)
      integer :: qqjsd(1), qqksd(1)
      integer :: ydasd(1), recd (1)

      integer :: var2(2)
      integer :: var5(5)

      character (len=MAX_LENGTH_VAR_NAME), parameter :: YDA_DNAM = 'yda_dim'

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Define_Netcdf_Out_Qqjk called by ', procID
      end if

!     ==========================
      call Define_Netcdf_Out_Gen  &
!     ==========================
     &  (ncid_qqjk, nhdf, numLon, numLat, numVert, hdfd, lond, latd, prsd,  &
     &   hdf_dim_name, lon_dim_name, lat_dim_name, prs_dim_name)

!     ------------------
!     Define dimensions.
!     ------------------

      call NcDef_dimension(ncid_qqjk, qqj_dim_name, numQjs, qqjsd(1))
!                                      ------------
      call NcDef_dimension(ncid_qqjk, qqk_dim_name, numQks, qqksd(1))
!                                      ------------
      call NcDef_dimension (ncid_qqjk, 'chr_dim', 80, chrd(1))
!                                       ---------
      if (do_qqjk_inchem) then
        call NcDef_dimension(ncid_qqjk, YDA_DNAM, numActiveChem, ydasd(1))
!                                        --------
      end if

      call NcDef_dimension(ncid_qqjk, rec_dim_name, NF_UNLIMITED, recd(1))
!                                      ------------

!     -----------------------------------------
!     Define variables and variable attributes.
!     -----------------------------------------

      call NcDef_variable (ncid_qqjk, qqj_dim_name, NF_FLOAT, 1, qqjsd, varid)
!                                     ------------

      call NcDef_var_attributes(ncid_qqjk, varid, 'long_name', 'Qqjs')
      call NcDef_var_attributes(ncid_qqjk, varid, 'units', 'unitless')
      call NcDef_var_attributes(ncid_qqjk, varid, 'coord_labels', 'qqj_labels')
      call NcDef_var_attributes(ncid_qqjk, varid, 'selection_category', 'NULL')

      call NcDef_variable (ncid_qqjk, qqk_dim_name, NF_FLOAT, 1, qqksd, varid)
!                                     ------------
      call NcDef_var_attributes(ncid_qqjk, varid, 'long_name', 'Qqks')
      call NcDef_var_attributes(ncid_qqjk, varid, 'units', 'unitless')
      call NcDef_var_attributes(ncid_qqjk, varid, 'coord_labels', 'qqk_labels')
      call NcDef_var_attributes(ncid_qqjk, varid, 'selection_category', 'NULL')

      if (do_qqjk_inchem) then
        call NcDef_variable (ncid_qqjk, YDA_DNAM, NF_FLOAT, 1, ydasd, varid)
!                                       --------
        call NcDef_var_attributes (ncid_qqjk, varid, 'long_name', 'Ydas')
      end if

      var2(:) = (/ hdfd(1), recd(1) /)
      call NcDef_variable (ncid_qqjk, hdr_var_name, NF_INT, 2, var2, varid)
!                                      ------------

      var2(:) = (/ chrd(1), qqjsd(1) /)
      call NcDef_variable (ncid_qqjk, 'qqj_labels', NF_CHAR, 2, var2, varid)
!                                      ------------
      call NcDef_var_attributes(ncid_qqjk, varid, 'long_name', 'qqj name')
      call NcDef_var_attributes(ncid_qqjk, varid, 'units', 'unitless')
      call NcDef_var_attributes(ncid_qqjk, varid, 'selection_category', 'NULL')

      var2(:) = (/ chrd(1), qqksd(1) /)
      call NcDef_variable (ncid_qqjk, 'qqk_labels', NF_CHAR, 2, var2, varid)
!                                      ------------
      call NcDef_var_attributes(ncid_qqjk, varid, 'long_name', 'qqk name')
      call NcDef_var_attributes(ncid_qqjk, varid, 'units', 'unitless')
      call NcDef_var_attributes(ncid_qqjk, varid, 'selection_category', 'NULL')

      var5(:) = (/ lond(1), latd(1), prsd(1), qqjsd(1), recd(1) /)
      call NcDef_variable (ncid_qqjk, qqj_var_name, NF_FLOAT, 5, var5, varid)
!                                      ------------

      var5(:) = (/ lond(1), latd(1), prsd(1), qqksd(1), recd(1) /)
      call NcDef_variable (ncid_qqjk, qqk_var_name, NF_FLOAT, 5, var5, varid)
!                                      ------------
      if (do_qqjk_inchem) then
        var5(:) = (/ lond(1), latd(1), prsd(1), ydasd(1), recd(1) /)
        call NcDef_variable (ncid_qqjk, 'yda', NF_FLOAT, 5, var5, varid)
!                                        -----
      end if

      call NcSetFill (ncid_qqjk, NF_NOFILL, omode)

      call NcEnd_def (ncid_qqjk)

      return

      end subroutine Define_Netcdf_Out_Qqjk
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Hdr_Qqjk
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Hdr_Qqjk  (gmiDomain, metFields, qj_labels, qk_labels)
!
!USES:
      use m_ncGeneralOpsOutput, only : Write_Netcdf_Hdr_Gen

      implicit none
!
! !INPUT PARAMETERS:
      character (len=MAX_LENGTH_LABELS), intent(in) :: qj_labels(:)
      character (len=MAX_LENGTH_LABELS), intent(in) :: qk_labels(:)
      type(t_gmiDomain), intent(in) :: gmiDomain
      type(t_metFields), intent(in) :: metFields
!
! !DESCRIPTION:
!   Creates some header information for the qj NetCDF output
!   file and writes it out.
!
! !LOCAL VARIABLES:
      integer :: ic
      integer :: cnt1d (1), cnt2d (2)
      integer :: strt1d(1), strt2d(2)

      real*8  :: prsdat(k1:k2)
      real*8  :: qqjsdat(numQjs)
      real*8  :: qqksdat(numQks)
      real*8  :: ydasdat(numActiveChem)
      real*8 , allocatable :: londeg(:), latdeg(:)
      real*8 , allocatable :: am(:), bm(:)
      real*8 :: pt
!EOP
!-----------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Hdr_Qqjk called by ', procID

      allocate(latdeg(ju1_gl:j2_gl))
      call Get_latdeg(gmiDomain, latdeg)

      allocate(londeg(i1_gl:i2_gl))
      call Get_londeg(gmiDomain, londeg)

!     =========================
      call Write_Netcdf_Hdr_Gen  &
!     =========================
     &  (ncid_qqjk, latdeg, londeg, pr_diag, procID, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   hdf_dim_name, lat_dim_name, lon_dim_name)

      allocate(am(k1:k2))
      allocate(bm(k1:k2))

      call Get_pt(metFields, pt)
      call Get_am(metFields, am)
      call Get_bm(metFields, bm)

      prsdat(:) = (am(:) * pt)  + (bm(:) * 1000.0d0)

      strt1d(1) = 1
      cnt1d (1) = numVert

      call Ncwr_1d (prsdat(kx1:kx2), ncid_qqjk, prs_dim_name, strt1d, cnt1d)

      do ic = 1, numQjs
        qqjsdat(ic) = ic
      end do

      strt1d(1) = 1
      cnt1d (1) = numQjs

      call Ncwr_1d (qqjsdat, ncid_qqjk, qqj_dim_name, strt1d, cnt1d)

      do ic = 1, numQks
        qqksdat(ic) = ic
      end do

      strt1d(1) = 1
      cnt1d (1) = numQks

      call Ncwr_1d (qqksdat, ncid_qqjk, qqk_dim_name, strt1d, cnt1d)

!     -----------
!     qqj labels.
!     -----------

      strt2d(:) = (/  1, 1 /)
      cnt2d (:) = (/ 80, numQjs /)

      call Ncwr_2d_Char(qj_labels(1:numQjs), ncid_qqjk, 'qqj_labels', strt2d, cnt2d)


!     -----------
!     qqk labels.
!     -----------

      strt2d(:) = (/  1, 1 /)
      cnt2d (:) = (/ 80, numQks /)

      call Ncwr_2d_Char(qk_labels(1:numQks), ncid_qqjk, 'qqk_labels', strt2d, cnt2d)


      if (do_qqjk_inchem) then

        do ic = 1, numActiveChem
          ydasdat(ic) = ic
        end do

        strt1d(1) = 1
        cnt1d (1) = numActiveChem

        call Ncwr_1d (ydasdat, ncid_qqjk, 'yda_dim', strt1d, cnt1d)

      end if

      return

      end subroutine Write_Netcdf_Hdr_Qqjk
!EOC
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Prep_Netcdf_Qqjk
!
! DESCRIPTION
!   This routine prepares the qj NetCDF output.
!
! ARGUMENTS
!   time_for_nc  : time for NetCDF output?
!
!-----------------------------------------------------------------------------

      subroutine Prep_Netcdf_Qqjk (time_for_nc, Chemistry)


      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      type(t_Chemistry), intent(in) :: Chemistry
      logical :: time_for_nc


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il

      integer, save :: counter = 0

      real*8  :: rcounter
      real*8, allocatable  :: qqjgmi(:,:,:,:)
      real*8, allocatable  :: qqkgmi(:,:,:,:)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Prep_Netcdf_Qqjk called by ', procID
      end if

!     =========================
      if (.not. do_qqjk_inchem) then
!     =========================

!       ============
        if (do_mean) then
!       ============

          if (counter == 0) then

!           ----------------
!           Reset mean sums.
!           ----------------

            qqj_mean(:,:,:,:) = 0.0d0
            qqk_mean(:,:,:,:) = 0.0d0

          end if

          Allocate (qqjgmi(i1:i2, ju1:j2, k1:k2, numQjo))
          call Get_qqjgmi(Chemistry, qqjgmi)

          Allocate (qqkgmi(i1:i2, ju1:j2, k1:k2, numQks))
          call Get_qqkgmi(Chemistry, qqkgmi)

          qqj_mean(:,:,:,:) = qqj_mean(:,:,:,:) + qqjgmi(:,:,kx1:kx2,:)
          qqk_mean(:,:,:,:) = qqk_mean(:,:,:,:) + qqkgmi(:,:,kx1:kx2,:)

          deallocate(qqjgmi)
          deallocate(qqkgmi)

          counter = counter + 1

          if (time_for_nc) then

!           ----------------
!           Calculate means.
!           ----------------

            rcounter = counter

            qqj_mean(:,:,:,:) = qqj_mean(:,:,:,:) / rcounter
            qqk_mean(:,:,:,:) = qqk_mean(:,:,:,:) / rcounter

            counter = 0

          end if

!       ======
        end if
!       ======

!     ======
      end if
!     ======

      return

      end subroutine Prep_Netcdf_Qqjk

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Write_Netcdf_Qqjk
!
! DESCRIPTION
!   This routine writes the qj NetCDF output.
!
! ARGUMENTS
!   None
!
!-----------------------------------------------------------------------------

      subroutine Write_Netcdf_Qqjk (nymd, nhms, gmi_sec)
!     
      implicit none
!     
! !INPUT PARAMETERS:
      integer, intent(in) :: nymd, nhms
      real*8 , intent(in) :: gmi_sec
!
! !DESCRIPTION:
! Writes the qqjk NetCDF output.
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_QQJ1_NC = 3900
      integer, parameter :: SG_QQK1_NC = 3901
      integer, parameter :: SG_YDA1_NC = 3902
!     
! !LOCAL VARIABLES:
      integer :: cnt2d (2)
      integer :: strt2d(2)
      integer :: hdr(NETCDF_HDF)
      real*8 , allocatable :: qqjGlob(:,:,:)
      real*8 , allocatable :: qqkGlob(:,:,:)
      real*8 , allocatable :: ydaGlob(:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Qqjk called by ', procID

      if (iAmRootProc) then

         strt2d(:) = (/ 1, rnum_out_qqjk /)

         cnt2d (:) = (/ NETCDF_HDF, 1 /)

         hdr(1) = Nint (gmi_sec)
         hdr(2) = nymd
         hdr(3) = nhms

         call Ncwr_2d_Int (hdr, ncid_qqjk, hdr_var_name, strt2d, cnt2d)

         allocate(qqjGlob(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))
         allocate(qqkGlob(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))
         if (do_qqjk_inchem) allocate(ydaGlob(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))

      end if

      call subDomain2Global (qqjGlob, qqj_nc, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &        i1, i2, ju1, j2, kx1, kx2, rootProc, procID, map1_u, numDomains, &
     &        SG_QQJ1_NC, commuWorld)

      call subDomain2Global (qqkGlob, qqk_nc, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &        i1, i2, ju1, j2, kx1, kx2, rootProc, procID, map1_u, numDomains, &
     &        SG_QQK1_NC, commuWorld)

      if (do_qqjk_inchem) then

         call subDomain2Global (ydaGlob, yda_nc, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &          i1, i2, ju1, j2, kx1, kx2, rootProc, procID, map1_u, numDomains, &
     &          SG_YDA1_NC, commuWorld)

      end if

      if (iAmRootProc) then
         deallocate(qqjGlob)
         deallocate(qqkGlob)
         if (do_qqjk_inchem) deallocate(ydaGlob)

         call Ncdo_Sync (ncid_qqjk)
      end if

      return

      end subroutine Write_Netcdf_Qqjk
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutput_Qqjk
!
! !INTERFACE:
!
      subroutine bufferOutput_Qqjk (Chemistry)
!
! !USES:
      use GmiBuffer_mod, only : bufferOutput4d, bufferOutput4d_Nomean
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_Chemistry), intent(in) :: Chemistry
!
! !DESCRIPTION:
! Buffers the qqjk data to the output file.
! The worker processors send a slice of the data to the root processor, 
! the root processor writes it out, they then send another slice and it 
! is written out, and so on.
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_QQJ2_NC = 3903
      integer, parameter :: SG_QQK2_NC = 3904
      integer, parameter :: SG_YDA2_NC = 3905
      character (len=MAX_LENGTH_VAR_NAME), parameter ::  YDA_VNAM = 'yda'
!
! !LOCAL VARIABLES:
      real*8, allocatable ::    yda(:,:,:,:)
      real*8, allocatable ::  qqjda(:,:,:,:)
      real*8, allocatable ::  qqkda(:,:,:,:)
      real*8, allocatable :: qqjgmi(:,:,:,:)
      real*8, allocatable :: qqkgmi(:,:,:,:)
!
!EOP
!----------------------------------------------------------------------------------------
!BOC

      if (numQjs > 0) then

        if (do_qqjk_inchem) then

           Allocate(yda  (i1:i2, ju1:j2, k1:k2, numActiveChem))
           Allocate(qqjda(i1:i2, ju1:j2, k1:k2, numQjs))
           Allocate(qqkda(i1:i2, ju1:j2, k1:k2, numQks))

           call Get_yda  (Chemistry, yda)
           call Get_qqjda(Chemistry, qqjda)
           call Get_qqkda(Chemistry, qqkda)
           
!         ===========================
          call bufferOutput4d_Nomean  &
!         ===========================
     &      (qqj_var_name, kx1, kx2, commuWorld, SG_QQJ2_NC,  &
     &       ncid_qqjk, numQjs, rnum_out_qqjk, qqjda, qqj_nc, &
     &       map1_u, numDomains, rootProc, procID, &
     &       i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)

!         ===========================
          call bufferOutput4d_Nomean  &
!         ===========================
     &      (qqk_var_name, kx1, kx2, commuWorld, SG_QQK2_NC,  &
     &       ncid_qqjk, numQks, rnum_out_qqjk, qqkda, qqk_nc, &
     &       map1_u, numDomains, rootProc, procID, &
     &       i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)

!         ===========================
          call bufferOutput4d_Nomean  &
!         ===========================
     &      (YDA_VNAM, kx1, kx2, commuWorld, SG_YDA2_NC,  &
     &       ncid_qqjk, numActiveChem, rnum_out_qqjk, yda, yda_nc, &
     &       map1_u, numDomains, rootProc, procID, &
     &       i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)

        else

           if (.not. do_mean) then
              Allocate(qqjgmi(i1:i2, ju1:j2, k1:k2, numQjo))
              Allocate(qqkgmi(i1:i2, ju1:j2, k1:k2, numQks))

              call Get_qqjgmi(Chemistry, qqjgmi)
              call Get_qqkgmi(Chemistry, qqkgmi)
           end if

!         ====================
          call bufferOutput4d  &
!         ====================
     &      (qqj_var_name, do_mean, kx1, kx2, commuWorld,  &
     &       SG_QQJ2_NC, ncid_qqjk, numQjs, rnum_out_qqjk, qqjgmi,  &
     &       qqj_mean, qqj_nc, map1_u, numDomains, rootProc, procID, &
     &       i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)

!         ====================
          call bufferOutput4d  &
!         ====================
     &      (qqk_var_name, do_mean, kx1, kx2, commuWorld,  &
     &       SG_QQK2_NC, ncid_qqjk, numQks, rnum_out_qqjk, qqkgmi,  &
     &       qqk_mean, qqk_nc, map1_u, numDomains, rootProc, procID, &
     &       i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)

           if (.not. do_mean) then
              deallocate(qqjgmi)
              deallocate(qqkgmi)
           end if
        end if
      end if

      return

      end subroutine bufferOutput_Qqjk
!EOC
!-----------------------------------------------------------------------------
end module GmiControlQqjk_mod
