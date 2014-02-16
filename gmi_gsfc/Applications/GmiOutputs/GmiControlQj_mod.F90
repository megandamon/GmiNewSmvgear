!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlQj_mod
!
! !INTERFACE:
!
module GmiControlQj_mod
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
      use GmiChemistryMethod_mod, only : t_Chemistry, Get_qjgmi, Get_num_qjo,  &
     &       Get_qj_labels
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      use GmiDomainDecomposition_mod, only : t_gmiDomain
      use GmiDomainDecomposition_mod, only : Get_rootProc, Get_procID, Get_iAmRootProc
      use GmiDomainDecomposition_mod, only : Get_numDomains
      use GmiDomainDecomposition_mod, only : Get_communicatorWorld
      use GmiDomainDecomposition_mod, only : Get_map1_u
      use GmiDomainDecomposition_mod, only : Get_latdeg, Get_londeg
      use GmiGrid_mod              , only : t_gmiGrid
      use GmiGrid_mod              , only : Get_i1, Get_i2, Get_ju1, Get_j2, Get_k1, Get_k2
      use GmiGrid_mod              , only : Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl
      use GmiGrid_mod              , only : Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jhi_gl
      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Get_hdf_dim_name,    &
     &       Get_lon_dim_name, Get_lat_dim_name
      use GmiDiagnosticsMethod_mod, only : Get_qjOutputFrequency, Get_hdr_var_name, Get_rec_dim_name
      use GmiDiagnosticsMethod_mod, only : Get_prs_dim_name, Get_qj_var_name, Get_qj_dim_name
      use GmiDiagnosticsMethod_mod, only : Get_do_mean, Get_k1r_gl, Get_k2r_gl
      use GmiDiagnosticsMethod_mod, only : Get_pr_qj, Get_pr_diag, Get_problem_name

      use GmiTimeControl_mod, only : t_GmiClock, GmiSplitDateTime, isOutTime, &
     &       Get_curGmiDate, Get_curGmiTime, Get_gmiSeconds, Get_gmiTimeStep, &
     &       Get_numTimeSteps

      use GmiMetFieldsControl_mod, only : t_metFields, Get_pt,    &
     &       Get_metdata_name, Get_mdt, Get_ai, Get_bi, Get_am, Get_bm
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: controlOutputQj, initializeOutputQj, finalizeOutputQj

#     include "gmi_phys_constants.h"
#     include "GmiParameters.h"
!
! !DESCRIPTION:
! Contains the necessary routines for the qj outputs. 
! Three routines are visible from the outside:
! %
! \begin{description}
! \item[initializeOutputQj:] called once in the initialization steps to
!      open the netCDF file, define variables in the file, write out header
!      information and allocate variables.
! \item[controlOutputQj:] called at the end of each model time step to
!      update variables and if needed for communications and writing out
!      data in the file.
! \item[finalizeOutputQj:] called at the end of the integration to
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
      real*8 , pointer, save :: qj_nc   (:,:,:) => null() ! photolysis rate constants (s^-1)
!     
      ! Variables manipulate by worker processors only         
      real*8 , pointer, save :: qj_mean (:,:,:,:) => null()
!
      ! netCDF file identifier
      integer,          save :: ncid_qj

      ! Counter for the number of records
      integer,          save :: rnum_out_qj

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

      integer,          save :: numQjOutrecs  ! # qjs for output file
      
      logical,          save :: iAmRootProc
      logical,          save :: pr_diag
      integer,          save :: procID, rootProc
      integer,          save :: numDomains
      integer,          save :: commuWorld
      integer, pointer, save :: map1_u       (:,:,:) => null()
      integer, pointer, save :: ewflag           (:) => null()
      integer, pointer, save :: nspole           (:) => null()
      character (len=MAX_LENGTH_VAR_NAME), save :: qj_dim_name, rec_dim_name, prs_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: qj_var_name, hdr_var_name, hdf_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: lon_dim_name, lat_dim_name
      logical           , save :: pr_qj, do_mean
      real*8            , save :: qjOutputFrequency
      real*8            , save :: mdt

      integer,          save :: nhdf
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
! !IROUTINE: initializeOutputQj
!
! !INTERFACE:
!
      subroutine initializeOutputQj(gmiGrid, gmiDomain, Diagnostics, &
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
      integer :: in1, k1r_gl, k2r_gl
      integer, allocatable :: ewflag_org(:)
      integer, allocatable :: nspole_org(:)
      character (len=MAX_LENGTH_FILE_NAME) :: fname
      character (len=MAX_LENGTH_FILE_NAME) :: problem_name
      integer  :: num_qjo
      character (len=MAX_LENGTH_LABELS) :: qj_labels(MAX_NUM_QJ)
!EOP
!-----------------------------------------------------------------------------
!BOC

      !#####################################################
      ! Initialization of variables used in all the routines
      !#####################################################

      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write (6,*) 'initializeOutputQj called by ', procID

      call Get_num_qjo(Chemistry, num_qjo)
      numQjOutrecs  = num_qjo

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
      call Get_qjOutputFrequency(Diagnostics, qjOutputFrequency)
      call Get_pr_qj       (Diagnostics, pr_qj )
      call Get_hdr_var_name(Diagnostics, hdr_var_name)
      call Get_rec_dim_name(Diagnostics, rec_dim_name)
      call Get_prs_dim_name(Diagnostics, prs_dim_name)
      call Get_qj_var_name (Diagnostics, qj_var_name)
      call Get_qj_dim_name (Diagnostics, qj_dim_name)
      call Get_hdf_dim_name (Diagnostics, hdf_dim_name)
      call Get_lat_dim_name (Diagnostics, lat_dim_name)
      call Get_lon_dim_name (Diagnostics, lon_dim_name)

      call Get_mdt (metFields, mdt)

      if (iAmRootProc) then

         !###########################
         ! Initialize the output file
         !###########################

         ! Determine the file name
         call Get_problem_name(Diagnostics, problem_name)
         call makeOutfileName (fname, '.qj.nc', problem_name)

         ! Create the file and assign a file identifier
         call Nccr_Wr (ncid_qj, fname)

         ! Define the variables in the file
         call Define_Netcdf_Out_Qj( )

         ! Write header data
         call Get_qj_labels(Chemistry, qj_labels)

         call Write_Netcdf_Hdr_Qj (gmiDomain, metFields, qj_labels)

         call Ncdo_Sync (ncid_qj)

         ! Initialize the counter for the number of records
         rnum_out_qj = 1
      end if

      !########################
      ! Allocation of variables
      !########################
      call allocateVariablesQj()

      return

      end subroutine initializeOutputQj
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: controlOutQj
!
! !INTERFACE:
!
      subroutine controlOutputQj (last_tstp, Chemistry, gmiClock)
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
      if (pr_diag) Write (6,*) 'controlOutputQj called by ', procID

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
     &        month_save, month, day, nhms, ndt, gmi_sec, qjOutputFrequency)

      month_save = month

      if (last_tstp .and. (.not. time_for_nc)) then
        if (Mod (num_time_steps, Nint (mdt) / ndt) == 0) then

          time_for_nc = .true.

        end if
      end if

      if (time_for_nc .or. do_mean) then
         !====================
         call Prep_Netcdf_Qj  (time_for_nc, Chemistry)
         !====================
      end if


!     ================
      if (time_for_nc) then
!     ================

!        =====================
         call Write_Netcdf_Qj (nymd, nhms, gmi_sec)
!        =====================

!        ======================
         call bufferOutput_Qj (Chemistry)
!        ======================

         if (iAmRootProc) then
             rnum_out_qj = rnum_out_qj + 1
         end if

!     ======
      end if
!     ======

      return

      end subroutine controlOutputQj
!EOC
!-----------------------------------------------------------------------------
!BOP
!     
! !IROUTINES: finalizeOutputQj
!
! !INTERFACE:
!
      subroutine finalizeOutputQj()
!
      implicit none
! 
! !DESCRIPTION:
! Deallocates variables necessary to produce qj outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'finalizeOutputQj called by ', procID

      if (pr_qj) then
     
         deallocate (qj_nc)

         if (do_mean) then
            deallocate (qj_mean)
         end if

         if (iAmRootProc) then
            call Nccl_Noerr (ncid_qj)
         end if
      end if

      return

      end subroutine finalizeOutputQj
!EOC
!-----------------------------------------------------------------------------
!BOP
!     
! !IROUTINES: allocateVariablesQj
!
! !INTERFACE:
!
      subroutine allocateVariablesQj()
!
      implicit none
! 
! !DESCRIPTION:
! Allocates variables necessary to produce qj outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'allocateVariablesQj called by ', procID

      if (pr_qj) then

         Allocate (qj_nc(i1:i2, ju1:j2, kx1:kx2))
         qj_nc = 0.0d0

         if (do_mean) then
            Allocate (qj_mean(i1:i2, ju1:j2, kx1:kx2, numQjOutrecs))
            qj_mean = 0.0d0
         end if

      end if

      return

      end subroutine allocateVariablesQj
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Define_Netcdf_Out_Qj
!
! !INTERFACE:
!
      subroutine Define_Netcdf_Out_Qj ()
!
! !USES:
      use m_ncGeneralOpsOutput, only : Define_Netcdf_Out_Gen
!
      implicit none
!
#     include "netcdf.inc"
!
! !DESCRIPTION:
! Makes the necessary definitions for the qj NetCDF output file.
!
! !LOCAL VARIABLES:
      integer :: ierr
      integer :: omode
      integer :: varid
      integer :: chrd(1)
      integer :: hdfd(1)
      integer :: lond(1), latd(1), prsd(1)
      integer :: qjsd(1), recd(1)
      integer :: var2(2)
      integer :: var5(5)
!EOC
!-----------------------------------------------------------------------------
!BOP
      if (pr_diag) Write (6,*) 'Define_Netcdf_Out_Qj called by ', procID

!     ==========================
      call Define_Netcdf_Out_Gen  &
!     ==========================
     &  (ncid_qj, nhdf, numLon, numLat, numVert, hdfd, lond, latd, prsd,  &
     &   hdf_dim_name, lon_dim_name, lat_dim_name, prs_dim_name)


!     ------------------
!     Define dimensions.
!     ------------------

      call NcDef_dimension(ncid_qj, qj_dim_name, numQjOutrecs, qjsd(1))
!                                    -----------
      call NcDef_dimension (ncid_qj, 'chr_dim', 80, chrd(1))
!                                     ---------
      call NcDef_dimension(ncid_qj, rec_dim_name, NF_UNLIMITED, recd(1))
!                                    ------------

!     -----------------------------------------
!     Define variables and variable attributes.
!     -----------------------------------------

      call NcDef_variable (ncid_qj, qj_dim_name, NF_FLOAT, 1, qjsd, varid)
!                                    -----------
      call NcDef_var_attributes (ncid_qj, varid, 'long_name', 'Qjs')
      call NcDef_var_attributes (ncid_qj, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_qj, varid, 'coord_labels', 'qj_labels')
      call NcDef_var_attributes (ncid_qj, varid, 'selection_category', 'NULL')

      var2(:) = (/ hdfd(1), recd(1) /)
      call NcDef_variable (ncid_qj, hdr_var_name, NF_INT, 2, var2, varid)
!                                    ------------

      var2(:) = (/ chrd(1), qjsd(1) /)
      call NcDef_variable (ncid_qj, 'qj_labels', NF_CHAR, 2, var2, varid)
!                                    -----------
      call NcDef_var_attributes (ncid_qj, varid, 'long_name', 'qj name')
      call NcDef_var_attributes (ncid_qj, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_qj, varid, 'selection_category', 'NULL')

      var5(:) = (/ lond(1), latd(1), prsd(1), qjsd(1), recd(1) /)
      call NcDef_variable (ncid_qj, qj_var_name, NF_FLOAT, 5, var5, varid)
!                                    -----------

      call NcSetFill (ncid_qj, NF_NOFILL, omode)

      call NcEnd_def (ncid_qj)

      return

      end subroutine Define_Netcdf_Out_Qj
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Hdr_Qj
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Hdr_Qj  (gmiDomain, metFields, qj_labels)
!
! !USES:
      use m_ncGeneralOpsOutput, only : Write_Netcdf_Hdr_Gen
!
      implicit none
!
! !INPUT PARAMETERS:
      character (len=MAX_LENGTH_LABELS), intent(in) :: qj_labels(:)
      type(t_gmiDomain), intent(in) :: gmiDomain
      type(t_metFields), intent(in) :: metFields
!
! !DESCRIPTION:
! Creates some header information for the qj NetCDF output file and writes 
!  it out.
!
! !LOCAL VARIABLES:
      integer :: ic
      integer :: cnt1d (1), cnt2d (2)
      integer :: strt1d(1), strt2d(2)
      real*8  :: prsdat(k1:k2)
      real*8, allocatable  :: prsdat_r(:)
      real*8  :: qjsdat(numQjOutrecs)
      real*8 , allocatable :: londeg(:), latdeg(:)
      real*8 , allocatable :: am(:), bm(:)
      real*8 :: pt
!EOP
!-----------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Hdr_Qj called by ', procID

      allocate(latdeg(ju1_gl:j2_gl))
      call Get_latdeg(gmiDomain, latdeg)

      allocate(londeg(i1_gl:i2_gl))
      call Get_londeg(gmiDomain, londeg)

!     =========================
      call Write_Netcdf_Hdr_Gen  &
!     =========================
     &  (ncid_qj, latdeg, londeg, pr_diag, procID, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   hdf_dim_name, lat_dim_name, lon_dim_name)

      allocate(am(k1:k2))
      allocate(bm(k1:k2))

      call Get_pt(metFields, pt)
      call Get_am(metFields, am)
      call Get_bm(metFields, bm)

      prsdat(:) = (am(:) * pt)  + (bm(:) * 1000.0d0)

      strt1d(1) = 1
      cnt1d (1) = numVert

      call Ncwr_1d (prsdat(kx1:kx2), ncid_qj, prs_dim_name, strt1d, cnt1d)

      do ic = 1, numQjOutrecs
        qjsdat(ic) = ic
      end do

      strt1d(1) = 1
      cnt1d (1) = numQjOutrecs

      call Ncwr_1d (qjsdat, ncid_qj, qj_dim_name, strt1d, cnt1d)

!     ----------
!     qj labels.
!     ----------

      strt2d(:) = (/  1, 1 /)
      cnt2d (:) = (/ 80, numQjOutrecs /)

      call Ncwr_2d_Char (qj_labels, ncid_qj, 'qj_labels', strt2d, cnt2d)

      return

      end subroutine Write_Netcdf_Hdr_Qj
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Prep_Netcdf_Qj
!
! !INTERFACE:
!
      subroutine Prep_Netcdf_Qj (time_for_nc, Chemistry)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical,           intent(in) :: time_for_nc ! time for NetCDF output?
      type(t_Chemistry), intent(in) :: Chemistry
!
! !DESCRIPTION:
! Prepares the qj NetCDF output.
!
! !LOCAL VARIABLES:
      integer :: il
      integer, save :: counter = 0
      real*8  :: rcounter
      real*8, allocatable  :: qjgmiLoc(:,:,:,:)
!EOC
!-----------------------------------------------------------------------------
!BOP
      if (pr_diag) Write (6,*) 'Prep_Netcdf_Qj called by ', procID

!     ============
      if (do_mean) then
!     ============

        if (counter == 0) then

!         ----------------
!         Reset mean sums.
!         ----------------

          qj_mean(:,:,:,:) = 0.0d0

        end if

!       ------------------------
!       Update mean sums.
!       ------------------------
        allocate(qjgmiLoc(i1:i2, ju1:j2, k1:k2, numQjOutrecs))
        call Get_qjgmi(Chemistry, qjgmiLoc)

        qj_mean(:,:,:,:) = qj_mean(:,:,:,:) + qjgmiLoc(:,:,kx1:kx2,:)

        deallocate(qjgmiLoc)

        counter = counter + 1

        if (time_for_nc) then

!         ----------------
!         Calculate means.
!         ----------------

          rcounter = counter

          qj_mean(:,:,:,:) = qj_mean(:,:,:,:) / rcounter

          counter = 0

        end if

      end if

      return

      end subroutine Prep_Netcdf_Qj
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Qj
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Qj (nymd, nhms, gmi_sec)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: nymd, nhms
      real*8 , intent(in) :: gmi_sec
!
! !DESCRIPTION:
! Writes the qj NetCDF output.
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_QJ1_NC = 2802
!
! !LOCAL VARIABLES:
      integer :: cnt2d (2)
      integer :: strt2d(2)
      integer :: hdr(NETCDF_HDF)
      real*8 , allocatable :: qjGlob(:,:,:)
!EOP
!-------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Qj called by ', procID

      if (iAmRootProc) then
         strt2d(:) = (/ 1, rnum_out_qj /)

         cnt2d (:) = (/ NETCDF_HDF, 1 /)

         hdr(1) = Nint (gmi_sec)
         hdr(2) = nymd
         hdr(3) = nhms

         call Ncwr_2d_Int (hdr, ncid_qj, hdr_var_name, strt2d, cnt2d)

         allocate(qjGlob(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))
      end if 

      call subDomain2Global (qjGlob, qj_nc, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &        i1, i2, ju1, j2, kx1, kx2, rootProc, procID, map1_u, numDomains, &
     &        SG_QJ1_NC, commuWorld)

      if (iAmRootProc) then
         deallocate(qjGlob)

         call Ncdo_Sync (ncid_qj)
      end if 

      return

      end subroutine Write_Netcdf_Qj
!
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutput_Qj
!
! !INTERFACE:
!
      subroutine bufferOutput_Qj (Chemistry)
!
! !USES:
      use GmiBuffer_mod, only : bufferOutput4d
!
      implicit none

! !INPUT PARAMETERS:
      type(t_Chemistry), intent(in) :: Chemistry
!
! !DESCRIPTION:
! Buffers the qj data to the output file.  The worker PEs send a slice of the 
! data to the root PE, the root PE writes it out, they then send another slice 
! and it is written out, and so on.
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_QJ2_NC = 2801
!
! !LOCAL VARIABLES:
      real*8, allocatable :: qjgmi(:,:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (numQjOutrecs > 0) then

         if (.not. do_mean) then
            Allocate (qjgmi(i1:i2, ju1:j2, kx1:kx2, numQjOutrecs))
            call Get_qjgmi(Chemistry, qjgmi)
         end if

!       ====================
        call bufferOutput4d  &
!       ====================
     &    (qj_var_name, do_mean, kx1, kx2, commuWorld,  &
     &     SG_QJ2_NC, ncid_qj, numQjOutrecs, rnum_out_qj, qjgmi,  &
     &     qj_mean, qj_nc, map1_u, numDomains, rootProc, procID, &
     &     i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)

         if (.not. do_mean) then
            deallocate (qjgmi)
         end if

      end if

      return

      end subroutine bufferOutput_Qj
!EOC
!-----------------------------------------------------------------------------
end module GmiControlQj_mod
