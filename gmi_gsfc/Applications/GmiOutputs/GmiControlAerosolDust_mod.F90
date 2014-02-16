!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlAerosolDust_mod
!
! !INTERFACE:
!
module GmiControlAerosolDust_mod
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
      use GmiChemistryMethod_mod    , only : t_Chemistry

      use GmiDomainDecomposition_mod, only : t_gmiDomain
      use GmiDomainDecomposition_mod, only : Get_rootProc, Get_procID, Get_iAmRootProc
      use GmiDomainDecomposition_mod, only : Get_numDomains
      use GmiDomainDecomposition_mod, only : Get_communicatorWorld
      use GmiDomainDecomposition_mod, only : Get_map1_u
      use GmiDomainDecomposition_mod, only : Get_mcorGlob, Get_latdeg, Get_londeg
      use GmiGrid_mod               , only : t_gmiGrid, Get_numSpecies
      use GmiGrid_mod               , only : Get_i1, Get_i2, Get_ju1, Get_j2, Get_k1, Get_k2
      use GmiGrid_mod               , only : Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl
      use GmiGrid_mod               , only : Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jhi_gl
      use GmiGrid_mod               , only : Get_ilo, Get_ihi, Get_julo, Get_jhi
      use GmiDiagnosticsMethod_mod, only : t_Diagnostics
      use GmiDiagnosticsMethod_mod, only : Get_hdr_var_name, Get_rec_dim_name
      use GmiDiagnosticsMethod_mod, only : Get_hdf_dim_name, Get_lat_dim_name, Get_lon_dim_name
      use GmiDiagnosticsMethod_mod, only : Get_prs_dim_name, Get_AerDust_var_name
      use GmiDiagnosticsMethod_mod, only : Get_do_mean, Get_k1r_gl, Get_k2r_gl
      use GmiDiagnosticsMethod_mod, only : Get_pr_AerDust, Get_pr_diag
      use GmiDiagnosticsMethod_mod, only : Get_outaerdust_name, Get_problem_name
      use GmiDiagnosticsMethod_mod, only : Get_pr_psf, Get_aerdustOutputFrequency

      use GmiTimeControl_mod, only : t_GmiClock, GmiSplitDateTime, isOutTime, &
     &       Get_curGmiDate, Get_curGmiTime, Get_gmiSeconds, Get_gmiTimeStep, &
     &       Get_numTimeSteps

      use GmiMetFieldsControl_mod, only : t_metFields, Get_pt, Get_mdt, Get_ai,&
     &       Get_bi, Get_am, Get_bm, Get_gridBoxHeight, Get_pctm2
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: controlOutputAerosolDust, initializeOutputAerosolDust
      public  :: finalizeOutputAerosolDust

#     include "gmi_phys_constants.h"
#     include "GmiParameters.h"
#     include "setkin_par.h"
#     include "gmi_AerDust_const.h"
!
! !DESCRIPTION:
! Contains the necessary routines for the species concentration outputs. 
! Three routines are visible from the outside:
! %
! \begin{description}
! \item[initializeOutputAerosolDust:] called once in the initialization steps to
!      open the netCDF file, define variables in the file, write out header
!      information and allocate variables.
! \item[controlOutputAerosolDust:] called at the end of each model time step to
!      update variables and if needed for communications and writing out
!      data in the file.
! \item[finalizeOutputAerosolDust:] called at the end of the integration to
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
      real*8 , pointer, save :: psf_AerDust_nc      (:,:) => null()
      real*8 , pointer, save :: optDepth        (:,:,:,:) => null()
      real*8 , pointer, save :: optDepth_nc       (:,:,:) => null()
      real*8 , pointer, save :: CM_AerDust_nc     (:,:,:) => null()
      real*8 , pointer, save :: grid_height_nc    (:,:,:) => null()
!     
      ! Variables manipulate by worker processors only         
      real*8 , pointer, save :: psf_AerDust_mean    (:,:) => null()

      character (len=MAX_LENGTH_LABELS) :: AerDust_labels  (num_AerDust)
      character (len=MAX_LENGTH_LABELS) :: CMAerDust_labels(num_CM_AerDust)
!
      ! netCDF file identifier
      integer,          save :: ncid_AerDust

      ! Counter for the number of records
      integer,          save :: rnum_out_AerDust

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

      logical,          save :: iAmRootProc
      logical,          save :: pr_diag
      integer,          save :: procID, rootProc
      integer,          save :: numDomains
      integer,          save :: commuWorld
      integer, pointer, save :: map1_u       (:,:,:) => null()

      integer,          save :: nhdf
      character (len=MAX_LENGTH_VAR_NAME), save :: rec_dim_name, prs_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: hdf_dim_name, lon_dim_name, lat_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: hdr_var_name, AerDust_var_name
      logical           , save :: pr_AerDust, do_mean, pr_psf
      real*8            , save :: aerdustOutputFrequency
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
! !IROUTINE: initializeOutputAerosolDust
!
! !INTERFACE:
!
      subroutine initializeOutputAerosolDust(gmiGrid, gmiDomain, Diagnostics, &
     &                     metFields)
!
! !USES:
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_gmiGrid    ), intent(in) :: gmiGrid
      type(t_gmiDomain  ), intent(in) :: gmiDomain
      type(t_metFields  ), intent(in) :: metFields
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
      character (len=MAX_LENGTH_FILE_NAME) :: outaerdust_name
!EOP
!-----------------------------------------------------------------------------
!BOC

      !#####################################################
      ! Initialization of variables used in all the routines
      !#####################################################

      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write (6,*) 'initializeOutputAerosolDust called by ', procID

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

      allocate(map1_u(2, 2, numDomains))
      call Get_map1_u (gmiDomain, map1_u)

      call Get_k1r_gl(Diagnostics, k1r_gl)
      call Get_k2r_gl(Diagnostics, k2r_gl)

      kx1     = k1r_gl
      kx2     = k2r_gl

      numLon  = i2_gl - i1_gl + 1
      numLat  = j2_gl - ju1_gl + 1
      numVert = kx2 - kx1 + 1

      nhdf = NETCDF_HDF

      call Get_do_mean     (Diagnostics, do_mean)
      call Get_aerdustOutputFrequency(Diagnostics, aerdustOutputFrequency)
      call Get_pr_AerDust  (Diagnostics, pr_AerDust)
      call Get_pr_psf      (Diagnostics, pr_psf    )
      call Get_hdf_dim_name (Diagnostics, hdf_dim_name)
      call Get_lat_dim_name (Diagnostics, lat_dim_name)
      call Get_lon_dim_name (Diagnostics, lon_dim_name)
      call Get_hdr_var_name(Diagnostics, hdr_var_name)
      call Get_rec_dim_name(Diagnostics, rec_dim_name)
      call Get_prs_dim_name(Diagnostics, prs_dim_name)
      call Get_AerDust_var_name(Diagnostics, AerDust_var_name)

      call Get_mdt(metFields, mdt)

      if (pr_AerDust) then

         if (iAmRootProc) then

            !###########################
            ! Initialize the output file
            !###########################

            ! Determine the file name
            call Get_problem_name   (Diagnostics, problem_name)
            call Get_outaerdust_name(Diagnostics, outaerdust_name)

            in1 = Len_Trim (outaerdust_name)
            call makeOutfileName (fname, '.' // outaerdust_name(1:in1) // '.nc', &
     &                            problem_name)

            ! Create the file and assign a file identifier
            call Nccr_Wr (ncid_AerDust, fname)

            ! Define the variables in the file
            call Define_Netcdf_Out_AerDust( )

            ! Write header data

            call Write_Netcdf_Hdr_AerDust (gmiDomain, metFields)

            call Ncdo_Sync (ncid_AerDust)

            ! Initialize the counter for the number of records
            rnum_out_AerDust = 1
         end if

         !########################
         ! Allocation of variables
         !########################
         call allocateVariablesAerosolDust()
      end if

      return

      end subroutine initializeOutputAerosolDust
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: controlOutAerosolDust
!
! !INTERFACE:
!
      subroutine controlOutputAerosolDust (last_tstp, Chemistry, gmiClock, &
     &                  metFields)
!
! !USES:
      use GmiTimeControl_mod, only : GmiSplitDateTime

      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: last_tstp
      type(t_gmiClock) , intent(in) :: gmiClock
      type(t_metFields) , intent(in) :: metFields
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_Chemistry) , intent(inOut) :: Chemistry
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
      if (pr_diag) Write (6,*) 'controlOutputAerosolDust called by ', procID

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
     &          month_save, month, day, nhms, ndt, gmi_sec, aerdustOutputFrequency)

      month_save = month

      if (last_tstp .and. (.not. time_for_nc)) then
        if (Mod (num_time_steps, Nint (mdt) / ndt) == 0) then

          time_for_nc = .true.

        end if
      end if

      if (time_for_nc .or. do_mean) then
         !====================
          call Prep_Netcdf_AerDust (time_for_nc, Chemistry, metFields)
          !====================
      endif

!     ================
      if (time_for_nc) then
!     ================

!        =====================
         call Write_Netcdf_AerDust (nymd, nhms, gmi_sec)
!        =====================

!        ======================
         call bufferOutput_AerDust ( )
!        ======================

        if (iAmRootProc) then
            rnum_out_AerDust = rnum_out_AerDust + 1
        end if

!     ======
      end if
!     ======

      return

      end subroutine controlOutputAerosolDust
!EOC
!-----------------------------------------------------------------------------
!BOP
!     
! !IROUTINES: finalizeOutputAerosolDust
!
! !INTERFACE:
!
      subroutine finalizeOutputAerosolDust()
!
      implicit none
! 
! !DESCRIPTION:
! Deallocates variables necessary to produce aerosol/dust outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'finalizeOutputAerosolDust called by ', procID

      !====================================================
      ! Deallocation on both the master and worker processors
      !====================================================

      deallocate (OptDepth)

      deallocate (OptDepth_nc)
      deallocate (CM_AerDust_nc)
      deallocate (grid_height_nc)
      deallocate (psf_AerDust_nc)

      if (do_mean) then
         if (pr_psf) then
            deallocate (psf_AerDust_mean)
         end if
      end if

      if (iAmRootProc) call Nccl_Noerr (ncid_AerDust)

      return

      end subroutine finalizeOutputAerosolDust
!EOC
!-----------------------------------------------------------------------------
!BOP
!     
! !IROUTINES: allocateVariablesAerosolDust
!
! !INTERFACE:
!
      subroutine allocateVariablesAerosolDust()
!
      implicit none
! 
! !DESCRIPTION:
! Allocates variables necessary to produce aerosol/dust outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'allocateVariablesAerosolDust called by ', procID

      Allocate (optDepth(i1:i2, ju1:j2, kx1:kx2, num_AerDust))

      Allocate (OptDepth_nc(i1:i2, ju1:j2, kx1:kx2))
      OptDepth_nc = 0.0d0
   
      Allocate (CM_AerDust_nc(i1:i2, ju1:j2,5))
      Allocate (grid_height_nc(i1:i2, ju1:j2,kx1:kx2))
      CM_AerDust_nc = 0.0d0
      grid_height_nc= 0.0d0

      Allocate (psf_AerDust_nc (i1:i2, ju1:j2))
      psf_AerDust_nc = 0.0d0

      if (do_mean) then

         if (pr_psf) then
            Allocate (psf_AerDust_mean (i1:i2, ju1:j2))
            psf_AerDust_mean = 0.0d0
         end if
      end if

      return

      end subroutine allocateVariablesAerosolDust
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Define_Netcdf_Out_AerDust
!
! DESCRIPTION
!   This routine makes the necessary definitions for the aersosol/dust
!   NetCDF output file.
!

      subroutine Define_Netcdf_Out_AerDust  ()
!
! !USES:
      use m_ncGeneralOpsOutput, only : Define_Netcdf_Out_Gen

      implicit none

#     include "netcdf.inc"

!     ----------------------
!     Variable declarations.
!     ----------------------
      character (len=MAX_LENGTH_VAR_NAME) :: prsp1_dim_name
      integer :: ierr , pos1
      integer :: omode
      integer :: varid

      integer :: chrd(1)
      integer :: hdfd(1)
      integer :: lond(1), latd(1), prsd(1), prsp1d(1)
      integer :: spcd(1), recd(1)
      integer :: spcd_2(1)

      integer :: var2(2), var3(3), var4(4), var5(5)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Define_Netcdf_Out_AerDust called by ', procID
      end if

!     ==========================
      call Define_Netcdf_Out_Gen  &
!     ==========================
     &  (ncid_AerDust, nhdf, numLon, numLat, numVert, hdfd, lond, latd, prsd,  &
     &   hdf_dim_name, lon_dim_name, lat_dim_name, prs_dim_name)


!     ------------------
!     Define dimensions.
!     ------------------

      prsp1_dim_name = prs_dim_name
      pos1 = Len_Trim (prs_dim_name) + 1
      prsp1_dim_name(pos1:pos1+2) = 'p1'

      call NcDef_dimension(ncid_AerDust, prsp1_dim_name, numVert+1, prsp1d(1))

      call NcDef_dimension(ncid_AerDust,'AerDust_dim', num_AerDust, spcd(1))
!                                         -----------
      call NcDef_dimension  &
     &        (ncid_AerDust,'CM_AerDust_dim', num_CM_AerDust, spcd_2(1))
!                           ----------------
      call NcDef_dimension (ncid_AerDust, 'chr_dim', 80, chrd(1))
!                                          ---------
      call NcDef_dimension(ncid_AerDust, rec_dim_name, NF_UNLIMITED, recd(1))
!                                         ------------
!     -----------------------------------------
!     Define variables and variable attributes.
!     -----------------------------------------

      var2(:) = (/ lond(1), latd(1) /)
      call NcDef_variable (ncid_AerDust, 'mcor', NF_FLOAT, 2, var2, varid)

      call NcDef_variable  &
     &        (ncid_AerDust, 'AerDust_dim', NF_FLOAT, 1, spcd, varid)
!                             -----------
      call NcDef_var_attributes (ncid_AerDust, varid, 'long_name', 'Aerosol Species index')
      call NcDef_var_attributes (ncid_AerDust, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_AerDust, varid, 'coord_labels', 'AerDust_labels')
      call NcDef_var_attributes (ncid_AerDust, varid, 'selection_category', 'NULL')

      call NcDef_variable  &
     &       (ncid_AerDust, 'CM_AerDust_dim', NF_FLOAT, 1, spcd_2, varid)
!                           ----------------
      call NcDef_var_attributes (ncid_AerDust, varid, 'long_name', 'Column Mass Species index')
      call NcDef_var_attributes (ncid_AerDust, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_AerDust, varid, 'coord_labels', 'CMAerDust_labels')
      call NcDef_var_attributes (ncid_AerDust, varid, 'selection_category', 'NULL')
!
      var2(:) = (/ hdfd(1), recd(1) /)
      call NcDef_variable (ncid_AerDust, hdr_var_name, NF_INT, 2, var2, varid)
!                                         ------------
      call NcDef_var_attributes (ncid_AerDust, varid, 'long_name', 'Header')
      call NcDef_var_attributes (ncid_AerDust, varid, 'units', 'gmi_sec, nymd, nhms')

      var2(:) = (/ chrd(1), spcd(1) /)
      call NcDef_variable  &
     &        (ncid_AerDust, 'AerDust_labels', NF_CHAR, 2, var2, varid)
!                            ---------------
      call NcDef_var_attributes (ncid_AerDust, varid, 'long_name', 'Aerosol/Dust name')
      call NcDef_var_attributes (ncid_AerDust, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_AerDust, varid, 'selection_category', 'NULL')

      call NcDef_variable (ncid_AerDust, 'am', NF_FLOAT, 1, prsd, varid)
!                                     ----
      call NcDef_var_attributes (ncid_AerDust, varid, 'long_name', 'Hybrid pressure term')
      call NcDef_var_attributes (ncid_AerDust, varid, 'units', 'unitless')

      call NcDef_variable (ncid_AerDust, 'bm', NF_FLOAT, 1, prsd, varid)
!                                     ----
      call NcDef_var_attributes (ncid_AerDust, varid, 'long_name', 'Hybrid sigma term')
      call NcDef_var_attributes (ncid_AerDust, varid, 'units', 'unitless')

      call NcDef_variable (ncid_AerDust, 'ai', NF_FLOAT, 1, prsp1d, varid)
!                                     ----
      call NcDef_var_attributes (ncid_AerDust, varid, 'long_name', 'Hybrid pressure edge term')
      call NcDef_var_attributes (ncid_AerDust, varid, 'units', 'unitless')

      call NcDef_variable (ncid_AerDust, 'bi', NF_FLOAT, 1, prsp1d, varid)
!                                     ----
      call NcDef_var_attributes (ncid_AerDust, varid, 'long_name', 'Hybrid sigma edge term')
      call NcDef_var_attributes (ncid_AerDust, varid, 'units', 'unitless')

      call NcDef_variable (ncid_AerDust, 'pt', NF_FLOAT, 0, prsp1d, varid)
!
! Monthly average aerosol/mineral dust data
!
      var5(:) = (/ lond(1), latd(1), prsd(1), spcd(1), recd(1) /)
      call NcDef_variable  &
     &         (ncid_AerDust, 'OptDepth', NF_FLOAT, 5, var5, varid)
!                            -----------
      var2(:) = (/ chrd(1), spcd_2(1) /)
      call NcDef_variable (ncid_AerDust, 'CMAerDust_labels', NF_CHAR, 2, var2, varid)
!                                          -----------------
      call NcDef_var_attributes (ncid_AerDust, varid, 'long_name', 'Column Mass label')
      call NcDef_var_attributes (ncid_AerDust, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_AerDust, varid, 'selection_category', 'NULL')
!
! Monthly average column mass for aerosol/dust
!
      var4(:) = (/ lond(1), latd(1), spcd_2(1), recd(1) /)
      call NcDef_variable  &
     &         (ncid_AerDust, 'ColumnMass_AerDust', NF_FLOAT, 4, var4, varid)
!                             --------------------
      call NcDef_var_attributes (ncid_AerDust, varid, 'long_name', 'Column Mass')
      call NcDef_var_attributes (ncid_AerDust, varid, 'units', 'g/m2')
!
! Monthly average grid_height
!
      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable  &
     &         (ncid_AerDust, 'Grid_Height', NF_FLOAT, 4, var4, varid)
!                              -----------
      call NcDef_var_attributes (ncid_AerDust, varid, 'long_name', 'Grid Height')
      call NcDef_var_attributes (ncid_AerDust, varid, 'units', 'm')

!    -----------------
      if (pr_psf) then
        var3(:) = (/ lond(1), latd(1), recd(1) /)
        call NcDef_variable (ncid_AerDust, 'psf', NF_FLOAT, 3, var3, varid)
!                                       -----
        call NcDef_var_attributes (ncid_AerDust, varid, 'long_name', 'Surface pressure')
        call NcDef_var_attributes (ncid_AerDust, varid, 'units', 'hPa')
      end if

!     -------------------------
!     Define global attributes.
!     -------------------------

      call NcDef_glob_attributes (ncid_AerDust, 'title',  &
     &       'Aerosol/dust optical depth and surface area diagnostic file')

      call NcSetFill (ncid_AerDust, NF_NOFILL, omode)

      call NcEnd_def (ncid_AerDust)

      return

      end subroutine Define_Netcdf_Out_AerDust
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Write_Netcdf_Hdr_AerDust
!
! DESCRIPTION
!   This routine creates some header information for the qj NetCDF output
!   file and writes it out.
!
!-----------------------------------------------------------------------------

      subroutine Write_Netcdf_Hdr_AerDust (gmiDomain, metFields)

      use m_ncGeneralOpsOutput, only : Write_Netcdf_Hdr_Gen

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      type(t_gmiDomain), intent(in) :: gmiDomain
      type(t_metFields), intent(in) :: metFields

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ic

      integer :: cnt1d (1), cnt2d (2)
      integer :: strt1d(1), strt2d(2)

      real*8  :: prsdat(k1:k2)
      real*8, allocatable  :: prsdat_r(:)
      real*8, allocatable  :: am_r(:)
      real*8, allocatable  :: bm_r(:)
      real*8, allocatable  :: ai_r(:)
      real*8, allocatable  :: bi_r(:)

      character(len=80)  :: spclab(num_AerDust)
      real*8             :: aerdat(num_AerDust)
!
      character(len=80)  :: spclab2(num_CM_AerDust)
      real*8             :: aerdat2(num_CM_AerDust)

      logical, save :: first = .true.

      real*8 , allocatable :: mcor(:,:)
      real*8 , allocatable :: londeg(:), latdeg(:)
      real*8 , allocatable :: ai(:), bi(:)
      real*8 , allocatable :: am(:), bm(:)
      real*8 :: pt
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) then
        Write (6,*) 'Write_Netcdf_Hdr_AerDust called by ', procID
      end if

      if (first) then
         call Set_AerDust_Labels()
         first = .false.
      endif

      allocate(latdeg(ju1_gl:j2_gl))
      call Get_latdeg(gmiDomain, latdeg)

      allocate(londeg(i1_gl:i2_gl))
      call Get_londeg(gmiDomain, londeg)

      allocate(mcor(i1_gl:i2_gl,ju1_gl:j2_gl))
      call Get_mcorGlob(gmiDomain, mcor)

!     =========================
      call Write_Netcdf_Hdr_Gen  &
!     =========================
     &  (ncid_AerDust, latdeg, londeg, pr_diag, procID, i1_gl, i2_gl, ju1_gl, j2_gl, &
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

      call Ncwr_1d (prsdat, ncid_AerDust, prs_dim_name, strt1d, cnt1d)

      call Ncwr_1d (am, ncid_AerDust, 'am', strt1d, cnt1d)
      call Ncwr_1d (bm, ncid_AerDust, 'bm', strt1d, cnt1d)

      cnt1d(1) = numVert+1

      call Ncwr_1d (ai, ncid_AerDust, 'ai', strt1d, cnt1d)
      call Ncwr_1d (bi, ncid_AerDust, 'bi', strt1d, cnt1d)

      call Ncwr_Scal (pt, ncid_AerDust, 'pt')

      strt2d(:) = (/ 1, 1 /)
      cnt2d (:) = (/ numLon, numLat /)

      call Ncwr_2d (mcor, ncid_AerDust, 'mcor', strt2d, cnt2d)

!     -------------------
!     aerosol/dust labels.
!     -------------------

      do ic = 1, num_AerDust
        aerdat(ic) = ic
      end do

      strt1d(1) = 1
      cnt1d (1) = num_AerDust

      call Ncwr_1d (aerdat, ncid_AerDust, 'AerDust_dim', strt1d, cnt1d)

      do ic = 1, num_AerDust
        spclab(ic) = AerDust_labels(ic)
      end do

      strt2d(:) = (/  1, 1 /)
      cnt2d (:) = (/ 80, num_AerDust /)

      call Ncwr_2d_Char(spclab, ncid_AerDust, 'AerDust_labels', strt2d, cnt2d)

!     -----------------------------------
!     aerosol/dust labels for column mass.
!     -----------------------------------

      do ic = 1, num_CM_AerDust
        aerdat2(ic) = ic
      end do

      strt1d(1) = 1
      cnt1d (1) = num_CM_AerDust

      call Ncwr_1d (aerdat2, ncid_AerDust, 'CM_AerDust_dim', strt1d, cnt1d)

      do ic = 1, num_CM_AerDust
        spclab2(ic) = CMAerDust_labels(ic)
      end do

      strt2d(:) = (/  1, 1 /)
      cnt2d (:) = (/ 80, num_CM_AerDust /)

      call Ncwr_2d_Char(spclab2, ncid_AerDust,'CMAerDust_labels', strt2d, cnt2d)

      return

      end subroutine Write_Netcdf_Hdr_AerDust

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Prep_Netcdf_AerDust
!
! DESCRIPTION
!   This routine prepares the aerosol/dust NetCDF output.
!
! ARGUMENTS
!   time_for_nc  : time for NetCDF output?
!
!-----------------------------------------------------------------------------

      subroutine Prep_Netcdf_AerDust (time_for_nc, Chemistry, metFields)
!
! !USES:
      use GmiChemistryMethod_mod, only : t_Chemistry
      use GmiChemistryMethod_mod, only : Get_dust, Get_wAersl, Get_dAersl
      use GmiChemistryMethod_mod, only : Get_optDepth, Set_optDepth

      implicit none


!     ----------------------
!     Argument declarations.
!     ----------------------

      logical :: time_for_nc
      type(t_Chemistry) :: Chemistry
      type(t_metFields) :: metFields

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer, save :: counter = 0

      real*8  :: rcounter

      integer :: ii, l
      real*8  :: dust  (i1:i2,ju1:j2,k1:k2,NSADdust)
      real*8  :: wAersl(i1:i2,ju1:j2,k1:k2,NSADaer)
      real*8  :: dAersl(i1:i2,ju1:j2,k1:k2,2      )
      real*8, allocatable :: optDepth0(:,:,:,:)
      real*8, allocatable :: gridBoxHeight(:,:,:)
      real*8, allocatable :: pctm2(:,:)
      real*8, allocatable :: opticalDepth(:,:,:,:)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Prep_Netcdf_AerDust called by ', procID
      end if

      if (counter == 0) then
         CM_AerDust_nc (:,:,:) = 0.0d0
         grid_height_nc(:,:,:) = 0.0d0
      end if

     allocate(gridBoxHeight(i1:i2,ju1:j2,k1:k2))
     call Get_gridBoxHeight(metFields, gridBoxHeight)

     allocate(pctm2(ilo:ihi, julo:jhi))
     call Get_pctm2(metFields, pctm2)

!     ============
      if (do_mean) then
!     ============

        if (counter == 0) then
           if (pr_psf) psf_AerDust_mean(:,:) = 0.0d0
        end if

        call Get_dust   (Chemistry, dust)
        call Get_wAersl (Chemistry, wAersl)
        call Get_dAersl (Chemistry, dAersl)

        counter = counter + 1

        if (pr_psf) &
     &     psf_AerDust_mean(i1:i2,ju1:j2) = psf_AerDust_mean(i1:i2,ju1:j2) + &
     &                                      pctm2(i1:i2,ju1:j2)

        grid_height_nc(:,:,:) = grid_height_nc(:,:,:) + gridBoxHeight(:,:,kx1:kx2)

        ! Dust           Column Mass (g/m2)
        do ii=1,NSADdust
           do l = k1, k2
              CM_AerDust_nc (:,:,1) = CM_AerDust_nc (:,:,1) + &
     &            DUST(:,:,l,ii) * gridBoxHeight(:,:,l) *1000.d0
           end do
        end do

#ifndef nonZeroInd_tracers
        do l = k1, k2
           ! Sulfate        Column Mass (g/m2)
           CM_AerDust_nc (:,:,2) = CM_AerDust_nc (:,:,2) + &
     &        WAERSL(:,:,l,1)  * gridBoxHeight(:,:,l) *1000.d0

           ! Black Carbon   Column Mass (g/m2)
           CM_AerDust_nc (:,:,3) = CM_AerDust_nc (:,:,3) + &
     &       (WAERSL(:,:,l,2) + DAERSL(:,:,l,1)) * &
     &        gridBoxHeight(:,:,l) *1000.d0

           ! Organic Carbon Column Mass (g/m2)
           CM_AerDust_nc (:,:,4) = CM_AerDust_nc (:,:,4) + &
     &       (WAERSL(:,:,l,3) + DAERSL(:,:,l,2)) * &
     &        gridBoxHeight(:,:,l) *1000.d0

           ! Sea Salt       Column Mass (g/m2)
           CM_AerDust_nc (:,:,5) = CM_AerDust_nc (:,:,5) + &
     &       (WAERSL(:,:,l,4) + WAERSL(:,:,l,5) + WAERSL(:,:,l,6) + WAERSL(:,:,l,7)) * &
     &       gridBoxHeight(:,:,l) *1000.d0
        end do
#endif

        if (time_for_nc) then

!         ----------------
!         Calculate means.
!         ----------------

          rcounter = counter

          if (pr_psf)  psf_AerDust_mean(:,:) = psf_AerDust_mean(:,:) / rcounter

          CM_AerDust_nc (:,:,:) = CM_AerDust_nc (:,:,:) / rcounter

          grid_height_nc(:,:,:) = grid_height_nc(:,:,:) / rcounter

          Allocate (opticalDepth(i1:i2, ju1:j2, k1:k2, num_AerDust))
          call Get_optDepth (Chemistry, opticalDepth)

          optDepth(:,:,:,:) = opticalDepth(:,:,kx1:kx2,:) / rcounter

          allocate(optDepth0(i1:i2,ju1:j2,k1:k2,num_AerDust))
          OptDepth0(:,:,:,:) = 0.0d0
          call Set_optDepth(Chemistry, optDepth0)

          deallocate(optDepth0)
          deallocate(opticalDepth)

          counter = 0

          if (pr_psf) psf_AerDust_nc(:,:) = psf_AerDust_mean(:,:)

        end if

      else

        if (time_for_nc) then

!         --------------------------------------------------
!         Fill output arrays with current "non-mean" values.
!         --------------------------------------------------

          if (pr_psf) psf_AerDust_nc(i1:i2,ju1:j2) = pctm2(i1:i2,ju1:j2)
         end if
   

      end if

     deallocate(gridBoxHeight)

      return

      end subroutine Prep_Netcdf_AerDust

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Write_Netcdf_AerDust
!
! DESCRIPTION
!   This routine writes the aerosol/dust variable data NetCDF output.
!
!
! ARGUMENTS
!   None
!
!-----------------------------------------------------------------------------

      subroutine Write_Netcdf_AerDust (nymd, nhms, gmi_sec)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: nymd, nhms
      real*8 , intent(in) :: gmi_sec
!
! !DEFINED PARAMETERS:
      integer, parameter :: PSF_MsgNum        = 100
      integer, parameter :: AerDust_MsgNum    = 101
      integer, parameter :: GridHeight_MsgNum = 102
!
! !LOCAL VARIABLES:
      integer :: cnt2d (2), cnt3d (3), cnt4d (4)
      integer :: strt2d(2), strt3d(3), strt4d(4)
      integer :: hdr(NETCDF_HDF)
      integer         :: klen
      real*8, allocatable :: grid_heightGlob(:,:,:)
      real*8, allocatable :: psf_AerDustGlob(:,:)
      real*8, allocatable :: CM_AerDustGlob(:,:,:)

!EOP
!------------------------------------------------------------------------------
!BOC
!
      if (pr_diag) Write (6,*) 'Write_Netcdf_AerDust called by ', procID

      if (iAmRootProc) then

         strt2d(:) = (/ 1, rnum_out_AerDust /)

         cnt2d (:) = (/ NETCDF_HDF, 1 /)

         hdr(1) = Nint (gmi_sec)
         hdr(2) = nymd
         hdr(3) = nhms

         call Ncwr_2d_Int (hdr, ncid_AerDust, hdr_var_name, strt2d, cnt2d)

         if (pr_psf) allocate(psf_AerDustGlob(i1_gl:i2_gl,ju1_gl:j2_gl))
         allocate(CM_AerDustGlob(i1_gl:i2_gl,ju1_gl:j2_gl,1:num_CM_AerDust))
         allocate(grid_heightGlob(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))
      end if

      call subDomain2Global (CM_AerDustGlob, CM_AerDust_nc, i1_gl, i2_gl, &
     &              ju1_gl, j2_gl, i1, i2, ju1, j2, 1, num_CM_AerDust,    &
     &              rootProc, procID, map1_u, numDomains, AerDust_MsgNum, &
     &              commuWorld)

      if (pr_psf) then
          call subDomain2Global (psf_AerDustGlob, psf_AerDust_nc, i1_gl, i2_gl,&
     &              ju1_gl, j2_gl, i1, i2, ju1, j2, rootProc, procID,   &
     &              map1_u, numDomains, PSF_MsgNum, commuWorld)
      end if

      call subDomain2Global (grid_heightGlob, grid_height_nc, i1_gl, i2_gl, &
     &              ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, rootProc,     &
     &              procID, map1_u, numDomains, AerDust_MsgNum, commuWorld)

      if (iAmRootProc) then
         strt4d(:) = (/ 1, 1, 1, rnum_out_AerDust /)
         cnt4d (:) = (/ numLon, numLat, 5, 1 /)

         call Ncwr_4d (CM_AerDustGlob, ncid_AerDust, 'ColumnMass_AerDust',  &
     &                    strt4d, cnt4d)

         cnt4d (:) = (/ numLon, numLat, numVert, 1 /)

         call Ncwr_4d (grid_heightGlob, ncid_AerDust, 'Grid_Height',  &
     &                    strt4d, cnt4d)

         strt3d(:) = (/ 1, 1, rnum_out_AerDust /)
         cnt3d (:) = (/ numLon, numLat, 1 /)

         if (pr_psf)  &
     &      call Ncwr_3d (psf_AerDustGlob, ncid_AerDust, 'psf', strt3d, cnt3d)

         if (pr_psf) deallocate(psf_AerDustGlob)
         deallocate(CM_AerDustGlob)
         deallocate(grid_heightGlob)

         call Ncdo_Sync (ncid_AerDust)

      end if

      return

      end subroutine Write_Netcdf_AerDust

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Set_AerDust_Labels
!
! !INTERFACE:
!
      subroutine Set_AerDust_Labels()
!
! !USES:
!
      implicit none
!
! !DESCRIPTION: This routine sets aerosol/dust labels for output netCDF file.
!
! !REVISION HISTORY:
!   February2005, Jules Kouatchou (Jules.Kouatchou.1@gsfc.nasa.gov)
!     Original code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
! Labels for aerosol and mineral dust
!
      AerDust_labels( 1) = 'Cloud optical depths (1000 nm)'
      AerDust_labels( 2) = 'Max Overlap Cld Frac'
      AerDust_labels( 3) = 'Random Overlap Cld Frac'
      AerDust_labels( 4) = 'Dust optical depths (400 nm)'
      AerDust_labels( 5) = 'Dust surface areas'
      AerDust_labels( 6) = 'Sulfate Optical Depth (400 nm)'
      AerDust_labels( 7) = 'Hygroscopic growth of SO4'
      AerDust_labels( 8) = 'Sulfate Surface Area'
      AerDust_labels( 9) = 'Black Carbon Optical Depth (400 nm)'
      AerDust_labels(10) = 'Hygroscopic growth of Black Carbon'
      AerDust_labels(11) = 'Black Carbon Surface Area'
      AerDust_labels(12) = 'Organic Carbon Optical Depth (400 nm)'
      AerDust_labels(13) = 'Hygroscopic growth of Organic Carbon'
      AerDust_labels(14) = 'Organic Carbon Surface Area'
      AerDust_labels(15) = 'Sea Salt (accum) Opt Depth (400 nm)'
      AerDust_labels(16) = 'Hygroscopic growth of Sea Salt (accum)'
      AerDust_labels(17) = 'Sea Salt (accum) Surface Area'
      AerDust_labels(18) = 'Sea Salt (coarse1) Opt Depth(400 nm)'
      AerDust_labels(19) = 'Hygroscopic growth of Sea Salt (coarse1)'
      AerDust_labels(20) = 'Sea Salt (coarse1) Surface Area'
      AerDust_labels(21) = 'Sea Salt (coarse2) Opt Depth(400 nm)'
      AerDust_labels(22) = 'Hygroscopic growth of Sea Salt (coarse2)'
      AerDust_labels(23) = 'Sea Salt (coarse2) Surface Area'
      AerDust_labels(24) = 'Sea Salt (coarse3) Opt Depth(400 nm)'
      AerDust_labels(25) = 'Hygroscopic growth of Sea Salt (coarse3)'
      AerDust_labels(26) = 'Sea Salt (coarse3) Surface Area'
!
! Labels for column mass
!
      CMAerDust_labels(1) = 'Dust Column Mass (g/m2)'
      CMAerDust_labels(2) = 'Sulfate Column Mass (g/m2)'
      CMAerDust_labels(3) = 'Black Carbon Column Mass (g/m2)'
      CMAerDust_labels(4) = 'Organic Carbon Column Mass (g/m2)'
      CMAerDust_labels(5) = 'Sea Salt Column Mass (g/m2)'

      return

      end subroutine Set_AerDust_Labels
!EOC
!-------------------------------------------------------------------------
!
! DESCRIPTION
!   This routine buffers the aerosol/dust data to the output file.
!   The Slaves send a slice of the data to the Master, the Master writes it out,
!   they then send another slice and it is written out, and so on.
!
! ARGUMENTS
!   None
!
!-----------------------------------------------------------------------------

      subroutine bufferOutput_AerDust ( )
!
! !USES:
      use GmiBuffer_mod, only : bufferOutput4d_Nomean
!
      implicit none
!
! !DEFINED PARAMETERS:
      integer, parameter :: OptDepth_MsgNum = 1005
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
      call bufferOutput4d_Nomean (AerDust_var_name, kx1, kx2, commuWorld, &
     &           OptDepth_MsgNum, ncid_AerDust, num_AerDust, rnum_out_AerDust, &
     &           optDepth, OptDepth_nc, map1_u, numDomains, rootProc, procID,  &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)

      return

      end subroutine bufferOutput_AerDust
!EOC
!-----------------------------------------------------------------------------

end module GmiControlAerosolDust_mod
