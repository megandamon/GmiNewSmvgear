!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlSAD_mod
!
! !INTERFACE:
!
module GmiControlSAD_mod
!
! !USES:
      use  GmiSub2Glob_mod      , only : subDomain2Global
      use  GmiMessagePassing_mod, only : synchronizeGroup
      use GmiFileOperations_mod , only : makeOutfileName
      use m_netcdf_io_close     , only : Nccl, Nccl_Noerr
      use m_netcdf_io_write, only : Ncwr_2d_Int, Ncwr_5d, Ncwr_3d, Ncwr_2d, &
     &        Ncwr_4d
      use m_netcdf_io_write     , only : Ncwr_1d, Ncwr_Scal, Ncwr_2d_Char
      use m_netcdf_io_create    , only : Ncdo_Sync, Nccr_Wr
      use m_netcdf_io_define    , only : NcDef_glob_attributes, NcDef_dimension
      use m_netcdf_io_define    , only : NcDef_var_attributes, NcDef_variable
      use m_netcdf_io_define    , only : NcSetFill, NcEnd_def
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

      use GmiChemistryMethod_mod  , only : t_Chemistry, Get_num_sad, Get_sad_opt
      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Get_hdf_dim_name,    &
     &       Get_lon_dim_name, Get_lat_dim_name, Get_sadOutputFrequency, Get_k1r_gl, &
     &       Get_hdr_var_name, Get_rec_dim_name, Get_prs_dim_name, Get_k2r_gl, &
     &       Get_sad_var_name, Get_sad_dim_name, Get_do_mean, Get_pr_sad,      &
     &       Get_pr_diag, Get_problem_name

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
      public  :: controlOutputSAD, initializeOutputSAD
      public  :: finalizeOutputSAD

#     include "gmi_phys_constants.h"
#     include "GmiParameters.h"
!
! !DESCRIPTION:
! Contains the necessary routines for the SAD outputs. 
! Three routines are visible from the outside:
! %
! \begin{description}
! \item[initializeOutputSAD:] called once in the initialization steps to
!      open the netCDF file, define variables in the file, write out header
!      information and allocate variables.
! \item[controlOutputSAD:] called at the end of each model time step to
!      update variables and if needed for communications and writing out
!      data in the file.
! \item[finalizeOutputSAD:] called at the end of the integration to
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
      real*8, pointer, save :: sad_nc     (:,:,:) => null() ! surface area densities (cm^2/cm^3)
      real*8, pointer, save :: h2oback_nc (:,:,:) => null() ! background      h2o  (mixing ratio)
      real*8, pointer, save :: h2ocond_nc (:,:,:) => null() ! condensed phase h2o  (mixing ratio)
      real*8, pointer, save :: hno3cond_nc(:,:,:) => null() ! condensed phase hno3 (mixing ratio)
      real*8, pointer, save :: hno3gas_nc (:,:,:) => null() ! gas       phase hno3 (mixing ratio)
      real*8, pointer, save :: reffice_nc (:,:,:) => null() ! effective radius of ICE aerosols  (cm)
      real*8, pointer, save :: reffsts_nc (:,:,:) => null() ! effective radius of STS aerosols  (cm)
      real*8, pointer, save :: vfall_nc   (:,:,:) => null() ! effective aerosol fall velocities (cm/s)

      ! Variables manipulate by worker processors only         
      real*8, pointer, save :: sad_mean   (:,:,:,:) => null()
      real*8, pointer, save :: h2oback_mean (:,:,:) => null()
      real*8, pointer, save :: h2ocond_mean (:,:,:) => null()
      real*8, pointer, save :: hno3cond_mean(:,:,:) => null()
      real*8, pointer, save :: hno3gas_mean (:,:,:) => null()
      real*8, pointer, save :: reffice_mean (:,:,:) => null()
      real*8, pointer, save :: reffsts_mean (:,:,:) => null()
      real*8, pointer, save :: vfall_mean   (:,:,:) => null()
!
      ! netCDF file identifier
      integer,          save :: ncid_sad

      ! Counter for the number of records
      integer,          save :: rnum_out_sad

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

      integer,          save :: numSAD
      integer,          save :: sad_opt
      
      logical,          save :: iAmRootProc
      logical,          save :: pr_diag
      integer,          save :: procID, rootProc
      integer,          save :: numDomains
      integer,          save :: commuWorld
      integer, pointer, save :: map1_u       (:,:,:) => null()
      character (len=MAX_LENGTH_VAR_NAME), save :: sad_dim_name, rec_dim_name, prs_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: hdf_dim_name, lon_dim_name, lat_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: sad_var_name, hdr_var_name
      logical           , save :: pr_sad, do_mean   
      real*8            , save :: sadOutputFrequency
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
! !IROUTINE: initializeOutputSAD
!
! !INTERFACE:
!
      subroutine initializeOutputSAD(Diagnostics, gmiGrid, gmiDomain, &
     &                     Chemistry, metFields)
!
! !USES:
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_gmiGrid  ), intent(in) :: gmiGrid
      type(t_gmiDomain), intent(in) :: gmiDomain
      type(t_Diagnostics), intent(in) :: Diagnostics
      type(t_Chemistry), intent(in) :: Chemistry
      type(t_metFields), intent(in) :: metFields
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
!EOP
!-----------------------------------------------------------------------------
!BOC

      !#####################################################
      ! Initialization of variables used in all the routines
      !#####################################################

      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write (6,*) 'initializeOutputSAD called by ', procID

      call Get_num_sad(Chemistry, numSAD)

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
      call Get_sadOutputFrequency(Diagnostics, sadOutputFrequency)
      call Get_pr_sad      (Diagnostics, pr_sad)
      call Get_hdr_var_name(Diagnostics, hdr_var_name)
      call Get_rec_dim_name(Diagnostics, rec_dim_name)
      call Get_prs_dim_name(Diagnostics, prs_dim_name)
      call Get_sad_var_name(Diagnostics, sad_var_name)
      call Get_sad_dim_name(Diagnostics, sad_dim_name)
      call Get_hdf_dim_name (Diagnostics, hdf_dim_name)
      call Get_lat_dim_name (Diagnostics, lat_dim_name)
      call Get_lon_dim_name (Diagnostics, lon_dim_name)

      call Get_mdt (metFields, mdt)

      call Get_sad_opt (Chemistry, sad_opt)

      if (pr_sad) then

         if (iAmRootProc) then

            !###########################
            ! Initialize the output file
            !###########################

            ! Determine the file name
            call Get_problem_name(Diagnostics, problem_name)
            call makeOutfileName (fname, '.sad.nc', problem_name)

            ! Create the file and assign a file identifier
            call Nccr_Wr (ncid_sad, fname)

            ! Define the variables in the file
            call Define_Netcdf_Out_SAD(Diagnostics)

            ! Write header data
            call Write_Netcdf_Hdr_SAD (gmiDomain, Diagnostics, metFields)

            call Ncdo_Sync (ncid_sad)

            ! Initialize the counter for the number of records
            rnum_out_sad = 1
         end if

         !########################
         ! Allocation of variables
         !########################
         call allocateVariablesSAD()
      end if

      return

      end subroutine initializeOutputSAD
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: controlOutSAD
!
! !INTERFACE:
!
      subroutine controlOutputSAD (last_tstp, Chemistry, gmiClock)
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
      if (pr_diag) Write (6,*) 'controlOutputSAD called by ', procID

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
     &      month_save, month, day, nhms, ndt, gmi_sec, sadOutputFrequency)

      month_save = month

      if (last_tstp .and. (.not. time_for_nc)) then
        if (Mod (num_time_steps, Nint (mdt) / ndt) == 0) then

          time_for_nc = .true.

        end if
      end if

      if (time_for_nc .or. do_mean) then
         !====================
          call Prep_Netcdf_SAD  (time_for_nc, Chemistry)
         !====================
      end if


!     ================
      if (time_for_nc) then
!     ================

!        =====================
         call Write_Netcdf_SAD (nymd, nhms, gmi_sec)
!        =====================

!       ======================
        call bufferOutput_SAD (Chemistry)
!       ======================

        if (iAmRootProc) then
            rnum_out_sad = rnum_out_sad + 1
        end if

!     ======
      end if
!     ======

      return

      end subroutine controlOutputSAD
!EOC
!-----------------------------------------------------------------------------
!BOP
!     
! !IROUTINES: finalizeOutputSAD
!
! !INTERFACE:
!
      subroutine finalizeOutputSAD()
!
      implicit none
! 
! !DESCRIPTION:
! Deallocates variables necessary to produce qj outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'finalizeOutputSAD called by ', procID

      if (pr_sad) then
     
         deallocate (sad_nc)

         if ((sad_opt == 1) .or. (sad_opt == 2)) then
            deallocate (hno3cond_nc)
            deallocate (hno3gas_nc )

            if (sad_opt == 2) then
               deallocate (h2oback_nc)
               deallocate (h2ocond_nc)
               deallocate (reffice_nc)
               deallocate (reffsts_nc)
               deallocate (vfall_nc  )
            end if
         end if

         if (do_mean) then
            deallocate (sad_mean)

            if ((sad_opt == 1) .or. (sad_opt == 2)) then
               deallocate (hno3cond_mean)
               deallocate (hno3gas_mean )

               if (sad_opt == 2) then
                  deallocate (h2oback_mean)
                  deallocate (h2ocond_mean)
                  deallocate (reffice_mean)
                  deallocate (reffsts_mean)
                  deallocate (vfall_mean  )
               end if
            end if
         end if

         if (iAmRootProc) call Nccl_Noerr (ncid_sad)
      end if

      return

      end subroutine finalizeOutputSAD
!EOC
!-----------------------------------------------------------------------------
!BOP
!     
! !IROUTINES: allocateVariablesSAD
!
! !INTERFACE:
!
      subroutine allocateVariablesSAD()
!
      implicit none
! 
! !DESCRIPTION:
! Allocates variables necessary to produce qj outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'allocateVariablesSAD called by ', procID

      if (pr_sad) then

         Allocate (sad_nc(i1:i2, ju1:j2, kx1:kx2))
         sad_nc = 0.0d0

         if ((sad_opt == 1) .or. (sad_opt == 2)) then
            Allocate (hno3cond_nc(i1:i2, ju1:j2, kx1:kx2))
            Allocate (hno3gas_nc (i1:i2, ju1:j2, kx1:kx2))
            hno3cond_nc = 0.0d0; hno3gas_nc = 0.0d0

            if (sad_opt == 2) then
               Allocate (h2oback_nc(i1:i2, ju1:j2, kx1:kx2))
               Allocate (h2ocond_nc(i1:i2, ju1:j2, kx1:kx2))
               Allocate (reffice_nc(i1:i2, ju1:j2, kx1:kx2))
               Allocate (reffsts_nc(i1:i2, ju1:j2, kx1:kx2))
               Allocate (vfall_nc  (i1:i2, ju1:j2, kx1:kx2))
               h2oback_nc = 0.0d0; h2ocond_nc = 0.0d0
               reffice_nc = 0.0d0; reffsts_nc = 0.0d0
               vfall_nc   = 0.0d0
            end if
         end if

         if (do_mean) then
            Allocate (sad_mean(i1:i2, ju1:j2, kx1:kx2, numSAD))
            sad_mean = 0.0d0

            if ((sad_opt == 1) .or. (sad_opt == 2)) then
               Allocate (hno3cond_mean(i1:i2, ju1:j2, kx1:kx2))
               Allocate (hno3gas_mean (i1:i2, ju1:j2, kx1:kx2))
               hno3cond_mean = 0.0d0; hno3gas_mean = 0.0d0
  
               if (sad_opt == 2) then
                  Allocate (h2oback_mean(i1:i2, ju1:j2, kx1:kx2))
                  Allocate (h2ocond_mean(i1:i2, ju1:j2, kx1:kx2))
                  Allocate (reffice_mean(i1:i2, ju1:j2, kx1:kx2))
                  Allocate (reffsts_mean(i1:i2, ju1:j2, kx1:kx2))
                  Allocate (vfall_mean  (i1:i2, ju1:j2, kx1:kx2))
                  h2oback_mean = 0.0d0; h2ocond_mean = 0.0d0
                  reffice_mean = 0.0d0; reffsts_mean = 0.0d0
                  vfall_mean   = 0.0d0
                end if
            end if
         end if

      end if

      return

      end subroutine allocateVariablesSAD
!EOC
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Define_Netcdf_Out_SAD
!
! DESCRIPTION
!   This routine makes the necessary definitions for the qj NetCDF output
!   file.

      subroutine Define_Netcdf_Out_SAD (Diagnostics)

      use m_ncGeneralOpsOutput, only : Define_Netcdf_Out_Gen

      implicit none

#     include "netcdf.inc"

!
! !INPUT PARAMETERS:
      type(t_Diagnostics), intent(in) :: Diagnostics

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ierr
      integer :: omode
      integer :: varid

      integer :: chrd(1)
      integer :: hdfd(1)
      integer :: lond(1), latd(1), prsd(1)
      integer :: ydasd(1), recd (1), sadd(1)

      integer :: var2(2), var4(4), var5(5)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Define_Netcdf_Out_SAD called by ', procID
      end if

!     ==========================
      call Define_Netcdf_Out_Gen  &
!     ==========================
     &  (ncid_sad, nhdf, numLon, numLat, numVert, hdfd, lond, latd, prsd,  &
     &   hdf_dim_name, lon_dim_name, lat_dim_name, prs_dim_name)

!     ------------------
!     Define dimensions.
!     ------------------

      call NcDef_dimension(ncid_sad, sad_dim_name, numSAD, sadd(1))
!                                     ------------
      call NcDef_dimension(ncid_sad, rec_dim_name, NF_UNLIMITED, recd(1))
!                                     ------------

!     -----------------------------------------
!     Define variables and variable attributes.
!     -----------------------------------------

      call NcDef_variable (ncid_sad, sad_dim_name, NF_FLOAT, 1, sadd, varid)
!                                     ------------
      call NcDef_var_attributes (ncid_sad, varid, 'long_name', 'Sad')

      var2(:) = (/ hdfd(1), recd(1) /)
      call NcDef_variable (ncid_sad, hdr_var_name, NF_INT, 2, var2, varid)
!                                    ------------
      var5(:) = (/ lond(1), latd(1), prsd(1), sadd(1), recd(1) /)
      call NcDef_variable (ncid_sad, sad_var_name, NF_FLOAT, 5, var5, varid)
!                                    ------------
      if ((sad_opt == 1) .or. (sad_opt == 2)) then

        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)

        call NcDef_variable (ncid_sad, 'hno3cond', NF_FLOAT, 4, var4, varid)
!                                       ----------
        call NcDef_variable (ncid_sad, 'hno3gas', NF_FLOAT, 4, var4, varid)
!                                       ---------
        if (sad_opt == 2) then

          call NcDef_variable (ncid_sad, 'h2oback', NF_FLOAT, 4, var4, varid)
!                                         ---------
          call NcDef_variable (ncid_sad, 'h2ocond', NF_FLOAT, 4, var4, varid)
!                                         ---------
          call NcDef_variable (ncid_sad, 'reffsts', NF_FLOAT, 4, var4, varid)
!                                         ---------
          call NcDef_variable (ncid_sad, 'reffice', NF_FLOAT, 4, var4, varid)
!                                         ---------
          call NcDef_variable (ncid_sad, 'vfall', NF_FLOAT, 4, var4, varid)

        end if

      end if

      call NcSetFill (ncid_sad, NF_NOFILL, omode)

      call NcEnd_def (ncid_sad)


      return

      end subroutine Define_Netcdf_Out_SAD
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Hdr_SAD
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Hdr_SAD  (gmiDomain, Diagnostics, metFields)
!
! !USES:
      use m_ncGeneralOpsOutput, only : Write_Netcdf_Hdr_Gen
!
      implicit none
!
! !INPUT PARAMETERS:
      type (t_Diagnostics), intent(in) :: Diagnostics
      type(t_gmiDomain), intent(in) :: gmiDomain
      type(t_metFields), intent(in) :: metFields
!
! !DESCRIPTION:
! Creates some header information for the SAD netCDF output
! file and writes it out.
!
! !LOCAL VARIABLES:
      integer :: ic
      integer :: cnt1d (1), cnt2d (2)
      integer :: strt1d(1), strt2d(2)
      real*8  :: prsdat(k1:k2)
      real*8  :: saddat(numSAD)
      real*8 , allocatable :: londeg(:), latdeg(:)
      real*8 , allocatable :: am(:), bm(:)
      real*8 :: pt
!
!EOP
!-----------------------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Hdr_SAD called by ', procID

      allocate(latdeg(ju1_gl:j2_gl))
      call Get_latdeg(gmiDomain, latdeg)

      allocate(londeg(i1_gl:i2_gl))
      call Get_londeg(gmiDomain, londeg)

!     =========================
      call Write_Netcdf_Hdr_Gen  &
!     =========================
     &  (ncid_sad, latdeg, londeg, pr_diag, procID, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   hdf_dim_name, lat_dim_name, lon_dim_name)

      allocate(am(k1:k2))
      allocate(bm(k1:k2))

      call Get_pt(metFields, pt)
      call Get_am(metFields, am)
      call Get_bm(metFields, bm)

      prsdat(:) = (am(:) * pt)  + (bm(:) * 1000.0d0)

      strt1d(1) = 1
      cnt1d (1) = numVert

      call Ncwr_1d (prsdat(kx1:kx2), ncid_sad, prs_dim_name, strt1d, cnt1d)

      do ic = 1, numSAD
        saddat(ic) = ic
      end do

      strt1d(1) = 1
      cnt1d (1) = numSAD

      call Ncwr_1d (saddat, ncid_sad, sad_dim_name, strt1d, cnt1d)

      return

      return

      end subroutine Write_Netcdf_Hdr_SAD
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Prep_Netcdf_SAD
!
! !INTERFACE:
!
      subroutine Prep_Netcdf_SAD (time_for_nc, Chemistry)
!
! !USES:
      use GmiChemistryMethod_mod, only : t_Chemistry
      use GmiChemistryMethod_mod, only : Get_sadgmi  , Get_vfall
      use GmiChemistryMethod_mod, only : Get_hno3cond, Get_hno3gas
      use GmiChemistryMethod_mod, only : Get_h2ocond , Get_h2oback
      use GmiChemistryMethod_mod, only : Get_reffice , Get_reffsts
!
      implicit none
!
! !INPUT PARAMETERS:
      logical          , intent(in) :: time_for_nc ! time for netCDF output?
      type(t_Chemistry), intent(in) :: Chemistry
!
! !DESCRIPTION:
! Prepares the SAD netCDF output.
!
! !LOCAL VARIABLES:
      integer, save :: counter = 0
      real*8  :: rcounter
      real*8 , allocatable :: sadgmi0  (:,:,:,:), h2ocond (:,:,:)
      real*8 , allocatable :: vfall    (:,:,:)  , h2oback (:,:,:)
      real*8 , allocatable :: hno3cond (:,:,:)  , hno3gas (:,:,:)
      real*8 , allocatable :: reffice  (:,:,:)  , reffsts (:,:,:)
!
!EOP
!-----------------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Prep_Netcdf_SAD called by ', procID

      if ((sad_opt == 1) .or. (sad_opt == 2)) then
         Allocate (hno3cond(i1:i2, ju1:j2, k1:k2))
         call Get_hno3cond(Chemistry, hno3cond)

         Allocate (hno3gas (i1:i2, ju1:j2, k1:k2))
         call Get_hno3gas (Chemistry, hno3gas )

         if (sad_opt == 2) then
            Allocate (h2oback(i1:i2, ju1:j2, k1:k2))
            call Get_h2oback(Chemistry, h2oback)

            Allocate (h2ocond(i1:i2, ju1:j2, k1:k2))
            call Get_h2ocond(Chemistry, h2ocond)

            Allocate (reffice(i1:i2, ju1:j2, k1:k2))
            call Get_reffice(Chemistry, reffice)

            Allocate (reffsts(i1:i2, ju1:j2, k1:k2))
            call Get_reffsts(Chemistry, reffsts)

            Allocate (vfall  (i1:i2, ju1:j2, k1:k2))
            call Get_vfall  (Chemistry, vfall  )
         end if
      end if

!     ============
      if (do_mean) then
!     ============

        if (counter == 0) then

!         ----------------
!         Reset mean sums.
!         ----------------

          sad_mean(:,:,:,:) = 0.0d0

          if ((sad_opt == 1) .or. (sad_opt == 2)) then
            hno3cond_mean(:,:,:) = 0.0d0
            hno3gas_mean (:,:,:) = 0.0d0

            if (sad_opt == 2) then
              h2oback_mean(:,:,:) = 0.0d0
              h2ocond_mean(:,:,:) = 0.0d0
              reffice_mean(:,:,:) = 0.0d0
              reffsts_mean(:,:,:) = 0.0d0
              vfall_mean  (:,:,:) = 0.0d0
            end if
          end if

        end if

!       -----------------
!       Update mean sums.
!       -----------------

        allocate(sadgmi0(i1:i2, ju1:j2, k1:k2, numSad))
        call Get_sadgmi(Chemistry, sadgmi0)

        sad_mean(:,:,:,:) = sad_mean(:,:,:,:) + sadgmi0(:,:,kx1:kx2,:)

        if ((sad_opt == 1) .or. (sad_opt == 2)) then
          hno3cond_mean(:,:,:) = hno3cond_mean(:,:,:) + hno3cond(:,:,kx1:kx2)
          hno3gas_mean (:,:,:) = hno3gas_mean (:,:,:) + hno3gas (:,:,kx1:kx2)

          if (sad_opt == 2) then
            h2oback_mean(:,:,:) = h2oback_mean(:,:,:) + h2oback(:,:,kx1:kx2)
            h2ocond_mean(:,:,:) = h2ocond_mean(:,:,:) + h2ocond(:,:,kx1:kx2)
            reffice_mean(:,:,:) = reffice_mean(:,:,:) + reffice(:,:,kx1:kx2)
            reffsts_mean(:,:,:) = reffsts_mean(:,:,:) + reffsts(:,:,kx1:kx2)
            vfall_mean  (:,:,:) = vfall_mean  (:,:,:) + vfall  (:,:,kx1:kx2)
          end if
        end if


        counter = counter + 1


        if (time_for_nc) then

!         ----------------
!         Calculate means.
!         ----------------

          rcounter = counter

          sad_mean(:,:,:,:) = sad_mean(:,:,:,:) / rcounter

          if ((sad_opt == 1) .or. (sad_opt == 2)) then
            hno3cond_mean(:,:,:) = hno3cond_mean(:,:,:) / rcounter
            hno3gas_mean (:,:,:) = hno3gas_mean (:,:,:) / rcounter

            if (sad_opt == 2) then
              h2oback_mean(:,:,:) = h2oback_mean(:,:,:) / rcounter
              h2ocond_mean(:,:,:) = h2ocond_mean(:,:,:) / rcounter
              reffice_mean(:,:,:) = reffice_mean(:,:,:) / rcounter
              reffsts_mean(:,:,:) = reffsts_mean(:,:,:) / rcounter
              vfall_mean  (:,:,:) = vfall_mean  (:,:,:) / rcounter
            end if
          end if

          counter = 0

!         -----------------------------------------------
!         Fill output arrays with calculated mean values.
!         -----------------------------------------------

          if ((sad_opt == 1) .or. (sad_opt == 2)) then
            hno3cond_nc(:,:,:) = hno3cond_mean(:,:,:)
            hno3gas_nc (:,:,:) = hno3gas_mean (:,:,:)

            if (sad_opt == 2) then
              h2oback_nc(:,:,:) = h2oback_mean(:,:,:)
              h2ocond_nc(:,:,:) = h2ocond_mean(:,:,:)
              reffice_nc(:,:,:) = reffice_mean(:,:,:)
              reffsts_nc(:,:,:) = reffsts_mean(:,:,:)
              vfall_nc  (:,:,:) = vfall_mean  (:,:,:)
            end if
          end if

        end if

!     ====
      else
!     ====

        if (time_for_nc) then

!         --------------------------------------------------
!         Fill output arrays with current "non-mean" values.
!         --------------------------------------------------

          if ((sad_opt == 1) .or. (sad_opt == 2)) then
            hno3cond_nc(:,:,:) = hno3cond(:,:,kx1:kx2)
            hno3gas_nc (:,:,:) = hno3gas (:,:,kx1:kx2)

            if (sad_opt == 2) then
              h2oback_nc(:,:,:) = h2oback(:,:,kx1:kx2)
              h2ocond_nc(:,:,:) = h2ocond(:,:,kx1:kx2)
              reffice_nc(:,:,:) = reffice(:,:,kx1:kx2)
              reffsts_nc(:,:,:) = reffsts(:,:,kx1:kx2)
              vfall_nc  (:,:,:) = vfall  (:,:,kx1:kx2)
            end if
          end if

        end if

      end if

      return

      end subroutine Prep_Netcdf_SAD
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_SAD
!
! !INTERFACE:
!

      subroutine Write_Netcdf_SAD (nymd, nhms, gmi_sec)
!     
      implicit none
!     
! !INPUT PARAMETERS:
      integer, intent(in) :: nymd, nhms
      real*8 , intent(in) :: gmi_sec
!
! !DESCRIPTION:
! Writes the SAD netCDF output.
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_SAD1_NC      = 264
      integer, parameter :: SG_H2OBACK_NC   = 310
      integer, parameter :: SG_H2OCOND_NC   = 312
      integer, parameter :: SG_HNO3COND_NC  = 314
      integer, parameter :: SG_HNO3GAS_NC   = 316
      integer, parameter :: SG_REFFICE_NC   = 326
      integer, parameter :: SG_REFFSTS_NC   = 328
      integer, parameter :: SG_VFALL_NC     = 330
!
! !LOCAL VARIABLES:
      integer :: cnt2d (2), cnt4d (4)
      integer :: strt2d(2), strt4d(4)
      integer :: hdr(NETCDF_HDF)
      real*8, allocatable :: hno3gasGlob (:,:,:)
      real*8, allocatable :: hno3condGlob(:,:,:)
      real*8, allocatable :: h2obackGlob (:,:,:)
      real*8, allocatable :: h2ocondGlob (:,:,:)
      real*8, allocatable :: reffstsGlob (:,:,:)
      real*8, allocatable :: refficeGlob (:,:,:)
      real*8, allocatable :: vfallGlob   (:,:,:)
!
!EOP
!--------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_SAD called by ', procID

      if (iAmRootProc) then

         strt2d(:) = (/ 1, rnum_out_sad /)
         strt4d(:) = (/ 1, 1, 1, rnum_out_sad /)

         cnt2d (:) = (/ NETCDF_HDF, 1 /)
         cnt4d (:) = (/ numLon, numLat, numVert  , 1 /)

         hdr(1) = Nint (gmi_sec)
         hdr(2) = nymd
         hdr(3) = nhms

         call Ncwr_2d_Int (hdr, ncid_sad, hdr_var_name, strt2d, cnt2d)

         if ((sad_opt == 1) .or. (sad_opt == 2)) then
            allocate(hno3gasGlob(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))
            allocate(hno3condGlob(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))

            if (sad_opt == 2) then
               allocate(h2obackGlob(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))
               allocate(h2ocondGlob(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))
               allocate(reffstsGlob(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))
               allocate(refficeGlob(i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))
               allocate(vfallGlob  (i1_gl:i2_gl,ju1_gl:j2_gl,kx1:kx2))
            end if
         end if
      end if

      if ((sad_opt == 1) .or. (sad_opt == 2)) then

!       ========================
        call subDomain2Global  &
!       ========================
     &    (hno3condGlob, hno3cond_nc, i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1,& 
     &    j2, kx1, kx2, rootProc, procID, map1_u,  numDomains,  &
     &     SG_HNO3COND_NC, commuWorld)

!       ========================
        call subDomain2Global  &
!       ========================
     &    (hno3gasGlob, hno3gas_nc, i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, &
     &     j2, kx1, kx2, rootProc, procID, map1_u,  numDomains,  &
     &     SG_HNO3GAS_NC, commuWorld)

        if (sad_opt == 2) then

!         ========================
          call subDomain2Global  &
!         ========================
     &      (h2obackGlob, h2oback_nc, i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1,& 
     &       j2, kx1, kx2, rootProc, procID, map1_u, numDomains, &
     &       SG_H2OBACK_NC, commuWorld)

!         ========================
          call subDomain2Global  &
!         ========================
     &      (h2ocondGlob, h2ocond_nc, i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1,& 
     &       j2, kx1, kx2, rootProc, procID, map1_u, numDomains, &
     &       SG_H2OCOND_NC, commuWorld)

!         ========================
          call subDomain2Global  &
!         ========================
     &      (reffstsGlob, reffsts_nc, i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1,& 
     &       j2, kx1, kx2, rootProc, procID, map1_u, numDomains, &
     &       SG_REFFSTS_NC, commuWorld)

!         ========================
          call subDomain2Global  &
!         ========================
     &      (refficeGlob, reffice_nc, i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1,& 
     &       j2, kx1, kx2, rootProc, procID, map1_u, numDomains, &
     &       SG_REFFICE_NC, commuWorld)

!         ========================
          call subDomain2Global  &
!         ========================
     &      (vfallGlob, vfall_nc, i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2,& 
     &       kx1, kx2, rootProc, procID, map1_u, numDomains, &
     &       SG_VFALL_NC, commuWorld)

        end if
      end if

      if (iAmRootProc) then
         if ((sad_opt == 1) .or. (sad_opt == 2)) then
            call Ncwr_4d (hno3condGlob, ncid_sad, 'hno3cond', strt4d, cnt4d)
            call Ncwr_4d (hno3gasGlob,  ncid_sad, 'hno3gas',  strt4d, cnt4d)

            if (sad_opt == 2) then
               call Ncwr_4d (h2obackGlob, ncid_sad, 'h2oback', strt4d, cnt4d)
               call Ncwr_4d (h2ocondGlob, ncid_sad, 'h2ocond', strt4d, cnt4d)
               call Ncwr_4d (reffstsGlob, ncid_sad, 'reffsts', strt4d, cnt4d)
               call Ncwr_4d (refficeGlob, ncid_sad, 'reffice', strt4d, cnt4d)
               call Ncwr_4d (vfallGlob  , ncid_sad, 'vfall',   strt4d, cnt4d)
            end if
         end if

         if ((sad_opt == 1) .or. (sad_opt == 2)) then
            deallocate(hno3gasGlob)
            deallocate(hno3condGlob)

            if (sad_opt == 2) then
               deallocate(h2obackGlob)
               deallocate(h2ocondGlob)
               deallocate(reffstsGlob)
               deallocate(refficeGlob)
               deallocate(vfallGlob  )
            end if
         end if

         call Ncdo_Sync (ncid_sad)
      end if

      return

      end subroutine Write_Netcdf_SAD
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutput_SAD
!
! !INTERFACE:
!
      subroutine bufferOutput_SAD (Chemistry)
!
! !USES:
      use GmiBuffer_mod, only : bufferOutput4d
      use GmiChemistryMethod_mod, only : Get_sadgmi
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_Chemistry), intent(in) :: Chemistry
!
! !DESCRIPTION:
! Buffers the SAD data to the output file.
! The worker processors send a slice of the data to the root processor, 
! the root processor writes it out, they then send another slice and it 
! is written out, and so on.
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_SAD2_NC      = 266
!
! !LOCAL VARIABLES:
      real*8, allocatable :: sadgmi(:,:,:,:)
!
!EOP
!----------------------------------------------------------------------------------------
!BOC

      if (numSAD > 0) then
         if (.not. do_mean) then
             Allocate(sadgmi(i1:i2, ju1:j2, k1:k2, numSAD))
             call Get_sadgmi(Chemistry, sadgmi)
         end if

!        ====================
         call bufferOutput4d  &
!        ====================
     &         (sad_var_name, do_mean, kx1, kx2, commuWorld,  &
     &          SG_SAD2_NC, ncid_sad, numSAD, rnum_out_sad, sadgmi,  &
     &          sad_mean, sad_nc, map1_u, numDomains, &
     &          rootProc, procID, i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)

         if (.not. do_mean) then
            if (iAmRootProc) then
               deallocate(sadgmi)
           end if
         end if
      end if

      return

      end subroutine bufferOutput_SAD
!EOC
!-----------------------------------------------------------------------------
end module GmiControlSAD_mod
