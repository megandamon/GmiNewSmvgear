!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlRestart_mod
!
! !INTERFACE:
!
module GmiControlRestart_mod
!
! !USES:
      use  GmiSub2Glob_mod      , only : subDomain2Global
      use  GmiMessagePassing_mod, only : synchronizeGroup
      use m_ncGeneralOpsOutput  , only: Is_Out_Freq_Time
      use GmiFileOperations_mod , only : makeOutfileName
      use m_netcdf_io_close     , only : Nccl, Nccl_Noerr
      use m_netcdf_io_write, only : Ncwr_2d_Int, Ncwr_5d, Ncwr_3d, Ncwr_2d, &
     &      Ncwr_1d_Int, Ncwr_4d
      use m_netcdf_io_write     , only : Ncwr_1d, Ncwr_Scal, Ncwr_2d_Char
      use m_netcdf_io_create    , only : Ncdo_Sync, Nccr_Wr
      use m_netcdf_io_define    , only : NcDef_glob_attributes, NcDef_dimension
      use m_netcdf_io_define    , only : NcDef_var_attributes, NcDef_variable
      use m_netcdf_io_define    , only : NcSetFill, NcEnd_def
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration, Get_const_var_name
      use GmiChemistryMethod_mod , only : t_Chemistry, Get_overheadO3col,      &
     &       Get_const_labels, Get_sad_opt
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_rootProc,        &
     &       Get_procID, Get_iAmRootProc, Get_numDomains, Get_map1_u,          &
     &       Get_communicatorWorld, Get_latdeg, Get_londeg

      use GmiGrid_mod, only : t_gmiGrid, Get_numSpecies, Get_i1, Get_i2,       &
     &       Get_ju1, Get_j2, Get_k1, Get_k2, Get_i1_gl, Get_i2_gl, Get_ju1_gl,&
     &       Get_j2_gl, Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jhi_gl
      use GmiPrintError_mod        , only : GmiPrintError
!
      use GmiTimeControl_mod, only : t_GmiClock, GmiSplitDateTime, isOutTime, &
     &       Get_curGmiDate, Get_curGmiTime, Get_gmiSeconds, Get_gmiTimeStep, &
     &       Get_numTimeSteps
      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Get_pr_diag,         &
     &       Get_problem_name, Get_hdr_var_name, Get_hdf_dim_name,             &
     &       Get_rec_dim_name, Get_lat_dim_name, Get_lon_dim_name,             &
     &       Get_prs_dim_name, Get_spc_dim_name, Get_do_overwrt_rst,           &
     &       Get_pr_rst_period, Get_restart_infile_name, Get_pr_restart,       &
     &       Get_pr_qqjk
      use GmiMetFieldsControl_mod, only : t_metFields, Get_pt, Get_mdt, &
     &       Get_met_opt, Get_do_cycle_met, Get_am, Get_bm, Get_m1num_recs,    &
     &       Get_met_num_infiles, Get_tmet1_1, Get_tmet1_2, Get_m1rnum_in,     &
     &       Get_met1_infile_num, Get_pctm2Glob

      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: controlOutputRestart, initializeOutputRestart
      public  :: finalizeOutputRestart
!
#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
!
! !DESCRIPTION:
! Contains the necessary routines to produce the restart output file. 
! Three routines are visible from the outside:
! %
! \begin{description}
! \item[initializeOutputRestart:] called once in the initialization steps to
!      open the netCDF file, define variables in the file, write out header
!      information and allocate variables.
! \item[controlOutputRestart:] called at the end of each model time step to
!      update variables and if needed for communications and writing out
!      data in the file.
! \item[finalizeOutputRestart:] called at the end of the integration to
!      deallocate all the variables and close the netCDF file.
! \end{description}
! %
! The above routines have as arguments derived types.
! %
! This module is self-contained. All the main routines are called by all
! the processors (master included) but each processor does not execute
! all portions of the code.
!
                                ! 3D const slice (mixing ratio)
      real*8 , pointer, save :: const_rst    (:,:,:) => null()
                                ! condensed phase h2o array (mixing ratio)
      real*8 , pointer, save :: h2ocond_rst  (:,:,:) => null()

      integer,          save :: ncid_rst     ! NetCDF       species output file id
      integer,          save :: rnum_out_rst ! next Netcdf output record to write

      ! Grid information of the variables to be written out
      integer,          save :: i1, i2     ! longitude
      integer,          save :: i1_gl, i2_gl
      integer,          save :: ilo_gl, ihi_gl
      integer,          save :: ju1        ! latitude
      integer,          save :: j2
      integer,          save :: ju1_gl, j2_gl
      integer,          save :: julo_gl, jhi_gl
      integer,          save :: k1, k2   ! vertical
      integer,          save :: numLon     ! number of longitudes
      integer,          save :: numLat     ! number of latitudes
      integer,          save :: numVert    ! number of vertical levels

      logical,          save :: iAmRootProc
      integer,          save :: procID, rootProc
      integer,          save :: numDomains, subdomain
      integer,          save :: commuWorld
      integer, pointer, save :: map1_u       (:,:,:) => null()

      integer,          save :: numSpecies

      integer,            save :: nhdf, sad_opt
      logical,            save :: pr_diag
      logical,            save :: do_overwrt_rst, pr_restart, pr_qqjk
      real*8,             save :: pr_rst_period
      real*8,             save :: mdt
      integer,            save :: met_opt, met_num_infiles
      logical,            save :: do_cycle_met
      character (len=MAX_LENGTH_VAR_NAME), save :: rec_dim_name, prs_dim_name, spc_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: hdr_var_name, hdf_dim_name, lat_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: lon_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: const_var_name
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
! !IROUTINE: initializeOutputRestart
!
! !INTERFACE:
!
      subroutine initializeOutputRestart(SpeciesConcentration, gmiGrid, &
     &                     gmiDomain, Diagnostics, Chemistry, metFields)
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
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !DESCRIPTION:
! The routines (1) opens a netCDF file, (2) defines variables in the file,
! (3) writes out header data in the file, and (4) allocates varaibles.
!
! !LOCAL VARIABLES:
      character (len=75)  :: err_msg
      character (len=MAX_LENGTH_FILE_NAME) :: fname
      character (len=MAX_LENGTH_FILE_NAME) :: problem_name
      character (len=MAX_LENGTH_FILE_NAME) :: restart_infile_name
      character (len=MAX_LENGTH_SPECIES_NAME), pointer :: const_labels(:)
!EOP
!-----------------------------------------------------------------------------
!BOC
      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)
      
      if (pr_diag) Write (6,*) 'initializeOutputRestart called by ', procID

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
      call Get_numSpecies (gmiGrid, numSpecies )

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

      numLon  = i2_gl - i1_gl  + 1
      numLat  = j2_gl - ju1_gl + 1
      numVert = k2    - k1     + 1

      nhdf = NETCDF_HDF

      call Get_pr_qqjk      (Diagnostics, pr_qqjk     )
      call Get_pr_restart   (Diagnostics, pr_restart  )
      call Get_hdf_dim_name (Diagnostics, hdf_dim_name)
      call Get_lat_dim_name (Diagnostics, lat_dim_name)
      call Get_lon_dim_name (Diagnostics, lon_dim_name)
      call Get_spc_dim_name (Diagnostics, spc_dim_name)
      call Get_rec_dim_name (Diagnostics, rec_dim_name)
      call Get_prs_dim_name (Diagnostics, prs_dim_name)
      call Get_pr_rst_period(Diagnostics, pr_rst_period)
      call Get_do_overwrt_rst(Diagnostics, do_overwrt_rst)

      call Get_const_var_name(SpeciesConcentration, const_var_name)

      call Get_sad_opt(Chemistry, sad_opt)

      call Get_mdt            (metFields, mdt)
      call Get_met_opt        (metFields, met_opt)
      call Get_do_cycle_met   (metFields, do_cycle_met)
      call Get_met_num_infiles(metFields, met_num_infiles)

      allocate(const_labels(numSpecies))

      if (iAmRootProc) then

         !###########################
         ! Initialize the output file
         !###########################

         ! Determine the file name
         call Get_problem_name(Diagnostics, problem_name)
         call makeOutfileName (fname, '.rst.nc', problem_name)

         call Get_restart_infile_name(Diagnostics, restart_infile_name)
         if (Trim (restart_infile_name) == Trim (fname)) then
            err_msg = 'restart_infile_name/fname problem in Init_Restart.'
            call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
         end if

         ! Create the file and assign a file identifier
         call Nccr_Wr (ncid_rst, fname)

         ! Define the variables in the file
         call Define_Netcdf_Out_Restart(const_var_name)

         ! Write header data
         call Get_const_labels(Chemistry, const_labels)

         call Write_Netcdf_Hdr_Restart (gmiDomain, metFields, const_labels)

         call Ncdo_Sync (ncid_rst)

         ! Initialize the counter for the number of records
         rnum_out_rst = 1
      end if

      !########################
      ! Allocation of variables
      !########################

      call allocateVariablesRestart()

      return

      end subroutine initializeOutputRestart
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: controlOutputRestart
!
! !INTERFACE:
!
      subroutine controlOutputRestart(last_tstp, &
     &                  Chemistry, SpeciesConcentration, gmiClock, metFields)
!
! !USES:
      use GmiChemistryMethod_mod, only : t_Chemistry
!
      implicit none
!
! !INPUT PARAMETERS:
      logical                     , intent(in) :: last_tstp
      type(t_GmiClock)            , intent(in) :: GmiClock
      type(t_metFields)           , intent(in) :: metFields
      type(t_Chemistry)           , intent(in) :: Chemistry
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !DESCRIPTION:
! This routine controls the restart file output. It is called at each time step.
!
! !LOCAL VARIABLES:
      logical :: time_for_rst
      integer :: day, month, idumyear, ndt, num_time_steps, nhms, nymd
      real*8  :: gmi_sec, tstep
      integer, save :: month_save = -999
      logical, save :: printed_on_this_day = .false.
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'controlOutputRestart called by ', procID

      call Get_curGmiDate(gmiClock, nymd)
      call Get_curGmitime(gmiClock, nhms)
      call Get_numTimeSteps(gmiClock, num_time_steps)
      call Get_gmiSeconds(gmiClock, gmi_sec)
      call Get_gmiTimeStep(gmiClock, tstep)

      ndt = Nint(tstep)
      
      call GmiSplitDateTime (nymd, idumyear, month, day)

      if (month_save == -999) month_save = month

      time_for_rst = .false.

      call isOutTime (time_for_rst, printed_on_this_day, &
     &          month_save, month, day, nhms, ndt, gmi_sec, pr_rst_period)

      month_save = month

      if (last_tstp .and. (.not. time_for_rst)) then
        if (Mod (num_time_steps, Nint (mdt) / ndt) == 0) then

!         ------------------------------------------------------------
!         Always update restart file after the last step if you are at
!         the end of a met record.
!         ------------------------------------------------------------

          time_for_rst = .true.
        end if
      end if


!     =================
      if (time_for_rst) then
!     =================

         call Prep_Netcdf_Restart (Chemistry)

         call Write_Netcdf_Restart (metFields, nymd, nhms, gmi_sec)

         call bufferOutput_Restart (const_var_name, SpeciesConcentration)

        if (procID == rootProc) then
           if(.not. do_overwrt_rst) rnum_out_rst = rnum_out_rst + 1
        end if

!     ======
      end if
!     ======

      return

      end subroutine controlOutputRestart
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finalizeOutputRestart
!
! !INTERFACE:
!
      subroutine finalizeOutputRestart()
!
      implicit none
!
! !DESCRIPTION:
! Deallocates variables necessary to produce restart outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if (pr_diag) Write(6,*) 'finalizeOutputRestart called by ', procID

      !==========================
      ! Deallocation of variables
      !==========================

      deallocate (const_rst)

      if (sad_opt == 2) then
         deallocate(h2ocond_rst)
      end if

      !======================
      ! Close the netCDF file
      !======================
      if (iAmRootProc) call Nccl_Noerr (ncid_rst)
      
      return
  
      end subroutine finalizeOutputRestart
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocateVariablesRestart
!
! !INTERFACE:
!
      subroutine allocateVariablesRestart()
!
      implicit none
!
! !DESCRIPTION:
! Allocates variables necessary to produce restart outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'allocateVariablesRestart called by ', procID

      Allocate(const_rst(i1:i2, ju1:j2, k1:k2))
      const_rst(:,:,:) = 0.0d0

      if (sad_opt == 2) then
         Allocate(h2ocond_rst(i1:i2, ju1:j2, k1:k2))
         h2ocond_rst(:,:,:) = 0.0d0
      end if

      return

      end subroutine allocateVariablesRestart
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Define_Netcdf_Out_Restart
!
      subroutine Define_Netcdf_Out_Restart  (const_var_name)
!         
! !USES:
      use m_ncGeneralOpsOutput, only: Define_Netcdf_Out_Gen
!       
      implicit none
!       
#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
      character (len=* ) , intent(in) :: const_var_name
!    
! !DESCRIPTION:
! This routine makes the necessary definitions for the restart
! netCDF output file.
!
! !LOCAL VARIABLES:
      integer :: ierr
      integer :: nchr, nspc
      integer :: omode
      integer :: varid

      integer :: chrd(1)
      integer :: hdfd(1)
      integer :: lond(1), latd(1), prsd(1)
      integer :: spcd(1), recd(1)

      integer :: var1(1), var2(2), var3(3)
      integer :: var4(4), var5(5)

!
!EOP
!-------------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Define_Netcdf_Out_Restart called by ', procID

      nchr = MAX_LENGTH_SPECIES_NAME
      nspc = numSpecies

!     ==========================
      call Define_Netcdf_Out_Gen  &
!     ==========================
     &  (ncid_rst, nhdf, numLon, numLat, numVert, hdfd, lond, latd, prsd,  &
     &   hdf_dim_name, lon_dim_name, lat_dim_name, prs_dim_name)

!     ------------------
!     Define dimensions.
!     ------------------

      call NcDef_dimension(ncid_rst, 'chr_dim', nchr, chrd(1))
!                                     ---------
      call NcDef_dimension(ncid_rst, spc_dim_name, nspc, spcd(1))
!                                     ------------
      call NcDef_dimension(ncid_rst, rec_dim_name, NF_UNLIMITED, recd(1))
!                                     ------------
!     -----------------------------------------
!     Define variables and variable attributes.
!     -----------------------------------------

      call NcDef_variable (ncid_rst, spc_dim_name, NF_FLOAT, 1, spcd, varid)
!                                    ------------
      call NcDef_var_attributes (ncid_rst, varid, 'long_name', 'Species index')
      call NcDef_var_attributes (ncid_rst, varid, 'units', 'unitless')
      call NcDef_var_attributes  &
     &          (ncid_rst, varid, 'coord_labels', 'const_labels')
      call NcDef_var_attributes (ncid_rst, varid, 'selection_category', 'NULL')

      var2(:) = (/ chrd(1), spcd(1) /)
      call NcDef_variable (ncid_rst, 'const_labels', NF_CHAR, 2, var2, varid)
!                                     --------------
      call NcDef_var_attributes (ncid_rst, varid, 'long_name', 'Species name')
      call NcDef_var_attributes (ncid_rst, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_rst, varid, 'selection_category', 'NULL')


      var1(:) = (/ recd(1) /)
      call NcDef_variable (ncid_rst, 'tbegin_days', NF_FLOAT, 1, var1, varid)
!                                     -------------
      call NcDef_var_attributes (ncid_rst, varid, 'long_name',  &
     &          'Starting time of the simulation')
      call NcDef_var_attributes (ncid_rst, varid, 'units', 'days')

      var1(:) = (/ recd(1) /)
      call NcDef_variable (ncid_rst, 'nymd', NF_INT, 1, var1, varid)
!                                     ------
      call NcDef_var_attributes (ncid_rst, varid, 'long_name',  &
     &          'Starting year/month/day')
      call NcDef_var_attributes (ncid_rst, varid, 'format', 'yyyymmdd')

      var1(:) = (/ recd(1) /)
      call NcDef_variable (ncid_rst, 'nhms', NF_INT, 1, var1, varid)
!                                     ------
      call NcDef_var_attributes (ncid_rst, varid, 'long_name',  &
     &          'Starting hour/minute/second')
      call NcDef_var_attributes (ncid_rst, varid, 'format', 'hhmmss')

      var1(:) = (/ recd(1) /)
      call NcDef_variable (ncid_rst, 'gmi_sec', NF_FLOAT, 1, var1, varid)
!                                     ---------
      call NcDef_var_attributes (ncid_rst, varid, 'long_name',  &
     &          'Total elapsed simulation time')
      call NcDef_var_attributes (ncid_rst, varid, 'units', 's')

      var1(:) = (/ recd(1) /)
      call NcDef_variable (ncid_rst, 'pr_qqjk_int', NF_INT, 1, var1, varid)
!                                     -------------
      call NcDef_var_attributes  (ncid_rst, varid, 'long_name', 'qqjk file write flag')
      call NcDef_var_attributes (ncid_rst, varid, 'format', '0|1')

      var1(:) = (/ recd(1) /)
      call NcDef_variable (ncid_rst, 'met_infile_num', NF_INT, 1, var1, varid)
!                                     ----------------
      call NcDef_var_attributes (ncid_rst, varid, 'long_name',  &
     &          'Met data input file index')
      call NcDef_var_attributes (ncid_rst, varid, 'units', 'unitless')

      var1(:) = (/ recd(1) /)
      call NcDef_variable (ncid_rst, 'mrnum_in', NF_INT, 1, var1, varid)
!                                     ----------
      call NcDef_var_attributes (ncid_rst, varid, 'long_name',  &
     &          'Met data input record index')
      call NcDef_var_attributes (ncid_rst, varid, 'units', 'unitless')

      var1(:) = (/ recd(1) /)
      call NcDef_variable (ncid_rst, 'tmet1', NF_FLOAT, 1, var1, varid)
!                                     -------
      call NcDef_var_attributes (ncid_rst, varid, 'long_name',  &
     &          'Time tag of last met data read')
      call NcDef_var_attributes (ncid_rst, varid, 'units', 's')

      var3(:) = (/ lond(1), latd(1), recd(1) /)
      call NcDef_variable (ncid_rst, 'pctm2', NF_FLOAT, 3, var3, varid)
!                                     -------

      if (sad_opt == 2) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_rst, 'h2ocond', NF_FLOAT, 4, var4, varid)
!                                       ---------
      end if

      var5(:) = (/ lond(1), latd(1), prsd(1), spcd(1), recd(1) /)
      call NcDef_variable (ncid_rst, const_var_name, NF_FLOAT, 5, var5, varid)
!                                     --------------
      call NcDef_var_attributes (ncid_rst, varid, 'long_name', 'Constituent')
      call NcDef_var_attributes (ncid_rst, varid, 'units', 'volume mixing ratio')

!     -------------------------
!     Define global attributes.
!     -------------------------

      call NcDef_glob_attributes (ncid_rst, 'title', 'Gmimod restart file')

      call NcSetFill (ncid_rst, NF_NOFILL, omode)

      call NcEnd_def (ncid_rst)

      return

      end subroutine Define_Netcdf_Out_Restart
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Hdr_Restart
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Hdr_Restart (gmiDomain, metFields, const_labels)
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
! This routine creates some header information for the restart netCDF output
! file and writes it out.
! 
! !LOCAL VARIABLES:
      integer :: ic
      integer :: cnt1d (1), cnt2d (2)
      integer :: strt1d(1), strt2d(2)
      real*8  :: prsdat(k1:k2)
      real*8  :: spcdat(numSpecies)
      real*8 , allocatable :: londeg(:), latdeg(:)
      real*8 , allocatable :: bm(:), am(:)
      real*8 :: pt
!
!EOP
!---------------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Hdr_Restart called by ', procID

      allocate(latdeg(ju1_gl:j2_gl))
      call Get_latdeg(gmiDomain, latdeg)

      allocate(londeg(i1_gl:i2_gl))
      call Get_londeg(gmiDomain, londeg)

!     =========================
      call Write_Netcdf_Hdr_Gen  &
!     =========================
     &  (ncid_rst, latdeg, londeg, pr_diag, procID, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   hdf_dim_name, lat_dim_name, lon_dim_name)

!     ---------
!     Pressure.
!     ---------

      call Get_pt(metFields, pt)

      allocate(am(k1:k2))
      call Get_am(metFields, am)

      allocate(bm(k1:k2))
      call Get_bm(metFields, bm)

      prsdat(:) = (am(:) * pt)  + (bm(:) * 1000.0d0)

      strt1d(1) = 1
      cnt1d (1) = numVert

      call Ncwr_1d (prsdat, ncid_rst, prs_dim_name, strt1d, cnt1d)

!     ----------------
!     Species numbers.
!     ----------------

      do ic = 1, numSpecies
        spcdat(ic) = ic
      end do

      strt1d(1) = 1
      cnt1d (1) = numSpecies

      call Ncwr_1d (spcdat, ncid_rst, spc_dim_name, strt1d, cnt1d)

!     ---------------
!     Species labels.
!     ---------------

      strt2d(:) = (/ 1, 1 /)
      cnt2d (:) = (/ MAX_LENGTH_SPECIES_NAME, numSpecies /)

      call Ncwr_2d_Char (const_labels, ncid_rst, 'const_labels', strt2d, cnt2d)

      return

      end subroutine Write_Netcdf_Hdr_Restart
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Prep_Netcdf_Restart
!
! !INTERFACE:
!
      subroutine Prep_Netcdf_Restart (Chemistry)
!
! !USES:
      use GmiChemistryMethod_mod, only : t_Chemistry, Get_h2ocond

      implicit none
!
! !INPUT PARAMETERS:
      type(t_Chemistry), intent(in) :: Chemistry
!
! !DESCRIPTION:
! Prepares the restart netCDF outputs.
!
!EOP
!-----------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Prep_Netcdf_Restart called by ', procID

      if (sad_opt == 2) then
         call Get_h2ocond(Chemistry, h2ocond_rst)
      end if

      return

      end subroutine Prep_Netcdf_Restart
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Restart
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Restart (metFields, nymd, nhms, gmi_sec)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: nymd, nhms
      real*8 , intent(in) :: gmi_sec
      type (t_metFields), intent(in) :: metFields
!
! !DESCRIPTION:
! Writes the restart netCDF outputs.
!
! DEFINED PARAMETERS:
      integer, parameter :: SG_H2OCOND_RST = 5555
!
! !LOCAL VARIABLES:
      integer :: file_num
      integer :: pr_qqjk_int
      integer :: rec_num
      integer :: iout(1)
      integer :: cnt1d (1), cnt3d (3), cnt4d (4)
      integer :: strt1d(1), strt3d(3), strt4d(4)
      real*8  :: rdays
      real*8  :: time1
      real*8  :: rout(1)
      integer :: m1rnum_in, m1num_recs, met1_infile_num
      real*8  :: tmet1_1, tmet1_2
      real*8  :: pctm2_rst(i1_gl:i2_gl, ju1_gl:j2_gl)
      real*8, allocatable :: h2ocond_rstGlob(:,:,:)
      real*8, allocatable :: pctm2Glob(:,:)
!
!EOP
!---------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Restart called by ', procID
   
      if (sad_opt == 2) then

         if (iAmRootProc) then
            allocate(h2ocond_rstGlob(i1_gl:i2_gl,ju1_gl:j2_gl,k1:k2))
            h2ocond_rstGlob(:,:,:) = 0.0d0
         end if

         call subDomain2Global(h2ocond_rstGlob, h2ocond_rst, i1_gl, i2_gl, &
     &          ju1_gl, j2_gl, i1, i2, ju1, j2, k1, k2,  &
     &          rootProc, procID, map1_u, numDomains, SG_H2OCOND_RST, commuWorld)
      end if

      if (iAmRootProc) then
         call Get_met_opt         (metFields, met_opt        )
         call Get_do_cycle_met    (metFields, do_cycle_met   )
         call Get_met_num_infiles (metFields, met_num_infiles   )

         call Get_tmet1_1         (metFields, tmet1_1        )
         call Get_tmet1_2         (metFields, tmet1_2        )
         call Get_m1rnum_in       (metFields, m1rnum_in      )
         call Get_m1num_recs      (metFields, m1num_recs     )
         call Get_met1_infile_num (metFields, met1_infile_num)

         allocate(pctm2Glob(ilo_gl:ihi_gl,julo_gl:jhi_gl))
         call Get_pctm2Glob(metFields, pctm2Glob)

         strt1d(:) = (/          rnum_out_rst /)
         strt3d(:) = (/    1, 1, rnum_out_rst /)
         strt4d(:) = (/ 1, 1, 1, rnum_out_rst /)

         cnt1d (:) = (/ 1 /)
         cnt3d (:) = (/ numLon, numLat, 1 /)
         cnt4d (:) = (/ numLon, numLat, numVert, 1 /)

         rdays = gmi_sec / SECPDY

         rout(1) = rdays
         call Ncwr_1d     (rout, ncid_rst, 'tbegin_days', strt1d, cnt1d)

         iout(1) = nymd
         call Ncwr_1d_Int (iout, ncid_rst, 'nymd',        strt1d, cnt1d)

         iout(1) = nhms
         call Ncwr_1d_Int (iout, ncid_rst, 'nhms',        strt1d, cnt1d)

         rout(1) = gmi_sec
         call Ncwr_1d     (rout, ncid_rst, 'gmi_sec',     strt1d, cnt1d)

         if (pr_qqjk) then
           pr_qqjk_int = 1
         else
           pr_qqjk_int = 0
         end if

         iout(1) = pr_qqjk_int
         call Ncwr_1d_Int (iout, ncid_rst, 'pr_qqjk_int', strt1d, cnt1d)

         if (met_opt /= 1) then

            if (m1rnum_in > m1num_recs) then
               file_num = met1_infile_num + 1
               if (do_cycle_met .and. (file_num > met_num_infiles)) then
                  file_num = 1
               end if
               rec_num  = 1
               time1    = tmet1_2
            else
               file_num = met1_infile_num
               rec_num  = m1rnum_in - 1
               time1    = tmet1_2
            end if

         else

           file_num = met1_infile_num
           rec_num  = m1rnum_in
           time1    = tmet1_1

         end if

         iout(1) = file_num
         call Ncwr_1d_Int (iout, ncid_rst, 'met_infile_num', strt1d, cnt1d)

         iout(1) = rec_num
         call Ncwr_1d_Int (iout, ncid_rst, 'mrnum_in',       strt1d, cnt1d)

         rout(1) = time1
         call Ncwr_1d (rout, ncid_rst, 'tmet1', strt1d, cnt1d)

         pctm2_rst(:,:) = pctm2Glob(i1_gl:i2_gl,ju1_gl:j2_gl)
   
         call Ncwr_3d (pctm2_rst, ncid_rst, 'pctm2', strt3d, cnt3d)
   
         if (sad_opt == 2) then
           call Ncwr_4d (h2ocond_rstGlob, ncid_rst, 'h2ocond', strt4d, cnt4d)
         end if

!        ==============
         call Ncdo_Sync (ncid_rst)
!        ==============

      end if

      return

      end subroutine Write_Netcdf_Restart
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutput_Restart
!
! !INTERFACE:
!
      subroutine bufferOutput_Restart (const_var_name, SpeciesConcentration)
!
! !USES:
      use m_netcdf_io_create  , only : Ncdo_Sync
      use GmiNcOutputSlice_mod, only : writeSlice4d

      implicit none
!
! !INPUT PARAMETERS:
      character (len=* ) , intent(in) :: const_var_name
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !DESCRIPTION:
!   This routine buffers the restart data to the output file.  
!   The worker processors send a slice of the data to the root processor, 
!   the root processor writes it out, they then send another slice and it 
!   is written out, and so on.
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_CONST2_RST = 2210
!
! !LOCAL VARIABLES:
      integer :: ic
      type (t_GmiArrayBundle), pointer :: concentration(:)
      real*8, allocatable :: const_rstGlob(:,:,:)
!
!EOP
!-----------------------------------------------------------------------------------
!BOC
      call Get_concentration(SpeciesConcentration, concentration)

      if (iAmRootProc) allocate(const_rstGlob(i1_gl:i2_gl,ju1_gl:j2_gl,k1:k2))

      do ic = 1, numSpecies
         if (iAmRootProc) const_rstGlob(:,:,:) = 0.0d0

         const_rst(:,:,:) = concentration(ic)%pArray3D(:,:,:)

         call subDomain2Global (const_rstGlob, const_rst, i1_gl, i2_gl, &
     &          ju1_gl, j2_gl, i1, i2, ju1, j2, k1, k2, rootProc, &
     &           procID, map1_u, numDomains, SG_CONST2_RST, commuWorld)

         if (iAmRootProc) then
            call writeSlice4d (i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2, &
     &                const_rstGlob, ncid_rst, const_var_name, ic, rnum_out_rst)
         end if

         call synchronizeGroup (commuWorld)

      end do

      if (iAmRootProc) then
         call Ncdo_Sync (ncid_rst)
         deallocate(const_rstGlob)
      end if

      return

      end subroutine bufferOutput_Restart
!EOC
!-----------------------------------------------------------------------------
end module GmiControlRestart_mod
