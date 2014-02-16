!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
      module GmiControlFlux_mod

      use  GmiSub2Glob_mod, only : subDomain2Global
      use GmiFileOperations_mod, only : makeOutfileName

      use m_netcdf_io_close, only : Nccl, Nccl_Noerr
      use m_netcdf_io_write, only : Ncwr_2d_Int, Ncwr_5d, Ncwr_3d, Ncwr_2d
      use m_netcdf_io_write, only : Ncwr_1d, Ncwr_Scal, Ncwr_2d_Char
      use m_netcdf_io_create, only : Ncdo_Sync, Nccr_Wr
      use m_netcdf_io_define, only : NcDef_glob_attributes, NcDef_dimension
      use m_netcdf_io_define, only : NcDef_var_attributes, NcDef_variable
      use m_netcdf_io_define, only : NcSetFill, NcEnd_def
      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Set_do_flux_reset,   &
     &       Get_pr_flux, Get_pr_psf_flux, Get_pr_const_flux, Get_flux_name,   &
     &       Get_flux_species, Get_flux_species_num, Get_flux_var_name,      &
     &       Get_hdf_dim_name, Get_lon_dim_name, Get_lat_dim_name,             &
     &       Get_fluxOutputFrequency, Get_hdr_var_name, Get_rec_dim_name,             &
     &       Get_prs_dim_name,  Get_k1r_gl, Get_k2r_gl, Get_pr_diag,           &
     &       Get_problem_name

      use GmiChemistryMethod_mod, only : t_Chemistry, Get_const_labels

      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_rootProc,        &
     &       Get_procID, Get_iAmRootProc,          &
     &       Get_numDomains, Get_communicatorWorld, Get_map1_u
      use GmiDomainDecomposition_mod, only : Get_mcor, Get_mcorGlob, &
     &       Get_latdeg, Get_londeg

      use GmiGrid_mod, only : t_gmiGrid, Get_numSpecies, Get_i1, Get_i2,       &
     &       Get_ju1, Get_j2, Get_k1, Get_k2, Get_i1_gl, Get_i2_gl, Get_ju1_gl,&
     &       Get_j2_gl, Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jhi_gl

      use GmiTimeControl_mod, only : t_GmiClock, GmiSplitDateTime, isOutTime, &
     &       Get_curGmiDate, Get_curGmiTime, Get_gmiSeconds, Get_gmiTimeStep, &
     &       Get_numTimeSteps

      use GmiMetFieldsControl_mod, only : t_metFields, Get_pt,  &
     &       Get_metdata_name, Get_mdt, Get_ai, Get_bi, Get_am, Get_bm

      use GmiAdvectionMethod_mod, only : t_Advection, resetFlux, &
     &       Get_flux_x, Set_flux_x, Get_flux_y, Set_flux_y, &
     &       Get_flux_z, Set_flux_z, Get_psf_flux, Set_psf_flux, &
     &       Get_count_flux, Set_count_flux, Get_air_mass_flux, Set_air_mass_flux

      implicit none

      private
      public  :: initializeOutputFlux, controlOutputFlux, finalizeOutputFlux

#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"

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

      ! netCDF file identifier
      integer,          save :: ncid_flux

      ! Counter for the number of records 
      integer,          save :: rnum_out_flux

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
      integer,          save :: numLon     ! number of local longitudes 
      integer,          save :: numLat     ! number of local latitudes
      integer,          save :: numVert    ! number of vertical levels
      integer,          save :: numLonGlob ! number of global longitudes 
      integer,          save :: numLatGlob ! number of global latitudes
      integer,          save :: numSpecies ! number of species

      logical,          save :: iAmRootProc
      integer,          save :: procID, rootProc
      integer,          save :: numDomains
      integer,          save :: commuWorld
      integer, pointer, save :: map1_u       (:,:,:) => null()
      real*8 , pointer, save :: mcor(:,:) => null()

      integer,          save :: nhdf
      logical,          save :: pr_diag

      real*8, pointer, save :: air_mass_flux_nc  (:,:,:,:) => null()
      real*8, pointer, save :: flux_x_nc (:,:,:,:) => null()
      real*8, pointer, save :: flux_y_nc (:,:,:,:) => null()
      real*8, pointer, save :: flux_z_nc (:,:,:,:) => null()
      real*8, pointer, save :: psf_flux_nc  (:,:) => null()

      character(len=80), save :: flux_name
      logical          , save :: pr_flux 
      logical          , save :: do_flux_reset
      logical          , save :: pr_const_flux
      logical          , save :: pr_psf_flux
      integer, pointer , save :: flux_species(:) => null()
      integer          , save :: flux_species_num
      real*8           , save :: fluxOutputFrequency
      character(len=MAX_LENGTH_VAR_NAME), save :: flux_var_name
      character (len=MAX_LENGTH_VAR_NAME), save :: rec_dim_name, prs_dim_name, hdr_var_name
      character (len=MAX_LENGTH_VAR_NAME), save :: hdf_dim_name, lon_dim_name, lat_dim_name

      real*8           , save :: mdt
!EOP
!-----------------------------------------------------------------------------
      contains
!-------------------------------------------------------------------------------
!BOP
! 
! !IROUTINE: initializeOutputFlux
! 
! !INTERFACE:
! 
      subroutine initializeOutputFlux(gmiGrid, gmiDomain, Diagnostics, &
     &                     Chemistry, metFields)
!     
! !USES:!
      implicit none
!     
! !INPUT PARAMETERS:
      type(t_gmiGrid    ), intent(in) :: gmiGrid
      type(t_metFields  ), intent(in) :: metFields
      type(t_gmiDomain  ), intent(in) :: gmiDomain
      type(t_Chemistry  ), intent(in) :: Chemistry  
      type(t_Diagnostics), intent(in) :: Diagnostics
!     
! !DESCRIPTION:
! The routines (1) opens a netCDF file, (2) defines variables in the file,
! (3) writes out header data in the file, and (4) allocates varaibles.
!     
! !LOCAL VARIABLES:     
      character (len=75)  :: err_msg       
      integer :: in1, k1r_gl, k2r_gl, ic
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
      
      if (pr_diag) Write (6,*) 'initializeOutputFlux called by ', procID
      
      call Get_pr_flux      (Diagnostics, pr_flux )

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

      call Get_k1r_gl      (Diagnostics, k1r_gl)
      call Get_k2r_gl      (Diagnostics, k2r_gl)

      kx1     = k1 ! k1r_gl
      kx2     = k2 ! k2r_gl

      numLon  = i2_gl - i1_gl + 1
      numLat  = j2_gl - ju1_gl + 1
      numVert = kx2 - kx1 + 1

      numLonGlob  = i2_gl - i1_gl + 1
      numLatGlob  = j2_gl - ju1_gl + 1

      nhdf = NETCDF_HDF

      call Get_pr_psf_flux  (Diagnostics, pr_psf_flux )
      call Get_pr_const_flux(Diagnostics, pr_const_flux )

      allocate(flux_species(numSpecies))
      call Get_flux_species (Diagnostics, flux_species)
      call Get_flux_species_num (Diagnostics, flux_species_num)
      call Get_fluxOutputFrequency (Diagnostics, fluxOutputFrequency)
      call Get_flux_var_name(Diagnostics, flux_var_name)

      call Get_hdr_var_name(Diagnostics, hdr_var_name)
      call Get_rec_dim_name(Diagnostics, rec_dim_name)
      call Get_prs_dim_name(Diagnostics, prs_dim_name)
      call Get_hdf_dim_name (Diagnostics, hdf_dim_name)
      call Get_lat_dim_name (Diagnostics, lat_dim_name)
      call Get_lon_dim_name (Diagnostics, lon_dim_name)

      allocate(mcor(i1:i2,ju1:j2))
      call Get_mcor(gmiDomain, mcor)

      call Get_mdt(metFields, mdt)

      if (iAmRootProc) then

         !###########################
         ! Initialize the output file
         !###########################

         ! Determine the file name
         call Get_flux_name   (Diagnostics, flux_name)
         call Get_problem_name(Diagnostics, problem_name)

         in1 = Len_Trim (flux_name)
         call makeOutfileName (fname, '.' // flux_name(1:in1) // '.nc', problem_name)

         ! Create the file and assign a file identifier
         call Nccr_Wr (ncid_flux, fname)

         ! Define the variables in the file
         call Define_Netcdf_Out_Flux( )

         ! Write header data

         call Write_Netcdf_Hdr_Flux (gmiDomain, Chemistry, metFields)

         call Ncdo_Sync (ncid_flux)

         ! Initialize the counter for the number of records
         rnum_out_flux = 1
      end if

      !########################
      ! Allocation of variables
      !########################
      call allocateOutputFlux()

      return

      end subroutine initializeOutputFlux
!EOC
!-----------------------------------------------------------------------------
!
! ROUTINE
!   controlOutputFlux
!
! DESCRIPTION
!   This routine controls the flux file output.
!
! ARGUMENTS
!   last_tstp : last time step?
!
!-----------------------------------------------------------------------------

      subroutine controlOutputFlux (last_tstp, Diagnostics, gmiClock, Advection)
!
! !USES:
      use GmiTimeControl_mod, only : GmiSplitDateTime
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: last_tstp
      type (t_gmiClock), intent(in) :: gmiClock
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Advection  ), intent(inOut) :: Advection  
      type (t_Diagnostics), intent(inOut) :: Diagnostics
!
! !LOACAL VARIABLES:
      logical :: time_for_flux
      integer :: day, month, idumyear, ndt, num_time_steps, nhms, nymd
      real*8  :: gmi_sec, tstep
      integer, save :: month_save = -999
      logical, save :: printed_on_this_day = .false.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'controlOutputFlux called by ', procID

      call Get_curGmiDate(gmiClock, nymd)
      call Get_curGmitime(gmiClock, nhms)
      call Get_numTimeSteps(gmiClock, num_time_steps)
      call Get_gmiSeconds(gmiClock, gmi_sec)
      call Get_gmiTimeStep(gmiClock, tstep)

      ndt = Nint(tstep)

      call GmiSplitDateTime (nymd, idumyear, month, day)

      if (month_save == -999) month_save = month

      time_for_flux = .false.

      call isOutTime (time_for_flux, printed_on_this_day, &
     &       month_save, month, day, nhms, ndt, gmi_sec, fluxOutputFrequency)

      month_save = month

      if (last_tstp .and. (.not. time_for_flux)) then
        if (Mod (num_time_steps, Nint (mdt) / ndt) == 0) then

!         ------------------------------------------------------------
!         Always update restart file after the last step if you are at
!         the end of a met record.
!         ------------------------------------------------------------

          time_for_flux = .true.

        endif
      endif

!     ================
      if (time_for_flux) then
!     ================        
        
         call Prep_Netcdf_Flux (Advection)

         call Write_Netcdf_Flux (nymd, nhms, gmi_sec)

         call bufferOutput_Flux (Advection)

         if (iAmRootProc) rnum_out_flux = rnum_out_flux + 1

         call Set_do_flux_reset(Diagnostics, .true.)

         call resetFlux(Advection, Diagnostics)
      endif

      return

      end subroutine controlOutputFlux
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINES: finalizeOutputFlux
!
! !INTERFACE:
!
      subroutine finalizeOutputFlux()
!
      implicit none
!
! !DESCRIPTION:
! Deallocates variables necessary to produce flux outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'finalizeOutputFlux called by ', procID

      deallocate (air_mass_flux_nc)

      if (pr_const_flux) then
         deallocate (flux_x_nc)
         deallocate (flux_y_nc)
         deallocate (flux_z_nc)
      endif

      if (pr_psf_flux) then
         deallocate (psf_flux_nc )
      endif

      if (iAmRootProc) call Nccl_Noerr (ncid_flux)

      return

      end subroutine finalizeOutputFlux
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocateOutputFlux
!   
! !INTERFACE:
!
      subroutine allocateOutputFlux ( )

!
      implicit none

!EOP
!-----------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'allocateOutputFlux called by ', procID

!... air mass flux
      Allocate (air_mass_flux_nc(i1:i2, ju1:j2, k1:k2, 3))
      air_mass_flux_nc = 0.0d0

      if (pr_const_flux) then
!... constituent E-W mass flux
         Allocate (flux_x_nc(i1:i2, ju1:j2, k1:k2, 1:flux_species_num))
         flux_x_nc = 0.0d0

!... constituent N-S mass flux
         Allocate (flux_y_nc(i1:i2, ju1:j2, k1:k2, 1:flux_species_num))
         flux_y_nc = 0.0d0

!... constituent vertical mass flux
         Allocate (flux_z_nc(i1:i2, ju1:j2, k1:k2, 1:flux_species_num))
         flux_z_nc = 0.0d0
      endif

      if (pr_psf_flux) then
         Allocate (psf_flux_nc (i1:i2, ju1:j2))
         psf_flux_nc = 0.0d0
      endif

      return

      end subroutine allocateOutputFlux
!EOC
!-----------------------------------------------------------------------------
!
! ROUTINE
!   bufferOutput_Flux
!
! DESCRIPTION
!   This routine buffers the flux data to the output file.  The Slaves send a
!   slice of the data to the Master, the Master writes it out, they then send
!   another slice and it is written out, and so on.
!
! ARGUMENTS
!   None
!
!-----------------------------------------------------------------------------

      subroutine bufferOutput_Flux (Advection)

      use GmiBuffer_mod, only : bufferOutput4d_Nomean

      implicit none
!
      type (t_Advection), intent(in) :: Advection
!
! ! DEFINED PARAMETERS:
      integer, parameter :: SG_FLUXm_NC     = 6450 
      integer, parameter :: SG_FLUXx_NC     = 6452 
      integer, parameter :: SG_FLUXy_NC     = 6454 
      integer, parameter :: SG_FLUXz_NC     = 6456 
      integer, parameter :: SG_FLUXx1_NC    = 6462 
      integer, parameter :: SG_FLUXy1_NC    = 6464 
      integer, parameter :: SG_FLUXz1_NC    = 6466 
      integer, parameter :: SG_PSF_FLUX_NC  = 6468 
      character (len=MAX_LENGTH_VAR_NAME), parameter :: CFlux_x_VNAM = 'const_flux_x'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: CFlux_y_VNAM = 'const_flux_y'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: CFlux_z_VNAM = 'const_flux_z'

      real*8, allocatable :: flux_x        (:,:,:,:) 
      real*8, allocatable :: flux_y        (:,:,:,:) 
      real*8, allocatable :: flux_z        (:,:,:,:) 

!     ----------------
!     Begin execution.
!     ----------------

      if (flux_species_num > 0 .and. pr_const_flux) then
!.sds.. X component

         allocate(flux_x(i1:i2,ju1:j2,k1:k2,flux_species_num))
         allocate(flux_y(i1:i2,ju1:j2,k1:k2,flux_species_num))
         allocate(flux_z(i1:i2,ju1:j2,k1:k2,flux_species_num))

         call Get_flux_x (Advection, flux_x)
         call Get_flux_y (Advection, flux_y)
         call Get_flux_z (Advection, flux_z)

         call bufferOutput4d_Nomean (CFlux_x_VNAM, k1, k2, commuWorld,  &
     &              SG_Fluxx1_NC, ncid_flux, flux_species_num, rnum_out_flux,  &
     &              flux_x, flux_x_nc, map1_u, numDomains, rootProc, procID, &
     &              i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)

!.sds.. Y component

         call bufferOutput4d_Nomean (CFlux_y_VNAM, k1, k2, commuWorld,  &
     &              SG_Fluxy1_NC, ncid_flux, flux_species_num, rnum_out_flux,  &
     &              flux_y, flux_y_nc, map1_u, numDomains, rootProc, procID, &
     &              i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)

!.sds.. Z component

         call bufferOutput4d_Nomean (CFlux_z_VNAM, k1, k2, commuWorld,  &
     &              SG_Fluxz1_NC, ncid_flux, flux_species_num, rnum_out_flux,  &
     &              flux_z, flux_z_nc, map1_u, numDomains, rootProc, procID, &
     &              i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)

         deallocate(flux_x)
         deallocate(flux_y)
         deallocate(flux_z)
      endif

      return

      end subroutine bufferOutput_Flux

!-----------------------------------------------------------------------------
!
! ROUTINE
!    Write_Netcdf_Flux
!
! DESCRIPTION
!   This routine collects flux data from the Slaves and sends it to the
!   Master for output.
!
! ARGUMENTS
!   None
!
!-----------------------------------------------------------------------------

      subroutine Write_Netcdf_Flux (nymd, nhms, gmi_sec)

      implicit none
! !INPUT PARAMETERS:
      integer, intent(in) :: nymd, nhms
      real*8 , intent(in) :: gmi_sec
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_FLUXm_NC     = 6450
      integer, parameter :: SG_FLUXx_NC     = 6452
      integer, parameter :: SG_FLUXy_NC     = 6454
      integer, parameter :: SG_FLUXz_NC     = 6456
      integer, parameter :: SG_FLUXx1_NC    = 6462 
      integer, parameter :: SG_FLUXy1_NC    = 6464
      integer, parameter :: SG_FLUXz1_NC    = 6466
      integer, parameter :: SG_PSF_FLUX_NC  = 6468
!
! !LOCAL VARIABLES:
      INTEGER :: Numb_Spec, ic
      integer :: cnt2d (2), cnt3d (3), cnt4d (4), cnt5d (5)
      integer :: strt2d(2), strt3d(3), strt4d(4), strt5d(5)

      integer :: hdr(NETCDF_HDF)
      real*8, allocatable :: air_mass_fluxGlob(:,:,:,:)
      real*8, allocatable :: flux_xGlob(:,:,:,:)
      real*8, allocatable :: flux_yGlob(:,:,:,:)
      real*8, allocatable :: flux_zGlob(:,:,:,:)
      real*8, allocatable :: psf_fluxGlob(:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Flux called by ', procID

      strt2d(:) = (/ 1, rnum_out_flux /)
      strt3d(:) = (/ 1, 1, rnum_out_flux /)

      cnt2d (:) = (/ NETCDF_HDF, 1 /)
      cnt3d (:) = (/ numLonGlob, numLatGlob, 1 /)

     if (iAmRootProc) then

        hdr(1) = Nint (gmi_sec)
        hdr(2) = nymd
        hdr(3) = nhms

        call Ncwr_2d_Int (hdr, ncid_flux, hdr_var_name, strt2d, cnt2d)

        cnt5d (:) = (/ numLonGlob, numLatGlob, numVert, 3, 1 /)
        strt5d(:) = (/ 1, 1, 1, 1, rnum_out_flux /)

!... do air mass flux
        allocate(air_mass_fluxGlob(i1_gl:i2_gl,ju1_gl:j2_gl,k1:k2,3))
     end if

      do ic=1,3
         call subDomain2Global(air_mass_fluxGlob(:,:,:,ic), &
     &           air_mass_flux_nc(:,:,:,ic), i1_gl, i2_gl, ju1_gl, j2_gl, i1, &
     &           i2, ju1, j2, k1, k2, rootProc, procID, map1_u, numDomains, &
     &           SG_Fluxm_NC, commuWorld)
      enddo

     if (iAmRootProc) then
        call Ncwr_5d (air_mass_fluxGlob, ncid_flux, 'amf', strt5d, cnt5d)
        deallocate(air_mass_fluxGlob)
     end if

!... do constituent fluxes
      if (pr_const_flux) then

         if (iAmRootProc) then
            allocate(flux_xGlob(i1_gl:i2_gl,ju1_gl:j2_gl,k1:k2,flux_species_num))
            allocate(flux_yGlob(i1_gl:i2_gl,ju1_gl:j2_gl,k1:k2,flux_species_num))
            allocate(flux_zGlob(i1_gl:i2_gl,ju1_gl:j2_gl,k1:k2,flux_species_num))
         end if

         do ic=1,flux_species_num

            call subDomain2Global(flux_xGlob(:,:,:,ic), flux_x_nc(:,:,:,ic),  &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, k1, k2, &
     &           rootProc, procID, map1_u, numDomains, SG_Fluxx_NC, commuWorld)

            call subDomain2Global(flux_yGlob(:,:,:,ic), flux_y_nc(:,:,:,ic),  &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, k1, k2, &
     &           rootProc, procID, map1_u, numDomains, SG_Fluxy_NC, commuWorld)

            call subDomain2Global(flux_zGlob(:,:,:,ic), flux_z_nc(:,:,:,ic),  &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, k1, k2, &
     &            rootProc, procID, map1_u, numDomains, SG_Fluxz_NC, commuWorld)

          enddo
          if (iAmRootProc) then
             cnt5d (:) = (/ numLonGlob, numLatGlob, numVert, flux_species_num, 1 /)
             strt5d(:) = (/ 1, 1, 1, 1, rnum_out_flux /)

!... and change flux_x_nc(:,:,:,:) to flux_x_nc(:,:,:,ic)
             call Ncwr_5d (flux_xGlob, ncid_flux, 'const_flux_x', strt5d, cnt5d)

             call Ncwr_5d (flux_yGlob, ncid_flux, 'const_flux_y', strt5d, cnt5d)

             call Ncwr_5d (flux_zGlob, ncid_flux, 'const_flux_z', strt5d, cnt5d)

             deallocate(flux_xGlob)
             deallocate(flux_yGlob)
             deallocate(flux_zGlob)
           end if

      endif


!... save consistent surface pressure
      if (pr_psf_flux) then

         if (iAmRootProc) then
            allocate(psf_fluxGlob(i1_gl:i2_gl, ju1_gl:j2_gl))
         end if

         call subDomain2Global (psf_fluxGlob, psf_flux_nc, i1_gl, i2_gl, &
     &           ju1_gl, j2_gl, i1, i2, ju1, j2, rootProc, procID, &
     &           map1_u, numDomains, SG_PSF_Flux_NC, commuWorld)

         if (iAmRootProc) then
            cnt3d (:) = (/ numLonGlob, numLatGlob, 1 /)
            strt3d(:) = (/ 1, 1, rnum_out_flux /)

            call Ncwr_3d (psf_fluxGlob, ncid_flux, 'Avg_sfc_press', strt3d, cnt3d)

            deallocate(psf_fluxGlob)
         end if
      endif

     if (iAmRootProc) call Ncdo_Sync (ncid_flux)

      return

      end subroutine Write_Netcdf_Flux
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Define_Netcdf_Out_Flux
!
! !INTERFACE:
!
      subroutine Define_Netcdf_Out_Flux  ()
!
! !USES:
      use m_ncGeneralOpsOutput, only: Define_Netcdf_Out_Gen
!
      implicit none
!
#     include "netcdf.inc"
!
! !DESCRIPTION:
! Makes the necessary definitions for the flux species netCDF output file.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_VAR_NAME) :: prsp1_dim_name
      integer :: ierr
      integer :: nchr1, nchr2, nstr, omode, pos1, varid
      integer :: chrd1(1), chrd2 (1), chrd3(1), strd (1)
      integer :: hdfd (1), lond (1), latd  (1)
      integer :: fluxd(1), prsd (1), prsp1d(1)
      integer :: recd (1), spcd  (1), threed(1), n3d
      integer :: spcd_2(1), spcd_3(1), spcd_4(1)
      integer :: var2(2), var3(3), var4(4), var5(5)
!
!EOC
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Define_Netcdf_Out_Flux called by ', procID

      nchr1 = MAX_LENGTH_SPECIES_NAME
      nchr2 = 50
      n3d   = 3
      nstr  =  1

!     ==========================
      call Define_Netcdf_Out_Gen  &
!     ==========================
     &  (ncid_flux, nhdf, numLonGlob, numLatGlob, numVert, hdfd, lond, latd, &
     &   prsd, hdf_dim_name, lon_dim_name, lat_dim_name, prs_dim_name)

!     ------------------
!     Define dimensions.
!     ------------------

      prsp1_dim_name = prs_dim_name
      pos1 = Len_Trim (prs_dim_name) + 1
      prsp1_dim_name(pos1:pos1+2) = 'p1'

      call NcDef_dimension(ncid_flux, prsp1_dim_name, numVert+1, prsp1d(1))
!                                        --------------

      call NcDef_dimension(ncid_flux, 'chr_dim1', nchr1, chrd1(1))
!                                        ----------
      call NcDef_dimension(ncid_flux, 'three', n3d, threed(1))
!                                        ----------

      if (flux_species_num > 0) then
         call NcDef_dimension (ncid_flux, 'flux_species', flux_species_num, fluxd(1))
      endif

      call NcDef_dimension(ncid_flux, rec_dim_name, NF_UNLIMITED, recd(1))
!                                       ------------

!     -----------------------------------------
!     Define variables and variable attributes.
!     -----------------------------------------
!
      call NcDef_variable (ncid_flux, 'am', NF_FLOAT, 1, prsd, varid)
!                                     ----
      call NcDef_var_attributes (ncid_flux, varid, 'long_name', 'Hybrid pressure term')
      call NcDef_var_attributes (ncid_flux, varid, 'units', 'unitless')

      call NcDef_variable (ncid_flux, 'bm', NF_FLOAT, 1, prsd, varid)
!                                     ----
      call NcDef_var_attributes (ncid_flux, varid, 'long_name', 'Hybrid sigma term ')
      call NcDef_var_attributes (ncid_flux, varid, 'units', 'unitless')

      call NcDef_variable (ncid_flux, 'ai', NF_FLOAT, 1, prsp1d, varid)
!                                     ----
      call NcDef_var_attributes (ncid_flux, varid, 'long_name', 'Hybrid pressure edge term')
      call NcDef_var_attributes (ncid_flux, varid, 'units', 'unitless')

      call NcDef_variable (ncid_flux, 'bi', NF_FLOAT, 1, prsp1d, varid)
!                                     ----
      call NcDef_var_attributes (ncid_flux, varid, 'long_name', 'Hybrid sigma edge term')
      call NcDef_var_attributes (ncid_flux, varid, 'units', 'unitless')

      call NcDef_variable (ncid_flux, 'pt', NF_FLOAT, 0, prsp1d, varid)
!                                     ----
      call NcDef_var_attributes (ncid_flux, varid, 'long_name', 'Hybrid coordinate constant for am/ai')
      call NcDef_var_attributes (ncid_flux, varid, 'units', 'unitless?')

!ccc Flux Species
!
      if (flux_species_num > 0) then
        call NcDef_variable (ncid_flux, 'flux_species', NF_FLOAT, 1, fluxd, varid)
!                                         --------------
        call NcDef_var_attributes (ncid_flux, varid, 'long_name', 'Flux Species index')
        call NcDef_var_attributes (ncid_flux, varid, 'units', 'unitless')
        call NcDef_var_attributes (ncid_flux, varid, 'coord_labels', 'flux_labels')
        call NcDef_var_attributes (ncid_flux, varid, 'selection_category', 'NULL')

        var2(:) = (/ chrd1(1), fluxd(1) /)
        call NcDef_variable (ncid_flux, 'flux_labels', NF_CHAR, 2, var2, varid)
!                                         -------------
        call NcDef_var_attributes (ncid_flux, varid, 'long_name', 'Flux Species name')
        call NcDef_var_attributes (ncid_flux, varid, 'units', 'unitless')
        call NcDef_var_attributes (ncid_flux, varid, 'selection_category', 'NULL')

!       ----------
!
        var2(:) = (/ lond(1), latd(1) /)
        call NcDef_variable (ncid_flux, 'mcor', NF_FLOAT, 2, var2, varid)
!                                    ------
        call NcDef_var_attributes (ncid_flux, varid, 'long_name', 'Grid box area')
        call NcDef_var_attributes (ncid_flux, varid, 'units', 'm^2')


! ================
!  Put Const Values
        if (pr_const_flux) then

          var5(:) = (/ lond(1), latd(1), prsd(1), fluxd(1), recd(1) /)
          call NcDef_variable (ncid_flux, 'const_flux_x', NF_FLOAT, 5, var5, varid)
!                                            ------------
          call NcDef_var_attributes (ncid_flux, varid, 'long_name','Constituent E-W flux')
          call NcDef_var_attributes (ncid_flux, varid, 'units', 'kg/fluxOutputFrequency')
          call NcDef_variable (ncid_flux, 'const_flux_y', NF_FLOAT, 5, var5, varid)
!                                            ------------
          call NcDef_var_attributes (ncid_flux, varid, 'long_name','Constituent N-S flux')
          call NcDef_var_attributes (ncid_flux, varid, 'units', 'kg/fluxOutputFrequency')
          call NcDef_variable (ncid_flux, 'const_flux_z', NF_FLOAT, 5, var5, varid)
!                                            ------------
          call NcDef_var_attributes (ncid_flux, varid, 'long_name','Constituent Vertical flux')
          call NcDef_var_attributes (ncid_flux, varid, 'units', 'kg/fluxOutputFrequency')
        endif

      endif
! ================
!    -----------------
      if (pr_psf_flux) then
        var3(:) = (/ lond(1), latd(1), recd(1) /)
        call NcDef_variable (ncid_flux, 'Avg_sfc_press', NF_FLOAT, 3, var3, varid)
!                                         -------------
        call NcDef_var_attributes (ncid_flux, varid, 'long_name', 'Averaged Surface pressure')
        call NcDef_var_attributes (ncid_flux, varid, 'units', 'hPa')
      endif
! ==================
      var5(:) = (/ lond(1), latd(1), prsd(1), threed(1), recd(1) /)
      call NcDef_variable (ncid_flux, 'amf', NF_FLOAT, 5, var5, varid)
!                                       ---
      call NcDef_var_attributes (ncid_flux, varid, 'long_name','Air Mass Flux')
      call NcDef_var_attributes (ncid_flux, varid, 'units', 'kg/fluxOutputFrequency')

      var2(:) = (/ hdfd(1), recd(1) /)
      call NcDef_variable (ncid_flux, hdr_var_name, NF_INT, 2, var2, varid)
!                                       ------------
      call NcDef_var_attributes (ncid_flux, varid, 'long_name', 'Header')
      call NcDef_var_attributes (ncid_flux, varid, 'units', 'gmi_sec, nymd, nhms')
!
!     -------------------------
!     Define global attributes.
!     -------------------------

      call NcDef_glob_attributes (ncid_flux, 'title',  &
     &     'Gmimod species flux file')

      call NcSetFill (ncid_flux, NF_NOFILL, omode)

      call NcEnd_def (ncid_flux)

      return

      end subroutine Define_Netcdf_Out_Flux
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Prep_Netcdf_Flux
!
! !INTERFACE:
!
      subroutine Prep_Netcdf_Flux (Advection)

      implicit none
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Advection), intent(inOut) :: Advection
!
! !DESCRIPTION:
!  Prepares the flux NetCDF output.
!
! !LOCAL VARIABLES:
      integer :: ic, icx
      integer :: ik, il, ikf
      real*8, allocatable :: air_mass_flux (:,:,:,:) 
      real*8, allocatable :: flux_x        (:,:,:,:) 
      real*8, allocatable :: flux_y        (:,:,:,:) 
      real*8, allocatable :: flux_z        (:,:,:,:) 
      real*8, allocatable :: psf_flux      (:,:)     
      integer :: count_flux
!
!EOC
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Prep_Netcdf_Flux called by ', procID

      call Get_count_flux (Advection, count_flux)
      
!       --------------------------------
!       Fill air_mass_flux_nc with current values.
!       --------------------------------
!... air mass fluxes
!... make units kg into box during integration time and flip

      allocate(air_mass_flux(i1:i2,ju1:j2,k1:k2,3))
      call Get_air_mass_flux (Advection, air_mass_flux)

      do il=1,3
         do ik=k1,k2
            ikf = k2-ik+1
            air_mass_flux_nc(:,:,ik,il) = air_mass_flux(:,:,ikf,il) *  &
     &                                  (PASPMB / GMI_G ) * mcor(:,:)
         enddo
      enddo

      deallocate(air_mass_flux)

!       ----------------------------------------------------
!       Flip constituent flux arrays in vertical for output.
!       ----------------------------------------------------
      if (pr_const_flux) then
         allocate(flux_x(i1:i2,ju1:j2,k1:k2,flux_species_num))
         allocate(flux_y(i1:i2,ju1:j2,k1:k2,flux_species_num))
         allocate(flux_z(i1:i2,ju1:j2,k1:k2,flux_species_num))

         call Get_flux_x (Advection, flux_x)
         call Get_flux_y (Advection, flux_y)
         call Get_flux_z (Advection, flux_z)

         flux_x(:,:,:,:) = flux_x(:,:,k2:k1:-1,:)
         flux_y(:,:,:,:) = flux_y(:,:,k2:k1:-1,:)
         flux_z(:,:,:,:) = flux_z(:,:,k2:k1:-1,:)

         call Set_flux_x (Advection, flux_x)
         call Set_flux_y (Advection, flux_y)
         call Set_flux_z (Advection, flux_z)

         deallocate(flux_x)
         deallocate(flux_y)
         deallocate(flux_z)
      endif

!       --------------------------------
!       Fill surface pressure concurent with flux arrays with current values.
!       --------------------------------
      if (pr_psf_flux) then
         allocate(psf_flux(i1:i2,ju1:j2))

         call Get_psf_flux (Advection, psf_flux)
         psf_flux_nc(:,:) = psf_flux(:,:) / count_flux

         deallocate(psf_flux)
      endif

      return

      end subroutine Prep_Netcdf_Flux
!EOC
!-----------------------------------------------------------------------------
!BOC
!
! ROUTINE
!   Write_Netcdf_Hdr_Flux
!
! DESCRIPTION
!   This routine creates some header information for the const NetCDF output
!   file and writes it out.
!
!-----------------------------------------------------------------------------

      subroutine Write_Netcdf_Hdr_Flux (gmiDomain, Chemistry, metFields)

      use m_ncGeneralOpsOutput, only : Write_Netcdf_Hdr_Gen

      implicit none

      type(t_gmiDomain), intent(in) :: gmiDomain
      type(t_Chemistry), intent(in) :: Chemistry
      type(t_metFields), intent(in) :: metFields

!     Variable declarations.
!     ----------------------

      character (len=MAX_LENGTH_SPECIES_NAME) :: spclab(flux_species_num)
      character (len=50) :: metnam(1)

      integer :: ic, icx

      integer :: cnt1d (1), cnt2d (2)
      integer :: strt1d(1), strt2d(2)

      real*8  :: prsdat(k1:k2)
      real*8  :: spcdat(flux_species_num)
      real*8 , allocatable :: londeg(:), latdeg(:)
      real*8 , allocatable :: mcorGlob(:,:)
      character (len=MAX_LENGTH_SPECIES_NAME) :: const_labels(numSpecies)
      real*8 , allocatable :: ai(:), bi(:)
      real*8 , allocatable :: am(:), bm(:)
      real*8 :: pt
!
!EOP
!------------------------------------------------------------------------------
!BOP

      if (pr_diag) Write (6,*) 'Write_Netcdf_Hdr_Flux called by ', procID

      allocate(latdeg(ju1_gl:j2_gl))
      call Get_latdeg(gmiDomain, latdeg)

      allocate(londeg(i1_gl:i2_gl))
      call Get_londeg(gmiDomain, londeg)

!     =========================
      call Write_Netcdf_Hdr_Gen  &
!     =========================
     &  (ncid_flux, latdeg, londeg, pr_diag, procID, i1_gl, i2_gl, ju1_gl, j2_gl,&
     &   hdf_dim_name, lat_dim_name, lon_dim_name)

      call Get_const_labels(Chemistry, const_labels)

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

      call Ncwr_1d (prsdat(kx1:kx2), ncid_flux, prs_dim_name, strt1d, cnt1d)
      call Ncwr_1d (am(kx1:kx2), ncid_flux, 'am', strt1d, cnt1d)
      call Ncwr_1d (bm(kx1:kx2), ncid_flux, 'bm', strt1d, cnt1d)

      cnt1d(1) = numvert + 1

      call Ncwr_1d (ai(kx1-1:kx2), ncid_flux, 'ai', strt1d, cnt1d)
      call Ncwr_1d (bi(kx1-1:kx2), ncid_flux, 'bi', strt1d, cnt1d)

      call Ncwr_Scal (pt, ncid_flux, 'pt')

!     ------------
!     Species map.
!     ------------

      icx = 0

      do ic = 1, numSpecies

        if (flux_species(ic) /= 0) then
          icx = icx + 1
          spcdat(icx) = ic
        endif

      enddo

      strt1d(1) = 1
      cnt1d (1) = flux_species_num

      if (flux_species_num > 0) then
        call Ncwr_1d (spcdat, ncid_flux, 'flux_species', strt1d, cnt1d)
      endif

!     ------------
!     Flux labels.
!     ------------

      icx = 0

      do ic = 1, numSpecies
        if (flux_species(ic) /= 0) then
          icx = icx + 1
          spclab(icx) = const_labels(ic)
        endif
      enddo

      strt2d(:) = (/  1, 1 /)
      cnt2d (:) = (/ MAX_LENGTH_SPECIES_NAME, flux_species_num /)

      if (flux_species_num > 0) then
         call Ncwr_2d_Char (spclab, ncid_flux, 'flux_labels', strt2d, cnt2d)
      endif

!     --------------
!     Grid box area.
!     --------------

      allocate(mcorGlob(i1_gl:i2_gl,ju1_gl:j2_gl))
      call Get_mcorGlob(gmiDomain, mcorGlob)

      strt2d(:) = (/ 1, 1 /)
      cnt2d (:) = (/ numLonGlob, numLatGlob /)

      call Ncwr_2d (mcorGlob, ncid_flux, 'mcor', strt2d, cnt2d)

      deallocate(mcorGlob)

      return

      end subroutine Write_Netcdf_Hdr_Flux
!EOC
!!-----------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: resetFlux
!!
!! !INTERFACE:
!!
!      subroutine resetFlux ( )
!!
!      implicit none
!!
!! !DESCRIPTION:
!! Resets the flux diagnostics.
!!
!! !LOCAL VARIABLES:
!      logical, save :: first = .true.
!!
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!      if (pr_diag) Write (6,*) 'resetFlux called by ', procID
!
!      if (first) then
!         first = .false.
!         do_flux_reset = .true.
!      endif
!
!!... reset accumulated quantities
!      if (pr_flux) then
!
!!... zero air mass flux
!         air_mass_flux = 0.0d0
!
!!... reset counter
!         count_flux = 0
!
!!... constituent component fluxes
!         if (pr_const_flux) then
!            flux_x = 0.0d0
!            flux_y = 0.0d0
!            flux_z = 0.0d0
!         endif
!
!!... constituent component fluxes
!         if (pr_psf_flux) then
!            psf_flux = 0.0d0
!         endif
!
!!... reset switch
!         do_flux_reset = .false.
!
!       endif
!
!      return
!
!      end subroutine resetFlux
!!EOC
!-----------------------------------------------------------------------------
   
      end module GmiControlFlux_mod
