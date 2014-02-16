!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlMetInput_mod
!
      module GmiControlMetInput_mod
!
! !USES:
      use ESMF_Mod, only : ESMF_Config, ESMF_MAXSTR, ESMF_ConfigGetAttribute
      use GmiESMFrcFileReading_mod, only : rcEsmfReadLogical
      use GmiPrintError_mod, only : GmiPrintError
      use GmiAdvectionMethod_mod, only : t_Advection, Set_num_adv_time_steps, &
     &       Get_do_var_adv_tstp, Get_advec_consrv_opt, Get_pmet2_opt, &
     &       Get_press_fix_opt, Get_j1p, Get_j2p
      use GmiSurfaceAlbedo_mod, only : setSurfaceAlbedo, setSurfaceAlbedoUV
      use GmiChemistryMethod_mod, only : Get_sasdir_data, Get_sasdif_data,     &
     &       t_Chemistry, Get_saldir_data, Get_saldif_data, Get_uvalbedo_data, &
     &       Get_uvalbedo_opt, Get_sfalbedo_opt

      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Get_pr_diag,    &
     &       Get_rd_restart

      use GmiTimeControl_mod, only : t_GmiClock, Get_curGmiDate

      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_procID
      use GmiGrid_mod, only : t_gmiGrid
      use GmiGrid_mod, only : Get_i1, Get_i2, Get_ju1, Get_jv1
      use GmiGrid_mod, only : Get_j2, Get_k1, Get_k2
      use GmiMetFieldsControl_mod, only : t_metFields, Get_surf_alb_uv, &
     &       Update_Met1, Control_Met1_File_Input, Get_sasdir, Get_sasdif, &
     &       Update_Met2, Control_Met2_File_Input, Get_saldir, Get_saldif, &
     &       Set_sasdir, Set_sasdif, Set_saldir, Set_saldif, Set_surf_alb_uv, &
     &       Get_met_opt

      implicit none

      private
      public  :: Control_Met1_Input
      public  :: Control_Met2_Input, Set_Restart_Partial

#     include "GmiParameters.h"
#     include "gmi_time_constants.h"
!
!EOP
!=============================================================================
      contains
!=============================================================================
!BOP
!
! !IROUTINE: Control_Met1_Input
!
! !INTERFACE:
!
      subroutine Control_Met1_Input (metFields,  Advection, &
     &                 gmiDomain, Diagnostics, first_tstp, new_met_rec)

      implicit none
!
! !INPUT PARAMETERS:
      logical :: first_tstp   ! first time step?
      logical :: new_met_rec  ! new met record?
      type (t_gmiDomain  ), intent(in) :: gmiDomain
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Advection), intent(inOut) :: Advection
      type (t_metFields), intent(inOut) :: metFields
!
! !DESCRIPTION:
! Controls and updates the met1 data (ps, u, v, & kel).
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      logical :: out_of_data
      integer :: procID, j1p, j2p
      logical :: pr_diag, do_var_adv_tstp, rd_restart
      integer :: advec_consrv_opt, pmet2_opt, press_fix_opt
      integer :: num_adv_time_steps
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID (gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write (6,*) 'Control_Met1_Input called by ', procID

      out_of_data = .false.

      if ((.not. first_tstp) .and. new_met_rec) then

         call Control_Met1_File_Input (metFields, out_of_data)

        if (out_of_data) then
          err_msg = 'No more met data in Control_Met1_Input.'
          call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
        end if

      end if

      call Get_j1p(Advection, j1p)
      call Get_j2p(Advection, j2p)
      call Get_pmet2_opt(Advection, pmet2_opt)
      call Get_press_fix_opt(Advection, press_fix_opt)
      call Get_do_var_adv_tstp(Advection, do_var_adv_tstp)
      call Get_advec_consrv_opt(Advection, advec_consrv_opt)

      call Get_rd_restart (Diagnostics, rd_restart)

      call Update_Met1 (metFields, gmiDomain, new_met_rec, &
     &          j1p, j2p, rd_restart, num_adv_time_steps, do_var_adv_tstp, &
     &          advec_consrv_opt, pmet2_opt, press_fix_opt)

      call Set_num_adv_time_steps(Advection, num_adv_time_steps)

      return

      end subroutine Control_Met1_Input
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Control_Met2_Input
!
! !INTERFACE:
!
      subroutine Control_Met2_Input (metFields, Chemistry, gmiDomain, gmiGrid, &
     &                  Diagnostics, gmiClock, first_tstp, new_met_rec)

      implicit none
!
! !INPUT PARAMETERS:
      logical :: first_tstp    ! first time step?
      logical :: new_met_rec   ! new met record?
      type (t_gmiGrid    ), intent(in) :: gmiGrid
      type (t_gmiClock   ), intent(in) :: gmiClock 
      type (t_gmiDomain  ), intent(in) :: gmiDomain
      type (t_Chemistry  ), intent(in) :: Chemistry
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT PARAMETERS:
      type (t_metFields), intent(inOut) :: metFields
!
! !DESCRIPTION:
! Controls and updates the met2 data (everything except ps, u, v, & kel).
!
! !LOCAL VARIABLES:
      character(len=ESMF_MAXSTR) :: err_msg, IAm
      integer :: rc, STATUS, met_opt
      integer :: procID, sfalbedo_opt, uvalbedo_opt, nymd
      integer :: i1, i2, ju1, j2, k1, k2
      logical :: pr_diag, do_wetdep
      logical :: out_of_data
      real*8, allocatable :: sasdir_data(:,:,:), sasdif_data(:,:,:)
      real*8, allocatable :: saldir_data(:,:,:), saldif_data(:,:,:)
      real*8, allocatable :: uvalbedo_data(:,:,:)
      real*8, allocatable :: sasdir(:,:), sasdif(:,:)
      real*8, allocatable :: saldir(:,:), saldif(:,:)
      real*8, allocatable :: surf_alb_uv(:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      IAm = "Control_Met2_Input"

      call Get_procID      (gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write (6,*) IAm, ' called by ', procID

      call Get_i1    (gmiGrid, i1    )
      call Get_i2    (gmiGrid, i2    )
      call Get_ju1   (gmiGrid, ju1   )
      call Get_j2    (gmiGrid, j2    )

      call Get_met_opt(metFields, met_opt)

      call Get_curGmiDate(gmiClock, nymd)

      out_of_data = .false.

      if ( met_opt == 3) then
        if ((.not. first_tstp) .and. new_met_rec) then

!       ============================
          call Control_Met2_File_Input  &
!       ============================
     &      (metFields, gmiDomain, new_met_rec, out_of_data)

          if (out_of_data) then
            err_msg = 'No more met data in Control_Met2_Input.'
            call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
          end if

        end if

!     ================
        call Update_Met2 (metFields, new_met_rec)
!     ================
      end if

      call Get_sfalbedo_opt(Chemistry, sfalbedo_opt)

      if ((sfalbedo_opt == 1) .or. (sfalbedo_opt == 2)) then
         allocate(sasdir_data(i1:i2, ju1:j2, MONTHS_PER_YEAR))
         call Get_sasdir_data(Chemistry, sasdir_data)

         allocate(sasdif_data(i1:i2, ju1:j2, MONTHS_PER_YEAR))
         call Get_sasdif_data(Chemistry, sasdif_data)

         allocate(saldir_data(i1:i2, ju1:j2, MONTHS_PER_YEAR))
         call Get_saldir_data(Chemistry, saldir_data)

         allocate(saldif_data(i1:i2, ju1:j2, MONTHS_PER_YEAR))
         call Get_saldif_data(Chemistry, saldif_data)

         allocate(sasdir(i1:i2, ju1:j2))
         allocate(sasdif(i1:i2, ju1:j2))
         allocate(saldir(i1:i2, ju1:j2))
         allocate(saldif(i1:i2, ju1:j2))

         call Get_sasdir(metFields, sasdir)
         call Get_sasdif(metFields, sasdif)
         call Get_saldir(metFields, saldir)
         call Get_saldif(metFields, saldif)

         call setSurfaceAlbedo(pr_diag, procID, sfalbedo_opt, nymd,         &
     &           sasdir, sasdif, saldir, saldif, sasdir_data, sasdif_data,  &
     &           saldir_data, saldif_data, i1, i2, ju1, j2)

         call Set_sasdir(metFields, sasdir)
         call Set_sasdif(metFields, sasdif)
         call Set_saldir(metFields, saldir)
         call Set_saldif(metFields, saldif)

         deallocate(sasdir_data)
         deallocate(sasdif_data)
         deallocate(saldir_data)
         deallocate(saldif_data)
      end if

      call Get_uvalbedo_opt(Chemistry, uvalbedo_opt)

      if ((uvalbedo_opt == 1) .or. (uvalbedo_opt == 2)) then
         allocate(uvalbedo_data(i1:i2, ju1:j2, MONTHS_PER_YEAR))
         call Get_uvalbedo_data(Chemistry, uvalbedo_data)

         allocate(surf_alb_uv(i1:i2, ju1:j2))
         call Get_surf_alb_uv(metFields, surf_alb_uv)

         call setSurfaceAlbedoUV (pr_diag, procID, nymd, uvalbedo_opt,         &
     &                uvalbedo_data, surf_alb_uv, i1, i2, ju1, j2)
      
         call Set_surf_alb_uv(metFields, surf_alb_uv)

         deallocate(uvalbedo_data)
      end if

      return

      end subroutine Control_Met2_Input
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Restart_Partial
!
! !INTERFACE:
!
      subroutine Set_Restart_Partial (metFields, Diagnostics, gmiGrid, gmiDomain)
!
! !USES:
      use GmiWrapMaster_mod, only : wrapMaster_2d
      use GmiCheckRange_mod, only : checkRange2d
      use m_netcdf_io_open , only : Ncop_Rd
      use m_netcdf_io_close , only : Nccl
      use m_netcdf_io_read , only : Ncrd_3d
      use m_netcdf_io_checks , only : Ncdoes_Var_Exist
      use GmiGrid_mod, only : t_gmiGrid, Get_gmi_nborder
      use GmiGrid_mod, only : Get_k1, Get_k2
      use GmiGrid_mod, only : Get_i1_gl, Get_ju1_gl
      use GmiGrid_mod, only : Get_i2_gl, Get_j2_gl
      use GmiGrid_mod, only : Get_numSpecies
      use GmiGrid_mod, only : Get_ilo_gl, Get_ihi_gl, Get_julo_gl
      use GmiGrid_mod, only : Get_jvlo_gl, Get_jhi_gl
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_procID
      use GmiDiagnosticsMethod_mod  , only : t_Diagnostics, Get_pr_diag, &
     &       Get_restart_inrec, Get_rd_restart, Get_restart_infile_name
      use GmiMetFieldsControl_mod, only : t_metFields, Set_pctm2Glob

      implicit none
!
! !INPUT PARAMETERS:
      type (t_gmiGrid    ), intent(in) :: gmiGrid
      type (t_gmiDomain  ), intent(in) :: gmiDomain
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_metFields), intent(inOut) :: metFields
!
! !DESCRIPTION:
! Reads in the met1 variable pctm2 from the restart file.
!
! !DECLARED PARAMETERS:
      character (len=MAX_LENGTH_VAR_NAME), parameter :: PCTM2_VNAM   = 'pctm2'
!
! !LOCAL VARIABLES:
      logical :: pr_diag, rd_restart
      integer :: i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2, procID
      integer :: ilo_gl, ihi_gl, julo_gl, jhi_gl, gmi_nborder
      integer :: il, ij, ik, ilong, ilat, ivert
      integer :: inc, jnc
      integer :: ncid, restart_inrec
      integer :: cnt3d (3), cnt4d (4)
      integer :: strt3d(3), strt4d(4)
      real*8, allocatable :: pctm2_3d_nc (:, :, :)
      real*8, allocatable :: pctm2Glob (:,:)
      character (len=MAX_LENGTH_FILE_NAME) :: restart_infile_name
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_rd_restart(Diagnostics, rd_restart)

      if (rd_restart) then
         call Get_procID (gmiDomain, procID)
         call Get_pr_diag(Diagnostics, pr_diag)

         if (pr_diag) Write (6,*) 'Set_Restart_Partial called by ', procID

         call Get_k1    (gmiGrid, k1    )
         call Get_k2    (gmiGrid, k2    )
         call Get_i1_gl (gmiGrid, i1_gl )
         call Get_i2_gl (gmiGrid, i2_gl )
         call Get_ju1_gl(gmiGrid, ju1_gl)
         call Get_j2_gl (gmiGrid, j2_gl )
         call Get_ilo_gl (gmiGrid, ilo_gl )
         call Get_ihi_gl (gmiGrid, ihi_gl )
         call Get_julo_gl(gmiGrid, julo_gl)
         call Get_jhi_gl (gmiGrid, jhi_gl )
         call Get_gmi_nborder (gmiGrid, gmi_nborder )

         allocate(pctm2Glob(ilo_gl:ihi_gl,julo_gl:jhi_gl))
         pctm2Glob = 0.0d0

         ilong = i2_gl - i1_gl  + 1
         ilat  = j2_gl - ju1_gl + 1
         ivert = k2    - k1     + 1

         call Get_restart_inrec(Diagnostics, restart_inrec)
         call Get_restart_infile_name(Diagnostics, restart_infile_name)

         call Ncop_Rd (ncid, restart_infile_name)

         if (Ncdoes_Var_Exist (ncid, PCTM2_VNAM)) then

            ! -------------------------------------------------------
            ! Get pctm2; note that old restart files may not have it.
            ! -------------------------------------------------------

            allocate(pctm2_3d_nc (ilong, ilat, 1))

            strt3d(:) = (/ 1, 1, restart_inrec /)
            cnt3d (:) = (/ ilong, ilat, 1 /)

            call Ncrd_3d (pctm2_3d_nc, ncid, PCTM2_VNAM, strt3d, cnt3d)

            do ij = ju1_gl, j2_gl
               jnc = ij - ju1_gl + 1
               do il = i1_gl, i2_gl
                  inc = il - i1_gl + 1
                  pctm2Glob(il,ij) = pctm2_3d_nc(inc,jnc,1)
               end do
            end do

            call wrapMaster_2d (pctm2Glob, i1_gl, i2_gl, ju1_gl, j2_gl, ilo_gl,&
     &           ihi_gl, julo_gl, jhi_gl, gmi_nborder)

            call checkRange2d ('pctm2', procID, i1_gl, ihi_gl, ju1_gl, jhi_gl, &
     &           pctm2Glob, 0.0d0, 5000.0d0)

            deallocate(pctm2_3d_nc)
         else

            pctm2Glob(:,:) = -999.0d0

         end if

         call Set_pctm2Glob(metFields, pctm2Glob)
         deallocate(pctm2Glob)
      end if

      return

      end subroutine Set_Restart_Partial
!EOC
!------------------------------------------------------------------------------
      end module GmiControlMetInput_mod
