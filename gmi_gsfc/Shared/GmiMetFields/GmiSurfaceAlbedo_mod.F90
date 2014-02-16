!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE:
!
      module GmiSurfaceAlbedo_mod
!
! !USES:
      use GmiTimeControl_mod, only : GmiSplitDateTime
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: setSurfaceAlbedo, setSurfaceAlbedoUV
!
! !AUTHOR:
! Jules Kouatchou, Jules.Kouatchou-1@nasa.gov
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setSurfaceAlbedo
!
! !INTERFACE:
!
      subroutine setSurfaceAlbedo(pr_diag, procID, sfalbedo_opt, nymd,         &
     &    sasdir, sasdif, saldir, saldif, sasdir_data, sasdif_data, &
     &    saldir_data, saldif_data, i1, i2, ju1, j2)
!
      implicit none
!
! !INPUT VARIABLES:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: i1, i2, ju1, j2
      integer, intent(in) :: procID, sfalbedo_opt, nymd
      real*8 , intent(in) :: sasdir_data(i1:,ju1:,1:)
      real*8 , intent(in) :: sasdif_data(i1:,ju1:,1:)
      real*8 , intent(in) :: saldir_data(i1:,ju1:,1:)
      real*8 , intent(in) :: saldif_data(i1:,ju1:,1:)
!
! !OUTPUT VARIABLES:
                              ! surface albedo for direct  light (uv/vis) (0-1)
      real*8 , intent(out) :: sasdir(i1:i2,ju1:j2) 
                              ! surface albedo for diffuse light (uv/vis) (0-1)
      real*8 , intent(out) :: sasdif(i1:i2,ju1:j2)
                              ! surface albedo for direct  light (nr IR)  (0-1)
      real*8 , intent(out) :: saldir(i1:i2,ju1:j2)
                              ! surface albedo for diffuse light (nr IR)  (0-1)
      real*8 , intent(out) :: saldif(i1:i2,ju1:j2)
!
! !DESCRIPTION:
! Sets four types of surface albedo for radiation routines.
!
! !LOCAL VARIABLES:
      logical, save :: first = .true.
      integer, save :: month_save = -999
      integer :: idumday, idumyear
      integer :: month
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'setSurfaceAlbedo called by ', procID

      if (first) then

         first = .false.

         if (sfalbedo_opt == 1) then
            sasdir(:,:) = sasdir_data(:,:,1)
            sasdif(:,:) = sasdif_data(:,:,1)
            saldir(:,:) = saldir_data(:,:,1)
            saldif(:,:) = saldif_data(:,:,1)
         end if

      end if


      if (sfalbedo_opt == 2) then

         call GmiSplitDateTime (nymd, idumyear, month, idumday)

         if (month /= month_save) then
            sasdir(:,:) = sasdir_data(:,:,month)
            sasdif(:,:) = sasdif_data(:,:,month)
            saldir(:,:) = saldir_data(:,:,month)
            saldif(:,:) = saldif_data(:,:,month)

            month_save = month
         end if
      end if

      return

      end subroutine setSurfaceAlbedo
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setSurfaceAlbedoUV
!
! !INTERFACE:
!
      subroutine setSurfaceAlbedoUV (pr_diag, procID, nymd, uvalbedo_opt,         &
     &                uvalbedo_data, surf_alb_uv,i1, i2, ju1, j2)
!
! !USES:
      use GmiTimeControl_mod, only : GmiSplitDateTime
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: i1, i2, ju1, j2
      integer, intent(in) :: uvalbedo_opt, procID
      integer, intent(in) :: nymd
      real*8 , intent(in) :: uvalbedo_data(i1:,ju1:,1:)
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: surf_alb_uv(i1:i2,ju1:j2)
!
! !DESCRIPTION:
! Sets the bulk surface albedo.
!
! !LOCAL VARIABLES:
      logical, save :: first = .true.
      integer :: idumday, idumyear
      integer :: month
      integer, save :: month_save = -999
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'setSurfaceAlbedoUV called by ', procID

      if (first) then
         first = .false.

         if (uvalbedo_opt == 1) then
            surf_alb_uv(:,:) = uvalbedo_data(:,:,1)
         end if
      end if

      if (uvalbedo_opt == 2) then
         call GmiSplitDateTime (nymd, idumyear, month, idumday)

         if (month /= month_save) then
            surf_alb_uv(:,:) = uvalbedo_data(:,:,month)
            month_save = month
         end if
      end if

      return

      end subroutine setSurfaceAlbedoUV
!EOC
!------------------------------------------------------------------------------
      end module GmiSurfaceAlbedo_mod
