!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiUpdateForcingBC_mod
!
      module GmiUpdateForcingBC_mod
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiTimeControl_mod, only : GmiSplitDateTime
!
      implicit none
!
      private
      public  :: updateForcingBC
!
#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
!
!EOP
!-------------------------------------------------------------------------
      contains
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: updateForcingBC
!
! !INTERFACE:
!
      subroutine updateForcingBC (concentration, nymd,       &
     &                 gmi_sec,  fbc_j1, fbc_j2, forc_bc_num, forc_bc_kmax,    &
     &                 forc_bc_kmin, forc_bc_opt, forc_bc_map, forc_bc_incrpyr,&
     &                 forc_bc_start_num, forc_bc_years, forc_bc_data, &
                       last_year, jlatmd, &
                       ju1_gl, j2_gl, i1, i2, ju1, j2, k1, k2, num_species)
!
      implicit none

#     include "gmi_forc_bc.h"
!
! !INPUT PARAMETERS:
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: ju1_gl, j2_gl
      integer, intent(in) :: num_species
      integer, intent(in) :: forc_bc_kmax, forc_bc_kmin
      integer, intent(in) :: forc_bc_num
      integer, intent(in) :: forc_bc_opt
      integer, intent(in) :: forc_bc_years
      integer, intent(in) :: fbc_j1, fbc_j2
      integer, intent(in) :: forc_bc_map (forc_bc_num)
      real*8 , intent(in) :: forc_bc_incrpyr
      integer, intent(in) :: nymd
      real*8 , intent(in) :: gmi_sec
      ! latitude of zone center in latitude direction (rad)
      real*8 , intent(in) :: forc_bc_data(:,:,:,:)
      integer, intent(in) :: jlatmd(j2_gl-ju1_gl+1)
!
! !INPUT/OUTPUT PARAMETERS:
      integer, intent(inOut) :: last_year
      integer, intent(inOut) :: forc_bc_start_num
      ! species concentration, known at zone centers (mixing ratio)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!
! !DESCRIPTION:
!  This routine updates the forcing boundary conditions.
!
! !LOCAL VARIABLES:
      integer :: idumday, month, year
      integer :: ic, ij, jj, jx
      real*8  :: fbc_emiss_fac, fbc_emission, rjx
! !REVISION HISTORY:
!  Original code
!EOP
!-------------------------------------------------------------------------
!BOC

      !====================
      if (forc_bc_opt == 1) then
      !====================
         call GmiSplitDateTime (nymd, year, month, idumday)

         do ic = 1, forc_bc_num
            do ij = ju1, j2
               jj = ij + 1 - ju1_gl
               concentration(forc_bc_map(ic))%pArray3D(:,ij,forc_bc_kmin:forc_bc_kmax) =  &
     &                      forc_bc_data(jlatmd(jj),month,forc_bc_start_num,ic)
            end do
         end do

      !=========================
      else if (forc_bc_opt == 2) then
      !=========================
         call GmiSplitDateTime (nymd, year, month, idumday)

         if (year /= last_year) then
            last_year         = year
            forc_bc_start_num = forc_bc_start_num + 1
         end if

         do ic = 1, forc_bc_num
            do ij = ju1, j2
               jj = ij + 1 - ju1_gl
               concentration(forc_bc_map(ic))%pArray3D(:,ij,forc_bc_kmin:forc_bc_kmax) =  &
     &                       forc_bc_data(jlatmd(jj),month,forc_bc_start_num,ic)
            end do
         end do

      !=========================
      else if (forc_bc_opt == 3) then
      !=========================
         fbc_emiss_fac = (gmi_sec / SECPYR) + (forc_bc_start_num - 1)
         fbc_emission  = fbc_emiss_fac * forc_bc_incrpyr * PPT_FAC

         do ic = 1, forc_bc_num
            do ij = ju1, j2
               jj = ij + 1 - ju1_gl
               if ((ij >= fbc_j1) .and. (ij <= fbc_j2)) then
                  concentration(forc_bc_map(ic))%pArray3D(:,ij,forc_bc_kmin:forc_bc_kmax)  &
     &                       = fbc_emission
               end if
            end do
         end do
      end if

      return

      end subroutine updateForcingBC
!EOC
!------------------------------------------------------------------------------

      end module GmiUpdateForcingBC_mod
