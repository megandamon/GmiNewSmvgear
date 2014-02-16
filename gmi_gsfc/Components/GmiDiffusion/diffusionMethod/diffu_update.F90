!=============================================================================
!
! $Id: diffu_update.F90,v 1.4 2006-09-21 17:01:02 kouatch Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   diffu_update.F
!
! ROUTINES
!   Update_Diffu
!   Update_PBL_Mixing
!   Adjust_Kzz
!
!=============================================================================

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update_Diffu
!
! DESCRIPTION
!   This routine updates diffusion.
!
! ARGUMENTS
!   tdt     : model time step (s)
!   vert_diffu_coef : scalar vertical diffusion coefficient (m^2/s)
!   pbl     : planetary boundary layer depth (m)
!   tropp   : tropopause pressure (mb)
!   kzz     : array of vertical diffusion coefficients (m^2/s)
!   press3c : atmospheric pressure at the center of each grid box (mb)
!   press3e : atmospheric pressure at the edge   of each grid box (mb)
!   pctm1   : surface pressure field at t1, known at zone centers (mb)
!   const   : species concentration (mixing ratio)
!
!-----------------------------------------------------------------------------

      subroutine Update_Diffu  &
     &  (tdt, vert_diffu_coef, pbl, tropp, kzz, press3c, press3e,  &
     &   pctm1, concentration, &
     &   pr_diag, loc_proc, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in   ) :: pr_diag
      integer, intent(in   ) :: loc_proc
      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in   ) :: num_species
      real*8 , intent(in   ) :: tdt
      real*8 , intent(in   ) :: vert_diffu_coef
      real*8 , intent(in   ) :: pbl    (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: tropp  (i1:i2, ju1:j2)
      real*8 , intent(inout) :: kzz    (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in   ) :: press3c(ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in   ) :: press3e(ilo:ihi, julo:jhi, k1-1:k2)
      real*8 , intent(in   ) :: pctm1  (ilo:ihi, julo:jhi)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Update_Diffu called by ', loc_proc
      end if

!     ===============
      call Adjust_Kzz  &
!     ===============
     &  (press3c, tropp, kzz, &
     &   pr_diag, loc_proc,  &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)

!     ==================
      call Do_Vert_Diffu  &
!     ==================
     &  (tdt, vert_diffu_coef, pbl, kzz, press3c, press3e,  &
     &   pctm1, concentration, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species)

      return

      end subroutine Update_Diffu 

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update_PBL_Mixing
!
! DESCRIPTION
!   This routine updates diffusion through a boundary layer mixing model.
!
! ARGUMENTS
!   tdt            : model time step (s)
!   pbl_mixing_tau : time scale for pbl mixing (s)
!   pbl            : planetary boundary layer depth (m)
!   grid_height    : height of each grid box (m)
!   mass           : array of air mass in each grid box (kg)
!   const          : species concentration (mixing ratio)
!
!-----------------------------------------------------------------------------

      subroutine Update_PBL_Mixing  &
     &  (tdt, pbl_mixing_tau, pbl, grid_height, mass, concentration, &
     &   pr_diag, loc_proc, &
     &   i1, i2, ju1, j2, k1, k2, num_species)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in   ) :: pr_diag
      integer, intent(in   ) :: loc_proc
      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in   ) :: num_species

      real*8 , intent(in   ) :: tdt
      real*8 , intent(in   ) :: pbl_mixing_tau
      real*8 , intent(in   ) :: pbl        (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: grid_height(i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in   ) :: mass       (i1:i2, ju1:j2, k1:k2)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ic, il, ij, ik, kpbl
      real*8  :: frac_pbl, height, mass_pbl, mixing_factor, well_mixed

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Update_PBL_Mixing called by ', loc_proc
      end if

!     write(*,*) 'Bigyani tdt and pbl_mixing_tau =', tdt, pbl_mixing_tau
      mixing_factor = Exp(-tdt / pbl_mixing_tau)

!     ---------------------------------------
!     Loop over all longitudes and latitudes.
!     ---------------------------------------

      do il = i1, i2
      do ij = ju1, j2

        mass_pbl = 0.0d0
        height = 0.0d0

!       -----------------------------------
!       Find the total mass within the pbl.
!       -----------------------------------

        kloop: do ik = k1, k2
          height = height + grid_height(il,ij,ik)
          if (pbl(il, ij) >= height) then
            mass_pbl = mass_pbl + mass(il, ij, ik)
          else

!           ---------------------------------------------------------------------
!           This grid box is fractionally in the PBL. Save part of its mass also,
!           save the k index, and exit the loop.
!           ---------------------------------------------------------------------

            frac_pbl = (pbl(il,ij) -  &
     &       (height - grid_height(il,ij,ik))) / grid_height(il,ij,ik)

            mass_pbl = mass_pbl + mass(il,ij,ik) * frac_pbl

            kpbl = ik

            exit kloop

          end if
        end do kloop

!       --------------------------------------------------------
!       Now if a real PBL was found (kpbl > k1) loop over the
!       species, calculate a well mixed
!       mixing ratio, and apply it to the grids within the pbl
!       and then to the grid box
!       which is fractionally within the pbl.
!       ---------------------------------------------------------

        if (kpbl > k1) then
        do ic = 1, num_species
          well_mixed = (Sum(concentration(ic)%pArray3D(il,ij,k1:kpbl-1) *  &
     &       mass(il,ij,k1:kpbl-1)) +  &
     &       frac_pbl * concentration(ic)%pArray3D(il,ij,kpbl) * mass(il,ij,kpbl)) /  &
     &       mass_pbl

          concentration(ic)%pArray3D(il,ij,1:kpbl-1) =  &
     &      concentration(ic)%pArray3D(il,ij,k1:kpbl-1) * mixing_factor +  &
     &      well_mixed * (1.0d0 - mixing_factor)

          concentration(ic)%pArray3D(il,ij,kpbl) = &
     &       concentration(ic)%pArray3D(il,ij,kpbl) + frac_pbl *  &
     &      (concentration(ic)%pArray3D(il,ij,kpbl)*(mixing_factor - 1.0d0) +  &
     &       well_mixed * (1.0d0 - mixing_factor))

        end do
        end if

      end do
      end do

      return

      end subroutine Update_PBL_Mixing

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Adjust_Kzz
!
! DESCRIPTION
!   This routine makes the kzz (vertical diffusion) zero above the
!   tropopause.
!
! ARGUMENTS
!   press3c : atmospheric pressure at the center of each grid box (mb)
!   tropp   : tropopause pressure (mb)
!   kzz     : array of vertical diffusion coefficients (m^2/s)
!
!-----------------------------------------------------------------------------

      subroutine Adjust_Kzz  &
     &  (press3c, tropp, kzz, &
     &   pr_diag, loc_proc,    &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in   ) :: pr_diag
      integer, intent(in   ) :: loc_proc
      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      real*8 , intent(in   ) :: press3c(ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in   ) :: tropp  (i1:i2, ju1:j2)
      real*8 , intent(inout) :: kzz    (i1:i2, ju1:j2, k1:k2)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Adjust_Kzz called by ', loc_proc
      end if

      if (Maxval (kzz) > 0.0d0) then
        where (press3c(i1:i2,ju1:j2,:) <  &
     &         (Spread (tropp(:,:), 3, k2) + 100.0d0))
          kzz(i1:i2,ju1:j2,:) = 0.0d0
        end where
      end if

      return

      end subroutine Adjust_Kzz
