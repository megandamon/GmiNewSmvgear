
!=============================================================================
!
! $Id: convec_dao2.F90,v 1.12 2013-07-31 14:23:48 ssteenro Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Shian-Jiann Lin, DAO)
!   jrt@llnl.gov
!
! FILE
!   convec_dao2.F
!
! ROUTINES
!   Do_Convec_Dao2
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Convec_Dao2
!
! DESCRIPTION
!   This is the cumulus transport routine. It takes wet convection into
!   account.
!
!             TOP
!              ^
!              |
!              |
!    ------ cldmas(k) -----------   ap(k), bp(k)
!
!        q(k), dtrn(k)   (layer-k)
!
!    ----- cldmas(k+1) ----------- ap(k+1), bp(k+1)
!              ^
!              |
!              |
!           SURFACE
!
!   Definition of the hybrid sigma-p coordinate:
!
!     Pressure at layer edges is defined as follows:
!
!       p(i,j,k) = (ap(k) * pt) + (bp(k) * ps(i,j))
!
!       bp(k1)   = 0.0d0
!       bp(k2+1) = 1.0d0
!
!     Pressure at model top is:
!
!       ptop     = ap(k1) * pt
!
!     For pure sigma system set:
!
!       pt       = ptop
!       ap(k)    = 1.0d0 for all k
!       bp(k)    = sige(k) (i.e., sigma at edges)
!       ps       = psfc - ptop, where psfc is the "true" surface pressure
!
! ARGUMENTS
!   tdt    : model time step (s)
!   pt     : pressure = (ai * pt) + (bi * psf) (mb)
!   ai     : pressure = (ai * pt) + (bi * psf), ai at zone interface
!   bi     : pressure = (ai * pt) + (bi * psf), bi at zone interface
!   cldmas : cloud mass flux,  should be positive definite (kg/m^2*s)
!   dtrn   : detrainment rate, should be positive definite
!            (DAO:kg/m^2*s, NCAR:s^-1)
!            (set dtrn to 0.0d0 if data not available)
!   ps     : psurf - ptop, where psurf is the "true" surface pressure (mb)
!   qq     : species concentration (mixing ratio)
!
!-----------------------------------------------------------------------------

      subroutine Do_Convec_Dao2  &
     &  (tdt, cldmas, dtrn, concentration, bmass, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, numSpecies)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      use GmiSpcConcentrationMethod_mod, only : isFixedConcentration

      implicit none

#     include "gmi_phys_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in   ) :: numSpecies
      real*8 , intent(in   ) :: tdt
      real*8 , intent(in   ) :: bmass (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in   ) :: cldmas(i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in   ) :: dtrn  (i1:i2, ju1:j2, k1:k2)
      type (t_GmiArrayBundle), intent(inout) :: concentration(numSpecies)

!     -----------------------
!     Parameter declarations.
!     -----------------------

!     -----------------------------------------------------
!     NSPLIT : internal time step for convective transport;
!              note that 300s is rather conservative;
!              450s works well with DAO's data
!     -----------------------------------------------------

      integer, parameter :: NSPLIT = 300

      real*8,  parameter :: TINY   = 1.0d-14

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij, ik, ic
      integer :: istep
      integer :: jn, js
!c    integer :: jnp, jump
      integer :: k2m1, ktop
      integer :: ns

      real*8  :: cmout
!      real*8  :: dapx, dbkx
      real*8  :: dq
      real*8  :: entrn
      real*8  :: fg
      real*8  :: qcc
      real*8  :: rns
      real*8  :: sdt
      real*8  :: tprod

      real*8  :: ap(k1:k2+1)
      real*8  :: bp(k1:k2+1)

      real*8  :: qb(i1:i2, ju1:j2)
      real*8  :: qc(i1:i2, ju1:j2)
      real*8  :: wk(i1:i2, ju1:j2)

      real*8  :: temp_cldmas(i1:i2, ju1:j2, k1:k2)
      real*8  :: temp_dtrn  (i1:i2, ju1:j2, k1:k2)

!     ----------------
!     Begin execution.
!     ----------------

      qb(:,:) = 0.0d0
      qc(:,:) = 0.0d0
      wk(:,:) = 0.0d0

!     ------------------------------
!     Flip variables in k dimension.
!     ------------------------------

      temp_cldmas(:,:,k1:k2) = cldmas(:,:,k2:k1:-1)
      temp_dtrn  (:,:,k1:k2) = dtrn  (:,:,k2:k1:-1)

      do ic = 1, numSpecies
         concentration(ic)%pArray3D(:,:,k1:k2) = concentration(ic)%pArray3D(:,:,k2:k1:-1)
      end do

!     --------------------------------
!     Define active convective region.
!     --------------------------------

      ktop = k1 + 1
      k2m1 = k2 - 1


!     ----------------------------------------------------
!     Polar regions are too cold to have moist convection.
!     ----------------------------------------------------

!c?   For now, just do all latitude zones.
!c    jnp  = ilat
!c    jump = (jnp - 1) / 20
!c    js   = jump + 1
!c    jn   = jnp - js + 1
      js   = ju1
      jn   = j2


      ns   = Nint (tdt) / NSPLIT
      ns   = Max (ns, 1)
      rns  = ns
      sdt  = tdt / rns



      icloop: do ic = 1, numSpecies

!       ====================================
        if (isFixedConcentration(ic)) cycle icloop
!       ====================================

        isteploop: do istep = 1, ns

          do ij = js, jn
            do il = i1, i2

!             --------------------------------------------------------------
!             Warning:  It is assumed here that the lowest two layers are
!             treated as a single "cloud base" layer in the original cumulus
!             parameterization code that produced the cloud mass flux
!             and the detrainment.
!             --------------------------------------------------------------

              if (temp_cldmas(il,ij,k2m1) > TINY) then

!               -------------------------------------------
!               Compute mean mixing ratio below cloud base.
!               -------------------------------------------

                qb(il,ij) =  &
     &            (concentration(ic)%pArray3D(il,ij,k2)   * bmass(il,ij,k2) +  &
     &             concentration(ic)%pArray3D(il,ij,k2m1) * bmass(il,ij,k2m1)) /  &
     &            (bmass(il,ij,k2) + bmass(il,ij,k2m1))

!               ---------------------------------------------
!               Compute mixing ratio inside the cloud:
!               qb is the first guess; qc is the final value.
!               ---------------------------------------------

                tprod     =  temp_cldmas(il,ij,k2m1) * sdt
                wk(il,ij) =  bmass(il,ij,k2) + bmass(il,ij,k2m1)

                qc(il,ij) =  &
     &            (wk(il,ij) * qb(il,ij) +  &
     &             tprod * concentration(ic)%pArray3D(il,ij,k2-2)) /  &
     &            (wk(il,ij) + tprod)

!               ------------------------------------------------------------
!               Compute net change in mixing ratio.
!               Changes below cloud base are proportional to background mass
!               (but not to make the mixing ratio in the PBL negative).
!               wk is the total mass to be transported out of the PBL.
!               ------------------------------------------------------------

                dq = qb(il,ij) - qc(il,ij)

                if ((dq > concentration(ic)%pArray3D(il,ij,k2)) .or.  &
     &              (dq > concentration(ic)%pArray3D(il,ij,k2m1))) then

!                 ---------------------------
!                 Complete mixing in the PBL.
!                 ---------------------------

                  concentration(ic)%pArray3D(il,ij,k2m1) = qc(il,ij)
                  concentration(ic)%pArray3D(il,ij,k2)   = qc(il,ij)

                else

                  concentration(ic)%pArray3D(il,ij,k2m1) = concentration(ic)%pArray3D(il,ij,k2m1) - dq
                  concentration(ic)%pArray3D(il,ij,k2)   = concentration(ic)%pArray3D(il,ij,k2)   - dq

                end if

              else

                qc(il,ij) = concentration(ic)%pArray3D(il,ij,k2-2)

              end if

            end do
          end do


!         ---------
!         Interior.
!         ---------

          do ik = k2 - 2, ktop, -1
            do ij = js, jn
              do il = i1, i2

                if (temp_cldmas(il,ij,ik+1) > TINY) then

                  qcc   = qc(il,ij)
                  cmout = temp_cldmas(il,ij,ik) + temp_dtrn(il,ij,ik)
                  entrn = cmout - temp_cldmas(il,ij,ik+1)

!                 ----------------------------------------------------
!                 Entrainment must be >= 0.0d0.  The condition that
!                 entrn > 0 is strong enough to ensure that cmout > 0.
!                 ----------------------------------------------------

                  if (entrn > 0.0d0) then

!                   ---------------------
!                   Modify qc in layer-k.
!                   ---------------------

                    qc(il,ij) =  &
     &                (temp_cldmas(il,ij,ik+1) * qcc +  &
     &                 entrn * concentration(ic)%pArray3D(il,ij,ik)) /  &
     &                cmout

                  end if

!                 ---------------------
!                 Update qq in layer-k.
!                 ---------------------

                  concentration(ic)%pArray3D(il,ij,ik) =  &
     &              concentration(ic)%pArray3D(il,ij,ik) +  &
     &              sdt *  &
     &              (temp_cldmas(il,ij,ik)   *  &
     &               (concentration(ic)%pArray3D(il,ij,ik-1) - qc(il,ij)) +  &
     &               temp_cldmas(il,ij,ik+1) *  &
     &               (qcc               - concentration(ic)%pArray3D(il,ij,ik))) /  &
     &               bmass(il,ij,ik)

                else

!                 ---------------------------------------------
!                 No cloud transport if cloud mass flux < TINY;
!                 change qc to qq.
!                 ---------------------------------------------

                  qc(il,ij) = concentration(ic)%pArray3D(il,ij,ik)

                end if

              end do
            end do
          end do

        end do isteploop

      end do icloop

      do ic = 1, numSpecies
         do ij = js, jn
           where (concentration(ic)%pArray3D(:,ij,:) < 1.0d-30)
             concentration(ic)%pArray3D(:,ij,:) = 1.0d-30
           end where
         end do
      end do

!     -----------------------------------
!     Flip variables back in k dimension.
!     -----------------------------------

      do ic = 1, numSpecies
         concentration(ic)%pArray3D(:,:,k1:k2) = concentration(ic)%pArray3D(:,:,k2:k1:-1)
      end do

      return

      end subroutine Do_Convec_Dao2
