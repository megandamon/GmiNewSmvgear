!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MUDULE: GmiUpdateSyntheticSpecies_mod
!
      module GmiUpdateSyntheticSpecies_mod
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: updateSyntheticSpecies
!
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
!
! !DESCRIPTION:
!  Routines for updating the synthetic species.
!
! !AUTHOR:
! John Tannahill, LLNL, jrt@llnl.gov
! Jules Kouatchou, Jules.Kouatchou@nasa.gov
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: updateSyntheticSpecies
!
! !INTERFACE:
!
      subroutine updateSyntheticSpecies  &
     &  (do_nodoz, latdeg, press3e, concentration, synoz_threshold, tdt, &
     &   mw, io3_num, isynoz_num, ihno3_num, &
     &   num_nox, num_noy, nox_map, noy_map, oz_eq_synoz_opt, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ju1_gl, j2_gl, &
     &   numSpecies, chem_mecha)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: ilo, ihi, julo, jhi, ju1_gl, j2_gl
      integer, intent(in) :: numSpecies
      integer, intent(in) :: oz_eq_synoz_opt
      integer, intent(in) :: io3_num, isynoz_num, ihno3_num
      integer, intent(in) :: num_nox, num_noy
      integer, intent(in) :: nox_map(1:), noy_map(1:)
      real*8 , intent(in) :: synoz_threshold
      real*8 , intent(in) :: tdt
      real*8 , intent(in) :: mw(numSpecies)
      logical, intent(in) :: do_nodoz
      real*8 , intent(in) :: latdeg (ju1_gl:j2_gl)
      real*8 , intent(in) :: press3e(ilo:ihi, julo:jhi, k1-1:k2)
      character (len=*), intent(in) :: chem_mecha
!
! !OUTPUT PARAMETERS:
      type (t_GmiArrayBundle), intent(inOut) :: concentration(numSpecies)
!
! !DESCRIPTION:
! Updates the synthetic species.
!
!EOP
!------------------------------------------------------------------------------
!BOC

!     ===============
      call Calc_Synoz  &
!     ===============
     &  (latdeg, press3e, concentration, tdt, oz_eq_synoz_opt, &
     &   mw, chem_mecha, io3_num, isynoz_num, numSpecies, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ju1_gl, j2_gl)


      if (do_nodoz) then

!       ===============
        call Calc_Nodoz  &
!       ===============
     &    (latdeg, press3e, concentration, synoz_threshold, tdt,  &
     &   mw, ihno3_num, isynoz_num, num_nox, num_noy, nox_map, noy_map, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ju1_gl, j2_gl, numSpecies)

      end if

      return

      end subroutine updateSyntheticSpecies
!EOC
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Calc_Synoz
!
! DESCRIPTION
!   This routine calculates Synoz emissions (mixing ratio/s) and removes
!   Synoz from the surface.
!
! ARGUMENTS
!   latdeg  : latitude (deg)
!   press3e : atmospheric pressure at the edge of each grid box (mb)
!   const   : species concentration, known at zone centers (mixing ratio)
!
!-----------------------------------------------------------------------------

      subroutine Calc_Synoz  &
     &  (latdeg, press3e, concentration, tdt, oz_eq_synoz_opt, &
     &   mw, chem_mecha, io3_num, isynoz_num, numSpecies, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ju1_gl, j2_gl)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      implicit none


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in   ) :: ilo, ihi, julo, jhi, ju1_gl, j2_gl
      integer, intent(in   ) :: numSpecies
      integer, intent(in   ) :: oz_eq_synoz_opt
      integer, intent(in   ) :: io3_num, isynoz_num
      real*8 , intent(in   ) :: mw(numSpecies)
      real*8 , intent(in   ) :: tdt
      real*8 , intent(in   ) :: latdeg (ju1_gl:j2_gl)
      real*8 , intent(in   ) :: press3e(ilo:ihi, julo:jhi, k1-1:k2)
      character (len=*), intent(in) :: chem_mecha
      type (t_GmiArrayBundle), intent(inOut) :: concentration(numSpecies)
!     -----------------------
!     Parameter declarations.
!     -----------------------

      integer, parameter ::  &
     &  KSTX_REM = 1,  &
     &  KEND_REM = 2

      real*8, parameter  ::  &
     &  LOLAT = -30.0d0,  &
     &  HILAT =  30.0d0

      real*8, parameter  ::  &
     &  LOPRS =  10.0d0,  &
     &  HIPRS =  70.0d0

      real*8, parameter  ::  &
     &  SYNOZ_RELAX_DAYS =   2.0d0,  &
     &  SYNOZ_SRCFLUX    = 550.0d0,    & ! Synoz source flux (tg/yr)
     &  SYNOZ_SRFVAl     =   2.5d-8


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij
      integer :: jsedge, jnedge
      integer :: kstx_emi, kend_emi

      logical, save :: doit     = .false.
      logical, save :: first    = .true.

      integer, save :: jstx_emi = -999
      integer, save :: jend_emi = -999

      real*8  :: rsecpdy, rsecpyr
      real*8  :: synoz_adjustment
      real*8  :: synoz_srcregion   ! relative size of Synoz source region to
                                   ! full atmosphere

      real*8, save :: synozadj_constant

      real*8  :: latdeg_edge(ju1_gl:j2_gl+1)  ! latitude at zone edge (deg)


!     ----------------
!     Begin execution.
!     ----------------

!     ==========
      if (first) then
!     ==========

        first = .false.

!       ---------------------------------------------------------
!       Determine latdeg_edge, but stay away from the Poles.
!       Synoz doesn't need latdeg_edge close to the Poles anyway.
!       latdeg_edge(ij) is the southern edge in latitude.
!       ---------------------------------------------------------

        latdeg_edge(ju1_gl)   = -90.0d0
        latdeg_edge(ju1_gl+1) = -90.0d0

        do ij = ju1_gl+2, j2_gl-1
          latdeg_edge(ij)     = (latdeg(ij-1) + latdeg(ij)) * 0.5d0
        end do

        latdeg_edge(j2_gl)    =  90.0d0
        latdeg_edge(j2_gl+1)  =  90.0d0


        if (Any ((latdeg_edge(ju1:j2)     >= LOLAT) .and.  &
     &           (latdeg_edge(ju1+1:j2+1) <= HILAT))) then

          doit = .true.

          jsedge   =  &
     &      Transfer  &
     &        (Minloc (latdeg_edge,  &
     &                 mask = (latdeg_edge >= LOLAT) .and.  &
     &                        (latdeg_edge <  HILAT)),  &
     &         1)

          jnedge   =  &
     &      Transfer  &
     &        (Maxloc (latdeg_edge,  &
     &                 mask = (latdeg_edge >  LOLAT) .and.  &
     &                        (latdeg_edge <= HILAT)),  &
     &         1)


          jstx_emi =  &
     &      Transfer  &
     &        (Minloc (latdeg_edge(ju1:j2),  &
     &                 mask = (latdeg_edge(ju1:j2)     >= LOLAT) .and.  &
     &                        (latdeg_edge(ju1+1:j2+1) <= HILAT)),  &
     &         1)

          jstx_emi = jstx_emi  + (ju1 - 1)


          jend_emi =  &
     &      Transfer  &
     &        (Maxloc (latdeg_edge(ju1:j2),  &
     &                 mask = (latdeg_edge(ju1:j2)     >= LOLAT) .and.  &
     &                        (latdeg_edge(ju1+1:j2+1) <= HILAT)),  &
     &         1)

          jend_emi = jend_emi + (ju1 - 1)


          synoz_srcregion =  &
     &      (Sin (latdeg_edge(jnedge) * RADPDEG) -  &
     &       Sin (latdeg_edge(jsedge) * RADPDEG)) *  &
     &      0.5d0 *  &
!c   &      ((HIPRS - LOPRS) / 993.407d0)
     &      ((HIPRS - LOPRS) / AVG_SRFPRS)


          rsecpyr = SECPYR

          synozadj_constant =  &
     &      tdt * SYNOZ_SRCFLUX *  &
     &      (MWTAIR / mw(isynoz_num)) /  &
     &      (rsecpyr * MASSATM * TGPKG * synoz_srcregion)

        end if

!     ======
      end if
!     ======


!     =========
      if (doit) then
!     =========

!       --------------------------
!       Calculate Synoz emissions.
!       --------------------------

        do ij = jstx_emi, jend_emi
          do il = i1, i2

            kstx_emi =  &
     &        Transfer  &
     &          (Maxloc (press3e(il,ij,k1:k2),  &
     &                   mask = (press3e(il,ij,0:k2-1) <= HIPRS) .and.  &
     &                          (press3e(il,ij,k1:k2)  >= LOPRS)),  &
     &           1)

            kstx_emi = kstx_emi + (k1 - 1)


            kend_emi =  &
     &        Transfer  &
     &          (Minloc (press3e(il,ij,k1:k2),  &
     &                   mask = (press3e(il,ij,0:k2-1) <= HIPRS) .and.  &
     &                          (press3e(il,ij,k1:k2)  >= LOPRS)),  &
     &           1)

            kend_emi = kend_emi + (k1 - 1)


            synoz_adjustment =  &
     &        synozadj_constant *  &
     &        ((HIPRS - LOPRS)  /  &
     &         (press3e(il,ij,kstx_emi-1) - press3e(il,ij,kend_emi)))


            concentration(isynoz_num)%pArray3D(il,ij,kstx_emi:kend_emi) =  &
     &        concentration(isynoz_num)%pArray3D(il,ij,kstx_emi:kend_emi) +  &
     &        synoz_adjustment

            IF ( oz_eq_synoz_opt /= 0 .and. chem_mecha /= 'strat_trop_aerosol'  &
     &           .and. chem_mecha /= 'strat_trop' ) THEN
               if (io3_num /= 0) then
                 if( (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol') .and. il.eq.i1 )  &
                    print *,'WARNING: ------ SYNOZ feeding back into O3 ------'
                 concentration(io3_num)%pArray3D(il,ij,kstx_emi:kend_emi) =  &
     &             concentration(io3_num)%pArray3D(il,ij,kstx_emi:kend_emi) +  &
     &             synoz_adjustment
               end if
            ENDIF

          end do
        end do

!     ======
      end if
!     ======


!     -----------------------------------
!     Calculate Synoz removal at surface.
!     -----------------------------------

      rsecpdy = SECPDY

      concentration(isynoz_num)%pArray3D(:,:,KSTX_REM:KEND_REM) =  &
     &  SYNOZ_SRFVAL +  &
     &  Max ((concentration(isynoz_num)%pArray3D(:,:,KSTX_REM:KEND_REM) - SYNOZ_SRFVAL),  &
     &       0.0d0) *  &
     &  Exp (-tdt / (SYNOZ_RELAX_DAYS * rsecpdy))


      return

      end subroutine Calc_Synoz


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Calc_Nodoz
!
! DESCRIPTION
!   This routine calculates Nodoz emissions (mixing ratio/s) and
!   repartitions NOx and NOy above the troposphere.
!
! ARGUMENTS
!   latdeg  : latitude (deg)
!   press3e : atmospheric pressure at the edge of each grid box (mb)
!   const   : species concentration, known at zone centers (mixing ratio)
!
!-----------------------------------------------------------------------------

      subroutine Calc_Nodoz  &
     &  (latdeg, press3e, concentration, synoz_threshold, tdt,  &
     &   mw, ihno3_num, isynoz_num, num_nox, num_noy, nox_map, noy_map, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ju1_gl, j2_gl, numSpecies)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in   ) :: ilo, ihi, julo, jhi, ju1_gl, j2_gl
      integer, intent(in   ) :: numSpecies
      integer, intent(in   ) :: ihno3_num, isynoz_num
      integer, intent(in   ) :: num_nox, num_noy
      integer, intent(in   ) :: nox_map(1:), noy_map(1:)
      real*8 , intent(in   ) :: synoz_threshold
      real*8 , intent(in   ) :: tdt
      real*8 , intent(in   ) :: mw(numSpecies)
      real*8 , intent(in   ) :: latdeg (ju1_gl:j2_gl)
      real*8 , intent(in   ) :: press3e(ilo:ihi, julo:jhi, k1-1:k2)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(numSpecies)


!     -----------------------
!     Parameter declarations.
!     -----------------------

      integer, parameter ::  &
     &  KEND_REM = 4,  &
     &  KSTX_REM = 1

      real*8, parameter ::  &
     &  HILAT =  30.0d0,  &
     &  LOLAT = -30.0d0

      real*8, parameter ::  &
     &  HIPRS  =  40.0d0,  &
     &  LOPRS  =   3.0d0,  &
     &  RLXPRS = 100.0d0

      real*8, parameter ::  &
     &  NODOZ_RELAX_DAYS = 3.0d0,  &
     &  NODOZ_SRCFLUX    = 2.295d0,    & ! Nodoz as HNO3 (tg/yr)
     &  NODOZ_SRFVAL     = 0.5d-9,  &
     &  RATIO            = 0.2d0,  &
     &  RATIO_RELAX_DAYS = 7.0d0


!     ----------------------
!     Variable declarations.
!     ----------------------

      logical, save :: doit  = .false.
      logical, save :: first = .true.

      integer :: il, ij, ic

!     --------------------------------------------------------------------
!     jnedge : global index of the northernmost latitude edge that lies on
!              or is south of HILAT; used to define the meridional extent
!              of the region of emission
!     jsedge : global index of the southernmost latitude edge that lies on
!              or is north of LOLAT; used to define the meridional extent
!              of the region of emission
!     --------------------------------------------------------------------

      integer :: jnedge, jsedge

!     -----------------------------------------------------------------
!     kend_emi : index of the highest (lowest pressure)  cell that lies
!                entirely between HIPRS and LOPRS
!     kstx_emi : index of the lowest  (highest pressure) cell that lies
!                entirely between HIPRS and LOPRS
!     -----------------------------------------------------------------

      integer :: kend_emi, kstx_emi

!     ----------------------------------------------------------------
!     kend_rlx : index in const of the highest (lowest pressure)  cell
!                with an upper edge below (higher pressure) RLXPRS
!     kstx_rlx : index in const of the lowest  (highest pressure) cell
!                with a Synoz value above SYNOZ_TROP (?)
!     ----------------------------------------------------------------

      integer :: kend_rlx, kstx_rlx

!     ------------------------------------------------------------------
!     jend_emi : index of the northernmost cell in the emission region
!     jstx_emi : index of the southernmost cell on the current processor
!                that is located entirely north of LOLAT
!     ------------------------------------------------------------------

      integer, save :: jstx_emi = -999
      integer, save :: jend_emi = -999

      real*8  :: nodoz_increment

!     --------------------------------------------------------------
!     nodoz_srcregion : relative size of Nodoz source region to full
!                       atmosphere
!     --------------------------------------------------------------

      real*8  :: nodoz_srcregion

      real*8  :: rsecpdy, rsecpyr

      real*8, save  :: nodoz_incconstant

!     -----------------------------------------
!     latdeg_edge : latitude at zone edge (deg)
!     -----------------------------------------

      real*8  :: latdeg_edge(ju1_gl:j2_gl+1)

      real*8  :: nox_adj  (k1:k2)
      real*8  :: nox_ratio(k1:k2)

      real*8  :: noy_adj  (k1:k2)

      real*8  :: sum_nox  (k1:k2)
      real*8  :: sum_noy  (k1:k2)


!     ----------------
!     Begin execution.
!     ----------------

!     ==========
      if (first) then
!     ==========

        first = .false.

!       ---------------------------------------------------------
!       Determine latdeg_edge, but stay away from the Poles.
!       Nodoz doesn't need latdeg_edge close to the Poles anyway.
!       latdeg_edge(ij) is the southern edge in latitude.
!       ---------------------------------------------------------

        latdeg_edge(ju1_gl)   = -90.0d0
        latdeg_edge(ju1_gl+1) = -90.0d0

        do ij = ju1_gl+2, j2_gl-1
          latdeg_edge(ij)     = (latdeg(ij-1) + latdeg(ij)) * 0.5d0
        end do

        latdeg_edge(j2_gl)    =  90.0d0
        latdeg_edge(j2_gl+1)  =  90.0d0


        if (Any ((latdeg_edge(ju1:j2)     >= LOLAT) .and.  &
     &           (latdeg_edge(ju1+1:j2+1) <= HILAT))) then

          doit = .true.

!         -----------------------------------------------------------------
!         The masks in the following two equations preselect only those
!         edges contained within LOLAT and HILAT.
!
!         Note that jsedge must be < jnedge for emission of Nodoz to occur.
!         -----------------------------------------------------------------

          jnedge =  &
     &      Transfer  &
     &        (Maxloc (latdeg_edge,  &
     &                 mask = (latdeg_edge >  LOLAT) .and.  &
     &                        (latdeg_edge <= HILAT)),  &
     &         1)

          jsedge =  &
     &      Transfer  &
     &        (Minloc (latdeg_edge,  &
     &                 mask = (latdeg_edge >= LOLAT) .and.  &
     &                        (latdeg_edge <  HILAT)),  &
     &         1)


!         -------------------------------------------------------------
!         The masks in the following two equations preselect only those
!         cells that are contained entirely within LOLAT and HILAT.
!         The result of the Minloc function is the smallest index,
!         counting from 1 in the subselected latdeg array, that
!         corresponds to a cell with a southern boundary lying on or
!         north of LOLAT and a northern boundary on or south of HILAT.
!         -------------------------------------------------------------

          jend_emi =  &
     &      Transfer  &
     &        (Maxloc (latdeg_edge(ju1:j2),  &
     &                 mask = (latdeg_edge(ju1:j2)     >= LOLAT) .and.  &
     &                        (latdeg_edge(ju1+1:j2+1) <= HILAT)),  &
     &         1)

          jstx_emi =  &
     &      Transfer  &
     &        (Minloc (latdeg_edge(ju1:j2),  &
     &                 mask = (latdeg_edge(ju1:j2)     >= LOLAT) .and.  &
     &                        (latdeg_edge(ju1+1:j2+1) <= HILAT)),  &
     &         1)

!         ------------------------------------------------------------
!         (ju1 - 1) below is the offset that must be added to jend_emi
!         and jstx_emi for indexing into the const array.
!         ------------------------------------------------------------

          jend_emi = jend_emi + (ju1 - 1)
          jstx_emi = jstx_emi + (ju1 - 1)


!         ------------------------------------------------------------
!         Note that AVG_SRFPRS below should be the true annual average
!         surface pressure for the global met field used.
!         ------------------------------------------------------------

          nodoz_srcregion =  &
     &      (Sin (latdeg_edge(jnedge) * RADPDEG) -  &
     &       Sin (latdeg_edge(jsedge) * RADPDEG)) *  &
     &      0.5d0 *  &
!c   &      ((HIPRS - LOPRS) / 993.407d0)
     &      ((HIPRS - LOPRS) / AVG_SRFPRS)


          rsecpyr = SECPYR

          nodoz_incconstant =  &
     &      tdt * NODOZ_SRCFLUX *  &
     &      (MWTAIR / mw(ihno3_num)) /  &
     &      (rsecpyr * MASSATM * TGPKG * nodoz_srcregion)

        end if

!     ======
      end if
!     ======


!     =========
      if (doit) then
!     =========

!       --------------------------
!       Calculate Nodoz emissions.
!       --------------------------

        do ij = jstx_emi, jend_emi
          do il = i1, i2

            kend_emi =  &
     &        Transfer  &
     &          (Minloc (press3e(il,ij,k1:k2),  &
     &                   mask = (press3e(il,ij,0:k2-1) <= HIPRS) .and.  &
     &                          (press3e(il,ij,k1:k2)  >= LOPRS)),  &
     &           1)

            kstx_emi =  &
     &        Transfer  &
     &          (Maxloc (press3e(il,ij,k1:k2),  &
     &                   mask = (press3e(il,ij,0:k2-1) <= HIPRS) .and.  &
     &                          (press3e(il,ij,k1:k2)  >= LOPRS)),  &
     &           1)

!           ----------------------------------------------------------
!           If, at some future time the parallel decomposition extends
!           to sending portions of the vertical column to different
!           processors, the k1 offset would need to be accounted for.
!           ----------------------------------------------------------

!c          kend_emi = kend_emi + (k1 - 1)
!c          kstx_emi = kstx_emi + (k1 - 1)


!           -----------------------------------------------------------
!           The high pressure edge index here (kstx_emi) is decremented
!           to account for the zero-based indexing of press3e.
!           -----------------------------------------------------------

            nodoz_increment =  &
     &        nodoz_incconstant *  &
     &        ((HIPRS - LOPRS)  /  &
     &         (press3e(il,ij,kstx_emi-1) - press3e(il,ij,kend_emi)))


!           -------------------------------
!           Actual Nodoz emission operator.
!           -------------------------------

            concentration(ihno3_num)%pArray3D(il,ij,kstx_emi:kend_emi) =  &
     &         concentration(ihno3_num)%pArray3D(il,ij,kstx_emi:kend_emi) +  &
     &        nodoz_increment

          end do
        end do

!     ======
      end if
!     ======

      if ((num_nox == 0) .and. (num_noy == 0)) then

!     -----------------------------------
!     Calculate Nodoz removal at surface.
!     -----------------------------------

      rsecpdy = SECPDY

       concentration(ihno3_num)%pArray3D(:,:,KSTX_REM:KEND_REM) =  &
     &  NODOZ_SRFVAL +  &
     &  Max (( concentration(ihno3_num)%pArray3D(:,:,KSTX_REM:KEND_REM) - NODOZ_SRFVAL),  &
     &       0.0d0) *  &
     &  Exp (-tdt / (NODOZ_RELAX_DAYS * rsecpdy))



      else

!     --------------------------
!     Calculate NOy partitioning.
!     --------------------------

      do ij = ju1, j2
        do il = i1, i2

          kend_rlx =  &
     &      Transfer  &
     &        (Minloc (press3e(il,ij,k1:k2),  &
     &                 mask = (press3e(il,ij,k1:k2) >= RLXPRS)),  &
     &         1)

          kstx_rlx =  &
     &      Transfer  &
     &        (Maxloc (press3e(il,ij,k1:k2),  &
     &                 mask = ( concentration(isynoz_num)%pArray3D(il,ij,k1:k2) >  &
     &                         synoz_threshold)),  &
     &         1)

!         ----------------------------------------------------------
!         If, at some future time the parallel decomposition extends
!         to sending portions of the vertical column to different
!         processors, the k1 offset would need to be accounted for.
!         ----------------------------------------------------------

!c        kend_rlx = kend_rlx + (k1 - 1)
!c        kstx_rlx = kstx_rlx + (k1 - 1)


!         ------------------------------------
!         Assemble NOx and NOy sums and ratio.
!         ------------------------------------

          sum_nox(:) = 0.0d0

          do ic = 1, num_nox

            sum_nox = sum_nox + concentration(nox_map(ic))%pArray3D(il,ij,:)

          end do

          sum_noy(:) = sum_nox(:)

          do ic = 1, num_noy

            sum_noy = sum_noy +  concentration(noy_map(ic))%pArray3D(il,ij,:)

          end do

          nox_ratio(:) = sum_nox(:) / sum_noy(:)


!         ----------------------------------------------------------
!         Between the Synoz tropopause and RLXPRS, relax the NOx/NOy
!         ratio to the imposed RATIO with a RATIO_RELAX_DAYS
!         relaxation time constant.
!         ----------------------------------------------------------

          rsecpdy    = SECPDY

          nox_adj(:) = 1.0d0

          if (kend_rlx >= kstx_rlx) then

            nox_adj(kstx_rlx:kend_rlx) =  &
     &        (RATIO +  &
     &         (nox_ratio(kstx_rlx:kend_rlx) - RATIO) *  &
     &         Exp(-tdt / (RATIO_RELAX_DAYS * rsecpdy))) /  &
     &        nox_ratio(kstx_rlx:kend_rlx)

          end if


!         --------------------------------------------------
!         Above RLXPRS (at lower pressures), force NOx / NOy
!         ratio, retaining subgroup partitioning.
!         --------------------------------------------------

          nox_adj(kend_rlx+1:k2) = RATIO / nox_ratio(kend_rlx+1:k2)


!         ------------------
!         Adjust NOx values.
!         ------------------

          do ic = 1, num_nox

             concentration(nox_map(ic))%pArray3D(il,ij,kstx_rlx:k2) =  &
     &        nox_adj(kstx_rlx:k2) *  &
     &         concentration(nox_map(ic))%pArray3D(il,ij,kstx_rlx:k2)

          end do


!         ------------------
!         Adjust NOy values.
!         ------------------

          noy_adj(:) = 1.0d0

          noy_adj(kstx_rlx:k2) =  &
     &      (1.0d0 - nox_adj  (kstx_rlx:k2) * nox_ratio(kstx_rlx:k2)) /  &
     &      (1.0d0 - nox_ratio(kstx_rlx:k2))

          do ic = 1, num_noy

             concentration(noy_map(ic))%pArray3D(il,ij,kstx_rlx:k2) =  &
     &        noy_adj(kstx_rlx:k2) *  &
     &         concentration(noy_map(ic))%pArray3D(il,ij,kstx_rlx:k2)

          end do

        end do
      end do

      end if

      return

      end subroutine Calc_Nodoz

!------------------------------------------------------------------------------

      end module GmiUpdateSyntheticSpecies_mod
