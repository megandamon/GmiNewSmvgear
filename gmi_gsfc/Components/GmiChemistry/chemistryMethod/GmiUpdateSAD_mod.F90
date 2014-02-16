!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiUpdateSAD_mod
!
! !INTERFACE:
!
      module GmiUpdateSAD_mod
!
! !USED:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiTimeControl_mod       , only : GmiSplitDateTime
      use GmiReduce_mod            , only : Gmi_Min1_Reduce
!
      implicit none
!
#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"
#     include "gmi_sad_constants.h"
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public :: updateSurfaceAreaDensities

! !AUTHOR:
! David B. Considine, NASA GSFC
!
! !REVISION HISTORY:
!  Original code written by David B. Considine
!  LLNL modifications by John Tannahill, jrt@llnl.gov
!  14Dec2007, included the routines in a module - Jules Kouatchou
!  14Dec2007, passed tropp and nameOfModel to Update_Sad2 - Jules Kouatchou
!
!EOP
!=============================================================================

      CONTAINS
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: updateSurfaceAreaDensities
!
! !INTERFACE:
!
      subroutine updateSurfaceAreaDensities (rateintv, tropp, press3c, press3e,  &
     &                 kel, concentration, ch4clim, h2oclim, hno3cond, hno3gas,  &
     &                 lbssad, sadgmi, h2oback, h2ocond, reffice, reffsts, vfall,&
     &                 dehydmin, dehyd_opt, h2oclim_opt, lbssad_opt, sad_opt,    &
     &                 ihno3cond_num, idehyd_num, ih2oaircr_num, ich4_num,       &
     &                 ihno3_num, ih2o_num, ih2ocond_num, nymd, pr_diag, procID, &
     &                 commuWorld, numDomains, num_species, num_sad,             &
     &                 lbssad_timpyr, h2oclim_timpyr, ju1_gl, j2_gl, ilo, ihi,   &
     &                 julo, jhi, i1, i2, ju1, j2, k1, k2, metdata_name_model,   &
     &                 chem_mecha)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer            , intent(in) :: ilo, ihi, julo, jhi
      integer            , intent(in) :: ju1_gl, j2_gl
      integer            , intent(in) :: i1, i2, ju1, j2, k1, k2
      integer            , intent(in) :: num_species, num_sad
      integer            , intent(in) :: lbssad_timpyr, h2oclim_timpyr
      integer            , intent(in) :: ih2o_num, ih2ocond_num, ihno3_num
      integer            , intent(in) :: ihno3cond_num, nymd
      integer            , intent(in) :: idehyd_num, ih2oaircr_num, ich4_num
      integer            , intent(in) :: procID, dehyd_opt
      integer            , intent(in) :: sad_opt, h2oclim_opt, lbssad_opt
      integer            , intent(in) :: numDomains, commuWorld
      logical            , intent(in) :: pr_diag
      real*8             , intent(in) :: rateintv
      real*8             , intent(in) :: kel    (ilo:ihi, julo:jhi, k1:k2)
      real*8             , intent(in) :: tropp      (i1:i2, ju1:j2)
      real*8             , intent(in) :: press3c(ilo:ihi, julo:jhi, k1:k2)
      real*8             , intent(in) :: press3e(ilo:ihi, julo:jhi, k1-1:k2)
      character (len=* ) , intent(in) :: chem_mecha
      character (len=MAX_LENGTH_MET_NAME) ,intent(in) :: metdata_name_model
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: dehydmin
      real*8  :: ch4clim (i1:i2,   ju1:j2,   k1:k2, h2oclim_timpyr)
      real*8  :: h2oclim (i1:i2,   ju1:j2,   k1:k2, h2oclim_timpyr)
      real*8  :: hno3cond(i1:i2,   ju1:j2,   k1:k2)
      real*8  :: hno3gas (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: lbssad  (i1:i2,   ju1:j2,   k1:k2, lbssad_timpyr)
      real*8  :: sadgmi  (i1:i2,   ju1:j2,   k1:k2, num_sad)
      real*8  :: h2oback (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: h2ocond (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: reffice (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: reffsts (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: vfall   (i1:i2,   ju1:j2,   k1:k2)

!
! !INPUT/OUTPUT PARAMETERS:
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)

!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
      integer :: idumday
      integer :: idumyear
      integer :: month
      real*8, allocatable  :: h2ocombo(:, :, :, :)
      real*8, allocatable  :: sadcombo(:, :, :, :)
      character(len=MAX_LENGTH_MET_NAME) :: nameOfModel
!
!EOP
!-----------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'updateSurfaceAreaDensities called by ', procID

      nameOfModel = ' '
      ! nameOfModel = TRIM(metdata_name_model) !<---- for running with GEOS5

      if (sad_opt == 1) then 
!
!... set aerosols (sadgmi array), hno3cond to zero and hno3gas to model gas phase HNO3
!
         call Update_Sad1 (ihno3_num, concentration, hno3cond, hno3gas, sadgmi, &
     &               pr_diag, procID, i1, i2, ju1, j2, k1, k2, num_sad,       &
     &               num_species)

      else if (sad_opt == 2) then
         if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then
            allocate(h2ocombo(i1:i2,  ju1:j2,  k1:k2,  h2oclim_timpyr))
            allocate(sadcombo(i1:i2,  ju1:j2,  k1:k2,  num_sad       ))

            ! Note: h2ocombo is not used when dehyd_opt = 0.
            ! ---------------------------------------------
            h2ocombo(:,:,:,:) = h2oclim(:,:,:,:)

            call GmiSplitDateTime (nymd, idumyear, month, idumday)
!
!... use model gas phase H2O in troposphere and H2O climotology above
!
            where (press3c(i1:i2,ju1:j2,:) > Spread (tropp(:,:), 3, k2))
               h2ocombo(:,:,:,month) = concentration(ih2o_num)%pArray3D(:,:,:)
            end where

            call Update_Sad2 (rateintv, tropp, press3c, press3e, kel,            &
     &                  concentration, ch4clim, h2ocombo, hno3cond, hno3gas,     &
     &                  lbssad, sadgmi, h2oback, h2ocond, reffice, reffsts,      &
     &                  vfall, dehydmin, dehyd_opt, h2oclim_opt, lbssad_opt,     &
     &                  ihno3cond_num, idehyd_num, ih2oaircr_num, ich4_num,      &
     &                  ihno3_num, ih2o_num, ih2ocond_num, nymd, pr_diag,        &
     &                  procID, ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2,      &
     &                  ju1, j2, k1, k2, lbssad_timpyr, h2oclim_timpyr, num_sad, &
     &                  num_species, numDomains, commuWorld, nameOfModel)

            sadcombo(:,:,:,:) = sadgmi(:,:,:,:)
!
            call Update_Sad3 (press3c, kel, concentration, lbssad, lbssad_opt,   &
     &                  sadcombo, pr_diag, procID, nymd, ih2o_num, ILBSSAD, i1, i2, ju1,  &
     &                  j2, k1, k2, ilo, ihi, julo, jhi, lbssad_timpyr, num_sad, &
     &                  num_species)
!
!... set lbssad array in troposphere back to what was calculated in Update_Sad2
!
            where (press3c(i1:i2,ju1:j2,:) > Spread (tropp(:,:), 3, k2))
               sadgmi(:,:,:,ILBSSAD) = sadcombo(:,:,:,ILBSSAD)
            end where

            deallocate(h2ocombo)
            deallocate(sadcombo)

         else

            call Update_Sad2 (rateintv, tropp, press3c, press3e, kel,            &
     &                  concentration, ch4clim, h2oclim, hno3cond, hno3gas,      &
     &                  lbssad, sadgmi, h2oback, h2ocond, reffice, reffsts,      &
     &                  vfall, dehydmin, dehyd_opt, h2oclim_opt, lbssad_opt,     &
     &                  ihno3cond_num, idehyd_num, ih2oaircr_num, ich4_num,      &
     &                  ihno3_num, ih2o_num, ih2ocond_num, nymd, pr_diag,        &
     &                  procID, ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2,      &
     &                  ju1, j2, k1, k2, lbssad_timpyr, h2oclim_timpyr, num_sad, &
     &                  num_species, numDomains, commuWorld, nameOfModel)

         endif
      else if (sad_opt == 3) then
!
!... adjust lbssad array with relative humidity
!
         call Update_Sad3 (press3c, kel, concentration, lbssad, lbssad_opt,      &
     &               sadgmi, pr_diag, procID, nymd, ih2o_num, ILBSSAD, i1, i2, ju1, j2,   &
     &               k1, k2, ilo, ihi, julo, jhi, lbssad_timpyr, num_sad,        &
     &               num_species)

      end if

      return

      end subroutine updateSurfaceAreaDensities
!EOC
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update_Sad1
!
! DESCRIPTION
!   This routine sets the stratospheric sulfuric acid aerosol surface area
!   densities.
!
! ARGUMENTS
!   ihno3_num : const array index for HNO3
!   const     : species concentration, known at zone centers (mixing ratio)
!   hno3cond  : condensed phase hno3 array (mixing ratio)
!   hno3gas   : gas       phase hno3 array (mixing ratio)
!   sadgmi    : surface area densities (cm^2/cm^3)
!
!-----------------------------------------------------------------------------

      subroutine Update_Sad1  &
     &  (ihno3_num, concentration, hno3cond, hno3gas, sadgmi, &
     &   pr_diag, procID, i1, i2, ju1, j2, k1, k2,  &
     &   num_sad, numSpecies)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID, i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: num_sad, numSpecies
      integer :: ihno3_num
      real*8  :: hno3cond(i1:i2, ju1:j2, k1:k2)
      real*8  :: hno3gas (i1:i2, ju1:j2, k1:k2)
      real*8  :: sadgmi  (i1:i2, ju1:j2, k1:k2, num_sad)
      type (t_GmiArrayBundle), intent(in) :: concentration(numSpecies)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Update_Sad1 called by ', procID
      end if


      sadgmi(:,:,:,:) = 0.0d0

      hno3cond(:,:,:) = 0.0d0

      hno3gas (:,:,:) = concentration(ihno3_num)%pArray3D(:,:,:)


      return

      end subroutine Update_Sad1


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update_Sad2
!
! DESCRIPTION
!   This routine sets the stratospheric sulfuric acid aerosol surface area
!   densities.
!
! ARGUMENTS
!   sadintv  : surface area density time step (s)
!   pres3c   : atmospheric pressure at the center of each grid box (mb)
!   pres3e   : atmospheric pressure at the edge   of each grid box (mb)
!   temp3    : temperature (degK)
!   const    : species concentration, known at zone centers (mixing ratio)
!   ch4clim  : array of ch4 climatology
!   h2oclim  : array of h2o climatology
!   hno3cond : condensed phase hno3 array (mixing ratio)
!   hno3gas  : gas       phase hno3 array (mixing ratio)
!   lbssad   : liquid binary sulfate background surface area density (cm^-1)
!   sadgmi   : surface area densities (cm^2/cm^3)
!   h2oback  : background      h2o array (mixing ratio)
!   h2ocond  : condensed phase h2o array (mixing ratio)
!   reffice  : effective radius of ICE aerosols  (cm)
!   reffsts  : effective radius of STS aerosols  (cm)
!   vfall    : effective aerosol fall velocities (cm/s)
!   dehyd_opt: 0=disabled (AGCM mode) 1=enabled (dehyd available)
!   h2oclim_opt : 0=disabled (Using transported CH4 and H20 instead) 1=enabled.
!
!-----------------------------------------------------------------------------

      subroutine Update_Sad2  &
     &  (sadintv, tropp, pres3c, pres3e, temp3, concentration, ch4clim, &
     &   h2oclim, hno3cond, hno3gas, lbssad, sadgmi, h2oback, h2ocond,  &
     &   reffice, reffsts, vfall, dehydmin, dehyd_opt, h2oclim_opt,     &
     &   lbssad_opt, ihno3cond_num, idehyd_num, ih2oaircr_num, ich4_num,&
     &   ihno3_num, ih2o_num, ih2ocond_num, nymd, pr_diag, procID, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi,&
     &   i1, i2, ju1, j2, k1, k2, lbssad_timpyr, h2oclim_timpyr, num_sad, &
     &   numSpecies, numDomains, commuWorld, nameOfModel)

      implicit none


!     ----------------------
!     Argument declarations.
!     ----------------------
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID, commuWorld
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: ju1_gl, j2_gl
      integer, intent(in) :: lbssad_timpyr, h2oclim_timpyr, num_sad, numSpecies
      integer, intent(in) :: nymd
      INTEGER, INTENT(IN) :: dehyd_opt, h2oclim_opt, lbssad_opt, ih2ocond_num
      integer, intent(in) :: numDomains
      integer, intent(in) :: ihno3_num, idehyd_num, ih2oaircr_num, ich4_num, ih2o_num
      integer, intent(in) :: ihno3cond_num
      real*8 , intent(out) :: dehydmin
      real*8  :: sadintv
      real*8,  intent(in) :: tropp(i1:i2, ju1:j2)
      real*8  :: pres3c  (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: pres3e  (ilo:ihi, julo:jhi, k1-1:k2)
      real*8  :: temp3   (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: ch4clim (i1:i2,   ju1:j2,   k1:k2, h2oclim_timpyr)
      real*8  :: h2oclim (i1:i2,   ju1:j2,   k1:k2, h2oclim_timpyr)
      real*8  :: hno3cond(i1:i2,   ju1:j2,   k1:k2)
      real*8  :: hno3gas (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: lbssad  (i1:i2,   ju1:j2,   k1:k2, lbssad_timpyr)
      real*8  :: sadgmi  (i1:i2,   ju1:j2,   k1:k2, num_sad)
      real*8  :: h2oback (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: h2ocond (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: reffice (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: reffsts (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: vfall   (i1:i2,   ju1:j2,   k1:k2)
      CHARACTER(LEN=*), INTENT(IN) :: nameOfModel
      type (t_GmiArrayBundle), intent(inout) :: concentration(numSpecies)


!     -----------------------
!     Parameter declarations.
!     -----------------------

      logical, parameter :: WR_ARGS = .false.

      integer, parameter :: WR_PROC =  1

      integer, parameter :: X1 = 36
      integer, parameter :: X2 = 36

      integer, parameter :: Y1 =  4
      integer, parameter :: Y2 =  4

      integer, parameter :: Z1 =  1
      integer, parameter :: Z2 =  1


!     ----------------------
!     Variable declarations.
!     ----------------------

      logical :: is_before

      integer :: idumday
      integer :: idumyear
      integer :: ik
      integer :: month

      real*8  :: fac

      real*8  :: dehyd   (i1:i2, ju1:j2, k1:k2)
      real*8  :: denssts (i1:i2, ju1:j2, k1:k2)
      real*8  :: dz      (i1:i2, ju1:j2, k1:k2)

      real*8  :: h2so4gas(i1:i2, ju1:j2, k1:k2)
      real*8  :: reffnat (i1:i2, ju1:j2, k1:k2)

      real*8  :: sadgmi_iice(i1:i2, ju1:j2, k1:k2)
      real*8  :: sadgmi_ilbs(i1:i2, ju1:j2, k1:k2)
      real*8  :: sadgmi_inat(i1:i2, ju1:j2, k1:k2)
      real*8  :: sadgmi_ists(i1:i2, ju1:j2, k1:k2)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Update_Sad2 called by ', procID
      end if


      denssts (:,:,:) = 0.0d0

      h2so4gas(:,:,:) = 0.0d0; reffnat(:,:,:) = 0.0d0


      call GmiSplitDateTime (nymd, idumyear, month, idumday)

      fac = ((GAS_CONST_J * GPKG) / (GMI_G * MWTAIR)) * CMPM


      hno3cond(:,:,:) = 0.0d0

      hno3gas(:,:,:) = concentration(ihno3_num)%pArray3D(:,:,:)
      IF (h2oclim_opt == 3) THEN
         hno3gas(:,:,:) = hno3gas(:,:,:)+concentration(ihno3cond_num)%pArray3D(:,:,:)
      END IF

      IF (dehyd_opt == 0) THEN
         dehyd(:,:,:) = 0.0d0
         h2ocond(:,:,:) = concentration(ih2ocond_num)%pArray3D(:,:,:)
      ELSE
         dehyd(:,:,:) = concentration(idehyd_num)%pArray3D(:,:,:)
      END IF

      IF (h2oclim_opt == 3) THEN
         h2oback(:,:,:) =  concentration(ih2o_num)%pArray3D(:,:,:) +  &
     &          concentration(ih2oaircr_num)%pArray3D(:,:,:) + h2ocond(:,:,:)
      ELSE
         h2oback(:,:,:) = concentration(ih2oaircr_num)%pArray3D(:,:,:) +  &
     &         h2oclim(:,:,:,month) +  2.0d0 *  &
     &  (concentration(ich4_num)%pArray3D(:,:,:) - ch4clim(:,:,:,month))
      END IF

      IF (lbssad_opt == 4) THEN
         sadgmi(:,:,:,ILBSSAD) = lbssad(:,:,:,1)
      ELSE
         sadgmi(:,:,:,ILBSSAD) = lbssad(:,:,:,month)
      END IF

      do ik = k1, k2

        dz(i1:i2,ju1:j2,ik) =  &
     &    fac * temp3(i1:i2,ju1:j2,ik) *  &
     &    Log (pres3e(i1:i2,ju1:j2,ik-1) / pres3e(i1:i2,ju1:j2,ik))

      end do


      sadgmi_iice(:,:,:) = sadgmi(:,:,:,IICESAD)
      sadgmi_ilbs(:,:,:) = sadgmi(:,:,:,ILBSSAD)
      sadgmi_inat(:,:,:) = sadgmi(:,:,:,INATSAD)
      sadgmi_ists(:,:,:) = sadgmi(:,:,:,ISTSSAD)

      if (WR_ARGS .and. (procID == WR_PROC)) then
        is_before = .true.
        call Write_Cond_Args  &
     &    (is_before, X1, X2, Y1, Y2, Z1, Z2,  &
     &     sadintv, pres3c, temp3, dehyd, dz,  &
     &     h2oback, hno3cond, hno3gas, sadgmi_ilbs,  &
     &     denssts, h2ocond, h2so4gas, sadgmi_iice,  &
     &     sadgmi_inat, sadgmi_ists, reffice, reffnat, reffsts, vfall, &
     &     pr_diag, procID, &
     &     ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)
      end if

!     =============
      call Condense  &
!     =============
     &  (sadintv, tropp, pres3c, temp3, dehyd, dz,  &
     &   h2oback, hno3cond, hno3gas, sadgmi_ilbs,  &
     &   denssts, h2ocond, h2so4gas, sadgmi_iice,  &
     &   sadgmi_inat, sadgmi_ists, reffice, reffnat, reffsts, vfall, &
     &   pr_diag, procID, ju1_gl, j2_gl, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2,&
     &   TRIM(nameOfModel))

      if (WR_ARGS .and. (procID == WR_PROC)) then
        is_before = .false.
        call Write_Cond_Args  &
     &    (is_before, X1, X2, Y1, Y2, Z1, Z2,  &
     &     sadintv, pres3c, temp3, dehyd, dz,  &
     &     h2oback, hno3cond, hno3gas, sadgmi_ilbs,  &
     &     denssts, h2ocond, h2so4gas, sadgmi_iice,  &
     &     sadgmi_inat, sadgmi_ists, reffice, reffnat, reffsts, vfall, &
     &     pr_diag, procID, &
     &     ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)
      end if

      sadgmi(:,:,:,IICESAD) = sadgmi_iice(:,:,:)
      sadgmi(:,:,:,ILBSSAD) = sadgmi_ilbs(:,:,:)
      sadgmi(:,:,:,INATSAD) = sadgmi_inat(:,:,:)
      sadgmi(:,:,:,ISTSSAD) = sadgmi_ists(:,:,:)


      call GmiSplitDateTime (nymd, idumyear, month, idumday)

!     For now, just zero out sootsad.
      sadgmi(:,:,:,ISOOTSAD) = 0.0d0

      IF(dehyd_opt == 0) THEN
        concentration(idehyd_num)%pArray3D(:,:,:) = 0.0d0
        concentration(ih2ocond_num)%pArray3D(:,:,:) = h2ocond(:,:,:)
      ELSE
        concentration(idehyd_num)%pArray3D(:,:,:) = dehyd  (:,:,:)
      END IF

      concentration(ihno3_num)%pArray3D(:,:,:)  = hno3gas(:,:,:)
      IF(h2oclim_opt == 3) THEN
        concentration(ihno3cond_num)%pArray3D(:,:,:)  = hno3cond(:,:,:)
      END IF

      IF(h2oclim_opt == 3) THEN
       concentration(ih2o_num)%pArray3D(:,:,:) =  &
     &  concentration(ih2o_num)%pArray3D(:,:,:) +  &
     &  concentration(ih2oaircr_num)%pArray3D  (:,:,:) -  &
     &  h2ocond(:,:,:)
      ELSE
       concentration(ih2o_num)%pArray3D(:,:,:) =  &
     &  h2oclim(:,:,:,month) +  &
     &  2.0d0 * (concentration(ich4_num)%pArray3D(:,:,:) - ch4clim(:,:,:,month)) +  &
     &  concentration(ih2oaircr_num)%pArray3D  (:,:,:) -  &
     &  h2ocond(:,:,:) -  &
     &  concentration(idehyd_num)%pArray3D  (:,:,:)
      END IF

      concentration(ih2o_num)%pArray3D(:,:,:) =  &
     &  Max (concentration(ih2o_num)%pArray3D(:,:,:), 0.01d-06)

!     ---------------------------------------------------------
!     Set dehydmin and do parallel communications if necessary.
!
!     Note: dehydmin is INTENT(OUT).  In GEOS-5 this value is
!     more easily generated in the GMI_GridCompMod run method.
!     ---------------------------------------------------------

      IF (nameOfModel == "GEOS-5" .OR. dehyd_opt == 0) RETURN

      dehydmin = Minval (concentration(idehyd_num)%pArray3D(:,:,:))

      dehydmin = Min (0.0d0, dehydmin)

!     =======================
      call Gmi_Min1_Reduce (dehydmin, numDomains, commuWorld)
!     =======================

      return

      end subroutine Update_Sad2

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update_Sad3
!
! DESCRIPTION
!   This routine sets the tropospheric sulfuric aerosol surface area
!   densities.
!
!.sds   Looks like it adjusts read in lbssad with local relative humidities
!
!
! ARGUMENTS
!   const    : species concentration, known at zone centers (mixing ratio)
!   hno3cond : condensed phase hno3 array (mixing ratio)
!   hno3gas  : gas       phase hno3 array (mixing ratio)
!   sadgmi   : surface area densities (cm^2/cm^3)
!
!-----------------------------------------------------------------------------

      subroutine Update_Sad3  &
     &  (pres3c, temp3, concentration, lbssad, lbssad_opt, sadgmi, &
     &  pr_diag, procID, nymd, ih2o_num, ILBSSAD, &
     &  i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, &
     &  lbssad_timpyr, num_sad, numSpecies)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: nymd
      integer, intent(in) :: ih2o_num, ILBSSAD
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in) :: lbssad_timpyr, lbssad_opt, num_sad, numSpecies
      real*8, intent(in)  :: pres3c(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: temp3 (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: lbssad(i1:i2,   ju1:j2,   k1:k2, lbssad_timpyr)
      real*8, intent(inout) :: sadgmi(i1:i2,   ju1:j2,   k1:k2, num_sad)
      type (t_GmiArrayBundle), intent(in) :: concentration(numSpecies)

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: idumday, idumyear
      integer :: month

      real*8  :: rh(i1:i2, ju1:j2, k1:k2)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Update_Sad3 called by ', procID
      end if

!     -----------------------------------------------------------
!     Calculate relative humidity from Seinfeld (1986) p. 181.
!     The first rh is the temperature dependent parameter a.
!     The second rh is the saturation vapor pressure of water.
!     The third rh is the actual relative humidity as a fraction.
!     Then make sure rh is between 0 and 1.
!     -----------------------------------------------------------

      rh(:,:,:) = 1.0d0 - (373.15d0 / temp3(i1:i2,ju1:j2,:))

      rh(:,:,:) =  &
     &  1013.25d0 * Exp (13.3185d0 * rh(:,:,:)    -  &
     &                    1.9760d0 * rh(:,:,:)**2 -  &
     &                    0.6445d0 * rh(:,:,:)**3 -  &
     &                    0.1299d0 * rh(:,:,:)**4)

      rh(:,:,:) =  &
     &  concentration(ih2o_num)%pArray3D(:,:,:) *  &
     &  pres3c(i1:i2,ju1:j2,:) / rh(:,:,:)

      rh(:,:,:) = Max (Min (rh(:,:,:), 1.0d0), 0.0d0)

!     -----------------------------------------------------------------
!     Choose month only when all 12 months of lbssad are available.
!     -----------------------------------------------------------------

      IF(lbssad_opt == 4) THEN
       month = 1
      ELSE
       CALL GmiSplitDateTime (nymd, idumyear, month, idumday)
      END IF

!     -----------------------------------------------------------------
!     Now calculate the change in surface area density due to humidity.
!     This is a fit to data from Chuang in the form suggested by Grant.
!     -----------------------------------------------------------------

      sadgmi(:,:,:,ILBSSAD) =  &
     &  lbssad(:,:,:,month) *  &
     &  (1.0d0 +  &
     &   1.25824d0 * (Exp (2.71811d0 * rh(:,:,:)**11.3214) - 1.0d0) +  &
     &   24.94300d0 * (1.0d0 - Exp (-0.0329662 * rh(:,:,:))))**2


      return

      end subroutine Update_Sad3
!EOC
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Write_Cond_Args
!
! DESCRIPTION
!   This routine writes out some debug info.
!
! ARGUMENTS
!   Not important.
!
!-----------------------------------------------------------------------------

      subroutine Write_Cond_Args  &
     &  (is_before, x1, x2, y1, y2, z1, z2,  &
     &   sadintv, pres3c, temp3, dehyd, dz,  &
     &   h2oback, hno3cond, hno3gas, lbssad,  &
     &   denssts, h2ocond, h2so4gas, icesad, natsad, stssad,  &
     &   reffice, reffnat, reffsts, vfall, &
     &   pr_diag, procID, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)

      use GmiASCIIoperations_mod, only : AsciiOpenWrite
      use GmiFlush_mod          , only : GmiFlush

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      logical :: is_before
      integer :: x1, x2
      integer :: y1, y2
      integer :: z1, z2
      real*8  :: sadintv
      real*8  :: pres3c  (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: temp3   (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: dehyd   (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: dz      (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: h2oback (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: hno3cond(i1:i2,   ju1:j2,   k1:k2)
      real*8  :: hno3gas (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: lbssad  (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: denssts (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: h2ocond (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: h2so4gas(i1:i2,   ju1:j2,   k1:k2)
      real*8  :: icesad  (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: natsad  (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: stssad  (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: reffice (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: reffnat (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: reffsts (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: vfall   (i1:i2,   ju1:j2,   k1:k2)


!     ----------------------
!     Variable declarations.
!     ----------------------

      logical, save :: first = .true.

      integer, save :: iaft  = -1
      integer, save :: ibef  = -1
      integer       :: iu
      integer       :: ix, iy, iz


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6, *) 'Write_Cond_Args called by ', procID
      end if


!     ==========
      if (first) then
!     ==========

        first = .false.

        call AsciiOpenWrite (ibef, 'cond_bef.out')
        call AsciiOpenWrite (iaft, 'cond_aft.out')

      end if


      if (is_before) then
        iu = ibef
        Write (iu, *) ' '
        Write (iu, *) '=================================='
        Write (iu, *) 'ARGUMENTS BEFORE CALL TO CONDENSE:'
        Write (iu, *) '=================================='
      else
        iu = iaft
        Write (iu, *) ' '
        Write (iu, *) '================================='
        Write (iu, *) 'ARGUMENTS AFTER CALL TO CONDENSE:'
        Write (iu, *) '================================='
      end if
      Write (iu, *) ' '


      Write (iu, *) 'lon1, lon2, lat1, lat2, alt1, alt2:'
      Write (iu, *) ' ', x1, x2, y1, y2, z1, z2
      Write (iu, *) ' '


      Write (iu, *) 'sadintv:  ', sadintv
      Write (iu, *) ' '


      Write (iu, *) 'lat, alt:  '
      do iz = z1, z2
        do iy = y1, y2
          Write (iu, *) ' ', iy, iz
        end do
      end do
      Write (iu, *) ' '


      Write (iu, *) 'lon, lat, alt, pres3c(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, pres3c(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '

      Write (iu, *) 'lon, lat, alt, temp3(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, temp3(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '


      Write (iu, *) 'lon, lat, alt, dehyd(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, dehyd(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '

      Write (iu, *) 'lon, lat, alt, dz(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, dz(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '

      Write (iu, *) 'lon, lat, alt, denssts(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
          Write (iu, *) ' ', ix, iy, iz, denssts(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '


      Write (iu, *) 'lon, lat, alt, h2oback(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, h2oback(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '

      Write (iu, *) 'lon, lat, alt, h2ocond(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, h2ocond(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '

      Write (iu, *) 'lon, lat, alt, hno3gas(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, hno3gas(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '

      Write (iu, *) 'lon, lat, alt, hno3cond(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, hno3cond(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '

      Write (iu, *) 'lon, lat, alt, h2so4gas(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, h2so4gas(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '


      Write (iu, *) 'lon, lat, alt, lbssad(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, lbssad(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '

      Write (iu, *) 'lon, lat, alt, stssad(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, stssad(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '

      Write (iu, *) 'lon, lat, alt, natsad(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, natsad(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '

      Write (iu, *) 'lon, lat, alt, icesad(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, icesad(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '


      Write (iu, *) 'lon, lat, alt, reffsts(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, reffsts(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '

      Write (iu, *) 'lon, lat, alt, reffnat(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, reffnat(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '

      Write (iu, *) 'lon, lat, alt, reffice(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, reffice(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '


      Write (iu, *) 'lon, lat, alt, vfall(lon,lat,alt):  '
      do iz = z1, z2
        do iy = y1, y2
          do ix = x1, x2
            Write (iu, *) ' ', ix, iy, iz, vfall(ix,iy,iz)
          end do
        end do
      end do
      Write (iu, *) ' '

      call GmiFlush (iu)

      return

      end subroutine Write_Cond_Args
!EOC
!------------------------------------------------------------------------------
      end module GmiUpdateSAD_mod
