!------------------------------------------------------------------------------
!BOP
! !MODULE: GmiThermalRateConstants_mod
!
      module GmiThermalRateConstants_mod
!
! !USES:
      use GmiTimeControl_mod, only : GmiSplitDateTime
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!
      implicit none

      private
!
! !PUBLIC MEMBER FUNCTIONS:
      public  :: calcThermalRateConstants, Accum_Qqjk

#     include "GmiParameters.h"
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
!   - January 24, 2005 * Jules Kouatchou
!     Modified the routine Update_Qk.
!     Added the variables "sadcol2", "radA" (surface area)
!     and "rhcol" (relative humidity column) that are arguments of the
!     routine Kcalc.
!   - September 15, 2005 * Jules Kouatchou
!     Modified the routine Update_Qk.
!     Added  "phot_opt" as argument and the condition
!     "if ((phot_opt == 3) .or. (phot_opt == 8))" before
!     "if (do_AerDust_Calc)". This was done to make sure that "sadcol2"
!     and "radA" are initialized only if fastJ or fastJX is employed.
!
!=============================================================================


!-----------------------------------------------------------------------------
!BOP
!
! !ROUTINE: calcThermalRateConstants
!
! !INTERFACE:
!
! ARGUMENTS
!   do_wetchem         : do wet chemistry?
!   metdata_name_org   : first  part of metdata_name, e.g., "NCAR"
!   metdata_name_model : second part of metdata_name, e.g., "MATCH"
!   num_time_steps     : number of time steps
!   ih2o_num           : const array index for water vapor
!   imgas_num          : const array index for air density
!   nymd               : current year/month/day (YYYYMMDD)
!   rxnr_adjust_map    : mapping of reaction rate adjustment number to
!                        reaction rate number
!   latdeg      : latitude of grid centers (deg)
!   pres3c      : atmospheric pressure at the center of each grid box (mb)
!   tropp       : tropopause pressure (mb)
!   temp3       : temperature (degK)
!   clwc        : cloud liquid water content, grid box average (g/m^3)
!   cmf         : convective mass flux   (kg/m^2*s)
!   sadgmi      : surface area densities (cm^2/cm^3)
!   qkgmi       : thermal rate constants (units vary)
!   const       : species concentration, known at zone centers
!                 (mixing ratio)
!   rxnr_adjust : array of reaction rate adjustment factors
!
!-----------------------------------------------------------------------------

      subroutine calcThermalRateConstants  &
     &  (do_wetchem, metdata_name_org, metdata_name_model,  &
     &   num_time_steps, ih2o_num, imgas_num, nymd, rxnr_adjust_map,  &
     &   latdeg, pres3c, tropp, temp3, clwc, cmf, sadgmi, qkgmi,  &
     &   concentration, rxnr_adjust, Eradius, Tarea, relativeHumidity, &
     &   do_AerDust_Calc, phot_opt, &
     &   pr_diag, loc_proc, num_rxnr_adjust, rxnr_adjust_timpyr, &
     &   ivert, num_sad, num_qks, num_molefrac, num_species, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, chem_mecha)


      implicit none

#     include "gmi_time_constants.h"
#     include "gmi_phys_constants.h"
#     include "setkin_par.h"
#     include "gmi_AerDust_const.h"

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: num_species
      integer, intent(in) :: num_sad, num_qks, num_molefrac
      integer, intent(in) :: ivert
      integer, intent(in) :: num_rxnr_adjust, rxnr_adjust_timpyr
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: ju1_gl, j2_gl
      logical :: do_wetchem
      character (len=*) :: chem_mecha
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_model
      integer :: num_time_steps
      integer :: ih2o_num
      integer :: imgas_num
      integer :: nymd
      integer, intent(in) :: rxnr_adjust_map(num_rxnr_adjust)
      real*8  :: latdeg(ju1_gl:j2_gl)
      real*8  :: pres3c(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: tropp (i1:i2,   ju1:j2)
      real*8  :: temp3 (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: cmf   (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: sadgmi(i1:i2,   ju1:j2,   k1:k2, num_sad)
      real*8  :: qkgmi (i1:i2,   ju1:j2,   k1:k2, num_qks)
      real*8  :: rxnr_adjust(i1:i2, ju1:j2, k1:k2, num_rxnr_adjust, rxnr_adjust_timpyr)
      real*8  :: clwc  (i1:i2,   ju1:j2,   k1:k2)
      logical :: do_AerDust_Calc
      integer :: phot_opt
      REAL*8 :: ERADIUS (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
      REAL*8 :: TAREA   (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
      REAL*8 :: relativeHumidity (i1:i2, ju1:j2, k1:k2)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)
!
! !DEFINED PARAMETERS:
      logical, parameter :: DO_QK_STATS = .false.
!
! !DESCRIPTION:
! Updates the thermal rate constants (i.e., qk's).
!
! !LOCAL VARIABLES:
      integer :: idumyear, idumday
      integer :: il, ij, ic, im, iq
      real*8  :: adcol   (k1:k2)
      real*8  :: prescol (k1:k2)
      real*8  :: tempcol (k1:k2)
      real*8  :: cmfcol  (k1:k2)
      real*8  :: lwccol  (k1:k2)
      real*8  :: watcol  (k1:k2)
      real*8  :: constcol(num_molefrac, k1:k2)
      real*8  :: qkcol   (num_qks,      k1:k2)
      real*8  :: sadcol  (num_sad,      k1:k2)
      real*8  :: sadcol2 (NSADdust+NSADaer,   k1:k2)
      real*8  :: radA    (NSADdust+NSADaer,   k1:k2)
      real*8  :: rhcol   (              k1:k2)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) &
         Write (6,*) 'calcThermalRateConstants (old: Update_Qk) called by ', loc_proc

      lwccol (:) = 0.0d0
      qkcol(:,:) = 0.0d0

      do ij = ju1, j2
         do il = i1, i2
            prescol(:) = pres3c(il,ij,:)
            tempcol(:) = temp3 (il,ij,:)

            if (do_wetchem) then
               if ((metdata_name_org  (1:4) == 'NCAR') .and.  &
     &             (metdata_name_model(1:4) == 'CCM3')) then
                  lwccol(:) = clwc(il,ij,:)
               else
                  cmfcol(:) = cmf(il,ij,:) * SECPHR * GMI_G / 100.0d0
                  watcol(:) = concentration(ih2o_num)%pArray3D(il,ij,:)

                  call Calc_Lwc (ij, latdeg, cmfcol, prescol, tempcol, watcol, &
     &                  lwccol, pr_diag, loc_proc, k1, k2, ju1_gl, j2_gl)
               end if
            end if

            !==================================================='
            ! Variables needed for gas/heterogeneous chemistry.'
            ! adcol   = air density (molec/cm3)'
            ! FRH     = relative humidity fraction (0-1)'
            ! radA    = effective radius of aerosol (cm)'
            ! sadcol2 = surface area of aerosol/volume of air (cm2/cm3)'
            !==================================================='

            adcol(:) = concentration(imgas_num)%pArray3D(il,ij,:)

            do ic = 1, num_molefrac
               constcol(ic,:) = concentration(ic)%pArray3D(il,ij,:) * adcol(:)
            end do

            do ic = 1, num_sad
               sadcol(ic,:) = sadgmi(il,ij,:,ic)
            end do

            do ic = 1, NSADdust+NSADaer
               sadcol2(ic,:) = 0.0d0
               radA   (ic,:) = 0.0d0
            end do
            rhcol(:)     = 0.0d0

            if ((chem_mecha == 'troposphere').or. (chem_mecha == 'strat_trop')  &
     &         .or. chem_mecha == 'strat_trop_aerosol' ) then
               if ((phot_opt == 3) .or. (phot_opt == 8)) then
                  if (do_AerDust_Calc) then
                     do ic = 1, NSADdust+NSADaer
                        sadcol2(ic,:) = TAREA  (il,ij,:,ic)
                        radA   (ic,:) = ERADIUS(il,ij,:,ic)
                     end do
                     rhcol(:)     = relativeHumidity(il,ij,:)
                  end if
               end if
            end if

            call Kcalc (ivert, sadcol, sadcol2, prescol, tropp(il,ij),  &
     &             tempcol, lwccol, adcol, constcol, qkcol, radA, rhcol)

            do iq = 1, num_qks
               qkgmi(il,ij,:,iq) = Max (0.0d0, qkcol(iq,:))
            end do
         end do
      end do

      if (rxnr_adjust_timpyr == MONTHS_PER_YEAR) then
         call GmiSplitDateTime (nymd, idumyear, im, idumday)
      else
         im = 1
      end if

      do iq = 1, num_rxnr_adjust
         qkgmi(:,:,:,rxnr_adjust_map(iq)) =  &
     &         qkgmi(:,:,:,rxnr_adjust_map(iq)) * rxnr_adjust(:,:,:,iq,im)
      end do

      if (DO_QK_STATS) then
         call Calc_Qk_Stats (qkgmi, num_qks, num_time_steps, &
     &              pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2)
      end if

      return

      end subroutine calcThermalRateConstants
!EOC
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Calc_Lwc
!
! DESCRIPTION
!   This routine calculates liquid water in each grid box in a single column.
!   It uses parameterizations of cloud water from Kiehl and cloud
!   fractions from Xu.
!
! ARGUMENTS
!   ij      : latitude index for this column
!   latdeg  : latitude of grid box center   (deg)
!   cmfcol  : convective mass flux used for cloud fraction calculation
!             (g/cm^2/hour ???)
!   prescol : atmospheric pressure at the center of each grid box (mb)
!   tempcol : temperature (degK)
!   watcol  : water vapor concentration, known at zone centers (mixing ratio)
!   lwccol  : liquid water in each grid box (g/m^3)
!
!-----------------------------------------------------------------------------

      subroutine Calc_Lwc  &
     &  (ij, latdeg, cmfcol, prescol, tempcol, watcol, lwccol, &
     &   pr_diag, loc_proc, k1, k2, ju1_gl, j2_gl)

      implicit none

#     include "gmi_phys_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: k1, k2
      integer, intent(in) :: ju1_gl, j2_gl
      integer :: ij
      real*8  :: latdeg (ju1_gl:j2_gl)
      real*8  :: cmfcol (k1:k2)
      real*8  :: prescol(k1:k2)
      real*8  :: tempcol(k1:k2)
      real*8  :: watcol (k1:k2)
      real*8  :: lwccol (k1:k2)


!     ----------------------
!     Variable declarations.
!     ----------------------

      real*8  :: c0 (k1:k2)
      real*8  :: c1 (k1:k2)
      real*8  :: c2 (k1:k2)

      real*8  :: cfc(k1:k2)
      real*8  :: cfs(k1:k2)

      real*8  :: rh (k1:k2)
      real*8  :: rhc(k1:k2)

      real*8  :: height_asl(k1:k2)  ! height above sea level (m)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Calc_Lwc called by ', loc_proc
      end if


      c0 (:) = 0.0d0; c1(:) = 0.0d0; c2(:) = 0.0d0

      rhc(:) = 0.0d0


      height_asl(:) =  &
     &  -29.3d0 * tempcol(:) * Log (prescol(:) / 1013.25d0)

      height_asl(:) = Max (height_asl(:), 0.0d0)


!     -------------------------------------------------------------------
!     Calculate cloud water from a parameterization of Kiehl found in
!     "Sensitivity of the Simulated Climate to a Diagnostic Formulation
!     for Cloud Liquid Water" by James Hack, Journal of Climate, Vol. 11,
!     July 1998, p 1499.
!     -------------------------------------------------------------------

      lwccol(:) =  &
     &  0.18d0 *  &
     &  Exp (-height_asl(:) /  &
     &       (1080.0d0 + (2000.0d0 * Cos (latdeg(ij) * RADPDEG)**2)))

!     ----------------------------------------------------------------
!     Calculate convective and stratiform cloud fractions from a
!     parameterization of Xu found in "Evaluation of Cloudiness
!     Parameterizations Using a Cumulus Ensemble Model" by Kuan-Man Xu
!     and Steven Krueger, Monthly Weather Review, Vol. 119 p 342.
!     ----------------------------------------------------------------

!     --------------------------------------------------------------
!     First calculate relative humidity from Seinfeld (1986) p. 181.
!     The first  rh is the temperature dependent parameter a.
!     The second rh is the saturation vapor pressure of water.
!     The third  rh is the actual relative humidity as a fraction.
!     Then make sure rh is between 0 and 1.
!     --------------------------------------------------------------

      rh(:) = 1.0d0 - (373.15d0 / tempcol(:))

      rh(:) =  &
     &  1013.25d0 * Exp (13.3185d0 * rh(:)    -  &
     &                    1.9760d0 * rh(:)**2 -  &
     &                    0.6445d0 * rh(:)**3 -  &
     &                    0.1299d0 * rh(:)**4)

      rh(:) = watcol(:) * prescol(:) / rh(:)

      rh(:) = Max (Min (rh(:), 1.0d0), 0.0d0)

!     -----------------------------------------------------------------------
!     Now set up the constants to be used in the stratiform parameterization.
!     -----------------------------------------------------------------------

      where ((height_asl(:) < 15000.0d0) .and.  &
     &       (height_asl(:) >  5900.0d0))
        c1 (:) =  0.012d0
        c2 (:) =  0.321d0
        rhc(:) =  0.445d0
      end where

      where ((height_asl(:) <= 5900.0d0) .and.  &
     &       (height_asl(:) >  2350.0d0))
        c1 (:) = -0.102d0
        c2 (:) =  7.473d0
        rhc(:) =  0.715d0
      end where

      where  (height_asl(:) <= 2350.0d0)
        c1 (:) =  0.144d0
        c2 (:) =  8.050d0
        rhc(:) =  0.780d0
      end where

!     ----------------------------------------
!     Calculate the stratiform cloud fraction.
!     ----------------------------------------

      where (rh(:) > rhc(:))

        cfs(:) = c1(:) * (rh(:) - rhc(:)) +  &
     &           c2(:) * (rh(:) - rhc(:))**2

      elsewhere

        cfs(:) = 0.0d0

      end where

!     --------------------------------------------------------------------
!     Now set up the constants to be used in the cumulus parameterization.
!     --------------------------------------------------------------------

      where ((height_asl(:) < 15000.0d0) .and.  &
     &       (height_asl(:) >  5900.0d0))
        c0(:) = 0.0131d0
        c1(:) = 0.0123d0
        c2(:) = 0.0030d0
      end where

      where ((height_asl(:) <= 5900.0d0) .and.  &
     &       (height_asl(:) >  2350.0d0))
        c0(:) = 0.0722d0
        c1(:) = 0.0760d0
        c2(:) = 0.0254d0
      end where

      where  (height_asl(:) <= 2350.0d0)
        c0(:) = 0.0337d0
        c1(:) = 0.0453d0
        c2(:) = 0.0237d0
      end where

!     -------------------------------------
!     Calculate the cumulus cloud fraction.
!     -------------------------------------

      where (cmfcol(:) > 0.01d0)

        cfc(:) = c0(:) +  &
     &           c1(:) *  Log10 (cmfcol(:)) +  &
     &           c2(:) * (Log10 (cmfcol(:)))**2

      elsewhere

        cfc(:) = 0.0d0

      end where

!     ----------------------------------------------------------
!     Adjust liquid water in each grid box in this column due to
!     cloud fractions.
!     ----------------------------------------------------------

      lwccol(:) = Min (cfc(:)+cfs(:), 1.0d0) * lwccol(:)

      return

      end subroutine Calc_Lwc

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Calc_Qk_Stats
!
! DESCRIPTION
!   This routine writes out some qk statistics (min/maxes for each dim4).
!
! ARGUMENTS
!   rate_array     : thermal rate constants (i.e., qk's)
!   dim4           : 4th dimension of rate_array
!   num_time_steps : number of time steps
!
!-----------------------------------------------------------------------------

      subroutine Calc_Qk_Stats  &
     &  (rate_array, dim4, num_time_steps, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2)

      use GmiFlush_mod, only : GmiFlush

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer :: dim4
      integer :: num_time_steps
      real*8  :: rate_array(i1:i2, ju1:j2, k1:k2, dim4)

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: iq

      real*8  :: minrate
      real*8  :: maxrate

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Calc_Qk_Stats called by ', loc_proc
      end if

      Write (77,*) 'QK mins/maxes follow =>'
      Write (77,*) '  Step #:  ', num_time_steps+1
      Write (77,*) ' '


      do iq = 1, dim4

        minrate = Minval (rate_array(:,:,:,iq))
        maxrate = Maxval (rate_array(:,:,:,iq))

        Write (77,900) iq, minrate, maxrate

      end do

      Write (77,*) ' '

      call GmiFlush (77)

 900  format (i4, ', ', e20.10, ', ', e20.10)

      return

      end subroutine Calc_Qk_Stats

!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Accum_Qqjk
! 
! !INTERFACE:
! 
      subroutine Accum_Qqjk  &
     &  (do_qqjk_reset, imgas_num, concentration, qjgmi, qkgmi, qqjgmi, qqkgmi, &
     &   num_molefrac, num_species, num_qks, num_qjs, num_qjo, &
     &   pr_diag, loc_proc, ilong, i1, i2, ju1, j2, k1, k2)
!    
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!     
      implicit none
      
!
! !INPUT PARAMETERS:
!!   do_qqjk_reset : reset qqjk accumulators (i.e., qqjk output occurred on
!!                   previous time step)?
!!   imgas_num     : const array index for air density
!!   const         : species concentration, known at zone centers
!!                   (mixing ratio)
!!   qjgmi         : photolysis rate constants (s^-1)
!!   qkgmi         : thermal    rate constants (units vary)
!!   qqjgmi        : rates of photolytic processes (molecules/cm^3*s)
!!   qqkgmi        : rates of thermal    processes (molecules/cm^3*s)

      integer, intent(in) :: num_molefrac, num_species, num_qks, num_qjs, num_qjo
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: ilong
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2

      logical :: do_qqjk_reset
      integer :: imgas_num
      real*8  :: qjgmi (i1:i2, ju1:j2, k1:k2, num_qjo)
      real*8  :: qkgmi (i1:i2, ju1:j2, k1:k2, num_qks)
      real*8  :: qqjgmi(i1:i2, ju1:j2, k1:k2, num_qjo)
      real*8  :: qqkgmi(i1:i2, ju1:j2, k1:k2, num_qks)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!
! !DESCRIPTION:
!   This routine accumulates the rates of the photolysis and thermal
!   processes.
!
! !LOCAL VARIABLES:
      !logical, save :: first = .true.
      integer :: ic, iq, il, ij, ik
      real*8  :: const2 (i1:i2, imgas_num-1)
      real*8  :: qjgmi2 (i1:i2, num_qjs)
      real*8  :: qqjgmi2(i1:i2, num_qjs)
      real*8  :: qkgmi2 (i1:i2, num_qks)
      real*8  :: qqkgmi2(i1:i2, num_qks)
!
! !REVISION HISTORY:
!   Initial code.
!EOP
!-------------------------------------------------------------------------
!BOC

      if (pr_diag) Write (6,*) 'Accum_Qqjk called by ', loc_proc

!      if (first) then
!        first = .false.

        do_qqjk_reset = .true.
!      end if

      if (do_qqjk_reset) then
        do_qqjk_reset = .false.

        qqjgmi(:,:,:,:) = 0.0d0
        qqkgmi(:,:,:,:) = 0.0d0
      end if

      const2 (:,:) = 0.0d0
      qjgmi2 (:,:) = 0.0d0
      qqjgmi2(:,:) = 0.0d0
      qkgmi2 (:,:) = 0.0d0
      qqkgmi2(:,:) = 0.0d0

      do ik = k1, k2
        do ij = ju1, j2

          do ic = 1, num_molefrac
            do il = i1, i2
              const2(il,ic) =  concentration(ic)%pArray3D(il,ij,ik) * &
     &                         concentration(imgas_num)%pArray3D(il,ij,ik)
            end do
          end do

          if ((imgas_num - 1) > num_molefrac) then
            do ic = num_molefrac+1, imgas_num-1
              do il = i1, i2
                const2(il,ic) = concentration(ic)%pArray3D(il,ij,ik)
              end do
            end do
          end if

          do iq = 1, num_qjs
            do il = i1, i2
              qjgmi2(il,iq) = qjgmi(il,ij,ik,iq)
            end do
          end do

          do iq = 1, num_qks
            do il = i1, i2
              qkgmi2(il,iq) = qkgmi(il,ij,ik,iq)
            end do
          end do

!         =====================
          call Calc_Rate_Setkin  &
!         =====================
     &      (ilong, num_qjs, num_qks, imgas_num-1,  &
     &       qkgmi2, qjgmi2, const2, qqkgmi2, qqjgmi2)

          do iq = 1, num_qjs
            do il = i1, i2
              qqjgmi(il,ij,ik,iq) =  &
     &          qqjgmi2(il,iq) / concentration(imgas_num)%pArray3D(il,ij,ik)
            end do
          end do

          do iq = 1, num_qks
            do il = i1, i2
              qqkgmi(il,ij,ik,iq) =  &
     &          qqkgmi2(il,iq) / concentration(imgas_num)%pArray3D(il,ij,ik)
            end do
          end do

        end do
      end do

      return

      end subroutine Accum_Qqjk
!EOC
!--------------------------------------------------------------------------------


      end module GmiThermalRateConstants_mod
