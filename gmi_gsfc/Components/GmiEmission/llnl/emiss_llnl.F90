!
!=============================================================================
!
! $Id: emiss_llnl.F90,v 1.16 2013-08-06 19:50:51 ssteenro Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   emiss_llnl.F
!
! ROUTINES
!   Add_Emiss_Llnl
!
! HISTORY
!   - December 8, 2005 * Bigyani Das
!     Changes are made for the aerocom run which are controlled by 2 parameters
!     do_aerocom and do_dust_emiss to add dust and sea salt emissions in
!     the model's  1st level. Originally sea salt emissions
!     were added in the emission fields with carbon emission,
!     now it is added with dust (dust and sea salt emissions are now
!     in the same files, before sea salt emissions were with carbon).
!  - October 2009 Steve Steenrod
!     Took out do_dust_emiss
!=============================================================================
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Add_Emiss_Llnl
!
! DESCRIPTION
!   This routine adds emissions to const.
!
!     1) Take emissions of kg/s and multiply by the time step to get total
!        kg of emission over time step.
!
!     2) Divide by mass of the zone to obtain emissions in terms of mixing
!        ratio; also multiply by the ratio of molecular weight of air to
!        molecular weight of the chemical emission to get volume mixing
!        ratio from mass mixing ratio.
!
!     3) Add emitted mixing ratio amount to existing mixing ratio of const.
!
! ARGUMENTS
!   pr_surf_emiss  : should the surface emissions be accumulated for output?
!   mcor           : surface area of each grid box (m^2)
!   surf_emiss_out : accumulated surface emissions for output (kg/m^2/time)
!   mass           : total mass of the atmosphere within each grid box (kg)
!   const          : species concentration, known at zone centers
!                    (mixing ratio)
!   emissionArray  : array of emissions (kg/s)
!   emiss_dust     : tbd
!   emiss_aero     :
!   pbl            :
!-----------------------------------------------------------------------------
!
      subroutine Add_Emiss_Llnl  &
     &  (do_aerocom,  &
     &   pr_surf_emiss, pr_emiss_3d, mcor, surf_emiss_out, emiss_3d_out,  &
     &   mass, concentration, emissionArray,  emiss_dust, emiss_aero_t, emiss_aero,  &
     &   pbl, gridBoxHeight, &
     &   IBOC, IBBC, INOC, IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3, ISSLT4, &
     &   IFSO2, INSO2, INDMS, IDUST1, IDUST2, IDUST3, IDUST4, IDUST5, &
     &   pr_diag, loc_proc, &
     &   chem_opt, emiss_aero_opt, emiss_dust_opt, &
     &   do_semiss_inchem, emiss_map, emiss_map_dust, emiss_map_aero, &
     &   ndust, naero, nymd, mw, tdt, &
     &   emiss_timpyr, num_emiss, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species)
!
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiTimeControl_mod  , only : GmiSplitDateTime
      use GmiSeaSaltMethod_mod, only : SourceSeaSalt
      use GmiSeaSaltMethod_mod, only : srcEmissSeaSalt
      use GmiDustMethod_mod   , only : SourceDust
      use GmiDustMethod_mod   , only : srcEmissDust
!
      implicit none
!
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
#     include "setkin_smv2par.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in   ) :: IBOC, IBBC, INOC, IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3, ISSLT4
      integer, intent(in   ) :: IFSO2, INSO2, INDMS, IDUST1, IDUST2, IDUST3, IDUST4, IDUST5
      logical, intent(in   ) :: pr_diag
      integer, intent(in   ) :: loc_proc
      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in   ) :: num_species, emiss_timpyr, num_emiss
      integer, intent(in   ) :: ndust, naero, nymd
      real*8 , intent(in   ) :: mw(num_species)
      real*8 , intent(in   ) :: tdt
      integer, intent(in   ) :: chem_opt, emiss_aero_opt, emiss_dust_opt
      integer, intent(in   ) :: emiss_map(num_emiss)
      integer, intent(inout) :: emiss_map_dust(num_species)
      integer, intent(inout) :: emiss_map_aero(num_species)
      logical, intent(in   ) :: do_semiss_inchem
      logical, intent(in   ) :: do_aerocom
      logical, intent(in   ) :: pr_surf_emiss, pr_emiss_3d
      real*8 , intent(in   ) :: mcor (i1:i2, ju1:j2)
      real*8 , intent(inout) :: surf_emiss_out(i1:i2, ju1:j2, num_species)
      real*8 , intent(inout) :: emiss_3d_out(i1:i2, ju1:j2, k1:k2, num_species)
      real*8 , intent(in   ) :: mass (i1:i2, ju1:j2, k1:k2)
!      real*8 , intent(in   ) :: emiss(i1:i2, ju1:j2, k1:k2, num_emiss)
      real*8 , intent(inout) :: emiss_dust(i1:i2, ju1:j2, ndust)
      real*8 , intent(in   ) :: emiss_aero_t(i1:i2, ju1:j2, naero, emiss_timpyr)
      real*8 , intent(out  ) :: emiss_aero(i1:i2, ju1:j2, naero)
      real*8 , intent(in   ) :: pbl (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: gridBoxHeight(i1:i2, ju1:j2, k1:k2)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)
      type (t_GmiArrayBundle), intent(in   ) :: emissionArray(num_emiss  )
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: ic, icx
      integer :: idumday, idumyear
      integer :: il, ij, ik
!      integer :: inum
      integer :: it
      integer :: kstrt
      integer :: month
!
      real*8  :: mw_fac
!
      real*8  :: mass_pbl(i1:i2, ju1:j2)
!
      real*8  :: emass(i1:i2, ju1:j2, k1:k2)
      real*8  :: bmass(i1:i2, ju1:j2, k1:k2)
!
      real*8  :: za(i1:i2, ju1:j2, k1:k2)
      real*8  :: zq(i1:i2, ju1:j2, k1-1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Add_Emiss_Llnl called by ', loc_proc
      end if
!
!
      emass(:,:,:) = 0.0d0
      bmass(:,:,:) = 0.0d0
!
!
      if (emiss_timpyr == MONTHS_PER_YEAR) then
        call GmiSplitDateTime  (nymd, idumyear, month, idumday)
        it = month
      else
        it = 1
      end if
!
!
      if (do_semiss_inchem) then
        kstrt = k1 + 1
      else
        kstrt = k1
      end if
!
!      if (emiss_aero_opt /= 0) then
!        emiss_aero(i1:i2,ju1:j2,1:naero) =  &
!     &    emiss_aero_t(i1:i2,ju1:j2,1:naero,it)
!      end if
!
      if (emiss_aero_opt == 1) then                        ! GMI    aerosol emissions
         emiss_aero(i1:i2,ju1:j2,1:naero) = emiss_aero_t(i1:i2,ju1:j2,1:naero,it)
      elseif (emiss_aero_opt == 2) then                    ! GOCART aerosol emissions
         emiss_aero(i1:i2,ju1:j2,1:5)     = emiss_aero_t(i1:i2,ju1:j2,1:5,it)
      end if
!
!      inum = 0
!
!     ================================
      SPCLOOP: do icx = 1, num_emiss
!     ================================
!
        ic = emiss_map(icx)
!
        if (ic > 0) then
!
!          inum = inum + 1
!
          if ((chem_opt == 8) .and.  &
     &        ((ic == IFSO2) .or. (ic == INSO2) .or.  &
     &         (ic == INDMS))) then
!           ------------------------------------------------
!           For sulfur chemistry, do emissions in chemistry.
!           ------------------------------------------------
!
!           =============
            cycle SPCLOOP
!           =============
          end if
!
!
          mw_fac = MWTAIR / mw(ic)
!
!.sds... if doing in-chem surface emissions, must be active species, otherwise emit here
          if ( ic .gt. SK_NACT .and. do_semiss_inchem) then
            where (emissionArray(icx)%pArray3D(:,:,k1) /= 0.0d0)
               emass(:,:,k1) = emissionArray(icx)%pArray3D(:,:,k1) * tdt
!
               emass(:,:,k1) = (emass(:,:,k1) / mass (:,:,k1)) * mw_fac
!
               concentration(ic)%pArray3D(:,:,k1) =  &
     &           concentration(ic)%pArray3D(:,:,k1) + emass(:,:,k1)
            end where
          endif
!.sds... if doing in-chem surface emissions, kstrt=2 otherwise = 1
          where (emissionArray(icx)%pArray3D(:,:,kstrt:k2) /= 0.0d0)
            emass(:,:,kstrt:k2) =  &
     &        emissionArray(icx)%pArray3D(:,:,kstrt:k2) * tdt
!
            emass(:,:,kstrt:k2) =  &
     &        (emass(:,:,kstrt:k2) /  &
     &         mass (:,:,kstrt:k2)) * mw_fac
!
            concentration(ic)%pArray3D(:,:,kstrt:k2) =  &
     &        concentration(ic)%pArray3D(:,:,kstrt:k2) +  &
     &        emass(:,:,kstrt:k2)
!
          end where
!
!
          if ((pr_surf_emiss) .and. (.not. do_semiss_inchem))  &
     &      surf_emiss_out(:,:,ic) = surf_emiss_out(:,:,ic) +  &
     &        (emissionArray(icx)%pArray3D(:,:,1) * tdt / mcor(:,:))
!
          if (pr_emiss_3d) then
!.sds          if ((pr_emiss_3d) .and. (.not. do_semiss_inchem)) then
            do ij = ju1, j2
              do il = i1, i2
                do ik = kstrt, k2
                  if (emissionArray(icx)%pArray3D(il,ij,ik) /= 0.0d0) then
                      emiss_3d_out(il,ij,ik,ic) = emiss_3d_out(il,ij,ik,ic) +  &
     &                   (emissionArray(icx)%pArray3D(il,ij,ik)*tdt/mcor(il,ij) )
                  end if
                end do
                if (kstrt /= k1) then
                  emiss_3d_out(il,ij,1,ic) = emiss_3d_out(il,ij,1,ic) +  &
     &                   (emissionArray(icx)%pArray3D(il,ij,1)*tdt/mcor(il,ij) )
                end if
              end do
            end do
          end if
!
!                                ============
!          if (inum == num_emiss) exit SPCLOOP
!                                ============
!
        end if
!
!
!
!     ==============
      end do SPCLOOP
!     ==============
!
!
!     =====================================================
      if ((emiss_aero_opt /= 0) .or. (emiss_dust_opt /= 0)) then
!     =====================================================
!
!c?
!       -------------------------------------------------------------------
!       Add dust emission and other aerosols to some lowest vertical levels
!       (emitted uniformly below the PBL).
!       -------------------------------------------------------------------
!
!       --------------------------------------------------------------------
!       Compute the height of full and half-sigma levels above ground level.
!       --------------------------------------------------------------------
!
        zq(i1:i2, ju1:j2, k1-1) = 0.0d0
!
!        do ik = k1, k2
!          zq(i1:i2,ju1:j2,ik) =  zq(i1:i2,ju1:j2,ik-1) + gridBoxHeight(i1:i2,ju1:j2,ik)
!        end do
!
        do ik = k1, k2
          zq(i1:i2,ju1:j2,ik) =  zq(i1:i2,ju1:j2,ik-1) + gridBoxHeight(i1:i2,ju1:j2,ik)
          za(i1:i2,ju1:j2,ik) = 0.5d0 * (zq(i1:i2,ju1:j2,ik) + zq(i1:i2,ju1:j2,ik-1))
        end do
!
!       ---------------------------------------------
!       Calculate total air mass (kg) within the PBL.
!       ---------------------------------------------
!
        mass_pbl(:,:) = mass(:,:,1)
!
        do ik = k1+1, k2
!
          where (za(:,:,ik) <= pbl(:,:))
            mass_pbl(:,:) = mass_pbl(:,:) + mass(:,:,ik)
          end where
!
        end do
!
!       ----------------------------------------------------------------
!       Add dust and aerosol emissions (kg/s) to concentrations; unit of
!       aerosol and dust concentrations is kg/(kg air).
!       ----------------------------------------------------------------
!
!!!!! GOCART Source for dust
        if (emiss_dust_opt == 2) then
           call SourceDust (mcor, tdt, nymd, i1, i2, ju1, j2, k1, k2)
           emiss_dust(i1:i2,ju1:j2,1:ndust) = srcEmissDust(i1:i2,ju1:j2,1:ndust)
           emiss_map_dust(1) = IDUST1
           emiss_map_dust(2) = IDUST2
           emiss_map_dust(3) = IDUST3
           emiss_map_dust(4) = IDUST4
           emiss_map_dust(5) = IDUST5
        end if
!!!!! GOCART Source for sea salt
        if (emiss_aero_opt == 2) then
           call SourceSeaSalt (mcor, tdt, nymd, i1, i2, ju1, j2, k1, k2)
           emiss_aero(i1:i2,ju1:j2,6:naero) = srcEmissSeaSalt(i1:i2,ju1:j2,1:4)
           emiss_map_aero(6) = ISSLT1
           emiss_map_aero(7) = ISSLT2
           emiss_map_aero(8) = ISSLT3
           emiss_map_aero(9) = ISSLT4
        end if
!!!!!
!
        do ij = ju1, j2
          do il = i1, i2
!
            IKLOOP2: do ik = k1, k2
!
              if ((za(il,ij,ik) > pbl(il,ij)) .and. (ik /= k1)) then
!               ============
                exit IKLOOP2
!               ============
              end if
!
!             -----
!             Dust.
!             -----
!
              if (emiss_dust_opt /= 0) then
                do icx = 1, ndust
!
                  ic = emiss_map_dust(icx)
!
! Bigyani's correction
!                 if (mass_pbl(il,ij) .le. 1.e-20)then
!                    mass_pbl(il,ij) = 1.e-10
!                 endif
! End Bigyani's correction
!
                  if (do_aerocom)then
                     if ((ic /= INDMS) .and. (ik == k1)) then
                       concentration(ic)%pArray3D(il,ij,ik) =  &
     &                   concentration(ic)%pArray3D(il,ij,ik) +  &
     &                     emiss_dust(il,ij,icx) * tdt/ mass(il,ij,ik)
                     endif
                  else
                     concentration(ic)%pArray3D(il,ij,ik) =  &
     &                 concentration(ic)%pArray3D(il,ij,ik) +  &
     &                 emiss_dust(il,ij,icx) * tdt / mass_pbl(il,ij)
                  end if
!
                end do
              end if
!
!             ------------------------------
!             Other aerosol (carbon & sslt).
!             ------------------------------
!
              if (emiss_aero_opt /= 0) then
                do icx = 1, naero
!
                  ic = emiss_map_aero(icx)
!
                  if ((ic == IBOC) .or. (ic == IBBC)) then
!
                    concentration(ic)%pArray3D(il,ij,ik) =  &
     &                concentration(ic)%pArray3D(il,ij,ik) +  &
     &                emiss_aero(il,ij,icx) * tdt / mass_pbl(il,ij)
!
                  else if (((ic == INOC)   .or. (ic == IFOC)   .or.  &
     &                      (ic == IFBC)   .or. (ic == ISSLT1) .or.  &
     &                      (ic == ISSLT2) .or. (ic == ISSLT3) .or.  &
     &                      (ic == ISSLT4)) .and.  &
     &                     (ik == k1)) then
!
                    concentration(ic)%pArray3D(il,ij,ik) =  &
     &                concentration(ic)%pArray3D(il,ij,ik) +  &
     &                emiss_aero(il,ij,icx) * tdt / mass(il,ij,ik)
!
                  end if
!
                end do
              end if
!
            end do IKLOOP2
!
          end do
        end do
!
!     ======
      end if
!     ======
!
!
      return
!
      end
!
!
