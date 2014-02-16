!=============================================================================
!
!
! CODE DEVELOPER
!   Adapted from Do_Convec_Ncar routine
!     by Stephen Steenrod, SSAI, Feb 2005
!
! FILE
!   convec_gmao.F
!
! ROUTINES
!   Do_Convec_Geos4
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Convec_GEOS4
!
! DESCRIPTION
!   This is the interface routine to Conv_Tran.  It formats the GEOS4 variables
!   to satisfy Conv_Tran.
!
! ARGUMENTS
!   metdata_name_org   : first  part of metdata_name, e.g., "GMAO"
!   metdata_name_model : second part of metdata_name, e.g., "GEOS4"
!   det_ent      : flag for doing detrainment then entrainment
!   do_downdraft : flag for doing downdrafts
!   do_wetdep    : flag for doing wet deposition
!   pr_wet_depos : flag for saving and then printing accumulated wet
!                  deposition
!   chem_opt     : chemistry option
!   ih2o2_num    : const array index for H2O2
!   ihno3_num    : const array index for HNO3
!   lwi_flags    : land water ice flags (1-water 2-land 3-ice)
!   tdt          : model time step  (s)
!   mw           : array of species' molecular weights (g/mol)
!   mcor         : horizontal area of each grid box    (m^2)
!   pbl          : planetary boundary layer thickness  (m)
!... Z-M is for Zheng-McFarland convection (deep)
!   zmdu         : detrainment mass flux from Z-M convective updraft (Pa/s)
!   zmeu         : entrainment mass flux into Z-M convective updraft (Pa/s)
!   zmed         : entrainment mass flux into Z-M convective downdraft (Pa/s)
!   zmmd         : Z-M convective mass flux in downdraft (Pa/s)
!   zmmu         : Z-M upward mass flux (Pa/s)
!... HK is for Hack convection (shallow)
!   hkdu         : detrainment mass flux from Hack convective updraft (Pa/s)
!   hkeu         : entrainment mass flux into Hack convective updraft (Pa/s)
!   hkmu         : Hack upward mass flux (Pa/s)
!
!   grid_height  : grid box height  (m)
!   mass         : mass of air in each grid box (kg)
!   wet_depos    : accumulated wet deposition   (kg/m^2)
!   kel          : temperature      (degK)
!   press3e      : atmospheric pressure at the edge of each grid box (mb)
!   const        : species concentration, known at zone centers (mixing ratio)
!-----------------------------------------------------------------------------
!-micro_aerosol  MICRO_AEROSOL
!   humidity     : specific humidity (g/kg)
!   press3c      : atmospheric pressure at the center of each grid box (mb)
!   REL_SCAV_EFF_new
!                : relative scavenging efficiency for aerosols (0-1)
!
!-----------------------------------------------------------------------------

      subroutine Do_Convec_GEOS4  &
     &  (chem_mecha, metdata_name_org, metdata_name_model, det_ent,  &
     &   do_downdraft, do_wetdep,  &
     &   pr_wet_depos, chem_opt, ih2o2_num, ihno3_num, lwi_flags, tdt,  &
     &   mw, mcor, pbl, zmdu, zmeu, zmed, zmmd, zmmu, hkdu, hkeu, hkmu, &
     &   grid_height, mass, wet_depos, kel, press3e, concentration,  &
#ifdef MICRO_AEROSOL
     &   humidity, press3c, REL_SCAV_EFF_new, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ilong, num_species)
#else
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ilong, num_species)
#endif

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiSpcConcentrationMethod_mod, only : isFixedConcentration
      use GmiPrintError_mod, only : GmiPrintError

      implicit  none

#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"

#ifdef MICRO_AEROSOL
#     include "gmi_micro_aerosol.h"
#elif GOCARTaerosol
#     include "gocart_aerosol.h"
#elif strat_trop_aerosol
#     include "gocart_aerosol.h"
#elif strat_trop
#     include "gocart_aerosol.h"
#else
#     include "gmi_aerosol.h"
#endif

!c?   Tight coupling to setkin?
#     include "setkin_par.h"
#     include "setkin_depos.h"
#ifdef MICRO_AEROSOL
#     include "gmi_time_constants.h"
#     include "umaerosol.h"
#endif

#ifdef nonZeroInd
      integer, parameter :: IFSO2_l = 1            
      integer, parameter :: INSO2_l = 1
#elif nonZeroInd_tracers
      integer, parameter :: IFSO2_l = 1            
      integer, parameter :: INSO2_l = 1
#else
      integer, parameter :: IFSO2_l = IFSO2        
      integer, parameter :: INSO2_l = INSO2
#endif


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in) :: ilong, num_species

      character (len=* ) :: chem_mecha
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_model
      logical :: det_ent
      logical, intent(in   ) :: do_downdraft
      logical, intent(in   ) :: do_wetdep
      logical, intent(in   ) :: pr_wet_depos
      integer, intent(in   ) :: chem_opt
      integer, intent(in   ) :: ih2o2_num
      integer, intent(in   ) :: ihno3_num
      integer, intent(in   ) :: lwi_flags(i1:i2, ju1:j2)
      real*8 , intent(in   ) :: tdt
      real*8 , intent(in   ) :: mw(num_species)
      real*8 , intent(in   ) :: mcor       (i1:i2,   ju1:j2)
      real*8 , intent(in   ) :: pbl        (i1:i2,   ju1:j2)
      real*8 , intent(in   ) :: zmdu       (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: zmeu       (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: zmed       (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: zmmd       (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: zmmu       (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: hkdu       (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: hkeu       (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: hkmu       (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: grid_height(i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: mass       (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(inout) :: wet_depos  (i1:i2,   ju1:j2,   num_species)
      real*8 , intent(in   ) :: kel        (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in   ) :: press3e    (ilo:ihi, julo:jhi, k1-1:k2)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)

#ifdef MICRO_AEROSOL
      real*8 , intent(in   ) :: humidity   (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: press3c    (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(out  ) :: REL_SCAV_EFF_new(i1:i2, ju1:j2, k1:k2, num_species)
#endif

!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter ::  &
     &  MBSTH = 1.0d-15      ! threshold below which we treat mass fluxes
                             ! as zero (mb/s)


!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      integer :: iku
      integer :: il, ij, ik, ic
      integer :: il2g                      ! gathered index to operate over
      integer :: itdt_conv
      integer :: num_zmconv_steps, num_hkconv_steps

      real*8  :: tdt_zmconv, tdt_hkconv

      real*8  :: xmbsth

      integer :: ideep(ilong)              ! gathering array
      integer :: pbli (ilong)              ! index of pbl height

      real*8  :: col_sum_depos   (i1:i2)   ! sum of column wet_deposition loss
      real*8  :: kloss           (i1:i2)   ! wet loss rate constant at all
                                           ! longitudes (s^-1)
      real*8  :: updraft_velocity(i1:i2)   ! velocity in convective updraft (m/s)

      real*8  :: dpi(ilong, k1:k2)     ! delta pressure between interfaces

      real*8  :: zmdui(ilong, k1:k2)   ! mass flux detraining from updraft (Pa/s)
      real*8  :: zmeui(ilong, k1:k2)   ! mass flux entraining into updraft (Pa/s)
      real*8  :: zmedi(ilong, k1:k2)   ! mass flux entraining into downdraft (Pa/s)
      real*8  :: zmmdi(ilong, k1:k2)   ! mass flux down (Pa/s)
      real*8  :: zmmui(ilong, k1:k2)   ! mass flux up (Pa/s)

      real*8  :: hkmui(ilong, k1:k2)   ! mass flux up (Pa/s)
      real*8  :: hkdui(ilong, k1:k2)   ! mass flux detraining from updraft (Pa/s)
      real*8  :: hkeui(ilong, k1:k2)   ! mass flux entraining into updraft (Pa/s)
      real*8  :: zero(ilong, k1:k2)    ! array of zeros for Hack downdrafts

      real*8  :: fracis(i1:i2, k1:k2, num_species)  ! insoluble fraction of tracer
      real*8  :: qq    (i1:i2, k1:k2, num_species)  ! tracer array including moisture

#ifdef MICRO_AEROSOL
!-micro_aerosol----begin------------------------------------------------------
!     --------------------
!     adding for umaerosol
!     --------------------

      real*8  :: relhume   (i1:i2, ju1:j2, k1:k2)
      real*8  :: h2osat    (i1:i2, ju1:j2, k1:k2)
      real*8  :: h2ogas    (i1:i2, ju1:j2, k1:k2)
      real*8  :: so4mfrac  (i1:i2, ju1:j2, k1:k2)
      real*8  :: so4dens   (i1:i2, ju1:j2, k1:k2)
      real*8  :: wetmas    (i1:i2, ju1:j2, k1:k2)
      real*8  :: so4aer    (i1:i2, ju1:j2, k1:k2, naer)
      real*8  :: so4radv   (i1:i2, ju1:j2, k1:k2, nso4)
#ifndef NONON
      real*8  :: so4non    (i1:i2, ju1:j2, k1:k2, nnon)
      real*8  :: aernon    (i1:i2, ju1:j2, k1:k2, nnon)
      real*8  :: xnonum    (i1:i2, ju1:j2, k1:k2, nnon)
      real*8  :: xso4non   (i1:i2, ju1:j2, k1:k2, nnon)
      real*8  :: fso4non   (i1:i2, ju1:j2, k1:k2, nnon)
      real*8  :: rso4min   (i1:i2, ju1:j2, k1:k2)
#endif
      real*8  :: xcldnum   (i1:i2, ju1:j2, k1:k2)
      real*8  :: so4cldcoag(i1:i2, ju1:j2, k1:k2)

!-micro_aerosol----end--------------------------------------------------------
#endif

!     ----------------
!     Begin execution.
!     ----------------

      ideep(:) = 0
      pbli (:) = 0

      kloss(:) = 0.0d0

      updraft_velocity(:) = 0.0d0

      dpi(:,:) = 0.0d0

      zmdui(:,:) = 0.0d0
      zmeui(:,:) = 0.0d0
      zmedi(:,:) = 0.0d0
      zmmui(:,:) = 0.0d0
      zmmdi(:,:) = 0.0d0

      hkdui(:,:) = 0.0d0
      hkeui(:,:) = 0.0d0
      hkmui(:,:) = 0.0d0
      zero(:,:) = 0.0d0

      fracis(:,:,:) = 1.0d0


      xmbsth = MBSTH


!cc!? For now just hardwire a few values for the radon/lead problem
!c    (chem_opt = 1).  Setkin will eventually provide this information
!c    for all species when doing a full chemistry calculation.

      if (chem_opt == 1) then

        aerosol(1) = 0
        aerosol(2) = 1

        hstar(1)   = 9.3d-3
        hstar(2)   = 0.0d0

        delh_298_over_r(1) = 2600.0d0
        delh_298_over_r(2) =    0.0d0

        retention_eff  (1) =    0.0d0
        retention_eff  (2) =    0.0d0

      else if (chem_opt == 6) then

        aerosol(1) = 15
        aerosol(2) = 15

        hstar(1)   = 0.0d0
        hstar(2)   = 0.0d0

        delh_298_over_r(1) =    0.0d0
        delh_298_over_r(2) =    0.0d0

        retention_eff  (1) =    0.0d0
        retention_eff  (2) =    0.0d0

      else if (chem_opt == 8) then

        hstar(IFSO2_l) = 600.0d0
        hstar(INSO2_l) = 600.0d0

      endif

!.sds. make like gocart
       if(IFSO2 .gt. 0) hstar(IFSO2_l) = 600.0d0
       if(INSO2 .gt. 0) hstar(INSO2_l) = 600.0d0
!.sds. 


!     -----------------------------------------------------------
!     Calculate any needed sub-cycling of the convection operators
!     by comparing the mass flux over the full time step to the
!     mass within the grid box.
!     -----------------------------------------------------------

      num_zmconv_steps =  &
     &  Maxval (tdt * zmmu(:,:,:) * Spread(mcor(:,:),3,k2-k1+1) /  &
     &          mass(:,:,:)) + 1.0d0

      tdt_zmconv       = tdt / real(num_zmconv_steps)

      num_hkconv_steps =  &
     &  Maxval (tdt * hkmu(:,:,:) * Spread (mcor(:,:), 3, k2-k1+1) /  &
     &          mass(:,:,:)) + 1.0d0

      tdt_hkconv       = tdt / real(num_hkconv_steps)

#ifdef MICRO_AEROSOL
!-micro_aerosol----begin------------------------------------------------------
!     if(chem_mecha == 'micro_aerosol') then
!     --------------------
!     adding for umaerosol
!     --------------------

!     compute relative humidity (0-1)

      relhume(:,:,:) = 1.0d0 - (373.15d0 / kel(i1:i2,ju1:j2,k1:k2))
      relhume(:,:,:) =  &
     &  1013.25d0 * Exp (13.3185d0 * relhume(:,:,:)    -  &
     &                    1.9760d0 * relhume(:,:,:)**2 -  &
     &                    0.6445d0 * relhume(:,:,:)**3 -  &
     &                    0.1299d0 * relhume(:,:,:)**4)

      relhume(:,:,:) =  &
     &  humidity(:,:,:) * MWTAIR / 18.0d0 /  &
     &  GPKG * press3c(i1:i2,ju1:j2,k1:k2) / relhume(:,:,:)

      relhume(:,:,:) = Max (Min (relhume(:,:,:), 0.95d0), 0.0d0)

!c    compute so4 aerosol size (m)

      so4aer(:,:,:,1) = concentration(ISO4M1)%pArray3D(:,:,:)
      so4aer(:,:,:,2) = concentration(ISO4N1)%pArray3D(:,:,:)
      so4aer(:,:,:,3) = concentration(ISO4M2)%pArray3D(:,:,:)
      so4aer(:,:,:,4) = concentration(ISO4N2)%pArray3D(:,:,:)

      do il = i1, i2
        do ij = ju1, j2
          do ik = k1, k2
            h2osat(il,ij,ik) = h2osat_f(kel(il,ij,ik))
            h2ogas(il,ij,ik) = relhume(il,ij,ik)*h2osat(il,ij,ik)

            so4mfrac(il,ij,ik) = so4mfrac_f(kel(il,ij,ik),  &
     &                                      h2osat(il,ij,ik),  &
     &                                      h2ogas(il,ij,ik))
            so4dens(il,ij,ik) = so4dens_f(kel(il,ij,ik),  &
     &                                    so4mfrac(il,ij,ik))
          end do
        end do
      end do

      do iso4 = 1, nso4

        iso4n = iso4 * nmomso4
        iso4m = iso4n - 1

        wetmas(:,:,:) = max(r2so4min, so4aer(:,:,:,iso4m)  &
     &                / max(epsilo,   so4aer(:,:,:,iso4n)))  &
     &                / so4mfrac(:,:,:)

        so4radv(:,:,:,iso4) = (r3q * wetmas(:,:,:)  &
     &                      / (GMI_PI * so4dens(:,:,:)))**r1td
        so4radv(:,:,:,iso4) = max(so4radvmin,  &
     &                        min(so4radv(:,:,:,iso4),so4radvmax))

      end do

#ifndef NONON
      do ic = 1, nnon
         ! from ISO4NOC to ISO4S4
         so4non(:,:,:,ic) = concentration(ISO4NOC+ic-1)%pArray3D(:,:,:)

         ! from INOC to ISSLT4
         aernon(:,:,:,ic) = concentration(INOC   +ic-1)%pArray3D(:,:,:)
      end do
!      so4non(:,:,:,1:nnon) = concentration(ISO4NOC:ISO4S4)%pArray3D(:,:,:)
!      aernon(:,:,:,1:nnon) = concentration(INOC:ISSLT4)%pArray3D(:,:,:)

!c    non-so4 aerosol radius (m) and mass (kg/particle)
      do inon = 1, nnon

        radvolm(:) = radgnon(:,inon)  &
     &             * exp(r3h*log(siggnon(:,inon))**2)
        radvnon(inon) = sum(fracnon(:,inon)*radvolm(:)**3)
        pmsnon(inon) = r4td*GMI_PI*rhonon(inon)  &
     &               * radvnon(inon)
        radvnon(inon) = radvnon(inon)**r1td

      end do

!c    non-so4 aerosol number concentration (#particles/kg air)
      do inon = 1, nnon
        xnonum(:,:,:,inon) = aernon(:,:,:,inon) / pmsnon(inon)
      end do

!c    so4 molecule # on non-so4 aerosol surface (#molec/particle)
      xso4non(:,:,:,1:nnon) = max(1.0d-30,so4non(:,:,:,1:nnon)) /  &
     &  max(r1,xnonum(:,:,:,1:nnon)) / so4min

!c    radius of each so4 molecule (m)
      rso4min(:,:,:) =  &
     &  (so4min * r3q /  &
     &  (GMI_PI*so4mfrac(:,:,:)*so4dens(:,:,:)))**r1td

!c    area fraction of so4 molecules on non-so4 aerosol surface
      do inon = 1, nnon

        fso4non(:,:,:,inon) =  &
     &    rso4min(:,:,:)**2.0d0 * xso4non(:,:,:,inon)  &
     &    / (4.0d0*(radvnon(inon)+rso4min(:,:,:))**2.0d0)

      end do
#endif

!c    compute cloud droplet number (cm^-3) to be the sum of
!c    mode 2 so4 number and non-so4 aerosol

      xcldnum(:,:,:) = so4aer(:,:,:,4)

#ifndef NONON
      do inon = 1, nnon

        if (inon <= 9) then                  ! 5 oc/bc, dust 1-4
          xcldnum(:,:,:) = xcldnum(:,:,:) + xnonum(:,:,:,inon)  &
     &                   * min(fso4non(:,:,:,inon)/fso4crit, r1)
        else                                 ! ss 1-4
          xcldnum(:,:,:) = xcldnum(:,:,:) + xnonum(:,:,:,inon)
        end if

      end do
#endif

!c    change unit from kg^-1 to cm^-3
      xcldnum(:,:,:) = min(r3000,max(r10,xcldnum(:,:,:) *  &
     &  (press3c(i1:i2,ju1:j2,k1:k2) * MWTAIR * MB2CGS) /  &
     &  (GPKG * kel(i1:i2,ju1:j2,k1:k2) * BOLTZMN_E * AVOGAD)))

!c    compute Brownian coagulation coefficient K12 (cm3 s-1)
!c    for the 1st so4 mode with cloud droplets
!c    fit from Table 12.3 and Fig.12.5, in Seinfeld and Pandis (1997)
!c    as log10(K12) = -1.8436 * log10(D(um)) - 9.21

      so4radv(:,:,:,1) = min(so4radv(:,:,:,1),1.0d-7)

      so4cldcoag(:,:,:) =  &
     &  10 ** (a11 * log10(so4radv(:,:,:,1)*2.0d+6) + a12)


!     ==============================
      ICLOOP: do ic = 1, num_species
!     ==============================

!c      calculate scavenging efficiency of aerosols

        if (isFixedConcentration(ic) .or.  &
     &     (aerosol(ic) <= 0)  .or.  &
     &     (aerosol(ic) >  NUM_AEROSOL)) then
!         ============
          cycle ICLOOP
!         ============
        end if

!c      SO4 mode 1 (assume cloud lifetime 4 hr)
        if (aerosol(ic) == 1) then

          REL_SCAV_EFF_new(:,:,:,ic) = min(r1,  &
     &      so4cldcoag(:,:,:) * 4.0d0 * SECPHR * xcldnum(:,:,:))

        else
#ifdef NONON
          REL_SCAV_EFF_new(:,:,:,ic) = REL_SCAV_EFF(aerosol(ic))
#else
          inon = aerosol(ic) - nso4

          if (inon>0 .and. inon<=9) then         ! 5 oc/bc, dust1-4

            REL_SCAV_EFF_new(:,:,:,ic) =  &
     &        min(fso4non(:,:,:,inon)/fso4crit, r1)

          else                             ! mode 2 so4, ss 1-4

            REL_SCAV_EFF_new(:,:,:,ic) = REL_SCAV_EFF(aerosol(ic))

          end if
#endif
        end if

!     =============
      end do ICLOOP
!     =============
!-micro_aerosol----end--------------------------------------------------------
!     endif     ! chem_mecha == 'micro_aerosol'
#endif


!     =======
      ijloop: do ij = ju1, j2
!     =======

!... initialize (flip) lat circle of constituents for convective trans
        do ic = 1, num_species
           qq(:,k2:k1:-1,ic) = concentration(ic)%pArray3D(:,ij,k1:k2)
        end do

!... do Zhang-McFarland (deep) convection first
        il2g = 0

        do il = i1, i2

!... Verify that there is convection in the current cell.
!          if (Maxval (zmmu(il,ij,:)) .gt. 0.0d0) then
          if (Maxval (zmmu(il,ij,:)) .gt. 0.0d0 .or.  &
     &        Minval (zmmd(il,ij,:)) .lt. 0.0d0) then

            il2g = il2g + 1

            ideep(il2g) = il

            dpi(il2g,:) = (press3e(il,ij,k2-1:k1-1:-1) -  &
     &         press3e(il,ij,k2:k1:-1)) * PASPMB

!... Zhang-McFarland convection - deep
            zmmui(il2g,:) = zmmu(il,ij,k2:k1:-1) * GMI_G

            zmeui(il2g,:) = zmeu(il,ij,k2:k1:-1) * GMI_G
!... should this be exactly consistant with zmmu and zmeu?
            zmdui(il2g,:) = zmdu(il,ij,k2:k1:-1) * GMI_G

            if(do_downdraft) then
              zmmdi(il2g,:) = zmmd(il,ij,k2:k1:-1) * GMI_G
              zmedi(il2g,:) = zmed(il,ij,k2:k1:-1) * GMI_G
            else
              zmmdi(il2g,:) = 0.0
              zmedi(il2g,:) = 0.0
            endif

          endif

        enddo
!=============
!       =======
        il2gifzm: if (il2g .gt. 0) then
!       =======

          if (do_wetdep) then

!           ----------------------------------------------------------------
!           Calculate the insoluble fraction at each level for each species.
!           ----------------------------------------------------------------

            where ((lwi_flags(i1:i2,ij) .eq. 1) .or.  &
     &             (lwi_flags(i1:i2,ij) .eq. 3))
              updraft_velocity(:) =  5.0d0
            elsewhere
              updraft_velocity(:) = 10.0d0
            end where

            do ik = k1, k2

              iku = k2 - ik + 1

              do ic = 1, num_species

#ifdef MICRO_AEROSOL
                if (ic == ihno3_num .or. ic == ISO4G) then
#else
                if (ic == ihno3_num) then
#endif

                  kloss(:) = 5.0d-3

                else if (aerosol(ic) .ge. 1) then


#ifdef MICRO_AEROSOL
                  kloss(:) = 5.0d-3  &
     &                     * REL_SCAV_EFF_new(:,ij,ik,ic)
#else
                  kloss(:) = 5.0d-3 * REL_SCAV_EFF(aerosol(ic))
#endif

                else if (hstar(ic) .gt. 0.0d0) then

!                 =======================
                  call Calc_Wet_Loss_Rate  &
!                 =======================
     &              (.true., ic, ih2o2_num, delh_298_over_r(ic),  &
     &               hstar(ic), retention_eff(ic), press3e(:,ij,ik),  &
     &               kel(:,ij,ik), kloss(:), i1, i2, ilo, ihi)

                  kloss(:) = 5.0d-3 * kloss(:)

                else

                  kloss(:) = 0.0d0

                endif

                fracis(:,iku,ic) =  &
     &            Exp (-kloss(:) * grid_height(:,ij,ik) /  &
     &                 updraft_velocity(:))

              enddo

            enddo

          else

           fracis(:,:,:) = 1.0d0

          endif


!         ----------------------------------------------------------------
!         Find the index of the top of the planetary boundary layer (pbl).
!         Convection will assume well mixed tracers below that level.
!         ----------------------------------------------------------------

          do il = 1, il2g

            pbli(il) = 0

            ikloopzm: do ik = k1, k2

              if (pbl(ideep(il),ij) .lt.  &
     &            Sum (grid_height(ideep(il),ij,k1:ik))) then
                pbli(il) = ik
                exit ikloopzm
              endif

            enddo ikloopzm

            if (pbli(il) .eq. 0) then

              err_msg = 'Could not find pbl in Do_Convec_GEOS4'
              call GmiPrintError (err_msg, .true., 2, ideep(il), ij,  &
     &                         1, pbl(ideep(il),ij), 0.0d0)

            endif

            pbli(il) = k2 - pbli(il)

          enddo

!... Zheng-McFarland (deep) convection
          do itdt_conv = 1, num_zmconv_steps

            call Conv_Tran  &
     &        (il2g, tdt_zmconv, xmbsth, ideep, pbli, zmdui, zmeui,  &
     &         zmmui, zmmdi, dpi, fracis, qq, i1, i2, k1, k2, ilong, num_species)

          enddo

!       =============
        endif il2gifzm
!       =============



!... do Hack (shallow) convection second

        il2g = 0

        do il = i1, i2

!... Verify that there is convection in the current cell.
          if (Maxval (hkmu(il,ij,:)) .gt. 0.0d0) then

            il2g = il2g + 1

            ideep(il2g) = il

            dpi(il2g,:) =  &
     &        (press3e(il,ij,k2-1:k1-1:-1) -  &
     &         press3e(il,ij,k2:k1:-1)) *  &
     &        PASPMB


!... compress Hack convection - shallow
            hkmui(il2g,:) = hkmu(il,ij,k2:k1:-1) * GMI_G

            hkeui(il2g,:) = hkeu(il,ij,k2:k1:-1) * GMI_G
!... should this be exactly consistant with hkmu and hkeu?
            hkdui(il2g,:) = hkdu(il,ij,k2:k1:-1) * GMI_G

          endif

        enddo
!=============
!... Now do Hack Convection - shallow
!       =======
        il2gifhk: if (il2g /= 0) then
!       =======
!c
!c... Not doing wet scaveging in the Hack convection, so set fracis = 1
!c
!          if (do_wetdep) then
!
!c           ----------------------------------------------------------------
!c           Calculate the insoluble fraction at each level for each species.
!c           ----------------------------------------------------------------
!
!            where ((lwi_flags(i1:i2,ij) == 1) .or.
!     &             (lwi_flags(i1:i2,ij) == 3))
!              updraft_velocity(:) =  5.0d0
!            elsewhere
!              updraft_velocity(:) = 10.0d0
!            end where
!
!            do ik = k1, k2
!
!              iku = k2 - ik + 1
!
!              do ic = 1, num_species
!
!                if (ic == ihno3_num) then
!
!                  kloss(:) = 5.0d-3
!
!                else if (aerosol(ic) >= 1) then
!
!                  kloss(:) = 5.0d-3 * REL_SCAV_EFF(aerosol(ic))
!
!                else if (hstar(ic) > 0.0d0) then
!
!c                 =======================
!                  call Calc_Wet_Loss_Rate
!c                 =======================
!     &              (.true., ic, ih2o2_num, delh_298_over_r(ic)
!     &               , hstar(ic), retention_eff(ic), press3e(:,ij,ik)
!     &               , kel(:,ij,ik), kloss(:), i1, i2, ilo, ihi)
!
!                  kloss(:) = 5.0d-3 * kloss(:)
!
!                else
!
!                  kloss(:) = 0.0d0
!
!                endif
!
!                fracis(:,iku,ic) =
!     &            Exp (-kloss(:) * grid_height(:,ij,ik) /
!     &                 updraft_velocity(:))
!
!              enddo
!
!            enddo
!
!          else
!
           fracis(:,:,:) = 1.0d0
!
!          endif


!         ----------------------------------------------------------------
!         Find the index of the top of the planetary boundary layer (pbl).
!         Convection will assume well mixed tracers below that level.
!         ----------------------------------------------------------------

          do il = 1, il2g

            pbli(il) = 0

            ikloophk: do ik = k1, k2

              if (pbl(ideep(il),ij) <  &
     &            Sum (grid_height(ideep(il),ij,k1:ik))) then
                pbli(il) = ik
                exit ikloophk
              endif

            enddo ikloophk

            if (pbli(il) == 0) then

              err_msg = 'Could not find pbl in Do_Convec_GEOS4.'
              call GmiPrintError (err_msg, .true., 2, ideep(il), ij,  &
     &                         1, pbl(ideep(il),ij), 0.0d0)

            endif

            pbli(il) = k2 - pbli(il)

          enddo
!         ==========
!... Hack convection - shallow
          do itdt_conv = 1, num_hkconv_steps
!... no downdrafts in Hack convection
            zero(:,:) = 0.0d0

            call Conv_Tran  &
     &        (il2g, tdt_hkconv, xmbsth, ideep, pbli, hkdui, hkeui,  &
     &         hkmui, zero, dpi, fracis, qq, i1, i2, k1, k2, ilong, num_species)

          enddo

!       =============
        endif il2gifhk
!       =============

        if (pr_wet_depos) then

!         ---------------------------------------------------
!         Calculate the wet deposition from the difference in
!         mixing ratio.
!         ---------------------------------------------------

          do ic = 1, num_species

            if ((hstar(ic) > 0.0d0) .or.  &
     &          (aerosol(ic) >= 1)  .or.  &
     &          (ic == ihno3_num)) then

              col_sum_depos(i1:i2) = 0.0d0

              do ik = k1, k2

                iku = k2 - ik + 1

                col_sum_depos(:) =  &
     &            col_sum_depos(:) +  &
     &            (concentration(ic)%pArray3D(:,ij,ik) - qq(:,iku,ic)) *  &
     &            mass(:,ij,ik)

              enddo

              where (col_sum_depos(i1:i2) > 0.0d0)  &
     &          wet_depos(i1:i2,ij,ic) =  &
     &            wet_depos(i1:i2,ij,ic) +  &
     &            col_sum_depos(i1:i2) * mw(ic) /  &
     &            MWTAIR / mcor(i1:i2,ij)

            endif

          enddo

        endif

!... put back into const (and flip)
        do ic = 1, num_species
           concentration(ic)%pArray3D(:,ij,k1:k2) = qq(:,k2:k1:-1,ic)
        end do

!     =============
      enddo ijloop
!     =============

      return

      end subroutine Do_Convec_GEOS4
