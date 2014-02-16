!=============================================================================
!
! $Id: convec_ncar.F90,v 1.21 2012-01-05 21:41:30 jkouatch Exp $
!
! CODE DEVELOPER
!   Original code from Phil Rasch, NCAR.
!   LLNL modifications:  Dan Bergmann
!                          dbergmann@llnl.gov
!                        John Tannahill
!                          jrt@llnl.gov
!                        Al Franz
!                          franz2@llnl.gov
!
! FILE
!   convec_ncar.F
!
! ROUTINES
!   Do_Convec_Ncar
!   Conv_Tran
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Convec_Ncar
!
! DESCRIPTION
!   This is the interface routine to Conv_Tran.  It formats the gem variables
!   to satisfy Conv_Tran.
!
! ARGUMENTS
!   metdata_name_org   : first  part of metdata_name, e.g., "NCAR"
!   metdata_name_model : second part of metdata_name, e.g., "MATCH"
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
!   cldmas       : convective mass flux in     updraft (kg/m^2/s)
!   dtrn         : detrainment rate (DAO:kg/m^2*s, NCAR:s^-1)
!   eu           : entrainment into convective updraft (s^-1)
!   ed           : entrainment into convective downdraft (s^-1)
!   md           : convective mass flux in downdraft (kg/m^2/s)
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
!-----------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------
      subroutine Do_Convec_Ncar  &
     &  (chem_mecha, metdata_name_org, metdata_name_model, det_ent, do_downdraft,  &
     &   do_wetdep,  &
     &   pr_wet_depos, chem_opt, ih2o2_num, ihno3_num, lwi_flags, tdt,  &
     &   mw, mcor, pbl, cldmas, dtrn, eu, ed, md, grid_height, mass, &
     &   wet_depos, kel, press3e, concentration,  &
#ifdef MICRO_AEROSOL
     &   humidity, press3c, REL_SCAV_EFF_new, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ilong, ivert, num_species)
#else
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ilong, ivert, num_species)
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



!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in) :: ilong, ivert, num_species

      character (len=* ) :: chem_mecha
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_model
      logical, intent(in   ) :: det_ent
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
      real*8 , intent(in   ) :: cldmas     (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: dtrn       (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: eu         (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: ed         (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: md         (i1:i2,   ju1:j2,   k1:k2)
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
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      integer :: iku
      integer :: il, ij, ik, ic
      integer :: il2g                      ! gathered index to operate over
#ifdef MICRO_AEROSOL
      integer :: inum
#endif

      integer :: itdt_conv
      integer :: num_conv_steps

      real*8  :: rnum_conv_steps
      real*8  :: tdt_conv

      real*8  :: xmbsth

      integer :: ideep(ilong)              ! gathering array
      integer :: pbli (ilong)              ! index of pbl height

      real*8  :: col_sum_depos   (i1:i2)   ! sum of column wet_deposition loss
      real*8  :: kloss           (i1:i2)   ! wet loss rate constant at all
                                           ! longitudes (s^-1)
      real*8  :: updraft_velocity(i1:i2)   ! velocity in convective updraft
                                           ! (m/s)

      real*8  :: dpi(ilong,  k1:k2)         ! delta pressure between interfaces
      real*8  :: dui(ilong,  k1:k2)         ! mass detraining from updraft
      real*8  :: eui(ilong,  k1:k2)         ! mass entraining into updraft
      real*8  :: mui(ilong,  k1:k2)         ! mass flux up
      real*8  :: mdi(ilong,  k1:k2)         ! mass flux down

      real*8  ::  &
     &  fracis(i1:i2, k1:k2, num_species)  ! insoluble fraction of tracer
      real*8  ::  &
     &  qq    (i1:i2, k1:k2, num_species)  ! tracer array including moisture

! mk
      real*8 :: nemui(ilong, ivert)     ! non-entraining mass flux up
      real*8 :: eneui(ilong, ivert)    ! mass entraining into non-entraining updraft
      real*8 :: dneui(ilong, ivert)    ! mass detraining from non-entraining updraft
      real*8 :: zero(ilong, ivert)     ! array of zeros for non-ent updraft calc

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

!----------------------------------------------------------------------------

!     ----------------
!     Begin execution.
!     ----------------

      ideep(:) = 0
      pbli (:) = 0

      kloss(:) = 0.0d0

      updraft_velocity(:) = 0.0d0

      dpi(:,:) = 0.0d0
      dui(:,:) = 0.0d0
      eui(:,:) = 0.0d0
      mui(:,:) = 0.0d0
      mdi(:,:) = 0.0d0

      fracis(:,:,:) = 0.0d0

      nemui(:,:) = 0.0d0
      eneui(:,:) = 0.0d0
      dneui(:,:) = 0.0d0
      zero(:,:) = 0.0d0


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

      end if

!.sds. make like gocart
       if(IFSO2 .gt. 0) hstar(IFSO2_l) = 600.0d0
       if(INSO2 .gt. 0) hstar(INSO2_l) = 600.0d0
!.sds. 


!     -----------------------------------------------------------
!     Calculate any needed sub-cycling of the convection operator
!     by comparing the mass flux over the full time step to the
!     mass within the grid box.
!     -----------------------------------------------------------

      num_conv_steps =  &
     &  Maxval (tdt * cldmas(:,:,:) * Spread (mcor(:,:), 3, k2-k1+1) /  &
     &          mass(:,:,:)) + 1.0d0

      rnum_conv_steps = num_conv_steps
      tdt_conv        = tdt / rnum_conv_steps

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

        il2g = 0

!       =======
        illoop: do il = i1, i2
!       =======

!         ----------------------------------------------------
!         Verify that there is convection in the current cell.
!         ----------------------------------------------------

          if (Maxval (cldmas(il,ij,:)) >= 0.0d0) then

            il2g = il2g + 1

            ideep(il2g) = il

            dpi(il2g,:) =  &
     &        (press3e(il,ij,k2-1:k1-1:-1) -  &
     &         press3e(il,ij,k2:k1:-1)) *  &
     &        PASPMB


            mui(il2g,:) = cldmas(il,ij,k2:k1:-1) * GMI_G

            if (metdata_name_org(1:3) == 'DAO' .or. &
     &          (metdata_name_org(1:4) == 'GMAO' .and.  &
     &           metdata_name_model(1:5) == 'GEOS5')  ) then

              eui(il2g,k1+1:k2-1) = Max (0.0d0,  &
     &               cldmas(il,ij,k2-1:k1+1:-1) -  &
     &               cldmas(il,ij,k2-2:k1  :-1) +  &
     &               dtrn  (il,ij,k2-1:k1+1:-1))

              eui(il2g,:) = eui(il2g,:) * GMI_G

              if (det_ent) dui(il2g,:) = dtrn(il,ij,k2:k1:-1) * GMI_G

              if (metdata_name_model(1:2) /= 'GS') then

                if(do_downdraft) mdi(il2g,:) = md(il,ij,k2:k1:-1)*GMI_G

              end if

            else

              eui(il2g,:) = eu    (il,ij,k2:k1:-1) * dpi(il2g,:)
!



              if (det_ent)  &
     &          dui(il2g,:) = dtrn(il,ij,k2:k1:-1) * dpi(il2g,:)

              if(do_downdraft) then

!               --------------------------------------------------------------
!               Special algorithm for calculating downdraft mass flux from the
!               GISS downdraft entrainment. Limited to grid boxes 1-7 and
!               limited to be less than or equal to updraft mass flux.
!               --------------------------------------------------------------

                if (metdata_name_org(1:4) == 'GISS') then
!mk
                  mdi(il2g,1:k2-7) = 0.0
!mk jan 21
!                 mdi(il2g,k2-6:k2) = ed(il,ij,1:7) * dpi(il2g,k2-6:k2)
!mk
                  mdi(il2g,k2-6:k2) = ed(il,ij,7:1:-1) * dpi(il2g,k2-6:k2)



                  do ik = k2-6, k2-1

                    mdi(il2g,ik) = mdi(il2g,ik-1) - 2.0d0 * mdi(il2g,ik)

                  end do

                  mdi(il2g,k2) = 0.0d0

                  where (abs(mdi(il2g,:)) > mui(il2g,:))  &
     &              mdi(il2g,:) = -mui(il2g,:)

                else

                  mdi(il2g,:) = md(il,ij,k2:k1:-1) * GMI_G

                end if

              end if

            end if

          end if

!       =============
        end do illoop
!       =============

!mk
!... parameterize deep convection in GISS fields
!
!... need to parameterize "deep" (non-entraining) convection from
!...  the existing GISS fields. 'eui' has been converted to top-down
!...  so need to work from level 5 from the bottom upward
!...  The method was suggested by UCI (Prather)
!
        if(metdata_name_org(1:4) == 'GISS' .and. il2g /= 0) then
!... start with 12% of value at level 5 of met field (empirical)
          nemui(1:il2g,k2-4) = 0.12d0*mui(1:il2g,k2-4)
!... build upward without entrainment
          do ik=k2-5,k1,-1
            nemui(1:il2g,ik) = Min(nemui(1:il2g,ik+1),mui(1:il2g,ik))
          enddo
!... build downward without detrainment
          do ik=k2-3,k2
            nemui(1:il2g,ik) = Min(nemui(1:il2g,ik-1),mui(1:il2g,ik))
          enddo

!... now remove this flux from the original total mass flux up field
!...  => mui now becomes the "shallow" mass flux up
          mui(1:il2g,:) = mui(1:il2g,:) - nemui(1:il2g,:)

!... calculate the entrainment into the "deep" updraft
          eneui(1:il2g,:) =  Max (0.0d0  &
     &      ,(nemui(1:il2g,1:k2-1) - nemui(1:il2g,2:k2)) )

!... redo "shallow" entrainment so total maintained
!...  => eui now becomes the "shallow" entrainment into updraft
          eui(1:il2g,:) =  eui(1:il2g,:) - eneui(1:il2g,:)

!... calculate the detrainment from the deep updraft
!mk
!         if (det_ent) then
!           dneui(1:il2g,:) =  Max (0.0d0
!    &       ,(nemui(1:il2g,2:k2) - nemui(1:il2g,1:k2-1)) )
!... redo "shallow" detrainment so total maintained
!...  => dui now becomes the "shallow" detrainment from updraft
!           dui(1:il2g,:) =  dui(1:il2g,:) - dneui(1:il2g,:)
!         endif

        endif


!       =======
        il2gif: if (il2g /= 0) then
!       =======

          if (do_wetdep) then

!           ----------------------------------------------------------------
!           Calculate the insoluble fraction at each level for each species.
!           ----------------------------------------------------------------

            where ((lwi_flags(i1:i2,ij) == 1) .or.  &
     &             (lwi_flags(i1:i2,ij) == 3))
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

                else if (aerosol(ic) >= 1) then


#ifdef MICRO_AEROSOL
                  kloss(:) = 5.0d-3  &
     &                     * REL_SCAV_EFF_new(:,ij,ik,ic)
#else
                  kloss(:) = 5.0d-3 * REL_SCAV_EFF(aerosol(ic))
#endif

                else if (hstar(ic) > 0.0d0) then

!                 =======================
                  call Calc_Wet_Loss_Rate  &
!                 =======================
     &              (.true., ic, ih2o2_num, delh_298_over_r(ic),  &
     &               hstar(ic), retention_eff(ic), press3e(:,ij,ik),  &
     &               kel(:,ij,ik), kloss(:), i1, i2, ilo, ihi)

                  kloss(:) = 5.0d-3 * kloss(:)

                else

                  kloss(:) = 0.0d0

                end if

                fracis(:,iku,ic) =  &
     &            Exp (-kloss(:) * grid_height(:,ij,ik) /  &
     &                 updraft_velocity(:))

              end do

            end do

          else

            fracis(:,:,:) = 1.0d0

          end if


!         ----------------------------------------------------------------
!         Find the index of the top of the planetary boundary layer (pbl).
!         Convection will assume well mixed tracers below that level.
!         ----------------------------------------------------------------

          do il = 1, il2g

            pbli(il) = 0

            ikloop: do ik = k1, k2

              if (pbl(ideep(il),ij) < Sum (grid_height(ideep(il),ij,k1:ik))) then
                pbli(il) = ik
                exit ikloop
              end if

            end do ikloop

            if (pbli(il) == 0) then

              err_msg = 'Could not find pbl in Do_Convec_Ncar.'
              call GmiPrintError (err_msg, .true., 2, ideep(il), ij,  &
     &                         1, pbl(ideep(il),ij), 0.0d0)

            end if

            pbli(il) = k2 - pbli(il)

          end do


!            write(*,*)'hi mk : num_conv_steps = ',num_conv_steps
!         ==========
          itdtcloop: do itdt_conv = 1, num_conv_steps
!         ==========

            do ic = 1, num_species
            qq(:,k2:k1:-1,ic) = concentration(ic)%pArray3D(:,ij,k1:k2)
            end do

!mk
!... call the "deep" convection first
            if(metdata_name_org(1:4) == 'GISS' ) then
!              ==============
               call Conv_Tran  &
!              ==============
     &          (il2g, tdt_conv, xmbsth, ideep, pbli, dneui, eneui,  &
     &           nemui, zero, dpi, fracis, qq, &
     &           i1, i2, k1, k2, ilong, num_species)
            endif

!... call the "shallow" convection for GISS or total for other met fields
!           ==============
            call Conv_Tran  &
!           ==============
     &        (il2g, tdt_conv, xmbsth, ideep, pbli, dui, eui,  &
     &         mui, mdi, dpi, fracis, qq, &
     &           i1, i2, k1, k2, ilong, num_species)


            if (pr_wet_depos) then

!             ---------------------------------------------------
!             Calculate the wet deposition from the difference in
!             mixing ratio.
!             ---------------------------------------------------

              do ic = 1, num_species

                if ((hstar(ic) > 0.0d0) .or.  &
     &              (aerosol(ic) >= 1)  .or.  &
     &              (ic == ihno3_num)) then

                  col_sum_depos(i1:i2) = 0.0d0

                  do ik = k1, k2

                    iku = k2 - ik + 1

                    col_sum_depos(:) =  &
     &                col_sum_depos(:) +  &
     &                (concentration(ic)%pArray3D(:,ij,ik) - qq(:,iku,ic)) *  &
     &                mass(:,ij,ik)

                  end do

                  where (col_sum_depos(i1:i2) > 0.0d0)  &
     &              wet_depos(i1:i2,ij,ic) =  &
     &                wet_depos(i1:i2,ij,ic) +  &
     &                col_sum_depos(i1:i2) * mw(ic) /  &
     &                MWTAIR / mcor(i1:i2,ij)

                end if

              end do

            end if

            do ic = 1, num_species
              concentration(ic)%pArray3D(:,ij,k1:k2) = qq(:,k2:k1:-1,ic)
            end do

!         ================
          end do itdtcloop
!         ================

!       =============
        end if il2gif
!       =============

!     =============
      end do ijloop
!     =============

      return

      end subroutine Do_Convec_Ncar


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Conv_Tran
!
! DESCRIPTION
!   This routine performs convective transport of trace species.
!
!   Note that we assume that the tracers are in a moist mixing ratio
!   (this should change soon).
!
! ARGUMENTS
!   INPUT:
!     il2g   : gathered max lon indices over which to operate
!     delt   : convection time step (s)
!     xmbsth : threshold below which we treat mass fluxes as zero (mb/s)
!     ideep  : gathering array
!     pbli   : index of planetary boundary layer
!     dui    : mass detraining from updraft
!     eui    : mass entraining into updraft
!     mui    : mass flux up
!     mdi    : mass flux down
!     dpi    : delta pressure between interfaces
!     fracis : insoluble fraction of tracer
!   INPUT/OUTPUT:
!     qq     : tracer array including moisture (mixing ratio)
!
!-----------------------------------------------------------------------------

      subroutine Conv_Tran  &
     &  (il2g, delt, xmbsth, ideep, pbli, dui, eui, mui, mdi, dpi,  &
     &   fracis, qq, i1, i2, k1, k2, ilong, num_species)

      use GmiSpcConcentrationMethod_mod, only : isFixedConcentration
      use GmiPrintError_mod, only : GmiPrintError

      implicit  none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, k1, k2, ilong, num_species
      integer :: il2g
      real*8  :: delt
      real*8  :: xmbsth
      integer :: ideep (ilong)
      integer :: pbli  (ilong)
      real*8, intent(in )  :: dui   (ilong, k1:k2)
      real*8, intent(in )  :: eui   (ilong, k1:k2)
      real*8, intent(in )  :: mui   (ilong, k1:k2)
      real*8, intent(in )  :: mdi   (ilong, k1:k2)
      real*8  :: dpi   (ilong, k1:k2)
      real*8  :: fracis(i1:i2, k1:k2, num_species)
      real*8  :: qq    (i1:i2, k1:k2, num_species)


!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter :: SMALL = 1.0d-36

!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      integer :: il, ik, ic
      integer :: km1, kp1

      real*8  :: avg_pbl    ! average mixing ratio in pbl
      real*8  :: cabv       ! mix ratio of constituent above
      real*8  :: cbel       ! mix ratio of constituent below
      real*8  :: cdifr      ! normalized diff. between cabv and cbel
      real*8  :: fluxin     ! flux coming into  each box at a k level
      real*8  :: fluxout    ! flux going out of each box at a k level
      real*8  :: maxc
      real*8  :: minc
      real*8  :: scav       ! mixing ratio of a scavenged tracer
      real*8  :: sqrt_fisg
      real*8  :: temp_conu  !.sds - new temp variable for updraft concen calc

      real*8  :: chat  (il2g, k1:k2)  ! mix ratio in env.      at interfaces
      real*8  :: conu  (il2g, k1:k2)  ! mix ratio in updraft   at interfaces
      real*8  :: cond  (il2g, k1:k2)  ! mix ratio in downdraft at interfaces
      real*8  :: dcondt(il2g, k1:k2)  ! gathered tend.  array
      real*8  :: edi   (il2g, k1:k2)  ! gathered downdraft entrainment
      real*8  :: fisg  (il2g, k1:k2)  ! gathered insoluble frac. of tracer
      real*8  :: xconst(il2g, k1:k2)  ! gathered tracer array


!     ----------------
!     Begin execution.
!     ----------------
!     =======
      icloop: do ic = 1, num_species
!     =======

!       ====================================
        if (isFixedConcentration(ic)) cycle icloop
!       ====================================

!       ------------------------------------------------
!       Gather up the constituent and set tend. to zero.
!       ------------------------------------------------

        do ik = k1, k2
          do il = 1, il2g

            fisg  (il,ik) = fracis(ideep(il),ik,ic)
            xconst(il,ik) = qq    (ideep(il),ik,ic)

          end do
        end do

!       -----------------------------------------
!       From now on work only with gathered data.
!       -----------------------------------------

!       ----------------------------------------------------
!       Interpolate environment tracer values to interfaces.
!       ----------------------------------------------------

!       ========
        ikloop1: do ik = k1, k2
!       ========

          km1 = Max (k1, ik-1)
          kp1 = Min (k2, ik+1)

          edi(:,ik) = mdi(:,ik) - mdi(:,kp1)

!         ========
          illoop1: do il = 1, il2g
!         ========

            minc = Min (xconst(il,km1), xconst(il,ik))
            maxc = Max (xconst(il,km1), xconst(il,ik))

            if (minc < 0.0d0) then

              cdifr = 0.0d0

            else

              cdifr = Abs (xconst(il,ik) - xconst(il,km1)) /  &
     &                Max (maxc, SMALL)

            end if


            if (cdifr > 1.0d-6) then

!             -------------------------------------------------------
!             The two layers differ significantly, so use a geometric
!             averaging procedure.
!             -------------------------------------------------------

              cabv = Max (xconst(il,km1), (maxc * 1.0d-12))
              cbel = Max (xconst(il,ik),  (maxc * 1.0d-12))

              chat(il,ik) =  &
     &          Log (cabv/cbel) / (cabv - cbel) * cabv * cbel

            else

!             ---------------------------------------------------
!             The two layers have only a small difference, so use
!             arithmetic mean.
!             ---------------------------------------------------

              chat(il,ik) = 0.5d0 * (xconst(il,ik) + xconst(il,km1))

            end if


!           -------------------------------------
!           Provisional up and down draft values.
!           -------------------------------------

            conu(il,ik) = chat(il,ik)
            cond(il,ik) = chat(il,km1)

!           ------------------
!           Provisional tends.
!           ------------------

            dcondt(il,ik) = 0.0d0

!         ==============
          end do illoop1
!         ==============

!       ==============
        end do ikloop1
!       ==============

        where (edi(:,:) < 0.0d0) edi(:,:) = 0.0d0

!       ========
        illoop2: do il = 1, il2g
!       ========

!         ----------------------------------------------------------------------
!         Find the mixing ratio in the downdraft, top of atmosphere down.
!         NOTE: mass flux   in downdraft (mdi) will be zero or negative.
!               entrainment in downdraft (edi) will be zero or positive.
!         ----------------------------------------------------------------------

!         =========
          ikloopdd: do ik = 2, k2
!         =========

            if     ((mdi(il,ik-1) < -xmbsth)  &
     &        .or.  (edi(il,ik)   >  xmbsth)) then

              cond(il,ik) = (cond(il,ik-1)*mdi(il,ik-1)  &
     &                     - xconst(il,ik)*edi(il,ik))  &
     &                    / (mdi(il,ik-1) - edi(il,ik))

            end if

!         ===============
          end do ikloopdd
!         ===============

!         ------------------------------------------------------
!         Calculate updrafts with scavenging from bottom to top.
!         Include the downdrafts.
!         ------------------------------------------------------

!         -------------------------------------
!         Do the bottom most levels in the pbl.
!         -------------------------------------

          if (Sum (dpi(il,k2:pbli(il):-1)) == 0.0d0)  then

            err_msg = 'Problem in Conv_Tran.'
            call GmiPrintError (err_msg, .true., 2, ideep(il), pbli(il),  &
     &                       il, dpi(il,k2), dpi(il,pbli(il)))

          endif

          avg_pbl =  &
     &      Sum (xconst(il,k2:pbli(il):-1) * dpi(il,k2:pbli(il):-1)) /  &
     &      Sum (dpi(il,k2:pbli(il):-1))

          sqrt_fisg = Sqrt (fisg(il,pbli(il)))

          scav    = avg_pbl * (1.0d0 -  sqrt_fisg)

          conu(il,pbli(il)) = avg_pbl * sqrt_fisg

          fluxin  = (mui(il,pbli(il)) + mdi(il,pbli(il))) *  &
     &          Min (chat(il,pbli(il)), xconst(il,pbli(il)-1))  &
     &            - (mdi(il,pbli(il)) * cond(il,pbli(il)))

          fluxout = mui(il,pbli(il)) * conu(il,pbli(il)) +  &
     &              mui(il,pbli(il)) * scav

          dcondt(il,k2) = (fluxin - fluxout) /  &
     &                    Sum (dpi(il,k2:pbli(il):-1))

          if (avg_pbl > 0.0d0) then

            do ik = k2, pbli(il), -1

              qq(ideep(il),ik,ic) =  &
     &          (avg_pbl + (dcondt(il,k2)*delt)) /  &
     &          avg_pbl * xconst(il,ik)

            end do

          end if


!         ---------------------------
!         Loop over all other levels.
!         ---------------------------

!         ========
          ikloop2: do ik = pbli(il)-1, 2, -1
!         ========

            kp1 = ik + 1
            km1 = ik - 1

!           -----------------------------
!           Find mixing ratio in updraft.
!           -----------------------------

!!!!!! Old convection
!
!            if ( (mui(il,kp1) - dui(il,ik) + eui(il,ik)) > xmbsth) then
!
!              sqrt_fisg   = Sqrt (fisg(il,ik))
!
!              scav        = ((mui(il,kp1)-dui(il,ik)) * conu(il,kp1)  *  &
!     &                       (1.0d0 - fisg(il,ik)) +  &
!     &                       eui(il,ik)  * xconst(il,ik) *  &
!     &                       (1.0d0 - sqrt_fisg)) /  &
!     &                      ((mui(il,kp1)-dui(il,ik)) + eui(il,ik))
!
!              conu(il,ik) = ((mui(il,kp1)-dui(il,ik)) * conu(il,kp1)  *  &
!     &                       fisg(il,ik) +  &
!     &                       eui(il,ik)  * xconst(il,ik) *  &
!     &                       sqrt_fisg) /  &
!     &                      ((mui(il,kp1)-dui(il,ik)) + eui(il,ik))
!
!            else
!
!              conu(il,ik) = xconst(il,ik)
!
!              scav        = 0.0d0
!
!            end if
!
!!!!!!! New convection
!
            if ( (mui(il,kp1) - dui(il,ik) + eui(il,ik)) > xmbsth) then

              sqrt_fisg   = Sqrt (fisg(il,ik))

!... mix updraft and entrainment, then detrain
!!Old Code with a bug
!              temp_conu  = ( mui(il,kp1) * conu(il,kp1) +  &
!     &                       eui(il,ik)  * xconst(il,ik) ) /  &
!     &                      (mui(il,kp1) + eui(il,ik))

!--Updated code:
!... decide which mixing in the updraft because of potential divide by 0.0

              if ( abs(mui(il,kp1) + eui(il,ik)) .ge. xmbsth) then
                !... mix updraft and entrainment, then detrain this concentration
                temp_conu  = ( mui(il,kp1) * conu(il,kp1) +  &
     &                         eui(il,ik)  * xconst(il,ik) ) /  &
     &                        (mui(il,kp1) + eui(il,ik))
                !... if updraft and entrainment flux sum to 0.0, then take average concen
              else
                temp_conu  = (xconst(il,ik) + conu(il,kp1)) / 2.0d0
              endif
!
              scav        = ((mui(il,kp1) * conu(il,kp1) -  &
     &                        dui(il,ik) * temp_conu ) *  &
     &                       (1.0d0 - fisg(il,ik)) +  &
     &                       eui(il,ik)  * xconst(il,ik) *  &
     &                       (1.0d0 - sqrt_fisg)) /  &
     &                      (mui(il,kp1) - dui(il,ik) + eui(il,ik))

              conu(il,ik) = ((mui(il,kp1) * conu(il,kp1) -  &
     &                        dui(il,ik) * temp_conu ) *  &
     &                       fisg(il,ik) +  &
     &                       eui(il,ik)  * xconst(il,ik) *  &
     &                       sqrt_fisg) /  &
     &                      (mui(il,kp1) - dui(il,ik) + eui(il,ik))
            else

              conu(il,ik) = xconst(il,ik)
              scav        = 0.0d0
            end if
!
!           -------------------------------------------------------
!           Calculate fluxes into and out of box.  With scavenging
!           included the net flux for the whole column is no longer
!           guaranteed to be zero.
!           Include the downdrafts.
!           -------------------------------------------------------

            fluxin  = mui(il,kp1) * conu(il,kp1) +  &
     &                (mui(il,ik)+mdi(il,ik))  &
     &                * Min (chat(il,ik), xconst(il,km1))  &
     &                - (mdi(il,ik) * cond(il,ik))

            fluxout = mui(il,ik)  * conu(il,ik) +  &
     &                (mui(il,kp1) + mdi(il,kp1))  &
     &                * Min (chat(il,kp1), xconst(il,ik)) +  &
     &                (mui(il,ik) + mui(il,kp1)) * 0.5d0 * scav  &
     &                - (mdi(il,kp1) * cond(il,kp1))

            dcondt(il,ik) = (fluxin - fluxout) / dpi(il,ik)

!           --------------------------------------------
!           Update and scatter data back to full arrays.
!           --------------------------------------------
            qq(ideep(il),ik,ic) = xconst(il,ik) + (dcondt(il,ik) * delt)

!         ==============
          end do ikloop2
!         ==============

!       ==============
        end do illoop2
!       ==============

!     =============
      end do icloop
!     =============


      return

      end



