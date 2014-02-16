
!=============================================================================
!
! $Id: convec_update.F90,v 1.7 2009-11-13 18:36:19 kouatch Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   convec_update.F
!
! ROUTINES
!   Update_Convec
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update_Convec
!
! DESCRIPTION
!   This routine updates convection.
!
! ARGUMENTS
!   metdata_name_org   : first  part of metdata_name, e.g., "NCAR"
!   metdata_name_model : second part of metdata_name, e.g., "MATCH"
!   det_ent      : flag for doing detrainment then entrainment
!   do_downdraft : flag for doing downdrafts
!   do_wetdep    : flag for doing wet deposition
!   pr_wet_depos : flag for saving and then printing accumulated wet
!                  deposition
!   chem_opt     : chemistry  option
!   convec_opt   : convection option
!   ih2o2_num    : concentration array index for H2O2
!   ihno3_num    : concentration array index for HNO3
!   lwi_flags    : land water ice flag (1-water 2-land 3-ice)
!   tdt          : model time step (s)
!   pt           : pressure = (ai * pt) + (bi * psf)   (mb)
!   ai           : pressure = (ai * pt) + (bi * psf), ai at zone interface
!   bi           : pressure = (ai * pt) + (bi * psf), bi at zone interface
!   mw           : array of species' molecular weights (g/mol)
!   mcor         : horizontal area of each grid box    (m^2)
!   pbl          : planetary boundary layer thickness  (m)
!   cmf          : convective mass flux   (kg/m^2*s)
!   dtrn         : detrainment rate (DAO:kg/m^2*s, NCAR:s^-1)
!   eu           : entrainment into convective updraft (s^-1)
!   ed           : entrainment into convective downdraft (s^-1)
!   md           : convective mass flux in downdraft (kg/m^2/s)
!   grid_height  : grid box height  (m)
!   mass         : mass of air in each grid box   (kg)
!   wet_depos    : accumulated wet deposition     (kg/m^2)
!   pctm2        : CTM surface pressure at t1+tdt (mb)
!   kel          : temperature (degK)
!   press3e      : atmospheric pressure at the edge of each grid box (mb)
!   concentration: species concentration, known at zone centers (mixing ratio)
!-micro_aerosol
!   humidity     : specific humidity (g/kg)
!   press3c      : atmospheric pressure at the center of each grid box (mb)
!   REL_SCAV_EFF_new
!                : relative scavenging efficiency for aerosols (0-1)
!
!-----------------------------------------------------------------------------

      subroutine Update_Convec  &
     &  (chem_mecha, metdata_name_org, metdata_name_model, det_ent, do_downdraft,  &
     &   do_old_ncar, do_wetdep,  &
     &   pr_wet_depos, chem_opt, convec_opt, ih2o2_num, ihno3_num,  &
     &   lwi_flags, tdt, mw, mcor, pbl, cmf, dtrn, eu,  &
     &   ed, md, zmdu, zmeu, zmed, zmmd, zmmu, hkdu, hkeu, hkmu,    &
     &   grid_height, mass, wet_depos, kel, press3e, concentration, bmass,  &
#ifdef MICRO_AEROSOL
     &   humidity, press3c, REL_SCAV_EFF_new, &
     &   pr_diag, loc_proc, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ilong, ivert, num_species)
#else
     &   pr_diag, loc_proc, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ilong, ivert, num_species)
#endif

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      implicit none

#     include "GmiParameters.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical           , intent(in   ) :: pr_diag
      integer           , intent(in   ) :: loc_proc
      integer           , intent(in   ) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer           , intent(in   ) :: ilong, ivert, num_species

      character (len=* ), intent(in   ) :: chem_mecha
      character (len=MAX_LENGTH_MET_NAME), intent(in   ) :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME), intent(in   ) :: metdata_name_model
      logical           , intent(in   ) :: det_ent
      logical           , intent(in   ) :: do_downdraft
      logical           , intent(in   ) :: do_old_ncar
      logical           , intent(in   ) :: do_wetdep
      logical           , intent(in   ) :: pr_wet_depos
      integer           , intent(in   ) :: chem_opt
      integer           , intent(in   ) :: convec_opt
      integer           , intent(in   ) :: ih2o2_num
      integer           , intent(in   ) :: ihno3_num
      integer           , intent(in   ) :: lwi_flags(i1:i2, ju1:j2)
      real*8            , intent(in   ) :: tdt
      real*8            , intent(in   ) :: mw(num_species)
      real*8            , intent(in   ) :: mcor       (i1:i2, ju1:j2)
      real*8            , intent(in   ) :: pbl        (i1:i2, ju1:j2)
      real*8            , intent(in   ) :: cmf        (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: dtrn       (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: eu         (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: ed         (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: md         (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: grid_height(i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: mass       (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: bmass      (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(inout) :: wet_depos  (i1:i2, ju1:j2, num_species)
      real*8            , intent(in   ) :: kel        (ilo:ihi, julo:jhi, k1:k2)
      real*8            , intent(in   ) :: press3e    (ilo:ihi, julo:jhi, k1-1:k2)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)
!... GEOS4 convective fluxes
      real*8            , intent(in   ) :: zmdu       (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: zmeu       (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: zmed       (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: zmmd       (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: zmmu       (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: hkdu       (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: hkeu       (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: hkmu       (i1:i2, ju1:j2, k1:k2)
#ifdef MICRO_AEROSOL
      real*8            , intent(in   ) :: humidity   (i1:i2, ju1:j2, k1:k2)
      real*8            , intent(in   ) :: press3c    (ilo:ihi, julo:jhi, k1:k2)
      real*8            , intent(out  ) :: REL_SCAV_EFF_new(i1:i2, ju1:j2, k1:k2, num_species)
#endif


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Update_Convec called by ', loc_proc
      end if


      if (convec_opt == 1) then

!       ===================
        call Do_Convec_Dao2  &
!       ===================
     &    (tdt, cmf, dtrn, concentration, bmass, &
           i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species)

      else if (convec_opt == 2) then

! mk

        if(do_old_ncar) then

         write(*,*) 'calling old convec ncar'
!       =======================
        call Do_Convec_Ncar_old  &
!       =======================
     &  (metdata_name_org, metdata_name_model, do_wetdep,  &
     &   pr_wet_depos, chem_opt, ih2o2_num, ihno3_num, lwi_flags, tdt,  &
     &   mw, mcor, pbl, cmf, dtrn, eu, grid_height, mass, wet_depos,  &
     &   kel, press3e, concentration, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ilong, ivert, num_species)

        else

!       ===================
        call Do_Convec_Ncar  &
!       ===================
     &  (chem_mecha, metdata_name_org, metdata_name_model, det_ent, do_downdraft,  &
     &   do_wetdep,  &
     &   pr_wet_depos, chem_opt, ih2o2_num, ihno3_num, lwi_flags, tdt,  &
     &   mw, mcor, pbl, cmf, dtrn, eu, ed, md, grid_height, mass, &
     &   wet_depos, kel, press3e, concentration,  &
#ifdef MICRO_AEROSOL
     &   humidity, press3c, REL_SCAV_EFF_new, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ilong, ivert, num_species)
#else
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ilong, ivert, num_species)
#endif

        end if

      else if (convec_opt == 3) then

!       ===================
        call Do_Convec_GEOS4  &
!       ===================
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

      end if


      return

      end

