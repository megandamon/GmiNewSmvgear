!=============================================================================
!
! $Id: convec_ncar_old.F90,v 1.15 2012-01-05 21:41:30 jkouatch Exp $
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
!   old_convec_ncar.F
!
! ROUTINES
!   Do_Convec_Ncar_old
!   Conv_Tran_old
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Convec_Ncar_old
!
! DESCRIPTION
!   This is the interface routine to Conv_Tran_old.  It formats the gem variables
!   to satisfy Conv_Tran_old.
!
! ARGUMENTS
!   metdata_name_org   : first  part of metdata_name, e.g., "NCAR"
!   metdata_name_model : second part of metdata_name, e.g., "MATCH"
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
!   cldmas       : convective mass flux in     updraft (kg/m^2*s)
!   dtrn         : detrainment rate (DAO:kg/m^2*s, NCAR:s^-1)
!   eu           : entrainment into convective updraft (s^-1)
!   grid_height  : grid box height  (m)
!   mass         : mass of air in each grid box (kg)
!   wet_depos    : accumulated wet deposition   (kg/m^2)
!   kel          : temperature      (degK)
!   press3e      : atmospheric pressure at the edge of each grid box (mb)
!   const        : species concentration, known at zone centers (mixing ratio)
!
!-----------------------------------------------------------------------------

      subroutine Do_Convec_Ncar_old  &
     &  (metdata_name_org, metdata_name_model, do_wetdep,  &
     &   pr_wet_depos, chem_opt, ih2o2_num, ihno3_num, lwi_flags, tdt,  &
     &   mw, mcor, pbl, cldmas, dtrn, eu, grid_height, mass, wet_depos,  &
     &   kel, press3e, concentration, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ilong, ivert, num_species)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiPrintError_mod, only : GmiPrintError

      implicit  none

#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"
#ifdef GOCARTaerosol
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
      integer, intent(in) :: ilong, ivert, num_species

      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_model
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
      real*8 , intent(in   ) :: grid_height(i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in   ) :: mass       (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(inout) :: wet_depos  (i1:i2,   ju1:j2,   num_species)
      real*8 , intent(in   ) :: kel        (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in   ) :: press3e    (ilo:ihi, julo:jhi, k1-1:k2)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)

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

      real*8  :: dpi(ilong, ivert)         ! delta pressure between interfaces
      real*8  :: eui(ilong, ivert)         ! mass entraining from updraft
      real*8  :: mui(ilong, ivert)         ! mass flux up

      real*8  ::  &
     &  fracis(i1:i2, k1:k2, num_species)  ! insoluble fraction of tracer
      real*8  ::  &
     &  qq    (i1:i2, k1:k2, num_species)  ! tracer array including moisture


!     ----------------
!     Begin execution.
!     ----------------

      ideep(:) = 0
      pbli (:) = 0

      kloss(:) = 0.0d0

      updraft_velocity(:) = 0.0d0

      dpi(:,:) = 0.0d0
      eui(:,:) = 0.0d0
      mui(:,:) = 0.0d0

      fracis(:,:,:) = 0.0d0


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
     &  Maxval (tdt * cldmas(:,:,:) * Spread (mcor(:,:), 3, ivert) /  &
     &          mass(:,:,:)) +  &
     &         1.0d0

      rnum_conv_steps = num_conv_steps
      tdt_conv        = tdt / rnum_conv_steps


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


            if ((metdata_name_org  (1:4) == 'NCAR') .and.  &
     &          (metdata_name_model(1:4) == 'CCM3')) then

              mui(il2g,:) = cldmas(il,ij,k2:k1:-1)
              eui(il2g,:) = eu    (il,ij,k2:k1:-1)

              do ik = k1+1, k2
                eui(il2g,ik) =  &
     &            Max (eui(il2g,ik), (mui(il2g,ik-1) - mui(il2g,ik)))
              end do

            else if (metdata_name_org(1:3) == 'DAO') then

              mui(il2g,:) = cldmas(il,ij,k2:k1:-1) * GMI_G

              eui(il2g,k1+1:k2-1) =  &
     &          Max (0.0d0,  &
     &               cldmas(il,ij,k2-1:k1+1:-1) -  &
     &               cldmas(il,ij,k2-2:k1  :-1) +  &
     &               dtrn  (il,ij,k2-1:k1+1:-1))

              eui(il2g,:) = eui(il2g,:) * GMI_G

            else

              mui(il2g,:) = cldmas(il,ij,k2:k1:-1) * GMI_G
              eui(il2g,:) = eu    (il,ij,k2:k1:-1) * dpi(il2g,:)

            end if

          end if

!       =============
        end do illoop
!       =============


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

                if (ic == ihno3_num) then

                  kloss(:) = 5.0d-3

                else if (aerosol(ic) >= 1) then

                  kloss(:) = 5.0d-3 * REL_SCAV_EFF(aerosol(ic))

                else if (hstar(ic) > 0.0d0) then

!                 =======================
                  call Calc_Wet_Loss_Rate  &
!                 =======================
     &              (.true., ic, ih2o2_num, delh_298_over_r(ic), hstar(ic),  &
     &               retention_eff(ic), press3e(:,ij,ik), kel(:,ij,ik),  &
     &               kloss(:), i1, i2, ilo, ihi)

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

              if (pbl(ideep(il),ij) <  &
     &            Sum (grid_height(ideep(il),ij,k1:ik))) then
                pbli(il) = ik
                exit ikloop
              end if

            end do ikloop

            if (pbli(il) == 0) then

              err_msg = 'Could not find pbl in Do_Convec_Ncar_old.'
              call GmiPrintError (err_msg, .true., 2, ideep(il), ij,  &
     &                         1, pbl(ideep(il),ij), 0.0d0)

            end if

            pbli(il) = k2 - pbli(il)

          end do


!         ==========
          itdtcloop: do itdt_conv = 1, num_conv_steps
!         ==========

            do ic = 1, num_species
               qq(:,k2:k1:-1,ic) = concentration(ic)%pArray3D(:,ij,k1:k2)
            end do

!           ==============
            call Conv_Tran_old  &
!           ==============
     &        (il2g, tdt_conv, xmbsth, ideep, pbli, eui,  &
     &         mui, dpi, fracis, qq, &
     &         i1, i2, k1, k2, ilong, ivert, num_species)


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

      end subroutine Do_Convec_Ncar_old


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Conv_Tran_old
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
!     eui    : mass entraining from updraft (s^-1)
!     mui    : mass flux up
!     dpi    : delta pressure between interfaces
!     fracis : insoluble fraction of tracer
!   INPUT/OUTPUT:
!     qq     : tracer array including moisture (mixing ratio)
!
!-----------------------------------------------------------------------------

      subroutine Conv_Tran_old  &
     &  (il2g, delt, xmbsth, ideep, pbli, eui, mui, dpi, fracis, qq, &
     &   i1, i2, k1, k2, ilong, ivert, num_species)

      use GmiSpcConcentrationMethod_mod, only : isFixedConcentration
      use GmiPrintError_mod, only : GmiPrintError

      implicit  none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, k1, k2, ilong, ivert, num_species
      integer :: il2g
      real*8  :: delt
      real*8  :: xmbsth
      integer :: ideep (ilong)
      integer :: pbli  (ilong)
      real*8  :: eui   (ilong, ivert)
      real*8  :: mui   (ilong, ivert)
      real*8  :: dpi   (ilong, ivert)
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

      real*8  :: chat  (il2g,ivert)  ! mix ratio in env.      at interfaces
      real*8  :: conu  (il2g,ivert)  ! mix ratio in updraft   at interfaces
      real*8  :: dcondt(il2g,ivert)  ! gathered tend.  array
      real*8  :: fisg  (il2g,ivert)  ! gathered insoluble frac. of tracer
      real*8  :: xconst(il2g,ivert)  ! gathered tracer array


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


!       ========
        illoop2: do il = 1, il2g
!       ========

!         ------------------------------------------------------
!         Calculate updrafts with scavenging from bottom to top.
!         ------------------------------------------------------

!         -------------------------------------
!         Do the bottom most levels in the pbl.
!         -------------------------------------

          if (Sum (dpi(il,k2:pbli(il):-1)) == 0.0d0)  then

            err_msg = 'Problem in Conv_Tran_old.'
            call GmiPrintError (err_msg, .true., 2, ideep(il), pbli(il),  &
     &                       0, 0.0d0, 0.0d0)

          end if

          avg_pbl =  &
     &      Sum (xconst(il,k2:pbli(il):-1) * dpi(il, k2:pbli(il):-1)) /  &
     &      Sum (dpi(il,k2:pbli(il):-1))

          sqrt_fisg = Sqrt (fisg(il,pbli(il)))

          scav    = avg_pbl * (1.0d0 -  sqrt_fisg)

          conu(il,pbli(il)) = avg_pbl * sqrt_fisg

          fluxin  = mui(il,pbli(il)) *  &
     &              Min (chat(il,pbli(il)), xconst(il,pbli(il)-1))

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

            if ((mui(il,kp1) + eui(il,ik)) > xmbsth) then

              sqrt_fisg   = Sqrt (fisg(il,ik))

              scav        = (mui(il,kp1) * conu(il,kp1)  *  &
     &                       (1.0d0 - fisg(il,ik)) +  &
     &                       eui(il,ik)  * xconst(il,ik) *  &
     &                       (1.0d0 - sqrt_fisg)) /  &
     &                      (mui(il,kp1) + eui(il,ik))

              conu(il,ik) = (mui(il,kp1) * conu(il,kp1)  *  &
     &                       fisg(il,ik) +  &
     &                       eui(il,ik)  * xconst(il,ik) *  &
     &                       sqrt_fisg) /  &
     &                      (mui(il,kp1) + eui(il,ik))

            else

              conu(il,ik) = xconst(il,ik)

              scav        = 0.0d0

            end if

!           -------------------------------------------------------
!           Calculate fluxes into and out of box.  With scavenging
!           included the net flux for the whole column is no longer
!           guaranteed to be zero.
!           -------------------------------------------------------

            fluxin  = mui(il,kp1) * conu(il,kp1) +  &
     &                mui(il,ik)  * Min (chat(il,ik), xconst(il,km1))

            fluxout = mui(il,ik)  * conu(il,ik) +  &
     &                mui(il,kp1) * Min (chat(il,kp1), xconst(il,ik)) +  &
     &                (mui(il,ik) + mui(il,kp1)) * 0.5d0 * scav

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



