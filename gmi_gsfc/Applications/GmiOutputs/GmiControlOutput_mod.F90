!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlOutput_mod
!
      module GmiControlOutput_mod
!
! !USES:
      use ftiming_dao
      use GmiTimeControl_mod, only : t_gmiClock
      use GmiDiagnosticsMethod_mod, only : t_Diagnostics
      use GmiDiagnosticsMethod_mod, only : Get_pr_qj, Get_pr_qk, Get_pr_qqjk,  &
     &       Get_pr_sad, Get_pr_tend, Get_pr_ascii, Get_pr_cloud, Get_pr_flux, &
     &       Get_pr_const, Get_pr_freq1, Get_pr_freq2, Get_pr_freq3,           &
     &       Get_pr_freq4, Get_pr_restart, Get_pr_netcdf, Get_pr_col_diag,     &
     &       Get_pr_AerDust, Get_pr_overpass1, Get_pr_overpass2,               &
     &       Get_pr_overpass3, Get_pr_overpass4, Get_do_ftiming, Get_pr_diag,  &
     &       Get_pr_time
      use GmiEmissionMethod_mod  , only : t_Emission
      use GmiChemistryMethod_mod , only : t_Chemistry, Get_phot_opt
      use GmiDepositionMethod_mod, only : t_Deposition
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_procID
      use GmiGrid_mod, only : t_gmiGrid
      use GmiMetFieldsControl_mod, only : t_metFields
      use GmiAdvectionMethod_mod , only : t_Advection
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: controlOutputFiles, initializeOutputFiles
!
#     include "gmi_time_constants.h"
#     include "GmiParameters.h"
!
! !AUTHOR:
! Jules Kouatchou, Jules.Kouatchou-1@nasa.gov
!
!EOP
!=============================================================================
      contains
!=============================================================================
!BOP
!
! !IROUTINE: initializeOutputFiles
!
! !INTERFACE:
!
      subroutine initializeOutputFiles(gmiGrid, gmiDomain, Diagnostics,        &
     &               SpeciesConcentration, Chemistry, Emission, metFields,     &
     &               chem_mecha)
!
! !USES:
      use GmiControlColumn_mod       , only : initializeOutputColumn
      use GmiControlRestart_mod      , only : initializeOutputRestart
      use GmiControlASCII_mod        , only : initializeOutputASCII
      use GmiControlFlux_mod         , only : initializeOutputFlux
      use GmiControlQj_mod           , only : initializeOutputQj
      use GmiControlQk_mod           , only : initializeOutputQk
      use GmiControlSAD_mod          , only : initializeOutputSAD
      use GmiControlQqjk_mod         , only : initializeOutputQqjk
      use GmiControlFreq1_mod        , only : initializeOutputFreq1
      use GmiControlFreq2_mod        , only : initializeOutputFreq2
      use GmiControlFreq3_mod        , only : initializeOutputFreq3
      use GmiControlFreq4_mod        , only : initializeOutputFreq4
      use GmiControlOverpass1_mod    , only : initializeOutputOverpass1
      use GmiControlOverpass2_mod    , only : initializeOutputOverpass2
      use GmiControlOverpass3_mod    , only : initializeOutputOverpass3
      use GmiControlOverpass4_mod    , only : initializeOutputOverpass4
      use GmiControlTendencies_mod   , only : initializeOutputTendencies
      use GmiControlAerosolDust_mod  , only : initializeOutputAerosolDust
      use GmiControlConcentration_mod, only : initializeOutputConcentration
      use GmiControlCloudModuleGT_mod, only : initializeOutputCloudGT
!
      implicit none
!
! !INPUT PARAMETERS:
      character(len=*)            , intent(in) :: chem_mecha
      type(t_gmiGrid             ), intent(in) :: gmiGrid
      type(t_gmiDomain           ), intent(in) :: gmiDomain
      type(t_Diagnostics         ), intent(in) :: Diagnostics
      type(t_metFields           ), intent(in) :: metFields
      type(t_Chemistry           ), intent(in) :: Chemistry
      type(t_Emission            ), intent(in) :: Emission
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !DESCRIPTION:
! Initializes all the GMI main components.
!
! !LOCAL VARIABLES:
      integer :: procID, phot_opt
      logical :: pr_diag
      logical :: pr_const, pr_sad, pr_qj, pr_qk, pr_qqjk, pr_cloud,            &
     &           pr_AerDust, pr_flux, pr_time, pr_overpass1, pr_overpass2,     &
     &           pr_overpass3, pr_overpass4, pr_restart, pr_freq1, pr_freq2,   &
     &           pr_freq3, pr_freq4, pr_ascii, pr_netcdf, pr_col_diag, pr_tend

!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID (gmiDomain, procID)
      call Get_pr_diag (Diagnostics, pr_diag)

      if (pr_diag) Write(6,*) ' initializeOutputFiles called by ', procID

      call Get_pr_qj       (Diagnostics, pr_qj   )
      call Get_pr_qk       (Diagnostics, pr_qk   )
      call Get_pr_sad      (Diagnostics, pr_sad  )
      call Get_pr_qqjk     (Diagnostics, pr_qqjk )
      call Get_pr_tend     (Diagnostics, pr_tend )
      call Get_pr_flux     (Diagnostics, pr_flux )
      call Get_pr_ascii    (Diagnostics, pr_ascii)
      call Get_pr_cloud    (Diagnostics, pr_cloud)
      call Get_pr_const    (Diagnostics, pr_const)
      call Get_pr_freq1    (Diagnostics, pr_freq1  )
      call Get_pr_freq2    (Diagnostics, pr_freq2  )
      call Get_pr_freq3    (Diagnostics, pr_freq3  )
      call Get_pr_freq4    (Diagnostics, pr_freq4  )
      call Get_pr_netcdf   (Diagnostics, pr_netcdf)
      call Get_pr_restart  (Diagnostics, pr_restart)
      call Get_pr_AerDust  (Diagnostics, pr_AerDust)
      call Get_pr_col_diag (Diagnostics, pr_col_diag)
      call Get_pr_overpass1(Diagnostics, pr_overpass1)
      call Get_pr_overpass2(Diagnostics, pr_overpass2)
      call Get_pr_overpass3(Diagnostics, pr_overpass3)
      call Get_pr_overpass4(Diagnostics, pr_overpass4)

      !----------------------
      ! Initialize ASCII file
      !----------------------

      if (pr_ascii) call initializeOutputASCII(gmiGrid, gmiDomain, Diagnostics,&
     &                                         Chemistry, speciesConcentration)

      !------------------------
      ! Initialize netCDF files
      !------------------------

      if (pr_col_diag .or. pr_netcdf .or. pr_restart) then

         if (pr_overpass1) call initializeOutputOverpass1(gmiGrid, gmiDomain,  &
     &                             Diagnostics, Chemistry, metFields)
         if (pr_overpass2) call initializeOutputOverpass2(gmiGrid, gmiDomain,  &
     &                             Diagnostics, Chemistry, metFields)
         if (pr_overpass3) call initializeOutputOverpass3(gmiGrid, gmiDomain,  &
     &                             Diagnostics, Chemistry, metFields)
         if (pr_overpass4) call initializeOutputOverpass4(gmiGrid, gmiDomain,  &
     &                             Diagnostics, Chemistry, metFields)

         if (pr_const) call initializeOutputConcentration(gmiGrid, gmiDomain,  &
     &                           Diagnostics, Chemistry, SpeciesConcentration, &
     &                           Emission, metFields)

         if (pr_AerDust .and. chem_mecha == 'gocart_aerosol') then 

            call initializeOutputAerosolDust(gmiGrid, gmiDomain, Diagnostics, metFields)
         end if

         if ((chem_mecha == 'troposphere').or.(chem_mecha == 'strat_trop')  &
     &         .or. chem_mecha == 'strat_trop_aerosol' ) then

            call Get_phot_opt(Chemistry, phot_opt)
            if (phot_opt == 3) then
               if (pr_AerDust) call initializeOutputAerosolDust(gmiGrid,     &
     &                                    gmiDomain, Diagnostics, metFields)
            end if
         end if

         if (pr_freq1) call initializeOutputFreq1(gmiGrid, gmiDomain,          &
     &                            Diagnostics, Chemistry, metFields, chem_mecha)
         if (pr_freq2) call initializeOutputFreq2(gmiGrid, gmiDomain,          &
     &                            Diagnostics, Chemistry, metFields, chem_mecha)
         if (pr_freq3) call initializeOutputFreq3(gmiGrid, gmiDomain,          &
     &                            Diagnostics, Chemistry, metFields, chem_mecha)
         if (pr_freq4) call initializeOutputFreq4(gmiGrid, gmiDomain,          &
     &                            Diagnostics, Chemistry, metFields, chem_mecha)

         if (pr_qj) call initializeOutputQj(gmiGrid, gmiDomain, Diagnostics,   &
     &                          Chemistry, metFields)

         if (pr_qk) call initializeOutputQk(gmiGrid, gmiDomain, Diagnostics,   &
     &                          Chemistry, metFields)

         if (pr_qqjk) call initializeOutputQqjk(gmiGrid, gmiDomain,            &
     &                          Diagnostics, Chemistry, metFields)

         if (pr_sad) call initializeOutputSAD(Diagnostics, gmiGrid, gmiDomain, &
     &                          Chemistry, metFields)

         if (pr_tend) call initializeOutputTendencies(gmiGrid, gmiDomain,      &
     &                          Diagnostics, Chemistry, metFields)

         if (pr_restart) call initializeOutputRestart(SpeciesConcentration,    &
     &                           gmiGrid, gmiDomain, Diagnostics, Chemistry,   &
     &                           metFields)

         if (pr_flux) call initializeOutputFlux (gmiGrid, gmiDomain,           &
     &                             Diagnostics, Chemistry, metFields)

         if (pr_col_diag) call initializeOutputColumn (gmiGrid, gmiDomain,     &
     &                             Diagnostics, Chemistry, metFields)

#ifdef GTmodule
         if (pr_cloud) call initializeOutputCloudGT(gmiGrid, gmiDomain, &
     &                                           Diagnostics, metFields)
#endif
      end if

      return

      end subroutine initializeOutputFiles
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: controlOutputFiles
!
! !INTERFACE:
!
      subroutine controlOutputFiles (last_tstp, Chemistry, Deposition, Emission, &
     &              SpeciesConcentration, Diagnostics, gmiDomain, gmiGrid,     &
     &              gmiClock, metFields, Advection, chem_mecha)
!
! !USES:
      use GmiControlColumn_mod   , only : controlOutputColumn
      use GmiControlRestart_mod  , only : controlOutputRestart
      use GmiControlASCII_mod    , only : controlOutputASCII
      use GmiControlAerosolDust_mod  , only : controlOutputAerosolDust
      use GmiControlQj_mod           , only : controlOutputQj
      use GmiControlQk_mod           , only : controlOutputQk
      use GmiControlSAD_mod         , only : controlOutputSAD
      use GmiControlQqjk_mod         , only : controlOutputQqjk
      use GmiControlTendencies_mod   , only : controlOutputTendencies
      use GmiControlConcentration_mod, only : controlOutputConcentration
      use GmiControlOverpass1_mod, only : controlOutputOverpass1
      use GmiControlOverpass2_mod, only : controlOutputOverpass2
      use GmiControlOverpass3_mod, only : controlOutputOverpass3
      use GmiControlOverpass4_mod, only : controlOutputOverpass4
      use GmiControlFreq1_mod    , only : controlOutputFreq1
      use GmiControlFreq2_mod    , only : controlOutputFreq2
      use GmiControlFreq3_mod    , only : controlOutputFreq3
      use GmiControlFreq4_mod    , only : controlOutputFreq4
      use GmiControlCloudModuleGT_mod  , only : controlOutputCloudGT
      use GmiControlFlux_mod     , only : controlOutputFlux
!
      implicit none
!
! !INPUT PARAMETERS:
      character(len=*)   , intent(in) :: chem_mecha
      logical            , intent(in) :: last_tstp ! last time step?
      type(t_metFields ) , intent(in) :: metFields
      type(t_gmiDomain ) , intent(in) :: gmiDomain
      type(t_gmiGrid   ) , intent(in) :: gmiGrid  
      type(t_gmiClock  ) , intent(in) :: gmiClock  
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_Chemistry   ), intent(inOut) :: Chemistry
      type(t_Emission    ), intent(inOut) :: Emission 
      type(t_Deposition  ), intent(inOut) :: Deposition
      type(t_Advection   ), intent(inOut) :: Advection 
      type (t_Diagnostics), intent(inOut) :: Diagnostics
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !DESCRIPTION:
! Controls the output.
!
! !LOCAL VARIABLES:
      integer :: procID, phot_opt
      logical :: do_ftiming, pr_diag
      logical :: pr_const, pr_sad, pr_qj, pr_qk, pr_qqjk, pr_cloud,            &
     &           pr_AerDust, pr_flux, pr_time, pr_overpass1, pr_overpass2,     &
     &           pr_overpass3, pr_overpass4, pr_restart, pr_freq1, pr_freq2,   &
     &           pr_freq3, pr_freq4, pr_ascii, pr_netcdf, pr_col_diag, pr_tend
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write (6,*) 'controlOutputFiles called by ', procID

      call Get_pr_qj       (Diagnostics, pr_qj)
      call Get_pr_qk       (Diagnostics, pr_qk)
      call Get_pr_sad      (Diagnostics, pr_sad)
      call Get_pr_qqjk     (Diagnostics, pr_qqjk)
      call Get_pr_tend     (Diagnostics, pr_tend)
      call Get_pr_flux     (Diagnostics, pr_flux)
      call Get_pr_time     (Diagnostics, pr_time)
      call Get_pr_ascii    (Diagnostics, pr_ascii)
      call Get_pr_const    (Diagnostics, pr_const)
      call Get_pr_cloud    (Diagnostics, pr_cloud)
      call Get_pr_freq1    (Diagnostics, pr_freq1)
      call Get_pr_freq2    (Diagnostics, pr_freq2)
      call Get_pr_freq3    (Diagnostics, pr_freq3)
      call Get_pr_freq4    (Diagnostics, pr_freq4)
      call Get_pr_netcdf   (Diagnostics, pr_netcdf)
      call Get_pr_AerDust  (Diagnostics, pr_AerDust)
      call Get_pr_restart  (Diagnostics, pr_restart)
      call Get_pr_col_diag (Diagnostics, pr_col_diag)
      call Get_pr_overpass1(Diagnostics, pr_overpass1)
      call Get_pr_overpass2(Diagnostics, pr_overpass2)
      call Get_pr_overpass3(Diagnostics, pr_overpass3)
      call Get_pr_overpass4(Diagnostics, pr_overpass4)

      call Get_do_ftiming        (Diagnostics, do_ftiming)

!     ==========================
      call controlScreenOutput (gmiClock, pr_time, pr_diag, procID)
!     ==========================

      if (pr_ascii) call controlOutputASCII(SpeciesConcentration, Chemistry, &
     &                                gmiClock, metFields)

      if (do_ftiming) call Ftiming_On  ('gmiWritingOutput')

      if (pr_netcdf) then
         if (pr_col_diag) call controlOutputColumn (gmiDomain,                 &
     &                                SpeciesConcentration, metFields, gmiClock)

         if (pr_const) call controlOutputConcentration (last_tstp, Chemistry,  &
     &                             Emission, Deposition, SpeciesConcentration, &
     &                             gmiClock, metFields)
    
         if (pr_qj) call controlOutputQj (last_tstp, Chemistry, gmiClock)
    
         if (pr_qk) call controlOutputQk (last_tstp, Chemistry, gmiClock)

         if (pr_sad) call controlOutputSAD (last_tstp, Chemistry, gmiClock)
    
         if (pr_qqjk) call controlOutputQqjk (last_tstp, Diagnostics,     &
     &                           Chemistry, gmiClock)
    
         if (pr_tend) call controlOutputTendencies(last_tstp,             &
     &                            SpeciesConcentration, gmiClock)

         if (pr_AerDust .and. chem_mecha == 'gocart_aerosol') then 
            call controlOutputAerosolDust (last_tstp, Chemistry, gmiClock, metFields)
         end if

         if ((chem_mecha == 'troposphere').or.(chem_mecha == 'strat_trop')  &
     &         .or. chem_mecha == 'strat_trop_aerosol' ) then
            call Get_phot_opt(Chemistry, phot_opt)
            if (phot_opt == 3) then
               if (pr_AerDust) then
                  call controlOutputAerosolDust (last_tstp, Chemistry, &
     &                        gmiClock, metFields)
               end if
            end if
         end if
!         
         if (pr_restart) call controlOutputRestart (last_tstp, Chemistry,      & 
     &                               SpeciesConcentration, gmiClock, metFields)

         if (pr_flux) call controlOutputFlux(last_tstp, Diagnostics, gmiClock, &
     &                            Advection)

         if (pr_freq1) call controlOutputFreq1 (last_tstp, Chemistry,          &
     &                             SpeciesConcentration, gmiClock, metFields)

         if (pr_freq2) call controlOutputFreq2 (last_tstp, Chemistry,          &
     &                          SpeciesConcentration, gmiClock, metFields)

         if (pr_freq3) call controlOutputFreq3 (last_tstp, Chemistry,          &
     &                             SpeciesConcentration, gmiClock, metFields)

         if (pr_freq4) call controlOutputFreq4 (last_tstp, Chemistry,          &
     &                             SpeciesConcentration, gmiClock, metFields)

#ifdef GTmodule
         if (pr_cloud) call controlOutputCloudGT(last_tstp, Chemistry,         &
     &                             gmiClock, metFields)
#endif

         if (pr_overpass1) call controlOutputOverpass1 (last_tstp, Chemistry,  &
              & SpeciesConcentration, Diagnostics, gmiClock, &
              & Emission, metFields)

         if (pr_overpass2) call controlOutputOverpass2 (last_tstp, Chemistry,  &
              &                   SpeciesConcentration, Diagnostics, gmiClock, &
              & Emission, metFields)

         if (pr_overpass3) call controlOutputOverpass3 (last_tstp, Chemistry,  &
              &                   SpeciesConcentration, Diagnostics, gmiClock, &
              & Emission, metFields)

         if (pr_overpass4) call controlOutputOverpass4 (last_tstp, Chemistry,  &
              &                   SpeciesConcentration, Diagnostics, gmiClock, &
              & Emission, metFields)

      end if

      if (do_ftiming) call Ftiming_Off  ('gmiWritingOutput')

      return

      end subroutine controlOutputFiles
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: controlScreenOutput
!
! !INTERFACE:
!
      subroutine controlScreenOutput (gmiClock, pr_time, pr_diag, procID)
!
! !USES:
      use GmiTimeControl_mod, only : Get_curGmiDate, Get_curGmiTime,           &
     &       Get_gmiSeconds, Get_gmiTimeStep
!
      implicit none

! !INPUT PARAMETERS:
      logical, intent(in) :: pr_time, pr_diag
      integer, intent(in) :: procID
      type (t_gmiClock), intent(in) :: gmiClock
!
! !DESCRIPTION:
! Controls the screen output.
!
! !LOCAL VARIABLES:
      character (len=8)  :: cdummy
      character (len=10) :: chms
      integer :: idays, nymd, nhms, ndt
      real*8  :: rsecpdy
      real*8  :: gmi_sec, tdt
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'controlScreenOutput called by ', procID

      call Get_curGmiDate (gmiClock, nymd   )
      call Get_curGmiTime (gmiClock, nhms   )
      call Get_gmiSeconds (gmiClock, gmi_sec)
      call Get_gmiTimeStep(gmiClock, tdt    )

      ndt = tdt

!     ---------------------
!     Write time to screen.
!     ---------------------

      call Date_And_Time (cdummy, chms)

      rsecpdy = SECPDY

      if (Mod (Nint (gmi_sec), Nint (rsecpdy)) < ndt) then

!       ------------------------------------
!       Print time at start of each new day.
!       ------------------------------------

        idays = Nint (gmi_sec / rsecpdy)

        IF (procID == 0) THEN
           Write (6,900) idays, nymd, nhms, chms(1:6)
        END IF

      else if (pr_time) then

!       --------------------------
!       Print time each time step.
!       --------------------------

        IF (procID == 0) THEN
           Write (6,910) Nint (gmi_sec), nymd, nhms, chms(1:6)
        END IF

      end if


 900  format (1x, '--day, ymd, hms, wc/hms:  ',  &
     &        i10, ', ', i8.8, ', ', i6.6, ', ', a6)

 910  format (1x, '**sec, ymd, hms, wc/hms:  ',  &
     &        i10, ', ', i8.8, ', ', i6.6, ', ', a6)

      return

      end subroutine controlScreenOutput
!EOC
!------------------------------------------------------------------------------
      end module GmiControlOutput_mod
