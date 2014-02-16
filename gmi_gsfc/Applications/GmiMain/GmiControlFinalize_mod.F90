#include "GmiESMF_ErrLog.h"
!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE:
!
      module GmiControlFinalize_mod
!
! !USES:
      use ESMF_Mod
      use GmiESMF_ErrorChecking_mod
      use GmiESMFderivedType_mod       , only : t_gmiESMF
      use Ftiming_Dao
      use GmiFileUnit_mod, only : GetFileUnitNumber

      use GmiPrintError_mod, only : GmiPrintError
      use m_netcdf_io_close, only : Nccl_Noerr
      use GmiMessagePassing_mod, only : synchronizeGroup, broadcastInteger
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_iAmRootProc,     &
     &       Get_communicatorWorld, Get_procID, Get_rootProc, Get_numDomains,  &
     &       Get_numLatDomains, Get_numLonDomains

      use GmiControlASCII_mod      , only : finalizeOutputASCII
      use GmiControlRestart_mod    , only : finalizeOutputRestart
      use GmiControlQj_mod         , only : finalizeOutputQj
      use GmiControlQk_mod         , only : finalizeOutputQk
      use GmiControlSAD_mod        , only : finalizeOutputSAD
      use GmiControlQqjk_mod       , only : finalizeOutputQqjk
      use GmiControlFlux_mod       , only : finalizeOutputFlux
      use GmiControlFreq1_mod      , only : finalizeOutputFreq1
      use GmiControlFreq2_mod      , only : finalizeOutputFreq2
      use GmiControlFreq3_mod      , only : finalizeOutputFreq3
      use GmiControlFreq4_mod      , only : finalizeOutputFreq4
      use GmiControlCloudModuleGT_mod, only : finalizeOutputCloudGT
      use GmiControlAerosolDust_mod, only : finalizeOutputAerosolDust

      use GmiControlOverpass1_mod    , only : finalizeOutputOverpass1
      use GmiControlOverpass2_mod    , only : finalizeOutputOverpass2
      use GmiControlOverpass3_mod    , only : finalizeOutputOverpass3
      use GmiControlOverpass4_mod    , only : finalizeOutputOverpass4
      use GmiControlTendencies_mod   , only : finalizeOutputTendencies

      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Get_do_ftiming,      &
     &       Get_pr_diag, Get_pr_const, Get_pr_restart, Get_pr_AerDust,        &
     &       Get_pr_freq1, Get_pr_freq2, Get_pr_freq3, Get_pr_freq4,           &
     &       Get_pr_qqjk, Get_pr_flux, Get_pr_tend, Get_pr_qj, Get_pr_qk,      &
     &       Get_pr_cloud, Get_pr_const, Get_pr_overpass1, Get_pr_overpass2,   &
     &       Get_pr_overpass3, Get_pr_overpass4, Get_pr_sad, Get_pr_ascii

      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: gmiControlFinalize
!
! !AUTHOR:
! Jules Kouatchou, Jules.Kouatchou-1@nasa.gov
!
!EOP
!=============================================================================
      contains
!=============================================================================
!!BOP
!
! !IROUTINE: gmiControlFinalize
!
! !INTERFACE:
!
      subroutine gmiControlFinalize (advCoreESMF, esmfClock, gmiDomain, Diagnostics)
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_gmiDomain  ), intent(in) :: gmiDomain
      type(t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Clock), intent(inOut) :: esmfClock
      type (t_gmiESMF), intent(inOut) :: advCoreESMF
!
! !DESCRIPTION:
! Finalizes all the components and closes all the files.
!
! !LOCAL VARIABLES:
      integer :: ftim_lun, STATUS, rc
      integer :: ierr
      integer :: commu_world
      logical :: do_ftiming, iAmRootProc
      integer :: procID, rootProc
      integer :: numDomains, numLatDomains, numLonDomains
      logical :: pr_diag, pr_restart, pr_AerDust, pr_freq1, pr_freq2, pr_freq3,&
                 pr_freq4, pr_qqjk, pr_flux, pr_tend, pr_qj, pr_qk, pr_cloud,  &
     &           pr_overpass1, pr_overpass2, pr_overpass3, pr_overpass4,       &
     &           pr_const, pr_sad, pr_ascii
      integer :: lun_array(1)
      character(len=ESMF_MAXSTR), parameter :: IAm = "gmiControlFinalize"
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID (gmiDomain, procID      )
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write (6,*) trim(Iam), ' called by ', procID

      call Get_communicatorWorld (gmiDomain, commu_world )

      ! ---------------------------
      ! Print out fine timing data.
      ! ---------------------------

      call Get_do_ftiming(Diagnostics, do_ftiming)

      if (do_ftiming) then

         call Get_rootProc          (gmiDomain, rootProc    )
         call Get_iAmRootProc       (gmiDomain, iAmRootProc )
         call Get_numDomains        (gmiDomain, numDomains  )
         call Get_numLonDomains     (gmiDomain, numLonDomains)
         call Get_numLatDomains     (gmiDomain, numLatDomains)

         call Ftiming_Off ('whole_GMI')

         if (iAmRootProc) then
            call GetFileUnitNumber (ftim_lun, ierr)
         end if

         lun_array(1) = ftim_lun
         call broadcastInteger (lun_array, 1, rootProc, commu_world)
         ftim_lun     = lun_array(1)

         call Ftiming_Prt (ftim_lun, numDomains, numLonDomains, &
     &                      numLatDomains, procID, commu_world)
      end if

      ! ------------------------------------------------------
      ! Synchronize all the processors before exiting.
      ! ------------------------------------------------------

      call synchronizeGroup (commu_world)

      ! -----------------------------
      ! Close the netCDF output files
      ! -----------------------------

      call Get_pr_qj(Diagnostics, pr_qj)
      call Get_pr_qk(Diagnostics, pr_qk)
      call Get_pr_sad(Diagnostics, pr_sad)
      call Get_pr_qqjk(Diagnostics, pr_qqjk)
      call Get_pr_tend(Diagnostics, pr_tend)
      call Get_pr_ascii(Diagnostics, pr_ascii)
      call Get_pr_freq1(Diagnostics, pr_freq1)
      call Get_pr_freq2(Diagnostics, pr_freq2)
      call Get_pr_freq3(Diagnostics, pr_freq3)
      call Get_pr_freq4(Diagnostics, pr_freq4)
      call Get_pr_cloud(Diagnostics, pr_cloud)
      call Get_pr_const(Diagnostics, pr_const)
      call Get_pr_AerDust(Diagnostics, pr_AerDust)
      call Get_pr_restart(Diagnostics, pr_restart)
      call Get_pr_overpass1(Diagnostics, pr_overpass1)
      call Get_pr_overpass2(Diagnostics, pr_overpass2)
      call Get_pr_overpass3(Diagnostics, pr_overpass3)
      call Get_pr_overpass4(Diagnostics, pr_overpass4)

      if (pr_qj)        call finalizeOutputQj()
      if (pr_qk)        call finalizeOutputQk()
      if (pr_sad)       call finalizeOutputSAD()
      if (pr_qqjk)      call finalizeOutputQqjk()
      if (pr_flux)      call finalizeOutputFlux()
      if (pr_tend)      call finalizeOutputTendencies()
      if (pr_ascii)     call finalizeOutputASCII()
      if (pr_freq1)     call finalizeOutputFreq1()
      if (pr_freq2)     call finalizeOutputFreq2()
      if (pr_freq3)     call finalizeOutputFreq3()
      if (pr_freq4)     call finalizeOutputFreq4()
      if (pr_restart)   call finalizeOutputRestart()
      if (pr_AerDust)   call finalizeOutputAerosolDust()
      if (pr_overpass1) call finalizeOutputOverpass1()
      if (pr_overpass2) call finalizeOutputOverpass2()
      if (pr_overpass3) call finalizeOutputOverpass3()
      if (pr_overpass4) call finalizeOutputOverpass4()

#ifdef GTmodule
      if (pr_cloud) call finalizeOutputCloudGT()
#endif

     !call ESMF_GridCompFinalize(advCoreESMF%compGridded, &
     !           importState = advCoreESMF%stateImp, & 
     !           exportState = advCoreESMF%stateExp, & 
     !           clock       = esmfClock,            &
     !           rc          = STATUS)
     !VERIFY_(STATUS)

      return

      end subroutine gmiControlFinalize
!EOC
!------------------------------------------------------------------------------
      end module GmiControlFinalize_mod
