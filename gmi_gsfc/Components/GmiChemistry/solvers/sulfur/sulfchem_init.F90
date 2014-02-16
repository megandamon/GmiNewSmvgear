!=============================================================================
!
! $Id: sulfchem_init.F90,v 1.4 2013-09-03 16:16:43 jkouatch Exp $
!
! CODE DEVELOPER
!   Xiaohong Liu
!
! FILE
!   sulfchem_init.F
!
! ROUTINES
!   Do_Sulf_Init
!
! ARGUMENTS
!   pr_diag        : print some diagnostic output to screen?
!   loc_proc       : local processor #
!   qj_infile_name : qj input file name
!   qj_var_name    : qj variable name
!
!=============================================================================

      subroutine Do_Sulf_Init (sulfSavedVars, pr_diag, loc_proc, &
                               nlat, qj_infile_name, qj_var_name)

      use GmiSolver_SavedVariables_mod, only : t_SulfSaved

      use m_netcdf_io_open , only : Ncop_Rd
      use m_netcdf_io_close, only : Nccl
      use m_netcdf_io_read , only : Ncrd_3d

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical :: pr_diag
      integer :: loc_proc, nlat
      character*(*) :: qj_infile_name
      character*(*) :: qj_var_name
      type(t_SulfSaved), intent(inOut) :: sulfSavedVars

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ncid_qj

      integer :: cnt3d (3)
      integer :: strt3d(3)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Do_Sulf_Init called by ', loc_proc
      end if

      sulfSavedVars%lon_dummy(1) = 0.0d0

      allocate(sulfSavedVars%qjh2o2_2d(nlat,5,12))

      call Ncop_Rd (ncid_qj, qj_infile_name)

      strt3d(:) = (/ 1, 1, 1/)

      cnt3d (:) = (/ nlat, 5, 12 /)

      call Ncrd_3d (sulfSavedVars%qjh2o2_2d, ncid_qj, qj_var_name, strt3d, cnt3d)

      call Nccl (ncid_qj)


      return

      end

