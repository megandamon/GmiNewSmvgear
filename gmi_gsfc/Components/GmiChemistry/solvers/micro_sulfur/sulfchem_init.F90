!=============================================================================
!
! $Id: sulfchem_init.F90,v 1.3 2011-08-09 22:12:58 mrdamon Exp $
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

      subroutine Do_Sulf_Init  &
     &  (pr_diag, loc_proc, nlat, qj_infile_name, qj_var_name)

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

      call Ncop_Rd (ncid_qj, qj_infile_name)


      strt3d(:) = (/ 1, 1, 1/)

      cnt3d (:) = (/ nlat, 5, 12 /)


      call Ncrd_3d (qjh2o2_2d, ncid_qj, qj_var_name, strt3d, cnt3d)


      call Nccl (ncid_qj)


      return

      end

