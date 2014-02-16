!=======================================================================
!
! FILE
!   setkin_surf_bc.h
!
! DESCRIPTION
!   This include file contains information about surface boundary
!   conditions.
!
!
!=======================================================================

!     -----------------------
!     Parameter declarations.
!     -----------------------

!     -------------------------------------------------------
!     NUM_SBC : number of surface bounday conditions to reset
!     K_SBC   : max k to which species will be adjusted
!     -------------------------------------------------------

      integer, parameter :: NUM_SBC = 0

      integer, parameter :: K_SBC   = 2

!     ------------------
!     Integer variables.
!     ------------------

      integer, save :: sbc_map(NUM_SBC) = &
     &  (/ 0 /)

!                                  --^--

