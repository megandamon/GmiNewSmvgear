  module ReadReactionRateAdjFac_mod

      use m_netcdf_io_open , only : Ncop_Rd
      use m_netcdf_io_close, only : Nccl
      use m_netcdf_io_read , only : Ncrd_5d
      use GmiCheckRange_mod, only : checkRange5d

    implicit none

    private
    public  :: readReactionRateAdjustmentFactors

!=============================================================================
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
!=============================================================================

    CONTAINS

!-----------------------------------------------------------------------------
!
! ROUTINE
!   readReactionRateAdjustmentFactors
!
! DESCRIPTION
!   This routine reads in a set of reaction rate adjustment factors data.
!   Recall that all NetCDF arrays always start at subscript 1, and this is
!   not always the case for Gmimod dimensions.
!
! ARGUMENTS
!   rxnr_adjust             : array of reaction rate adjustment factors
!   rxnr_adjust_infile_name : reaction rate adjustment factor input file name
!   rxnr_adjust_var_name    : NetCDF reaction rate adjustment factor variable
!                             name
!
!-----------------------------------------------------------------------------

      subroutine readReactionRateAdjustmentFactors  &
     &  (rxnr_adjust, rxnr_adjust_infile_name, rxnr_adjust_var_name, &
     &   num_rxnr_adjust, rxnr_adjust_timpyr, pr_diag, loc_proc, &
     &   i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, loc_proc, i1_gl, ju1_gl
      integer, intent(in) :: num_rxnr_adjust, rxnr_adjust_timpyr
      real*8 :: rxnr_adjust(i1:i2, ju1:j2, k1:k2,  &
     &                      num_rxnr_adjust, rxnr_adjust_timpyr)
      character (len=*) :: rxnr_adjust_infile_name
      character (len=*) :: rxnr_adjust_var_name

      integer :: il, ij, ik, ic, it, inb, jnb, ncid_rra
      integer :: cnt5d (5), strt5d(5)
      real*8  :: rxnr_adjust_tmp(i2-i1+1, j2-ju1+1, k2-k1+1, 1, 1)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'readReactionRateAdjustmentFactors called by ', loc_proc
      end if

      strt5d(1) = i1  -  i1_gl + 1
      strt5d(2) = ju1 - ju1_gl + 1
      strt5d(3) = k1

      cnt5d (:) = (/ i2-i1+1, j2-ju1+1, k2-k1+1, 1, 1 /)

      call Ncop_Rd (ncid_rra, rxnr_adjust_infile_name)

      do it = 1, rxnr_adjust_timpyr
         strt5d(5) = it

         do ic = 1, num_rxnr_adjust
            strt5d(4) = ic

            call Ncrd_5d (rxnr_adjust_tmp, ncid_rra, rxnr_adjust_var_name,  &
     &                     strt5d, cnt5d)

            do ik = k1, k2
               do ij = ju1, j2
                  jnb = ij - ju1 + 1
                  do il = i1, i2
                     inb = il - i1 + 1
                     rxnr_adjust(il,ij,ik,ic,it) = rxnr_adjust_tmp(inb,jnb,ik,1,1)
                  end do
               end do
            end do
         end do
      end do

      call Nccl (ncid_rra)

!      call Check_Range_5d  &
      call checkRange5d  &
     &  ('rxnr_adjust', loc_proc, i1, i2, ju1, j2, k1, k2, 1,  &
     &   num_rxnr_adjust, 1, rxnr_adjust_timpyr, rxnr_adjust,  &
     &   0.0d0, 1.0d2)


      return

      end subroutine readReactionRateAdjustmentFactors

  end module ReadReactionRateAdjFac_mod
