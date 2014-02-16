    module ReadFixedConcentration_mod

      use m_netcdf_io_open , only : Ncop_Rd
      use m_netcdf_io_close, only : Nccl
      use m_netcdf_io_read , only : Ncrd_5d
      use GmiCheckRange_mod, only : checkRange5d

    implicit none

    private
    public   :: readFixedConcentration

    CONTAINS

!-----------------------------------------------------------------------------
!
! ROUTINE
!   readFixedConcentration
!
! DESCRIPTION
!   This routine reads in a set of fixed cosnt data.  Recall that all NetCDF
!   arrays always start at subscript 1, and this is not always the case for
!   Gmimod dimensions.
!
! ARGUMENTS
!   fixed_const    : array of fixed consts (mixing ratio)
!   fixed_const_infile_name : fixed const input file name
!   const_var_name : NetCDF const variable name
!
!-----------------------------------------------------------------------------

      subroutine readFixedConcentration  &
     &  (fixed_const, fixed_const_infile_name, const_var_name, &
     &   num_fixed_const, fixed_const_timpyr, pr_diag, loc_proc, &
     &   ilong, ilat, ivert, i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

     integer, intent(in) :: num_fixed_const, fixed_const_timpyr
     integer, intent(in) :: i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl
     integer, intent(in) :: ilat, ilong, ivert, loc_proc
     logical, intent(in) :: pr_diag
      real*8  :: fixed_const(i1:i2, ju1:j2, k1:k2,  &
     &                       num_fixed_const, fixed_const_timpyr)
      character (len=*), intent(in) :: fixed_const_infile_name
      character (len=*), intent(in) :: const_var_name

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij, ik, ic, it
      integer :: inb, jnb
      integer :: ncid_fc

      integer :: cnt5d (5)
      integer :: strt5d(5)

      real*8  :: fconst_tmp(ilong, ilat, ivert, 1, 1)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'readFixedConcentration called by ', loc_proc
      end if

      strt5d(1) = i1  -  i1_gl + 1
      strt5d(2) = ju1 - ju1_gl + 1
      strt5d(3) = k1

      cnt5d (:) = (/ ilong, ilat, ivert, 1, 1 /)

      call Ncop_Rd (ncid_fc, fixed_const_infile_name)

      do it = 1, fixed_const_timpyr

        strt5d(5) = it

        do ic = 1, num_fixed_const

          strt5d(4) = ic

          call Ncrd_5d (fconst_tmp, ncid_fc, const_var_name, strt5d, cnt5d)

          do ik = k1, k2
            do ij = ju1, j2
              jnb = ij - ju1 + 1
              do il = i1, i2
                inb = il - i1 + 1

                fixed_const(il,ij,ik,ic,it) = fconst_tmp(inb,jnb,ik,1,1)

              end do
            end do
          end do

        end do

      end do

      call Nccl (ncid_fc)

!      call Check_Range_5d  &
      call checkRange5d  &
     &  ('fixed_const', loc_proc, i1, i2, ju1, j2, k1, k2, 1,  &
     &   num_fixed_const, 1, fixed_const_timpyr, fixed_const,  &
     &   0.0d0, 1.0d20)

      return

      end subroutine readFixedConcentration
    end module ReadFixedConcentration_mod
