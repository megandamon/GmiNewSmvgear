

    module ReadLiqBinarySulfate_mod

    use m_netcdf_io_open , only : Ncop_Rd
    use m_netcdf_io_read , only : Ncrd_4d, Ncrd_3d
    use m_netcdf_io_close, only : Nccl
    use GmiCheckRange_mod, only : checkRange4d

    implicit none

    private
    public  :: ReadLiqBinarySulfate

#   include "GmiParameters.h"

    CONTAINS

!=============================================================================
!
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   ReadLiqBinarySulfate
!
!-----------------------------------------------------------------------------

      subroutine ReadLiqBinarySulfate (lbssad, lbssad_init_val, lbssad_opt, &
                 lbssad_infile_name, loc_proc, pr_diag, &
                 i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl, lbssad_timpyr)
!
! !USES:
!
      implicit none
!
! !INPUT PARAMETERS:
!
      integer          , intent(in)   :: i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl
      integer          , intent(in)   :: lbssad_timpyr
      integer          , intent(in)   :: lbssad_opt, loc_proc
      real*8           , intent(in)   :: lbssad_init_val
      logical          , intent(in)   :: pr_diag
      character (len=*), intent(in)   :: lbssad_infile_name
!
! !OUTPUT PARAMETERS:
!
      real*8 , intent(out)  :: lbssad(i1:i2, ju1:j2, k1:k2, lbssad_timpyr)
!              ! liquid binary sulfate background surface area density (cm^-1) 

! !DESCRIPTION:
!  This routine sets lbssad, the liquid binary sulfate background surface
!  area density (cm^-1).
!
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!

      if (pr_diag) then
        Write (6,*) 'ReadLiqBinarySulfate called by ', loc_proc
      end if


!     ====================
      if (lbssad_opt == 1) then
!     ====================

        lbssad(i1:i2,ju1:j2,:,:) = lbssad_init_val


!     =========================
      else if (lbssad_opt == 2) then
!     =========================

        call ReadLbssad3 (lbssad, lbssad_infile_name, &
              i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl, lbssad_timpyr)


!     =========================
      else if (lbssad_opt == 3) then
!     =========================

        call ReadLbssad2 (lbssad, lbssad_infile_name, &
              i1, i2, ju1, j2, k1, k2, ju1_gl, lbssad_timpyr)


!     ======
      end if
!     ======

      call checkRange4d  &
     &  ('lbssad', loc_proc, i1, i2, ju1, j2, k1, k2,  &
     &   1, lbssad_timpyr, lbssad, 0.0d0, 1.0d0)

      return

      end subroutine ReadLiqBinarySulfate


!-----------------------------------------------------------------------------
!
! ROUTINE
!   ReadLbssad2
!
! ARGUMENTS
!   lbssad             : liquid binary sulfate background surface area
!                        density (cm^-1)
!   lbssad_infile_name : lbssad input file name
!
!-----------------------------------------------------------------------------

      subroutine ReadLbssad2 (lbssad, lbssad_infile_name, &
                  i1, i2, ju1, j2, k1, k2, ju1_gl, lbssad_timpyr)
!
      implicit none
!
! !INPUT PARAMETERS:
!
      integer          , intent(in)   :: i1, i2, ju1, j2, k1, k2, ju1_gl
      integer          , intent(in)   :: lbssad_timpyr
      character (len=*), intent(in)   :: lbssad_infile_name
!
! !OUTPUT PARAMETERS:
!
      real*8 , intent(out)  :: lbssad(i1:i2, ju1:j2, k1:k2, lbssad_timpyr)

!     -----------------------
!     Parameter declarations.
!     -----------------------

      character (len=MAX_LENGTH_VAR_NAME), parameter :: LBSSAD_VNAM = 'lbssad'

! !DESCRIPTION:
!   This routine reads in a set of lbssad data, the liquid binary sulfate
!   background surface area density (cm^-1).  Recall that all NetCDF arrays
!   always start at subscript 1, and this is not always the case for Gmimod
!   dimensions.
!
! !LOCAL VARIABLES:
!
      integer :: il, ij, ik, it, jnb, ncid_lb
      integer :: cnt3d (3), strt3d(3)
      real*8, allocatable  :: lbssad1(:, :, :)
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
!... read in only zonal avg fields
      strt3d(1) = ju1 - ju1_gl + 1
      strt3d(2) = k1

      cnt3d (:) = (/ j2-ju1+1, k2-k1+1, 1 /)

      call Ncop_Rd (ncid_lb, lbssad_infile_name)

!... read in zonal average aerosol files
      allocate (lbssad1( j2-ju1+1, k2-k1+1, 1 ))
      do it = 1, lbssad_timpyr
        strt3d(3) = it
        call Ncrd_3d (lbssad1, ncid_lb, LBSSAD_VNAM, strt3d, cnt3d)
!... put in 3d array
        do ik = k1, k2
          do ij = ju1, j2
            jnb = ij - ju1 + 1
            do il = i1, i2
              lbssad(il,ij,ik,it) = lbssad1(jnb,ik,1)
            end do
          end do
        end do
      end do
      deallocate(lbssad1)

      call Nccl (ncid_lb)

      return

      end subroutine ReadLbssad2


!-----------------------------------------------------------------------------
!
! ROUTINE
!   ReadLbssad3
!
! ARGUMENTS
!   lbssad             : liquid binary sulfate background surface area
!                        density (cm^-1)
!   lbssad_infile_name : lbssad input file name
!
!-----------------------------------------------------------------------------

      subroutine ReadLbssad3 (lbssad, lbssad_infile_name, &
                  i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl, lbssad_timpyr)
!
      implicit none
!
! !INPUT PARAMETERS:
!
      integer          , intent(in)   :: i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl
      integer          , intent(in)   :: lbssad_timpyr
      character (len=*), intent(in)   :: lbssad_infile_name
!
! !OUTPUT PARAMETERS:
!
      real*8 , intent(out)  :: lbssad(i1:i2, ju1:j2, k1:k2, lbssad_timpyr)

!     -----------------------
!     Parameter declarations.
!     -----------------------

      character (len=MAX_LENGTH_VAR_NAME), parameter :: LBSSAD_VNAM = 'lbssad'
!
! !DESCRIPTION:
!  This routine reads in a set of lbssad data from full 3-d file,
!  the liquid binary sulfatebackground surface area density (cm^-1).
!  Recall that all NetCDF arrays always start at subscript 1,
!  and this is not always the case for Gmimod dimensions.
!
! !LOCAL VARIABLES:
!
      integer :: il, ij, ik, it, inb, jnb, ncid_lb
      integer :: cnt4d (4), strt4d(4)
      real*8, allocatable  :: lbssad1(:, :, :, :)
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!

      strt4d(1) = i1  -  i1_gl + 1
      strt4d(2) = ju1 - ju1_gl + 1
      strt4d(3) = k1

      cnt4d (:) = (/ i2-i1+1, j2-ju1+1, k2-k1+1, 1 /)


      call Ncop_Rd (ncid_lb, lbssad_infile_name)

!... read in full 3d aerosol files
      allocate (lbssad1(i2-i1+1, j2-ju1+1, k2-k1+1, 1 ))
      do it = 1, lbssad_timpyr
        strt4d(4) = it
        call Ncrd_4d (lbssad1, ncid_lb, LBSSAD_VNAM, strt4d, cnt4d)
!... put in 4d array
        do ik = k1, k2
          do ij = ju1, j2
            jnb = ij - ju1 + 1
            do il = i1, i2
              inb = il - i1 + 1
              lbssad(il,ij,ik,it) = lbssad1(inb,jnb,ik,1)
            end do
          end do
        end do
      end do
      deallocate(lbssad1)
!... close netCDF file
      call Nccl (ncid_lb)

      return

      end subroutine ReadLbssad3

    end module ReadLiqBinarySulfate_mod
