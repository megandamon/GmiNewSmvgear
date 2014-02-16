!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: ReadNOxScalingFactor_mod
!
! !INTERFACE:
!
module ReadNOxScalingFactor_mod
!
! !USES:
  use m_netcdf_io_open , only : Ncop_Rd
  use m_netcdf_io_close, only : Nccl
  use m_netcdf_io_read , only : Ncrd_3d, Ncrd_3d

  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
   private
   public  :: readNOffScalingFactor
   public  :: readNObbScalingFactor

#  include "GmiParameters.h"
!
! !DESCRIPTION:
! Routines to read in (every hour) scaling factors for fossil fuel
! NO and biomass burning NO emissions.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  19 Sept2007 Initial code.
!EOP
!-------------------------------------------------------------------------
    contains
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readNOffScalingFactor
!
! !INTERFACE:
!
  subroutine readNOffScalingFactor(scFacNOff, scFactorNOff_infile_name, &
                         i1, i2, ju1, j2,  i1_gl, ju1_gl)

  implicit none
!
! !INPUT PARAMETERS:
  integer, intent(in) :: i1, i2, ju1, j2
  integer, intent(in) :: i1_gl, ju1_gl
  character (len=MAX_LENGTH_FILE_NAME), intent(in) :: scFactorNOff_infile_name
!
! !OUTPUT PARAMETERS:
  real*8 , intent(out) :: scFacNOff(i1:i2, ju1:j2, 24)
!
! !DESCRIPTION:
!  Reads in scaling factor for anthropogenic NOx emissions.
!
! !LOCAL VARIABLES:
  integer :: il, ij, inb, jnb
  integer :: ncid
  integer :: count3d(3), start3d(3)
  real*8  :: temp_todn(i2-i1+1,j2-ju1+1, 24)
!
! !REVISION HISTORY:
!  19Sep2007
!EOP
!---------------------------------------------------------------------------
!BOC
  start3d(:) = (/i1-i1_gl+1, ju1-ju1_gl+1, 1  /)
  count3d(:) = (/i2-i1+1   , j2-ju1+1    , 24 /)

  call Ncop_Rd (ncid, scFactorNOff_infile_name)

  call Ncrd_3d (temp_todn, ncid, 'todn', start3d, count3d)

  do ij = ju1, j2
     jnb = ij - ju1 + 1
     do il = i1, i2
        inb = il - i1 + 1
        scFacNOff(il,ij,:) = temp_todn(inb,jnb,:)
     end do
  end do

  call Nccl (ncid)

  return

  end subroutine readNOffScalingFactor
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readNObbScalingFactor
!
! !INTERFACE:
!
  subroutine readNObbScalingFactor(scFacNObb, scFactorNObb_infile_name, &
                         i1, i2, ju1, j2,  i1_gl, ju1_gl)

  implicit none
!
! !INPUT PARAMETERS:
  integer, intent(in) :: i1, i2, ju1, j2
  integer, intent(in) :: i1_gl, ju1_gl
  character (len=MAX_LENGTH_FILE_NAME), intent(in) :: scFactorNObb_infile_name
!
! !OUTPUT PARAMETERS:
  real*8 , intent(out) :: scFacNObb(i1:i2, ju1:j2, 24)
!
! !DESCRIPTION:
!  Reads in scaling factor for anthropogenic NOx emissions.
!
! !LOCAL VARIABLES:
  integer :: il, ij, inb, jnb
  integer :: ncid
  integer :: count3d(3), start3d(3)
  real*8  :: temp_todn(i2-i1+1,j2-ju1+1, 24)
!
! !REVISION HISTORY:
!  19Sep2007
!EOP
!---------------------------------------------------------------------------
!BOC
  start3d(:) = (/i1-i1_gl+1, ju1-ju1_gl+1, 1  /)
  count3d(:) = (/i2-i1+1   , j2-ju1+1    , 24 /)

  call Ncop_Rd (ncid, scFactorNObb_infile_name)

  call Ncrd_3d (temp_todn, ncid, 'todn', start3d, count3d)

  do ij = ju1, j2
     jnb = ij - ju1 + 1
     do il = i1, i2
        inb = il - i1 + 1
        scFacNObb(il,ij,:) = temp_todn(inb,jnb,:)
     end do
  end do

  call Nccl (ncid)

  return

  end subroutine readNObbScalingFactor
!EOC
!-------------------------------------------------------------------------
end module ReadNOxScalingFactor_mod
