!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
       module GmiReadGocartSourceFiles_mod

    use m_netcdf_io_open , only : Ncop_Rd
    use m_netcdf_io_read , only : Ncrd_4d, Ncrd_2d
    use m_netcdf_io_close, only : Nccl

       implicit none

       private
       public   :: GmiReadGocartSourceFiles

#     include "GmiParameters.h"

       contains

!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GmiReadGocartSourceFiles
!
! !INTERFACE:
!
  subroutine GmiReadGocartSourceFiles           &
                    (GOCARTerod_infile_name, GOCARTerod_mod_infile_name, &
                     GOCARTocean_infile_name, mcor, latdeg, londeg, &
                     i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl)
!
! !USES:
!
  use GmiDustMethod_mod    , only : ndcls, nDustBin, ndsrc, erod, erod_mod
  use GmiDustMethod_mod    , only : Allocate_erod, Allocate_erod_mod
  use GmiDustMethod_mod    , only : Allocate_srcEmissDust, Update_erod
  use GmiWaterMethod_mod   , only : water, Allocate_water, Update_water
  use GmiSeaSaltMethod_mod , only : Allocate_srcEmissSeaSalt

  implicit none

! !INPUT PARAMETERS:
  integer            , intent(in) :: i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl
  real*8             , intent(in) :: mcor(i1:i2,ju1:j2)
  real*8             , intent(in) :: latdeg(ju1_gl:j2_gl), londeg(i1_gl:i2_gl)
  character (len=MAX_LENGTH_FILE_NAME), intent(in) :: GOCARTerod_infile_name 
  character (len=MAX_LENGTH_FILE_NAME), intent(in) :: GOCARTerod_mod_infile_name 
  character (len=MAX_LENGTH_FILE_NAME), intent(in) :: GOCARTocean_infile_name    
!
! !DESCRIPTION: This routine reads in netCDF files for GOCART dust and sea
!  salt calculations. All the worker processors do the reading.
!
! !LOCAL VARIABLES:
  character (len=75) :: err_msg

  integer :: jlen, ilen
  integer :: ncid_ocean, ncid_erod, ncid_erod_mod
  integer :: ic, it, il, ij, ik, inb, jnb

  integer :: count2d(2), count4d(4)
  integer :: start2d(2), start4d(4)
  real*8, allocatable :: temp2d(:,:)
  real*8, allocatable :: temp4d(:,:,:,:)
!
! !REVISION HISTORY:
!   February2006, Jules Kouatchou (Jules.Kouatchou.1@gsfc.nasa.gov)
!     Original code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

   call Allocate_srcEmissDust   (i1, i2, ju1, j2)
   call Allocate_srcEmissSeaSalt(i1, i2, ju1, j2)
   call Allocate_erod_mod       (i1, i2, ju1, j2)
   call Allocate_water          (i1, i2, ju1, j2)
   call Allocate_erod           (i1, i2, ju1, j2)

   ilen = i2 - i1  + 1
   jlen = j2 - ju1 + 1

   allocate(temp2d(ilen,jlen))
   start2d(:) = (/ i1 - i1_gl + 1, ju1 - ju1_gl + 1 /)
   count2d (:) = (/ ilen, jlen /)

   call Ncop_Rd (ncid_erod_mod, GOCARTerod_mod_infile_name)
   call Ncrd_2d (temp2d, ncid_erod_mod, 'EROD_MOD', start2d, count2d )
   call Nccl (ncid_erod_mod)

   do ij = ju1, j2
      jnb = ij - ju1 + 1
      do il = i1, i2
         inb = il - i1 + 1
         erod_mod(il,ij) = temp2d(inb,jnb)
      end do
   end do

   call Ncop_Rd (ncid_ocean, GOCARTocean_infile_name)
   call Ncrd_2d (temp2d, ncid_ocean, 'WATER', start2d, count2d )
   call Nccl (ncid_ocean)

   do ij = ju1, j2
      jnb = ij - ju1 + 1
      do il = i1, i2
         inb = il - i1 + 1
         water(il,ij) = temp2d(inb,jnb)
      end do
   end do

   allocate(temp4d(ilen,jlen,ndcls, ndsrc))
   start4d(:) = (/ i1 - i1_gl + 1, ju1 - ju1_gl + 1, 1, 1 /)
   count4d (:) = (/ ilen, jlen, ndcls, ndsrc /)

   call Ncop_Rd (ncid_erod, GOCARTerod_infile_name)
   call Ncrd_4d (temp4d, ncid_erod, 'EROD', start4d, count4d )
   call Nccl (ncid_erod)

   do ij = ju1, j2
      jnb = ij - ju1 + 1
      do il = i1, i2
         inb = il - i1 + 1
         erod(il,ij,1:ndcls, 1:ndsrc) = temp4d(inb,jnb, 1:ndcls, 1:ndsrc)
      end do
   end do

   
   call Update_water(mcor, i1, i2, ju1, j2)
   call Update_erod (latdeg, londeg, mcor, i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl)

   deallocate(temp2d)
   deallocate(temp4d)

   return

   end subroutine GmiReadGocartSourceFiles
!EOC
!-------------------------------------------------------------------------
       end module GmiReadGocartSourceFiles_mod
