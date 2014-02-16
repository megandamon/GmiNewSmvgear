!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: ReadAerosolDust_mod
!
! !INTERFACE:
!
       module ReadAerosolDust_mod

       use m_netcdf_io_open , only : Ncop_Rd
       use m_netcdf_io_close, only : Nccl
       use m_netcdf_io_read , only : Ncrd_5d
!
       implicit none
!
       private
       public   :: ReadAerosolDust
!
! !DESCRIPTION:
!  Routine to read a specified record of a netCDF file containing 
!  aerosol and mineral dust data. The reading is done by all the
!  worker processors.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  12Jan2007 Initial code.
!
!EOP
!-------------------------------------------------------------------------

       CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ReadAerosolDust
!
! !INTERFACE:
!
      subroutine ReadAerosolDust (dust, wAersl, dAersl, AerDust_infile_name, &
                  curRecord, i1,i2, ju1, j2, k1, k2, i1_gl, ju1_gl)
!
! !USES:
!
      implicit none

#     include "GmiParameters.h"
#     include "setkin_par.h"
!
! !INPUT PARAMETERS:
      integer            , intent(in) :: i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl
      integer            , intent(in) :: curRecord
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: AerDust_infile_name
      real*8             , intent(out):: dust  (i1:i2, ju1:j2, k1:k2, nSADdust)
      real*8             , intent(out):: wAersl(i1:i2, ju1:j2, k1:k2, nSADaer )
      real*8             , intent(out):: dAersl(i1:i2, ju1:j2, k1:k2, 2       )
!
! !DESCRIPTION: This routine reads in a record (for a given day or month) the 
!               netCDF aerosol/dust file (monthly or daily) and populates
!               the corresponding aerosol/dust variables.
!               This routine is called by the worker processors only.
!               Each worker processor reads its portion of the data.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      character (len=12) :: aero_var_name='aerofield'
      integer :: ij, il, ii, j_len, i_len, k_len, ik, jnb, inb
      integer :: ncid_aero, nferr, var_id_aero
      integer :: cnt5d (5)
      integer :: strt5d(5)
      real*8, allocatable :: tmp_aerofield(:,:,:,:)
!
! !REVISION HISTORY:
!   February2005, Jules Kouatchou (Jules.Kouatchou.1@gsfc.nasa.gov)
!     Original code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifndef nonZeroInd_tracers
      i_len = i2 - i1  + 1
      j_len = j2 - ju1 + 1
      k_len = k2 - k1  + 1

      strt5d(1) = i1  -  i1_gl + 1
      strt5d(2) = ju1 - ju1_gl + 1
      strt5d(3) = k1
      strt5d(4) = 1
      strt5d(5) = curRecord

      cnt5d(:) = (/ i_len, j_len, k_len, 16, 1 /)

      allocate(tmp_aerofield(i_len, j_len, k_len, 16))
      tmp_aerofield = 0.0d0

      ! Open the netCDF file
      call Ncop_Rd (ncid_aero, AerDust_infile_name)

      ! Read the data
      call Ncrd_5d (tmp_aerofield, ncid_aero, 'aerofield', strt5d, cnt5d)

      ! Close the file
      call Nccl (ncid_aero)

      !-------------------------------------------------------------------
      ! Load tmp_aerofield (i.e., aerosol/dust data that was just read in)
      ! into proper location in aerosol/dust arrays.
      !-------------------------------------------------------------------

      do ik = k1, k2
         do ij = ju1, j2
            jnb = ij - ju1 + 1
            do il = i1, i2
               inb = il - i1 + 1

               ! Dust
               do ii=1,NSADdust
                  dust(il,ij,ik,ii) = tmp_aerofield(inb,jnb,ik,ii)
               end do

               ! Sulfate SO4
               wAersl(il,ij,ik,1) =   tmp_aerofield(inb,jnb,ik,16)

               ! Hydrophobic BC
               dAersl(il,ij,ik,1) =   tmp_aerofield(inb,jnb,ik, 8)

               ! Hydrophilic BC
               wAersl(il,ij,ik,2) =   tmp_aerofield(inb,jnb,ik, 9)

               ! Hydrophobic OC
               dAersl(il,ij,ik,2) =   tmp_aerofield(inb,jnb,ik,10)
         
               ! Hydrophilic OC
               wAersl(il,ij,ik,3) =   tmp_aerofield(inb,jnb,ik,11)

               ! Sea Salt (accum)
               wAersl(il,ij,ik,4) =   tmp_aerofield(inb,jnb,ik,12)
!.sds... expanded number of sea salt bins carried to 7
                ! Sea Salt (coarse1)
               wAersl(il,ij,ik,5) =   tmp_aerofield(inb,jnb,ik,13)
               ! Sea Salt (coarse2)
               wAersl(il,ij,ik,6) =   tmp_aerofield(inb,jnb,ik,14)
               ! Sea Salt (coarse3)
               wAersl(il,ij,ik,7) =   tmp_aerofield(inb,jnb,ik,15)

            enddo
         enddo
      enddo

#endif
      return

      end subroutine ReadAerosolDust
!EOC
!-------------------------------------------------------------------------

   end module ReadAerosolDust_mod
