!------------------------------------------------------------------------------
! NASA GSFC - SIVO 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ReadUValbedoData_mod
!
      module ReadUValbedoData_mod
!
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: readUValbedoData
!
! !AUTHOR:
! Jules Kouatchou, Jules.Kouatchou-1@nasa.gov
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readUValbedoData
!
! !INTERFACE:
!
      subroutine readUValbedoData (uvalbedo_data, uvalbedo_infile_name,      &
     &             i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl, pr_diag, procID)
!
! !USES:
      use GmiASCIIoperations_mod, only : AsciiOpenRead
      use GmiPrintError_mod, only : GmiPrintError

      implicit none

#     include "GmiParameters.h"
#     include "gmi_time_constants.h"
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl
      character(len=*), intent(in) :: uvalbedo_infile_name
!
! !OUTPUT PARAMETERS:
      real*8, intent(out) :: uvalbedo_data(i1:i2, ju1:j2, MONTHS_PER_YEAR)
!
! !DESCRIPTION:
! Reads in the monthly UV albedo data.
! The each processor reads first the global data and then extracts its local data.
!
! !DEFINED PARAMETERS:
      integer, parameter :: NUM_PER_FULLDATA_LINE = 7
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      character (len=80) :: dumline
      integer :: fulldata_lines_per_im
      integer :: il, ij
      integer :: iline
      integer :: im, ix
      integer :: inum_left, inum_read
      integer :: lun, ilong, ilat
      integer :: num_extra
      real*8  :: uvalbedo(NUM_PER_FULLDATA_LINE)
      real*8, allocatable :: loc_uvalbedo_data(:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*)'readUValbedoData called by ', procID

      ilong = i2_gl - i1_gl  + 1
      ilat  = j2_gl - ju1_gl + 1

      allocate(loc_uvalbedo_data(ilong, ilat, MONTHS_PER_YEAR))
      loc_uvalbedo_data(:,:,:) = 0.0d0

      call AsciiOpenRead (lun, uvalbedo_infile_name)

      fulldata_lines_per_im = ((ilong * ilat) / NUM_PER_FULLDATA_LINE)

      num_extra = (ilong * ilat) - (NUM_PER_FULLDATA_LINE * fulldata_lines_per_im)

      do im = 1, MONTHS_PER_YEAR

         Read (lun, '(a80)') dumline

         iline     = 0
         inum_read = 0
         inum_left = 0

         do ij = 1, ilat
            do il = 1, ilong
               if (inum_left == 0) then
                  uvalbedo(:) = -1.0d0

                  iline = iline + 1

                  if (iline <= fulldata_lines_per_im) then
                     inum_read = NUM_PER_FULLDATA_LINE
                   else if (num_extra /= 0) then
                      inum_read = num_extra
                   else
                      err_msg = 'Problem in Read_Uvalbedo_Data.'
                      call GmiPrintError (err_msg, .true., 2, iline, num_extra,  &
     &                           0, 0.0d0, 0.0d0)
                   end if

                   Read (lun,*) (uvalbedo(ix), ix=1,inum_read)

                   inum_left = inum_read
                end if

                loc_uvalbedo_data(il,ij,im) = uvalbedo(inum_read-inum_left+1)

                inum_left = inum_left - 1
             end do
          end do
      end do

      uvalbedo_data(i1:i2, ju1:j2, :)  = loc_uvalbedo_data(i1:i2, ju1:j2, :)

      Close (lun)

      return

      end subroutine readUValbedoData
!EOC
!-----------------------------------------------------------------------------
      end module ReadUValbedoData_mod
