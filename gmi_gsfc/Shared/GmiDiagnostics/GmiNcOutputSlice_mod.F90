!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiNcOutputSlice_mod
!
! !INTERFACE:
!
module GmiNcOutputSlice_mod
!
! !USES:
  use m_netcdf_io_write , only : Ncwr_5d, Ncwr_6d
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: prepSlice4d
  public  :: prepSlice4d_Nomean
  public  :: prepSlice4d_Map
  public  :: prepSlice4d_Map_Nomean
  public  :: prepSlice4d_Emiss
  public  :: writeSlice4d
  public  :: writeSlice5d
!
! !DESCRIPTION:
!  Provides routines for writing out slices of arrays.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------------
!
! ROUTINE
!   prepSlice4d
!
! DESCRIPTION
!   This routine prepares one 3D slice of a 4D array for NetCDF output.
!
! ARGUMENTS
!   do_mean    : should means or current values be put in periodic output
!                files?
!   k1x, k2x   : altitude  dimensions for norm_array, mean_array, & slice_nc
!   snum_out   : 4th dimension slice index
!   dim4       : 4th dimension of norm_array & mean_array
!   norm_array : normal 4D array
!   mean_array : mean   4D array
!   slice_nc   : slice of 3D data extracted
!
!-----------------------------------------------------------------------------

      subroutine prepSlice4d  &
     &  (do_mean, i1x, i2x, j1x, j2x, k1x, k2x, snum_out, dim4, norm_array, mean_array,  &
     &   slice_nc)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------
    
      integer, intent(in) :: i1x, i2x, j1x, j2x

      logical :: do_mean
      integer :: k1x, k2x
      integer :: snum_out
      integer :: dim4
      real*8  :: norm_array(i1x:i2x, j1x:j2x, k1x:k2x, dim4)
      real*8  :: mean_array(i1x:i2x, j1x:j2x, k1x:k2x, dim4)
      real*8  :: slice_nc  (i1x:i2x, j1x:j2x, k1x:k2x)
!EOP
!------------------------------------------------------------------------------
!BOC
      if (do_mean) then
        slice_nc(:,:,:) = mean_array(:,:,:,snum_out)
      else
        slice_nc(:,:,:) = norm_array(:,:,:,snum_out)
      end if

      return

      end subroutine prepSlice4d
!EOC
!------------------------------------------------------------------------------
!BOP
!
! ROUTINE
!   prepSlice4d_Nomean
!
! DESCRIPTION
!   This routine prepares one 3D slice of a snapshot of a 4D array for NetCDF
!   output (no mean values).
!
! ARGUMENTS
!   k1x, k2x   : altitude  dimensions for norm_array & slice_nc
!   snum_out   : 4th dimension slice index
!   dim4       : 4th dimension of norm_array
!   norm_array : normal 4D array
!   slice_nc   : slice of 3D data extracted
!
!-----------------------------------------------------------------------------

      subroutine prepSlice4d_Nomean  &
     &  (i1x, i2x, j1x, j2x, k1x, k2x, snum_out, dim4, norm_array, slice_nc)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1x, i2x, j1x, j2x
      integer :: k1x, k2x
      integer :: snum_out
      integer :: dim4
      real*8  :: norm_array(i1x:i2x, j1x:j2x, k1x:k2x, dim4)
      real*8  :: slice_nc  (i1x:i2x, j1x:j2x, k1x:k2x)
!EOP
!------------------------------------------------------------------------------
!BOC

      slice_nc(i1x:i2x,j1x:j2x,:) =  &
     &  norm_array(i1x:i2x,j1x:j2x,:,snum_out)

      return

      end subroutine prepSlice4d_Nomean
!EOC
!------------------------------------------------------------------------------
!BOP
!
! ROUTINE
!   prepSlice4d_Map
!
! DESCRIPTION
!   This routine prepares one 3D slice of a 4D mapped array for NetCDF
!   output.
!
! ARGUMENTS
!   do_mean     : should means or current values be put in periodic output
!                 files?
!   k1x, k2x    : altitude  dimensions for norm_array, mean_array, & slice_nc
!   snum_out    : 4th dimension slice index
!   dim4        : 4th dimension of norm_array & mean_array
!   num_outrecs : number of species that will be output
!   outrec_map  : mapping of species output number to species input number
!   norm_array  : normal 4D array
!   mean_array  : mean   4D array
!   slice_nc    : slice of 3D data extracted
!
!-----------------------------------------------------------------------------

      subroutine prepSlice4d_Map  &
     &  (do_mean, i1x, i2x, j1x, j2x, k1x, k2x, snum_out, dim4, num_outrecs, outrec_map,  &
     &   norm_array, mean_array, slice_nc)

      implicit none


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1x, i2x, j1x, j2x
      logical :: do_mean
      integer :: k1x, k2x
      integer :: snum_out
      integer :: dim4
      integer :: num_outrecs
      integer :: outrec_map(num_outrecs)
      real*8  :: norm_array(i1x:i2x, j1x:j2x, k1x:k2x, dim4)
      real*8  :: mean_array(i1x:i2x, j1x:j2x, k1x:k2x, num_outrecs)
      real*8  :: slice_nc  (i1x:i2x, j1x:j2x, k1x:k2x)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ic
!EOP
!------------------------------------------------------------------------------
!BOC

      if (do_mean) then

!       ------------------------------------------------------------------
!       No "unmapping" needed here as any mean data was only collected for
!       the slices that will be output.
!       ------------------------------------------------------------------

        slice_nc(:,:,:) = mean_array(:,:,:,snum_out)

      else

        ic = outrec_map(snum_out)

        slice_nc(:,:,:) = norm_array(:,:,:,ic)

      end if

      return

      end subroutine prepSlice4d_Map
!EOC
!------------------------------------------------------------------------------
!BOP
!
! ROUTINE
!   prepSlice4d_Map_Nomean
!
! DESCRIPTION
!   This routine prepares one 3D slice of a snapshot of a 4D mapped array for
!   NetCDF output (no mean values).
!
! ARGUMENTS
!   k1x, k2x    : altitude  dimensions for norm_array & slice_nc
!   snum_out    : 4th dimension slice index
!   dim4        : 4th dimension of norm_array
!   num_outrecs : number of species that will be output
!   outrec_map  : mapping of species output number to species input number
!   norm_array  : normal 4D array
!   slice_nc    : slice of 3D data extracted
!
!-----------------------------------------------------------------------------

      subroutine prepSlice4d_Map_Nomean  &
     &  (i1x, i2x, j1x, j2x, k1x, k2x, snum_out, dim4, num_outrecs, &
     &   outrec_map, norm_array, slice_nc)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1x, i2x, j1x, j2x
      integer :: k1x, k2x
      integer :: snum_out
      integer :: dim4
      integer :: num_outrecs
      integer :: outrec_map(num_outrecs)
      real*8  :: norm_array(i1x:i2x, j1x:j2x, k1x:k2x, dim4)
      real*8  :: slice_nc  (i1x:i2x, j1x:j2x, k1x:k2x)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ic
!EOP
!------------------------------------------------------------------------------
!BOC

      ic = outrec_map(snum_out)
      slice_nc(:,:,:) = norm_array(:,:,:,ic)

      return

      end subroutine prepSlice4d_Map_Nomean
!EOC
!------------------------------------------------------------------------------
!BOP
!
! ROUTINE
!   prepSlice4d_Emiss
!
! DESCRIPTION
!   This routine prepares one 3D slice of a 4D mapped array for NetCDF
!   output.
!
! ARGUMENTS
!   do_mean     : should means or current values be put in periodic output
!                 files?
!   k1x, k2x    : altitude  dimensions for norm_array, mean_array, & slice_nc
!   snum_out    : 4th dimension slice index
!   dim4        : 4th dimension of norm_array & mean_array
!   num_outrecs : number of species that will be output
!   outrec_map  : mapping of species output number to species input number
!   norm_array  : normal 4D array
!   mean_array  : mean   4D array
!   slice_nc    : slice of 3D data extracted
!
!-----------------------------------------------------------------------------

      subroutine prepSlice4d_Emiss  &
     &  (i1x, i2x, j1x, j2x, k1x, k2x, snum_out, dim4, num_outrecs, outrec_map,  &
     &   norm_array, slice_nc)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1x, i2x, j1x, j2x
      integer :: k1x, k2x
      integer :: snum_out
      integer :: dim4
      integer :: num_outrecs
      integer :: outrec_map(num_outrecs)
      real*8  :: norm_array(i1x:i2x, j1x:j2x, k1x:k2x, dim4)
      real*8  :: slice_nc  (i1x:i2x, j1x:j2x, k1x:k2x)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ic
!EOP
!------------------------------------------------------------------------------
!BOC

!       ------------------------------------------------------------------
!       No "unmapping" needed here as any mean data was only collected for
!       the slices that will be output.
!       ------------------------------------------------------------------

      ic = outrec_map(snum_out)
      slice_nc(:,:,:) = norm_array(:,:,:,ic)

      return

      end subroutine prepSlice4d_Emiss
!EOC
!------------------------------------------------------------------------------
!BOP
!
! ROUTINE
!   writeSlice4d
!
! DESCRIPTION
!   This routine outputs one 3D slice of a 4D array in NetCDF format.
!
! ARGUMENTS
!   k1x, k2x : altitude  dimensions for slice_nc
!   slice_nc : 3D slice to output
!   ncid     : NetCDF file id
!   varname  : NetCDF variable name
!   snum_out : 4th dimension slice  index
!   rnum_out : 5th dimension record index (i.e., NetCDF record to write to)
!
!-----------------------------------------------------------------------------

      subroutine writeSlice4d  &
     &  (i1x, i2x, j1x, j2x, k1x, k2x, slice_nc, ncid, varname, snum_out, rnum_out)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1x, i2x, j1x, j2x
      integer :: k1x, k2x
      real*8  :: slice_nc(i1x:i2x, j1x:j2x, k1x:k2x)
      integer :: ncid
      character (len=*) :: varname
      integer :: snum_out
      integer :: rnum_out

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ilon, ilat, ivert

      integer :: cnt5d (5)
      integer :: strt5d(5)
!EOP
!------------------------------------------------------------------------------
!BOC

      ilon  = i2x - i1x + 1
      ilat  = j2x - j1x + 1
      ivert = k2x - k1x + 1

      strt5d(:) = (/ 1, 1, 1, snum_out, rnum_out /)

      cnt5d (:) = (/ ilon, ilat, ivert, 1, 1 /)

!     ============
      call Ncwr_5d  (slice_nc(i1x:i2x, j1x:j2x, k1x:k2x), ncid, varname, strt5d, cnt5d)
!     ============

      return

      end subroutine writeSlice4d
!EOC
!------------------------------------------------------------------------------
!BOP
!
! ROUTINE
!   writeSlice5d
!
! DESCRIPTION
!   This routine outputs one 3D slice of a 5D array in NetCDF format.
!
! ARGUMENTS
!   k1x, k2x : altitude  dimensions for slice_nc
!   slice_nc : 3D slice to output
!   ncid     : NetCDF file id
!   varname  : NetCDF variable name
!   snum_out : 4th dimension slice  index
!   snum_out5: 5th dimension slice  index
!   rnum_out : 6th dimension record index (i.e., NetCDF record to write to)
!
!-----------------------------------------------------------------------------

      subroutine writeSlice5d  &
     &  (i1x, i2x, j1x, j2x, k1x, k2x, slice_nc, ncid, varname, snum_out, snum_out5,  &
     &   rnum_out)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1x, i2x, j1x, j2x
      integer :: k1x, k2x
      real*8  :: slice_nc(i1x:i2x, j1x:j2x, k1x:k2x)
      integer :: ncid
      character (len=*) :: varname
      integer :: snum_out
      integer :: snum_out5
      integer :: rnum_out

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ilon, ilat, ivert

      integer :: cnt6d (6)
      integer :: strt6d(6)
!EOP
!------------------------------------------------------------------------------
!BOC

      ilon  = i2x - i1x + 1
      ilat  = j2x - j1x + 1
      ivert = k2x - k1x + 1

      strt6d(:) = (/ 1, 1, 1, snum_out, snum_out5, rnum_out /)

      cnt6d (:) = (/ ilon, ilat, ivert, 1, 1, 1 /)

!     ============
      call Ncwr_6d (slice_nc(i1x:i2x, j1x:j2x, k1x:k2x), ncid, varname, strt6d, cnt6d)
!     ============

      return

      end subroutine writeSlice5d
!EOC
!------------------------------------------------------------------------------
end module GmiNcOutputSlice_mod
