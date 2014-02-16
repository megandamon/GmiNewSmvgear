!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: ReadWaterClimatology_mod
!
! !INTERFACE:
!
  module ReadWaterClimatology_mod
!
    use m_netcdf_io_open , only : Ncop_Rd
    use m_netcdf_io_read , only : Ncrd_4d
    use m_netcdf_io_close, only : Nccl
    use GmiCheckRange_mod, only : checkRange4d

      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  readWaterClimatology, readH2ocond

#     include "GmiParameters.h"
!
! !DESCRIPTION:
!  Contains routine to set the water climatology.
!
! !AUTHOR:
!   John Tannahill, LLNL, jrt@llnl.gov
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readWaterClimatology
!
! !INTERFACE:
!
      subroutine readWaterClimatology (h2oclim, ch4clim, &
                 h2oclim_init_val, ch4clim_init_val, h2oclim_opt, &
                 h2oclim_infile_name, loc_proc, pr_diag, &
                 i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl, h2oclim_timpyr)
!
! !USES:
!
      implicit none
!
! !INPUT PARAMETERS:
!   
      integer          , intent(in)   :: i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl
      integer          , intent(in)   :: h2oclim_timpyr
      integer          , intent(in)   :: h2oclim_opt, loc_proc
      real*8           , intent(in)   :: h2oclim_init_val, ch4clim_init_val
      logical          , intent(in)   :: pr_diag
      character (len=*), intent(in)   :: h2oclim_infile_name
!
! !OUTPUT PARAMETERS:
!
      real*8 , intent(out)  :: h2oclim(i1:i2, ju1:j2, k1:k2, h2oclim_timpyr)
!                            ! array of h2o climatology
      real*8 , intent(out)  :: ch4clim(i1:i2, ju1:j2, k1:k2, h2oclim_timpyr)
!                            ! array of ch4 climatology
!
! !DESCRIPTION:
!   This routine sets the water climatology.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_FILE_NAME) :: err_msg
      integer             :: ierr
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
        Write (6,*) 'ReadWaterClimatology called by ', loc_proc
      end if

!     =====================
      if (h2oclim_opt == 1) then
!     =====================

        h2oclim(i1:i2,ju1:j2,:,:) = h2oclim_init_val
        ch4clim(i1:i2,ju1:j2,:,:) = ch4clim_init_val

!     ==========================
      else if (h2oclim_opt == 2) then
!     ==========================

        call ReadH2oclim (h2oclim, ch4clim, h2oclim_infile_name, &
                 i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl, h2oclim_timpyr)

!     ======
      end if
!     ======

      !call Check_Range_4d  &
      call checkRange4d  &
     &  ('h2oclim', loc_proc, i1, i2, ju1, j2, k1, k2, 1,  &
     &   h2oclim_timpyr, h2oclim, 0.0d0, 1.0d-01)

      !call Check_Range_4d  &
      call checkRange4d  &
     &  ('ch4clim', loc_proc, i1, i2, ju1, j2, k1, k2, 1,  &
     &   h2oclim_timpyr, ch4clim, 0.0d0, 1.0d-04)

      return

      end subroutine ReadWaterClimatology

!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadH2oclim
!
! !INTERFACE:
!
      subroutine ReadH2oclim  &
                 (h2oclim, ch4clim, h2oclim_infile_name, &
                  i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl, h2oclim_timpyr)
!
      implicit none
!
! !INPUT PARAMETERS:
!
      integer          , intent(in)   :: i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl
      integer          , intent(in)   :: h2oclim_timpyr
      character (len=*), intent(in)   :: h2oclim_infile_name
!
! !OUTPUT PARAMETERS:
!
      real*8 , intent(out)  :: h2oclim(i1:i2, ju1:j2, k1:k2, h2oclim_timpyr)
!                            ! array of h2o climatology
      real*8 , intent(out)  :: ch4clim(i1:i2, ju1:j2, k1:k2, h2oclim_timpyr)
!                            ! array of ch4 climatology

!     -----------------------
!     Parameter declarations.
!     -----------------------

      character (len=MAX_LENGTH_VAR_NAME), parameter :: H2OCLIM_VNAM = 'h2o_clim'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: CH4CLIM_VNAM = 'ch4_clim'

!
! !DESCRIPTION:
!  This routine reads in a set of water climatology data and an associated
!  set of methane (ch4) climatology. Recall that all NetCDF arrays always 
!  start at subscript 1, and this is not always the case for GMI dimensions.
!
! !LOCAL VARIABLES:
!
      integer :: il, ij, ik, it,inb, jnb, ncid_hc
      integer :: cnt4d (4), strt4d(4)

      real*8  :: h2oclim1(i2-i1+1, j2-ju1+1, k2-k1+1, 1)
      real*8  :: ch4clim1(i2-i1+1, j2-ju1+1, k2-k1+1, 1)
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

      call Ncop_Rd (ncid_hc, h2oclim_infile_name)

      do it = 1, h2oclim_timpyr

        strt4d(4) = it

        call Ncrd_4d (h2oclim1, ncid_hc, H2OCLIM_VNAM, strt4d, cnt4d)
        call Ncrd_4d (ch4clim1, ncid_hc, CH4CLIM_VNAM, strt4d, cnt4d)

        do ik = k1, k2
          do ij = ju1, j2
            jnb = ij - ju1 + 1
            do il = i1, i2
              inb = il - i1 + 1

              h2oclim(il,ij,ik,it) = h2oclim1(inb,jnb,ik,1)
              ch4clim(il,ij,ik,it) = ch4clim1(inb,jnb,ik,1)

            end do
          end do
        end do

      end do

      call Nccl (ncid_hc)

      return

      end subroutine ReadH2oclim
!EOC
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadH2ocond
!
! !INTERFACE:
!
      subroutine ReadH2ocond  (h2ocond, restart_inrec, restart_infile_name, &
     &                i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl)
!
      implicit none
!
! !INPUT PARAMETERS:
!
      integer          , intent(in) :: i1, i2, ju1, j2, k1, k2
      integer          , intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer          , intent(in) :: restart_inrec
      character (len=*), intent(in) :: restart_infile_name
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out)  :: h2ocond(i1:i2,ju1:j2,k1:k2)
!
! !DEFINED PARAMETERS:
      character (len=MAX_LENGTH_VAR_NAME), parameter :: H2OCOND_VNAM = 'h2ocond'
!
! !DESCRIPTION:
!  Reads the h2ocond array from the restart file.
!
! !LOCAL VARIABLES:
!
      integer :: il, ij, ik, inc, jnc, ncid
      integer :: ilong, ilat, ivert
      integer :: cnt4d (4), strt4d(4)
      real*8, allocatable :: h2ocond4d_nc(:,:,:,:)
      real*8, allocatable :: h2ocondGlob (:,:,:)
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
      ilong = i2_gl - i1_gl  + 1
      ilat  = j2_gl - ju1_gl + 1
      ivert = k2    - k1     + 1

      allocate(h2ocond4d_nc(ilong, ilat, ivert, 1))
      allocate(h2ocondGlob (i1_gl:i2_gl, ju1_gl:j2_gl, k1:k2))

      call Ncop_Rd (ncid, restart_infile_name)

      strt4d(:) = (/ 1, 1, 1, restart_inrec /)
      cnt4d (:) = (/ ilong, ilat, ivert, 1 /)

      call Ncrd_4d (h2ocond4d_nc, ncid, H2OCOND_VNAM, strt4d, cnt4d)

      do ik = k1, k2
         do ij = ju1_gl, j2_gl
            jnc = ij - ju1_gl + 1
            do il = i1_gl, i2_gl
               inc = il - i1_gl + 1
               h2ocondGlob(il,ij,ik) = h2ocond4d_nc(inc,jnc,ik,1)
            end do
         end do
      end do

      call Nccl (ncid)

      h2ocond(i1:i2,ju1:j2,k1:k2) = h2ocondGlob(i1:i2,ju1:j2,k1:k2)

      deallocate(h2ocond4d_nc)
      deallocate(h2ocondGlob)

      return

      end subroutine readH2ocond
!EOC
!------------------------------------------------------------------------

  end module ReadWaterClimatology_mod
