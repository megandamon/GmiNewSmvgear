!-------------------------------------------------------------------------
!BOP
!
! !MODULE: ReadForcedBC_mod
!
! !INTERFACE:
!
    module ReadForcedBC_mod
!
    implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: readForcedBcData
!
! !DESCRIPTION:
! Routine to set or read forced boundary condition data.
!
! !AUTHOR:
! John Tannahill (jrt@llnl.gov) and Jules Kouatchou (kouatchou@gsfc.nasa.gov)
!
! !REVISION HISTORY:
!
!EOP
!-------------------------------------------------------------------------

  CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readForcedBcData
!
! !INTERFACE:
!
      subroutine readForcedBcData (pr_diag, loc_proc, forc_bc_opt, &
                            forc_bc_years, forc_bc_num, forc_bc_init_val, &
                            forc_bc_infile_name, forc_bc_data)
!
! !USES:
!
      use GmiASCIIoperations_mod, only : AsciiOpenRead
!
      implicit none
!
#     include "GmiParameters.h"
#     include "gmi_forc_bc.h"
!
! !INPUT PARAMETERS:
      logical          , intent(in) :: pr_diag
      integer          , intent(in) :: loc_proc
      integer          , intent(in) :: forc_bc_opt, forc_bc_years, forc_bc_num
      real*8           , intent(in) :: forc_bc_init_val
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: forc_bc_infile_name
! !OUTPUT PARAMETERS:
      REAL*8  , intent(inOut) :: forc_bc_data(1:FBC_LATDIM,1:FBC_MONDIM,1:forc_bc_years,1:forc_bc_num)
!
! !DESCRIPTION:
! This routine sets/reads the forcing boundary condition values.
!
! !LOCAL VARIABLES:
      integer :: ic, ij, im, iy
      integer :: lun
!
!-----------------------------------------------------------------------------
!BOC

      if (pr_diag) then
        Write (6,*) 'readForcedBcData called by ', loc_proc
      end if


!     =====================
      if (forc_bc_opt == 1) then
!     =====================

          forc_bc_data(1:FBC_LATDIM,1:FBC_MONDIM,1:forc_bc_years,1:forc_bc_num) =  &
     &         forc_bc_init_val

!     ==========================
      else if (forc_bc_opt == 2) then
!     ==========================
         call AsciiOpenRead (lun, forc_bc_infile_name)

         do ic = 1, forc_bc_num
            do iy = 1, forc_bc_years
               do im = 1, FBC_MONDIM
                  Read (lun, 900) (forc_bc_data(ij,im,iy,ic), ij = 1, FBC_LATDIM)
               end do
            end do
         end do

 900     format (9x, 9e12.2/2x, 10e12.2)

         Close (lun)

!     ======
      end if
!     ======

      return

      end subroutine readForcedBcData
!EOC
!---------------------------------------------------------------------------

    end module ReadForcedBC_mod
