!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: ReadLossFrequency_mod
!
! !INTERFACE:
!
  module ReadLossFrequency_mod
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  readLossFrequency
!
! !DESCRIPTION:
!  Contains routine to set the loss frequency.
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


!-----------------------------------------------------------------------------
!
! ROUTINE
!   readLossFrequency
!
! ARGUMENTS
!   loss_freq : loss_rates (s^-1)
!
!-----------------------------------------------------------------------------

      subroutine readLossFrequency (loss_freq, lossData, loss_data_infile_name, &
                 loss_init_val, loss_freq_opt, loc_proc, pr_diag, &
                 i1, i2, ju1, j2, k1, k2, numSpecies, KDIM, JDIM, MDIM)

      implicit none
!
! !INPUT PARAMETERS:
!
      integer          , intent(in)   :: i1, i2, ju1, j2, k1, k2, numSpecies
      integer          , intent(in)   :: KDIM, JDIM, MDIM
      integer          , intent(in)   :: loss_freq_opt, loc_proc
      real*8           , intent(in)   :: loss_init_val
      logical          , intent(in)   :: pr_diag
      character (len=*), intent(in)   :: loss_data_infile_name
!
! !OUTPUT PARAMETERS:
!
      real*8 , intent(out)  :: loss_freq(i1:i2, ju1:j2, k1:k2, numSpecies)
      real*8 , intent(out)  :: lossData(KDIM, JDIM, MDIM, 1)
!
! !DESCRIPTION:
!  This routine either sets the loss frequency or reads in some loss
!  data that will be used to set the loss frequency later.
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

      if (pr_diag) then
        Write (6,*) 'ReadLossFrequency called by ', loc_proc
      end if


!     =======================
      if (loss_freq_opt == 1) then
!     =======================

        loss_freq(i1:i2,ju1:j2,:,:) = 0.0d0
        loss_freq(i1:i2,ju1:j2,:,1) = loss_init_val

!     ============================
      else if (loss_freq_opt == 2) then
!     ============================

        call ReadLossData (lossData, loss_data_infile_name, KDIM, JDIM, MDIM)

!     ======
      end if
!     ======

      return

      end subroutine readLossFrequency
!EOC
!-----------------------------------------------------------------------------
!
! ROUTINE
!   readLossData
!
! DESCRIPTION
!   This routine reads in initial loss data.
!
! ARGUMENTS
!   None
!
!-----------------------------------------------------------------------------

      subroutine readLossData (lossData,loss_data_infile_name, KDIM, JDIM, MDIM )

      use GmiASCIIoperations_mod, only : AsciiOpenRead

      implicit none

      integer          , intent(in)   :: KDIM, JDIM, MDIM
      character (len=*), intent(in)   :: loss_data_infile_name

      real*8           , intent(out)  :: lossData(KDIM, JDIM, MDIM, 1)

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ij, ik, im
      integer :: lun


!     ----------------
!     Begin execution.
!     ----------------

      call AsciiOpenRead (lun, loss_data_infile_name)

      do im = 1, MDIM
        do ij = 1, JDIM

          Read (lun, 900) (lossData(ik,ij,im,1), ik = KDIM, 1, -1)

        end do
      end do

 900  format (20x, 6e10.3/(8e10.3))

      Close (lun)

      return

      end subroutine readLossData

  end module ReadLossFrequency_mod
