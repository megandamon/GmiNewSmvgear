!-----------------------------------------------------------------------
!                     NASA/GSFC - SIVO ASTG Code 610.3                 !
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m_FastjMIEparameters
!
! !INTERFACE:
      module m_FastjMIEparameters
!
! !PUBLIC DATA MEMBERS:
!
      private
      public  :: N_
      public  :: M_ 

      integer, parameter ::   N_ = 501
                                 ! no. of levels in Mie scattering arrays
                                 ! = 2*NC+1 = 4*LPAR + 5 + 2*sum(JADDLV)
      integer, parameter ::   M_ =   4
                                 ! no. of Gauss points used, must 
                                 ! = 4 in fast_JX (no option)
!
! !DESCRIPTION:
! 
!
! !AUTHOR:
! Michael Prather
!
! !HISTORY:
!
!EOP
!-----------------------------------------------------------------------

      end module m_FastjMIEparameters
