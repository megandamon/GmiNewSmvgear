!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: ReadLightningData_mod
!
! !INTERFACE:
!
     module ReadLightningData_mod
!
! !USES:
     use m_netcdf_io_open  , only : Ncop_Rd
     use m_netcdf_io_close , only : Nccl
     use m_netcdf_io_read  , only : Ncrd_1d, Ncrd_2d, Ncrd_3d
     use GmiTimeControl_mod, only : GmiSplitDateTime, GetDaysFromJanuary1
!
     implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
     private
     public  :: readLightningCoeffData

#     include "GmiParameters.h"
#     include "gmi_time_constants.h"

! !DESCRIPTION:
! Reads and stores local and global ratios for lightning routine
!
! !AUTHOR:
! Megan Damon
! Dale Allen

! !REVISION HISTORY:
!  Initial code.
!
!EOP
!------------------------------------------------------------------------------
     contains
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readLightningCoeffData 
!
! !INTERFACE:
!
      subroutine readLightningCoeffData &
     &           (lightCoeff_infile_name, globalCoeff, localCoeff, &
     &            midLatAdj, lightThreshold, ik0, pr_diag, procID, curMonth, &
                  curYear, i1, i2, ju1, j2, i1_gl, ju1_gl)

      implicit none

! !INPUT PARAMETERS:
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: lightCoeff_infile_name
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID, curMonth, curYear
      integer, intent(in) :: i1, i2, ju1, j2
      integer, intent(in) :: i1_gl, ju1_gl

!
! !OUTPUT PARAMETERS:
                              ! global ratios for the current month
      real*8 , intent(out) :: globalCoeff
                              ! local ratios for the current month
      real*8 , intent(out) :: localCoeff(i1:i2, ju1:j2)
                              ! mid latitude adjustment
      real*8 , intent(out) :: midLatAdj(i1:i2, ju1:j2)
                              ! threshold for lightning
      real*8 , intent(out) :: lightThreshold
                              ! desired layer for 3-d Convective mass flux
      integer , intent(out) :: ik0
      

! !DESCRIPTION: 
!
! !LOCAL VARIABLES:
      integer :: count4d(4), start4d(4), count1d (2)
      integer :: count2d(2), start2d(2)
      integer :: ncid
      integer :: ilength, jlength, il, ij, inb, jnb
      character (len=MAX_LENGTH_VAR_NAME), parameter ::          &
                             localCoeffVarname = 'ratio_local',  &
                              midLatAdjVarName = 'midlat',       &
                            globalCoeffVarname = 'ratio_global', &
                              thresholdVarname = 'threshold',    &
                                    ik0Varname = 'ik0'
      
      real*8 :: tempVar(i1:i2, ju1:j2)
      real*8 :: tempScalar(1)

!
! !AUTHOR:
!  
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if (pr_diag) then
         Write (6, *) 'readLightningCoeffData called by ', procID
         Write (6, *), 'Ready to read month ', curMonth, ' and year ', curYear
      end if

      IF (procID == 0) THEN
      Write (6, *), 'Ready to read month ', curMonth, ' and year ', curYear
      END IF
      start4d(1) = i1  -  i1_gl + 1
      start4d(2) = ju1 - ju1_gl + 1
      start4d(3) = curMonth
      start4d(4) = curYear


      ilength = i2 - i1 + 1
      jlength = j2 - ju1 + 1

      count4d(:) = (/ ilength, jlength, 1, 1 /)

      start2d(1) = start4d(1)
      start2d(2) = start4d(2)
      count2d(:) = (/ ilength, jlength /)

      count1d(1) = 1

      call Ncop_Rd (ncid, lightCoeff_infile_name)

      ! Reading local data based on the current month data
      call Ncrd_3d (localCoeff, ncid, localCoeffVarName, start4d, count4d)

      ! Reading mid latitude adjustment
      call Ncrd_2d (midLatAdj, ncid, midLatAdjVarName, start2d, count2d)
      
     ! Reading global coefficient for lightning
     ! redefine start2d and count2d
      start2d(1) = curMonth
      start2d(2) = curYear
      count2d(1) = 1
      count2d(2) = 1
      call Ncrd_1d (tempScalar, ncid, globalCoeffVarName, start2d, count2d)

      globalCoeff = tempScalar(1)
      
      ! Reading threshold
      call Ncrd_1d (tempScalar, ncid, thresholdVarName, (/1/), count1d)

      lightThreshold = tempScalar(1)

      ! Reading desired layer for 3-D Convective mass flux
      call Ncrd_1d (tempScalar, ncid, ik0Varname, (/1/), count1d)

      ik0 = int (tempScalar(1))

      call Nccl (ncid)

      return

      end subroutine readLightningCoeffData
!EOC

      end module ReadLightningData_mod
