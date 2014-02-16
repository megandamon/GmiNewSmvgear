    module ReadVegLaiData_mod

    use GmiASCIIoperations_mod, only : AsciiOpenRead

    implicit none

!=============================================================================
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
!   Original code from:
!     Harvard tropospheric emissions module for 3D applications;
!       by Yuhang Wang, Gerry Gardner, and Prof. Daniel Jacob
!       of Harvard University (Release V1.0)
!
    private
    public  ::  readVegetationData
    public  ::  readLeafAreaIndexData

#     include "GmiParameters.h"
#     include "gmi_emiss_constants.h"
#     include "gmi_time_constants.h"


!=============================================================================

   CONTAINS

!-----------------------------------------------------------------------------
!BOC
!
! ROUTINE
!   readVegetationData

      subroutine readVegetationData  &
     &  (veg_infile_name, ireg, iland, iuse, pr_diag, loc_proc, i1, i2, ju1, j2)


      implicit none
!
! !INPUT PARAMETERS:
      logical          , intent(in) :: pr_diag
      integer          , intent(in) :: loc_proc
      integer          , intent(in) :: i1, i2, ju1, j2
      character (len=*), intent(in) :: veg_infile_name
!                                      vegetation type input file name
!
! !OUTPUT PARAMETERS:
      integer, intent(out) :: ireg (i1:i2, ju1:j2)
!                             number of land types in a grid square
      integer, intent(out) :: iland(i1:i2, ju1:j2, NTYPE)
!                             land type id in grid square for ireg land types
      integer, intent(out) :: iuse (i1:i2, ju1:j2, NTYPE)
!                             fraction of grid box area occupied by land type (mil^-1?)
!
! !DESCRIPTION:
!   This routine reads in the vegetation data.
!
! !LOCAL VARIABLES:
      integer :: idummy
      integer :: il, ij, it, iost
      integer :: lun
      integer :: ireg_temp 
      integer :: iland_temp(NTYPE)
      integer :: iuse_temp (NTYPE)
      integer :: i_grid, j_grid
!EOP
!-------------------------------------------------------------------------
!BOC
      if (pr_diag) then
        Write (6,*) 'readVegetationData called by ', loc_proc
      end if

      call AsciiOpenRead (lun, veg_infile_name)

      READ_LOOP: do

         Read (lun, 905, iostat = iost) &
     &        i_grid, j_grid, ireg_temp,  &
     &        (iland_temp(it), it=1,ireg_temp),  &
     &        (iuse_temp (it), it=1,ireg_temp)

         if (iost < 0) exit READ_LOOP

         if (((i_grid >= i1 ) .and. (i_grid <= i2)) .and.  &
     &       ((j_grid >= ju1) .and. (j_grid <= j2))) then

            ireg (i_grid,j_grid)             = ireg_temp
            iland(i_grid,j_grid,1:ireg_temp) = iland_temp(1:ireg_temp)
            iuse (i_grid,j_grid,1:ireg_temp) = iuse_temp (1:ireg_temp)

         end if

      end do READ_LOOP

 905  format (33i4)

      Close (lun)

      return

      end subroutine readVegetationData
!EOC
!-----------------------------------------------------------------------------
!BOP
      subroutine readLeafAreaIndexData (lai_infile_name, ireg, xlai, xlai2, &
     &           pr_diag, loc_proc, nymd, i1, i2, ju1, j2, lai_day_save, &
                 firstReadLaiData)

      use GmiTimeControl_mod, only : GmiSplitDateTime

      implicit none
!
! !INPUT PARAMETERS:
      logical          , intent(in) :: pr_diag
      integer          , intent(in) :: loc_proc
      integer          , intent(in) :: i1, i2, ju1, j2
      integer          , intent(in) :: nymd
      character (len=*), intent(in) :: lai_infile_name
      integer          , intent(in) :: ireg(i1:i2, ju1:j2)
!
! !INPUT/OUTPUT PARAMETERS:
      logical          , intent(inOut) :: firstReadLaiData
      integer          , intent(inOut) :: lai_day_save 
      real*8           , intent(inOut) :: xlai (i1:i2, ju1:j2, NTYPE)
      real*8           , intent(inOut) :: xlai2(i1:i2, ju1:j2, NTYPE)
!
! !DESCRIPTION:
! This routine reads in the leaf area index file.
!
! !LOCAL VARIABLES:
      integer       :: iday
      integer       :: mdummy, ydummy
!EOP
!-------------------------------------------------------------------------
!BOC
      if (pr_diag) then
        Write (6,*) 'readLeafAreaIndexData called by ', loc_proc
      end if

      call GmiSplitDateTime (nymd, ydummy, mdummy, iday)

      if (iday /= lai_day_save) then

        lai_day_save = iday

!       ==================
        call Read_Lai_Data  &
!       ==================
     &    (lai_infile_name, nymd, ireg, xlai, xlai2, &
     &     pr_diag, loc_proc, i1, i2, ju1, j2, firstReadLaiData)

      end if

      return

      end subroutine readLeafAreaIndexData
!EOC
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Read_Lai_Data
!
! DESCRIPTION
!   This routine reads in the leaf area index (lai) if this a new month.
!
!   All 12 months of data are in one file and they are separated by a line
!   that has a zero and then a month number.
!
! ARGUMENTS
!   lai_infile_name : leaf area index input file name
!   nymd  : year/month/day (YYYYMMDD)
!   ireg  : number of land types in a grid square
!   xlai  : leaf area index of land type for month #1
!   xlai2 : leaf area index of land type for month #2
!
!-----------------------------------------------------------------------------

      subroutine Read_Lai_Data  &
     &  (lai_infile_name, nymd, ireg, xlai, xlai2, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, firstReadLaiData)

       use GmiTimeControl_mod, only : GmiSplitDateTime
       use GmiTimeControl_mod, only : GetDaysFromJanuary1

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2
      character (len=MAX_LENGTH_FILE_NAME) :: lai_infile_name
      integer :: nymd
      integer :: ireg (i1:i2, ju1:j2)
      real*8  :: xlai (i1:i2, ju1:j2, NTYPE)
      real*8  :: xlai2(i1:i2, ju1:j2, NTYPE)
      logical, intent(inOut) :: firstReadLaiData


!     -----------------------
!     Parameter declarations.
!     -----------------------

!c?   leap years?
      integer, parameter :: STARTDAY(MONTHS_PER_YEAR+1) =  &
!c?  &  (/  15,  45,  74, 105, 135, 166,
     &  (/  15,  46,  74, 105, 135, 166,  &
     &     196, 227, 258, 288, 319, 349, 380 /)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ddummy, ydummy
      integer :: il, ij, ik
      integer :: imonth
      integer :: imul, itd
      integer :: jday
      integer :: month1, month2
      integer :: monpyr

      real*8  :: rimul, ritd


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Read_Lai_Data called by ', loc_proc
      end if


      monpyr = MONTHS_PER_YEAR


!     ====================
      call GmiSplitDateTime  &
!     ====================
     &  (nymd, ydummy, imonth, ddummy)

!     =================
      call GetDaysFromJanuary1  &
!     =================
     &  (jday, nymd)

      if (jday >= STARTDAY(imonth)) then
        month1 = imonth
      else
        month1 = imonth - 1
        if (month1 == 0) month1 = monpyr
      end if


!     ==========================================
      if (firstReadLaiData .or. (jday == STARTDAY(month1))) then
!     ==========================================

        month2 = month1 + 1
        if (month2 == (monpyr+1)) month2 = 1

!       ========================
        call Read_Lai_Two_Months  &
!       ========================
     &    (lai_infile_name, month1, month2, ireg, xlai, xlai2, &
     &     pr_diag, loc_proc, i1, i2, ju1, j2)

        itd  = STARTDAY(month1+1) - STARTDAY(month1)
        ritd = itd

        do ij = ju1, j2
          do il = i1, i2
            do ik = 1, ireg(il,ij)

              xlai2(il,ij,ik) =  &
     &          (xlai2(il,ij,ik) - xlai(il,ij,ik)) / ritd

            end do
          end do
        end do

!       ==========
        if (firstReadLaiData) then
!       ==========

          firstReadLaiData = .false.

          if (jday >= STARTDAY(1)) then
            imul = jday - STARTDAY(month1)
          else
            imul = DAYS_PER_YEAR - STARTDAY(monpyr) +  &
     &             jday
          end if

          rimul = imul

!c?       cache misses?
          do ij = ju1, j2
            do il = i1, i2
              do ik = 1, ireg(il,ij)

                xlai(il,ij,ik)  =  &
     &            xlai(il,ij,ik) + (xlai2(il,ij,ik) * rimul)

              end do
            end do
          end do

        end if

!     ====
      else
!     ====

!c?     cache misses?
        do ij = ju1, j2
          do il = i1, i2
            do ik = 1, ireg(il,ij)

              xlai(il,ij,ik) =  &
     &          xlai(il,ij,ik) + xlai2(il,ij,ik)

            end do
          end do
        end do

!     ======
      end if
!     ======


      return

      end subroutine Read_Lai_Data


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Read_Lai_Two_Months
!
! DESCRIPTION
!   This routine reads in the leaf area index (lai) for two months.
!
!   All 12 months of data are in one file and they are separated by a line
!   that has a zero and then a month number.
!
! ARGUMENTS
!   lai_infile_name : leaf area index input file name
!   month1 : month #1
!   month2 : month #2
!   ireg   : number of land types in a grid square
!   xlai   : leaf area index of land type for month #1
!   xlai2  : leaf area index of land type for month #2
!
!-----------------------------------------------------------------------------

      subroutine Read_Lai_Two_Months  &
     &  (lai_infile_name, month1, month2, ireg, xlai, xlai2, &
         pr_diag, loc_proc, i1, i2, ju1, j2)

      use GmiPrintError_mod, only : GmiPrintError
      use GmiASCIIoperations_mod, only : AsciiOpenRead

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2
      character (len=MAX_LENGTH_FILE_NAME) :: lai_infile_name
      integer :: month1
      integer :: month2
      integer :: ireg (i1:i2, ju1:j2)
      real*8  :: xlai (i1:i2, ju1:j2, NTYPE)
      real*8  :: xlai2(i1:i2, ju1:j2, NTYPE)


!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      logical :: found_month1 = .false.
      logical :: found_month2 = .false.

      integer :: il, ij, ik, it
      integer :: index_type
      integer :: iost
      integer :: lun
      integer :: month_found

      real*8  :: xlai_tmp(NTYPE)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Read_Lai_Two_Months called by ', loc_proc
      end if

      found_month1 = .false.
      found_month2 = .false.

!c?   cache misses?
      do ij = ju1, j2
        do il = i1, i2
          do ik = 1, ireg(il,ij)

            xlai (il,ij,ik) = 0.0d0
            xlai2(il,ij,ik) = 0.0d0

          end do
        end do
      end do


!     ================
      call AsciiOpenRead  &
!     ================
     &  (lun, lai_infile_name)


      READ_LOOP: do

        Read (lun, 900, iostat = iost)  &
     &    il, ij, index_type, (xlai_tmp(it), it=1,index_type)


        if (((iost < 0) .and.  &
     &       (.not. (found_month1 .and. found_month2))) .or.  &
     &      (iost > 0)) then
          err_msg = 'Problem in Read_Lai_Data.'
          call GmiPrintError (err_msg, .true., 1, iost, 0, 0, 0.0d0, 0.0d0)
        end if


        if ((il == 0) .and.  &
     &      (.not. (found_month1 .and. found_month2))) then

          if (ij == month1) then

            found_month1 = .true.
            month_found  = month1

          else if (ij == month2) then

            found_month2 = .true.
            month_found  = month2

          end if

        else

          if ((iost < 0) .or.  &
     &        ((il == 0) .and. found_month1 .and. found_month2)) then

!           ==============
            exit READ_LOOP
!           ==============

          else if (((il >= i1 ) .and. (il <= i2)) .and.  &
     &             ((ij >= ju1) .and. (ij <= j2))) then

            if (month_found == month1) then

              xlai (il,ij,1:index_type) = xlai_tmp(1:index_type)

            else if (month_found == month2) then

              xlai2(il,ij,1:index_type) = xlai_tmp(1:index_type)

            end if

          end if

        end if

      end do READ_LOOP


      Close (lun)

 900  format (3i3, 20f5.1)

      return

      end subroutine Read_Lai_Two_Months

!-------------------------------------------------------------------

    end module ReadVegLaiData_mod

