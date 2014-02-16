!---------------------------------------------------
! HISTORY
!  - July 1, 2004 - Jules Kouatchou
!      o Remove MPI statements when MSG_OPTION is
!        set to MSG_NONE.
!---------------------------------------------------

!     ==================
      module Ftiming_Dao
!     ==================

      use GmiPrintError_mod, only : GmiPrintError
      use GmiFlush_mod     , only : GmiFlush

      implicit none


      integer, parameter, private :: NBLKS = 100


      character (len=24), private :: blkname(NBLKS)

      integer, private      :: tblk

      type tms
        private
        real*8 :: usr, sys
      end type tms

      type (tms), private   :: accum(NBLKS)
      type (tms), private   :: last (NBLKS)

      real*8, private       :: us_tmp1(NBLKS, 2)
      real*8, private       :: us_tmp2(NBLKS, 2)

      integer               :: ierwtime

      real*8, external      :: Mpi_Wtime

      contains


!     -----------------------
      subroutine Ftiming_Init ( )
!     -----------------------

      implicit none


      integer :: nn
      real*8  :: wclk


      tblk = 0

      do nn = 1, NBLKS

        accum(nn)%usr = 0.0d0
        accum(nn)%sys = 0.0d0

        last (nn)%usr = 0.0d0
        last (nn)%sys = 0.0d0

      end do


!     ------------------------------------------
!     To reduce the overhead for the first call.
!     ------------------------------------------

      wclk = Mpi_Wtime (ierwtime)


      return

      end subroutine Ftiming_Init


!     ---------------------
      subroutine Ftiming_On (blk_name)
!     ---------------------

      implicit none


      character (len=*)  :: blk_name


      character (len=24) :: ctmp

      integer :: iblk
      integer :: ii

      real*8  :: wclk


      ctmp = Trim (blk_name)

      iblk = 0

      do ii = 1, tblk

        if (ctmp == blkname(ii)) then
          iblk = ii
        end if

      end do

      if (iblk == 0) then

        tblk = tblk + 1
        iblk = tblk

        blkname(iblk) = Trim (blk_name)

      end if


      wclk = Mpi_Wtime (ierwtime)

      last(iblk)%usr = wclk
      last(iblk)%sys = 0.0d0


      return

      end subroutine Ftiming_On


!     ----------------------
      subroutine Ftiming_Off (blk_name)
!     ----------------------

      implicit none


      character (len=*)  :: blk_name


      character (len=24) :: ctmp

      integer :: iblk
      integer :: ii

      real*8  :: wclk


      ctmp = Trim (blk_name)

      iblk = 0

      do ii = 1, tblk

        if (ctmp == blkname(ii)) then
          iblk = ii
        end if

      end do

      if (iblk == 0) then
         Write (6,*) 'Stopping in Ftiming_Off in ', ctmp
         call GmiPrintError ('Ftiming_Off', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if


      wclk = Mpi_Wtime (ierwtime)

      accum(iblk)%usr = accum(iblk)%usr + wclk - last(iblk)%usr
      accum(iblk)%sys = 0.0d0

      last(iblk)%usr  = wclk
      last(iblk)%sys  = 0.0d0

      return

      end subroutine Ftiming_Off


!     ------------------------
      subroutine Ftiming_Reset (blk_name)
!     ------------------------

      implicit none


      character (len=*)  :: blk_name


      character (len=24) :: ctmp

      integer :: iblk
      integer :: ii


      ctmp = Trim (blk_name)

      iblk = 0

      do ii = 1, tblk

        if (ctmp == blkname(ii)) then
          iblk = ii
        end if

      end do

      if (iblk == 0) then
        Write (6,*) 'Stopping in Ftiming_Reset in ', ctmp
        call GmiPrintError ('Ftiming_Reset', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      accum(iblk)%usr = 0.0d0
      accum(iblk)%sys = 0.0d0

      last (iblk)%usr = 0.0d0
      last (iblk)%sys = 0.0d0


      return

      end subroutine Ftiming_Reset


!     ----------------------
      subroutine Ftiming_Prt &
!     ----------------------
     & (lu, numWorkerProcs, numLonProcs, numLatProcs, procID, commuWorld)

      implicit none

#     include "mpif.h"

      integer :: lu
      integer, intent(in) :: numWorkerProcs, numLonProcs, numLatProcs
      integer, intent(in) :: procID, commuWorld


      character (len=8)  :: cdate
      character (len=10) :: ctime
      character (len=23) :: cfilename

      integer :: ierr
      integer :: ii, jj, kk
      integer :: ijk
      integer :: ind
      integer :: nn

      real*8  :: cijk, sijk
      real*8  :: onpi, onpij, onpj

      real*8, allocatable  :: ctavg (:)
      real*8, allocatable  :: ctmax (:)
      real*8, allocatable  :: ctmin (:)

      real*8, allocatable  :: ctavgj(:, :)
      real*8, allocatable  :: ctmaxj(:, :)
      real*8, allocatable  :: ctminj(:, :)

      real*8, allocatable  :: ctavgi(:, :)
      real*8, allocatable  :: ctmaxi(:, :)
      real*8, allocatable  :: ctmini(:, :)

      real*8, allocatable  :: us_glob1(:, :)
      real*8, allocatable  :: us_glob2(:, :)

      allocate (ctavg(tblk))
      allocate (ctmax(tblk))
      allocate (ctmin(tblk))

      allocate (ctavgj(tblk, numLonProcs))
      allocate (ctmaxj(tblk, numLonProcs))
      allocate (ctminj(tblk, numLonProcs))

      allocate (ctavgi(tblk, numLatProcs))
      allocate (ctmaxi(tblk, numLatProcs))
      allocate (ctmini(tblk, numLatProcs))

      allocate (us_glob1(tblk, numWorkerProcs))
      allocate (us_glob2(tblk, numWorkerProcs))

      us_glob1(:,:) = 0.0d0
      us_glob2(:,:) = 0.0d0

      do nn = 1, tblk
        us_tmp1(nn,1) = accum(nn)%usr
        us_tmp1(nn,2) = accum(nn)%sys

        us_glob1(nn,procID+1) = us_tmp1(nn,1) + us_tmp1(nn,2)
      end do

!     ==================
      call Mpi_Allreduce  &
!     ==================
     &  (us_glob1, us_glob2, tblk*numWorkerProcs, MPI_REAL8, MPI_SUM,  &
     &   commuWorld, ierr)


      onpi  = 1.0d0 / numLonProcs
      onpj  = 1.0d0 / numLatProcs
      onpij = 1.0d0 / numWorkerProcs

!      do nn = 1, NBLKS
      do nn = 1, tblk

        ctavg(nn) = us_glob2(nn,1)
        ctmax(nn) = us_glob2(nn,1)
        ctmin(nn) = us_glob2(nn,1)

        do kk = 2, numWorkerProcs

          ctavg(nn) = ctavg(nn) + us_glob2(nn,kk)
          ctmax(nn) = Max (ctmax(nn), us_glob2(nn,kk))
          ctmin(nn) = Min (ctmin(nn), us_glob2(nn,kk))

        end do


        ctavg(nn) = onpij * ctavg(nn)


        do jj = 1, numLatProcs

          ind = (jj - 1) * numLonProcs + 1

          ctavgi(nn,jj) = us_glob2(nn,ind)
          ctmaxi(nn,jj) = us_glob2(nn,ind)
          ctmini(nn,jj) = us_glob2(nn,ind)

          do ii = 2, numLonProcs

            ind = (jj - 1) * numLonProcs + ii

            ctavgi(nn,jj) = ctavgi(nn,jj) + us_glob2(nn,ind)
            ctmaxi(nn,jj) = Max (ctmaxi(nn,jj), us_glob2(nn,ind))
            ctmini(nn,jj) = Min (ctmini(nn,jj), us_glob2(nn,ind))

          end do

          ctavgi(nn,jj) = onpi * ctavgi(nn,jj)

        end do


        do ii = 1, numLonProcs

          ctavgj(nn,ii) = us_glob2(nn,ii)
          ctmaxj(nn,ii) = us_glob2(nn,ii)
          ctminj(nn,ii) = us_glob2(nn,ii)

          do jj = 2, numLatProcs

            ind = (jj - 1) * numLonProcs + ii

            ctavgj(nn,ii) = ctavgj(nn,ii) + us_glob2(nn,ind)
            ctmaxj(nn,ii) = Max (ctmaxj(nn,ii), us_glob2(nn,ind))
            ctminj(nn,ii) = Min (ctminj(nn,ii), us_glob2(nn,ind))

          end do

          ctavgj(nn,ii) = onpj * ctavgj(nn,ii)

        end do

      end do


!     ==============
      call GmiFlush (6)
!     ==============

!     ================
      call Mpi_Barrier (commuWorld, ierr)
!     ================

      if (procID == 0) then

!       ==================
        call Date_And_Time (cdate, ctime)
!       ==================

        cfilename( 1:8)  = 'ftiming_'
        cfilename( 9:16) = cdate(1:8)
        cfilename(17:17) = '_'
        cfilename(18:23) = ctime(1:6)

        if (lu /= 6) Open (lu, file = cfilename)

        Write (lu,*)
        Write (lu,*)
        Write (lu,*) '-------------------------------------------------'
        Write (lu,*) '           Beginning of Timing Statistics'
        Write (lu,*) '-------------------------------------------------'
        Write (lu,*)
        Write (lu,*)
        Write (lu,*) 'Timing statistics are below, on a per-process basis for'
        Write (lu,*) 'various sections of code. When viewed over a very short'
        Write (lu,*) 'period of time, these statistics are useful for analyzing'
        Write (lu,*) 'load imbalances. However, for longer periods of time, the'
        Write (lu,*) 'load distribution will vary, and simply comparing the'
        Write (lu,*) 'times among the various tasks is no longer valid.'
        Write (lu,*) 'Instead, one must compare the MPI barrier time to the'
        Write (lu,*) 'overall gmi_step time to get a sense of the effect of'
        Write (lu,*) 'load imbalances.'
        Write (lu,*)
        Write (lu,*) 'Several MPI barrier calls have been added to the time'
        Write (lu,*) 'stepping routine, most notably one at the end of the'
        Write (lu,*) 'routine, which is likely to have effect due to the load'
        Write (lu,*) 'imbalances of the Chemistry. These barriers might or'
        Write (lu,*) 'might not affect performance depending on barrier'
        Write (lu,*) 'synchronization outside of the time stepping routine.'
        Write (lu,*)
        Write (lu,*) 'The detail below contains data with respect to al MPI'
        Write (lu,*) 'tasks as well as latitudinal variations at fixed'
        Write (lu,*) 'longitude and longitudinal variations at fixed latitude.'
        Write (lu,*)
        Write (lu,*) 'For questions, contact Art Mirin, mirin@llnl.gov.'
        Write (lu,*)
        Write (lu,*)  &
     &    '  ---------------------------------------------------------',  &
     &    '--------'
        Write (lu,*)  &
     &    '     Block                       Min Time    Max Time    ',  &
     &    'Avg Time'
        Write (lu,*)  &
     &    '  ---------------------------------------------------------',  &
     &    '--------'

        do nn = 1, tblk
          Write (lu,900) blkname(nn), ctmin(nn), ctmax(nn), ctavg(nn)
        end do

 900    format (3x, a24, 3x, 3f12.4)

        Write (lu,*)
        Write (lu,*)  &
     &    '  ---------------------------------------------------------',  &
     &    '--------'
        Write (lu,*)  &
     &    '                    Longitudinal Accounting                 '
        Write (lu,*)  &
     &    '  ---------------------------------------------------------',  &
     &    '--------'

        do ii = 1, numLonProcs

          Write (lu,*)
          Write (lu,910) ii
          Write (lu,*)

          do nn = 1, tblk

            Write (lu,900)  &
     &        blkname(nn), ctminj(nn,ii), ctmaxj(nn,ii), ctavgj(nn,ii)

          end do

        end do

 910    format ('  i = ', i5, '                            ',  &
     &          'Min         Max         Avg')

        Write (lu,*)
        Write (lu,*)  &
     &    '  ---------------------------------------------------------',  &
     &    '--------'
        Write (lu,*)  &
     &    '                    Latitudinal Accounting                  '
        Write (lu,*)  &
     &    '  ---------------------------------------------------------',  &
     &    '--------'

        do jj = 1, numLatProcs

          Write (lu,*)
          Write (lu,920) jj
          Write (lu,*)

          do nn = 1, tblk

            Write (lu,900)  &
     &        blkname(nn), ctmini(nn,jj), ctmaxi(nn,jj), ctavgi(nn,jj)

          end do

        end do

 920    format ('  j = ', i5, '                            ',  &
     &          'Min         Max         Avg')

        Write (lu,*)
        Write (lu,*)  &
     &    '  ---------------------------------------------------------',  &
     &    '------------------'
        Write (lu,*)  &
     &    '                    Per-Process Detail                      '
        Write (lu,*)  &
     &    '  ---------------------------------------------------------',  &
     &    '------------------'

        do nn = 1, tblk
          Write (lu,990) blkname(nn), (us_glob2(nn,kk),kk=1,numWorkerProcs)
        end do

 990    format (3x, a24, 3x, 4f12.4, /, (30x, 4f12.4))

        Write (lu,*)
        Write (lu,*)
        Write (lu,*) '-------------------------------------------------'
        Write (lu,*) '           End of Timing Statistics'
        Write (lu,*) '-------------------------------------------------'
        Write (lu,*)
        Write (lu,*)

!       ==============
        call GmiFlush (lu)
!       ==============

        if (lu /= 6) Close (lu)

      else if (lu == 6) then

!c?
!       --------------------------------------------------------------------
!       Unnecessary calculations in order to allow time for unit 6 to flush.
!       --------------------------------------------------------------------

        ijk  = 1000

        cijk = 1.0d0
        sijk = 1.0d0 / ijk

        do kk = 1, ijk
          do jj = 1, ijk
            do ii = 1, ijk
              cijk = cijk * Sin (sijk*ii*jj*kk)
            end do
          end do
        end do

!       ====================
        call Ftiming_Nothing (cijk)
!       ====================

      end if

      deallocate (ctavg)
      deallocate (ctmax)
      deallocate (ctmin)

      deallocate (ctavgj)
      deallocate (ctmaxj)
      deallocate (ctminj)

      deallocate (ctavgi)
      deallocate (ctmaxi)
      deallocate (ctmini)

      deallocate (us_glob1)
      deallocate (us_glob2)

!     ================
      call Mpi_Barrier  &
!     ================
     &  (commuWorld, ierr)

      return

      end subroutine Ftiming_Prt


!     --------------------------
      subroutine Ftiming_Nothing (aa)
!     --------------------------

      implicit none


      real*8  :: aa


      real*8  :: bb


      bb = Exp (aa)


      if (bb < 0.0d0) Write (6,*) bb


      return

      end subroutine Ftiming_Nothing


!     ======================
      end module Ftiming_Dao
!     ======================

