!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiReduce_mod
!
! !INTERFACE:
!
      module GmiReduce_mod
!
! !USES:
      use GmiMessagePassing_mod, only : sendInteger, sendReal8
      use GmiMessagePassing_mod, only : receiveInteger, receiveReal8, writeMpiError
      use GmiPrintError_mod, only : GmiPrintError
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private 
      public   :: Gmi_Maxpij_Reduce
      public   :: Gmi_Minpij_Reduce
      public   :: Gmi_Sum_Reduce
      public ::   Gmi_Max_Reduce      , Gmi_Max1_Reduce
      public ::   Gmi_Min_Reduce      , Gmi_Min1_Reduce
      public ::   Gmi_Sum_Pole_Reduce
      public ::   Gmi_Sum1_Pole_Reduce
      
#     include "mpif.h"
#     include "gem_msg_numbers.h"
!
! !DESCRIPTION:
! Routines for performing ''reduce'' operations.
!
! !AUTHOR:
!  John Tannahill, LLNL , jrt@llnl.gov
!  Jules Kouatchou, NASA GSFC, Jules.Kouatchou-1@nasa.gov
!
!EOP
!-----------------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gmi_Maxij_Reduce
!
! !INTERFACE:
!
      subroutine Gmi_Maxpij_Reduce  &
     &  (asize, msgnum, maxi, maxj, maxv, &
     &   commuWorld, rootProc, procID, numDomains, ivert)
!
      implicit none
!
! !INPUT PARAMETER:
      integer, intent(in) :: asize
      integer, intent(in) :: commuWorld
      integer, intent(in) :: numDomains, ivert
      integer, intent(in) :: rootProc, procID
      integer, intent(in) :: msgnum(3)
!
! !INPUT/OUTPUT PARAMETER:
      integer, intent(inOut) :: maxi(asize)
      integer, intent(inOut) :: maxj(asize)
      real*8 , intent(inOut) :: maxv(asize)
!
! !DESCRIPTION:
! Visits all the processors and reduces and broadcasts
! a set of real maximums and their corresponding ij locations to the
! root processor.
!
! !LOCAL VARIABLES:
      integer :: domain
      integer :: gprx, lprx
      integer :: nn
      integer :: loc_maxi(asize)
      integer :: loc_maxj(asize)
      real*8  :: loc_maxv(asize)
!EOP
!------------------------------------------------------------------------
!BOC
      gprx = rootProc

      if (procID /= gprx) then

         call sendInteger  (gprx, asize, msgnum(1), maxi, commuWorld)
         call sendInteger  (gprx, asize, msgnum(2), maxj, commuWorld)
         call sendReal8    (gprx, asize, msgnum(3), maxv, commuWorld)
 
      else

         loc_maxi(:) = maxi(:)
         loc_maxj(:) = maxj(:)
         loc_maxv(:) = maxv(:)

         do domain = 2, numDomains

            lprx = domain-1

            call receiveInteger  (lprx, asize, msgnum(1), maxi, commuWorld)
            call receiveInteger  (lprx, asize, msgnum(2), maxj, commuWorld)
            call receiveReal8    (lprx, asize, msgnum(3), maxv, commuWorld)

            do nn = 1, ivert

               if (maxv(nn) > loc_maxv(nn)) then
                  loc_maxi(nn) = maxi(nn)
                  loc_maxj(nn) = maxj(nn)
                  loc_maxv(nn) = maxv(nn)
               end if

            end do

         end do

         maxi(:) = loc_maxi(:)
         maxj(:) = loc_maxj(:)
         maxv(:) = loc_maxv(:)

      end if

      return

      end subroutine Gmi_Maxpij_Reduce
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gmi_Minpij
!
! !INTERFACE:
!
      subroutine Gmi_Minpij_Reduce  &
     &  (asize, msgnum, mini, minj, minv, &
     &   commuWorld, rootProc, procID, numDomains, ivert)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: commuWorld
      integer, intent(in) :: numDomains, ivert
      integer, intent(in) :: rootProc, procID
      integer, intent(in) :: msgnum(3)
      integer, intent(in) :: asize
!
! !INPUT/OUTPUT PARAMETERS:
      integer, intent(inOut) :: mini(asize)
      integer, intent(inOut) :: minj(asize)
      real*8 , intent(inOut) :: minv(asize)
!
! !DESCRIPTION:
!   This routine visits all the processors and reduces and broadcasts
!   a set of real minimums and their corresponding ij locations to the
!   Root.
!
! !LOCAL VARIABLES:
      integer :: domain
      integer :: gprx, lprx
      integer :: nn
      integer :: loc_mini(asize)
      integer :: loc_minj(asize)
      real*8  :: loc_minv(asize)
!EOP
!------------------------------------------------------------------------
!BOC
      gprx  = rootProc

      if (procID /= gprx) then

         call sendInteger  (gprx, asize, msgnum(1), mini, commuWorld)
         call sendInteger  (gprx, asize, msgnum(2), minj, commuWorld)
         call sendReal8    (gprx, asize, msgnum(3), minv, commuWorld)

      else

         loc_mini(:) = mini(:)
         loc_minj(:) = minj(:)
         loc_minv(:) = minv(:)

         do domain = 2, numDomains

            lprx = domain - 1

            call receiveInteger  (lprx, asize, msgnum(1), mini, commuWorld)
            call receiveInteger  (lprx, asize, msgnum(2), minj, commuWorld)
            call receiveReal8    (lprx, asize, msgnum(3), minv, commuWorld)

            do nn = 1, ivert

               if (minv(nn) < loc_minv(nn)) then
                  loc_mini(nn) = mini(nn)
                  loc_minj(nn) = minj(nn)
                  loc_minv(nn) = minv(nn)
               end if

            end do
         end do

         mini(:) = loc_mini(:)
         minj(:) = loc_minj(:)
         minv(:) = loc_minv(:)

      end if

      return

      end subroutine Gmi_Minpij_Reduce
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gmi_Sum_Reduce
!
! !INTERFACE:
!
      subroutine Gmi_Sum_Reduce (dim1, dim2, sum, commuWorld, rootProc, procID)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: commuWorld
      integer, intent(in) :: rootProc, procID
      integer, intent(in) :: dim1
      integer, intent(in) :: dim2
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: sum(dim1:dim2)
!
! !DESCRIPTION:
!   Visits all the processors and reduces and broadcasts a set of sums 
!   to the Root.
!
! !LOCAL VARIABLES:
      integer :: count
      integer :: ierr
      real*8  :: global_sum(dim1:dim2)
!EOP
!------------------------------------------------------------------------
!BOC
      count = dim2 - dim1 + 1

!     ===============
      call MPI_Reduce  &
!     ===============
     &  (sum, global_sum, count, MPI_REAL8,  &
     &   MPI_SUM, rootProc, commuWorld, ierr)

      if (ierr /= MPI_SUCCESS) then
        call writeMpiError (commuWorld, .true., ierr)
      end if

      sum(:) = global_sum(:)

      return

      end subroutine Gmi_Sum_Reduce
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gmi_Max_Reduce
!
! !INTERFACE:
!
      subroutine Gmi_Max_Reduce (dim1, dim2, rmax, numDomains, commuWorld)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: numDomains
      integer, intent(in) :: commuWorld
      integer, intent(in) :: dim1
      integer, intent(in) :: dim2
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: rmax(dim1:dim2)
!
! !DESCRIPTION:
!   This routine visits all the processors and reduces and broadcasts
!   a set of maximum values.
!
! !LOCAL VARIABLES:
      integer :: count
      integer :: ierr
      real*8  :: global_rmax(dim1:dim2)
      real*8  :: local_rmax (dim1:dim2)
!
!EOP
!------------------------------------------------------------------------
!BOC
      if (numDomains /= 1) then

        local_rmax(:) = rmax(:)

        count = dim2 - dim1 + 1

!       ==================
        call MPI_Allreduce  &
!       ==================
     &    (local_rmax, global_rmax, count, MPI_REAL8,  &
     &     MPI_MAX, commuWorld, ierr)


        if (ierr /= MPI_SUCCESS) then
          call writeMpiError (commuWorld, .true., ierr)
        end if

        rmax(:) = global_rmax(:)

      end if

      return

      end  subroutine Gmi_Max_Reduce
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gmi_Max1_Reduce
!
! !INTERFACE:
!
      subroutine Gmi_Max1_Reduce (rmax, numDomains, commuWorld)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: numDomains
      integer, intent(in) :: commuWorld
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: rmax
!
! !DESCRIPTION:
!   This routine visits all the processors and reduces and broadcasts
!   a single maximum value.
!
! !LOCAL VARIABLES:
      integer :: ierr
      real*8  :: global_rmax, local_rmax
!
!EOP
!------------------------------------------------------------------------
!BOC

      if (numDomains /= 1) then

        local_rmax = rmax

!       ==================
        call MPI_Allreduce  &
!       ==================
     &    (local_rmax, global_rmax, 1, MPI_REAL8,  &
     &     MPI_MAX, commuWorld, ierr)


        if (ierr /= MPI_SUCCESS) then
          call writeMpiError (commuWorld, .true., ierr)
        end if

        rmax = global_rmax

      end if

      return

      end subroutine Gmi_Max1_Reduce
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Gmi_Min_Reduce
!
! !INTERFACE:
!
      subroutine Gmi_Min_Reduce  &
     &  (dim1, dim2, rmin, numDomains, commuWorld)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: numDomains
      integer, intent(in) :: commuWorld
      integer, intent(in) :: dim1
      integer, intent(in) :: dim2
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: rmin(dim1:dim2)
!
! !DESCRIPTION:
!   This routine visits all the processors and reduces and broadcasts
!   a set of minimum values.
!
! !LOCAL VARIABLES:
      integer :: count
      integer :: ierr
      real*8  :: global_rmin(dim1:dim2)
      real*8  :: local_rmin (dim1:dim2)
!EOP
!------------------------------------------------------------------------
!BOC
      if (numDomains /= 1) then

        local_rmin(:) = rmin(:)

        count = dim2 - dim1 + 1

!       ==================
        call MPI_Allreduce  &
!       ==================
     &    (local_rmin, global_rmin, count, MPI_REAL8,  &
     &     MPI_MIN, commuWorld, ierr)

        if (ierr /= MPI_SUCCESS) then
          call writeMpiError (commuWorld, .true., ierr)
        end if

        rmin(:) = global_rmin(:)

      end if

      return

      end subroutine Gmi_Min_Reduce
!EOC
!-----------------------------------------------------------------------------
!BOC
!
! !IROUTINE: Gmi_Min1_Reduce
!
! !INTERFACE:
!
      subroutine Gmi_Min1_Reduce (rmin, numDomains, commuWorld)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: numDomains
      integer, intent(in) :: commuWorld
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: rmin
!
! !DESCRIPTION:
!   This routine visits all the processors and reduces and broadcasts
!   a single minimum value.
!
! !LOCAL VARIABLES:
      integer :: ierr
      real*8  :: global_rmin, local_rmin
!EOP
!------------------------------------------------------------------------
!BOC
      if (numDomains /= 1) then

        local_rmin = rmin

!       ==================
        call MPI_Allreduce  &
!       ==================
     &    (local_rmin, global_rmin, 1, MPI_REAL8,  &
     &     MPI_MIN, commuWorld, ierr)

        if (ierr /= MPI_SUCCESS) then
          call writeMpiError (commuWorld, .true., ierr)
        end if

        rmin = global_rmin

      end if

      return

      end subroutine Gmi_Min1_Reduce
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gmi_Sum_Pole_Reduce
!
! !INTERFACE:
!
      subroutine Gmi_Sum_Pole_Reduce  &
     &  (dim1, dim2, rsum, ju1, j2, ju1_gl, j2_gl, &
     &  commu_npole, commu_spole)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: ju1, j2, ju1_gl, j2_gl
      integer, intent(in) :: commu_npole, commu_spole
      integer, intent(in) :: dim1
      integer, intent(in) :: dim2
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: rsum(dim1:dim2)
!
! !DESCRIPTION:
!   This routine visits all the processors at either the North or South
!   Pole and reduces and broadcasts an array of sums.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      integer :: commu
      integer :: count
      integer :: ierr
      real*8  :: global_rsum(dim1:dim2)
      real*8  :: local_rsum (dim1:dim2)
!EOP
!------------------------------------------------------------------------
!BOC
      local_rsum(:) = rsum(:)
      count = dim2 - dim1 + 1

      if (((ju1 == ju1_gl) .and. (j2 == j2_gl)) .or.  &
     &    ((ju1 /= ju1_gl) .and. (j2 /= j2_gl))) then

        err_msg = 'Pole problem in Gmi_Sum_Pole_Reduce.'
        call GmiPrintError (err_msg, .true., 2, ju1, j2, 0, 0.0d0, 0.0d0)

      else if (ju1 == ju1_gl) then

        commu = commu_spole

      else if (j2 == j2_gl) then

        commu = commu_npole

      end if

!     ==================
      call MPI_Allreduce  &
!     ==================
     &  (local_rsum, global_rsum, count, MPI_REAL8,  &
     &   MPI_SUM, commu, ierr)

      if (ierr /= MPI_SUCCESS) then
        call writeMpiError (commu, .true., ierr)
      end if

      rsum(:) = global_rsum(:)

      return

      end subroutine Gmi_Sum_Pole_Reduce
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gmi_Sum1_Pole_Reduce
!
! !INTERFACE:
!
      subroutine Gmi_Sum1_Pole_Reduce (rsum, &
     &         numLonDomains, ju1, j2, ju1_gl, j2_gl, &
     &         commu_npole, commu_spole)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: numLonDomains
      integer, intent(in) :: ju1, j2, ju1_gl, j2_gl
      integer, intent(in) :: commu_npole, commu_spole
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: rsum
!
! !DESCRIPTION:
!   This routine visits all the processors at either the North or South
!   Pole and reduces and broadcasts a single sum.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      integer :: commu
      integer :: ierr
      real*8  :: global_rsum, local_rsum
!EOP
!------------------------------------------------------------------------
!BOC

      if (numLonDomains /= 1) then

        local_rsum = rsum

        if (((ju1 == ju1_gl) .and. (j2 == j2_gl)) .or.  &
     &      ((ju1 /= ju1_gl) .and. (j2 /= j2_gl))) then

          err_msg = 'Pole problem in Gmi_Sum1_Pole_Reduce.'
          call GmiPrintError (err_msg, .true., 2, ju1, j2, 0, 0.0d0, 0.0d0)

        else if (ju1 == ju1_gl) then

          commu = commu_spole

        else if (j2 == j2_gl) then

          commu = commu_npole

        end if

!       ==================
        call MPI_Allreduce  &
!       ==================
     &    (local_rsum, global_rsum, 1, MPI_REAL8,  &
     &     MPI_SUM, commu, ierr)


        if (ierr /= MPI_SUCCESS) then
          call writeMpiError (commu, .true., ierr)
        end if

        rsum = global_rsum

      end if

      return

      end subroutine Gmi_Sum1_Pole_Reduce
!EOC
!------------------------------------------------------------------------

      end module GmiReduce_mod
