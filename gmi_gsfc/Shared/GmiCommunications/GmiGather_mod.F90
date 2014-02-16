!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiGather_mod
!
! !INTERFACE:
!
      module GmiGather_mod
!
! !USES:
      use GmiPrintError_mod, only : GmiPrintError
      use GmiMessagePassing_mod, only :  writeMpiError
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public   :: Gmi_Pole_Allgather
!
#     include "mpif.h"
!
! !AUTHOR:
! John Tannahill, LLNL, jrt@llnl.gov
! Jules Kouatchou, NASA/GSFC, Jules.Kouatchou-1@nasa.gov
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gmi_Pole_Allgather
!
! !INTERFACE:
!
      subroutine Gmi_Pole_Allgather  &
     &  (jst, jend, send_array, recv_array, &
     &  mapi_all, commu_spole, commu_npole, numDomains, numLonDomains, &
     &  i1, i2, ju1, j2, k1, k2, ivert, ilo, ihi, julo, jhi, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: jst
                             ! starting latitude dimension of global recv_array
      integer, intent(in) :: jend
                             ! ending   latitude dimension of global recv_array
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ivert
      integer, intent(in) :: ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: numDomains, numLonDomains
      integer, intent(in) :: commu_spole, commu_npole
      integer, intent(in) :: mapi_all(2,numDomains)
      real*8 , intent(in) :: send_array(ilo:ihi, julo:jhi, k1:k2)
                             ! local array that contains the needed info on a processor
                             ! (local longitude dimensions)
      real*8, intent(out)  :: recv_array(i1_gl:i2_gl, jst:jend, k1:k2)
!      real*8, intent(out)  :: recv_array(k1:k2, jst:jend, i1_gl:i2_gl)
                             ! global array to fill in (global longitude dimensions)
!
! !DESCRIPTION:
!   This routine takes a 3D local array at a Pole and collects the data from
!   each processor.  It then puts all of this data on each processor in a
!   longitudinally global array.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      integer :: commu
      integer :: count
      integer :: ierr
      integer :: il, ij, ik, id, ilcnt, ic, idim, jdim, kdim
      integer :: scale
      integer :: istrts(numLonDomains)
      integer :: counts(0:numLonDomains-1)
      integer :: displs(0:numLonDomains-1)
      real*8  :: rterms((i2_gl-i1_gl+1) * (jend-jst+1) * IVERT)
      real*8  :: sterms((i2   -i1   +1) * (jend-jst+1) * IVERT)
!EOP
!--------------------------------------------------------------------------
!BOC
!
      counts(:) = 0
      displs(:) = 0
      rterms(:) = -99999.0d0
      sterms(:) = -99999.0d0
!
      scale = (jend - jst + 1) * IVERT
!
      count = 0
!
      do ik = 1, ivert
        do ij = jst, jend
          do il = i1, i2
!
            count = count + 1
!
!            send_array_rev_nb(ik,ij,il) = send_array(il,ij,ik)
!
            sterms(count) = send_array(il,ij,ik)
!
          end do
        end do
      end do
!
!
!     ==========================================
      if (((ju1 == ju1_gl) .and. (j2 == j2_gl)) .or.  &
     &    ((ju1 /= ju1_gl) .and. (j2 /= j2_gl))) then
!     ==========================================
!
        err_msg = 'Pole problem in Gmi_Pole_Allgather.'
!
        call GmiPrintError (err_msg, .true., 2, ju1, j2, 0, 0.0d0, 0.0d0)
!
!     =======================
      else if (ju1 == ju1_gl) then
!     =======================
!
        commu = commu_spole
!
        do il = 0, numLonDomains - 1
!
          counts(il) = (mapi_all(2,numDomains-numLonDomains+il+1) -  &
     &                  mapi_all(1,numDomains-numLonDomains+il+1) + 1) *  &
     &                 scale
!
          if (il == 0) then
            displs(0) = 0
          else
            displs(il) = displs(il-1) + counts(il-1)
          end if
!
        end do
!
!     =====================
      else if (j2 == j2_gl) then
!     =====================
!
        commu = commu_npole
!
        do il = 0, numLonDomains - 1
!
          counts(il) = (mapi_all(2,il+1) -  &
     &                  mapi_all(1,il+1) + 1) *  &
     &                 scale
!
          if (il == 0) then
            displs(0) = 0
          else
            displs(il) = displs(il-1) + counts(il-1)
          end if
!
        end do
!
      end if
!
!     ===================
      call MPI_Allgatherv  &
!     ===================
     &  (sterms, count,          MPI_REAL8,  &
     &   rterms, counts, displs, MPI_REAL8, commu, ierr)
!
      if (ierr /= MPI_SUCCESS) then
        call writeMpiError (commu, .true., ierr)
      end if
!
      do id = 1, numLonDomains
        istrts(id) = (i2-i1+1)*(id-1)
      enddo
!
      ilcnt = 0
      do id = 1, numLonDomains
        do ik = 1, ivert
          do ij = jst, jend
            do il = 1, i2-i1+1
              ilcnt = ilcnt+1
              recv_array(il+istrts(id),ij,ik) = rterms(ilcnt)
            end do
          end do
        end do
      end do
!
      return
!
      end subroutine Gmi_Pole_Allgather
!EOC
!----------------------------------------------------------------------------
      end module GmiGather_mod
!
