!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiSub2Glob_mod
!
! !INTERFACE:
!
      module GmiSub2Glob_mod
!
! !USES:
      use GmiSubDomainsBC_mod, only : subDomainsBc
      use GmiMessagePassing_mod, only : sendReal8, receiveReal8
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: subDomain2Global

      interface subDomain2Global
         module procedure subDomain2Global3D
         module procedure subDomain2Global2D
      end interface
!
! !DESCRIPTION:
! Routine for communications from worker processors to the root processor.
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
! !IROUTINE: subDomain2Global3D
!
! !INTERFACE:
!
      subroutine subDomain2Global3D (globArray, locArray, ilo_gl, ihi_gl, jlo_gl, &
     &              jhi_gl, ilo, ihi, jlo, jhi, klo, khi, &
     &                 rootProc, procID, map, numDomains, msg_id, msg_comm)
!
      implicit none
!
! !INPUT PARAMETERS:
               ! the global dimensions when called by the root proc.
      integer, intent(in) :: ilo_gl, ihi_gl, jlo_gl, jhi_gl
               ! the local dimensions when called by all the procs.
      integer, intent(in) :: ilo, ihi, jlo, jhi
      integer, intent(in) :: klo, khi
               ! total number of subdomains (same as number of processors)
      integer, intent(in) :: numDomains
      integer, intent(in) :: procID   ! processor id
      integer, intent(in) :: rootProc ! root processor id
      integer, intent(in) :: msg_id   ! message id
      integer, intent(in) :: msg_comm ! message passing communicator
               ! a map defining the subdomains; for each subdomain the first
               ! dimension contains the low index and high index into the
               ! global domain for the longitudinal direction, the second
               ! dimension is the same for the latitudinal direction; the
               ! third dimension runs over subdomain number; it is important
               ! that the inclusion/omission of the border zones be consistent
               ! with the above dimensions and the array;
               ! a geometric interpretation of this array is:
               !   (1,1,numDomains) = southwest corner
               !   (2,1,numDomains) = southeast corner
               !   (1,2,numDomains) = northwest corner
               !   (2,2,numDomains) = northeast corner
      integer, intent(in) :: map(2, 2, numDomains)
               ! Local array only covering the subdomain assigned to each
               ! processor.
      real*8 , intent(in) :: locArray(ilo:ihi,jlo:jhi,klo:khi)
!
! !INPUT/OUTPUT PARAMETERS:
               ! Global array only allocated by the root processor.
      real*8 , intent(inOut) :: globArray(ilo_gl:ihi_gl,jlo_gl:jhi_gl,klo:khi)
!
! !DESCRIPTION:
!   Takes a real array defined on a subdomain and aggregates the data into a 
!   larger array defined on the global domain. This is accomplished by sending 
!   the data from the subdomain processors to the root processor in a series 
!   of messages.
!
! !LOCAL VARIABLES:
      integer :: lpx, domainx, msize
      integer :: ilo_x, ihi_x, jlo_x, jhi_x
      integer :: imax
      real*8, allocatable  :: local(:, :, :)
!EOP
!------------------------------------------------------------------------
!BOC
      if (procID == rootProc) then

         ilo_x = map(1,1,procID+1)
         ihi_x = map(2,1,procID+1)
         jlo_x = map(1,2,procID+1)
         jhi_x = map(2,2,procID+1)
       
         globArray(ilo_x:ihi_x,jlo_x:jhi_x, klo:khi) =  &
     &                locArray(ilo_x:ihi_x,jlo_x:jhi_x, klo:khi)

         ! Have global processor loop through the subdomains.

         do domainx = 2, numDomains
            lpx = domainx - 1

            ! Global processor receives a message, sorts 
            ! it into global array.
            ilo_x = map(1,1,domainx)
            ihi_x = map(2,1,domainx)
            jlo_x = map(1,2,domainx)
            jhi_x = map(2,2,domainx)

            allocate(local(ilo_x:ihi_x,jlo_x:jhi_x, klo:khi))

            ! Receive a message from the subdomain processors.

            msize = (ihi_x - ilo_x + 1) * (jhi_x - jlo_x + 1) * (khi - klo + 1)

            call receiveReal8 (lpx, msize, msg_id, local, msg_comm)

            ! Aggregate the local data into the global domain.

            globArray(ilo_x:ihi_x,jlo_x:jhi_x,klo:khi) = &
     &                       local(ilo_x:ihi_x,jlo_x:jhi_x,klo:khi)

            deallocate(local)
         end do

      else

         ! Send message to global processor.

         ilo_x = map(1,1,procID+1)
         ihi_x = map(2,1,procID+1)
         jlo_x = map(1,2,procID+1)
         jhi_x = map(2,2,procID+1)

         msize = (ihi_x - ilo_x + 1) * (jhi_x - jlo_x + 1) * (khi - klo + 1)

         call sendReal8 (rootProc, msize, msg_id, locArray, msg_comm)

      end if

      return

      end subroutine subDomain2Global3D
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subDomain2Global2D
!
! !INTERFACE:
!
      subroutine subDomain2Global2D (globArray, locArray, &
     &              ilo_gl, ihi_gl, jlo_gl, jhi_gl, ilo, ihi, jlo, jhi, &
     &                 rootProc, procID, map, numDomains, msg_id, msg_comm)
!
      implicit none
!
! !INPUT PARAMETERS:
               ! the global dimensions when called by the root proc.
      integer, intent(in) :: ilo_gl, ihi_gl, jlo_gl, jhi_gl
               ! the local dimensions when called by all the procs.
      integer, intent(in) :: ilo, ihi, jlo, jhi
               ! total number of subdomains (same as number of processors)
      integer, intent(in) :: numDomains
      integer, intent(in) :: procID   ! processor id
      integer, intent(in) :: rootProc ! root processor id
      integer, intent(in) :: msg_id   ! message id
      integer, intent(in) :: msg_comm ! message passing communicator
               ! a map defining the subdomains; for each subdomain the first
               ! dimension contains the low index and high index into the
               ! global domain for the longitudinal direction, the second
               ! dimension is the same for the latitudinal direction; the
               ! third dimension runs over subdomain number; it is important
               ! that the inclusion/omission of the border zones be consistent
               ! with the above dimensions and the array;
               ! a geometric interpretation of this array is:
               !   (1,1,numDomains) = southwest corner
               !   (2,1,numDomains) = southeast corner
               !   (1,2,numDomains) = northwest corner
               !   (2,2,numDomains) = northeast corner
      integer, intent(in) :: map(2, 2, numDomains)
               ! Local array only covering the subdomain assigned to each
               ! processor.
      real*8 , intent(in) :: locArray(ilo:ihi,jlo:jhi)
!
! !INPUT/OUTPUT PARAMETERS:
               ! Global array only allocated by the root processor.
      real*8 , intent(inOut) :: globArray(ilo_gl:ihi_gl,jlo_gl:jhi_gl)
!
! !DESCRIPTION:
!   Takes a real array defined on a subdomain and aggregates the data into a 
!   larger array defined on the global domain. This is accomplished by sending 
!   the data from the subdomain processors to the root processor in a series 
!   of messages.
!
! !LOCAL VARIABLES:
      integer :: lpx, domainx, msize
      integer :: ilo_x, ihi_x, jlo_x, jhi_x
      integer :: imax
      real*8, allocatable  :: local(:, :)
!EOP
!------------------------------------------------------------------------
!BOC
      if (procID == rootProc) then

         ilo_x = map(1,1,rootProc+1)
         ihi_x = map(2,1,rootProc+1)
         jlo_x = map(1,2,rootProc+1)
         jhi_x = map(2,2,rootProc+1)
       
         globArray(ilo_x:ihi_x,jlo_x:jhi_x) = locArray(ilo_x:ihi_x,jlo_x:jhi_x)

         ! Have global processor loop through the subdomains.

         do domainx = 2, numDomains
            lpx = domainx -1

            ! Global processor receives a message, sorts 
            ! it into global array.
            ilo_x = map(1,1,domainx)
            ihi_x = map(2,1,domainx)
            jlo_x = map(1,2,domainx)
            jhi_x = map(2,2,domainx)

            allocate(local(ilo_x:ihi_x,jlo_x:jhi_x))

            ! Receive a message from the subdomain processors.

            msize = (ihi_x - ilo_x + 1) * (jhi_x - jlo_x + 1)

            call receiveReal8 (lpx, msize, msg_id, local, msg_comm)

            ! Aggregate the local data into the global domain.

            globArray(ilo_x:ihi_x,jlo_x:jhi_x) = &
     &                       local(ilo_x:ihi_x,jlo_x:jhi_x)

            deallocate(local)
         end do

      else

         ! Send message to global processor.

         ilo_x = map(1,1,procID+1)
         ihi_x = map(2,1,procID+1)
         jlo_x = map(1,2,procID+1)
         jhi_x = map(2,2,procID+1)

         msize = (ihi_x - ilo_x + 1) * (jhi_x - jlo_x + 1)

         call sendReal8 (rootProc, msize, msg_id, locArray, msg_comm)

      end if

      return

      end subroutine subDomain2Global2D
!EOC
!-----------------------------------------------------------------------------
      end module GmiSub2Glob_mod

