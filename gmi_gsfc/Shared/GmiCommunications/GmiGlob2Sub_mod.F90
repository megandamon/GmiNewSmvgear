!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiGlob2Sub_mod
!
! !INTERFACE:
!
      module GmiGlob2Sub_mod
!
! !USES:
      use GmiMessagePassing_mod, only : sendReal8, receiveReal8
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: global2Subdomain, global2Subdomain2D
!
!
! !DESCRIPTION:
! Routine for communications from the root processor to worker processors.
!
! !AUTHOR:
!  Jules Kouatchou, NASA GSFC, Jules.Kouatchou-1@nasa.gov
!
!EOP
!-----------------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: global2Subdomain
!
! !INTERFACE:
!
      subroutine global2Subdomain (globArray, locArray, ilo_gl, ihi_gl, &
     &                 jlo_gl, jhi_gl, ilo, ihi, jlo, jhi, klo, khi,  &
     &                 rootProc, procID, map, numDomains, msg_comm,  msg_id)
!
      implicit none
!
               ! the global dimensions when called by the root proc.;
      integer, intent(in) :: ilo_gl, ihi_gl, jlo_gl, jhi_gl
               ! the local dimensions when called by all procs.
      integer, intent(in) :: ilo, ihi, jlo, jhi
      integer, intent(in) :: klo, khi
               ! total number of subdomains (same as number of processors)
      integer, intent(in) :: numDomains
      integer, intent(in) :: procID ! processor id
      integer, intent(in) :: rootProc ! root processor id
      integer, intent(in) :: msg_id ! message identifier
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
               ! Global array only allocated by the root processor.
      real*8 , intent(in) :: globArray(ilo_gl:ihi_gl, jlo_gl:jhi_gl,klo:khi)
!
! !INPUT/OUTPUT PARAMETERS:
               ! Local array only covering local domain on each processor.
      real*8 , intent(inOut) :: locArray(ilo:ihi, jlo:jhi,klo:khi)
!
! !DESCRIPTION:
! Takes a real array defined on the global domain and segregates the data 
! into smaller arrays defined by the subdomain decomposition. These arrays are 
! then sent by the global and received by the subdomain processors in the form 
! of a message.
!
! !LOCAL VARIABLES:
      integer :: domainx, lpx, msize
      integer :: ilo_x, ihi_x, jlo_x, jhi_x
      real*8, allocatable  :: local(:, :, :)
!
!EOP
!------------------------------------------------------------------------
!BOC

!     -----------------------
!     Global processor tasks.
!     -----------------------

      if (procID == rootProc) then

         ilo_x = map(1,1,rootProc+1)
         ihi_x = map(2,1,rootProc+1)
         jlo_x = map(1,2,rootProc+1)
         jhi_x = map(2,2,rootProc+1)
         
         locArray(ilo_x:ihi_x,jlo_x:jhi_x, klo: khi) =  &
     &                globArray(ilo_x:ihi_x,jlo_x:jhi_x, klo: khi)

         !----------------------------
         ! Loop through the subdomains.
         !----------------------------

         do domainx = 2, numDomains

            ! ---------------------------------------------------------------
            ! Extract the relevant portion of the global array into a local
            ! scratch array of the proper dimensions for each subdomain.
            !---------------------------------------------------------------

            lpx = domainx - 1

            ilo_x = map(1,1,domainx)
            ihi_x = map(2,1,domainx)
            jlo_x = map(1,2,domainx)
            jhi_x = map(2,2,domainx)

            allocate(local(ilo_x:ihi_x,jlo_x:jhi_x, klo: khi))

            local(ilo_x:ihi_x,jlo_x:jhi_x, klo: khi) =  &
     &                globArray(ilo_x:ihi_x,jlo_x:jhi_x, klo: khi)

            ! -----------------------------------------------
            ! Send the message from the root processor to the 
            ! subdomain processor.
            ! -----------------------------------------------

            msize = (ihi_x - ilo_x + 1) * (jhi_x - jlo_x + 1) *  &
     &                 (khi - klo + 1)

            call sendReal8 (lpx, msize, msg_id, local, msg_comm)

            deallocate(local)

          end do

      else

         ! --------------------------
         ! Subdomain processor tasks.
         ! --------------------------

          ilo_x = map(1,1,procID + 1)
          ihi_x = map(2,1,procID + 1)
          jlo_x = map(1,2,procID + 1)
          jhi_x = map(2,2,procID + 1)

         ! ----------------------------------------------
         ! Receive the message from the global processor.
         ! ----------------------------------------------

         msize = (ihi_x - ilo_x + 1) * (jhi_x - jlo_x + 1) * (khi - klo + 1)
         call receiveReal8 (rootProc, msize, msg_id, locArray, msg_comm)

      end if

      return

      end subroutine global2Subdomain
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: global2Subdomain2D
!
! !INTERFACE:
!
      subroutine global2Subdomain2D(globArray, locArray, &
     &                 ilo_gl, ihi_gl, jlo_gl, jhi_gl, ilo, ihi, jlo, jhi, &
     &                 rootProc, procID, map, numDomains, msg_comm,  msg_id)
!
      implicit none
!
      integer, intent(in) :: ilo, ihi, jlo, jhi
      integer, intent(in) :: ilo_gl, ihi_gl, jlo_gl, jhi_gl
               ! total number of subdomains (same as number of processors)
      integer, intent(in) :: numDomains
      integer, intent(in) :: procID ! processor id
      integer, intent(in) :: rootProc ! root processor id
      integer, intent(in) :: msg_id ! message identifier
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
               ! Global array only allocated by the root processor.
      real*8 , intent(in) :: globArray(ilo_gl:ihi_gl,jlo_gl:jhi_gl)
!
! !INPUT/OUTPUT PARAMETERS:
               ! Local array only covering local domain on each processor.
      real*8 , intent(inOut) :: locArray(ilo:ihi,jlo:jhi)
!
! !DESCRIPTION:
! Takes a real array defined on the global domain and segregates the data 
! into smaller arrays defined by the subdomain decomposition. These arrays are 
! then sent by the global and received by the subdomain processors in the form 
! of a message.
!
! !LOCAL VARIABLES:
      integer :: domainx, lpx, msize
      integer :: ilo_x, ihi_x, jlo_x, jhi_x
      real*8, allocatable  :: local(:, :)
!
!EOP
!------------------------------------------------------------------------
!BOC

!     -----------------------
!     Global processor tasks.
!     -----------------------

      if (procID == rootProc) then

         ilo_x = map(1,1,rootProc+1)
         ihi_x = map(2,1,rootProc+1)
         jlo_x = map(1,2,rootProc+1)
         jhi_x = map(2,2,rootProc+1)
         
         locArray(ilo_x:ihi_x,jlo_x:jhi_x) = globArray(ilo_x:ihi_x,jlo_x:jhi_x)

         !----------------------------
         ! Loop through the subdomains.
         !----------------------------

         do domainx = 2, numDomains

            ! ---------------------------------------------------------------
            ! Extract the relevant portion of the global array into a local
            ! scratch array of the proper dimensions for each subdomain; 
            !---------------------------------------------------------------

            lpx = domainx - 1

            ilo_x = map(1,1,domainx)
            ihi_x = map(2,1,domainx)
            jlo_x = map(1,2,domainx)
            jhi_x = map(2,2,domainx)

            allocate(local(ilo_x:ihi_x,jlo_x:jhi_x))

            local(ilo_x:ihi_x,jlo_x:jhi_x) =  &
     &                globArray(ilo_x:ihi_x,jlo_x:jhi_x)

            ! -----------------------------------------------
            ! Send the message from the root processor to the 
            ! subdomain processor.
            ! -----------------------------------------------

            msize = (ihi_x - ilo_x + 1) * (jhi_x - jlo_x + 1) 

            call sendReal8 (lpx, msize, msg_id, local, msg_comm)

            deallocate(local)

         end do

      else

         ! --------------------------
         ! Subdomain processor tasks.
         ! --------------------------

          ilo_x = map(1,1,procID + 1)
          ihi_x = map(2,1,procID + 1)
          jlo_x = map(1,2,procID + 1)
          jhi_x = map(2,2,procID + 1)

         ! ----------------------------------------------
         ! Receive the message from the global processor.
         ! ----------------------------------------------

         msize = (ihi_x - ilo_x + 1) * (jhi_x - jlo_x + 1)
         call receiveReal8 (rootProc, msize, msg_id, locArray, msg_comm)

      end if

      return

      end subroutine global2Subdomain2D
!EOC
!-----------------------------------------------------------------------------
      end module GmiGlob2Sub_mod
