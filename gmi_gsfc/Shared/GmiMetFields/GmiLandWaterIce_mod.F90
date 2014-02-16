!=============================================================================
! NASA GSFC - Code 610.3
!=============================================================================
!BOP
!
! !MODULE:
!
      module GmiLandWaterIce_mod
!
! !USES:
      use GmiPrintError_mod, only : GmiPrintError
      use GmiGrid_mod, only : t_gmiGrid, Get_i1, Get_i2, Get_ju1, Get_j2, &
     &       Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl
!
      implicit none

      private 
      public   :: setLWIflags, setCMIflags, computeGlobCMIflags1

#     include "gmi_lwi_data_64.h"
#     include "gmi_lwi_data_91.h"

!=============================================================================
      contains
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   setLWIFlags
!
! DESCRIPTION
!   Get the LWI flag data supplied by NASA-Goddard.
!
!   These flags show whether over land, water, or ice:
!
!     lwi_flags(i,j) = 1 for water
!     lwi_flags(i,j) = 2 for land
!     lwi_flags(i,j) = 3 for ice
!
! ARGUMENTS
!   lwi_flags : array of flags that indicate land, water, or ice
!
!-----------------------------------------------------------------------------

      subroutine setLWIFlags (lwi_flags, &
     &             i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl
      integer :: lwi_flags(i1:i2, ju1:j2)

!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      integer :: il, ij
      integer :: sumsq

!     ----------------
!     Begin execution.
!     ----------------

!     ==========================================
      if (((i2_gl - i1_gl)  == (L64I2 - L64I1)) .and.  &
     &    ((j2_gl - ju1_gl) == (L64J2 - L64J1))) then
!     ==========================================

        lwi_flags(i1:i2, ju1:j2) = lwi_data_64(i1:i2, ju1:j2)

!     ===============================================
      else if (((i2_gl - i1_gl)  == (L91I2 - L91I1)) .and.  &
     &         ((j2_gl - ju1_gl) == (L91J2 - L91J1))) then
!     ===============================================

        lwi_flags(i1:i2, ju1:j2) = lwi_data_91(i1:i2, ju1:j2)

!     ======================================
      else if (((i2_gl - i1_gl  + 1) == 72) .and.  &
     &         ((j2_gl - ju1_gl + 1) == 46)) then
!     ======================================

        do ij = ju1, j2
          do il = i1, i2

            if (ij < j2_gl) then

              sumsq = lwi_data_91((il*2-1),(ij*2-1))**2 +  &
     &                lwi_data_91((il*2-1),(ij*2  ))**2 +  &
     &                lwi_data_91((il*2  ),(ij*2-1))**2 +  &
     &                lwi_data_91((il*2  ),(ij*2  ))**2

            else

              sumsq = lwi_data_91((il*2-1),(ij*2-1))**2 +  &
     &                lwi_data_91((il*2  ),(ij*2-1))**2 +  &
     &                lwi_data_91((il*2-1),(ij*2-1))**2 +  &
     &                lwi_data_91((il*2  ),(ij*2-1))**2

            end if

            if      ((sumsq ==  4) .or. (sumsq ==  7) .or.  &
     &               (sumsq == 10) .or. (sumsq == 12) .or.  &
     &               (sumsq == 15)) then

              lwi_flags(il,ij) = 1

            else if ((sumsq == 13) .or. (sumsq == 16) .or.  &
     &               (sumsq == 18) .or. (sumsq == 21) .or.  &
     &               (sumsq == 26)) then

              lwi_flags(il,ij) = 2

            else if ((sumsq == 20) .or. (sumsq == 23) .or.  &
     &               (sumsq == 28) .or. (sumsq == 31) .or.  &
     &               (sumsq == 36)) then

              lwi_flags(il,ij) = 3

            else

              err_msg = 'Error finding 4x5 lwi data in Get_Lwi_Flags.'
              call GmiPrintError  &
     &          (err_msg, .true., 2, ju1_gl, j2_gl, 0, 0.0d0, 0.0d0)

            end if

          end do
        end do


!     ====
      else
!     ====

        err_msg = 'LWI indexing problem in Get_Lwi_Flags.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, ju1_gl, j2_gl, 0, 0.0d0, 0.0d0)

      end if


      return

      end subroutine setLWIFlags

!================================================================================

!-------------------------------------------------------------------
!c
!c ROUTINE
!c   setCMIflags
!c
!c DESCRIPTION
!c   Get the LWI flag data supplied by NASA-Goddard.
!c   Determine if location is marine, continental, or ice.
!c
!c
!c     Dale's routine implemented
!c     cmi_flags(i,j) = 1 for marine location
!c     cmi_flags(i,j) = 2 for continental location
!c     cmi_flags(i,j) = 3 for ice
!c
!c ARGUMENTS
!c   cmi_flags : array of flags that indicate continental, marine, or ice
!c
!c-----------------------------------------------------------------------------

      subroutine setCMIflags (cmi_flags, &
     &            i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl
      integer :: cmi_flags(i1:i2, ju1:j2)

!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      integer :: il, ij
      integer :: sumsq

!     ----------------
!     Begin execution.
!     ----------------

!     ==========================================
      if (((i2_gl - i1_gl)  == (L64I2 - L64I1)) .and.  &
     &    ((j2_gl - ju1_gl) == (L64J2 - L64J1))) then
!     ==========================================

        cmi_flags(i1:i2, ju1:j2) = lwi_data_64(i1:i2, ju1:j2)

!     ===============================================
      else if (((i2_gl - i1_gl)  == (L91I2 - L91I1)) .and.  &
     &         ((j2_gl - ju1_gl) == (L91J2 - L91J1))) then
!     ===============================================

        cmi_flags(i1:i2, ju1:j2) = lwi_data_91(i1:i2, ju1:j2)

!     ======================================
      else if (((i2_gl - i1_gl  + 1) == 72) .and.  &
     &         ((j2_gl - ju1_gl + 1) == 46)) then
!     ======================================

        do ij = ju1, j2
          do il = i1, i2

           ! if (ij < j2) then
            if (ij < j2_gl) then

              sumsq = lwi_data_91((il*2-1),(ij*2-1))**2 +  &
     &                lwi_data_91((il*2-1),(ij*2  ))**2 +  &
     &                lwi_data_91((il*2  ),(ij*2-1))**2 +  &
     &                lwi_data_91((il*2  ),(ij*2  ))**2

            else

              sumsq = lwi_data_91((il*2-1),(ij*2-1))**2 +  &
     &                lwi_data_91((il*2  ),(ij*2-1))**2 +  &
     &                lwi_data_91((il*2-1),(ij*2-1))**2 +  &
     &                lwi_data_91((il*2  ),(ij*2-1))**2

      if(il.eq.58.and.ij.eq.19) then
         write(6,*) 'in else: il, ij = ', il, ij

         write(6,*)'(il*2-1),(ij*2-1) = ',(il*2-1),(ij*2-1)
         write(6,*)'(il*2  ),(ij*2-1) = ',(il*2  ),(ij*2-1)
         write(6,*)'(il*2-1),(ij*2-1) = ',(il*2-1),(ij*2-1)
         write(6,*)'(il*2  ),(ij*2-1) = ',(il*2  ),(ij*2-1)

         write(6,*)'lwi_data_91((il*2-1),(ij*2-1)) = ',  &
     &              lwi_data_91((il*2-1),(ij*2-1))
         write(6,*)'lwi_data_91((il*2  ),(ij*2-1)) = ',  &
     &              lwi_data_91((il*2  ),(ij*2-1))
         write(6,*)'lwi_data_91((il*2-1),(ij*2-1)) = ',  &
     &              lwi_data_91((il*2-1),(ij*2-1))
         write(6,*)'lwi_data_91((il*2  ),(ij*2-1)) = ',  &
     &              lwi_data_91((il*2-1),(ij*2-1))
      endif

            end if

            select case(sumsq)
            case (4,12)
                cmi_flags(il,ij) = 1
            case (7,10,13,15,16,18,21,26)
                cmi_flags(il,ij) = 2
            case(20,23,28,31,36)
                cmi_flags(il,ij) = 3
            case default
                print *, 'Invalid lwi flags'
            end select

      if(il.eq.58.and.ij.eq.19) then
         write(6,*)'Get_cmi_flags: cmi_flags(58,19) = ',  &
     &             cmi_flags(58,19)
      endif

          end do
        end do

!     ====
      else
!     ====

        err_msg = 'CMI indexing problem in Get_cmi_flags.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, ju1_gl, j2_gl, 0, 0.0d0, 0.0d0)

      endif


      return

      end subroutine setCMIflags 
!EOC
!================================================================================
!BOP
!
! !IROUTINE: computeGlobCMIflags1
!
! !INTERFACE:
!
      subroutine computeGlobCMIflags1(gmiDomain, pr_diag, cmi_flags, cmi_flags1, &
     &                  i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl)
!
! !USES:
      use GmiGlob2Sub_mod          , only : global2Subdomain
      use GmiSub2Glob_mod          , only : subDomain2Global
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_rootProc, Get_map1_u
      use GmiDomainDecomposition_mod, only : Get_procID
      use GmiDomainDecomposition_mod, only : Get_numDomains
      use GmiDomainDecomposition_mod, only : Get_communicatorWorld, Get_iAmRootProc

      implicit none
!
! !INPUT PARAMETERS:
      logical           , intent(in) :: pr_diag
      integer           , intent(in) :: i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl
      type (t_gmiDomain), intent(in) :: gmiDomain
!
! !INPUT/OUTPUT PARAMETERS:
      integer          , intent(inOut) :: cmi_flags(i1:i2,ju1:j2)
      integer          , intent(inOut) :: cmi_flags1(i1:i2,ju1:j2)
!
! !DESCRIPTION:
! This routine is for the Dale's lightning parameterization.
! It is called by both the root and worker processors.
!
! !LOCAL VARIABLES:
      logical :: iAmRootProc
      integer :: commuWorld, procID, numDomains, rootProc
      integer, allocatable :: map1_u(:,:,:)
      integer, allocatable :: cmi_flags0_root (:,:)

      integer, allocatable :: glob_cmi_flags (:,:)
      real*8 , allocatable :: glob_cmiFlags0Real (:,:,:)
      integer, allocatable :: glob_cmi1Flags      (:,:)
      real*8 , allocatable :: glob_cmi_flags_real (:,:,:)
      real*8 , allocatable :: glob_cmi1FlagsReal  (:,:,:)

      integer, allocatable :: loc_cmi_flags (:,:)
      real*8 , allocatable :: loc_cmiFlags0Real (:,:,:)
      integer, allocatable :: loc_cmi1Flags      (:,:)
      real*8 , allocatable :: loc_cmi_flags_real (:,:,:)
      real*8 , allocatable :: loc_cmi1FlagsReal  (:,:,:)

      integer, parameter :: MSGnum1 = 3001
      integer, parameter :: MSGnum2 = 3002
      integer, parameter :: MSGnum3 = 3003
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID           (gmiDomain, procID)

      if (pr_diag) then
         Write(6,*) "computeGlobCMIflags1 called by", procID
      end if

      call Get_rootProc         (gmiDomain, rootProc)
      call Get_iAmRootProc      (gmiDomain, iAmRootProc)
      call Get_numDomains       (gmiDomain, numDomains)
      call Get_communicatorWorld(gmiDomain, commuWorld)

      allocate(map1_u(2,2,numDomains))
      call Get_map1_u(gmiDomain, map1_u)

      if (iAmRootProc) then
         ! allocate the real version of cmi_flags
         allocate (Glob_cmi_flags_real (i1_gl:i2_gl, ju1_gl:j2_gl, 1:1))
         Glob_cmi_flags_real (:,:,:) = real (0)

         allocate (Glob_cmi_flags (i1_gl:i2_gl, ju1_gl:j2_gl))
      end if

      ! allocate the real version of cmi_flags
      allocate (Loc_cmi_flags_real (i1:i2, ju1:j2, 1:1))
      Loc_cmi_flags_real (:,:,:) = real (0)

      allocate (loc_cmi_flags(i1:i2, ju1:j2))

      loc_cmi_flags      (:,:) = cmi_flags (:,:)
      loc_cmi_flags_real (:,:,1) = real (loc_cmi_flags(:,:))

      ! all PEs should call this routine
      ! worker PEs send their copy of cmi_flagsGlob to the root PE

!     ========================
      call subDomain2Global  &
!     ========================  
     &    (Glob_cmi_flags_real, Loc_cmi_flags_real, i1_gl, i2_gl, ju1_gl, &
     &     j2_gl, i1, i2, ju1, j2,  1, 1, &
     &     rootProc, procID, map1_u, numDomains, MSGnum1, commuWorld)
         
      if (iAmRootProc) then
         glob_cmi_flags (:,:) = int (Glob_cmi_flags_real (:,:,1))

         allocate (glob_cmi1FlagsReal (i1_gl:i2_gl, ju1_gl:j2_gl, 1:1))
         allocate (glob_cmi1Flags (i1_gl:i2_gl, ju1_gl:j2_gl))
         allocate (glob_cmiFlags0Real (i1_gl:i2_gl, ju1_gl:j2_gl, 1:1))
         glob_cmi1FlagsReal (:,:,:) = real (0)
         glob_cmi1Flags     (:,:) = int (0)
         glob_cmiFlags0Real (:,:,:) = real (0)

         allocate (cmi_flags0_root (i1_gl:i2_gl,ju1_gl:j2_gl))
         cmi_flags0_root = calcCMI1flags (glob_cmi_flags  , i1_gl, i2_gl, ju1_gl, j2_gl)
         glob_cmi1Flags  = calcCMI1flags (cmi_flags0_root, i1_gl, i2_gl, ju1_gl, j2_gl)
         glob_cmi1FlagsReal (:,:,1) = real (glob_cmi1Flags (:,:))
         glob_cmiFlags0Real (:,:,1) = real (cmi_flags0_root (:,:))
      end if

      allocate (loc_cmi1FlagsReal (i1:i2, ju1:j2, 1:1))
      allocate (loc_cmi1Flags (i1:i2, ju1:j2))
      allocate (loc_cmiFlags0Real (i1:i2, ju1:j2, 1:1))
      loc_cmi1FlagsReal (:,:,:) = real (0)
      loc_cmi1Flags     (:,:) = int (0)
      loc_cmiFlags0Real (:,:,:) = real (0)
      
      ! the master distributes cmi1Flags to the workers
      call global2Subdomain (glob_cmi1FlagsReal, loc_cmi1FlagsReal, i1_gl, &
           & i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, 1, 1,  &
           &   rootProc, procID, map1_u, numDomains, commuWorld, MSGnum2)
      
      ! the master distributes cmi_flags0 to the workers
      call global2Subdomain (glob_cmiFlags0Real, loc_cmiFlags0Real, i1_gl, &
           & i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, 1, 1,  &
              &   rootProc, procID, map1_u, numDomains, commuWorld, MSGnum3)
      
      ! Worker processors set the cmi1Flags and cmi_flags (0)
      ! in the Emission derived type

      cmi_flags  (:,:) = int (loc_cmiFlags0Real (:,:,1))
      cmi_flags1 (:,:) = int (loc_cmi1FlagsReal (:,:,1))
      
      if (iAmRootProc) then
         ! the master deallocates it's copy of cmi_flags
         deallocate (cmi_flags0_root)
         deallocate (glob_cmi1Flags)
         deallocate (glob_cmi1FlagsReal)
         deallocate (glob_cmi_flags_real)
         deallocate (glob_cmiFlags0Real)
         deallocate (glob_cmi_flags)
      endif
      deallocate (map1_u)
      deallocate (loc_cmi1Flags)
      deallocate (loc_cmi1FlagsReal)
      deallocate (loc_cmi_flags_real)
      deallocate (loc_cmiFlags0Real)
      deallocate (loc_cmi_flags)

      return

      end subroutine computeGlobCMIflags1
!EOC
!-------------------------------------------------------------------------
!BOP

  function calcCmi1Flags (cmiFlags, i1, i2, ju1, j2) result (cmi1Flags)

    ! arguments
    integer, intent (in) :: cmiFlags (i1:i2, ju1:j2)
    integer, intent (in) :: i1, i2, ju1, j2

    ! local variables
    ! adjusted coastal-marine-ice flags
    integer :: cmi1Flags (i1:i2, ju1:j2)
!
! !DESCRIPTION:
!  This function should only be called if running on one processor
!  or by the root processor.

    ! dx1,dx2,dy1, and dy2 are work arrays used in
    ! determining if grid points are "coastal".
    integer :: dx1 (i1:i2, ju1:j2), dx2(i1:i2,ju1:j2)
    integer :: dy1 (i1:i2, ju1:j2), dy2(i1:i2,ju1:j2)

    ! initializations
    cmi1Flags (:,:) = 0
    dx1 (:,:) = 0
    dx2 (:,:) = 0
    dy1 (:,:) = 0
    dy2 (:,:) = 0

    ! start with original coastal-marine-ice values
    ! but set all ice values to marine
    cmi1Flags = cmiFlags
    where (cmi1Flags == 3) cmi1Flags = 1

    ! get the differences between adjacent grid boxes
    dx1 = cshift (cmi1Flags, 1,  1) - cmi1Flags
    dx2 = cshift (cmi1Flags, -1, 1) - cmi1Flags
    dy1 = cshift (cmi1Flags, 1,  2) - cmi1Flags
    dy2 = cshift (cmi1Flags, -1, 2) - cmi1Flags

    ! find where there is a difference and set the
    ! cell to a coastal one
    where ((dx1 .ne. 0) .or. (dx2 .ne. 0) .or. (dy1 .ne. 0) .or. (dy2 .ne. 0)) &
         cmi1Flags = 2

    ! set all cells that haven't been adjusted to coastal back ice
    where ((cmi1Flags .ne. 2) .and. (cmiFlags .eq. 3)) &
         cmi1Flags = 3

  end function calcCmi1Flags

!EOC
!================================================================================
      end module GmiLandWaterIce_mod
