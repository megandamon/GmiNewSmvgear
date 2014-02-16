!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiWrapMaster_mod
!
! !INTERFACE:
!
   module GmiWrapMaster_mod
!
! !USES:
!
   implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
   private
   public  wrapMaster_2d
   public  wrapMaster_3du
   public  wrapMaster_3dv
   public  wrapMaster_4d
!
! !DESCRIPTION:
!  Provides subroutines to wrap 2D/3D/4D ghost zone values for the
!  master processor. They operate on an entire array. The subroutines
!  are used in the GMI code.
!
! !AUTHOR: 
!  Jules Kouatchou, NASA/GSFC, Jules.Kouatchou-1@nasa.gov
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: wrapMaster_2d
!
! !INTERFACE:
!
      subroutine wrapMaster_2d (var_2d, i1, i2, ju1, j2, ilo, ihi, &
                                julo, jhi, gmi_nborder)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in)   :: i1, i2, ju1, j2, ilo, ihi, julo, jhi
      integer, intent(in)   :: gmi_nborder
!
! !INPUT/OUTPUT PARAMETERS:
!!    var_2d : the specified 2D array
      real*8, intent(inout) :: var_2d(ilo:ihi, julo:jhi)
!
! !DESCRIPTION:
!  Wraps the 2D ghost zone values for the master processor only.
!  It operates on an entire array.
!
! !LOCAL VARIABLES:
      integer :: il
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      do il = 1, gmi_nborder
         var_2d(i1-il,ju1:j2) = var_2d(i2-il+1,ju1:j2)
         var_2d(i2+il,ju1:j2) = var_2d(i1+il-1,ju1:j2)
      end do

      return

      end subroutine wrapMaster_2d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: wrapMaster_3du
!
! !INTERFACE:
!
      subroutine wrapMaster_3du (var_3du, i1, i2, ju1, j2, k1, k2,    &
                                   ilo, ihi, julo, jhi, gmi_nborder)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in)   :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in)   :: gmi_nborder
!
! !INPUT/OUTPUT PARAMETERS:
!!    var_3du : the specified 3D/julo array
      real*8, intent(inout)  :: var_3du(ilo:ihi, julo:jhi, k1:k2)
!
! !DESCRIPTION:
!  Wraps the 3D/julo ghost zone values for the master processor only.
!  It operates on an entire array.
!
! !LOCAL VARIABLES:
      integer :: il
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      do il = 1, gmi_nborder
         var_3du(i1-il,ju1:j2,:) = var_3du(i2-il+1,ju1:j2,:)
         var_3du(i2+il,ju1:j2,:) = var_3du(i1+il-1,ju1:j2,:)
      end do
 
      return

      end subroutine wrapMaster_3du
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: wrapMaster_3dv
!
! !INTERFACE:
!
      subroutine wrapMaster_3dv (var_3dv, i1, i2, jv1, j2, k1, k2,    &
                                   ilo, ihi, jvlo, jhi, gmi_nborder)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in)   :: i1, i2, jv1, j2, k1, k2, ilo, ihi, jvlo, jhi
      integer, intent(in)   :: gmi_nborder
!
! !INPUT/OUTPUT PARAMETERS:
!!    var_3dv : the specified 3D/jvlo array
      real*8, intent(inout)  :: var_3dv(ilo:ihi, jvlo:jhi, k1:k2)
!
! !DESCRIPTION:
!  Wraps the 3D/jvlo ghost zone values for the master processor only.
!  It operates on an entire array.
!
! !LOCAL VARIABLES:
      integer :: il
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      do il = 1, gmi_nborder
         var_3dv(i1-il,jv1:j2,:) = var_3dv(i2-il+1,jv1:j2,:)
         var_3dv(i2+il,jv1:j2,:) = var_3dv(i1+il-1,jv1:j2,:)
      end do
 
      return

      end subroutine wrapMaster_3dv
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: wrapMaster_4d
!
! !INTERFACE:
!
      subroutine wrapMaster_4d (var_4d, num_d4, i1, i2, ju1, j2, k1, k2, &
                                        ilo, ihi, julo, jhi, gmi_nborder)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in)  :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in)  :: gmi_nborder
      integer, intent(in)  :: num_d4 ! number of elements in 4th dimension
!
! !INPUT/OUTPUT PARAMETERS:
!!    var_4d : the specified 4D array
      real*8,  intent(inout)  :: var_4d(ilo:ihi, julo:jhi, k1:k2, num_d4)
!
! !DESCRIPTION:
!  Wraps the 4D ghost zone values for the master processor only.
!  It operates on an entire array.
!
! !LOCAL VARIABLES:
      integer :: il
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      do il = 1, gmi_nborder
         var_4d(i1-il,ju1:j2,:,:) = var_4d(i2-il+1,ju1:j2,:,:)
         var_4d(i2+il,ju1:j2,:,:) = var_4d(i1+il-1,ju1:j2,:,:)
      end do
 
      return

      end subroutine wrapMaster_4d
!EOC
!-------------------------------------------------------------------------
end module GmiWrapMaster_mod
