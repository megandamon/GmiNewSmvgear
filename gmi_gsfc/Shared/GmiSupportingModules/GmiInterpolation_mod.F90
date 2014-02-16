    module GmiInterpolation_mod

    implicit none

!=============================================================================
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original Interp* code from Ricky Rood, DAO)
!   jrt@llnl.gov
!
    private
    public :: Interp
    public :: Interp_Linear
    public :: Interp_Bilinear
    public :: Interp_Trilinear
!
!=============================================================================

      CONTAINS

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Interp
!
! DESCRIPTION
!   Linear interpolation is performed on an array of input values:
!     xin(1) < xout(j) < xin(n) or
!     xin(1) > xout(j) > xin(n)
!
! ARGUMENTS
!   xin  : input levels
!   yin  : input values
!   nn   : number of input levels
!   xout : output levels
!   yout : output values
!   mm   : number of output levels
!
!-----------------------------------------------------------------------------

      subroutine Interp  &
     &  (xin, yin, nn, xout, yout, mm)

      implicit none


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: nn
      integer :: mm
      real*8  :: xin (nn)
      real*8  :: yin (nn)
      real*8  :: xout(mm)
      real*8  :: yout(mm)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ij, jj
      integer :: start


!     ----------------
!     Begin execution.
!     ----------------

      start = 2


!     =====================
      if (xin(1) < xin(nn)) then  ! increasing
!     =====================

        do ij = 1, mm

          if (xout(ij) <= xin(1)) then

            yout(ij) = yin(1)

          else if (xout(ij) >= xin(nn)) then

            yout(ij) = yin(nn)

          else

            jjloop1: do jj = start, nn

              if (xout(ij) <= xin(jj)) then

                yout(ij) =  &
     &            yin(jj-1) +  &
     &            ((yin(jj) - yin(jj-1)) / (xin(jj) - xin(jj-1)) *  &
     &             (xout(ij) - xin(jj-1)))

                start = jj

!               ============
                exit jjloop1
!               ============

              end if

            end do jjloop1

          end if

        end do


!     ==========================
      else if (xin(1) > xin(nn)) then  ! decreasing
!     ==========================

        do ij = 1, mm

          if (xout(ij) >= xin(1)) then

            yout(ij) = yin(1)

          else if (xout(ij) <= xin(nn)) then

            yout(ij) = yin(nn)

          else

            jjloop2: do jj = start, nn

              if (xout(ij) >= xin(jj)) then

                yout(ij) =  &
     &            yin(jj-1) +  &
     &            ((yin(jj) - yin(jj-1)) / (xin(jj) - xin(jj-1)) *  &
     &             (xout(ij) - xin(jj-1)))

                start = jj

!               ============
                exit jjloop2
!               ============

              end if

            end do jjloop2

          end if

        end do

!     ====
      else
!     ====

        yout(:) = xin(1)

      end if


      return

      end subroutine Interp


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Interp_Linear
!
! DESCRIPTION
!   This routine does a simple linear interpolation (see equation below).
!
! ARGUMENTS
!   a1 : value of function at point 1 (input)
!   ax : value of function at point 3 (output)
!   a2 : value of function at point 2 (input)
!   b1 : coordinate value  at point 1 (input)
!   b2 : coordinate value  at point 2 (input)
!   b3 : coordinate value  at point 3 (input)
!
!-----------------------------------------------------------------------------

      subroutine Interp_Linear  &
     &  (a1, ax, a2, b1, b2, b3)

      implicit none


!     ----------------------
!     Argument declarations.
!     ----------------------

      real*8  :: a1, ax, a2
      real*8  :: b1, b2, b3


!     ----------------
!     Begin execution.
!     ----------------

      ax = a1 + (((b2 - b1) / (b3 - b1)) * (a2 - a1))

      return

      end subroutine Interp_Linear

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Interp_Bilinear
!
! DESCRIPTION
!   This routine does a simple bi-linear interpolation from a 2D table
!   with monotonically changing coordinates.  The coordinate system is
!   assumed to be Cartesian.
!   From "Numerical Recipes In Fortran, 2nd ed." by William H. Press, et al.
!
! ARGUMENTS
!   a1      : value of coordinate 1 (input)
!   a2      : value of coordinate 2 (input)
!   vinterp : value of function at coordinate(a1,a2) (output)
!   v1      : coordinate values of dimension 1 (input)
!   v2      : coordinate values of dimension 2 (input)
!   l1      : length of v1 (input)
!   l2      : length of v2 (input)
!   data    : 2D table of data for interpolation. (input)
!-----------------------------------------------------------------------------

      subroutine Interp_Bilinear  &
     &  (a1, a2, vinterp, v1, v2, l1, l2, data)

      implicit none


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: l1, l2

      real*8  :: a1, a2
      real*8  :: vinterp
      real*8  :: v1(l1), v2(l2)
      real*8  :: data(l1, l2)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: i1, i2

      real*8  :: t, u


!     ----------------
!     Begin execution.
!     ----------------


!     ----------------------------------------------------
!     Find the indices in the table where v1>a1 and v2>a2.
!     ----------------------------------------------------
      i1 = Minval (Minloc (abs(v1-a1)))
      if(a1.ge.v1(i1).and.a1.le.v1(i1+1)) i1 = i1+1
      if(a1.le.v1(i1).and.a1.ge.v1(i1+1)) i1 = i1+1
      i1 = max(min(i1,size(v1)),2)

      i2 = Minval (Minloc (abs(v2-a2)))
      if(a2.ge.v2(i2).and.a2.le.v2(i2+1)) i2 = i2+1
      if(a2.le.v2(i2).and.a2.ge.v2(i2+1)) i2 = i2+1
      i2 = max(min(i2,size(v2)),2)

!     ----------------------------------------------------------------
!     Interpolate the data to find the value vinterp at point (a1,a2).
!     ----------------------------------------------------------------

      t = (a1 - v1(i1-1)) / (v1(i1) - v1(i1-1))
      u = (a2 - v2(i2-1)) / (v2(i2) - v2(i2-1))

      t = max(min(t,1.0d0),0.0d0)
      u = max(min(u,1.0d0),0.0d0)

      vinterp = (1-t) * (1-u) * data(i1-1,i2-1) +  &
     &             t  * (1-u) * data(i1  ,i2-1) +  &
     &             t  *    u  * data(i1  ,i2  ) +  &
     &          (1-t) *    u  * data(i1-1,i2  )

      return

      end subroutine Interp_Bilinear

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Interp_Trilinear
!
! DESCRIPTION
!   This routine does a simple tri-linear interpolation from a 3D table
!   with monotonically changing coordinates.  The coordinate system is
!   assumed to be Cartesian.
!   From "Numerical Recipes In Fortran, 2nd ed." by William H. Press, et al.
!
! ARGUMENTS
!   a1      : value of coordinate 1 (input)
!   a2      : value of coordinate 2 (input)
!   a3      : value of coordinate 3 (input)
!   vinterp : value of function at coordinate(a1,a2,a3) (output)
!   v1      : coordinate values of dimension 1 (input)
!   v2      : coordinate values of dimension 2 (input)
!   v3      : coordinate values of dimension 3 (input)
!   l1      : length of v1 (input)
!   l2      : length of v2 (input)
!   l3      : length of v3 (input)
!   data    : 3D table of data for interpolation. (input)
!-----------------------------------------------------------------------------

   subroutine Interp_Trilinear  &
     (a1, a2, a3, vinterp, v1, v2, v3, l1, l2, l3, data)


   implicit none


!     ----------------------
!     Argument declarations.
!     ----------------------

   integer :: l1, l2, l3

   real*8  :: a1, a2, a3
   real*8  :: vinterp
   real*8  :: v1(l1), v2(l2), v3(l3)
   real*8  :: data(l1, l2, l3)


!     ----------------------
!     Variable declarations.
!     ----------------------

   integer :: i1, i2, i3

   real*8  :: t, u, v


!     ----------------
!     Begin execution.
!     ----------------


!     ----------------------------------------------------
!     Find the indices in the table where v1>a1 and v2>a2.
!     ----------------------------------------------------
   i1 = Minval (Minloc (abs(v1-a1)))
   i1 = max(min(i1,size(v1)-1),2)
   if(v1(1).lt.v1(l1)) then
      if(a1.ge.v1(i1).and.a1.le.v1(i1+1)) i1 = i1+1
    else
      if(a1.le.v1(i1).and.a1.ge.v1(i1+1)) i1 = i1+1
    endif
   i1 = max(min(i1,size(v1)),2)

   i2 = Minval (Minloc (abs(v2-a2)))
   i2 = max(min(i2,size(v2)-1),2)
   if(v2(1).lt.v2(l2)) then
      if(a2.ge.v2(i2).and.a2.le.v2(i2+1)) i2 = i2+1
    else
      if(a2.le.v2(i2).and.a2.ge.v2(i2+1)) i2 = i2+1
    endif
   i2 = max(min(i2,size(v2)),2)

   i3 = Minval (Minloc (abs(v3-a3)))
   i3 = max(min(i3,size(v3)-1),2)
   if(v3(1).lt.v3(l3)) then
      if((a3.ge.v3(i3).and.a3.le.v3(i3+1))) i3 = i3+1
    else
      if((a3.le.v3(i3).and.a3.ge.v3(i3+1))) i3 = i3+1
    endif
   i3 = max(min(i3,size(v3)),2)

!     ----------------------------------------------------------------
!     Interpolate the data to find the value vinterp at point (a1,a2).
!     ----------------------------------------------------------------

   t = (a1 - v1(i1-1)) / (v1(i1) - v1(i1-1))
   u = (a2 - v2(i2-1)) / (v2(i2) - v2(i2-1))
   v = (a3 - v3(i3-1)) / (v3(i3) - v3(i3-1))

   t = max(min(t,1.0d0),0.0d0)
   u = max(min(u,1.0d0),0.0d0)
   v = max(min(v,1.0d0),0.0d0)

   if(i1-1.lt.1 .or. i1.gt.l1) print *,'I1',i1,l1
   if(i2-1.lt.1 .or. i2.gt.l2) print *,'I2',i2,l2
   if(i3-1.lt.1 .or. i3.gt.l3) print *,'I3',i3,l3,size(v3)


   vinterp = (1-v) * (1-u) * (1-t) * data(i1-1,i2-1,i3-1) +  &
             (1-v) * (1-u) *    t  * data(i1  ,i2-1,i3-1) +  &
             (1-v) *    u  * (1-t) * data(i1-1,i2  ,i3-1) +  &
             (1-v) *    u  *    t  * data(i1  ,i2  ,i3-1) +  &
                v  * (1-u) * (1-t) * data(i1-1,i2-1,i3  ) +  &
                v  * (1-u) *    t  * data(i1  ,i2-1,i3  ) +  &
                v  *    u  * (1-t) * data(i1-1,i2  ,i3  ) +  &
                v  *    u  *    t  * data(i1  ,i2  ,i3  )

   return

   end subroutine Interp_Trilinear

!-------------------------------------------------------------------------------

    end module GmiInterpolation_mod
