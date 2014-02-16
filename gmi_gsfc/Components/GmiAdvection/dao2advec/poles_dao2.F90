!
!=============================================================================
!
! $Id: poles_dao2.F90,v 1.10 2013-07-31 15:00:47 ssteenro Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Shian-Jiann Lin, DAO)
!   jrt@llnl.gov
!
! FILE
!   poles_dao2.F
!
! ROUTINES
!   Do_Cross_Terms_Pole_I2d2
!   Do_Fyppm_Pole_I2d2
!   Do_Yadv_Pole_I2d2
!   Do_Ymist_Pole1_I2d2
!   Do_Ymist_Pole2_I2d2
!   Do_Divergence_Pole_Sum
!   Do_Yadv_Pole_Sum
!   Do_Ytp_Pole_Sum
!
!=============================================================================
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Cross_Terms_Pole_I2d2
!
! DESCRIPTION
!   This routine sets "va" at the Poles.
!
! ARGUMENTS
!   cry : Courant number in N-S direction
!   va  : average of Courant numbers from ij and ij+1
!
!-----------------------------------------------------------------------------
!
      subroutine Do_Cross_Terms_Pole_I2d2  &
     &  (cry, va, &
     &   numLonDomains, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   mapi_all, commu_npole, commu_spole, numDomains, ivert)
!
      use GmiGather_mod, only : Gmi_Pole_Allgather
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numLonDomains, numDomains
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, j1p
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ivert
      integer, intent(in) :: commu_spole, commu_npole
      integer, intent(in) :: mapi_all(2,numDomains)
!
      real*8  :: cry(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: va (ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: i2d2
      integer :: il, ik
!
      real*8  :: crysp_gl(i1_gl:i2_gl, ju1+1:ju1+1, k1:k2)
      real*8  :: crynp_gl(i1_gl:i2_gl,    j2:j2,    k1:k2)
!
      real*8  :: vasp_gl (i1_gl:i2_gl, ju1+1:ju1+1, k1:k2)
      real*8  :: vanp_gl (i1_gl:i2_gl,  j2:j2,  k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      i2d2 = i2_gl / 2
!
!
!     ====================
      if (j1p == ju1_gl+1) then
!     ====================
!
!       ---------------------------------------------
!       Polar Cap NOT Enlarged:
!       Get cross terms for N-S horizontal advection.
!       ---------------------------------------------
!
!       ==================
        if (ju1 == ju1_gl) then
!       ==================
!
          if (numLonDomains == 1) then
!
            do il = i1, i2d2
!
              va(il,ju1+1,:) = 0.5d0 * (cry(il,ju1+1,:) - cry(il+i2d2,ju1+1,:))
!
              va(il+i2d2,ju1+1,:) = -va(il,ju1+1,:)
!
            end do
!
          else
!
            crysp_gl(:,:,:) = 0.0d0
            vasp_gl (:,:,:) = 0.0d0
!
!           ==========================
            call Gmi_Pole_Allgather  &
!           ==========================
     &        (ju1+1, ju1+1, cry, crysp_gl, &
     &  mapi_all, commu_spole, commu_npole, numDomains, numLonDomains, &
     &  i1, i2, ju1, j2, k1, k2, ivert, ilo, ihi, julo, jhi, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl)
!
            do ik = k1, k2
              do il = i1_gl, i2d2
!
                vasp_gl(il,ju1+1,ik) =  &
     &            0.5d0 * (crysp_gl(il,ju1+1,ik) - crysp_gl(il+i2d2,ju1+1,ik))
!
                vasp_gl(il+i2d2,ju1+1,ik) = -vasp_gl(il,ju1+1,ik)
!
              end do
            end do
!
            va(i1:i2,ju1+1,:) = vasp_gl(i1:i2,ju1+1,:)
!
          end if
!
!       ======
        end if
!       ======
!
!
!       ================
        if (j2 == j2_gl) then
!       ================
!
          if (numLonDomains == 1) then
!
            do il = i1, i2d2
!
              va(il,j2,:) =  &
     &          0.5d0 * (cry(il,j2,:) - cry(il+i2d2,j2,:))
!
              va(il+i2d2,j2,:) = -va(il,j2,:)
!
            end do
!
          else
!
            crynp_gl(:,:,:) = 0.0d0
            vanp_gl (:,:,:) = 0.0d0
!
!           ==========================
            call Gmi_Pole_Allgather  &
!           ==========================
     &        (j2, j2, cry, crynp_gl, &
     &  mapi_all, commu_spole, commu_npole, numDomains, numLonDomains, &
     &  i1, i2, ju1, j2, k1, k2, ivert, ilo, ihi, julo, jhi, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl)
!
            do ik = k1, k2
              do il = i1_gl, i2d2
!
                vanp_gl(il,j2,ik) =  &
     &            0.5d0 * (crynp_gl(il,j2,ik) - crynp_gl(il+i2d2,j2,ik))
!
                vanp_gl(il+i2d2,j2,ik) = -vanp_gl(il,j2,ik)
!
              end do
            end do
!
            va(i1:i2,j2,:) = vanp_gl(i1:i2,j2,:)
!
          end if
!
!       ======
        end if
!       ======
!
!     ======
      end if
!     ======
!
!
      return
!
      end
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Fyppm_Pole_I2d2
!
! DESCRIPTION
!   This routine sets "al" & "ar" at the Poles.
!
! ARGUMENTS
!   al : left  edge value of the test parabola
!   ar : right edge value of the test parabola
!
!-----------------------------------------------------------------------------
!
      subroutine Do_Fyppm_Pole_I2d2  &
     &  (al, ar, &
     &   numLonDomains, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   mapi_all, commu_npole, commu_spole, numDomains, ivert)
!
      use GmiGather_mod, only : Gmi_Pole_Allgather
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numLonDomains, numDomains
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ivert
      integer, intent(in) :: commu_npole, commu_spole
      integer, intent(in) :: mapi_all(2, numDomains)
!
      real*8  :: al(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: ar(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: i2d2
      integer :: il, ik, ilos
!
      real*8  :: alsp_gl(i1_gl:i2_gl, ju1_gl+1:ju1_gl+1, k1:k2)
      real*8  :: arnp_gl(i1_gl:i2_gl,  j2_gl-1:j2_gl-1,  k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      i2d2 = i2_gl / 2
!
!
!     ==================
      if (ju1 == ju1_gl) then
!     ==================
!
        if (numLonDomains == 1) then
!
          do il = i1, i2d2
            al(il,     ju1_gl,:) = al(il+i2d2,ju1_gl+1,:)
            al(il+i2d2,ju1_gl,:) = al(il,     ju1_gl+1,:)
          end do
!
        else
!
          alsp_gl(:,:,:) = 0.0d0
!
!
!         ==========================
          call Gmi_Pole_Allgather  &
!         ==========================
     &      (ju1_gl+1, ju1_gl+1, al, alsp_gl, &
     &  mapi_all, commu_spole, commu_npole, numDomains, numLonDomains, &
     &  i1, i2, ju1, j2, k1, k2, ivert, ilo, ihi, julo, jhi, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl)
!
!          do il = i1_gl, i2d2
!            alsp_gl(il     ,ju1,:) = alsp_gl(il+i2d2,ju1+1,:)
!            alsp_gl(il+i2d2,ju1,:) = alsp_gl(il     ,ju1+1,:)
!          end do
!
          do ik = k1, k2
            do il = i1, i2
              ilos = il+i2d2
              if(ilos.gt.i2_gl) ilos = ilos-(i2_gl-i1_gl+1)
              al(il,ju1_gl,ik) = alsp_gl(ilos,ju1_gl+1,ik)
            end do
          end do
!
        end if
!
!     ======
      end if
!     ======
!
!
!     ================
      if (j2 == j2_gl) then
!     ================
!
        if (numLonDomains == 1) then
!
          do il = i1, i2d2
            ar(il,     j2,:)  = ar(il+i2d2,j2-1,:)
            ar(il+i2d2,j2,:)  = ar(il,     j2-1,:)
          end do
!
        else
!
          arnp_gl(:,:,:) = 0.0d0
!
!         ==========================
          call Gmi_Pole_Allgather  &
!         ==========================
     &      (j2_gl-1, j2_gl-1, ar, arnp_gl, &
     &  mapi_all, commu_spole, commu_npole, numDomains, numLonDomains, &
     &  i1, i2, ju1, j2, k1, k2, ivert, ilo, ihi, julo, jhi, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl)
!
!          do il = i1_gl, i2d2
!            arnp_gl(il     ,j2,:) = arnp_gl(il+i2d2,j2-1,:)
!            arnp_gl(il+i2d2,j2,:) = arnp_gl(il     ,j2-1,:)
!          end do
!
          do ik = k1, k2
            do il = i1, i2
              ilos = il+i2d2
              if(ilos.gt.i2_gl) ilos = ilos-(i2_gl-i1_gl+1)
              ar(il,j2,ik) = arnp_gl(ilos,j2-1,ik)
            end do
          end do
!
        end if
!
!     ======
      end if
!     ======
!
!
      return
!
      end
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Yadv_Pole_I2d2
!
! DESCRIPTION
!   This routine sets "qquwk" at the Poles.
!
! ARGUMENTS
!   qqu   : concentration contribution from E-W advection (mixing ratio)
!   qquwk : qqu working array (mixing ratio)
!
!-----------------------------------------------------------------------------
!
      subroutine Do_Yadv_Pole_I2d2 (qqu, qquwk, &
     &   numLonDomains, gmi_nborder, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   mapi_all, commu_npole, commu_spole, numDomains, ivert)
!
      use GmiGather_mod, only : Gmi_Pole_Allgather
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numLonDomains, gmi_nborder, numDomains
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, j1p
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ivert
      integer, intent(in) :: commu_npole, commu_spole
      integer, intent(in) :: mapi_all(2, numDomains)
!
      real*8, intent(in)  :: qqu  (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: qquwk(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: i2d2
      integer :: il, ij, ik
      integer :: ina, inb, inaa
!
      real*8  :: qqusp_gl  (i1_gl:i2_gl, ju1+1:ju1+gmi_nborder, k1:k2)
      real*8  :: qqunp_gl  (i1_gl:i2_gl, j2-gmi_nborder:j2-1,   k1:k2)
!
      real*8  :: qquwksp_gl(i1_gl:i2_gl, julo:ju1-1, k1:k2)
      real*8  :: qquwknp_gl(i1_gl:i2_gl, j2+1:jhi,   k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      i2d2 = i2_gl / 2
!
!
!     ====================
      if (j1p == ju1_gl+1) then
!     ====================
!
!       -----------------------
!       Polar Cap NOT Enlarged.
!       -----------------------
!
!       ==================
        if (ju1 == ju1_gl) then
!       ==================
!
          if (numLonDomains == 1) then
!
            do il = i1, i2d2
              do inb = 1, gmi_nborder
!
                qquwk(il,     ju1-inb,:) = qqu(il+i2d2,ju1+inb,:)
                qquwk(il+i2d2,ju1-inb,:) = qqu(il,     ju1+inb,:)
!
              end do
            end do
!
          else
!
            qqusp_gl  (:,:,:) = 0.0d0
            qquwksp_gl(:,:,:) = 0.0d0
!
!           ==========================
            call Gmi_Pole_Allgather  &
!           ==========================
     &        (ju1_gl+1, ju1_gl+gmi_nborder, qqu, qqusp_gl, &
     &  mapi_all, commu_spole, commu_npole, numDomains, numLonDomains, &
     &  i1, i2, ju1, j2, k1, k2, ivert, ilo, ihi, julo, jhi, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl)
!
            ij = ju1_gl
!
            do ik = k1, k2
              do il = i1_gl, i2d2
                do inb = 1, gmi_nborder
!
                  qquwksp_gl(il     ,ij-inb,ik) = qqusp_gl(il+i2d2,ij+inb,ik)
                  qquwksp_gl(il+i2d2,ij-inb,ik) = qqusp_gl(il     ,ij+inb,ik)
!
                end do
              end do
            end do
!
            do inb = 1, gmi_nborder
              qquwk(i1:i2,ij-inb,:) = qquwksp_gl(i1:i2,ij-inb,:)
            end do
!
          end if
!
!       ======
        end if
!       ======
!
!
!       ================
        if (j2 == j2_gl) then
!       ================
!
          if (numLonDomains == 1) then
!
            do il = i1, i2d2
              do inb = 1, gmi_nborder
!
                qquwk(il,     j2+inb,:) = qqu(il+i2d2,j2-inb,:)
                qquwk(il+i2d2,j2+inb,:) = qqu(il,     j2-inb,:)
!
              end do
            end do
!
          else
!
            qqunp_gl  (:,:,:) = 0.0d0
            qquwknp_gl(:,:,:) = 0.0d0
!
!           ==========================
            call Gmi_Pole_Allgather  &
!           ==========================
     &        (j2_gl-gmi_nborder, j2_gl-1, qqu, qqunp_gl, &
     &  mapi_all, commu_spole, commu_npole, numDomains, numLonDomains, &
     &  i1, i2, ju1, j2, k1, k2, ivert, ilo, ihi, julo, jhi, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl)
!
            ij = j2_gl
!
            do ik = k1, k2
              do il = i1_gl, i2d2
                do inb = 1, gmi_nborder
!
                  qquwknp_gl(il     ,ij+inb,ik) = qqunp_gl(il+i2d2,ij-inb,ik)
                  qquwknp_gl(il+i2d2,ij+inb,ik) = qqunp_gl(il     ,ij-inb,ik)
!
                end do
              end do
            end do
!
            do inb = 1, gmi_nborder
              qquwk(i1:i2,ij+inb,:) = qquwknp_gl(i1:i2,ij+inb,:)
            end do
!
          end if
!
!       ======
        end if
!       ======
!
!     ======
      end if
!     ======
!
!
      return
!
      end
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Ymist_Pole1_I2d2
!
! DESCRIPTION
!   This routine sets "dcy" at the Poles.
!
! ARGUMENTS
!   dcy : slope of concentration distribution in N-S direction (mixing ratio)
!   qqu : concentration contribution from E-W advection (mixing ratio)
!
!-----------------------------------------------------------------------------
!
      subroutine Do_Ymist_Pole1_I2d2 (dcy, qqu, &
     &   numLonDomains, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   mapi_all, commu_npole, commu_spole, numDomains, ivert)
!
      use GmiGather_mod, only : Gmi_Pole_Allgather
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numLonDomains, numDomains
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ivert
      integer, intent(in) :: commu_npole, commu_spole
      integer, intent(in) :: mapi_all(2, numDomains)
!
      real*8, intent(out) :: dcy(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: qqu(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: i2d2
      integer :: il, ik
!
      real*8  :: pmax, pmin
      real*8  :: r24
      real*8  :: tmp
!
      real*8  :: dcysp_gl(i1_gl:i2_gl, ju1+1:ju1+1, k1:k2)
      real*8  :: dcynp_gl(i1_gl:i2_gl,  j2-1:j2-1,  k1:k2)
!
      real*8  :: qqusp_gl(i1_gl:i2_gl, ju1:ju1+3, k1:k2)
      real*8  :: qqunp_gl(i1_gl:i2_gl,  j2-3:j2,  k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      i2d2 = i2_gl / 2
!
      r24  = 1.0d0 / 24.0d0
!
!
!     ==================
      if (ju1 == ju1_gl) then
!     ==================
!
!       ==================
        if (numLonDomains == 1) then
!       ==================
!
          do ik = k1, k2
            do il = i1, i2d2
!
              tmp  =  &
     &          ((8.0d0 * (qqu(il,ju1+2,ik) - qqu(il,ju1,ik))) +  &
     &           qqu(il+i2d2,ju1+1,ik) - qqu(il,ju1+3,ik)) *  &
     &          r24
!
              pmax = Max (qqu(il,ju1,ik), qqu(il,ju1+1,ik),  &
     &                    qqu(il,ju1+2,ik)) -  &
     &               qqu(il,ju1+1,ik)
!
              pmin = qqu(il,ju1+1,ik) -  &
     &               Min (qqu(il,ju1,ik), qqu(il,ju1+1,ik),  &
     &                    qqu(il,ju1+2,ik))
!
              dcy(il,ju1+1,ik) =  &
     &          Sign (Min (Abs (tmp), pmin, pmax), tmp)
!
            end do
          end do
!
          do ik = k1, k2
            do il = i1 + i2d2, i2
!
              tmp  =  &
     &          ((8.0d0 * (qqu(il,ju1+2,ik) - qqu(il,ju1,ik))) +  &
     &           qqu(il-i2d2,ju1+1,ik) - qqu(il,ju1+3,ik)) *  &
     &          r24
!
              pmax = Max (qqu(il,ju1,ik), qqu(il,ju1+1,ik),  &
     &                    qqu(il,ju1+2,ik)) -  &
     &               qqu(il,ju1+1,ik)
!
              pmin = qqu(il,ju1+1,ik) -  &
     &               Min (qqu(il,ju1,ik), qqu(il,ju1+1,ik),  &
     &                    qqu(il,ju1+2,ik))
!
              dcy(il,ju1+1,ik) =  &
     &          Sign (Min (Abs (tmp), pmin, pmax), tmp)
!
            end do
          end do
!
!       ====
        else
!       ====
!
          dcysp_gl(:,:,:) = 0.0d0
          qqusp_gl(:,:,:) = 0.0d0
!
!         ==========================
          call Gmi_Pole_Allgather  &
!         ==========================
     &      (ju1_gl, ju1_gl+3, qqu, qqusp_gl, &
     &  mapi_all, commu_spole, commu_npole, numDomains, numLonDomains, &
     &  i1, i2, ju1, j2, k1, k2, ivert, ilo, ihi, julo, jhi, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl)
!
          do ik = k1, k2
            do il = i1_gl, i2d2
!
              tmp  =  &
     &          ((8.0d0 *  &
     &            (qqusp_gl(il     ,ju1+2,ik) - qqusp_gl(il,ju1  ,ik))) +  &
     &             qqusp_gl(il+i2d2,ju1+1,ik) - qqusp_gl(il,ju1+3,ik)) *  &
     &          r24
!
              pmax = Max (qqusp_gl(il,ju1,ik), qqusp_gl(il,ju1+1,ik),  &
     &                    qqusp_gl(il,ju1+2,ik)) - qqusp_gl(il,ju1+1,ik)
!
              pmin = qqusp_gl(il,ju1+1,ik) -  &
     &               Min (qqusp_gl(il,ju1,ik), qqusp_gl(il,ju1+1,ik),  &
     &                    qqusp_gl(il,ju1+2,ik))
!
              dcysp_gl(il,ju1+1,ik) =  &
     &          Sign (Min (Abs (tmp), pmin, pmax), tmp)
!
            end do
          end do
!
          do ik = k1, k2
            do il = i1_gl + i2d2, i2_gl
!
              tmp  =  &
     &          ((8.0d0 *  &
     &            (qqusp_gl(il     ,ju1+2,ik) - qqusp_gl(il,ju1  ,ik))) +  &
     &             qqusp_gl(il-i2d2,ju1+1,ik) - qqusp_gl(il,ju1+3,ik)) *  &
     &          r24
!
              pmax = Max (qqusp_gl(il,ju1,ik), qqusp_gl(il,ju1+1,ik),  &
     &                    qqusp_gl(il,ju1+2,ik)) - qqusp_gl(il,ju1+1,ik)
!
              pmin = qqusp_gl(il,ju1+1,ik) -  &
     &               Min (qqusp_gl(il,ju1,ik), qqusp_gl(il,ju1+1,ik),  &
     &                    qqusp_gl(il,ju1+2,ik))
!
              dcysp_gl(il,ju1+1,ik) =  &
     &          Sign (Min (Abs (tmp), pmin, pmax), tmp)
!
            end do
          end do
!
          dcy(i1:i2,ju1+1,:) = dcysp_gl(i1:i2,ju1+1,:)
!
!       ======
        end if
!       ======
!
!     ======
      end if
!     ======
!
!
!     ================
      if (j2 == j2_gl) then
!     ================
!
!       ==================
        if (numLonDomains == 1) then
!       ==================
!
          do ik = k1, k2
            do il = i1, i2d2
!
              tmp  =  &
     &          ((8.0d0 * (qqu(il,j2,ik) - qqu(il,j2-2,ik))) +  &
     &           qqu(il,j2-3,ik) - qqu(il+i2d2,j2-1,ik)) *  &
     &          r24
!
              pmax = Max (qqu(il,j2-2,ik), qqu(il,j2-1,ik),  &
     &                    qqu(il,j2,ik)) -  &
     &               qqu(il,j2-1,ik)
!
              pmin = qqu(il,j2-1,ik) -  &
     &               Min (qqu(il,j2-2,ik), qqu(il,j2-1,ik),  &
     &                    qqu(il,j2,ik))
!
              dcy(il,j2-1,ik) =  &
     &          Sign (Min (Abs (tmp), pmin, pmax), tmp)
!
            end do
          end do
!
          do ik = k1, k2
            do il = i1 + i2d2, i2
!
              tmp  =  &
     &          ((8.0d0 * (qqu(il,j2,ik) - qqu(il,j2-2,ik))) +  &
     &           qqu(il,j2-3,ik) - qqu(il-i2d2,j2-1,ik)) *  &
     &          r24
!
              pmax = Max (qqu(il,j2-2,ik), qqu(il,j2-1,ik),  &
     &                    qqu(il,j2,ik)) -  &
     &               qqu(il,j2-1,ik)
!
              pmin = qqu(il,j2-1,ik) -  &
     &               Min (qqu(il,j2-2,ik), qqu(il,j2-1,ik),  &
     &                    qqu(il,j2,ik))
!
              dcy(il,j2-1,ik) =  &
     &          Sign (Min (Abs (tmp), pmin, pmax), tmp)
!
            end do
          end do
!
!       ====
        else
!       ====
!
          dcynp_gl(:,:,:) = 0.0d0
          qqunp_gl(:,:,:) = 0.0d0
!
!         ==========================
          call Gmi_Pole_Allgather  &
!         ==========================
     &      (j2_gl-3, j2_gl, qqu, qqunp_gl, &
     &  mapi_all, commu_spole, commu_npole, numDomains, numLonDomains, &
     &  i1, i2, ju1, j2, k1, k2, ivert, ilo, ihi, julo, jhi, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl)
!
          do ik = k1, k2
            do il = i1_gl, i2d2
!
              tmp  =  &
     &          ((8.0d0 *  &
     &            (qqunp_gl(il,j2,ik) - qqunp_gl(il,j2-2,ik))) +  &
     &             qqunp_gl(il,j2-3,ik) - qqunp_gl(il+i2d2,j2-1,ik)) *  &
     &          r24
!
              pmax = Max (qqunp_gl(il,j2-2,ik), qqunp_gl(il,j2-1,ik),  &
     &                    qqunp_gl(il,j2,ik)) -  &
     &               qqunp_gl(il,j2-1,ik)
!
              pmin = qqunp_gl(il,j2-1,ik) -  &
     &               Min (qqunp_gl(il,j2-2,ik), qqunp_gl(il,j2-1,ik),  &
     &                    qqunp_gl(il,j2,ik))
!
              dcynp_gl(il,j2-1,ik) =  &
     &          Sign (Min (Abs (tmp), pmin, pmax), tmp)
!
            end do
          end do
!
          do ik = k1, k2
            do il = i1_gl + i2d2, i2_gl
!
              tmp  =  &
     &          ((8.0d0 *  &
     &            (qqunp_gl(il,j2,ik) - qqunp_gl(il,j2-2,ik))) +  &
     &             qqunp_gl(il,j2-3,ik) - qqunp_gl(il-i2d2,j2-1,ik)) *  &
     &          r24
!
              pmax = Max (qqunp_gl(il,j2-2,ik), qqunp_gl(il,j2-1,ik),  &
     &                    qqunp_gl(il,j2,ik)) -  &
     &               qqunp_gl(il,j2-1,ik)
!
              pmin = qqunp_gl(il,j2-1,ik) -  &
     &               Min (qqunp_gl(il,j2-2,ik), qqunp_gl(il,j2-1,ik),  &
     &                    qqunp_gl(il,j2,ik))
!
              dcynp_gl(il,j2-1,ik) =  &
     &          Sign (Min (Abs (tmp), pmin, pmax), tmp)
!
            end do
          end do
!
          dcy(i1:i2,j2-1,:) = dcynp_gl(i1:i2,j2-1,:)
!
!       ======
        end if
!       ======
!
!     ======
      end if
!     ======
!
!
      return
!
      end
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Ymist_Pole2_I2d2
!
! DESCRIPTION
!   This routine sets "dcy" at the Poles.
!
! ARGUMENTS
!   dcy : slope of concentration distribution in N-S direction (mixing ratio)
!   qqu : concentration contribution from E-W advection (mixing ratio)
!
!-----------------------------------------------------------------------------
!
      subroutine Do_Ymist_Pole2_I2d2 (dcy, qqu, &
     &   numLonDomains, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   mapi_all, commu_npole, commu_spole, numDomains, ivert)
!
      use GmiGather_mod, only : Gmi_Pole_Allgather
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numLonDomains, numDomains
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, j1p
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ivert
      integer, intent(in) :: commu_npole, commu_spole
      integer, intent(in) :: mapi_all(2, numDomains)
!
      real*8  :: dcy(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: qqu(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: i2d2
      integer :: il, ik
!
      real*8  :: pmax, pmin
      real*8  :: tmp
!
      real*8  :: dcysp_gl(i1_gl:i2_gl, ju1:ju1, k1:k2)
      real*8  :: dcynp_gl(i1_gl:i2_gl,  j2:j2,  k1:k2)
!
      real*8  :: qqusp_gl(i1_gl:i2_gl, ju1:ju1+1, k1:k2)
      real*8  :: qqunp_gl(i1_gl:i2_gl,  j2-1:j2,  k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      i2d2 = i2_gl / 2
!
!
!     ==================
      if (ju1 == ju1_gl) then
!     ==================
!
!       ==================
        if (numLonDomains == 1) then
!       ==================
!
          if (j1p /= ju1_gl+1) then
!
            dcy(i1:i2,ju1,:) = 0.0d0
!
          else
!
!           -----------------------------------------------
!           Determine slope in South Polar cap for scalars.
!           -----------------------------------------------
!
            do ik = k1, k2
              do il = i1, i2d2
!
                tmp  =  &
     &            0.25d0 *  &
     &            (qqu(il,ju1+1,ik) - qqu(il+i2d2,ju1+1,ik))
!
                pmax =  &
     &            Max (qqu(il,ju1+1,ik), qqu(il,ju1,ik),  &
     &                 qqu(il+i2d2,ju1+1,ik)) -  &
     &            qqu(il,ju1,ik)
!
                pmin =  &
     &            qqu(il,ju1,ik) -  &
     &            Min (qqu(il,ju1+1,ik), qqu(il,ju1,ik),  &
     &                 qqu(il+i2d2,ju1+1,ik))
!
                dcy(il,ju1,ik) =  &
     &            Sign (Min (Abs (tmp), pmax, pmin), tmp)
!
              end do
            end do
!
            do il = i1 + i2d2, i2
              dcy(il,ju1,:) = -dcy(il-i2d2,ju1,:)
            end do
!
          end if
!
!       ====
        else
!       ====
!
          if (j1p /= ju1_gl+1) then
!
            dcy(i1:i2,ju1,:) = 0.0d0
!
          else
!
            dcysp_gl(:,:,:) = 0.0d0
            qqusp_gl(:,:,:) = 0.0d0
!
!           ==========================
            call Gmi_Pole_Allgather  &
!           ==========================
     &        (ju1_gl, ju1_gl+1, qqu, qqusp_gl, &
     &  mapi_all, commu_spole, commu_npole, numDomains, numLonDomains, &
     &  i1, i2, ju1, j2, k1, k2, ivert, ilo, ihi, julo, jhi, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl)
!
!           -----------------------------------------------
!           Determine slope in South Polar cap for scalars.
!           -----------------------------------------------
!
            do ik = k1, k2
              do il = i1_gl, i2d2
!
                tmp  =  &
     &            0.25d0 *  &
     &            (qqusp_gl(il,ju1+1,ik) - qqusp_gl(il+i2d2,ju1+1,ik))
!
                pmax =  &
     &            Max (qqusp_gl(il,ju1+1,ik), qqusp_gl(il,ju1,ik),  &
     &                 qqusp_gl(il+i2d2,ju1+1,ik)) -  &
     &            qqusp_gl(il,ju1,ik)
!
                pmin =  &
     &            qqusp_gl(il,ju1,ik) -  &
     &            Min (qqusp_gl(il,ju1+1,ik), qqusp_gl(il,ju1,ik),  &
     &                 qqusp_gl(il+i2d2,ju1+1,ik))
!
                dcysp_gl(il,ju1,ik) =  &
     &            Sign (Min (Abs (tmp), pmax, pmin), tmp)
!
              end do
            end do
!
            do il = i1_gl + i2d2, i2_gl
              dcysp_gl(il,ju1,:) = -dcysp_gl(il-i2d2,ju1,:)
            end do
!
            dcy(i1:i2,ju1,:) = dcysp_gl(i1:i2,ju1,:)
!
          end if
!
!       ======
        end if
!       ======
!
!     ======
      end if
!     ======
!
!
!     ================
      if (j2 == j2_gl) then
!     ================
!
!       ==================
        if (numLonDomains == 1) then
!       ==================
!
          if (j1p /= ju1_gl+1) then
!
            dcy(i1:i2,j2,:) = 0.0d0
!
          else
!
!           -----------------------------------------------
!           Determine slope in North Polar cap for scalars.
!           -----------------------------------------------
!
            do ik = k1, k2
              do il = i1, i2d2
!
                tmp  =  &
     &            0.25d0 *  &
     &            (qqu(il+i2d2,j2-1,ik) - qqu(il,j2-1,ik))
!
                pmax =  &
     &            Max (qqu(il+i2d2,j2-1,ik), qqu(il,j2,ik),  &
     &                 qqu(il,j2-1,ik)) -  &
     &            qqu(il,j2,ik)
!
                pmin =  &
     &            qqu(il,j2,ik) -  &
     &            Min (qqu(il+i2d2,j2-1,ik), qqu(il,j2,ik),  &
     &                 qqu(il,j2-1,ik))
!
                dcy(il,j2,ik) =  &
     &            Sign (Min (Abs (tmp), pmax, pmin), tmp)
!
              end do
            end do
!
            do il = i1 + i2d2, i2
              dcy(il,j2,:) = -dcy(il-i2d2,j2,:)
            end do
!
          end if
!
!       ====
        else
!       ====
!
          if (j1p /= ju1_gl+1) then
!
            dcy(i1:i2,j2,:) = 0.0d0
!
          else
!
            dcynp_gl(:,:,:) = 0.0d0
            qqunp_gl(:,:,:) = 0.0d0
!
!           ==========================
            call Gmi_Pole_Allgather  &
!           ==========================
     &        (j2_gl-1, j2_gl, qqu, qqunp_gl, &
     &  mapi_all, commu_spole, commu_npole, numDomains, numLonDomains, &
     &  i1, i2, ju1, j2, k1, k2, ivert, ilo, ihi, julo, jhi, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl)
!
!           -----------------------------------------------
!           Determine slope in North Polar cap for scalars.
!           -----------------------------------------------
!
            do ik = k1, k2
              do il = i1_gl, i2d2
!
                tmp  =  &
     &            0.25d0 *  &
     &            (qqunp_gl(il+i2d2,j2-1,ik) - qqunp_gl(il,j2-1,ik))
!
                pmax =  &
     &            Max (qqunp_gl(il+i2d2,j2-1,ik), qqunp_gl(il,j2,ik),  &
     &                 qqunp_gl(il,j2-1,ik)) -  &
     &            qqunp_gl(il,j2,ik)
!
                pmin =  &
     &            qqunp_gl(il,j2,ik) -  &
     &            Min (qqunp_gl(il+i2d2,j2-1,ik), qqunp_gl(il,j2,ik),  &
     &                 qqunp_gl(il,j2-1,ik))
!
                dcynp_gl(il,j2,ik) =  &
     &            Sign (Min (Abs (tmp), pmax, pmin), tmp)
!
              end do
            end do
!
            do il = i1_gl + i2d2, i2_gl
              dcynp_gl(il,j2,:) = -dcynp_gl(il-i2d2,j2,:)
            end do
!
            dcy(i1:i2,j2,:) = dcynp_gl(i1:i2,j2,:)
!
          end if
!
!       ======
        end if
!       ======
!
!     ======
      end if
!     ======
!
!
      return
!
      end
!
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Yadv_Pole_Sum
!
! DESCRIPTION
!   This routine sets the cross term due to N-S advection at the Poles.
!
! ARGUMENTS
!   ady : cross term due to N-S advection (mixing ratio)
!
!-----------------------------------------------------------------------------
!
      subroutine Do_Yadv_Pole_Sum (ady, &
     &   numLonDomains, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole)
!
      use GmiReduce_mod, only : Gmi_Sum_Pole_Reduce
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numLonDomains
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, j1p
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: commu_npole, commu_spole
!
      real*8  :: ady(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ik
!
      real*8  :: sumnp(k1:k2)
      real*8  :: sumsp(k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
!     ==================
      if (ju1 == ju1_gl) then
!     ==================
!
        do ik = k1, k2
!
          sumsp(ik) = 0.0d0
!
          if (j1p /= ju1_gl+1) then
!
            do il = i1, i2
!
              sumsp(ik) = sumsp(ik) + ady(il,ju1+1,ik)
!
            end do
!
          else
!
            do il = i1, i2
!
              sumsp(ik) = sumsp(ik) + ady(il,ju1,ik)
!
            end do
!
          end if
!
        end do
!
        if (numLonDomains /= 1) then
!         ===========================
          call Gmi_Sum_Pole_Reduce  &
!         ===========================
     &      (k1, k2, sumsp, ju1, j2, ju1_gl, j2_gl, &
     &       commu_npole, commu_spole)
        end if
!
        do ik = k1, k2
!
          sumsp(ik) = sumsp(ik) / i2_gl
!
          if (j1p /= ju1_gl+1) then
!
            do il = i1, i2
!
              ady(il,ju1+1,ik) = sumsp(ik)
              ady(il,ju1,ik)   = sumsp(ik)
!
            end do
!
          else
!
            do il = i1, i2
!
              ady(il,ju1,ik) = sumsp(ik)
!
            end do
!
          end if
!
        end do
!
!     ======
      end if
!     ======
!
!
!     ================
      if (j2 == j2_gl) then
!     ================
!
        do ik = k1, k2
!
          sumnp(ik) = 0.0d0
!
          if (j1p /= ju1_gl+1) then
!
            do il = i1, i2
!
              sumnp(ik) = sumnp(ik) + ady(il,j2-1,ik)
!
            end do
!
          else
!
            do il = i1, i2
!
              sumnp(ik) = sumnp(ik) + ady(il,j2,ik)
!
            end do
!
          end if
!
        end do
!
        if (numLonDomains /= 1) then
!         ===========================
          call Gmi_Sum_Pole_Reduce  &
!         ===========================
     &      (k1, k2, sumnp, ju1, j2, ju1_gl, j2_gl, &
     &       commu_npole, commu_spole)
        end if
!
        do ik = k1, k2
!
          sumnp(ik) = sumnp(ik) / i2_gl
!
          if (j1p /= ju1_gl+1) then
!
            do il = i1, i2
!
              ady(il,j2-1,ik) = sumnp(ik)
              ady(il,j2,ik)   = sumnp(ik)
!
            end do
!
          else
!
            do il = i1, i2
!
              ady(il,j2,ik) = sumnp(ik)
!
            end do
!
          end if
!
        end do
!
!     ======
      end if
!     ======
!
!
      return
!
      end
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Ytp_Pole_Sum
!
! DESCRIPTION
!   This routine sets "dq1" at the Poles.
!
! ARGUMENTS
!   geofac_pc : special geometrical factor (geofac) for Polar cap
!   dq1       : species density (mb)
!   qqv       : concentration contribution from N-S advection (mixing ratio)
!
!-----------------------------------------------------------------------------
!
      subroutine Do_Ytp_Pole_Sum  &
     &  (geofac_pc, dq1, qqv, fy, &
     &   numLonDomains, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, j2p, ivert, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &       commu_npole, commu_spole)
!
      use GmiReduce_mod, only : Gmi_Sum_Pole_Reduce
      use GmiBroadcast_mod   , only : Gmi_Pole_Broadcast
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numLonDomains, commu_npole, commu_spole
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, j1p, j2p
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ivert
!
      real*8, intent(in)    :: geofac_pc
      real*8, intent(inout) :: dq1(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)    :: qqv(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(inout) :: fy (ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ik
!
      real*8  :: ri2
!
      real*8  :: dq_np(k1:k2)
      real*8  :: dq_sp(k1:k2)
      real*8  :: dqik (k1:k2)
      real*8  :: sumnp(k1:k2)
      real*8  :: sumsp(k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      ri2 = i2_gl
!
      dqik(:) = 0.0d0
!
!... Do south pole
!     ==================
      if (ju1 == ju1_gl) then
!     ==================
!
!... Integrate N-S flux around polar cap lat circle for each level
        do ik = k1, k2
          sumsp(ik) = 0.0d0
!
          do il = i1, i2
            sumsp(ik) = sumsp(ik) + qqv(il,j1p,ik)
          enddo
!
        enddo
!
!... wrap in E-W direction
        if (i1 == i1_gl) then
          dqik(:) = dq1(i1,ju1,:)
        endif
!
!... do if using more than one processor
        if (numLonDomains /= 1) then
!         ===========================
          call Gmi_Sum_Pole_Reduce  &
!         ===========================
     &      (k1, k2, sumsp, ju1, j2, ju1_gl, j2_gl, &
     &       commu_npole, commu_spole)
!
!         ==========================
          call Gmi_Pole_Broadcast  &
!         ==========================
     &      (ivert, dqik, ju1, j2, ju1_gl, j2_gl, &
     &       commu_npole, commu_spole)
        endif
!
!... normalize and set inside polar cap
        do ik = k1, k2
!
          dq_sp(ik) = dqik(ik) - (sumsp(ik) / ri2 * geofac_pc)
!
          do il = i1, i2
            dq1(il,ju1,ik) = dq_sp(ik)
!... save polar flux
            fy(il,ju1,ik) = - (sumsp(ik) / ri2 * geofac_pc)
          enddo
!
          if (j1p /= ju1_gl+1) then
            do il = i1, i2
              dq1(il,ju1+1,ik) = dq_sp(ik)
!... save polar flux
              fy(il,ju1+1,ik) = - (sumsp(ik) / ri2 * geofac_pc)
            enddo
!
          endif
!
        enddo
!
!     ======
      endif
!     ======
!
!
!... Do north pole
!     ================
      if (j2 == j2_gl) then
!     ================
!
        do ik = k1, k2
          sumnp(ik) = 0.0d0
!
          do il = i1, i2
            sumnp(ik) = sumnp(ik) + qqv(il,j2p+1,ik)
          enddo
!
        enddo
!
        if (i1 == i1_gl) then
          dqik(:) = dq1(i1,j2,:)
        endif
!
        if (numLonDomains /= 1) then
!         ===========================
          call Gmi_Sum_Pole_Reduce  &
!         ===========================
     &      (k1, k2, sumnp, ju1, j2, ju1_gl, j2_gl, &
     &       commu_npole, commu_spole)
!
!         ==========================
          call Gmi_Pole_Broadcast  &
!         ==========================
     &      (ivert, dqik, ju1, j2, ju1_gl, j2_gl, &
     &       commu_npole, commu_spole)
        endif
!
        do ik = k1, k2
!
          dq_np(ik) = dqik(ik) + (sumnp(ik) / ri2 * geofac_pc)
!
          do il = i1, i2
            dq1(il,j2,ik) = dq_np(ik)
!... save polar flux
            fy(il,j2+1,ik) = (sumnp(ik) / ri2* geofac_pc)
          enddo
!
          if (j1p /= ju1_gl+1) then
!
            do il = i1, i2
              dq1(il,j2-1,ik) = dq_np(ik)
!... save polar flux
              fy(il,j2,ik) = (sumnp(ik) / ri2* geofac_pc)
            enddo
!
          endif
!
        enddo
!
!     ======
      endif
!     ======
!
!
      return
!
      end
!
!
