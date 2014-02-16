!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiProdLossDiagnostics_mod
!
      module GmiProdLossDiagnostics_mod
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration,       &
     &       Get_prod, Get_loss, Set_prod, Set_loss, Get_concentration,       &
     &       Get_net_cum_tend, Set_net_cum_tend

      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Get_pr_tend,        &
     &       Get_pr_ascii5, Get_num_tend_outrecs, Get_tend_outrec_map
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: doDiagnosticsBefore, doDiagnosticsAfter
      public  :: doDiagnosticsBefore_BC, doDiagnosticsAfter_BC

#     include "gmi_diag_constants_llnl.h"
#     include "GmiParameters.h"
!
                       ! mass difference before and after an operator
      real*8, pointer, save :: delta_mass (:) => null()
                       ! mass before operator
      real*8, pointer, save :: mass_before(:) => null()
                       ! mass after  operator
      real*8, pointer, save :: mass_after (:) => null()
      logical, save :: pr_ascii5, pr_tend
      integer, save :: num_tend_outrecs
      integer, save :: tend_outrec_map (MAX_NUM_CONST_GIO)
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doDiagnosticsBefore_BC
!
! !INTERFACE:
!
      subroutine doDiagnosticsBefore_BC(concentration, prod, loss, mass,  &
     &                 operator, i1, i2, ju1, j2, k1, k2, numSpecies)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer             , intent(in) :: operator ! operator that P/L diagnostics
                                                   ! are being saved for 
      integer             , intent(in) :: i1, i2, ju1, j2, k1, k2
      integer             , intent(in) :: numSpecies
      real*8              , intent(in) :: mass(i1:i2,ju1:j2,k1:k2)
      type (t_GmiArrayBundle), intent(in) :: concentration(numSpecies)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: prod(2, NUM_OPERATORS, 1:numSpecies)
      real*8 , intent(inOut) :: loss(2, NUM_OPERATORS, 1:numSpecies)
!
! !DESCRIPTION:
! Does any necessary diagnostics prior to Forced BC operator being called.
!
! !LOCAL VARIABLES:
      logical, save :: first = .true.
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_ascii5) then
         call saveProdLoss (operator, BEFORE, prod, loss, mass,             &
     &            concentration, numSpecies, i1, i2, ju1, j2, k1, k2)
      end if

      return

      end subroutine doDiagnosticsBefore_BC
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doDiagnosticsAfter_BC
!
! !INTERFACE:
!
      subroutine doDiagnosticsAfter_BC(concentration, prod, loss, mass, &
     &                operator, i1, i2, ju1, j2, k1, k2, numSpecies)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer          , intent(in) :: operator ! operator that P/L diagnostics
                                                !  are being saved for
      integer          , intent(in) :: i1, i2, ju1, j2, k1, k2
      integer          , intent(in) :: numSpecies
      real*8           , intent(in) :: mass(i1:i2,ju1:j2,k1:k2)
      type (t_GmiArrayBundle), intent(in) :: concentration(numSpecies)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: prod(2, NUM_OPERATORS, 1:numSpecies)
      real*8 , intent(inOut) :: loss(2, NUM_OPERATORS, 1:numSpecies)
!
! !DESCRIPTION:
! Does any necessary diagnostics after the Forced BC operator has been called.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_ascii5) then
         call saveProdLoss (operator, AFTER, prod, loss, mass,              &
     &            concentration, numSpecies, i1, i2, ju1, j2, k1, k2)
      end if

      return

      end subroutine doDiagnosticsAfter_BC
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doDiagnosticsBefore
!
! !INTERFACE:
!
      subroutine doDiagnosticsBefore(Diagnostics, SpeciesConcentration, mass,  &
     &   operator, i1, i2, ju1, j2, k1, k2, numSpecies)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer             , intent(in) :: operator ! operator that P/L diagnostics
                                                   ! are being saved for 
      integer             , intent(in) :: i1, i2, ju1, j2, k1, k2
      integer             , intent(in) :: numSpecies
      real*8              , intent(in) :: mass(i1:i2,ju1:j2,k1:k2)
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !DESCRIPTION:
! Does any necessary diagnostics prior to an operator being  called.
!
! !LOCAL VARIABLES:
      logical, save :: first = .true.
      real*8 , allocatable :: prod(:,:,:)
      real*8 , allocatable :: loss(:,:,:)
      type (t_GmiArrayBundle), pointer :: net_cum_tend(:,:)
      type (t_GmiArrayBundle), pointer :: concentration(:)
!EOP
!------------------------------------------------------------------------------
!BOC
      if (first) then
         first = .false.
         Allocate (delta_mass (1:numSpecies))
         Allocate (mass_after (1:numSpecies))
         Allocate (mass_before(1:numSpecies))
         delta_mass  = 0.0d0
         mass_before = 0.0d0
         mass_after  = 0.0d0

         call Get_pr_tend(Diagnostics, pr_tend)
         call Get_pr_ascii5(Diagnostics, pr_ascii5)
         call Get_tend_outrec_map(Diagnostics, tend_outrec_map)
         call Get_num_tend_outrecs(Diagnostics, num_tend_outrecs)
      end if

      if (pr_ascii5 .or. pr_tend) then
         call Get_concentration(SpeciesConcentration, concentration)

         if (pr_ascii5) then

            Allocate (loss(2, NUM_OPERATORS, 1:numSpecies))
            Allocate (prod(2, NUM_OPERATORS, 1:numSpecies))

            call Get_prod(SpeciesConcentration, prod)
            call Get_loss(SpeciesConcentration, loss)

            call saveProdLoss (operator, BEFORE, prod, loss, mass,             &
     &               concentration, numSpecies, i1, i2, ju1, j2, k1, k2)

            call Set_prod(SpeciesConcentration, prod)
            call Set_loss(SpeciesConcentration, loss)

            deallocate(prod)
            deallocate(loss)
         end if

         if (pr_tend) then
            call Get_net_cum_tend(SpeciesConcentration, net_cum_tend)

            call saveTendencies (operator, BEFORE, concentration, mass,       &
     &               net_cum_tend, i1, i2, ju1, j2, k1, k2, numSpecies)

            call Set_net_cum_tend(SpeciesConcentration, net_cum_tend)
         end if
      end if


      return

      end subroutine doDiagnosticsBefore
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doDiagnosticsAfter
!
! !INTERFACE:
!
      subroutine doDiagnosticsAfter(SpeciesConcentration, mass,        &
     &   operator, i1, i2, ju1, j2, k1, k2, numSpecies)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer          , intent(in) :: operator ! operator that P/L diagnostics
                                                !  are being saved for
      integer          , intent(in) :: i1, i2, ju1, j2, k1, k2
      integer          , intent(in) :: numSpecies
      real*8           , intent(in) :: mass(i1:i2,ju1:j2,k1:k2)
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !DESCRIPTION:
! Does any necessary diagnostics after an operator has been called.
!
! !LOCAL VARIABLES:
      real*8 , allocatable :: prod(:,:,:)
      real*8 , allocatable :: loss(:,:,:)
      type (t_GmiArrayBundle), pointer :: net_cum_tend(:,:)
      type (t_GmiArrayBundle), pointer :: concentration(:)
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_ascii5 .or. pr_tend) then
         call Get_concentration(SpeciesConcentration, concentration)

         if (pr_tend) then
            call Get_net_cum_tend(SpeciesConcentration, net_cum_tend)

            call saveTendencies (operator, AFTER, concentration, mass,        &
     &              net_cum_tend, i1, i2, ju1, j2, k1, k2, numSpecies)

            call Set_net_cum_tend(SpeciesConcentration, net_cum_tend)
         end if

         if (pr_ascii5) then

            Allocate (loss(2, NUM_OPERATORS, 1:numSpecies))
            Allocate (prod(2, NUM_OPERATORS, 1:numSpecies))

            call Get_prod(SpeciesConcentration, prod)
            call Get_loss(SpeciesConcentration, loss)

            call saveProdLoss (operator, AFTER, prod, loss, mass,              &
     &               concentration, numSpecies, i1, i2, ju1, j2, k1, k2)

            call Set_prod(SpeciesConcentration, prod)
            call Set_loss(SpeciesConcentration, loss)

            deallocate(prod)
            deallocate(loss)
         end if
      end if

      return

      end subroutine doDiagnosticsAfter
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: saveProdLoss
!
! !INTERFACE:
!
      subroutine saveProdLoss (operator, time_relation, prod, loss, mass,      &
     &               concentration, numSpecies, i1, i2, ju1, j2, k1, k2)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer                , intent(in) :: operator
      integer                , intent(in) :: time_relation
      integer                , intent(in) :: numSpecies, i1, i2, ju1, j2, k1, k2
      real*8                 , intent(in) :: mass(i1:i2,ju1:j2,k1:k2)
      type (t_GmiArrayBundle), intent(in) :: concentration(numSpecies)
!
! !INPUT/OUTPUT PARAMETERS
      real*8, intent(inOut) :: loss(2, NUM_OPERATORS, 1:numSpecies)
      real*8, intent(inOut) :: prod(2, NUM_OPERATORS, 1:numSpecies)
!
! !DESCRIPTION:
! Saves the production and loss diagnostics either before or after an operator.
! \begin{verbatim}
!   Usage example =>
!
!       call saveProdLoss (EMISS_OP, BEFORE)
!
!       =================
!       call Update_Emiss (arg1, arg2, ...)
!       =================
!
!       call saveProdLoss (EMISS_OP, AFTER)
! \end{verbatim}
!
! !LOCAL VARIABLES:
      logical, save :: first = .true.
      integer :: ic
!
!EOP
!------------------------------------------------------------------------------
!BOC
!     ==========
      if (first) then
!     ==========

        first = .false.

        prod(CUMULATIVE,:,:) = 0.0d0
        loss(CUMULATIVE,:,:) = 0.0d0

      end if

!     ============================
      if (time_relation == BEFORE) then
!     ============================

        do ic = 1, numSpecies
          mass_before(ic) = Sum (mass(:,:,:) * concentration(ic)%pArray3D(:,:,:))
        end do

      end if

!     ===========================
      if (time_relation == AFTER) then
!     ===========================

        do ic = 1, numSpecies

          prod(ONE_STEP,operator,ic) = 0.0d0
          loss(ONE_STEP,operator,ic) = 0.0d0

          mass_after(ic) = Sum (mass (:,:,:) * concentration(ic)%pArray3D(:,:,:))

          delta_mass(ic) =  mass_after(ic) - mass_before(ic)

          if (delta_mass(ic) > 0.0d0) then

            prod(ONE_STEP,operator,ic) =  &
     &        prod(ONE_STEP,operator,ic) + delta_mass(ic)

          else

            loss(ONE_STEP,operator,ic) =  &
     &        loss(ONE_STEP,operator,ic) - delta_mass(ic)

          end if

          prod(CUMULATIVE,operator,ic) =  &
     &      prod(CUMULATIVE,operator,ic) + prod(ONE_STEP,  operator,ic)

          loss(CUMULATIVE,operator,ic) =  &
     &      loss(CUMULATIVE,operator,ic) + loss(ONE_STEP,  operator,ic)

        end do

      end if

      return

      end subroutine saveProdLoss
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: saveTendencies
!
      subroutine saveTendencies (operator, time_relation, concentration, mass, &
     &               net_cum_tend, i1, i2, ju1, j2, k1, k2, numSpecies)
!
      implicit none
!
! !INPUT PARAMETERS:
      ! operator that net cumulative P/L diagnostics are being saved for
      integer                , intent(in) :: operator
      ! set to either BEFORE or AFTER the call to the operator
      integer                , intent(in) :: time_relation
      integer                , intent(in) :: i1, i2, ju1, j2, k1, k2
      integer                , intent(in) :: numSpecies
      real*8                 , intent(in) :: mass(i1:i2,ju1:j2,k1:k2)
      type (t_GmiArrayBundle), intent(in) :: concentration(numSpecies)
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_GmiArrayBundle), intent(inOut) :: net_cum_tend(num_tend_outrecs,NUM_OPERATORS)
!
! !DESCRIPTION:
!   Saves the net cumulative tendency diagnostics for netCDF
!   output, either before or after an operator.
!
! !LOCAL VARIABLES:
      logical, save :: first = .true.
      integer :: ic, icx
      real*8, allocatable       :: mass_afterx (:,:,:,:)
      real*8, allocatable, save :: mass_beforex(:,:,:,:)
!EOP
!------------------------------------------------------------------------------
!BOC

!     ==========
      if (first) then
!     ==========

        first = .false.

        do ic = 1, NUM_OPERATORS
           do icx = 1, num_tend_outrecs
              net_cum_tend(icx,ic)%pArray3D(:,:,:) = 0.0d0
           end do
        end do

        Allocate (mass_beforex(i1:i2, ju1:j2, k1:k2, num_tend_outrecs))
        mass_beforex = 0.0d0

      end if

!     ============================
      if (time_relation == BEFORE) then
!     ============================
         do icx = 1, num_tend_outrecs
            ic = tend_outrec_map(icx)
            mass_beforex(:,:,:,icx) =  &
     &           mass (:,:,:) * concentration(ic)%pArray3D(:,:,:)
         end do
      end if

!     ===========================
      if (time_relation == AFTER) then
!     ===========================
         Allocate (mass_afterx(i1:i2, ju1:j2, k1:k2, num_tend_outrecs))
         do icx = 1, num_tend_outrecs
            ic = tend_outrec_map(icx)
             mass_afterx(:,:,:,icx) =  &
     &           mass (:,:,:) * concentration(ic)%pArray3D(:,:,:)
             net_cum_tend(icx,operator)%pArray3D(:,:,:) =  &
     &           net_cum_tend(icx,operator)%pArray3D(:,:,:) +  &
     &           (mass_afterx(:,:,:,icx) - mass_beforex(:,:,:,icx))
         end do
         deallocate(mass_afterx)
      end if

      return

      end subroutine saveTendencies
!EOC
!------------------------------------------------------------------------------
      end module GmiProdLossDiagnostics_mod
