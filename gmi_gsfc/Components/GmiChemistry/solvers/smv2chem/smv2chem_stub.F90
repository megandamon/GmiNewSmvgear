
!
!     $Id: smv2chem_stub.F90,v 1.3 2011-08-09 22:12:58 mrdamon Exp $
!

      subroutine Do_Smv2_Init ( )

      use GmiPrintError_mod, only : GmiPrintError

      implicit none


      Write (6,*)
      Write (6,*)  &
     &  'Should not be in this stubbed-out version of Do_Smv2_Init.'
      Write (6,*)

      call GmiPrintError ('  ', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

      return

      end


      subroutine Do_Smv2_Solver ( )

      use GmiPrintError_mod, only : GmiPrintError

      implicit none


      Write (6,*)
      Write (6,*)  &
     &  'Should not be in this stubbed-out version of Do_Smv2_Solver.'
      Write (6,*)

      call GmiPrintError ('  ', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

      return

      end

