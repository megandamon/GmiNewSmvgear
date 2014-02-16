! $Id: FV_AdvectionEx.F90,v 1.2 2011-08-09 22:12:59 mrdamon Exp $

! *********************************************************************
! *****                      Main Program                          ****
! *****         Finite-Volume Advection Core (Lin/Rood)            ****
! *****         Driven with simple initial conditions              ****
! *********************************************************************

#define I_AM_MAIN

#include "MAPL_Generic.h"

Program FV_AdvectionEx

!=============================================================================
!BOP
!
! !PROGRAM: FV_AdvectionEx
!
! !DESCRIPTION: This is the ``cap'' or main program for the FV Advection
!   (FVadvcore) test example. This initializes a static wind field over 
!   the globe, defines initial values for one tracer quantity, then 
!   executes as many time steps as specified by the duration in CAP.rc
!
!   In this simple application, the cap is
!   very thin, relying on MAPL\_Cap for all of the work.
!   The behavior of the application is controlled through three
!   resource (or configuration) files. The MAPL\_Cap opens 
!   configuration files for itself and its two children. These have
!   the default names: CAP.rc, ROOT.rc, and HISTORY.rc. They must be present
!   in the run directory at run time. The name of the CAP's own
!   resource file is fixed as "CAP.rc", since this is the resource from
!   which the application "boots-up." The other two may be renamed from 
!   file names specified in the CAP.rc.
!
!   {\bf Grids:}\newline  We will use the simplest grid creation strategy
!     supported by MAPL\_Cap: that the cap fully defines the grid and 
!     passes it to the entire hierarchy below it. The grid is a regular lat-lon
!     grid with the poles at grid centers).  Tracer values are defined at 
!     cell centers while winds, Courant numbers and mass-fluxes are defined
!     on the cell walls (on a so-called C-grid) This grid is created 
!     automatically by the CAP, using parameters specified in the CAP.rc. 
!     \newline
!
!   {\bf Time:}\newline  The cap contains the time do loop. It creates the 
!      application's ``universal'' clock and advances it each time through 
!      the loop by an interval referred to as the application's ``heartbeat.''
!      Before this loop, it initializes the entire hierarchy, including
!      the history component, reading any necessary restarts. After the
!      time loop is completed, it finalizes all components, writing any
!      requested checkpoint files.
!
!      The clock may be advanced at either the top 
!      or bottom of the loop depending on the desired behavior of the 
!      children. This is controlled from the configuration. All of this
!      functionality is implemented in the generic MAPL\_Cap and can be 
!      controlled from the CAP.rc. The clock created by MAPL\_Cap uses
!      a Gregorian calendar.
!
!{\bf Resources in the CAP.rc:}
! 
!\makebox[1.1 in][l]{\bf Name } 
!\makebox[3.5 in][l]{\bf Description} 
!\makebox[1.0 in][l]{\bf Units       } 
!\makebox[1.0 in][l]{\bf Default   } 
! 
!\makebox[1.1 in][l]{\bf            CF\_FILE: } 
!\makebox[3.5 in][l]{\bf             Name of ROOT's config file  } 
!\makebox[1.0 in][l]{\bf            none  } 
!\makebox[1.0 in][l]{\bf            "ROOT.rc" } 
!\newline 
!\makebox[1.1 in][l]{\bf            CF\_FILE: } 
!\makebox[3.5 in][l]{\bf             Name of HISTORY's config file  } 
!\makebox[1.0 in][l]{\bf            none  } 
!\makebox[1.0 in][l]{\bf            'HISTORY.rc' } 
!\newline 
!\makebox[1.1 in][l]{\bf            TICK\_FIRST: } 
!\makebox[3.5 in][l]{\bf             Determines when clock is advanced  } 
!\makebox[1.0 in][l]{\bf            1 or 0  } 
!\makebox[1.0 in][l]{\bf            none } 
!\newline 
!\makebox[1.1 in][l]{\bf            BEG\_YY: } 
!\makebox[3.5 in][l]{\bf             Beginning year (integer) } 
!\makebox[1.0 in][l]{\bf            year  } 
!\makebox[1.0 in][l]{\bf            1 } 
!\newline 
!\makebox[1.1 in][l]{\bf            BEG\_MM: } 
!\makebox[3.5 in][l]{\bf             Beginning month (integer 1-12) } 
!\makebox[1.0 in][l]{\bf            month  } 
!\makebox[1.0 in][l]{\bf            1 } 
!\newline 
!\makebox[1.1 in][l]{\bf            BEG\_DD: } 
!\makebox[3.5 in][l]{\bf             Beginning day of month (integer 1-31) } 
!\makebox[1.0 in][l]{\bf            day  } 
!\makebox[1.0 in][l]{\bf            1 } 
!\newline 
!\makebox[1.1 in][l]{\bf            BEG\_H: } 
!\makebox[3.5 in][l]{\bf             Beginning hour of day (integer 0-23) } 
!\makebox[1.0 in][l]{\bf            hour  } 
!\makebox[1.0 in][l]{\bf            0 } 
!\newline 
!\makebox[1.1 in][l]{\bf            BEG\_M: } 
!\makebox[3.5 in][l]{\bf             Beginning minute (integer 0-59) } 
!\makebox[1.0 in][l]{\bf            minute  } 
!\makebox[1.0 in][l]{\bf            0 } 
!\newline 
!\makebox[1.1 in][l]{\bf            BEG\_S: } 
!\makebox[3.5 in][l]{\bf             Beginning second (integer 0-59) } 
!\makebox[1.0 in][l]{\bf            second  } 
!\makebox[1.0 in][l]{\bf            0 } 
!\newline 
!\makebox[1.1 in][l]{\bf            END\_YY: } 
!\makebox[3.5 in][l]{\bf             Ending year (integer) } 
!\makebox[1.0 in][l]{\bf            year  } 
!\makebox[1.0 in][l]{\bf            1 } 
!\newline 
!\makebox[1.1 in][l]{\bf            END\_MM: } 
!\makebox[3.5 in][l]{\bf             Ending month (integer 1-12) } 
!\makebox[1.0 in][l]{\bf            month  } 
!\makebox[1.0 in][l]{\bf            1 } 
!\newline 
!\makebox[1.1 in][l]{\bf            END\_DD: } 
!\makebox[3.5 in][l]{\bf             Ending day of month (integer 1-31) } 
!\makebox[1.0 in][l]{\bf            day  } 
!\makebox[1.0 in][l]{\bf            1 } 
!\newline 
!\makebox[1.1 in][l]{\bf            END\_H: } 
!\makebox[3.5 in][l]{\bf             Ending hour of day (integer 0-23) } 
!\makebox[1.0 in][l]{\bf            hour  } 
!\makebox[1.0 in][l]{\bf            0 } 
!\newline 
!\makebox[1.1 in][l]{\bf            END\_M: } 
!\makebox[3.5 in][l]{\bf             Ending minute (integer 0-59) } 
!\makebox[1.0 in][l]{\bf            minute  } 
!\makebox[1.0 in][l]{\bf            0 } 
!\newline 
!\makebox[1.1 in][l]{\bf            END\_S: } 
!\makebox[3.5 in][l]{\bf             Ending second (integer 0-59) } 
!\makebox[1.0 in][l]{\bf            second  } 
!\makebox[1.0 in][l]{\bf            0 } 
!\newline 
!\makebox[1.1 in][l]{\bf            RUN\_DT: } 
!\makebox[3.5 in][l]{\bf            App Clock Interval (the Heartbeat) } 
!\makebox[1.0 in][l]{\bf            seconds  } 
!\makebox[1.0 in][l]{\bf            none } 
!\newline  
!\makebox[1.1 in][l]{\bf            LATLON: } 
!\makebox[3.5 in][l]{\bf             1 -> regular lat-lon; 0 -> custom grid } 
!\makebox[1.0 in][l]{\bf            0 or 1  } 
!\makebox[1.0 in][l]{\bf            1 } 
!\newline 
!\makebox[1.1 in][l]{\bf            NX: } 
!\makebox[3.5 in][l]{\bf             Processing elements in 1st dimension } 
!\makebox[1.0 in][l]{\bf            none  } 
!\makebox[1.0 in][l]{\bf            1 } 
!\newline 
!\makebox[1.1 in][l]{\bf            NY: } 
!\makebox[3.5 in][l]{\bf             Processing elements in 2nd dimension } 
!\makebox[1.0 in][l]{\bf            none  } 
!\makebox[1.0 in][l]{\bf            1 } 
!\newline 
!\makebox[1.1 in][l]{\bf            IM\_WORLD: } 
!\makebox[3.5 in][l]{\bf             Grid size in 1st dimension } 
!\makebox[1.0 in][l]{\bf            none  } 
!\makebox[1.0 in][l]{\bf            none } 
!\newline 
!\makebox[1.1 in][l]{\bf            JM\_WORLD: } 
!\makebox[3.5 in][l]{\bf             Grid size in 2nd dimension } 
!\makebox[1.0 in][l]{\bf            none  } 
!\makebox[1.0 in][l]{\bf            none } 
!\newline 
!\makebox[1.1 in][l]{\bf            LM: } 
!\makebox[3.5 in][l]{\bf             Grid size in 3rd dimension } 
!\makebox[1.0 in][l]{\bf            none  } 
!\makebox[1.0 in][l]{\bf            1 } 
!\newline 
!\makebox[1.1 in][l]{\bf            GRIDNAME: } 
!\makebox[3.5 in][l]{\bf             Optional grid name } 
!\makebox[1.0 in][l]{\bf            none  } 
!\makebox[1.0 in][l]{\bf            'APPGRID' } 
!\newline 
!\makebox[1.1 in][l]{\bf            IMS: } 
!\makebox[3.5 in][l]{\bf             gridpoints in each PE along 1st dimension } 
!\makebox[1.0 in][l]{\bf            none  } 
!\makebox[1.0 in][l]{\bf            IMS } 
!\newline 
!\makebox[1.1 in][l]{\bf            JMS: } 
!\makebox[3.5 in][l]{\bf             gridpoints in each PE along 2nd dimension } 
!\makebox[1.0 in][l]{\bf            none  } 
!\makebox[1.0 in][l]{\bf            JMS } 
!\newline 
!\makebox[1.1 in][l]{\bf            POLEEDGE: } 
!\makebox[3.5 in][l]{\bf             1->gridedge at pole; 0->gridpoint at pole } 
!\makebox[1.0 in][l]{\bf            0 or 1  } 
!\makebox[1.0 in][l]{\bf            0        } 
!\newline 
!\makebox[1.1 in][l]{\bf            LON0: } 
!\makebox[3.5 in][l]{\bf             Longituce of center of first gridbox } 
!\makebox[1.0 in][l]{\bf            degrees  } 
!\makebox[1.0 in][l]{\bf            -90. } 
!\newline 
!
! !USES:

   use MAPL_Mod
   use FVadvcore_GridCompMod, only:  ROOT_SetServices => SetServices

   implicit none

!EOP

!EOC

   integer           :: STATUS
   character(len=18) :: Iam="FV_AdvectionEx"

   call MAPL_CAP(ROOT_SetServices, rc=STATUS)
   VERIFY_(STATUS)

   call exit(0)

end Program FV_AdvectionEx

!EOC