!-----------------------------------------------------------------------
!
      module globals
!
! Fortran90 interface to PIC library constants
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: december 14, 2005
!
      implicit none
      public
!
      integer, parameter :: LINEAR = 1, QUADRATIC = 2
      integer, parameter :: STANDARD = 1, LOOKAHEAD = 2, VECTOR = 3
      integer, parameter :: PERIODIC_2D = 1, DIRICHLET_2D = 2
      integer, parameter :: DIRICHLET_PERIODIC_2D = 3, VACUUM_2D = 4
      integer, parameter :: VACUUM_3D = 5, NEUMANN_2D = 6
!
      end module globals
