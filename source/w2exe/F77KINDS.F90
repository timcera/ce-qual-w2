MODULE F77KINDS
!
! This module was generated automatically by MKKIND.F90
!    (supplied with the plusFORT toolkit from Polyhedron Software).
! It maps common non-standard Fortran 77 types to Fortran 90 kinds.
! Because different Fortran 90 compilers use different kind
! numbers, MKKIND.F90 should be compiled and run with the
! compiler you intend to use.
!
   INTEGER,PARAMETER :: &
     I1KIND  = 1,       & ! INTEGER*1
     I2KIND  = 2,       & ! INTEGER*2
     I4KIND  = 4,       & ! INTEGER*4
     L1KIND  = 1,       & ! LOGICAL*1
     L2KIND  = 2,       & ! LOGICAL*2
     L4KIND  = 4,       & ! LOGICAL*4
     R4KIND  = 4,       & ! REAL*4
     R8KIND  = 8,       & ! REAL*8
     DPKIND  = 8,       & ! DOUBLE PRECISION
     CX8KIND = 4,       & ! COMPLEX*8
     CX16KIND= 8          ! COMPLEX*16
END MODULE F77KINDS
