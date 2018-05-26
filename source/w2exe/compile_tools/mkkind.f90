PROGRAM MKKIND
!
!  Copyright (c) 1996 - Polyhedron Software Ltd.  All rights reserved.
!
!  This program writes the source code of a Fortran 90 module (normally
!  called F77KINDS) which declares Fortran 90 kinds corresponding to some
!  common but non-standard Fortran 77 types (such as INTEGER*2).
!
!  When SPAG translates Fortran 77 to Fortran 90, it assumes that the
!  module F77KINDS exists, and translates
!
!         INTEGER*2 to INTEGER(I2KIND)       etc.
!
!  I2KIND is a PARAMETER in F77KINDS containing the KIND which
!  corresponds to INTEGER*2.
!
!  Because kind numbers vary between Fortran 90 implementations, the F77KINDS
!  module has to be tailored to the compiler in use.  To do this, simply
!  compile and run MKKIND, using the target compiler, and it will produce a
!  version of F77KINDS for that compiler.  Using this scheme localises the
!  system dependence to the F77KINDS module.
!
!  Note that not all compilers accept all of the non-standard Fortran 77
!  types, and you may have to modify this program to allow for that fact.
!  For example, if LOGICAL*2 is not a supported type, you would have to
!  comment out two lines, the declaration of L2, and the line that writes
!  the parameter definition of L2KIND.
!
   INTEGER*1 i1
   INTEGER*2 i2
   INTEGER*4 i4

   LOGICAL*1 l1
   LOGICAL*2 l2
   LOGICAL*4 l4

   REAL*4 r4
   REAL*8 r8
   DOUBLE PRECISION d

   COMPLEX*8 c8
   COMPLEX*16 c16

   CHARACTER(8),PARAMETER :: FMT = '(A,I2,A)'

   OPEN(11,FILE='F77KINDS.F90')
   WRITE(11,FMT) 'MODULE F77KINDS'
   WRITE(11,FMT) '!'
   WRITE(11,FMT) '! This module was generated automatically by MKKIND.F90'
   WRITE(11,FMT) '!    (supplied with the plusFORT toolkit from Polyhedron &
                   &Software).'
   WRITE(11,FMT) '! It maps common non-standard Fortran 77 types to Fortran &
                   &90 kinds.'
   WRITE(11,FMT) '! Because different Fortran 90 compilers use different kind'
   WRITE(11,FMT) '! numbers, MKKIND.F90 should be compiled and run with the'
   WRITE(11,FMT) '! compiler you intend to use.'
   WRITE(11,FMT) '!'
   WRITE(11,FMT) '   INTEGER,PARAMETER :: &'
   WRITE(11,FMT) '     I1KIND  =' , KIND(i1)  , ',       & ! INTEGER*1'
   WRITE(11,FMT) '     I2KIND  =' , KIND(i2)  , ',       & ! INTEGER*2'
   WRITE(11,FMT) '     I4KIND  =' , KIND(i4)  , ',       & ! INTEGER*4'
   WRITE(11,FMT) '     L1KIND  =' , KIND(l1)  , ',       & ! LOGICAL*1'
   WRITE(11,FMT) '     L2KIND  =' , KIND(l2)  , ',       & ! LOGICAL*2'
   WRITE(11,FMT) '     L4KIND  =' , KIND(l4)  , ',       & ! LOGICAL*4'
   WRITE(11,FMT) '     R4KIND  =' , KIND(r4)  , ',       & ! REAL*4'
   WRITE(11,FMT) '     R8KIND  =' , KIND(r8)  , ',       & ! REAL*8'
   WRITE(11,FMT) '     DPKIND  =' , KIND(d)   , ',       & ! DOUBLE PRECISION'
   WRITE(11,FMT) '     CX8KIND =' , KIND(c8)  , ',       & ! COMPLEX*8'
   WRITE(11,FMT) '     CX16KIND=' , KIND(c16) , '          ! COMPLEX*16'
   WRITE(11,FMT) 'END MODULE F77KINDS'
   CLOSE(11)
END PROGRAM MKKIND
