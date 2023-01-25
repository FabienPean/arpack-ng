!\BeginDoc
!
!\Name: pdznorm2
!
! Message Passing Layer: BLACS
!
!\Description:
!
!\Usage:
!  call pdznorm2 ( COMM, N, X, INC )
!
!\Arguments
!  COMM    BLACS Communicator for the processor grid.  (INPUT)
!
!\SCCS Information:
! FILE: norm2.F   SID: 1.2   DATE OF SID: 3/6/96
!
!-----------------------------------------------------------------------
!
      Double precision function pdznorm2 ( comm, n, x, inc )
!
!     %------------------------------%
!     | BLACS Variables and Routines |
!     %------------------------------%
!
      integer    comm
      external   dgsum2d, dgamx2d
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer      n, inc
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Complex*16
     &             x(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      Double precision
     &             max, buf, zero
      parameter    ( zero = 0.0 )
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    abs, sqrt
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Double precision
     &             dznrm2
      External     dznrm2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      pdznorm2 = dznrm2( n, x, inc)
!
      max = pdznorm2
      call dgamx2d( comm, 'All', ' ', 1, 1, max, 1, ra, ca,
     &              -1, -1, -1 )
      if ( max .eq. zero ) then
         pdznorm2 = zero
      else
         pdznorm2 = (pdznorm2/max)**2.0
         call dgsum2d( comm, 'All', ' ', 1, 1, pdznorm2, 1, -1, -1 )
         pdznorm2 = max * sqrt(abs(pdznorm2))
      endif
!
!     %-----------------%
!     | End of pdznorm2 |
!     %-----------------%
!
      return
      end
