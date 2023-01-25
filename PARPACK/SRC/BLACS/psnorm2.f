!\BeginDoc
!
!\Name: psnorm2
!
! Message Passing Layer: BLACS
!
!\Description:
!
!\Usage:
!  call psnorm2 ( COMM, N, X, INC )
!
!\Arguments
!  COMM    BLACS Communicator for the processor grid.  (INPUT)
!
!\SCCS Information:
! FILE: norm2.F   SID: 1.2   DATE OF SID: 2/22/96
!
!-----------------------------------------------------------------------
!
      Real function psnorm2 ( comm, n, x, inc )
!
!     %------------------------------%
!     | BLACS Variables and Routines |
!     %------------------------------%
!
      integer    comm
      external   sgsum2d, sgamx2d
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
      Real
     &             x(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      Real
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
      Real
     &             snrm2
      External     snrm2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      psnorm2 = snrm2( n, x, inc)
!
      max = psnorm2
      call sgamx2d( comm, 'All', ' ', 1, 1, max, 1, ra, ca,
     &              -1, -1, -1 )
      if ( max .eq. zero ) then
         psnorm2 = zero
      else
         psnorm2 = (psnorm2/max)**2.0
         call sgsum2d( comm, 'All', ' ', 1, 1, psnorm2, 1, -1, -1 )
         psnorm2 = max * sqrt(abs(psnorm2))
      endif
!
!     %----------------%
!     | End of psnorm2 |
!     %----------------%
!
      return
      end
