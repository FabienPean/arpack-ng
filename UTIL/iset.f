!
!-----------------------------------------------------------------------
!
!     Only work with increment equal to 1 right now.
!
      subroutine iset (n, value, array, inc)
!
      integer    n, value, inc
      integer    array(*)
!
      do 10 i = 1, n
         array(i) = value
   10 continue
!
      return
      end
