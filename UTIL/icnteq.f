!
!-----------------------------------------------------------------------
!
!     Count the number of elements equal to a specified integer value.
!
      integer function icnteq (n, array, value)
!
      integer    n, value
      integer    array(*)
!
      k = 0
      do 10 i = 1, n
         if (array(i) .eq. value) k = k + 1
   10 continue
      icnteq = k
!
      return
      end
