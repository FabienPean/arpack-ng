!\BeginDoc
!
!\Name: csortc
!
!\Description:
!  Sorts the Complex array in X into the order
!  specified by WHICH and optionally applies the permutation to the
!  Real  array Y.
!
!\Usage:
!  call csortc
!     ( WHICH, APPLY, N, X, Y )
!
!\Arguments
!  WHICH   Character*2.  (Input)
!          'LM' -> sort X into increasing order of magnitude.
!          'SM' -> sort X into decreasing order of magnitude.
!          'LR' -> sort X with real(X) in increasing algebraic order
!          'SR' -> sort X with real(X) in decreasing algebraic order
!          'LI' -> sort X with imag(X) in increasing algebraic order
!          'SI' -> sort X with imag(X) in decreasing algebraic order
!
!  APPLY   Logical.  (Input)
!          APPLY = .TRUE.  -> apply the sorted order to array Y.
!          APPLY = .FALSE. -> do not apply the sorted order to array Y.
!
!  N       Integer.  (INPUT)
!          Size of the arrays.
!
!  X       Complex array of length N.  (INPUT/OUTPUT)
!          This is the array to be sorted.
!
!  Y       Complex array of length N.  (INPUT/OUTPUT)
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Routines called:
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!     Adapted from the sort routine in LANSO.
!
!\SCCS Information: @(#)
! FILE: sortc.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine csortc (which, apply, n, x, y)
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character*2 which
      logical    apply
      integer    n
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Complex
     &           x(0:n-1), y(0:n-1)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i, igap, j
      Complex
     &           temp
      Real
     &           temp1, temp2
!
!     %--------------------%
!     | External functions |
!     %--------------------%
!
      Real
     &           slapy2
!
!     %--------------------%
!     | Intrinsic Functions |
!     %--------------------%
       Intrinsic
     &           real, aimag
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      igap = n / 2
!
      if (which .eq. 'LM') then
!
!        %--------------------------------------------%
!        | Sort X into increasing order of magnitude. |
!        %--------------------------------------------%
!
   10    continue
         if (igap .eq. 0) go to 9000
!
         do 30 i = igap, n-1
            j = i-igap
   20       continue
!
            if (j.lt.0) go to 30
!
            temp1 = slapy2(real(x(j)),aimag(x(j)))
            temp2 = slapy2(real(x(j+igap)),aimag(x(j+igap)))
!
            if (temp1.gt.temp2) then
                temp = x(j)
                x(j) = x(j+igap)
                x(j+igap) = temp
!
                if (apply) then
                    temp = y(j)
                    y(j) = y(j+igap)
                    y(j+igap) = temp
                end if
            else
                go to 30
            end if
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
!
      else if (which .eq. 'SM') then
!
!        %--------------------------------------------%
!        | Sort X into decreasing order of magnitude. |
!        %--------------------------------------------%
!
   40    continue
         if (igap .eq. 0) go to 9000
!
         do 60 i = igap, n-1
            j = i-igap
   50       continue
!
            if (j .lt. 0) go to 60
!
            temp1 = slapy2(real(x(j)),aimag(x(j)))
            temp2 = slapy2(real(x(j+igap)),aimag(x(j+igap)))
!
            if (temp1.lt.temp2) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
!
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
!
      else if (which .eq. 'LR') then
!
!        %------------------------------------------------%
!        | Sort XREAL into increasing order of algebraic. |
!        %------------------------------------------------%
!
   70    continue
         if (igap .eq. 0) go to 9000
!
         do 90 i = igap, n-1
            j = i-igap
   80       continue
!
            if (j.lt.0) go to 90
!
            if (real(x(j)).gt.real(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
!
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
!
      else if (which .eq. 'SR') then
!
!        %------------------------------------------------%
!        | Sort XREAL into decreasing order of algebraic. |
!        %------------------------------------------------%
!
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
!
            if (j.lt.0) go to 120
!
            if (real(x(j)).lt.real(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
!
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
!
      else if (which .eq. 'LI') then
!
!        %--------------------------------------------%
!        | Sort XIMAG into increasing algebraic order |
!        %--------------------------------------------%
!
  130    continue
         if (igap .eq. 0) go to 9000
         do 150 i = igap, n-1
            j = i-igap
  140       continue
!
            if (j.lt.0) go to 150
!
            if (aimag(x(j)).gt.aimag(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
!
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 150
            endif
            j = j-igap
            go to 140
  150    continue
         igap = igap / 2
         go to 130
!
      else if (which .eq. 'SI') then
!
!        %---------------------------------------------%
!        | Sort XIMAG into decreasing algebraic order  |
!        %---------------------------------------------%
!
  160    continue
         if (igap .eq. 0) go to 9000
         do 180 i = igap, n-1
            j = i-igap
  170       continue
!
            if (j.lt.0) go to 180
!
            if (aimag(x(j)).lt.aimag(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
!
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 180
            endif
            j = j-igap
            go to 170
  180    continue
         igap = igap / 2
         go to 160
      end if
!
 9000 continue
      return
!
!     %---------------%
!     | End of csortc |
!     %---------------%
!
      end
