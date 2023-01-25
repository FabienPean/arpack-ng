      program dsdrv2
!
!     Program to illustrate the idea of reverse communication
!     in shift and invert mode for a standard symmetric eigenvalue
!     problem.  The following program uses the two LAPACK subroutines
!     dgttrf.f and dgttrs.f to factor and solve a tridiagonal system of
!     equations.
!
!     We implement example two of ex-sym.doc in DOCUMENTS directory
!
!\Example-2
!     ... Suppose we want to solve A*x = lambda*x in shift-invert mode,
!         where A is derived from the central difference discretization
!         of the 1-dimensional Laplacian on [0,1]  with zero Dirichlet
!         boundary condition.
!     ... OP = (inv[A - sigma*I]) and  B = I.
!     ... Use mode 3 of DSAUPD.
!
!\BeginLib
!
!\Routines called:
!     dsaupd  ARPACK reverse communication interface routine.
!     dseupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     dgttrf  LAPACK tridiagonal factorization routine.
!     dgttrs  LAPACK tridiagonal solve routine.
!     daxpy   daxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     av      Matrix vector multiplication routine that computes A*x.
!
!\Author
!     Richard Lehoucq
!     Danny Sorensen
!     Chao Yang
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: sdrv2.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!----------------------------------------------------------------------
!
!     %-----------------------------%
!     | Define leading dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!
      integer          maxn, maxnev, maxncv, ldv
      parameter        (maxn=256, maxnev=10, maxncv=25,
     &                 ldv=maxn )
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      Double precision
     &                 v(ldv,maxncv), workl(maxncv*(maxncv+8)),
     &                 workd(3*maxn), d(maxncv,2), resid(maxn),
     &                 ad(maxn), adl(maxn), adu(maxn), adu2(maxn),
     &                 ax(maxn)
      logical          select(maxncv)
      integer          iparam(11), ipntr(11), ipiv(maxn)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character        bmat*1, which*2
      integer          ido, n, nev, ncv, lworkl, info, j, ierr,
     &                 nconv, maxitr, ishfts, mode
      logical          rvec
      Double precision
     &                 sigma, tol, h2
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Double precision
     &                 zero, one, two
      parameter        (zero = 0.0D+0, one = 1.0D+0,
     &                  two = 2.0D+0)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Double precision
     &                 dnrm2
      external         daxpy, dnrm2, dgttrf, dgttrs
!
!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
      intrinsic        abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %----------------------------------------------------%
!     | The number N is the dimension of the matrix.  A    |
!     | standard eigenvalue problem is solved (BMAT = 'I'. |
!     | NEV is the number of eigenvalues (closest to       |
!     | SIGMA) to be approximated.  Since the shift-invert |
!     | mode is used, WHICH is set to 'LM'.  The user can  |
!     | modify NEV, NCV, SIGMA to solve problems of        |
!     | different sizes, and to get different parts of the |
!     | spectrum.  However, The following conditions must  |
!     | be satisfied:                                      |
!     |                   N <= MAXN,                       |
!     |                 NEV <= MAXNEV,                     |
!     |             NEV + 1 <= NCV <= MAXNCV               |
!     %----------------------------------------------------%
!
      n = 100
      nev = 4
      ncv = 10
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SDRV2: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SDRV2: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SDRV2: NCV is greater than MAXNCV '
         go to 9000
      end if
!
      bmat = 'I'
      which = 'LM'
      sigma = zero
!
!     %--------------------------------------------------%
!     | The work array WORKL is used in DSAUPD as        |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in DSAUPD to start the Arnoldi         |
!     | iteration.                                       |
!     %--------------------------------------------------%
!
      lworkl = ncv*(ncv+8)
      tol = zero
      ido = 0
      info = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 3 of DSAUPD is used     |
!     | (IPARAM(7) = 3).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode   = 3
!
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
!
!     %-----------------------------------------------------%
!     | Call LAPACK routine to factor (A-SIGMA*I), where A  |
!     | is the 1-d Laplacian.                               |
!     %-----------------------------------------------------%
!
      h2 = one / dble((n+1)*(n+1))
      do 20 j=1,n
         ad(j) = two / h2 - sigma
         adl(j) = -one / h2
 20   continue
      call dcopy (n, adl, 1, adu, 1)
      call dgttrf (n, adl, ad, adu, adu2, ipiv, ierr)
      if (ierr .ne. 0) then
         print *, ' '
         print *, ' Error with _gttrf in SDRV2.'
         print *, ' '
         go to 9000
      end if
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call dsaupd ( ido, bmat, n, which, nev, tol, resid,
     &                 ncv, v, ldv, iparam, ipntr, workd, workl,
     &                 lworkl, info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %----------------------------------------%
!           | Perform y <-- OP*x = inv[A-sigma*I]*x. |
!           | The user only need the linear system   |
!           | solver here that takes workd(ipntr(1)) |
!           | as input, and returns the result to    |
!           | workd(ipntr(2)).                       |
!           %----------------------------------------%
!
            call dcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
!
            call dgttrs ('Notranspose', n, 1, adl, ad, adu, adu2, ipiv,
     &                   workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print *, ' '
               print *, ' Error with _gttrs in _SDRV2. '
               print *, ' '
               go to 9000
            end if
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         end if
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %----------------------------%
!        | Error message.  Check the  |
!        | documentation in DSAUPD    |
!        %----------------------------%
!
         print *, ' '
         print *, ' Error with _saupd, info = ',info
         print *, ' Check documentation of _saupd '
         print *, ' '
!
      else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
         rvec = .true.
!
         call dseupd ( rvec, 'All', select, d, v, ldv, sigma,
     &        bmat, n, which, nev, tol, resid, ncv, v, ldv,
     &        iparam, ipntr, workd, workl, lworkl, ierr )
!
!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%

         if ( ierr .ne. 0 ) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DSEUPD. |
!           %------------------------------------%
!
            print *, ' '
            print *, ' Error with _seupd, info = ', ierr
            print *, ' Check the documentation of _seupd '
            print *, ' '
!
         else
!
            nconv =  iparam(5)
            do 30 j=1, nconv
!
!              %---------------------------%
!              | Compute the residual norm |
!              |                           |
!              |   ||  A*x - lambda*x ||   |
!              |                           |
!              | for the NCONV accurately  |
!              | computed eigenvalues and  |
!              | eigenvectors.  (iparam(5) |
!              | indicates how many are    |
!              | accurate to the requested |
!              | tolerance)                |
!              %---------------------------%
!
               call av(n, v(1,j), ax)
               call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
               d(j,2) = dnrm2(n, ax, 1)
               d(j,2) = d(j,2) / abs(d(j,1))
!
 30         continue
!
!           %-------------------------------%
!           | Display computed residuals    |
!           %-------------------------------%
!
            call dmout(6, nconv, 2, d, maxncv, -6,
     &           'Ritz values and relative residuals')
         end if
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' '
            print *, ' No shifts could be applied during implicit',
     &               ' Arnoldi update, try increasing NCV.'
            print *, ' '
         end if
!
         print *, ' '
         print *, ' _SDRV2 '
         print *, ' ====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated',
     &            ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ',
     &              nconv
         print *, ' The number of Implicit Arnoldi update',
     &            ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
!
      end if
!
!     %---------------------------%
!     | Done with program dsdrv2. |
!     %---------------------------%
!
 9000 continue
!
      end
!
!------------------------------------------------------------------------
!     Matrix vector subroutine
!     where the matrix is the 1 dimensional discrete Laplacian on
!     the interval [0,1] with zero Dirichlet boundary condition.
!
      subroutine av (n, v, w)
      integer           n, j
      Double precision
     &                  v(n), w(n), one, two, h2
      parameter         (one = 1.0D+0,  two = 2.0D+0)

!
      w(1) =  two*v(1) - v(2)
      do 100 j = 2,n-1
         w(j) = - v(j-1) + two*v(j) - v(j+1)
  100 continue
      j = n
      w(j) = - v(j-1) + two*v(j)
!
!     Scale the vector w by (1 / h^2).
!
      h2 = one / dble((n+1)*(n+1))
      call dscal(n, one/h2, w, 1)
      return
      end
