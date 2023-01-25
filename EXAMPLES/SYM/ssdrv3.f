      program ssdrv3
!
!     Program to illustrate the idea of reverse communication in
!     inverse mode for a generalized symmetric eigenvalue problem.
!     The following program uses the two LAPACK subroutines sgttrf.f
!     and sgttrs.f to factor and solve a tridiagonal system of equations.
!
!     We implement example three of ex-sym.doc in DOCUMENTS directory
!
!\Example-3
!     ... Suppose we want to solve A*x = lambda*M*x in inverse mode,
!         where A and M are obtained by the finite element of the
!         1-dimensional discrete Laplacian
!                             d^2u / dx^2
!         on the interval [0,1] with zero Dirichlet boundary condition
!         using piecewise linear elements.
!
!     ... OP = inv[M]*A  and  B = M.
!
!     ... Use mode 2 of SSAUPD.
!
!\BeginLib
!
!\Routines called:
!     ssaupd  ARPACK reverse communication interface routine.
!     sseupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     sgttrf  LAPACK tridiagonal factorization routine.
!     sgttrs  LAPACK tridiagonal solve routine.
!     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     sscal   Level 1 BLAS that scales a vector by a scalar.
!     scopy   Level 1 BLAS that copies one vector to another.
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!     av      Matrix vector multiplication routine that computes A*x.
!     mv      Matrix vector multiplication routine that computes M*x.
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
! FILE: sdrv3.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!--------------------------------------------------------------------------
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
      Real
     &                 v(ldv,maxncv), workl(maxncv*(maxncv+8)),
     &                 workd(3*maxn), d(maxncv,2), resid(maxn),
     &                 ad(maxn), adl(maxn), adu(maxn), adu2(maxn),
     &                 ax(maxn), mx(maxn)
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
      Real
     &                 sigma, r1, r2, tol, h
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real
     &                 zero, one, four, six
      parameter        ( zero = 0.0E+0 , one = 1.0E+0 ,
     &                   four = 4.0E+0 , six = 6.0E+0  )
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Real
     &                 snrm2
      external         saxpy, scopy, sscal, snrm2, sgttrf, sgttrs
!
!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
      intrinsic        abs
!
!     %-----------------------%
!     | Executable statements |
!     %-----------------------%
!
!     %----------------------------------------------------%
!     | The number N is the dimension of the matrix. A     |
!     | generalized eigenvalue problem is solved (BMAT =   |
!     | 'G'.) NEV is the number of eigenvalues to be       |
!     | approximated.  The user can modify NEV, NCV, WHICH |
!     | to solve problems of different sizes, and to get   |
!     | different parts of the spectrum.  However, The     |
!     | following conditions must be satisfied:            |
!     |                     N <= MAXN,                     |
!     |                   NEV <= MAXNEV,                   |
!     |               NEV + 1 <= NCV <= MAXNCV             |
!     %----------------------------------------------------%
!
      n = 100
      nev = 4
      ncv = 10
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SDRV3: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SDRV3: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SDRV3: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'G'
      which = 'LM'
!
!     %--------------------------------------------------%
!     | The work array WORKL is used in SSAUPD as        |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in SSAUPD to start the Arnoldi         |
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
!     | iterations allowed.  Mode 2 of SSAUPD is used     |
!     | (IPARAM(7) = 2).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in SSAUPD.                          |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode   = 2
!
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
!
!     %------------------------------------------------%
!     | Call LAPACK routine to factor the mass matrix. |
!     | The mass matrix is the tridiagonal matrix      |
!     | arising from using piecewise linear finite     |
!     | elements on the interval [0, 1].               |
!     %------------------------------------------------%
!
      h = one / real (n+1)
!
      r1 = (four / six) * h
      r2 = (one / six) * h
      do 20 j=1,n
         ad(j) =  r1
         adl(j) = r2
 20   continue
      call scopy (n, adl, 1, adu, 1)
      call sgttrf (n, adl, ad, adu, adu2, ipiv, ierr)
      if (ierr .ne. 0) then
         print *, ' '
         print *, ' Error with _gttrf in _SDRV3. '
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
!        | Repeatedly call the routine SSAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call ssaupd ( ido, bmat, n, which, nev, tol, resid,
     &                 ncv, v, ldv, iparam, ipntr, workd, workl,
     &                 lworkl, info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %--------------------------------------%
!           | Perform  y <--- OP*x = inv[M]*A*x.   |
!           | The user should supply his/her own   |
!           | matrix vector multiplication (A*x)   |
!           | routine and a linear system solver   |
!           | here.  The matrix vector             |
!           | multiplication routine takes         |
!           | workd(ipntr(1)) as the input vector. |
!           | The final result is returned to      |
!           | workd(ipntr(2)). The result of A*x   |
!           | overwrites workd(ipntr(1)).          |
!           %--------------------------------------%
!
            call av (n, workd(ipntr(1)), workd(ipntr(2)))
            call scopy(n, workd(ipntr(2)),  1, workd(ipntr(1)), 1)
            call sgttrs ('Notranspose', n, 1, adl, ad, adu, adu2, ipiv,
     &                  workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print *, ' '
               print *, ' Error with _gttrs in _SDRV3.'
               print *, ' '
               go to 9000
            end if
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call SSAUPD again. |
!           %-----------------------------------------%
!
            go to 10

!
         else if (ido .eq. 2) then
!
!           %-----------------------------------------%
!           |         Perform  y <--- M*x.            |
!           | Need the matrix vector multiplication   |
!           | routine here that takes workd(ipntr(1)) |
!           | as the input and returns the result to  |
!           | workd(ipntr(2)).                        |
!           %-----------------------------------------%
!
            call mv (n, workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call SSAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         end if
!
!
!
!     %-----------------------------------------%
!     | Either we have convergence, or there is |
!     | an error.                               |
!     %-----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in SSAUPD  |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _saupd, info = ',info
         print *, ' Check the documentation of _saupd '
         print *, ' '
!
      else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using SSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
         rvec = .true.
!
         call sseupd ( rvec, 'All', select, d, v, ldv, sigma,
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
!
         if (ierr .ne. 0) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of SSEUPD. |
!           %------------------------------------%
!
            print *, ' '
            print *, ' Error with _seupd, info = ', ierr
            print *, ' Check the documentation of _seupd'
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
!              |  ||  A*x - lambda*M*x ||  |
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
               call mv(n, v(1,j), mx)
               call saxpy(n, -d(j,1), mx, 1, ax, 1)
               d(j,2) = snrm2(n, ax, 1)
               d(j,2) = d(j,2) / abs(d(j,1))
!
 30         continue
!
!           %-----------------------------%
!           | Display computed residuals. |
!           %-----------------------------%
!
            call smout(6, nconv, 2, d, maxncv, -6,
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
         print *, ' _SDRV3 '
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
!     | Done with program ssdrv3. |
!     %---------------------------%
!
 9000 continue
!
      end
!
!------------------------------------------------------------------------
!     Matrix vector subroutine
!     where the matrix is the 1 dimensional mass matrix
!     on the interval [0,1].
!
      subroutine mv (n, v, w)
      integer           n, j
      Real
     &                  v(n),w(n), one, four, six, h
      parameter         (one = 1.0E+0 , four = 4.0E+0 ,
     &                   six = 6.0E+0 )
!
      w(1) = four*v(1) + v(2)
      do 100 j = 2,n-1
         w(j) = v(j-1) + four*v(j) + v(j+1)
  100 continue
      j = n
      w(j) = v(j-1) + four*v(j)
!
!     Scale the vector w by h.
!
      h = one / (real (n+1)*six)
      call sscal(n, h, w, 1)
      return
      end
!
!--------------------------------------------------------------------
!     matrix vector subroutine
!
!     The matrix used is the stiffness matrix obtained from the finite
!     element discretization of the 1-dimensional discrete Laplacian
!     on the interval [0,1] with zero Dirichlet boundary condition using
!     piecewise linear elements.
!
      subroutine av (n, v, w)
      integer           n, j
      Real
     &                  v(n),w(n), two, one, h
      parameter         ( one = 1.0E+0 , two = 2.0E+0  )
!
      w(1) =  two*v(1) - v(2)
      do 100 j = 2,n-1
         w(j) = - v(j-1) + two*v(j) - v(j+1)
  100 continue
      j = n
      w(j) = - v(j-1) + two*v(j)
!
!     Scale the vector w by (1 / h).
!
      h = one / real (n+1)
      call sscal(n, one/h, w, 1)
      return
      end
