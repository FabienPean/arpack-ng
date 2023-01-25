      program sndrv2
!
!     Simple program to illustrate the idea of reverse communication
!     in shift-invert mode for a standard nonsymmetric eigenvalue problem.
!
!     We implement example two of ex-nonsym.doc in DOCUMENTS directory
!
!\Example-2
!     ... Suppose we want to solve A*x = lambda*x in shift-invert mode,
!         where A is derived from the centered difference discretization
!         of the 1-dimensional convection-diffusion operator
!                          (d^2u / dx^2) + rho*(du/dx)
!         on the interval [0,1] with zero Dirichlet boundary condition.
!
!     ... The shift sigma is a real number.
!
!     ... OP = inv[A-sigma*I] and  B = I.
!
!     ... Use mode 3 of SNAUPD.
!
!\BeginLib
!
!\Routines called:
!     snaupd  ARPACK reverse communication interface routine.
!     sneupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     sgttrf  LAPACK tridiagonal factorization routine.
!     sgttrs  LAPACK tridiagonal solve routine.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     scopy   Level 1 BLAS that copies one vector to another.
!     sdot    Level 1 BLAS that computes the dot product of two vectors.
!     snrm2   Level 1 BLAS that computes the norm of a vector.
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
! FILE: ndrv2.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
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
      integer           maxn, maxnev, maxncv, ldv
      parameter         (maxn=256, maxnev=10, maxncv=25,
     &                   ldv=maxn )
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      integer           iparam(11), ipntr(14), ipiv(maxn)
      logical           select(maxncv)
      Real
     &                  ax(maxn), d(maxncv,3), resid(maxn),
     &                  v(ldv, maxncv), workd(3*maxn),
     &                  workev(3*maxncv),
     &                  workl(3*maxncv*maxncv+6*maxncv),
     &                  dd(maxn), dl(maxn), du(maxn),
     &                  du2(maxn)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, j,
     &                  ierr, nconv, maxitr, ishfts, mode
      Real
     &                  tol, h, s,
     &                  sigmar, sigmai, s1, s2, s3
      logical           first, rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real
     &                   one, zero, two, rho
      common            /convct/ rho
      parameter         (one = 1.0E+0, zero = 0.0E+0,
     &                   two = 2.0E+0)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Real
     &                  sdot, snrm2, slapy2
      external          sgttrf, sgttrs, sdot, snrm2, slapy2
!
!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
      intrinsic         abs
!
!     %-----------------------%
!     | Executable statements |
!     %-----------------------%
!
!     %--------------------------------------------------%
!     | The number N is the dimension of the matrix.  A  |
!     | standard eigenvalue problem is solved (BMAT =    |
!     | 'I').  NEV is the number of eigenvalues (closest |
!     | to the shift SIGMAR) to be approximated.  Since  |
!     | the shift-invert mode is used, WHICH is set to   |
!     | 'LM'.  The user can modify NEV, NCV, SIGMAR to   |
!     | solve problems of different sizes, and to get    |
!     | different parts of the spectrum.  However, The   |
!     | following conditions must be satisfied:          |
!     |                 N <= MAXN,                       |
!     |               NEV <= MAXNEV,                     |
!     |           NEV + 2 <= NCV <= MAXNCV               |
!     %--------------------------------------------------%
!
      n     = 100
      nev   = 4
      ncv   = 20
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV2: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV2: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV2: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'I'
      which = 'LM'
      sigmar = 1.0E+0
      sigmai = 0.0E+0
!
!     %----------------------------------------------------%
!     | Construct C = A - SIGMA*I in real arithmetic, and  |
!     | factor C in real arithmetic using LAPACK           |
!     | subroutine sgttrf. The matrix A is chosen to be    |
!     | the tridiagonal matrix derived from standard       |
!     | central difference of the 1-d convection diffusion |
!     | operator u" + rho*u' on the interval [0, 1] with   |
!     | zero Dirichlet boundary condition.                 |
!     %----------------------------------------------------%
!
      rho = 1.0E+1
      h = one / real(n+1)
      s = rho*h / two
!
      s1 = -one-s
      s2 = two - sigmar
      s3 = -one+s
!
      do 10 j = 1, n-1
         dl(j) = s1
         dd(j) = s2
         du(j) = s3
  10  continue
      dd(n) = s2
!
      call sgttrf(n, dl, dd, du, du2, ipiv, ierr)
      if ( ierr .ne. 0 ) then
         print*, ' '
         print*, ' ERROR with _gttrf in _NDRV2.'
         print*, ' '
         go to 9000
      end if
!
!     %-----------------------------------------------------%
!     | The work array WORKL is used in SNAUPD as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.  The parameter TOL determines    |
!     | the stopping criterion. If TOL<=0, machine          |
!     | precision is used.  The variable IDO is used for    |
!     | reverse communication, and is initially set to 0.   |
!     | Setting INFO=0 indicates that a random vector is    |
!     | generated in SNAUPD to start the Arnoldi iteration. |
!     %-----------------------------------------------------%
!
      lworkl = 3*ncv**2+6*ncv
      tol    = zero
      ido    = 0
      info   = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed. Mode 3 of SNAUPD is used      |
!     | (IPARAM(7) = 3).  All these options can be        |
!     | changed by the user. For details see the          |
!     | documentation in SNAUPD.                          |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode   = 3

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
 20   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine SNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call snaupd ( ido, bmat, n, which, nev, tol, resid,
     &        ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,
     &        info )
!
         if ( ido .eq. -1 .or. ido .eq. 1) then
!
!           %-------------------------------------------%
!           | Perform  y <--- OP*x = inv[A-SIGMA*I]*x   |
!           | The user should supply his/her own linear |
!           | system solver here that takes             |
!           | workd(ipntr(1)) as the input, and returns |
!           | the result to workd(ipntr(2)).            |
!           %-------------------------------------------%
!
            call scopy( n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
!
            call sgttrs('N', n, 1, dl, dd, du, du2, ipiv,
     &                  workd(ipntr(2)), n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _gttrs in _NDRV2.'
               print*, ' '
               go to 9000
            end if
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call SNAUPD again. |
!           %-----------------------------------------%
!
            go to 20
!
         end if
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
!        | documentation in SNAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation in _naupd.'
         print *, ' '
!
      else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using SNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
         rvec = .true.
!
         call sneupd ( rvec, 'A', select, d, d(1,2), v, ldv,
     &        sigmar, sigmai, workev, bmat, n, which, nev, tol,
     &        resid, ncv, v, ldv, iparam, ipntr, workd,
     &        workl, lworkl, ierr )
!
!        %-----------------------------------------------%
!        | The real part of the eigenvalue is returned   |
!        | in the first column of the two dimensional    |
!        | array D, and the imaginary part is returned   |
!        | in the second column of D.  The corresponding |
!        | eigenvectors are returned in the first NEV    |
!        | columns of the two dimensional array V if     |
!        | requested.  Otherwise, an orthogonal basis    |
!        | for the invariant subspace corresponding to   |
!        | the eigenvalues in D is returned in V.        |
!        %-----------------------------------------------%
!
         if ( ierr .ne. 0 ) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of SNEUPD. |
!           %------------------------------------%
!
            print *, ' '
            print *, ' Error with _neupd, info = ', ierr
            print *, ' Check the documentation of _neupd. '
            print *, ' '
!
         else
!
            first  = .true.
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
               if (d(j,2) .eq. zero)  then
!
!                 %--------------------%
!                 | Ritz value is real |
!                 %--------------------%
!
                  call av(n, v(1,j), ax)
                  call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  d(j,3) = snrm2(n, ax, 1)
                  d(j,3) = d(j,3) / abs(d(j,1))
!
               else if (first) then
!
!                 %------------------------%
!                 | Ritz value is complex  |
!                 | Residual of one Ritz   |
!                 | value of the conjugate |
!                 | pair is computed.      |
!                 %------------------------%
!
                  call av(n, v(1,j), ax)
                  call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  call saxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                  d(j,3) = snrm2(n, ax, 1)
                  call av(n, v(1,j+1), ax)
                  call saxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                  call saxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                  d(j,3) = slapy2( d(j,3), snrm2(n, ax, 1) )
                  d(j+1,3) = d(j,3)
                  first = .false.
               else
                  first = .true.
               end if
!
 30         continue
!
!           %-----------------------------%
!           | Display computed residuals. |
!           %-----------------------------%
!
            call smout(6, nconv, 3, d, maxncv, -6,
     &           'Ritz values (Real,Imag) and relative residuals')
!
         end if
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' '
             print *, ' No shifts could be applied during implicit',
     &                ' Arnoldi update, try increasing NCV.'
             print *, ' '
         end if
!
         print *, ' '
         print *, ' _NDRV2 '
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
!     | Done with program sndrv2. |
!     %---------------------------%
!
 9000 continue
!
      end
!
!-------------------------------------------------------------------
!
!     matrix vector multiplication subroutine
!
      subroutine av (n, v, w)
      integer           n, j
      Real
     &                   v(n), w(n), rho, two, one, dd, dl, du, s, h
      common            /convct/ rho
      parameter         (one = 1.0E+0, two = 2.0E+0 )
!
!     Compute the matrix vector multiplication y<---A*x
!     where A is a n by n nonsymmetric tridiagonal matrix derived from
!     the central difference discretization of the 1-dimensional
!     convection diffusion operator on the interval [0,1] with
!     zero Dirichlet boundary condition.
!
!
      h = one / real(n+1)
      s = rho *h / two
      dd = two
      dl = -one - s
      du = -one + s
!
      w(1) =  dd*v(1) + du*v(2)
      do 10 j = 2,n-1
         w(j) = dl*v(j-1) + dd*v(j) + du*v(j+1)
 10   continue
      w(n) =  dl*v(n-1) + dd*v(n)
      return
      end
