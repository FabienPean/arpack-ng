      program sndrv3
!
!     Simple program to illustrate the idea of reverse communication
!     in inverse mode for a generalized nonsymmetric eigenvalue problem.
!
!     We implement example three of ex-nonsym.doc in DOCUMENTS directory
!
!\Example-3
!     ... Suppose we want to solve A*x = lambda*B*x in inverse mode,
!         where A and B are derived from the finite element discretization
!         of the 1-dimensional convection-diffusion operator
!                           (d^2u / dx^2) + rho*(du/dx)
!         on the interval [0,1] with zero Dirichlet boundary condition
!         using linear elements.
!
!     ... So OP = inv[M]*A  and  B = M.
!
!     ... Use mode 2 of SNAUPD.
!
!\BeginLib
!
!\Routines called:
!     snaupd  ARPACK reverse communication interface routine.
!     sneupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     spttrf  LAPACK symmetric positive definite tridiagonal factorization
!             routine.
!     spttrs  LAPACK symmetric positive definite tridiagonal solve routine.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
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
! FILE: ndrv3.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
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
      integer           iparam(11), ipntr(14)
      logical           select(maxncv)
      Real
     &                  ax(maxn), mx(maxn), d(maxncv, 3), resid(maxn),
     &                  v(ldv,maxncv), workd(3*maxn),
     &                  workev(3*maxncv),
     &                  workl(3*maxncv*maxncv+6*maxncv),
     &                  md(maxn), me(maxn-1)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, ierr, j,
     &                  nconv, maxitr, ishfts, mode
      Real
     &                  tol, sigmar, sigmai, h
      logical           first, rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real
     &                  zero, one
      parameter         (zero = 0.0E+0, one = 1.0E+0)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
      Real
     &                  snrm2, slapy2
      external          saxpy, snrm2, spttrf, spttrs, slapy2
!
!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
      intrinsic         abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %----------------------------------------------------%
!     | The number N is the dimension of the matrix.  A    |
!     | generalized eigenvalue problem is solved (BMAT =   |
!     | 'G').  NEV is the number of eigenvalues to be      |
!     | approximated.  The user can modify NEV, NCV, WHICH |
!     | to solve problems of different sizes, and to get   |
!     | different parts of the spectrum.  However, The     |
!     | following conditions must be satisfied:            |
!     |                    N <= MAXN,                      |
!     |                  NEV <= MAXNEV,                    |
!     |              NEV + 2 <= NCV <= MAXNCV              |
!     %----------------------------------------------------%
!
      n     = 100
      nev   = 4
      ncv   = 20
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV3: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV3: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV3: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'G'
      which = 'LM'
!
!     %------------------------------------------------%
!     | M is the mass matrix formed by using piecewise |
!     | linear elements on [0,1].                      |
!     %------------------------------------------------%
!
      h = one / real(n+1)
      do 20 j = 1, n-1
         md(j) = 4.0E+0*h
         me(j) = one*h
  20  continue
      md(n) = 4.0E+0*h
!
      call spttrf(n, md, me, ierr)
      if ( ierr .ne. 0 ) then
         print*, ' '
         print*, ' ERROR with _pttrf. '
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
      tol    = 0.0
      ido    = 0
      info   = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 2 of SNAUPD is used     |
!     | (IPARAM(7) = 2).  All these options can be        |
!     | changed by the user. For details, see the         |
!     | documentation in SNAUPD.                          |
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
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine SNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call snaupd ( ido, bmat, n, which, nev, tol, resid,
     &                 ncv, v, ldv, iparam, ipntr, workd,
     &                 workl, lworkl, info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %----------------------------------------%
!           | Perform  y <--- OP*x = inv[M]*A*x      |
!           | The user should supply his/her own     |
!           | matrix vector routine and a linear     |
!           | system solver.  The matrix-vector      |
!           | subroutine should take workd(ipntr(1)) |
!           | as input, and the final result should  |
!           | be returned to workd(ipntr(2)).        |
!           %----------------------------------------%
!
            call av (n, workd(ipntr(1)), workd(ipntr(2)))
            call spttrs(n, 1, md, me, workd(ipntr(2)), n,
     &                  ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _pttrs. '
               print*, ' '
               go to 9000
            end if
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call SNAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         else if ( ido .eq. 2) then
!
!           %-------------------------------------%
!           |        Perform  y <--- M*x          |
!           | The matrix vector multiplication    |
!           | routine should take workd(ipntr(1)) |
!           | as input and return the result to   |
!           | workd(ipntr(2)).                    |
!           %-------------------------------------%
!
            call mv (n, workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call SNAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         end if
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
!        | Error message. Check the |
!        | documentation in SNAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation of _naupd.'
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
         call sneupd ( rvec, 'A', select, d, d(1,2), v, ldv,
     &        sigmar, sigmai, workev, bmat, n, which, nev, tol,
     &        resid, ncv, v, ldv, iparam, ipntr, workd,
     &        workl, lworkl, ierr )
!
!        %-----------------------------------------------%
!        | The real part of the eigenvalue is returned   |
!        | in the first column of the two dimensional    |
!        | array D, and the IMAGINARY part is returned   |
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
            print *, ' Check the documentation of _neupd'
            print *, ' '
!
         else
!
            first = .true.
            nconv = iparam(5)
            do 30 j=1, iparam(5)
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
               if (d(j,2) .eq. zero)  then
!
!                 %--------------------%
!                 | Ritz value is real |
!                 %--------------------%
!
                  call av(n, v(1,j), ax)
                  call mv(n, v(1,j), mx)
                  call saxpy(n, -d(j,1), mx, 1, ax, 1)
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
                  call mv(n, v(1,j), mx)
                  call saxpy(n, -d(j,1), mx, 1, ax, 1)
                  call mv(n, v(1,j+1), mx)
                  call saxpy(n, d(j,2), mx, 1, ax, 1)
                  d(j,3) = snrm2(n, ax, 1)**2
                  call av(n, v(1,j+1), ax)
                  call mv(n, v(1,j+1), mx)
                  call saxpy(n, -d(j,1), mx, 1, ax, 1)
                  call mv(n, v(1,j), mx)
                  call saxpy(n, -d(j,2), mx, 1, ax, 1)
                  d(j,3) = slapy2( d(j,3), snrm2(n, ax, 1) )
                  d(j,3) = d(j,3) / slapy2(d(j,1),d(j,2))
                  d(j+1,3) = d(j,3)
                  first = .false.
               else
                  first = .true.
               end if
!
  30        continue
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
         print *, ' _NDRV3 '
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
!     | Done with program sndrv3. |
!     %---------------------------%
!
 9000 continue
!
      end
!
!==========================================================================
!
!     matrix vector multiplication subroutine
!
      subroutine av (n, v, w)
      integer           n, j
      Real
     &                  v(n), w(n), one, two, dd, dl, du, s, h, rho
      parameter         ( rho = 1.0E+1, one = 1.0E+0,
     &                    two = 2.0E+0)
!
!     Compute the matrix vector multiplication y<---A*x
!     where A is stiffness matrix obtained from the finite element
!     discretization of the 1-dimensional convection diffusion operator
!                           d^2u/dx^2 + rho*(du/dx)
!     on the interval [0,1] with zero Dirichlet boundary condition using
!     linear elements.
!
      h = one / real(n+1)
      s = rho / two
      dd = two / h
      dl = -one/h - s
      du = -one/h + s
!
      w(1) =  dd*v(1) + du*v(2)
      do 10 j = 2,n-1
         w(j) = dl*v(j-1) + dd*v(j) + du*v(j+1)
 10   continue
      w(n) =  dl*v(n-1) + dd*v(n)
      return
      end
!------------------------------------------------------------------------
      subroutine mv (n, v, w)
      integer           n, j
      Real
     &                  v(n), w(n), one, four, h
      parameter         ( one = 1.0E+0, four = 4.0E+0)
!
!     Compute the matrix vector multiplication y<---M*x
!     where M is the mass matrix formed by using piecewise linear
!     elements on [0,1].
!
      w(1) =  four*v(1) + one*v(2)
      do 10 j = 2,n-1
         w(j) = one*v(j-1) + four*v(j) + one*v(j+1)
 10   continue
      w(n) =  one*v(n-1) + four*v(n)
!
      h = one / real(n+1)
      call sscal(n, h, w, 1)
      return
      end
