      program cndrv4
!
!     Simple program to illustrate the idea of reverse communication
!     in shift and invert mode for a generalized complex nonsymmetric
!     eigenvalue problem.
!
!     We implement example four of ex-complex.doc in DOCUMENTS directory
!
!\Example-4
!     ... Suppose we want to solve A*x = lambda*B*x in shift-invert mode,
!         where A and B are derived from a finite element discretization
!         of a 1-dimensional convection-diffusion operator
!                         (d^2u/dx^2) + rho*(du/dx)
!         on the interval [0,1] with zero boundary condition using
!         piecewise linear elements.
!
!     ... where the shift sigma is a complex number.
!
!     ... OP = inv[A-SIGMA*M]*M  and  B = M.
!
!     ... Use mode 3 of CNAUPD.
!
!\BeginLib
!
!\Routines called:
!     cnaupd  ARPACK reverse communication interface routine.
!     cneupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     cgttrf  LAPACK tridiagonal factorization routine.
!     cgttrs  LAPACK tridiagonal solve routine.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     caxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     ccopy   Level 1 BLAS that copies one vector to another.
!     scnrm2  Level 1 BLAS that computes the norm of a complex vector.
!     av      Matrix vector multiplication routine that computes A*x.
!     mv      Matrix vector multiplication routine that computes M*x.
!
!\Author
!     Danny Sorensen
!     Richard Lehoucq
!     Chao Yang
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: ndrv4.F   SID: 2.4   DATE OF SID: 10/18/00   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!-----------------------------------------------------------------------
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
      Complex
     &                  ax(maxn), mx(maxn), d(maxncv),
     &                  v(ldv,maxncv), workd(3*maxn), resid(maxn),
     &                  workev(2*maxncv),
     &                  workl(3*maxncv*maxncv+5*maxncv),
     &                  dd(maxn), dl(maxn), du(maxn),
     &                  du2(maxn)
      Real
     &                  rwork(maxn), rd(maxncv,3)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, j, ierr,
     &                  nconv, maxitr, ishfts, mode
      Complex
     &                  rho, h, s,
     &                  sigma, s1, s2, s3
      common            /convct/ rho
!
      Real
     &                  tol
      logical           rvec
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Real
     &                  scnrm2, slapy2
      external          scnrm2, caxpy, ccopy, cgttrf, cgttrs,
     &                  slapy2
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Complex
     &                   one, zero, two, four, six
      parameter         (one = (1.0E+0, 0.0E+0) ,
     &                   zero = (0.0E+0, 0.0E+0) ,
     &                   two = (2.0E+0, 0.0E+0) ,
     &                   four = (4.0E+0, 0.0E+0) ,
     &                   six = (6.0E+0, 0.0E+0) )
!
!     %-----------------------%
!     | Executable statements |
!     %-----------------------%
!
!     %----------------------------------------------------%
!     | The number N is the dimension of the matrix.  A    |
!     | generalized eigenvalue problem is solved (BMAT =   |
!     | 'G').  NEV is the number of eigenvalues (closest   |
!     | to SIGMAR) to be approximated.  Since the          |
!     | shift-invert mode is used,  WHICH is set to 'LM'.  |
!     | The user can modify NEV, NCV, SIGMA to solve       |
!     | problems of different sizes, and to get different  |
!     | parts of the spectrum.  However, The following     |
!     | conditions must be satisfied:                      |
!     |                     N <= MAXN,                     |
!     |                   NEV <= MAXNEV,                   |
!     |               NEV + 2 <= NCV <= MAXNCV             |
!     %----------------------------------------------------%
!
      n     = 100
      nev   = 4
      ncv   = 20
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV4: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV4: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV4: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'G'
      which = 'LM'
      sigma = one
!
!     %--------------------------------------------------%
!     | Construct C = A - SIGMA*M in COMPLEX arithmetic. |
!     | Factor C in COMPLEX arithmetic (using LAPACK     |
!     | subroutine cgttrf). The matrix A is chosen to be |
!     | the tridiagonal matrix derived from the standard |
!     | central difference discretization of the 1-d     |
!     | convection-diffusion operator u``+ rho*u` on the |
!     | interval [0, 1] with zero Dirichlet boundary     |
!     | condition.  The matrix M is chosen to be the     |
!     | symmetric tridiagonal matrix with 4.0 on the     |
!     | diagonal and 1.0 on the off-diagonals.           |
!     %--------------------------------------------------%
!
      rho = (1.0E+1, 0.0E+0)
      h = one / cmplx(n+1)
      s = rho / two
!
      s1 = -one/h - s - sigma*h/six
      s2 = two/h  - four*sigma*h/six
      s3 = -one/h + s - sigma*h/six
!
      do 10 j = 1, n-1
	 dl(j) = s1
	 dd(j) = s2
	 du(j) = s3
  10  continue
      dd(n) = s2
!
      call cgttrf(n, dl, dd, du, du2, ipiv, ierr)
      if ( ierr .ne. 0 ) then
         print*, ' '
         print*, ' ERROR with _gttrf in _NDRV4.'
         print*, ' '
         go to 9000
      end if
!
!     %-----------------------------------------------------%
!     | The work array WORKL is used in CNAUPD as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.  The parameter TOL determines    |
!     | the stopping criterion. If TOL<=0, machine          |
!     | precision is used.  The variable IDO is used for    |
!     | reverse communication, and is initially set to 0.   |
!     | Setting INFO=0 indicates that a random vector is    |
!     | generated in CNAUPD to start the Arnoldi iteration. |
!     %-----------------------------------------------------%
!
      lworkl = 3*ncv**2+5*ncv
      tol    = 0.0
      ido    = 0
      info   = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed. Mode 3 of CNAUPD is used      |
!     | (IPARAM(7) = 3).  All these options can be        |
!     | changed by the user. For details see the          |
!     | documentation in CNAUPD.                          |
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
!     %------------------------------------------%
!     | M A I N   L O O P(Reverse communication) |
!     %------------------------------------------%
!
 20   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine CNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call cnaupd ( ido, bmat, n, which, nev, tol, resid, ncv,
     &        v, ldv, iparam, ipntr, workd, workl, lworkl,
     &        rwork, info )

!
         if (ido .eq. -1) then
!
!           %-------------------------------------------%
!           | Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x |
!           | to force starting vector into the range   |
!           | of OP.   The user should supply his/her   |
!           | own matrix vector multiplication routine  |
!           | and a linear system solver.  The matrix   |
!           | vector multiplication routine should take |
!           | workd(ipntr(1)) as the input. The final   |
!           | result should be returned to              |
!           | workd(ipntr(2)).                          |
!           %-------------------------------------------%
!
            call mv (n, workd(ipntr(1)), workd(ipntr(2)))
            call cgttrs('N', n, 1, dl, dd, du, du2, ipiv,
     &                  workd(ipntr(2)), n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _gttrs in _NDRV4.'
               print*, ' '
               go to 9000
            end if
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call CNAUPD again. |
!           %-----------------------------------------%
!
            go to 20
!
         else if ( ido .eq. 1) then
!
!           %-----------------------------------------%
!           | Perform y <-- OP*x = inv[A-sigma*M]*M*x |
!           | M*x has been saved in workd(ipntr(3)).  |
!           | The user only need the linear system    |
!           | solver here that takes workd(ipntr(3))  |
!           | as input, and returns the result to     |
!           | workd(ipntr(2)).                        |
!           %-----------------------------------------%
!
            call ccopy( n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
            call cgttrs ('N', n, 1, dl, dd, du, du2, ipiv,
     &                   workd(ipntr(2)), n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _gttrs in _NDRV4.'
               print*, ' '
               go to 9000
            end if
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call CNAUPD again. |
!           %-----------------------------------------%
!
            go to 20
!
         else if ( ido .eq. 2) then
!
!           %---------------------------------------------%
!           |          Perform  y <--- M*x                |
!           | Need matrix vector multiplication routine   |
!           | here that takes workd(ipntr(1)) as input    |
!           | and returns the result to workd(ipntr(2)).  |
!           %---------------------------------------------%
!
  	    call mv (n, workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call CNAUPD again. |
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
!        %----------------------------%
!        |  Error message, check the  |
!        |  documentation in CNAUPD   |
!        %----------------------------%
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
!        | Post-Process using CNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
         rvec = .true.
!
         call cneupd (rvec, 'A', select, d, v, ldv, sigma,
     &        workev, bmat, n, which, nev, tol, resid, ncv, v,
     &        ldv, iparam, ipntr, workd, workl, lworkl, rwork,
     &        ierr)
!
!        %----------------------------------------------%
!        | Eigenvalues are returned in the one          |
!        | dimensional array D.  The corresponding      |
!        | eigenvectors are returned in the first NCONV |
!        | (=IPARAM(5)) columns of the two dimensional  |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         if ( ierr .ne. 0) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of CNEUPD. |
!           %------------------------------------%
!
             print *, ' '
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
!
         else
!
             nconv = iparam(5)
             do 80 j=1, nconv
!
                call av(n, v(1,j), ax)
                call mv(n, v(1,j), mx)
                call caxpy(n, -d(j), mx, 1, ax, 1)
                rd(j,1) = real (d(j))
                rd(j,2) = aimag(d(j))
                rd(j,3) = scnrm2(n, ax, 1)
                rd(j,3) = rd(j,3) / slapy2(rd(j,1),rd(j,2))
  80         continue
!
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
!
             call smout(6, nconv, 3, rd, maxncv, -6,
     &            'Ritz values (Real, Imag) and direct residuals')
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
         print *, '_NDRV4 '
         print *, '====== '
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
 9000 continue
!
      end
!
!==========================================================================
!
!     matrix vector multiplication subroutine
!
      subroutine mv (n, v, w)
      integer           n, j
      Complex
     &                  v(n), w(n), one, four, six, h
      parameter         (one = (1.0E+0, 0.0E+0) ,
     &                   four = (4.0E+0, 0.0E+0) ,
     &                   six = (6.0E+0, 0.0E+0) )
!
!     Compute the matrix vector multiplication y<---M*x
!     where M is a n by n symmetric tridiagonal matrix with 4 on the
!     diagonal, 1 on the subdiagonal and superdiagonal.
!
      w(1) =  ( four*v(1) + one*v(2) ) / six
      do 40 j = 2,n-1
         w(j) = ( one*v(j-1) + four*v(j) + one*v(j+1) ) / six
 40   continue
      w(n) =  ( one*v(n-1) + four*v(n) ) / six
!
      h = one / cmplx(n+1)
      call cscal(n, h, w, 1)
      return
      end
!------------------------------------------------------------------
      subroutine av (n, v, w)
      integer           n, j
      Complex
     &                  v(n), w(n), one, two, dd, dl, du, s, h, rho
      parameter         (one = (1.0E+0, 0.0E+0) ,
     &                   two = (2.0E+0, 0.0E+0) )
      common            /convct/ rho
!
      h = one / cmplx(n+1)
      s = rho / two
      dd = two / h
      dl = -one/h - s
      du = -one/h + s
!
      w(1) =  dd*v(1) + du*v(2)
      do 40 j = 2,n-1
         w(j) = dl*v(j-1) + dd*v(j) + du*v(j+1)
 40   continue
      w(n) =  dl*v(n-1) + dd*v(n)
      return
      end

