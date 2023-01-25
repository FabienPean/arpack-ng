      program bug_142_gen
!
!     Avoid taking the initial vector in the range of OP after a restart
!     (generalized case)
!
!     Example program to illustrate the idea of reverse communication
!     for a standard nonsymmetric eigenvalue problem.
!
!     We implement example one of ex-nonsym.doc in DOCUMENTS directory
!
!\Example-1
!     ... Suppose we want to solve A*x = lambda*x in regular mode,
!         where A is obtained from the standard central difference
!         discretization of the convection-diffusion operator 
!                 (Laplacian u) + rho*(du / dx)
!         on the unit square [0,1]x[0,1] with zero Dirichlet boundary 
!         condition.
!
!     ... OP = A  and  B = I.
!
!     ... Assume "call av (nx,x,y)" computes y = A*x.c
!
!     ... Use mode 1 of DNAUPD.
!
!\BeginLib
!
!\Routines called:
!     dnaupd  ARPACK reverse communication interface routine.
!     dneupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     av      Matrix vector multiplication routine that computes A*x.
!     tv      Matrix vector multiplication routine that computes T*x, 
!             where T is a tridiagonal matrix.  It is used in routine
!             av.
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
! FILE: ndrv1.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!---------------------------------------------------------------------------
!
!     %-----------------------------%
!     | Define maximum dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!
      integer           maxn, maxnev, maxncv, ldv
      parameter         (maxn=256, maxnev=12, maxncv=30, ldv=maxn)
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      integer           iparam(11), ipntr(14)
      logical           select(maxncv)
      Double precision
     &                  ax(maxn), d(maxncv,3), resid(maxn), 
     &                  v(ldv,maxncv), workd(3*maxn), 
     &                  workev(3*maxncv), 
     &                  workl(3*maxncv*maxncv+6*maxncv), a(maxn, maxn)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character         bmat*1, which*2
      integer           ido, n, nx, nev, ncv, lworkl, info, j,
     &                  ierr, nconv, maxitr, ishfts, mode
      Double precision
     &                  tol, sigmar, sigmai
      logical           first, rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Double precision
     &                  zero
      parameter         (zero = 0.0D+0)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Double precision
     &                  dlapy2, dnrm2
      external          dlapy2, dnrm2, daxpy
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
!     %--------------------------------------------------%
!     | The number NX is the number of interior points   |
!     | in the discretization of the 2-dimensional       |
!     | convection-diffusion operator on the unit        |
!     | square with zero Dirichlet boundary condition.   | 
!     | The number N(=NX*NX) is the dimension of the     |
!     | matrix.  A standard eigenvalue problem is        |
!     | solved (BMAT = 'I').  NEV is the number of       |
!     | eigenvalues to be approximated.  The user can    |
!     | modify NX, NEV, NCV, WHICH to solve problems of  |
!     | different sizes, and to get different parts of   |
!     | the spectrum.  However, The following            |
!     | conditions must be satisfied:                    |
!     |                   N <= MAXN                      |
!     |                 NEV <= MAXNEV                    |
!     |           NEV + 2 <= NCV <= MAXNCV               | 
!     %--------------------------------------------------% 
!
      nx    = 10 
      n     = 11
      nev   = 1
      ncv   = 11 
      do i = 1,n
         do j = 1,n
            a(i,j) = 0.15d0/11
         end do
      end do
      do j = 2,n
         a(1,j) = a(1,j) + 0.85d0
      end do
      do i = 2,n
         a(i,1) = (1-a(1,1))/10
      end do

      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV1: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'G'
      which = 'LM'
!
!     %-----------------------------------------------------%
!     | The work array WORKL is used in DNAUPD as           |  
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.  The parameter TOL determines    |
!     | the stopping criterion. If TOL<=0, machine          |
!     | precision is used.  The variable IDO is used for    |
!     | reverse communication, and is initially set to 0.   |
!     | Setting INFO=0 indicates that a random vector is    |
!     | generated in DNAUPD to start the Arnoldi iteration. | 
!     %-----------------------------------------------------%
!
      lworkl  = 3*ncv**2+6*ncv 
      tol    = zero 
      ido    = 0
      info   = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | DNAUPD.                                           |
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
!        | Repeatedly call the routine DNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call dnaupd ( ido, bmat, n, which, nev, tol, resid, 
     &        ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, 
     &        info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- OP*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and return the matrix vector      |
!           | product to workd(ipntr(2)).               | 
!           %-------------------------------------------%
!
            call dgemv('N',n,n,1.0d0,a,maxn,
     &           workd(ipntr(1)),1,0.0d0,workd(ipntr(2)),1)
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DNAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         else if (ido .eq. 2) then

            call dcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)

            go to 10

         end if 
! 
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in DNAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation of _naupd'
         print *, ' '
         stop 1
!
      else 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
         rvec = .true.
!
         call dneupd ( rvec, 'A', select, d, d(1,2), v, ldv, 
     &        sigmar, sigmai, workev, bmat, n, which, nev, tol, 
     &        resid, ncv, v, ldv, iparam, ipntr, workd, workl,
     &        lworkl, ierr )
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
         if ( ierr .ne. 0) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DNEUPD. |
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
             nconv  = iparam(5)
             do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!
                if (d(j,2) .eq. zero)  then
!
!                  %--------------------%
!                  | Ritz value is real |
!                  %--------------------%
!
                   call dgemv('N',n,n,1.0d0,a,maxn,
     &                  v(1,j),1,0.0d0,ax,1)
                   call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                   d(j,3) = dnrm2(n, ax, 1)
                   d(j,3) = d(j,3) / abs(d(j,1))
!
                else if (first) then
!
!                  %------------------------%
!                  | Ritz value is complex. |
!                  | Residual of one Ritz   |
!                  | value of the conjugate |
!                  | pair is computed.      | 
!                  %------------------------%
!        
                   call dgemv('N',n,n,1.0d0,a,maxn,
     &                  v(1,j),1,0.0d0,ax,1)
                   call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                   call daxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                   d(j,3) = dnrm2(n, ax, 1)
                   call dgemv('N',n,n,1.0d0,a,maxn,
     &                  v(1,j+1),1,0.0d0,ax,1)
                   call daxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                   call daxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                   d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
                   d(j,3) = d(j,3) / dlapy2(d(j,1),d(j,2))
                   d(j+1,3) = d(j,3)
                   first = .false.
                else
                   first = .true.
                end if
!
 20          continue
!
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
!
             call dmout(6, nconv, 3, d, maxncv, -6,
     &            'Ritz values (Real,Imag) and relative residuals')
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
         print *, ' _NDRV1 '
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
!     | Done with program dndrv1. |
!     %---------------------------%
!
 9000 continue
!
      end
! 
