      program psndrv3
!
!     Message Passing Layer: MPI
!
!     Simple program to illustrate the idea of reverse communication
!     in inverse mode for a generalized nonsymmetric eigenvalue problem.
!
!     We implement example three of ex-nonsym.doc in DOCUMENTS directory
!
!\Example-3
!     ... Suppose we want to solve A*x = lambda*B*x in inverse mode,
!         where A is derived from the 1-dimensional convection-diffusion
!         operator on the interval [0,1] with zero boundary condition,
!         and M is the tridiagonal matrix with 4 on the diagonal and 1
!         on the subdiagonals.
!     ... So OP = inv[M]*A  and  B = M.
!     ... Use mode 2 of PSNAUPD.
!
!\BeginLib
!
!\Routines called:
!     psnaupd Parallel ARPACK reverse communication interface routine.
!     psneupd Parallel ARPACK routine that returns Ritz values and (optionally)
!              Ritz vectors.
!     spttrf  LAPACK symmetric positive definite tridiagonal factorization
!             routine.
!     spttrs  LAPACK symmetric positive definite tridiagonal solve routine.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     psnorm2 Parallel version of Level 1 BLAS that computes the norm of a vector.
!     av      Parallel Matrix vector multiplication routine that computes A*x.
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
!\Parallel Modifications
!     Kristi Maschhoff
!
!\Revision history:
!     Starting Point: Serial Code FILE: ndrv3.F   SID: 2.2
!
!\SCCS Information:
! FILE: ndrv3.F   SID: 1.1   DATE OF SID: 8/13/96   RELEASE: 1
!
!\Remarks
!     1. None
!
!\EndLib
!--------------------------------------------------------------------------
!
      include 'mpif.h'
      include 'debug.h'
      include 'stat.h'

!     %-------------------------------%
!     | MPI INTERFACE                 |
!     | ILP64 is not supported by MPI |
!     | integer*4 must be imposed in  |
!     | all calls involving MPI.      |
!     |                               |
!     | Use ierr for MPI calls.       |
!     %-------------------------------%

      integer*4         comm, myid, nprocs, rc, ierr

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
     &                  md(maxn), me(maxn-1), temp(maxn), temp_buf(maxn)
!
!     %------------------------------------%
!     | Local Scalars                      |
!     |                                    |
!     | Use info if ILP64 can be supported |
!     | (call to BLAS, LAPACK, ARPACK).    |
!     %------------------------------------%
!
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, nloc, j,
     &                  nconv, maxitr, ishfts, mode, blk
      Real
     &                  tol, sigmar, sigmai
      logical           first, rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real
     &                  zero, one
      parameter         (zero = 0.0, one = 1.0)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
      Real
     &                  psnorm2, slapy2
      external          saxpy, psnorm2, spttrf, spttrs, slapy2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      call MPI_INIT( ierr )
      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
!
      ndigit = -3
      logfil = 6
      mnaupd = 1
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
!
!     %--------------------------------------%
!     | Set up distribution of data to nodes |
!     %--------------------------------------%
!
      nloc = (n / nprocs )
      blk = nloc
      if ( mod(n, nprocs) .gt. 0 ) then
         if ( myid .eq. nprocs-1 ) nloc = nloc + mod(n, nprocs)
*      if ( mod(n, nprocs) .gt. myid ) nloc = nloc + 1
      endif
!
      if ( nloc .gt. maxn ) then
         print *, ' ERROR with _NDRV3: NLOC is greater than MAXN '
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
!     %-----------------------------------------------------%
!     | The matrix M is chosen to be the symmetric tri-     |
!     | diagonal matrix with 4 on the diagonal and 1 on the |
!     | off diagonals. It is factored by LAPACK subroutine  |
!     | spttrf.                                             |
!     %-----------------------------------------------------%
!
      do 20 j = 1, n-1
         md(j) = 4.0
         me(j) = one
  20  continue
      md(n) = 4.0*one
!
      call spttrf(n, md, me, info)
      if ( info .ne. 0 ) then
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
         call psnaupd( comm, ido, bmat, nloc, which, nev, tol, resid,
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
            call av (comm, nloc, n, workd(ipntr(1)), workd(ipntr(2)))
!======== Hack for Linear system ======= ccc
            call sscal(n, zero, temp, 1)
            call sscal(n, zero, temp_buf, 1)
            do 15 j=1,nloc
               temp_buf(myid*blk + j) = workd(ipntr(2) + j - 1)
   15       continue
            call MPI_ALLREDUCE( temp_buf, temp, n,
     &            MPI_REAL, MPI_SUM, comm, ierr )
            call spttrs(n, 1, md, me, temp, n,
     &                  info)
            if ( info .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _pttrs. '
               print*, ' '
               go to 9000
            end if
            do 16 j=1,nloc
               workd(ipntr(2) + j - 1 ) = temp(myid*blk + j)
   16       continue
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
            call mv (comm, nloc, workd(ipntr(1)), workd(ipntr(2)))
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
!        %---------------------------%
!        | Error message. Check the  |
!        | documentation in PSNAUPD. |
!        %---------------------------%
!
         if ( myid .eq. 0 ) then
            print *, ' '
            print *, ' Error with _naupd, info = ', info
            print *, ' Check the documentation of _naupd.'
            print *, ' '
         endif
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
         call psneupd ( comm, rvec, 'A', select, d, d(1,2), v, ldv,
     &        sigmar, sigmai, workev, bmat, nloc, which, nev, tol,
     &        resid, ncv, v, ldv, iparam, ipntr, workd,
     &        workl, lworkl, info )
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
         if ( info .ne. 0 ) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of SNEUPD. |
!           %------------------------------------%
!
            if ( myid .eq. 0 ) then
               print *, ' '
               print *, ' Error with _neupd, info = ', info
               print *, ' Check the documentation of _neupd'
               print *, ' '
            endif
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
                  call av(comm, nloc, n, v(1,j), ax)
                  call mv(comm, nloc, v(1,j), mx)
                  call saxpy(nloc, -d(j,1), mx, 1, ax, 1)
                  d(j,3) = psnorm2(comm, nloc, ax, 1)
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
                  call av(comm, nloc, n, v(1,j), ax)
                  call mv(comm, nloc, v(1,j), mx)
                  call saxpy(nloc, -d(j,1), mx, 1, ax, 1)
                  call mv(comm, nloc, v(1,j+1), mx)
                  call saxpy(nloc, d(j,2), mx, 1, ax, 1)
                  d(j,3) = psnorm2(comm, nloc, ax, 1)**2
                  call av(comm, nloc, n, v(1,j+1), ax)
                  call mv(comm, nloc, v(1,j+1), mx)
                  call saxpy(nloc, -d(j,1), mx, 1, ax, 1)
                  call mv(comm, nloc, v(1,j), mx)
                  call saxpy(nloc, -d(j,2), mx, 1, ax, 1)
                  d(j,3) = slapy2( d(j,3), psnorm2(comm,nloc,ax,1) )
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
            call psmout(comm, 6, nconv, 3, d, maxncv, -6,
     &           'Ritz values (Real,Imag) and direct residuals')
!
         end if
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
         if (myid .eq. 0)then
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' '
            print *, ' No shifts could be applied during implicit
     &                 Arnoldi update, try increasing NCV.'
            print *, ' '
         end if
!
         print *, ' '
         print *, '_NDRV3 '
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of processors is ', nprocs
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
         endif
      end if
!
!     %----------------------------%
!     | Done with program psndrv3. |
!     %----------------------------%
!
 9000 continue
!
!     %-------------------------%
!     | Release resources MPI |
!     %-------------------------%
!
      call MPI_FINALIZE(rc)
!
      end
!
!==========================================================================
!
!     parallel matrix vector multiplication subroutine
!
!     Compute the matrix vector multiplication y<---A*x
!     where A is a n by n nonsymmetric tridiagonal matrix derived
!     from the central difference discretization of the 1-dimensional
!     convection diffusion operator on the interval [0,1] with
!     zero Dirichlet boundary condition.
!
      subroutine av (comm, nloc, n, v, w)
!
!     .. MPI Declarations ...
      include           'mpif.h'
      integer*4         comm, nprocs, myid, ierr,
     &                  status(MPI_STATUS_SIZE)
!
      integer           nloc, n, j, next, prev
      Real
     &                  v(nloc), w(nloc), one, two, dd, dl, du,
     &                  s, h, rho, mv_buf
      parameter         ( rho = 10.0, one = 1.0,
     &                    two = 2.0)
!
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
      h = one / real(n+1)
      s = rho*h / two
      dd = two
      dl = -one - s
      du = -one + s
!
      w(1) =  dd*v(1) + du*v(2)
      do 10 j = 2,nloc-1
         w(j) = dl*v(j-1) + dd*v(j) + du*v(j+1)
 10   continue
      w(nloc) =  dl*v(nloc-1) + dd*v(nloc)
!
      next = myid + 1
      prev = myid - 1
      if ( myid .lt. nprocs-1 ) then
         call mpi_send( v(nloc), 1, MPI_REAL,
     &                  next, myid+1, comm, ierr )
      endif
      if ( myid .gt. 0 ) then
         call mpi_recv( mv_buf, 1, MPI_REAL, prev, myid,
     &                  comm, status, ierr )
         w(1) = w(1) + dl*mv_buf
      endif
!
      if ( myid .gt. 0 ) then
         call mpi_send( v(1), 1, MPI_REAL,
     &                  prev, myid-1, comm, ierr )
      endif
      if ( myid .lt. nprocs-1 ) then
         call mpi_recv( mv_buf, 1, MPI_REAL, next, myid,
     &                  comm, status, ierr )
         w(nloc) = w(nloc) + du*mv_buf
      endif
!
      return
      end
!------------------------------------------------------------------------
!
!     Compute the matrix vector multiplication y<---M*x
!     where M is a n by n tridiagonal matrix with 4 on the
!     diagonal, 1 on the subdiagonal and the superdiagonal.
!
      subroutine mv (comm, nloc, v, w)
!
!     .. MPI Declarations ...
      include           'mpif.h'
      integer*4         comm, nprocs, myid, ierr,
     &                  status(MPI_STATUS_SIZE)
!
      integer           nloc, j, next, prev
      Real
     &                  v(nloc), w(nloc), one, four, mv_buf
      parameter         ( one = 1.0, four = 4.0)
!
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
!
      w(1) =  four*v(1) + one*v(2)
      do 10 j = 2,nloc-1
         w(j) = one*v(j-1) + four*v(j) + one*v(j+1)
 10   continue
      w(nloc) =  one*v(nloc-1) + four*v(nloc)
!
      next = myid + 1
      prev = myid - 1
      if ( myid .lt. nprocs-1 ) then
         call mpi_send( v(nloc), 1, MPI_REAL,
     &                  next, myid+1, comm, ierr )
      endif
      if ( myid .gt. 0 ) then
         call mpi_recv( mv_buf, 1, MPI_REAL, prev, myid,
     &                  comm, status, ierr )
         w(1) = w(1) + mv_buf
      endif
!
      if ( myid .gt. 0 ) then
         call mpi_send( v(1), 1, MPI_REAL,
     &                  prev, myid-1, comm, ierr )
      endif
      if ( myid .lt. nprocs-1 ) then
         call mpi_recv( mv_buf, 1, MPI_REAL, next, myid,
     &                  comm, status, ierr )
         w(nloc) = w(nloc) + mv_buf
      endif
!
      return
      end
!------------------------------------------------------------
      subroutine mv2 (comm, n, v, w)
      integer           n, j, comm
      Real
     &                  v(n), w(n)
      do 10 j=1,n
         w(j) = v(j)
 10   continue
!
      return
      end
