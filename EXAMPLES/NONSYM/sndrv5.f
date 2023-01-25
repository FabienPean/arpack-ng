      program sndrv5
!
!     Simple program to illustrate the idea of reverse communication
!     in shift-invert mode for a generalized nonsymmetric eigenvalue problem.
!
!     We implement example five of ex-nonsym.doc in DOCUMENTS directory
!
!\Example-5
!
!     ... Suppose we want to solve A*x = lambda*B*x in shift-invert mode
!         The matrix A is the tridiagonal matrix with 2 on the diagonal,
!         -2 on the subdiagonal and 3 on the superdiagonal.  The matrix M
!         is the tridiagonal matrix with 4 on the diagonal and 1 on the
!         off-diagonals.
!     ... The shift sigma is a complex number (sigmar, sigmai).
!     ... OP = Real_Part{inv[A-(SIGMAR,SIGMAI)*M]*M and  B = M.
!     ... Use mode 3 of SNAUPD.
!
!\BeginLib
!
!\Routines called:
!     snaupd  ARPACK reverse communication interface routine.
!     sneupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     cgttrf  LAPACK complex matrix factorization routine.
!     cgttrs  LAPACK complex linear system solve routine.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     sdot    Level 1 BLAS that computes the dot product of two vectors.
!     snrm2   Level 1 BLAS that computes the norm of a vector
!     av      Matrix vector subroutine that computes A*x.
!     mv      Matrix vector subroutine that computes M*x.
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
! FILE: ndrv5.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!-------------------------------------------------------------------------
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
     &                  ax(maxn), mx(maxn), d(maxncv,3), resid(maxn),
     &                  v(ldv,maxncv), workd(3*maxn),
     &                  workev(3*maxncv),
     &                  workl(3*maxncv*maxncv+6*maxncv)
      Complex
     &                  cdd(maxn), cdl(maxn), cdu(maxn),
     &                  cdu2(maxn), ctemp(maxn)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, ierr, j,
     &                  nconv, maxitr, ishfts, mode
      Real
     &                  tol, numr, numi, denr, deni, sigmar, sigmai
      Complex
     &                  c1, c2, c3
      logical           first, rvec
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      external          cgttrf, cgttrs
      Real
     &                  sdot, snrm2, slapy2
      external          sdot, snrm2, slapy2
!
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real
     &                   zero
      parameter         (zero = 0.0E+0)
!
!     %--------------------%
!     | Intrinsic Function |
!     %--------------------%
!
      intrinsic          real, cmplx, abs
!
!     %-----------------------%
!     | Executable statements |
!     %-----------------------%
!
!     %----------------------------------------------------%
!     | The number N is the dimension of the matrix.  A    |
!     | generalized eigenvalue problem is solved (BMAT =   |
!     | 'G').  NEV is the number of eigenvalues (closest   |
!     | to the shift (SIGMAR,SIGMAI)) to be approximated.  |
!     | Since the shift-invert mode is used, WHICH is set  |
!     | to 'LM'.  The user can modify NEV, NCV, SIGMAR,    |
!     | SIGMAI to solve problems of different sizes, and   |
!     | to get different parts of the spectrum. However,   |
!     | The following conditions must be satisfied:        |
!     |                     N <= MAXN,                     |
!     |                   NEV <= MAXNEV,                   |
!     |               NEV + 2 <= NCV <= MAXNCV             |
!     %----------------------------------------------------%
!
      n     = 100
      nev   = 4
      ncv   = 20
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV5: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV5: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV5: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'G'
      which = 'LM'
      sigmar = 4.0E-1
      sigmai = 6.0E-1
!
!     %---------------------------------------------------%
!     | Construct C = A - (SIGMAR,SIGMAI)*M in complex    |
!     | arithmetic, and factor C in complex arithmetic    |
!     | (using LAPACK subroutine cgttrf). The matrix A is |
!     | chosen to be the tridiagonal matrix with -2 on    |
!     | the subdiagonal, 2 on the diagonal and 3 on the   |
!     | superdiagonal. The matrix M is chosen to the      |
!     | symmetric tridiagonal matrix with 4 on the        |
!     | diagonal and 1 on the off-diagonals.              |
!     %---------------------------------------------------%
!
      c1 = cmplx(-2.0E+0-sigmar, -sigmai)
      c2 = cmplx( 2.0E+0-4.0E+0*sigmar, -4.0E+0*sigmai)
      c3 = cmplx( 3.0E+0-sigmar, -sigmai)
!
      do 10 j = 1, n-1
         cdl(j) = c1
         cdd(j) = c2
         cdu(j) = c3
  10  continue
      cdd(n) = c2
!
      call cgttrf(n, cdl, cdd, cdu, cdu2, ipiv, ierr)
      if ( ierr .ne. 0 ) then
         print*, ' '
         print*, ' ERROR with _gttrf in _NDRV5.'
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
!     | This program uses exact shift with respect to     |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 3 of SNAUPD is used     |
!     | (IPARAM(7) = 3).  All these options can be        |
!     | changed by the user. For details, see the         |
!     | documentation in SNAUPD.                          |
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
         if (ido .eq. -1) then
!
!           %-------------------------------------------------------%
!           |                       Perform                         |
!           | y <--- OP*x = Real_Part{inv[A-(SIGMAR,SIGMAI)*M]*M*x} |
!           | to force starting vector into the range of OP. The    |
!           | user should supply his/her own matrix vector          |
!           | multiplication routine and a complex linear system    |
!           | solver.  The matrix vector multiplication routine     |
!           | should take workd(ipntr(1)) as the input. The final   |
!           | result (a real vector) should be returned to          |
!           | workd(ipntr(2)).                                      |
!           %-------------------------------------------------------%
!
            call mv (n, workd(ipntr(1)), workd(ipntr(2)))
            do 30 j = 1, n
               ctemp(j) = cmplx(workd(ipntr(2)+j-1))
  30        continue
!
            call cgttrs('N', n, 1, cdl, cdd, cdu, cdu2, ipiv,
     &                  ctemp, n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _gttrs in _NDRV5.'
               print*, ' '
               go to 9000
            end if
            do  40 j = 1, n
               workd(ipntr(2)+j-1) = real(ctemp(j))
  40        continue
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call SNAUPD again. |
!           %-----------------------------------------%
!
            go to 20
!
         else if ( ido .eq. 1) then
!
!           %-------------------------------------------------------%
!           |                        Perform                        |
!           | y <--- OP*x = Real_Part{inv[A-(SIGMAR,SIGMAI)*M]*M*x} |
!           | M*x has been saved in workd(ipntr(3)). The user only  |
!           | needs the complex linear system solver here that      |
!           | takes complex[workd(ipntr(3))] as input, and returns  |
!           | the result to workd(ipntr(2)).                        |
!           %-------------------------------------------------------%
!
            do 50 j = 1,n
               ctemp(j) = cmplx(workd(ipntr(3)+j-1))
  50        continue
            call cgttrs ('N', n, 1, cdl, cdd, cdu, cdu2, ipiv,
     &                   ctemp, n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _gttrs in _NDRV5.'
               print*, ' '
               go to 9000
            end if
            do  60 j = 1, n
               workd(ipntr(2)+j-1) = real(ctemp(j))
  60        continue
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call SNAUPD again. |
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
!           | L O O P   B A C K to call SNAUPD again. |
!           %-----------------------------------------%
!
            go to 20
!
         end if
!
!
!     %------------------------------------------%
!     | Either we have convergence, or there is  |
!     | an error.                                |
!     %------------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in SNAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _naupd info = ',info
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
         if ( ierr .ne. 0) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of SNEUPD. |
!           %------------------------------------%
!
            print *, ' '
            print *, ' Error with _neupd = ', ierr
            print *, ' Check the documentation of _neupd. '
            print *, ' '
!
         else
!
            first = .true.
            nconv =  iparam(5)
            do 70 j=1,nconv
!
!              %-------------------------------------%
!              | Use Rayleigh Quotient to recover    |
!              | eigenvalues of the original problem.|
!              %-------------------------------------%
!
               if ( d(j,2) .eq. zero ) then
!
!                 %---------------------------%
!                 |    Eigenvalue is real.    |
!                 | Compute d = x'(Ax)/x'(Mx).|
!                 %---------------------------%
!
                  call av(n, v(1,j), ax )
                  numr = sdot(n, v(1,j), 1, ax, 1)
                  call mv(n, v(1,j), ax )
                  denr = sdot(n, v(1,j), 1, ax, 1)
                  d(j,1) =  numr / denr
!
               else if (first) then
!
!                 %------------------------%
!                 | Eigenvalue is complex. |
!                 | Compute the first one  |
!                 | of the conjugate pair. |
!                 %------------------------%
!
!                 %----------------%
!                 | Compute x'(Ax) |
!                 %----------------%

                  call av(n, v(1,j), ax )
                  numr = sdot(n, v(1,j), 1, ax, 1)
                  numi = sdot(n, v(1,j+1), 1, ax, 1)
                  call av(n, v(1,j+1), ax)
                  numr = numr + sdot(n,v(1,j+1),1,ax,1)
                  numi = -numi + sdot(n,v(1,j),1,ax,1)
!
!                 %----------------%
!                 | Compute x'(Mx) |
!                 %----------------%
!
                  call mv(n, v(1,j), ax )
                  denr = sdot(n, v(1,j), 1, ax, 1)
                  deni = sdot(n, v(1,j+1), 1, ax, 1)
                  call mv(n, v(1,j+1), ax)
                  denr = denr + sdot(n,v(1,j+1),1,ax,1)
                  deni = -deni + sdot(n,v(1,j),1, ax,1)
!
!                 %----------------%
!                 | d=x'(Ax)/x'(Mx)|
!                 %----------------%
!
                  d(j,1) = (numr*denr+numi*deni) /
     &                      slapy2(denr, deni)
                  d(j,2) = (numi*denr-numr*deni) /
     &                      slapy2(denr, deni)
                  first = .false.
!
               else
!
!                 %------------------------------%
!                 | Get the second eigenvalue of |
!                 | the conjugate pair by taking |
!                 | the conjugate of the last    |
!                 | eigenvalue computed.         |
!                 %------------------------------%
!
                  d(j,1) = d(j-1,1)
                  d(j,2) = -d(j-1,2)
                  first = .true.
!
               end if
!
  70        continue
!
!           %---------------------------%
!           | Compute the residual norm |
!           |                           |
!           |   ||  A*x - lambda*x ||   |
!           |                           |
!           | for the NCONV accurately  |
!           | computed eigenvalues and  |
!           | eigenvectors.  (iparam(5) |
!           | indicates how many are    |
!           | accurate to the requested |
!           | tolerance)                |
!           %---------------------------%
!
            first  = .true.
            do 80 j=1, nconv
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
                  d(j,3) = snrm2(n, ax, 1)
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
  80        continue
!
!
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
!
             call smout(6, nconv, 3, d, maxncv, -6,
     &            'Ritz values (Real,Imag) and relative residuals')
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
         print *, ' _NDRV5 '
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
!     | Done with program sndrv5. |
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
      subroutine mv (n, v, w)
      integer           n, j
      Real
     &                  v(n), w(n), one, four
      parameter         (one = 1.0E+0, four = 4.0E+0)
!
!     Compute the matrix vector multiplication y<---M*x
!     where M is a n by n symmetric tridiagonal matrix with 4 on the
!     diagonal, 1 on the subdiagonal and superdiagonal.
!
      w(1) =  four*v(1) + one*v(2)
      do 10 j = 2,n-1
         w(j) = one*v(j-1) + four*v(j) + one*v(j+1)
 10   continue
      w(n) =  one*v(n-1) + four*v(n)
      return
      end
!------------------------------------------------------------------
      subroutine av (n, v, w)
      integer           n, j
      Real
     &                  v(n), w(n), three, two
      parameter         (three = 3.0E+0, two = 2.0E+0)
!
!     Compute the matrix vector multiplication y<---A*x
!     where M is a n by n symmetric tridiagonal matrix with 2 on the
!     diagonal, -2 on the subdiagonal and 3 on the superdiagonal.
!
      w(1) =  two*v(1) + three*v(2)
      do 10 j = 2,n-1
         w(j) = -two*v(j-1) + two*v(j) + three*v(j+1)
 10   continue
      w(n) =  -two*v(n-1) + two*v(n)
      return
      end

