
!Calculates matrix element for i-th row and j-th column.
pure function Aelem(i, j, param) Result(res)
  USE ModMtrx
  Integer(4), intent(in) :: i, j !row and column indices
  Type(Mtrx), intent(in) :: param !Arbitrary array of parameters
  Double precision :: res !Output value
  !Just uses parameters array as the input matrix
  res = param%d(i,j)
end

subroutine ExampleA()
  USE ModAppr
  Type(Mtrx) U, V, A, C, CA, AR, E, Ahat, A1, AB, ABin, R, Q
  Type(IntVec) per1, per2, per3, peri
  Type(Vector) cin, S
  Integer(4) n, k, maxsteps, maxswaps, swapsmade, i
  Double precision dsecnd, time
  Logical dp
  
  !WELCOMING AND INITIALIZATION
  print *, ''
  print *, 'Welcome to the low-rank approximation example!'
  print *, ''
  print *, 'We are going to use different methods to'
  print *, 'compute fast low-rank approximations and compare them.'
  print *, 'If you get no error here, you are fine!'
  print *, ''
  
  !Desired size
  n = 1000
  !Desired rank
  k = 10
  !Maximum number of steps for maxvol and maxvol-rect
  maxsteps = 2
  !Maximum number of row and column swaps for maxvol and maxvol-rect
  maxswaps = 4*k
  
  write(*,'(A,I0,A,I0,A)') ' We create ', n ,' by ', n, ' matrix'
  print *, 'with random singular vectors (now generating)...'
  call U%random(n)
  call V%random(n)
  call S%init(n)
  print *, 'Done!'
  print *, ''
  
  write(*,'(A,I0,A)') ' We set the first ', k, ' singular values'
  print *, 'to 100 and others to 1...'
  do i = 1, k
    S%d(i) = 100.0d0
  end do
  do i = k+1, n
    S%d(i) = 1.0d0
  end do
  A = U*(S .dot. V)
  print *, 'Done!'
  print *, ''
  
  write(*,'(A,I0,A)') ' We seek rank ', k, ' approximation, so Frobenius norm error of SVD is'
  print *, 'SVD error:', sqrt(dble(n-k))
  print *, ''
  
  !ACA
  print *, 'We start by trying adaptive cross approximation (ACA)'
  !We initialize the perm1 and perm2, which store permutations of rows and columns
  !We apply this permutations to put the maximum volume submatrix in the top left corner.
  !Some algorithms do it automatically.
  call per1%perm(n)
  call per2%perm(n)
  
  !Let's calculate the time
  time = dsecnd()
  
  !Find dominant k by k submatrix \hat A and construct C \hat A^{-1} R approximation
  !In case of ACA per1 is also used to select "test" rows. We use only 1
  call ACA(Aelem, A, n, n, k, 1, C, AR, per1, iNs_ = per2, jNs_ = per3)
  AR = .T.AR
  
  !Can also be used in combination with maxvol2 (uncomment to do it)
  !call per2%extend(n)
  !call per3%extend(n)
  !call C%permrows(per2, 1)
  !call AR%permcols(per3, 1)
  !C = .T.C
  !Ahat = C%subcols(k)
  !C = Ahat%rtsolve(C)
  !C = .T.C
  !Ahat = AR%subcols(k)
  !AR = Ahat%rtsolve(AR)
  !call maxvol2(Aelem, 2*k, 2*k, per2, per3, A, C, AR, .true.)
  
  time = dsecnd() - time
  print *, 'ACA time:', time
  
  E = A - C*AR
  print *, 'ACA error:', E%fnorm()
  print *, ''
  
  !MAXVOL
  print *, 'Next we perform MAXVOL approximation'
  call per1%perm(n)
  call per2%perm(n)
  
  !Let's calculate the time
  time = dsecnd()
  
  !Find dominant k by k submatrix \hat A and construct C \hat A^{-1} R approximation
  call maxvol(Aelem, n, n, k, per1, per2, A, C, AR, maxsteps, maxswaps, CA, .true.)
  
  time = dsecnd() - time
  print *, 'MAXVOL time:', time
  
  E = A - C*AR
  print *, 'MAXVOL error:', E%fnorm()
  print *, ''
  
  time = dsecnd()
  
  !MAXVOL2
  print *, 'Then let us try Householder-based MAXVOL2'
  print *, 'to construct fast CUR.'
  write(*,'(A,I0,A,I0)') ' We add rows and columns up to ', k, '*2 = ', k*2
  print *, 'Remember, that we search in rows and columns from MAXVOL'
  
  !Increase submatrix size to 2k by 2k
  call maxvol2(Aelem, 2*k, 2*k, per1, per2, A, CA, AR, .true.)
  
  time = dsecnd() - time
  print *, 'FAST CUR time:', time
  
  E = A - CA*AR
  print *, 'FAST CUR error:', E%fnorm()
  print *, ''
  
  !MAXVOL-PROJ
  print *, 'Nobody needs maxvol-rect separately, so let us use MAXVOL-PROJ'
  print *, 'We discard previous rows and columns'
  print *, '(To illustrate that maxvol-proj can work without initialization)'
  print *, 'And try to construct approximation from random start'
  call per1%perm(n)
  call per2%perm(n)
  
  time = dsecnd()
  
  !Large projective volume submatrix search
  !We set accuracy type to 0, since we expect relative error larger than 10^-9
  call maxvolproj(Aelem, n, n, k, 2*k, 2*k, per1, per2, A, C, AR, maxsteps, maxswaps, acc_type = 0)
  
  time = dsecnd() - time
  print *, 'MAXVOL-PROJ time:', time
  
  E = A - C*AR
  print *, 'MAXVOL-PROJ error:', E%fnorm()
  print *, ''
  
  !RRQR with DOMINANT_R
  print *, 'Now we construct Strong Rank Revealing QR with Dominant-R'
  !Let's reinitialize column permutations
  call per2%perm(n)
  !and rows
  call AR%init(k, n)
  
  time = dsecnd()
  
  !We'll need identity permuatation
  call peri%perm(n)
  !And some dummy permutation
  call per3%perm(n)
  
  !We will work in the copy of A
  call A1%copy(A)
  !We will use pre-maxvol to decrease time and increase accuracy
  !ABin and cin will be passed to Dominant-R to remove initialization of AB and c
  !Again, we use accuracy type equal to 0. Using dp one can check that error less than 10^-9 cant be reached with k columns
  call A1%premaxvol(k, per2, ABin, cin, dp = dp, acc_type = 0)
  !Let's limit ourselves to 2*k swaps
  call A1%dominantr(2, k, n, per2, swapsmade, 2*k, Ahat, AB, ABin, cin, dp = dp, acc_type = 0)
  
  !Uncomment to construct incomplete QR instead of CC^+A.
  !In case acc_type = 0, but dp turned out to be .true., it is better to use QR.
  !Ahat = .T.Ahat
  !We do not need Q; cin and Q are dummy variables here
  !call Ahat%halfLQ(R, cin, Q)
  !R = .T.R
  !C = Acols(Aelem, N, k, peri, per2, A)
  !Q = C*R
  do i = 1, k
    AR%d(i,i) = 1.0d0
  end do
  AR%d(:,k+1:n) = AB%d(:,1:n-k)
  !AR = R%rtsolve(AR)
  call AR%permcols(per2, 2)
  
  Ahat = A1%subcols(k)
  
  time = dsecnd() - time
  print *, 'DOMINANT-R RRQR time:', time
  
  E = A - Ahat*AR
  !Uncomment to construct incomplete QR instead of CC^+A.
  !E = A - Q*AR
  print *, 'DOMINANT-R RRQR error:', E%fnorm()
  
  !Quasioptimal approximation
  print *, 'Finally, we construct close to optimal column approximation'
  print *, 'First, let us use SVD vectors (given in advance)'
  print *, 'We will now generate A with exponentially decreasing singular values 2^-k'
  do i = 1, n
    S%d(i) = 2.0d0**(-i)
  end do
  A = U*(S .dot. V)
  print *, 'SVD error:', sqrt(sum(S%d(k+1:)**2))
  print *, ''
  !Since SVD is known, we use the version, which is given A - AV^TV in order to avoid SVD calculation
  E = S%subarray(n, k+1) .dot. V%subrows(n, k+1)
  !Select first k right singular vectors
  V = V%subrows(k)
  !Let's reinitialize column permutations
  call per2%perm(n)
  
  time = dsecnd()
  
  !To use with arbitrary V, calculate E as A - A*V^T*V directly
  !This only selects the columns, which are written in per2
  call BestColsOrth(V, E, per2)
  
  time = dsecnd() - time
  
  print *, 'Column selection time:', time
  print *, ''
  print *, 'Now compute the error of approximation with the selected columns'
  C = Acols(Aelem,n,k,peri,per2,A)
  call C%qr(Q,R)
  E = A - Q*(Q .Td. A)
  print *, 'Column approximation error:', E%fnorm()
  print *, ''
  print *, 'Now let us select rows for C \hat A^-1 R approximation'
  
  Q = .T.Q
  E = .T.E
  call per1%perm(n)
  
  time = dsecnd()
  
  !This only selects the rows, which are written in per1
  !Q and E are reused from column approximation
  call BestColsOrth(Q, E, per1)
  
  time = dsecnd() - time
  
  print *, 'Row selection time:', time
  print *, ''
  print *, 'Now compute the error of the skeleton approximation'
  AR = Arows(Aelem,k,n,per1,peri,A)
  !Select a submatrix
  Ahat = Arows(Aelem,k,k,per1,per2,A)
  !This can be done a little faster, if the resulting Q from BestColsOrth is saved
  AR = Ahat .Id. AR
  E = A - C*AR
  print *, 'Skeleton approximation error:', E%fnorm()
  print *, ''
  
  print *, 'Now let us show how to reach almost the same in O(Nr^2) without SVD'
  print *, 'I will use ACA as an initial approximation: it can be replaced with projective volume for better accuracy'
  
  call per1%perm(n)
  call per2%perm(n)
  call per3%perm(n)
  
  !Let's calculate the time
  time = dsecnd()
  
  call ACA(Aelem, A, n, n, k+10, 1, C, AR, per1, iNs_ = per2, jNs_ = per3)
  AR = .T.AR
  
  call AR%lq(R,Q)
  Ahat = .T.(C*R)
  call per1%perm(n)
  !Typical error of ACA is usually no more than 10 times larger than SVD error (Frobenius norm should be used instead of spectral)
  !Here I substitute the bound, in practice it should be estimated (see inside ACA for ways to do that)
  !Worst-case scenario is spectral norm of error * sqrt(n-k-5), but in practice Frobenius norm should be enough
  !However, for some reason, zero regularization works even better; I will leave the equation however
  call BestCols(Ahat, per1, k, 0*10*2.0d0**(-(k+11)))
  CA = Arows(Aelem, k, N, per1, peri, A)
  
  time = dsecnd() - time
  print *, 'ACA-based rows time:', time
  E = A - (C * (AR .dI. CA)) * CA
  print *, 'ACA-based row error:', E%fnorm()
  
  time = dsecnd()
  
  call CA%lq(R,V)
  AR = AR - ((AR .dT. V) * V)
  call C%qr(Q,R)
  Ahat = R*AR
  call per2%perm(n)
  call BestColsOrth(V, Ahat, per2, reg=0*10*2.0d0**(-(k+11)))
  Ahat = Acols(Aelem, k, k, peri, per2, CA)
  C = Acols(Aelem, N, k, peri, per2, A)
  time = dsecnd() - time
  print *, 'ACA-based CUR time:', time
  E = A - C*(Ahat .Id. CA)
  print *, 'ACA-based CUR error:', E%fnorm()
end
