      subroutine errjac(n,x,fjac,ldfjac,nprob)
      integer n,ldfjac,nprob
      real x(n),fjac(ldfjac,n)
c     **********
c
c     subroutine errjac
c
c     this subroutine is derived from vecjac which defines the
c     jacobian matrices of fourteen test functions. the problem
c     dimensions are as described in the prologue comments of vecfcn.
c     various errors are deliberately introduced to provide a test
c     for chkder.
c
c     the subroutine statement is
c
c       subroutine errjac(n,x,fjac,ldfjac,nprob)
c
c     where
c
c       n is a positive integer variable.
c
c       x is an array of length n.
c
c       fjac is an n by n array. on output fjac contains the
c         jacobian matrix, with various errors deliberately
c         introduced, of the nprob function evaluated at x.
c
c       ldfjac is a positive integer variable not less than n
c         which specifies the leading dimension of the array fjac.
c
c       nprob is a positive integer variable which defines the
c         number of the problem. nprob must not exceed 14.
c
c     subprograms called
c
c       fortran-supplied ... atan,cos,exp,amin1,sin,sqrt,
c                            max0,min0
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,ivar,j,k,k1,k2,ml,mu
      real c1,c3,c4,c5,c6,c9,eight,fiftn,five,four,h,hundrd,one,prod,
     *     six,sum,sum1,sum2,temp,temp1,temp2,temp3,temp4,ten,three,
     *     ti,tj,tk,tpi,twenty,two,zero
      real float
      data zero,one,two,three,four,five,six,eight,ten,fiftn,twenty,
     *     hundrd
     *     /0.0e0,1.0e0,2.0e0,3.0e0,4.0e0,5.0e0,6.0e0,8.0e0,1.0e1,
     *      1.5e1,2.0e1,1.0e2/
      data c1,c3,c4,c5,c6,c9 /1.0e4,2.0e2,2.02e1,1.98e1,1.8e2,2.9e1/
      float(ivar) = ivar
c
c     jacobian routine selector.
c
      go to (10,20,50,60,90,100,200,230,290,320,350,380,420,450),
     *      nprob
c
c     rosenbrock function with sign reversal affecting element (1,1).
c
   10 continue
      fjac(1,1) = one
      fjac(1,2) = zero
      fjac(2,1) = -twenty*x(1)
      fjac(2,2) = ten
      go to 490
c
c     powell singular function with sign reversal affecting element
c     (3,3).
c
   20 continue
      do 40 k = 1, 4
         do 30 j = 1, 4
            fjac(k,j) = zero
   30       continue
   40    continue
      fjac(1,1) = one
      fjac(1,2) = ten
      fjac(2,3) = sqrt(five)
      fjac(2,4) = -fjac(2,3)
      fjac(3,2) = two*(x(2) - two*x(3))
      fjac(3,3) = two*fjac(3,2)
      fjac(4,1) = two*sqrt(ten)*(x(1) - x(4))
      fjac(4,4) = -fjac(4,1)
      go to 490
c
c     powell badly scaled function with the sign of the jacobian
c     reversed.
c
   50 continue
      fjac(1,1) = -c1*x(2)
      fjac(1,2) = -c1*x(1)
      fjac(2,1) = exp(-x(1))
      fjac(2,2) = exp(-x(2))
      go to 490
c
c     wood function without error.
c
   60 continue
      do 80 k = 1, 4
         do 70 j = 1, 4
            fjac(k,j) = zero
   70       continue
   80    continue
      temp1 = x(2) - three*x(1)**2
      temp2 = x(4) - three*x(3)**2
      fjac(1,1) = -c3*temp1 + one
      fjac(1,2) = -c3*x(1)
      fjac(2,1) = -two*c3*x(1)
      fjac(2,2) = c3 + c4
      fjac(2,4) = c5
      fjac(3,3) = -c6*temp2 + one
      fjac(3,4) = -c6*x(3)
      fjac(4,2) = c5
      fjac(4,3) = -two*c6*x(3)
      fjac(4,4) = c6 + c4
      go to 490
c
c     helical valley function with multiplicative error affecting
c     elements (2,1) and (2,2).
c
   90 continue
      tpi = eight*atan(one)
      temp = x(1)**2 + x(2)**2
      temp1 = tpi*temp
      temp2 = sqrt(temp)
      fjac(1,1) = hundrd*x(2)/temp1
      fjac(1,2) = -hundrd*x(1)/temp1
      fjac(1,3) = ten
      fjac(2,1) = five*x(1)/temp2
      fjac(2,2) = five*x(2)/temp2
      fjac(2,3) = zero
      fjac(3,1) = zero
      fjac(3,2) = zero
      fjac(3,3) = one
      go to 490
c
c     watson function with sign reversals affecting the computation of
c     temp1.
c
  100 continue
      do 120 k = 1, n
         do 110 j = k, n
            fjac(k,j) = zero
  110       continue
  120    continue
      do 170 i = 1, 29
         ti = float(i)/c9
         sum1 = zero
         temp = one
         do 130 j = 2, n
            sum1 = sum1 + float(j-1)*temp*x(j)
            temp = ti*temp
  130       continue
         sum2 = zero
         temp = one
         do 140 j = 1, n
            sum2 = sum2 + temp*x(j)
            temp = ti*temp
  140       continue
         temp1 = two*(sum1 + sum2**2 + one)
         temp2 = two*sum2
         temp = ti**2
         tk = one
         do 160 k = 1, n
            tj = tk
            do 150 j = k, n
               fjac(k,j) = fjac(k,j)
     *                     + tj
     *                       *((float(k-1)/ti - temp2)
     *                         *(float(j-1)/ti - temp2) - temp1)
               tj = ti*tj
  150          continue
            tk = temp*tk
  160       continue
  170    continue
      fjac(1,1) = fjac(1,1) + six*x(1)**2 - two*x(2) + three
      fjac(1,2) = fjac(1,2) - two*x(1)
      fjac(2,2) = fjac(2,2) + one
      do 190 k = 1, n
         do 180 j = k, n
            fjac(j,k) = fjac(k,j)
  180       continue
  190    continue
      go to 490
c
c     chebyquad function with jacobian twice correct size.
c
  200 continue
      tk = one/float(n)
      do 220 j = 1, n
         temp1 = one
         temp2 = two*x(j) - one
         temp = two*temp2
         temp3 = zero
         temp4 = two
         do 210 k = 1, n
            fjac(k,j) = two*tk*temp4
            ti = four*temp2 + temp*temp4 - temp3
            temp3 = temp4
            temp4 = ti
            ti = temp*temp2 - temp1
            temp1 = temp2
            temp2 = ti
  210       continue
  220    continue
      go to 490
c
c     brown almost-linear function without error.
c
  230 continue
      prod = one
      do 250 j = 1, n
         prod = x(j)*prod
         do 240 k = 1, n
            fjac(k,j) = one
  240       continue
         fjac(j,j) = two
  250    continue
      do 280 j = 1, n
         temp = x(j)
         if (temp .ne. zero) go to 270
         temp = one
         prod = one
         do 260 k = 1, n
            if (k .ne. j) prod = x(k)*prod
  260       continue
  270    continue
         fjac(n,j) = prod/temp
  280    continue
      go to 490
c
c     discrete boundary value function with multiplicative error
c     affecting the jacobian diagonal.
c
  290 continue
      h = one/float(n+1)
      do 310 k = 1, n
         temp = three*(x(k) + float(k)*h + one)**2
         do 300 j = 1, n
            fjac(k,j) = zero
  300       continue
         fjac(k,k) = four + temp*h**2
         if (k .ne. 1) fjac(k,k-1) = -one
         if (k .ne. n) fjac(k,k+1) = -one
  310    continue
      go to 490
c
c     discrete integral equation function with sign error affecting
c     the jacobian diagonal.
c
  320 continue
      h = one/float(n+1)
      do 340 k = 1, n
         tk = float(k)*h
         do 330 j = 1, n
            tj = float(j)*h
            temp = three*(x(j) + tj + one)**2
            fjac(k,j) = h*amin1(tj*(one-tk),tk*(one-tj))*temp/two
  330       continue
         fjac(k,k) = fjac(k,k) - one
  340    continue
      go to 490
c
c     trigonometric function with sign errors affecting the
c     offdiagonal elements of the jacobian.
c
  350 continue
      do 370 j = 1, n
         temp = sin(x(j))
         do 360 k = 1, n
            fjac(k,j) = -temp
  360       continue
         fjac(j,j) = float(j+1)*temp - cos(x(j))
  370    continue
      go to 490
c
c     variably dimensioned function with operation error affecting
c     the upper triangular elements of the jacobian.
c
  380 continue
      sum = zero
      do 390 j = 1, n
         sum = sum + float(j)*(x(j) - one)
  390    continue
      temp = one + six*sum**2
      do 410 k = 1, n
         do 400 j = k, n
            fjac(k,j) = float(k*j)/temp
            fjac(j,k) = fjac(k,j)
  400       continue
         fjac(k,k) = fjac(k,k) + one
  410    continue
      go to 490
c
c     broyden tridiagonal function without error.
c
  420 continue
      do 440 k = 1, n
         do 430 j = 1, n
            fjac(k,j) = zero
  430       continue
         fjac(k,k) = three - four*x(k)
         if (k .ne. 1) fjac(k,k-1) = -one
         if (k .ne. n) fjac(k,k+1) = -two
  440    continue
      go to 490
c
c     broyden banded function with sign error affecting the jacobian
c     diagonal.
c
  450 continue
      ml = 5
      mu = 1
      do 480 k = 1, n
         do 460 j = 1, n
            fjac(k,j) = zero
  460       continue
         k1 = max0(1,k-ml)
         k2 = min0(k+mu,n)
         do 470 j = k1, k2
            if (j .ne. k) fjac(k,j) = -(one + two*x(j))
  470       continue
         fjac(k,k) = two - fiftn*x(k)**2
  480    continue
  490 continue
      return
c
c     last card of subroutine errjac.
c
      end

