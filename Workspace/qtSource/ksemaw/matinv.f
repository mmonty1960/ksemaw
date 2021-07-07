C subroutine per il calcolo della matrice inversa con le subroutine
C fornite da Numerical Recipe
C  n = dimension of the matrix to be inverted
C  np = max dimension of matrices used in called subroutine
C  as(n,n) = matrix to be inverted
C  b(n,n) = inverted matrix
      SUBROUTINE MATINV(n,np,as,b)
      INTEGER n,indx(np)
      REAL as(np,np),a(np,np),y(np,np),b(np,np)

***** copy as-matrix into a-matrix that will be lost
      do i=1,n
       do j=1,n
         a(i,j)=as(i,j)
       end do
      end do

*****Set up identity matrix.     
      do i=1,n
       do j=1,n
        y(i,j)=0.
       end do
       y(i,i)=1.
      end do
      call ludcmp(a,n,np,indx,d) ! Decompose the matrix just once.
      do j=1,n ! Find inverse by columns.
        call lubksb(a,n,np,indx,y(1,j))
      end do
C Note that FORTRAN stores two-dimensional matrices by column,
C so y(1,j) is the address of the jth column of y.

***** save results into b-matrix
      do i=1,n
       do j=1,n
         b(i,j)=y(i,j)
       end do
      end do
      
      return
      end



      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(np)
      REAL a(np,np),b(np)
*** Solves the set of n linear equations A * X = B. Here a is input,
*** not as the matrix A but rather as its LU decomposition, determined
*** by the routine ludcmp. indx is input as the permutation vector
*** returned by ludcmp. b(1:n) is input as the right-hand side vector B,
*** and returns with the solution vector X. a, n, np, and indx are not
*** modified by this routine and can be left in place for successive
*** calls with different right-hand sides b. This routine takes into 
*** account the possibility that b will begin with many zero elements,
*** so it is e cient for use in matrix inversion.
      INTEGER i,ii,j,ll
      REAL sum
      ii=0 
c When ii is set to a positive value, it will become the index
C of the first nonvanishing element of b. We now do
C the forward substitution, equation (2.3.6). The only new
c wrinkle is to unscramble the permutation as we go.
      do i=1,n
       ll=indx(i)
       sum=b(ll)
       b(ll)=b(i)
       if (ii.ne.0) then
        do j=ii,i-1
         sum=sum-a(i,j)*b(j)
        end do
       else if (sum.ne.0.) then
        ii=i ! A nonzero element was encountered, so from now on we
c              will have to do the sums in the loop above.
       endif
       b(i)=sum
      end do
      do i=n,1,-1 ! Now we do the backsubstitution, equation (2.3.7).
       sum=b(i)
       do j=i+1,n
        sum=sum-a(i,j)*b(j)
       end do
       b(i)=sum/a(i,i) ! Store a component of the solution vector X.
      end do
      return !All done!
      END



      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(np),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
*** From "Numerical recipe"
*** Largest expected n, and a small number.
*** Given a matrix a(1:n,1:n), with physical dimension np by np,
*** this routine replaces it by the LU decomposition of a rowwise
*** permutation of itself. a and n are input.a is output, arranged
*** as in equation (2.3.14) above; indx(1:n) is an output vector
*** that records the row permutation effected by the partial pivoting;
*** d is output as +-1 depending on whether the number of row interchanges
*** was even or odd,respectively. This routine is used in combination with
*** lubksb to solve linear equations or invert a matrix.
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)   ! vv stores the implicit scaling of
c                                      each row.
      d=1.                          ! No row interchanges yet.
      do i=1,n ! Loop over rows to get the implicit scaling information.
       aamax=0.
       do j=1,n
        if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
       end do
       if (aamax.eq.0.) then
c        pause
        write(*,*) 'singular matrix in ludcmp!!!!'
c                   ! No nonzero largest element
       end if
       vv(i)=1./aamax    !Save the scaling.
      end do
      do j=1,n      ! This is the loop over columns of Crout's method.
       do i=1,j-1   ! This is equation (2.3.12) except for i=j.
        sum=a(i,j)
        do k=1,i-1
         sum=sum-a(i,k)*a(k,j)
        end do
        a(i,j)=sum
       end do
       aamax=0.   ! Initialize for the search for largest pivot element.
       do i=j,n   ! This is i=j of equation (2.3.12) and i=j+1...N of 
c                   equation (2.3.13).
       sum=a(i,j)
       do k=1,j-1
        sum=sum-a(i,k)*a(k,j)
       end do
       a(i,j)=sum
       dum=vv(i)*abs(sum)      ! Figure of merit for the pivot.
       if (dum.ge.aamax) then  ! Is it better than the best so far?
        imax=i
        aamax=dum
       end if
       end do
       if (j.ne.imax) then      ! Do we need to interchange rows?
        do k=1,n                ! Yes, do so...
         dum=a(imax,k)
         a(imax,k)=a(j,k)
         a(j,k)=dum
        end do
        d=-d                   ! ...and change the parity of d
        vv(imax)=vv(j)         ! Also interchange the scale factor.
       end if
       indx(j)=imax
       if(a(j,j).eq.0.) a(j,j)=TINY
c If the pivot element is zero the matrix is singular (at least to the
c precision of the algorithm).
c For some applications on singular matrices, it is desirable to substitute TINY
c for zero.
       if(j.ne.n) then        ! Now,finally, divide by the pivot element.
        dum=1./a(j,j)
        do i=j+1,n
          a(i,j)=a(i,j)*dum
        end do
       end if
      end do         ! Go back for the next column in the reduction.
      return
      END
