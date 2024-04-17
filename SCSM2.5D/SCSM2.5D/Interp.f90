!C----------------------------------------------------------------------C
!C                                                                      C
!C     This subroutine calculates the coeffiecents b(i)=a(1,i), c(i)    C
!C     =A(2,i) and c(i)=A(3,i) of the cubic-spline interpolations:       C
!C       f(x)= y(i)+b(i)*(x-x(i))+c(i)*(x-x(i))**2+d(i)*(x-x(i))**3,    C
!C     in the intervals [x(i),x(i+1)], (i=1,2,...,n-1).                 C
!C                                                                      C
!C     Entries:                                                         C
!C               (1) x(n),f(n)............sampling points of x & f(x);  C
!C                                                                      C
!C     Returns:                                                         C
!C               (1) a(3,n).....three coefficients b(i),c(i) and d(i);  C
!C                                                                      C
!C----------------------------------------------------------------------C
      SUBROUTINE SPLINE_1D_A(N,X,F,A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),F(N),A(3,N)
      
      N1=N-1
      IF(N<2)THEN  ! CHECK INPUT
      WRITE(*,*)'   POINTS ARE NOT ENOUGH (N<2) !'
      GO TO 10
      ENDIF
      
      IF(N<3)THEN ! LINEAR INTERPOLATION
      A(1,1)=(F(2)-F(1))/(X(2)-X(1))
      A(2,1)=0.D0
      A(3,1)=0.D0
      A(1,2)=A(1,1)
      A(2,2)=0.D0
      A(3,2)=0.D0
      GO TO 10
      ENDIF

      ! STEP 1: PREPARATION

      A(3,1)=X(2)-X(1)
      A(2,2)=(F(2)-F(1))/A(3,1)
      DO I=2,N1
      A(3,I)=X(I+1)-X(I)
      A(1,I)=2.D0*(A(3,I-1)+A(3,I))
      A(2,I+1)=(F(I+1)-F(I))/A(3,I)
      A(2,I)=A(2,I+1)-A(2,I)
      END DO
      
      ! STEP 2: END CONDITIONS 

      A(1,1)=-A(3,1)
      A(1,N)=-A(3,N-1)
      A(2,1)=0.D0
      A(2,N)=0.D0
      IF(N/=3) THEN
      A(2,1)=A(2,3)/(X(4)-X(2))-A(2,2)/(X(3)-X(1))
      A(2,N)=A(2,N-1)/(X(N)-X(N-2))-A(2,N-2)/(X(N-1)-X(N-3))
      A(2,1)=A(2,1)*A(3,1)**2/(X(4)-X(1))
      A(2,N)=-A(2,N)*A(3,N-1)**2/(X(N)-X(N-3))
      END IF
      
      ! STEP 3: FORWARD ELIMINATION 

      DO I=2,N
      H=A(3,I-1)/A(1,I-1)
      A(1,I)=A(1,I)-H*A(3,I-1)
      A(2,I)=A(2,I)-H*A(2,I-1)
      END DO

      ! STEP 4: BACK SUBSTITUTION

      A(2,N)=A(2,N)/A(1,N)
      DO J=1,N1
      I=N-J
      A(2,I)=(A(2,I)-A(3,I)*A(2,I+1))/A(1,I)
      END DO

      ! STEP 5: COMPUTE SPLINE COEFFICIENTS: A(1,N),A(2,N),A(3,N)

      A(1,N)=(F(N)-F(N1))/A(3,N1)+A(3,N1)*(A(2,N1)+2.D0*A(2,N))
      DO I = 1, N1
      A(1,I)=(F(I+1)-F(I))/A(3,I)-A(3,I)*(A(2,I+1)+2.D0*A(2,I))
      A(3,I)=(A(2,I+1)-A(2,I))/A(3,I)
      A(2,I)=3.D0*A(2,I)
      END DO
      A(2,N)=3.D0*A(2,N)
      A(3,N)=A(3,N-1)
   
   10 RETURN
      END
      
!C----------------------------------------------------------------------C
!C                                                                      C
!C     This subroutine calculates the interpolated value fi,1st & 2nd   C
!C     derivatives (fx,fxx) of cubic-spline interpolation function:     C
!C     f(x)= y(i)+b(i)*(x-x(i))+c(i)*(x-x(i))**2+d(i)*(x-x(i))**3, in   C
!C     the intervals [x(i),x(i+1)], (i=1,2,...,n-1).                    C
!C                                                                      C
!C     Entries:                                                         C
!C               (1) xi...........                    estimated point;  C
!C               (2) x(n),y(n)............sampling points of x & f(x);  C
!c               (3) a(3,n)................cofficients:b(n),c(n),d(n);  C
!C                                                                      C
!C     Returns:                                                         C
!C               (1) fi.................function value at xi:fi=f(xi);  C
!C               (2) fx.................the 1st derivative: fx=f'(xi);  C
!C               (3) fxx...............the 2nd derivative:fxx=f''(xi);  C
!C                                                                      C
!----------------------------------------------------------------------C
      SUBROUTINE SPLINE_1D_F(XI,N,X,F,A,FI,FX,FXX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),F(N),A(3,N)
      
      IF(XI <= X(1)) THEN
      FI = F(1)
      FX = A(1,1)
      FXX= 0.D0
      GO TO 10
      END IF

      IF(XI >= X(N)) THEN
      FI = F(N)
      FX = A(1,N)
      FXX= 0.D0
      GO TO 10
      END IF

      I = 1
      J = N+1
      DO WHILE (J > I+1)
      K = (I+J)/2
      IF(XI < X(K)) THEN
      J=K
      ELSE
      I=K
      END IF
      END DO

      DX = XI - X(I)
      FI = F(I) + DX*(A(1,I) + DX*(A(2,I) + DX*A(3,I)))
      FX =A(1,I)+DX*(2.D0*A(2,I)+3.D0*A(3,I)*DX)
      FXX=2.D0*A(2,I)+6.D0*A(3,I)*DX
      
  10  RETURN
      END
      
      
!C----------------------------------------------------------------------C
!C                                                                      C
!C     calculate the Lagrange functions, the 1st & 2nd derivatives at   C
!C     givem point Xi.  .                                               C
!C                                                                      C
!C     Entries:                                                         C
!C                                                                      C
!C            X(N)..............................N points for Lagrange;  C
!C            Xi...................................interpolated point;  C
!C            IORD......=0, 1 or 2 for Lagrange, 1st & 2nd derivative;  C
!C                                                                      C
!C     Returns:                                                         C
!C                                                                      C
!C            DLI(N)..........Lagrange, 1st or 2nd derivative vectors;  C
!C                                                                      C
!C----------------------------------------------------------------------C
      SUBROUTINE CDLI(XI,N,X,IORD,DLI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),DLI(N),ID(:),IQ(:,:),IP(:)
      ALLOCATABLE ID,IQ,IP

!---- for li(x)/li(xi) ------------
      IF(IORD.EQ.0)THEN          
      DO I=1,N
      F1=1.D0
      F2=1.D0
      DO 2 J=1,N
      IF(J.EQ.I)GO TO 2
      F1=F1*(  XI-X(J))
      F2=F2*(X(I)-X(J))
   2  CONTINUE
      DLI(I)=F1/F2
      ENDDO
      ENDIF
      
!---- for 1st & 2nd derivatives ---
      IF(IORD.NE.0)THEN
      ALLOCATE(ID(N-1),IQ(N-1,N-1),IP(N-1))
                
      DO I=1,N
      !for li(xi)
      F1=1.D0
      II=0
      DO 3 J=1,N
      IF(J.EQ.I)GO TO 3
      F1=F1*(X(I)-X(J))
      II=II+1
      ID(II)=J
   3  CONTINUE
      
      !for li'(x)
      S=0.D0
      DO 6 K=1,II
      I1=ID(K)      
      F2=1.D0
      JJ=0 
      DO 4 J=1,II
      J1=ID(J)
      IF(J1.EQ.I1)GO TO 4
      F2=F2*(XI-X(J1))
      JJ=JJ+1
      IQ(K,JJ)=J1
   4  CONTINUE
      IP(K)=JJ   
      S=S+F2
   6  CONTINUE
      
      IF(IORD.EQ.1)GO TO 14    
      !for li''(x)
      S=0.D0
      DO 12 K=1,II
      DO 10 J=1,IP(K)
      I1=IQ(K,J)
      F2=1.D0
      DO 8 L=1,IP(K)
      J1=IQ(K,L)
      IF(J1.EQ.I1)GO TO 8
      F2=F2*(XI-X(J1))
   8  CONTINUE
      S=S+F2
  10  CONTINUE
  12  CONTINUE        
  14  DLI(I)=S/F1
      ENDDO
      DEALLOCATE (ID,IQ,IP)
      ENDIF
             
      RETURN
      END