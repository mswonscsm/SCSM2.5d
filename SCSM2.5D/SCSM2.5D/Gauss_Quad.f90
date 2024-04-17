Subroutine FFT_SR(Nky,ik,dkyc,wky,RS,R,Weight)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
Complex*16 :: Weight

Dimension tot_ky(NKY),y(NKY),SR0(NKY)
Complex*16,Dimension(:) :: SR1(NKY)
PI=4.d0*atan(1.d0)

do i=1,Nky
    Tot_ky(i)=dble(i-1)/dble(nky-1)*dkyc
enddo
dy=2*pi/(Tot_ky(2)*dble(nky))

do i=1,nky
    y(i)=dble(i-1)*dy
enddo

do i=1,nky
    RY=sqrt(r**2+y(i)**2)
    !SR0(i)=3.d0*pi/(4.d0*(pi**2-6)*RS**3)*(1.d0+cos(pi*RY/RS))
    SR0(i)=3.d0*pi/(4.d0*(pi**2-6)*y(nky)**3)*(1.d0+cos(pi*RY/y(nky)))
enddo

SR1(:)=0.d0
do i=1,nky
    do k=1,nky
        SR1(i)=SR1(i)+dcmplx(SR0(k)*dcos(Tot_ky(i)*y(k)),-SR0(k)*dsin(Tot_ky(i)*y(k)))
    enddo
enddo
Weight=SR1(ik)
!write(91,*)real(SR1(1:NKY))
!write(92,*)imag(SR1(1:NKY))
end subroutine

!-------------------------------------------------------C
!                                                       C
!     This suboutine is to calculate the points Z(*)    C
!     weights W(*) of Gauss-legendre-Lobbato integ-     C
!     ration formula.                                   C
!                                                       C
!     Entry:                                            C
!            NP.....................number of points    C
!                 (NP-1 is the degree of polynormial)   C
!                                                       C
!     Returns:                                          C
!            Z(*)..........points in the range [-1,1],  C
!            W(*).............weights for integration.  C
!                                                       C
!                                                       C
!     this subroutine involves the following routines:  C
!                                                       C
!            ZWGJD(), JACOBF(), ENDW1(), ENDW2(),       C
!            GAMMAF(), JACG(), PNORMJ().                C
!                                                       C
!-------------------------------------------------------C
      SUBROUTINE GLL(NP,Z,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(*),W(*)
      ALPHA = 0.0
      BETA  = 0.0
      N     = NP-1
      NM1   = N-1
      ONE   = 1.D0
      TWO   = 2.D0

      IF (NP.LE.1) THEN
      WRITE (6,*) 'Minimum number of Gauss-Lobatto points is 2'
      STOP
      ENDIF
      IF ((ALPHA.LE.-ONE).OR.(BETA.LE.-ONE)) THEN
      WRITE (6,*) 'Alpha and Beta must be greater than -1'
      STOP
      ENDIF

      IF (NM1.GT.0) THEN
      ALPG  = ALPHA+ONE
      BETG  = BETA+ONE
      CALL ZWGJD (Z(2),W(2),NM1,ALPG,BETG)
      ENDIF
      Z(1)  = -ONE
      Z(NP) =  ONE
      DO  110 I=2,NP-1
      W(I) = W(I)/(ONE-Z(I)**2)
  110 CONTINUE
      CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z(1))
      W(1)  = ENDW1 (N,ALPHA,BETA)/(TWO*PD)
      CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z(NP))
      W(NP) = ENDW2 (N,ALPHA,BETA)/(TWO*PD)

      RETURN
      END

      SUBROUTINE ZWGJD (Z,W,NP,ALPHA,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(*),W(*)

      N     = NP-1
      DN    = DBLE(FLOAT(N))
      ONE   = 1.D0
      TWO   = 2.D0
      APB   = ALPHA+BETA

      IF (NP.LE.0) THEN
      WRITE (6,*) 'Minimum number of Gauss points is 1'
      STOP
      ENDIF
      IF ((ALPHA.LE.-ONE).OR.(BETA.LE.-ONE)) THEN
      WRITE (6,*) 'Alpha and Beta must be greater than -1'
      STOP
      ENDIF

      IF (NP.EQ.1) THEN
      Z(1) = (BETA-ALPHA)/(APB+TWO)
      W(1) = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE)/GAMMAF(APB+TWO)&
                * TWO**(APB+ONE)
      RETURN
      ENDIF

      CALL JACG (Z,NP,ALPHA,BETA)

      NP1   = N+1
      NP2   = N+2
      DNP1  = DBLE(FLOAT(NP1))
      DNP2  = DBLE(FLOAT(NP2))
      FAC1  = DNP1+ALPHA+BETA+ONE
      FAC2  = FAC1+DNP1
      FAC3  = FAC2+ONE
      FNORM = PNORMJ(NP1,ALPHA,BETA)
      RCOEF = (FNORM*FAC2*FAC3)/(TWO*FAC1*DNP2)
      DO 100 I=1,NP
      CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,NP2,ALPHA,BETA,Z(I))
      W(I) = -RCOEF/(P*PDM1)
 100  CONTINUE
      RETURN
      END

      SUBROUTINE JACOBF (POLY,PDER,POLYM1,PDERM1,POLYM2,PDERM2,&
                         N,ALP,BET,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      APB  = ALP+BET
      POLY = 1.D0
      PDER = 0.D0
      IF (N .EQ. 0) RETURN
      POLYL = POLY
      PDERL = PDER
      POLY  = (ALP-BET+(APB+2.D0)*X)/2.D0
      PDER  = (APB+2.D0)/2.D0
      IF (N .EQ. 1) RETURN
      DO 20 K=2,N
      DK = DBLE(FLOAT(K))
      A1 = 2.D0*DK*(DK+APB)*(2.D0*DK+APB-2.D0)
      A2 = (2.D0*DK+APB-1.D0)*(ALP**2-BET**2)
      B3 = (2.D0*DK+APB-2.D0)
      A3 = B3*(B3+1.D0)*(B3+2.D0)
      A4 = 2.D0*(DK+ALP-1.D0)*(DK+BET-1.D0)*(2.D0*DK+APB)
      POLYN  = ((A2+A3*X)*POLY-A4*POLYL)/A1
      PDERN  = ((A2+A3*X)*PDER-A4*PDERL+A3*POLY)/A1
      PSAVE  = POLYL
      PDSAVE = PDERL
      POLYL  = POLY
      POLY   = POLYN
      PDERL  = PDER
      PDER   = PDERN
  20  CONTINUE
      POLYM1 = POLYL
      PDERM1 = PDERL
      POLYM2 = PSAVE
      PDERM2 = PDSAVE
      RETURN
      END

      DOUBLE PRECISION FUNCTION ENDW1 (N,ALPHA,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ZERO  = 0.D0
      ONE   = 1.D0
      TWO   = 2.D0
      THREE = 3.D0
      FOUR  = 4.D0
      APB   = ALPHA+BETA
      IF (N.EQ.0) THEN
      ENDW1 = ZERO
      RETURN
      ENDIF
      F1   = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE)
      F1   = F1*(APB+TWO)*TWO**(APB+TWO)/TWO
      IF (N.EQ.1) THEN
      ENDW1 = F1
      RETURN
      ENDIF
      FINT1 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE)
      FINT1 = FINT1*TWO**(APB+TWO)
      FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR)
      FINT2 = FINT2*TWO**(APB+THREE)
      F2    = (-TWO*(BETA+TWO)*FINT1 + (APB+FOUR)*FINT2)&
              * (APB+THREE)/FOUR
      IF (N.EQ.2) THEN
      ENDW1 = F2
      RETURN
      ENDIF
      DO 100 I=3,N
      DI   = DBLE(FLOAT(I-1))
      ABN  = ALPHA+BETA+DI
      ABNN = ABN+DI
      A1   = -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE))
      A2   =  (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO))
      A3   =  (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE))
      F3   =  -(A2*F2+A1*F1)/A3
      F1   = F2
      F2   = F3
 100  CONTINUE
      ENDW1  = F3
      RETURN
      END

      DOUBLE PRECISION FUNCTION ENDW2 (N,ALPHA,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ZERO  = 0.D0
      ONE   = 1.D0
      TWO   = 2.D0
      THREE = 3.D0
      FOUR  = 4.D0
      APB   = ALPHA+BETA
      IF (N.EQ.0) THEN
      ENDW2 = ZERO
      RETURN
      ENDIF
      F1   = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE)
      F1   = F1*(APB+TWO)*TWO**(APB+TWO)/TWO
      IF (N.EQ.1) THEN
      ENDW2 = F1
      RETURN
      ENDIF
      FINT1 = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE)
      FINT1 = FINT1*TWO**(APB+TWO)
      FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR)
      FINT2 = FINT2*TWO**(APB+THREE)
      F2    = (TWO*(ALPHA+TWO)*FINT1 - (APB+FOUR)*FINT2)&
              * (APB+THREE)/FOUR
      IF (N.EQ.2) THEN
      ENDW2 = F2
      RETURN
      ENDIF
      DO 100 I=3,N
      DI   = DBLE(FLOAT(I-1))
      ABN  = ALPHA+BETA+DI
      ABNN = ABN+DI
      A1   =  -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE))
      A2   =  (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO))
      A3   =  (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE))
      F3   =  -(A2*F2+A1*F1)/A3
      F1   = F2
      F2   = F3
 100  CONTINUE
      ENDW2  = F3
      RETURN
      END

      DOUBLE PRECISION FUNCTION GAMMAF (X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ZERO = 0.0D0
      HALF = 0.5D0
      ONE  = 1.0D0
      TWO  = 2.0D0
      FOUR = 4.0D0
      PI   = FOUR*DATAN(ONE)
      GAMMAF = ONE
      IF (X.EQ.-HALF) GAMMAF = -TWO*DSQRT(PI)
      IF (X.EQ. HALF) GAMMAF =  DSQRT(PI)
      IF (X.EQ. ONE ) GAMMAF =  ONE
      IF (X.EQ. TWO ) GAMMAF =  ONE
      IF (X.EQ. 1.5D0) GAMMAF =  DSQRT(PI)/2.D0
      IF (X.EQ. 2.5D0) GAMMAF =  1.5D0*DSQRT(PI)/2.D0
      IF (X.EQ. 3.5D0) GAMMAF =  2.5D0*1.5D0*DSQRT(PI)/2.D0
      IF (X.EQ. 3.D0 ) GAMMAF =  2.D0
      IF (X.EQ. 4.D0 ) GAMMAF = 6.D0
      IF (X.EQ. 5.D0 ) GAMMAF = 24.D0
      IF (X.EQ. 6.D0 ) GAMMAF = 120.D0
      RETURN
      END
      
      SUBROUTINE JACG (Z,NP,ALPHA,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION Z(*)
      KSTOP=10
      EPS=1.0D-12
      N   = NP-1
      DTH = 4.D0*DATAN(1.D0)/(2.D0*DBLE(FLOAT(N))+2.D0)
      DO 40 J=1,NP
      IF (J.EQ.1) THEN
      X = DCOS((2.D0*(DBLE(FLOAT(J))-1.D0)+1.D0)*DTH)
      ELSE
      X1 = DCOS((2.D0*(DBLE(FLOAT(J))-1.D0)+1.D0)*DTH)
      X2 = XLAST
      X  = (X1+X2)/2.D0
      ENDIF
      DO 30 K=1,KSTOP
      CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,NP,ALPHA,BETA,X)
      RECSUM = 0.D0
      JM = J-1
      DO 29 I=1,JM
      RECSUM = RECSUM+1.D0/(X-Z(NP-I+1))
 29   CONTINUE
      DELX = -P/(PD-RECSUM*P)
      X    = X+DELX
      IF (ABS(DELX) .LT. EPS) GOTO 31
 30   CONTINUE
 31   CONTINUE
      Z(NP-J+1) = X
      XLAST     = X
 40   CONTINUE
      DO 200 I=1,NP
      XMIN = 2.D0
      DO 100 J=I,NP
      IF (Z(J).LT.XMIN) THEN
      XMIN = Z(J)
      JMIN = J
      ENDIF
 100  CONTINUE
      IF (JMIN.NE.I) THEN
      SWAP = Z(I)
      Z(I) = Z(JMIN)
      Z(JMIN) = SWAP
      ENDIF
 200  CONTINUE
      RETURN
      END
      
      DOUBLE PRECISION FUNCTION PNORMJ (N,ALPHA,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ONE   = 1.D0
      TWO   = 2.D0
      DN    = DBLE(FLOAT(N))
      CONST = ALPHA+BETA+ONE
      IF (N.LE.1) THEN
      PROD   = GAMMAF(DN+ALPHA)*GAMMAF(DN+BETA)
      PROD   = PROD/(GAMMAF(DN)*GAMMAF(DN+ALPHA+BETA))
      PNORMJ = PROD * TWO**CONST/(TWO*DN+CONST)
      RETURN
      ENDIF
      PROD  = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE)
      PROD  = PROD/(TWO*(ONE+CONST)*GAMMAF(CONST+ONE))
      PROD  = PROD*(ONE+ALPHA)*(TWO+ALPHA)
      PROD  = PROD*(ONE+BETA)*(TWO+BETA)
      DO 100 I=3,N
      DINDX = DBLE(FLOAT(I))
      FRAC  = (DINDX+ALPHA)*(DINDX+BETA)/(DINDX*(DINDX+ALPHA+BETA))
      PROD  = PROD*FRAC
 100  CONTINUE
      PNORMJ = PROD * TWO**CONST/(TWO*DN+CONST)
      RETURN
      END