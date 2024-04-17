
!-------------------
!
!       PML
!
!------------------------   
subroutine PML(ID,NPT,dt,RT,R1,f0,VMAX,DMAX,DY,DY_G,PML1)

  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  
  DIMENSION MNO(NPT), R1(NPT)
  dimension DY_G(13,NPT), PML1(13,NPT)  !time&line(DYX)
  dimension DY(13)!think about dimension
  !---------------------------------------------------  
  dimension Yn(4),Y12(13)
  dimension DLI(4),Tn(4)
   
    !---------------
 PI=4.d0*atan(1.d0)
 D=RT
  
  a0=PI*f0;
  R=1.d-8
  d0=abs(log10(1.d0/R))*dble(DMAX)*VMAX/(2.d0*D)  !VMAX
  !====================for Test
  
  !====================================
  
    d1=d0*(R1(ID)/D)**DMAX;  
    c1=1.d0
  
    a1=a0*(1-R1(ID)/D)**DMAX; 
    !----------------------
    AAX=EXP(-(a1+d1/c1)*dt)
    BX=d1*dt/2.d0
  !=====================
  Y12=(0,0);
  !-----------------
  
   !Y12(1:4)=(DY(1:4)+DY_G(1:4,ID))/2.d0
   !PML1(1:4)=AAX*PML1(1:4)-BX*Y12(1:4)
    PML1(1:13,ID)=(AAX*PML1(1:13,ID)-BX*(DY(1:13)+AAX*DY_G(1:13,ID)))
  
  DY_G(1:13,ID)=DY(1:13)
  
  return
end subroutine
!---------------------------------------------
SUBROUTINE MS_DF2(I,K,I1,K1,N0,NX,NZ,XP,ZP,&
            NORDX,NORDZ,Mf0,IX,IDX,DEX,IZ,IDZ,DEZ,DFS)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  
  DIMENSION XP(*),ZP(3,*),ASX(NORDX),ASZ(NORDZ)
  DIMENSION IDX(NORDX),IDZ(NORDZ)
  DIMENSION DEX(NORDX),DEZ(NORDZ)
  DIMENSION DFS(2,NORDZ,4)
      
  Dimension NB(NORDX,NORDZ)
  DIMENSION DLX_LC(:),DLZ_LC(:),ID(:,:)
  dimension XE_GL(NORDX),XE_LC(NORDX),ZE_GL(NORDZ,NORDZ),DLX_LC2(NORDX),ID_GL(NORDX)
  dimension DFX(NORDX*NORDZ,NORDX*NORDZ),DFZ(NORDX*NORDZ,NORDX*NORDZ)
  Dimension AASX(NORDX),AASZ(NORDZ)
      
      ALLOCATABLE DLX_LC,DLZ_LC,ID
DEX=0.d0;DEZ=0.d0
IX=0;IZ=0;IDX=0;IDZ=0
DFS=0.d0
  
      NNX=(NX-1)*(NORDX-1)+1
      NNZ=(NZ-1)*(NORDZ-1)+1
      NP=NORDX*NORDZ
      
      I2=(I-1)*(NORDX-1)+I1
      K2=(K-1)*(NORDX-1)+K1
      
      I0=(I-1)*(NORDX-1)+1
      K0=(K-1)*(NORDZ-1)+1
      K_1=(K-2)*(NORDZ-1)+1
      
      IP=I0+I1-1; KP=K0+K1-1  !global X&Z Point value
      M0=(IP-1)*NNZ+KP    !should be same as N0: M0=N0
      !if(N0==34231)write(*,*)K0,K2
    
   !-----------------------------------------------   
      call ITYPE(IP,NORDX,NNX,ISTART,IEND,KIX)
      call ITYPE(KP,NORDZ,NNZ,KSTART,KEND,KIZ)
    
    if(MF0==1)then
        KIZ=KIZ+int(NORDZ/2)
        KSTART=KSTART-int(NORDZ/2)
        KEND=KEND-int(NORDZ/2)
        
        !print *, KIZ,KSTART,KEND
    endif
    
    if(MF0==2)then
        KIZ=KIZ-int(NORDZ/2)
        KSTART=KSTART+int(NORDZ/2)
        KEND=KEND+int(NORDZ/2)
    endif
    
    NN=0
    do II=1,NORDX
      do JJ=1,NORDZ
        NB(II,JJ)=((ISTART-1)+II-1)*NNZ+(KSTART-1)+JJ
        NN=NN+1
      enddo
      !if(MF0==1)print *, NB(II,1:NORDZ1)
    enddo
    
    do II=1,NORDX
        ID_GL(II)=((ISTART-1)+II-1)*NNZ+1
    enddo
    !write(*,*)ID_GL(1:NORDX)
    
      if(NB(KIX,KIZ)/=N0)then
      write(*,*)MF0,'ID not matched:',N0
      stop
      endif
  !----------------------------------------
    I3=ISTART
    do IX1=1,NORDX
      XE_GL(IX1)=XP(I3)
      I3=I3+1
      
      K3=KSTART
      do IZ1=1,NORDZ
        ZE_GL(IX1,IZ1)=ZP(1,NB(IX1,IZ1))
        K3=K3+1
      enddo
      !if(MF0==1)print *, ZE_GL(IX1,1:NORDZ)
    enddo
    
  !--------------------------------------
  ALLOCATE (DLX_LC(NORDX),DLZ_LC(NORDZ))
  
  XM=(XE_GL(NORDX)+XE_GL(1))/2.d0
  XD=(XE_GL(NORDX)-XE_GL(1))/2.d0
  AASX(1:NORDX)=(XE_GL(1:NORDX)-XM)/XD
  !if(NORDX/=NORDX1)XD=XD/2.d0
!=========================================
!-
!stop
IX=NORDX+NORDZ;IZ=NORDZ;
DFX=0.d0;DFZ=0.d0
No1=0
!write(*,*)AASX(1:NORDX1)
do II=1,NORDX
  call CDLI(AASX(II),NORDX,AASX,1,DLX_LC)
  
  jj0=1; jj1=NORDZ; NORDZ1=NORDZ
  if(MF0==1.or.MF0==3)then
    JJ0=4
    NORDZ1=NORDZ-3
  endif
  
  if(MF0==2)then
    JJ0=1; jj1=2;
    NORDZ1=NORDZ-3
  endif
  
  ZM=(ZE_GL(II,jj1)+ZE_GL(II,jj0))/2.d0
  ZD=(ZE_GL(II,jj1)-ZE_GL(II,jj0))/2.d0
  AASZ(1:NORDZ)=0.0
  AASZ(JJ0:jj1)=(ZE_GL(II,JJ0:jj1)-ZM)/ZD
  DZLC=ZD
  No=0
  !if(MF0==-1)print *, ZD,AASZ(JJ0:NORDZ)
  !--------------------------------------
  do JJ=JJ0,jj1
    !IP=(II-1)*NORDZ+JJ
    !IDX=0;IDZ=0
    
    call CDLI(AASZ(JJ),NORDZ1,AASZ(JJ0:jj1),1,DLZ_LC(JJ0:jj1))
    
    !if(N0==34885)write(91,*)XD,ZD
    NIP=(KIX-1)*NORDZ+KIZ
    !DXZ=0.d0
    !if(abs(ZP(2,NB2(IP)))>1.d-9)DXZ=ZP(2,NB2(IP))
    DXZ=ZP(2,NB(II,JJ))
    
    if(JJ==KIZ)then
    DFS(1,II,1)=XD*ZD
    DFS(1,II,2)=ZD/(XD*ZD)
    DFS(1,II,3)=-DXZ*XD/(XD*ZD)
    DFS(1,II,4)=XD/(XD*ZD)
    endif
    !if(N0==76127.and.NB(II,JJ)==76127)write(96,*)MF0,II,JJ,XD,ZD,DFS(1,II,1:4)
    if(II==KIX)then
    DFS(2,JJ,1)=XD*ZD
    DFS(2,JJ,2)=ZD/(XD*ZD)
    DFS(2,JJ,3)=-DXZ*XD/(XD*ZD)
    DFS(2,JJ,4)=XD/(XD*ZD)
    !write(*,*)JJ,XD,ZD
    endif
    !if(N0==76127.and.NB(II,JJ)==76127)write(96,*)MF0,II,JJ,XD,ZD,DFS(2,JJ,1:4)
    if(II==KIX.and.JJ==KIZ)then
    DEX(1:NORDX)=DLX_LC(1:NORDX)
    
    DEZ(JJ0:jj1)=DLZ_LC(JJ0:jj1)
    IX=NORDX; IZ=NORDZ
    IDX(1:NORDX)=NB(1:NORDX,KIZ)
    IDZ(1:NORDZ)=NB(KIX,1:NORDZ)
    !write(*,*)IDZ(1:NORDZ)
    No=1
    No1=No
    !goto 74
    endif
  enddo
enddo

if(No1==0)then
write(*,*)'MS_DF2: Not matched points'
stop
endif
!if(MF0==1)DEZ(:)=-DEZ(:)
!do ii=1,NORDX
!DEX(ii)=DEX(ii)*DFS(1,II,2)
!DEZ(ii)=DEZ(ii)*DFS(2,II,4)
!enddo
!write(*,*)'MS_DF2:  ',IDX(1:IX)
!write(*,*)'MS_DF2:  ',DEX(1:IX)

!=======================================================================
!----------------------      
      RETURN
      END
      
!-----------------------------------
!   ITYPE
!----------------------------------           
SUBROUTINE ITYPE(I,NORDX,NNX,ISTART,IEND,KI)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
   
    
    IMOD=mod(i-1,NORDX-1)
    IINT=i-IMOD
    
    IF(I==1) THEN
        MYTYPE=1         !left boundary
    ELSEIF (I==NNX) THEN
        MYTYPE=2         !right block
    ELSEIF (IMOD==0) THEN
        MYTYPE=3         !block boundary
    ELSE
        MYTYPE=4
    ENDIF 
    
    
    
    mid=int((NORDX-1)/2 )   
    Select case (MYTYPE)
    case(1)
        ISTART=i
        IEND=i+(NORDX-1)
        KI=1
    case(2)
        ISTART=i-(NORDX-1)
        IEND=i
        KI=NORDX
    case(3)
        ISTART=i-mid
        IEND=i+mid
        KI=mid+1
    case(4)
        ISTART=IINT
        IEND=IINT+(NORDX-1)
        KI=IMOD+1
    End select
    
    RETURN
END SUBROUTINE
!------------------------------
