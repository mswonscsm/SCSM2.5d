program Main
USE OMP_LIB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXDA=150000)
      CHARACTER LAB*80,FNAME*9
      integer :: ncore, myid, iTimes1, iTimes2, rate  
      
      DIMENSION XTO(:,:),ZTO(:,:),X(:),NZL(:)
      dimension XP(:),ZP(:,:)

      Dimension IDG(:)
      DIMENSION PM(:,:),PT1(:,:),PT11(:),PT2(:),PT_TEMP(:,:)
      dimension INDP(:),VMAX(:),VMIN(:),KT(:),MT(:),ITTI(:),ZTA(:)
      dimension XSR(:),ZSR(:)
      dimension FKY(:),WKY(:), FKP(:),WKP(:),Nani0(:),Nani(:)
      dimension tau_ep0(:,:,:), tau_sig0(:,:,:)
      dimension tau_ep1(:,:,:),tau_sig1(:,:,:)
      dimension tau_ep(:,:,:),tau_sig(:,:,:),tau_TTI(:,:,:)
      
      dimension FTS(:,:),FTS1(:,:)
      dimension FTG(:,:,:),FTG1(:,:,:),FST(:,:,:)
      dimension FGT(:,:,:,:)
      dimension Relax_C(:)
      !Integer, External :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
      
      ALLOCATABLE XTO,ZTO,IDG,VMAX,VMIN
      ALLOCATABLE X,NZL,KT,MT
      ALLOCATABLE XP,ZP,PM,PT1,PT2,INDP,PT11,ITTI,ZTA,PT_TEMP
      Allocatable XSR,ZSR
      allocatable FKY,WKY,FKP,WKP
      allocatable tau_ep1,tau_sig1,tau_ep0,tau_sig0,tau_ep,tau_sig,tau_TTI,Relax_C, Nani0,Nani
      allocatable FTS,FTS1, FTG,FTG1,FST,FGT
      
      WRITE(*,*)' '
      WRITE(*,*)'******************************************************'
      WRITE(*,*)'*   2.5D Subdomain Chebyshev Spectral Method (SCSM)  *'
      WRITE(*,*)'*   for seismic wavefield modelling in anisotripic   *'
      WRITE(*,*)'*   media (velocity-stress version).                 *'
      WRITE(*,*)'******************************************************'
        
      
      OPEN(1,FILE='2.5Dseis_SCSM.inp',STATUS='UNKNOWN')
      
      OPEN(88,FILE='Xgrid.out',STATUS='UNKNOWN')
      OPEN(89,FILE='Zgrid.out',STATUS='UNKNOWN')

      OPEN(91,FILE='Test1.out',STATUS='UNKNOWN')
      OPEN(92,FILE='Test2.out',STATUS='UNKNOWN')
      OPEN(93,FILE='Test3.out',STATUS='UNKNOWN')
      OPEN(94,FILE='Test4.out',STATUS='UNKNOWN')
      OPEN(95,FILE='Test5.out',STATUS='UNKNOWN')
      OPEN(96,FILE='Test6.out',STATUS='UNKNOWN')
      !-------------------------------

!==========================================
!
!  1) Input Data
!
!=========================================
    WRITE(*,*)' 1.  Input Data'
!C----------------------------------------------------------------------C
!C     (0) input subdomain size & Point Num                             C
!C----------------------------------------------------------------------C
      WRITE(*,*)'     1.0. NT,Source Data'
      READ(1,*)LAB
      READ(1,*)MODD,NT,f0,t0

!C----------------------------------------------------------------------C
!C     (1) input subdomain size & Point Num                             C
!C----------------------------------------------------------------------C
      WRITE(*,*)'     1.1. Input DXZ,NORDXZ'
      READ(1,*)LAB
      READ(1,*)DX,DZ,NORDX,NORDZ

!C----------------------------------------------------------------------C
!C     (2) input 2-D grid interfaces
!C----------------------------------------------------------------------C
      !input the interfaces
      WRITE(*,*)'     1.2. Input free-srface & interfaces'
      READ(1,*)LAB
      READ(1,*)NSF0,NSP,NFS

      !add interfaces for extension
      NSF=NSF0+2
      if(NFS==1)NSF=NSF0+1

      ALLOCATE (XTO(NSF,NSP),ZTO(NSF,NSP))

      XMIN= 1.D+6
      XMAX=-1.D+6
      ZMIN= 1.D+6
      ZMAX=-1.D+6

      DO I=1,NSF0
        READ(1,*)LAB
        DO J=1,NSP
          READ(1,*)XTO(I+1,J),ZTO(I+1,J)
          IF(XTO(I+1,J).LT.XMIN)XMIN=XTO(I+1,J)
          IF(XTO(I+1,J).GT.XMAX)XMAX=XTO(I+1,J)
          IF(ZTO(I+1,J).LT.ZMIN)ZMIN=ZTO(I+1,J)
          IF(ZTO(I+1,J).GT.ZMAX)ZMAX=ZTO(I+1,J)
        ENDDO
      ENDDO
    

!--------------------------------------
!   (3) Input Model Parameters
!---------------------------------------
    WRITE(*,*)'     1.3. Input Model'
    NMD=NSF0-1
    allocate(PT1(NMD,7),INDP(NMD),ITTI(NMD),ZTA(NMD))
    
    write(*,*)'No. Medium:',NMD
    READ(1,*)LAB
    READ(1,*)IEXT

    do i=1,NMD
      PT1(i,:)=0.d0
      READ(1,*)LAB
      READ(1,*)INDP(i)
      READ(1,*)PT1(i,1:INDP(i)+1)
      READ(1,*)ITTI(i)
      READ(1,*)ZTA(i)
      
    enddo

!--------------------------------------
!   (4) Source & Receiver
!---------------------------------------
    write(*,*)'     1.4. Input S&R vector'

    Read(1,*)LAB
    Read(1,*)NGV,MX,MZ
    
!--------------------------------------
    WRITE(*,*)'     1.5. Input Absorbing index & Parameter '
      READ(1,*)LAB
      READ(1,*)MAB,DMAX,Efact

!-----------------------------------
    Write(*,*)'     1.6. Input S & R'
    read(1,*)LAB
    READ(1,*)NSR
    allocate(XSR(NSR),ZSR(NSR))
    
    do i=1,NSR
        read(1,*)XSR(i),ZSR(i)
    enddo
!-------------------------------------
!-End of Input
  CLOSE(1)
  WRITE(*,*)'   End of Input data'
  WRITE(*,*)' ----------------------------------------------'
!===========================================================================

!==========================================
!---------------------------------------
!  2) Calculating if proper Grid-Model
!------------------------------------
!=========================================
  WRITE(*,*)' 2.  Grid<->Model Calculation'
    allocate(VMAX(NMD),VMIN(NMD))
    VTX=0.d0;VTN=1.d6

    do i=1,NMD
    if(INDP(i)==1)then
      VMAX(i)=Sqrt(PT1(i,2)/PT1(i,1))
      VMIN(i)=VMAX(i)
    endif

    if(INDP(i)==2)then
      VMAX(i)=sqrt((PT1(i,2)+2*PT1(i,3))/PT1(i,1))
      VMIN(i)=sqrt((PT1(i,3))/PT1(i,1))
    endif
    
    if(INDP(i)>=3)then
      Cmin=1.d15
      Cmax=-1.d15
      do j=1,INDP(i)
        if(abs(PT1(i,1+j))>1.d-6)then
        if(PT1(i,1+j)<Cmin)Cmin=abs(PT1(i,1+j))
        if(PT1(i,1+j)>Cmax)Cmax=abs(PT1(i,1+j))
        endif
      enddo
      
      VMIN(i)=sqrt(Cmin/PT1(i,1))
      VMAX(i)=sqrt(Cmax/PT1(i,1))
    endif

  write(*,*)'   VMAX=',VMAX(i),'  VMIN=',VMIN(i)
  if(VTX<VMAX(i))VTX=VMAX(i)
  if(VTN>VMIN(i))VTN=VMIN(i)
  !------------------------------------

    WLMAX=VMAX(i)/f0;  WLMIN=VMIN(i)/f0 !Wavelength
  write(*,*)'   WLMAX=',WLMAX,'  WLMIN=',WLMIN

    if(WLMIN<1.0d0*int(DX))then
      write(*,*)'   Subdomain Size Too large'
      write(*,*)'   Recommandation: DX <',WLMIN/2
      goto 10
    endif

    if(WLMAX>XMAX)then
      write(*,*)'   !!Grid Size Too small!!'
      write(*,*)'   !!Recommandation: Tot/2 >',2*WLMAX
      goto 10
    endif

    enddo
    !--------------------------
    WLMAX=VTX/f0;
    PI=4.d0*atan(1.d0)
    if(MAB==1)then
    NEXTD=INT(3.d0*WLMAX/DX+1) !ABL
    else
    NEXTD=INT(3.d0*WLMAX/DX+1)
    endif
    
      write(*,*)' NEXTD=',NEXTD
    WRITE(*,*)'   End of Calculation'
    WRITE(*,*)' --------------------------------'
!===============================================
!---------------------------------------------

!------------------------------------------
!
! 3)  Generating Grid
!
!-----------------------------------
    NX=INT((XMAX-XMIN)/DX)+1+2*NEXTD
    NZ=INT((ZMAX-ZMIN)/DZ)+1+(1+(1-NFS))*NEXTD
    NNX=(NX-1)*(NORDX-1)+1
    NNZ=(NZ-1)*(NORDZ-1)+1
    NPT=NNX*NNZ

    Allocate(X(NX),XP(NNX),ZP(3,3*NPT))
    Allocate(NZL(NSF-1))

    WRITE(*,*)' 3. Generating Grid'
CALL  GRID_2D(DX,DZ,NEXTD,NORDX,NORDZ,NSF,NSP,NFS,XTO,ZTO,XMIN,XMAX,&
    ZMIN,ZMAX,NX,NZ,NZL,XP,ZP,NNX,NNZ)

    deallocate(X)

    NPT=NNX*NNZ
    write(*,*)'   NX=',NX,' NZ=',NZ
    write(*,*)'   NNX=',NNX,' NNZ=',NNZ

    RM=1.d6
    do i=1,NNX
      do j=1,NNZ-1
        RK=ZP(1,(i-1)*NNZ+j+1)-ZP(1,(i-1)*NNZ+j)
        if(RM>RK)RM=RK
      enddo
    enddo
    write(*,*)' MIN Distance:',RM
    DTR=0.35d0*RM/VTX
  !===============
  !-DT
    dt=2.d-4
  !===============
  
  !===============
  !-Kyc
    !dkyc=pi/RM
    dkyc=1.5d0*2.d0*PI*f0/VTN;
    write(*,*)' Dkyc=',dkyc
  !===============
     if(dt>DTR)then
      write(*,*)' !!Smaller dt needed <',DTR
      goto 10
    endif

    WRITE(*,*)'   Grid ends'
    WRITE(*,*)' -----------------------------'
    
!--------------------
!-----------------------------------------
!
! 4)  Model: Density and Elastic moduli
!
!------------------------------------
      WRITE(*,*)' 4. Generating Model'
      allocate(PT11(7),PT2(21),PM(22,NNX*NNZ),MT(NMD))
      PM(1:22,:)=0.d0

!------------------------
    do i=1,NMD
      if(INDP(i)>1)then
        mno=1;
        if(INDP(i)==2)write(*,*)'   ISO_S'
        
        if(INDP(i)==4)write(*,*)'   VTI1_S'
        if(INDP(i)==3)write(*,*)'   VTI2_S'
        if(INDP(i)==5)write(*,*)'   VTI_S'
        
        if(INDP(i)==6)write(*,*)'   ORT'
        if(ITTI(i)==1)write(*,*)'    TTI'
      else
        mno=3;
        write(*,*)'   Water'
      endif
      MT(i)=mno
    enddo

  !================================
  allocate(KT(NMD))

    NMZ=NZL(1)
    do k=1,NMD
      NMZ=NMZ+NZL(k+1)
      KT(k)=NMZ*(NORDZ-1)+1
    enddo
  
  write(*,*)'Number of Subdomains in Layers:',NZL(:)
  write(*,*)'Nth point of Interfaces:',KT(1:NMD)
  

  do i=1,NNX
  do j=1,NNZ
    if(j<=KT(1))mn=1
    if(j>KT(NMD))mn=NMD

    if(NMD>1)then
      do k=1,NMD-1
        if(KT(k)<j.and.j<=KT(K+1))mn=k+1
      enddo
    endif

    !------------------------
    ID=(i-1)*NNZ+j
    PM(1,ID)=PT1(mn,1)

    PT11=0.d0;  PT2=0.d0

    PT11(1:INDP(mn))=PT1(mn,2:INDP(mn)+1)
    IDP=INDP(mn)
    
    !TTI case 
    ZTA0=0.0
    if(ITTI(mn)==1)ZTA0=tan(pi*ZTA(mn)/180.0)
    call TTI(INDP(mn),PT11,ITTI(mn),ZTA0,PT2)

    PM(2:22,ID)=PT2(1:21)
  enddo
  enddo
  
 
! Write orinigal model 
   ! do ii=1,NNX
   !     write(91,*)PM(1,1+(ii-1)*NNZ:ii*NNZ)
   !     write(92,*)PM(2,1+(ii-1)*NNZ:ii*NNZ)
   !     write(93,*)PM(4,1+(ii-1)*NNZ:ii*NNZ)
   !     write(94,*)PM(13,1+(ii-1)*NNZ:ii*NNZ)
   !     write(95,*)PM(17,1+(ii-1)*NNZ:ii*NNZ)
   !     write(96,*)PM(22,1+(ii-1)*NNZ:ii*NNZ)
   ! enddo
      
!stop
!--------------------------------------------
!Read Relx time input data layer by layer
Print *, 'Read Relax Times'
I_VISCO=0; Nsls=0; IREC=0 ; 
OPEN(2,FILE='relaxation_time.inp',STATUS='UNKNOWN')
read(2,*)LAB
!Update method 0: ODE, 1: Recursive convolution
Read(2,*)I_VISCO, IREC, Nhete
if(I_VISCO==1)then
    Allocate(Nani0(Nhete))
    READ(2,*)Nani0(1), Nsls
    Allocate(tau_ep0(Nhete,Nani0(1),Nsls), tau_sig0(Nhete,Nani0(1),Nsls)) 
    tau_ep0(:,:,:)=0.d0; tau_sig0(:,:,:)=0.d0
do i=1,Nhete
    if(I>1)READ(2,*)Nani0(I), Nsls
    print *,Nani0(I), Nsls
do j=1,Nani0(I)
    READ(2,*)tau_ep0(i,j,1:Nsls)
    READ(2,*)tau_sig0(i,j,1:Nsls)
    print *,tau_ep0(i,j,1:Nsls)
    print *,tau_sig0(i,j,1:Nsls)
enddo
enddo
close(2)

Nani_max=maxval(Nani0(1:Nhete))
endif


Allocate(tau_ep1(Nani_max,NSLS,NPT), tau_sig1(Nani_max,NSLS,NPT))
Allocate(tau_ep(NPT,NSLS,21), tau_sig(NPT,NSLS,21),Nani(NPT))
if(I_VISCO==1)then
  do i=1,NNX
  do j=1,NNZ
    if(j<=KT(1))mn=1
    if(j>KT(NMD))mn=NMD

    if(NMD>1)then
      do k=1,NMD-1
        if(KT(k)<j.and.j<=KT(K+1))mn=k+1
      enddo
    endif

    !------------------------
    ID=(i-1)*NNZ+j
    !PM(1,ID)=PT1(mn,1)
    Nani(ID)=Nani0(mn)
    tau_ep1(1:Nani_max,1:NSLS,ID)=tau_ep0(mn,1:Nani_max,1:Nsls)
    tau_sig1(1:Nani_max,1:NSLS,ID)=tau_sig0(mn,1:Nani_max,1:Nsls)
    !-----------------------
  enddo
  enddo
  
    do ii=1,npt
    if (Nani(ii)==1)then !Viscoacoustic
    TAU_EP(ii,1:NSLS,1)=TAU_EP1(1,1:NSLS,ii)
    TAU_SIG(ii,1:NSLS,1)=TAU_SIG1(1,1:NSLS,ii)
    
    else!ISOTROPIC
    do jj=1,NANI(ii)
    TAU_EP(ii,1:NSLS,jj)=TAU_EP1(jj,1:NSLS,ii);
    TAU_SIG(ii,1:NSLS,jj)=TAU_SIG1(jj,1:NSLS,ii);
    enddo
    endif 
    
    enddo    

endif
Print *, 'END Read Relax Times'

!-TTI for relax
allocate(tau_TTI(NPT,NSLS,21),RELAX_C(6))
tau_TTI(1:NPT,1:NSLS,1:21)=0.d0
  do i=1,NNX
  do j=1,NNZ
    if(j<=KT(1))mn=1
    if(j>KT(NMD))mn=NMD

    if(NMD>1)then
      do k=1,NMD-1
        if(KT(k)<j.and.j<=KT(K+1))mn=k+1
      enddo
    endif

    !------------------------
    ID=(i-1)*NNZ+j
    
    
    if(I_VISCO==1)then
    IDP=INDP(mn)
    do k=1,IDP
    PT11(k)=1.d0+sum(TAU_EP(ID,1:NSLS,k)/TAU_SIG(ID,1:NSLS,k)-1.d0)/dble(NSLS)
    Relax_C(k)=PT1(mn,k+1)/PT11(k)
    enddo
    
    
    do k=1,NSLS
    PT11=0.d0;  PT2=0.d0
    PT11(1:IDP)=(1.d0-TAU_EP(ID,k,1:IDP)/TAU_SIG(ID,k,1:IDP))/TAU_SIG(ID,k,1:IDP)*RELAX_C(1:IDP)
    !TTI case 
    if(ITTI(mn)==1)then
        ZTA0=tan(pi*ZTA(mn)/180.0)
    else
        ZTA0=0.d0
    endif
    
    call TTI(IDP,PT11,ITTI(mn),ZTA0,PT2)

    tau_TTI(ID,k,1:21)=PT2(1:21)
    enddo
    endif
    !-----------------------
  enddo
  enddo

!stop

!-Marmousi Model input
IF(IEXT==1)then
allocate(PT_TEMP(NPT,INDP(1)+1))

call EXTERNAL_GRID_MODEL(NX,NZ,DX,NORDX,NORDZ,NNX,NNZ,NEXTD,NFS,INDP(1),ITTI(1),PT_TEMP,PM)

if(I_VISCO==1)then

call EXTERNAL_GRID_Relax_time(NX,NZ,DX,NORDX,NORDZ,NNX,NNZ,NEXTD,NFS,&
                        ITTI(1),PT_TEMP,nani,nsls,tau_ep,tau_sig,tau_TTI)
                        
endif
print *, 'EXTERNAL GRID ENDs'
endif

!IF(IEXT==2)then
!allocate(PT_TEMP(NPT,INDP(1)+1))

!call EXTERNAL_GRID_MODEL2(NNX,NNZ,INDP(1),ITTI(1),PT_TEMP,PM)

!if(I_VISCO==1)then

!call EXTERNAL_GRID_Relax_time(NX,NZ,DX,NORDX,NORDZ,NNX,NNZ,NEXTD,NFS,&
!                        ITTI(1),PT_TEMP,nani,nsls,tau_ep,tau_sig,tau_TTI)
!endif
!print *, 'EXTERNAL GRID ENDs'
!endif

!stop
! Write orinigal model 
!    do ii=1,NNX
!        write(91,*)PM(1,1+(ii-1)*NNZ:ii*NNZ)
!        write(92,*)PM(2,1+(ii-1)*NNZ:ii*NNZ)
!        write(93,*)PM(4,1+(ii-1)*NNZ:ii*NNZ)
!        write(94,*)PM(13,1+(ii-1)*NNZ:ii*NNZ)
!        write(95,*)PM(17,1+(ii-1)*NNZ:ii*NNZ)
!        write(96,*)PM(22,1+(ii-1)*NNZ:ii*NNZ)
!    enddo
  !==================================
  !   Source & Receivers ID
  allocate(IDG(NSR))
  EX=NEXTD*DX; EZ=NEXTD*DZ
  do k=1,NSR
    XEX=XSR(k); ZEX=ZSR(k)+EZ
    !write(*,*)XEX,ZEX
    x1=XP(1)
    do i=1,NX-1
        if (abs(XEX-X1) <= DX/2.0+1.d-3) then
            i1=i
            exit
        endif
        x1=x1+DX
    enddo
    
    if (abs(XEX-X1) <= 0.1) then
     MX1=1!-int(nordx/2)+1
    else
     MX1=int(nordx/2)+1
    endif
    !print *, k, XEX,x1, MX1
    z1=ZP(1,1)
    do j=1,NZ
         if(abs(ZEX-Z1) <= DZ) then
            j1=j
            exit
        endif
        z1=z1+DZ
   
   enddo
   
    if (abs(ZEX-Z1) <= 0.1) then
     MZ1=1!-int(nordx/2)+1
    else
     MZ1=-int(nordz/2)+1
    endif
  if(k==1)then  !source
  IDG(k)=((i1-1)*(NORDX-1)+MX1-1)*NNZ+(j1-1)*(NORDZ-1)+MZ1
  !IDG(k)=((i1-1)*(NORDX-1)+MX1-1)*NNZ+NNZ
  else  !receiver
  IDG(k)=((i1-1)*(NORDX-1)+MX1-1)*NNZ+(j1-1)*(NORDZ-1)+MZ1
  !IDG(k)=((i1-1)*(NORDX-1)+MX1-1)*NNZ+NNZ
    endif
enddo
  !stop
  !==================================

    Deallocate(PT1,PT2,PT11)
    WRITE(*,*)'   Model ends'
    WRITE(*,*)' -----------------------------'


!stop
!------------------------------------------
!
!         C_DF: Main Calculation
!   (Differential multiplication& T-step
!
!------------------------------------------
      ALLOCATE(FTG(NSR,13,NT),FTG1(NSR,13,NT),FTS(4,NT),FTS1(4,NT))
      Allocate(FGT(8,NSR,13,NT),FST(8,4,NT))
      
      FTG=0.0;FTS=0.0
      FGT=0.0;FST=0.0
      WRITE(*,*)' 5. C_DF Going in'
    
    call SYSTEM_CLOCK(count_rate=rate)
    call SYSTEM_CLOCK(iTimes1)
  
  !=2D case
if(MODD==0)then
    Nky=1
  CALL C_DF(MODD,NT,Nky,1,dkyc,0.d0,f0,t0,dt,NX,NZ,NORDX,NORDZ,npt,NFS,IDG,NGV,&
      DX,DZ,NEXTD,XP,ZP,NMD,KT,MT,PM,&
      MAB,DMAX,CMAX,Efact,VTX,VTN,FTG1,NSR,I_VISCO, IREC,Nani, Nsls,tau_ep, tau_sig,tau_TTI)

else

  Nky=1+int(1.2*dkyc*(vtx*(nt-1)*dt+Wlmax)/(2*pi))
  dky=dkyc/dble(NKY-1)
  
  
  allocate(FKY(NKY),WKY(NKY),FKP(NKY),WKP(NKY))
  CALL GLL(NKY,FKP,WKP)
    II=0
    DO I=1,NKY
      II=II+1
      FKY(II)=0.5*Dkyc*(FKP(I)+1.)
      WKY(II)=0.5*Dkyc*WKP(I)
    ENDDO
      
  write(*,*)'Nky=',Nky,'    Dky=',dky

!=2.5D case
ncore=10
call omp_set_num_threads(ncore)
!$omp Parallel DEFAULT(SHARED) PRIVATE(wky1,FKY1,myid,FTG1)
  myid=OMP_GET_THREAD_NUM()
    
!$OMP DO
  do ik=1,NKY

    
    !$OMP CRITICAL
    FTG1(1:NSR,1:13,1:NT)=0.d0
    
    !Gauss-Quadrature
    FKY1=FKY(ik);
    Wky1=WKY(ik);
    write(*,*)ik,'/',NKY,'|',FKY1,'/',DKYC,'|','  Thread ID:',myid
    !$OMP END CRITICAL
    
CALL C_DF(MODD,NT,Nky,ik,dkyc,FKY1,f0,t0,dt,NX,NZ,NORDX,NORDZ,npt,NFS,IDG,NGV,&
DX,DZ,NEXTD,XP,ZP,NMD,KT,MT,PM,&
MAB,DMAX,CMAX,Efact,VTX,VTN,FTG1,NSR,I_VISCO, IREC,Nani, Nsls,tau_ep, tau_sig,tau_TTI)

    !$OMP CRITICAL
    FTG=FTG+WKY1/PI*FTG1*cos(0.d0*FKY1);
    FTG1(1:NSR,1:13,1:NT)=0.d0
    !$OMP END CRITICAL
    
  enddo;
!$omp enddo


!$OMP END parallel
!==========================================
endif

call SYSTEM_CLOCK(iTimes2)

    
  !=====================================
    WRITE(*,*)'C_DF ends'
      !--------------------------
      !-for output out of C_DF

do ii=1,NSR
    if(MODD==1)then
        call output_GP(ii,NSR,NT,dt,FTG)
    else
        call output_GP(ii,NSR,NT,dt,FTG1)
    endif

enddo
      !do ii=1,NNX
      ! write(91,*)ZP(1,1+(ii-1)*NNZ:ii*NNZ)
       ! write(91,*)ZP(2,1+(ii-1)*NNZ:ii*NNZ)
      !  write(93,*)PM(1,1+(ii-1)*NNZ:ii*NNZ)
      !  write(94,*)PM(2,1+(ii-1)*NNZ:ii*NNZ)
      !enddo
      !-----------------------
      deALLOCATE(FTG,FTG1,FTS,FTS1)
      Deallocate(PM,XP,ZP)
      Deallocate (XTO,ZTO)
      deallocate (NZL)
    !call CPU_TIME(FINISH)
    
    write(*,*)'TIME= ', real(iTimes2-iTimes1)/real(rate)/60.d0,'min'
10      write(*,*)'Program ends'
ENDprogram

!============================================
subroutine output_GP_TEST(i,FKY,NSR,NT,dt,FTG1)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
character(len=10) :: char_i1,char_i2
dimension FTG1(NSR,13,NT)

double precision,Dimension(:) :: T(NT)

do n=1,nt
    T(n)=dble(n-1)*dt
enddo

write(char_i1, '(I10)') I
write(char_i2, '(F6.3)') FKY
open(86, file='rec_real_'//trim(adjustl(char_i1))//'_'//trim(adjustl(char_i2))//'.out')
write(86,*)T(1:NT)
do ii=1,13
write(86,*)real(FTG1(i,ii,1:NT))
enddo
close(86)


end subroutine

!-------------------------------------
subroutine output_GP(i,NSR,NT,dt,FTG)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
character(len=5) :: char_i
dimension FTG(NSR,13,NT)

double precision,Dimension(:) :: T(NT)

do n=1,nt
    T(n)=dble(n-1)*dt
enddo

write(char_i, '(I5)') i
open(86, file='rec_real_'//trim(adjustl(char_i))//'.out')
write(86,*)T(1:NT)
do ii=1,13
write(86,*)real(FTG(i,ii,1:NT))
enddo
close(86)

end subroutine
