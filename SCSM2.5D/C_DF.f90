SUBROUTINE C_DF(MODD,NTD,Nky,ik,dkyc,wky,f0,t0,dt,NX,NZ,NORDX,NORDZ,NPT,NFS,IDG,NGV,DX,DZ,NEXTD,XP,ZP,NMD,KT,MT,PM0,&
                  MAB,DMAX,CMAX,Efact,VMAX,VMIN,FTG,NSR,I_VISCO, IREC,Nani, Nsls,tau_ep, tau_sig,tau_TTI)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    CHARACTER LAB*80
!-INPUT
  DIMENSION XP((NX-1)*(NORDX-1)+1),ZP(3,NPT)
  DImension KT(NMD),IDG(NSR),MT(NMD),Nani(NPT)
  dimension PM0(22,NPT)
  dimension tau_ep(NPT,nsls,21),tau_sig(NPT,nsls,21),tau_TTI(NPT,nsls,21)
  Dimension FTG(NSR,13,NTD)
!==================================================================

!-Dimension
INTEGER, DIMENSION(:),allocatable :: MNOG(:)
INTEGER, DIMENSION(:,:),allocatable:: IDXG(:,:),IDZG(:,:), IDZ_WS(:,:)
DIMENSION DEXG(:,:),DEZG(:,:),DEZ_WS(:,:), TEST(:,:)
DIMENSION DFSG(:,:,:,:),DFS_WS(:,:,:,:)
DIMENSION IDX(NORDX),DEX(NORDX),IDZ(NORDZ),DEZ(NORDZ)
DIMENSION DFS(2,NORDX,4)
dimension PM(22,NPT), GF(13,(NX-1)*(NORDX-1)+1)
DIMENSION PD(22), WK(NORDX)
DIMENSION ETA_G(:),RX0(:),RZ0(:)

Dimension Yhat_l(:,:,:)

DIMENSION  AT_Memory(NPT,NSLS,6,6),TSRC_ABC(NPT,NSLS,4)
DIMENSION FORCE1(3),DX_Y(13),DZ_Y(13),DY_Y(13),DX_Y_WS(13),DZ_Y_WS(13),DY_Y_WS(13)
DIMENSION V_AG(:,:,:),ST_AG(:,:,:,:), P_AG(:,:,:,:)

doubleprecision,dimension(:), allocatable :: P,P_N
doubleprecision,dimension(:,:),allocatable:: V,ST1, ST2, V_N, ST1_N, ST2_N,WV,WV_N

doubleprecision,dimension(:,:),allocatable :: P_T_G
doubleprecision,dimension(:,:,:),allocatable :: V_T_G, ST_T_G, WV_T_G 

doubleprecision,DIMENSION(:,:,:) :: var_memory(:,:,:),Var_memory_T(:,:,:),var_memory_WS(:,:,:)
doubleprecision,DIMENSION(:,:,:,:) :: Var_memory_T2(:,:,:,:)
doubleprecision,DIMENSION(:,:) :: V_kl(:,:), V_KL_WS(:,:)
doubleprecision,DIMENSION(:,:) :: DX_Y_G(:,:),DZ_Y_G(:,:), DY_Y_G(:,:)
doubleprecision,DIMENSION(:,:),allocatable:: PMLX_G, PMLZ_G

!-Allocatable
allocatable DEXG, DEZG, DFSG, DEZ_WS,DFS_WS
allocatable V_AG,ST_AG,P_AG
allocatable ETA_G,RX0,RZ0
Allocatable TEST
allocatable Yhat_l
allocatable var_memory,var_memory_T,V_kl,var_memory_t2,var_memory_WS,V_KL_WS

allocatable DX_Y_G,DZ_Y_G,DY_Y_G
!--Calculate Grid Size
    PI=4.d0*atan(1.d0)
    NNX=(NX-1)*(NORDX-1)+1
    NNZ=(NZ-1)*(NORDZ-1)+1
    NPT=NNX*NNZ
    NP=NORDX*NORDZ
    !print *,NPT
!-Allocate
allocate(MNOG(NPT))
allocate(IDXG(NPT,NORDX),DEXG(NPT,NORDX),IDZG(NPT,NORDZ),DEZG(NPT,NORDZ))
allocate(DFSG(NPT,2,NORDX,4))
allocate(IDZ_WS(NNX,NORDZ),DEZ_WS(NNX,NORDZ),DFS_WS(NNX,2,NORDX,4))
allocate(ST_AG(NPT,3,3,6),V_AG(NPT,6,6))
allocate(P_AG(NPT,3,4,4))

Allocate(V(3,NPT),ST1(3,NPT),ST2(3,NPT),V_N(3,NPT),ST1_N(3,NPT),ST2_N(3,NPT))
Allocate(P(NPT), WV(3,NPT), P_N(NPT),WV_N(3,NPT))
Allocate(V_T_G(3,5,NPT),ST_T_G(6,5,NPT), P_T_G(5,NPT), WV_T_G(3,5,NPT))
Allocate(ETA_G(NPT),RX0(NPT),RZ0(NPT))
Allocate(TEST(6,NPT))
ALLOCATE(DX_Y_G(13,NPT),DZ_Y_G(13,NPT),PMLX_G(13,NPT),PMLZ_G(13,NPT))
ALLOCATE(DY_Y_G(13,NPT))
!---------

PM(1:22,1:NPT)=PM0(1:22,1:NPT)
!Update method 0: ODE, 1: Recursive convolution
if(I_VISCO==1.and.Nsls>0)then
    if(MODD==0)print *, 'Viscoacoustic & Viscoelastic case.'
    Allocate(Yhat_l(NPT,Nsls,21))
    
    if(MODD==0)print *, 'Standard Linear Solid USE.'
    if(MODD==0)print *, 'NSLS=',Nsls
    Allocate(var_memory(NPT,Nsls,7) ,var_memory_T(NPT,Nsls,7),V_kl(NPT,6),var_memory_T2(2,NPT,Nsls,7))
    Allocate(V_kl_WS(NNX,6),var_memory_WS(NNX,Nsls,7))
    

Yhat_l=0.d0; 
do ii=1,npt
if (Nani(ii)==1)then !Viscoacoustic
    Yhat_l(ii,:,1)=(tau_ep(ii,:,1)/tau_sig(ii,:,1)-1.d0)/Nsls
else if(Nani(ii)==2)then    !ISOTROPIC
    !Lambda
    Yhat_l(ii,:,1:Nani(ii))=(tau_ep(ii,:,1:Nani(ii))/tau_sig(ii,:,1:Nani(ii))-1.d0)/dble(Nsls)
endif 
enddo

!stop

var_memory(:,:,:)=0.d0;var_memory_T(:,:,:)=0.d0; V_kl(:,:)=0.d0; var_memory_T2(:,:,:,:)=0.d0;
V_KL_WS(:,:)=0.d0; var_memory_WS(:,:,:)=0.d0;
endif

!print *, Yhat_l
RS=1.25*DX+1.d-6    !Source distribution radius
!-Abs distance & factors
    RT=dfloat(NEXTD)*DX!*dsqrt(2.d0)
    if(MAB==2)then
        EMAX=1.d0*2.d0*PI*f0;
        Amax=log(Efact)
        Amax2=log(5.d0)
    endif
    
!Matrix Init
MNOG(:)=0; ETA_G(:)=0.d0
IDXG(:,:)=0; DEXG(:,:)=0.d0; IDZG(:,:)=0; DEZG(:,:)=0.d0
DFSG(:,:,:,:)=0.d0; 
FTG(:,:,:)=0;

P_T_G=0; WV_T_G=0; V_T_G=0; ST_T_G=0
JWS=10**6
!-snapshot
ifile=1

I_Gauss_distribution=1 !-Source distribution
RX0(:)=0.d0; RZ0(:)=0.d0
DX_Y_G(:,:)=0.d0; DZ_Y_G(:,:)=0.d0; DY_Y_G(:,:)=0.d0
PMLX_G(:,:)=0.d0; PMLZ_G(:,:)=0.d0; 

!===============================
!-Before Chebyshev, DF & Matrix
if(MODD==0)write(*,*)'C_DF: Before T-step'
DO 14 I=1,NX-1
DO 15 K=1,NZ-1
NO=(I-1)*(NORDX-1)*NNZ+(K-1)*(NORDZ-1)+1
NORD1=NORDX-1; NORD3=NORDZ-1
IF(I.EQ.NX-1)NORD1=NORDX
IF(K.EQ.NZ-1)NORD3=NORDZ

!--------------------------
    DO 16 I1=1,NORD1
    DO 17 K1=1,NORD3
    ITX=(I-1)*(NORDX-1)+I1
    ITZ=(K-1)*(NORDZ-1)+K1
    !-------------------------
    ID=NO+(I1-1)*NNZ+(K1-1)    !global index
    IT=(I1-1)*NORDZ+K1         !local  index
    
    !Model & layer define
    call MODEL_CDF(NFS,ITZ,NMD,KT,MT,MNO,K_layer,MFS) 
    
    KM_visco=K_layer
    MNOG(ID)=MNO
    !-Identify the SR
    if(ID==IDG(1))then
      
      xs=XP(ITX); zs=zp(1,ID)
      if(wky<=abs(1.d-9).and.wky>=abs(0.d0))then
      write(*,*)'----------------------------'
      Write(*,*)'   NF X-axis:',I,I1,'->',ITX,XP(ITX)
      Write(*,*)'   NF Z-axis:',K,K1,'->',ITZ,ZP(1,ID)
      Write(*,*)'   Source Medium:',mno,ID
      write(*,*)'----------------------------'
      endif
      endif
     
    do nr=2,NSR
      if(wky<=abs(1.d-9).and.wky>=abs(0.d0))then
      if(ID==IDG(nr))then
      write(*,*)'----------------------------'
      Write(*,*)'   R X-axis:',I,I1,'->',ITX,XP(ITX)
      Write(*,*)'   R Z-axis:',K,K1,'->',ITZ,ZP(1,ID)
      Write(*,*)'   Receiver Medium:',mno,ID
      write(*,*)'----------------------------'
      endif
      endif
    enddo
    
    !-Chebyshev Differentiation
    MF0=0
    KFS=NNZ
    IF(MNO==1.and.NFS==1.and.ITZ==KFS)then
        MF0=3
    endif
    
    if(MNO==2)MF0=1

    call MS_DF2(I,K,I1,K1,ID,NX,NZ,XP,ZP,&
            NORDX,NORDZ,MF0,IX,IDX,DEX,IZ,IDZ,DEZ,DFS)
    IDXG(ID,1:NORDX)=IDX(1:NORDX);DEXG(ID,1:NORDX)=DEX(1:NORDX)
    IDZG(ID,1:NORDZ)=IDZ(1:NORDZ);DEZG(ID,1:NORDZ)=DEZ(1:NORDZ)
    DFSG(ID,1,1:NORDX,1:4)=DFS(1,1:NORDX,1:4)
    DFSG(ID,2,1:NORDZ,1:4)=DFS(2,1:NORDZ,1:4)
    
    
    if(MNO==2)then
        MF0=2
        call MS_DF2(I,K,I1,K1,ID,NX,NZ,XP,ZP,&
            NORDX,NORDZ,MF0,IX,IDX,DEX,IZ,IDZ,DEZ,DFS)
        
        IDZ_WS(ITX,1:NORDZ)=IDZ(1:NORDZ); DEZ_WS(ITX,1:NORDZ)=DEZ(1:NORDZ)
        DFS_WS(ITX,1,1:NORDX,1:4)=DFS(1,1:NORDX,1:4)
        DFS_WS(ITX,2,1:NORDZ,1:4)=DFS(2,1:NORDZ,1:4)
        
    endif
    
    !-Absorbing distance
    NAB1=0;NAB2=0
    IPL=I;KPL=K;IPL1=I1;KPL1=K1
    RX=0.d0;RZ=0.d0;R=0.d0  

    if(MAB>0)then
    If(I<=NEXTD)then
    IPL=NEXTD+1;IPL1=1; NAB1=1
    IDPX=(IPL-1)*(NORDX-1)+IPL1
    RX=sqrt((XP(IDPX)-XP(ITX))**2)
    else if(I>=NX-NEXTD) then
    IPL=NX-NEXTD;IPL1=1;  NAB1=1
    IDPX=(IPL-1)*(NORDX-1)+IPL1
    RX=sqrt((XP(IDPX)-XP(ITX))**2)
    end if
    
    If(K<=NEXTD)then
    KPL=NEXTD+1;KPL1=1; NAB2=1
    IDPZ=((IPL-1)*(NORDX-1)+(IPL1-1))*NNZ+(KPL-1)*(NORDZ-1)+KPL1
    RZ=sqrt((ZP(1,IDPZ)-ZP(1,ID))**2)
    else if(NFS==0.and.K>=NZ-NEXTD)then
    KPL=NZ-NEXTD;KPL1=1;  NAB2=1
    IDPZ=((IPL-1)*(NORDX-1)+(IPL1-1))*NNZ+(KPL-1)*(NORDZ-1)+KPL1
    RZ=sqrt((ZP(1,IDPZ)-ZP(1,ID))**2)
    end if
    endif

    !---------------------------
    RX0(ID)=RX; RZ0(ID)=RZ
    R=SQRT(RX**2+RZ**2)
    if(R>RT)R=RT
    !test(1,ID)=R
    PD(1:22)=PM(1:22,ID)
    test(1,ID)=PD(1)
    test(2,ID)=PD(2)
    EV=0.0;ETA=0.0
    PD1=PD(1)
   
    !-SRM Elastic modulic adjust
    if(MAB==2)then
    if(NAB1==1.or.NAB2==1)then
    !---Apply moduli and ETA*V for RSM
        alpha=Amax*(R/RT)**2.0
    
    PD(2:22)=PD(2:22)*EXP(-alpha)!SRM
    if(I_VISCO==1)tau_TTI(ID,1:NSLS,1:21)=tau_TTI(ID,1:NSLS,1:21)*EXP(-alpha)
    ETA_G(ID)=EMAX*(exp((R/(RT))**2.0)-1.d0)!/PD(1)  !P eta
    
    !----------------------------
  !==============================
    endif
    endif
    
    PM(1:22,ID)=PD(1:22)

!-displaying conversion relax -> unrelax for viscoelastic
if(I_VISCO==1)then

if(Nani(ID)==1)then
    PD(2)=PM(2,ID)/(1+sum(Yhat_L(ID,:,1)))    
else if(Nani(ID)==2)then    
    PD(2)=PM(3,ID)/(1+sum(Yhat_L(ID,:,1)))+2*PM(22,ID)/(1+sum(Yhat_L(ID,:,2)))
    PD(22)=PM(22,ID)/(1+sum(Yhat_L(ID,:,2)))
else if(Nani(ID)==3)then    !VTI2
    PD(2)=PM(2,ID)/(1+sum(Yhat_L(ID,:,1)))
    PD(20)=PM(20,ID)/(1+sum(Yhat_L(ID,:,2)))
    PD(22)=PM(22,ID)/(1+sum(Yhat_L(ID,:,3)))
    
    PD(3)=PD(2)-2*PD(22)
    PD(4)=PD(2)-2*PD(20)
end if 
if(id==IDG(1))then
    if(MNO==3)then
    if(MODD==0)print *, 'unrelaxed P:',sqrt(PM(2,ID)/PM(1,ID)) 
    if(MODD==0)print *, 'relaxed P:',sqrt(PD(2)/PD(1))
   
    else if(MNO==1)then
    if(Nani(ID)==2)then
        if(MODD==0)print *, 'unrelaxed P:',sqrt(PM(2,ID)/PM(1,ID)) 
        if(MODD==0)print *, 'relaxed P:',sqrt((PD(2))/PD(1))
    
        if(MODD==0)print *, 'unrelaxed S:',sqrt(PM(22,ID)/PM(1,ID)) 
        if(MODD==0)print *, 'relaxed S:',sqrt(PD(22)/PD(1))
    endif
    endif
endif
endif

PD(1:22)=PM(1:22,ID)
    !-Assign Matrix in diff medium 
    IX=0
    do ii=1,NORDX
        if(IDX(ii)==ID)IX=ii
    enddo
    if(IX==0)then
        
        write(*,*)'C_DF: Before Chebyshev IX not found'
        stop
    endif
    do JJ=1,NORDZ
        if(IDZ(JJ)==ID)IZ=JJ
    enddo
    if(IZ==0)then
        write(*,*)'C_DF: Before Chebyshev IZ not found'
        stop
    endif
    
    !-Solid & Acoustic matrix
if(mno==1.or.mno==2)then
        call ELASTIC_A(ID,NPT,NORDX,IX,IZ,WKy,DFS,PD,ST_AG,V_AG)
        if(I_VISCO==1)then
        call Memory_Matrix_VTI(ID,IREC,NPT,Nani(ID),NSLS,PD,tau_sig,tau_ep,tau_TTI,AT_Memory)
        endif
endif

if (mno==3.or.mno==2)then
    !if(mno==2)write(91,*)PD(1:2)
        call Acoustic_A(ID,NPT,NORDX,IX,IZ,WKy,DFS,PD,P_AG)
endif

if(I_VISCO==1.and.IREC==1)then
        call TSRC_Matrix(ID,IREC,NPT,NSLS,dt,PD,tau_sig,tau_ep,TSRC_ABC,1)        
endif
!-------------------------------
    17  continue
    16  continue
15  continue
14  continue

if(MODD==0)write(*,*)'C_DF: Before T-step ends'
!=====================

!------
!-Init
P=0.d0; WV=0.d0; P_N=0.d0;WV_N=0.d0
V=0.d0; ST1=0.d0; ST2=0.d0; V_N=0.d0; ST1_N=0.d0; ST2_N=0.d0;
!======================================================================
!-T-stepping
if(MODD==0)write(*,*)'C_DF: T-step starts'
do NT=1,NTD
  T=dfloat(NT-1)*dt
  SE1=(PI*f0)**2;  
  SR1=(1.d0-2.d0*SE1*(T-t0)**2)*exp(-SE1*(T-t0)**2) 
  SR2=-4*SE1*(t-t0)*exp(-se1*(t-t0)**2)+SR1*(-2*se1*(t-t0))
  V_N=0.0;ST1_N=0.0;ST2_N=0.0;P_N=0.0;WV_N=0.0
 
DO 20 I=1,NX-1
DO 21 K=1,NZ-1
  NO=(I-1)*(NORDX-1)*NNZ+(K-1)*(NORDZ-1)+1
      !---
  NORD1=NORDX-1; NORD3=NORDZ-1
  IF(I.EQ.NX-1)NORD1=NORDX
  IF(K.EQ.NZ-1)NORD3=NORDZ
 
!=====================================
!----------------------
!3)---- Chebyshev points in (I,K) -----
!----------------------
  DO 22 I1=1,NORD1
  DO 23 K1=1,NORD3
  !-------------------
  !-Skip the b.d to set 0
    ITX=(I-1)*(NORDX-1)+I1;ITZ=(K-1)*(NORDZ-1)+K1
    !-------------------------
    ID=NO+(I1-1)*NNZ+(K1-1)    !global index
    IT=(I1-1)*NORDZ+K1         !local  index
    
    !-Identify FS
    
    call Assign_Y(ID,MNOG,NPT,NORDX,IDXG,DEXG,NORDZ,IDZG,DEZG,DFSG,wky,WV,P,V,ST1,ST2,DX_Y,DY_Y,DZ_Y)
    
    if(MNOG(ID)==2)call Assign_Y_WS(ID,NPT,NNX,ITX,NORDX,IDXG,DEXG,NORDZ,&
        IDZ_WS,DEZ_WS,DFS_WS,wky,WV,P,V,ST1,ST2,DX_Y_WS,DY_Y_WS,DZ_Y_WS)
    !---------------
    !-PML part
    if(MAB==1)then
    if(ABS(RX0(ID))>1.d-6)then
        call PML(ID,NPT,dt,RT,RX0,f0,VMAX,DMAX,DX_Y,DX_Y_G,PMLX_G)
        DX_Y(1:13)=DX_Y(1:13)+PMLX_G(1:13,ID)
    endif
    if(ABS(RZ0(ID))>1.d-6)then
        call PML(ID,NPT,dt,RT,RZ0,f0,VMAX,DMAX,DZ_Y,DZ_Y_G,PMLZ_G)
        DZ_Y(1:13)=DZ_Y(1:13)+PMLZ_G(1:13,ID)
    endif
    endif
    
    
    if(MNOG(ID)==1.or.MNOG(ID)==2)then
    if(I_VISCO==1)then
    !-Memory var calculation
    if(irec==0)then
        !AT_Memory_L(1:NSLS,1:6,1:6)=AT_MEMORY(ID,1:NSLS,1:6,1:6)
       
        call Memory_var_T3(ID,IREC,NPT,NSLS,dt,DX_Y,DY_Y,DZ_Y,NORDX,IDXG,NORDZ,&
            IDZG,DFSG,tau_sig,tau_ep,AT_memory,V_kl,Var_Memory,Var_Memory_T)
        
    else
        !1st order
        call TSRC_1(ID,IREC,NPT,NSLS,dt,DX_Y,DY_Y,DZ_Y,NORDX,IDXG,NORDZ,&
            IDZG,DFSG,tau_sig,tau_ep,TSRC_ABC,AT_memory,V_kl,Var_Memory,Var_Memory_T)
        
        !2nd order
        !call TSRC_2(ID,IREC,NPT,NSLS,dt,DX_Y,DY_Y,DZ_Y,NORDX,IDXG,NORDZ,&
        !    IDZG,DFSG,tau_sig,tau_ep,TSRC_ABC,AT_memory,V_kl,Var_Memory,Var_Memory_T2)
    endif
    endif
        
        !-Elastic calculation
    call ELASTIC_T2(ID,NPT,V,ST1,ST2,DX_Y,DY_Y,DZ_Y,V_AG,ST_AG,V_T_G,ST_T_G)
       if(MAB==2)V_T_G(1:3,1,ID)=V_T_G(1:3,1,ID)-ETA_G(ID)*V(1:3,ID)
        
    if(I_VISCO==1)then
        
        if(irec==1)then
            do iv=1,6
            ST_T_G(iv,1,ID)=ST_T_G(iv,1,ID)+sum(Var_MEMORY(ID,1:NSLS,iv))
            enddo
        else
            do iv=1,6
            ST_T_G(iv,1,ID)=ST_T_G(iv,1,ID)+(sum(Var_MEMORY(ID,1:NSLS,iv))&
                                +sum(Var_MEMORY_T(ID,1:NSLS,iv)))/2.d0
            enddo
            Var_MEMORY_T(ID,1:NSLS,1:7)=Var_MEMORY(ID,1:NSLS,1:7)
        endif
        
        endif
      
      if(NFS==1.and.ITZ==KFS)then
       call FS_T(ID,NPT,ZP,NNX,ITX,NORDX,IDXG,DEXG,NORDZ,IDZG,DEZG,DFSG,&
                    DX_Y,DY_Y,DZ_Y,PM,V,ST1,ST2,V_T_G,ST_T_G)
       endif
    endif
    
    !-------------------------------
    if(MNOG(ID)==3)then
    if(I_VISCO==1)then
        !-Memory var calculation
    if(irec==0)then
        call Memory_var_W_T2(ID,IREC,NPT,NSLS,wky,dt,DX_Y,DY_Y,DZ_Y,NORDX,&
        IDXG,NORDZ,IDZG,DFSG,PM,tau_sig,tau_ep,V_kl,Var_Memory,Var_Memory_T)
    else
        call Memory_var_TSRC_T(ID,IREC,NPT,NSLS,wky,dt,DX_Y,DY_Y,DZ_Y,NORDX,IDXG,NORDZ,&
        IDZG,DFSG,PM,tau_sig,tau_ep,TSRC_ABC,V_kl,Var_Memory)
    endif
    endif
        
        call Acoustic_T2(ID,NPT,P,WV,DX_Y,DY_Y,DZ_Y,P_AG,P_T_G,WV_T_G)
        
        if(MAB==2)WV_T_G(1:3,1,ID)=WV_T_G(1:3,1,ID)-ETA_G(ID)*WV(1:3,ID)
        
    if(I_VISCO==1)then
        if(irec==1)then
            P_T_G(1,ID)=P_T_G(1,ID)-sum(Var_MEMORY(ID,1:NSLS,1))
        else
            P_T_G(1,ID)=P_T_G(1,ID)-(sum(Var_MEMORY(ID,1:NSLS,1))+sum(Var_MEMORY(ID,1:NSLS,2)))/2.d0
            Var_MEMORY(ID,1:NSLS,1)=Var_MEMORY(ID,1:NSLS,2)
        endif
    endif
        
    endif
    
    !-------------------------------
if(MNOG(ID)==2)then
    if(I_VISCO==1)then
        !-Memory var calculation
    if(irec==0)then
       call Memory_var_WS_T(ID,ITX,IREC,NNX,NPT,NSLS,wky,dt,ZP,DX_Y,DY_Y,DZ_Y,&
       DX_Y_WS,DY_Y_WS,DZ_Y_WS,NORDX,IDXG,NORDZ,IDZG,DFSG,PM,tau_sig,tau_ep,V_kl_WS,Var_Memory_WS)
        
    else
     call Memory_var_TSRC_WS_T(ID,ITX,IREC,NNX,NPT,NSLS,wky,dt,ZP,DX_Y,DY_Y,DZ_Y,&
        DX_Y_WS,DY_Y_WS,DZ_Y_WS,NORDX,IDXG,NORDZ,IDZG,DFSG,PM,tau_sig,tau_ep,TSRC_ABC,V_kl_WS,Var_Memory_WS)
        
    endif
    endif
    
    call Acoustic_T2(ID,NPT,P,WV,DX_Y_WS,DY_Y_WS,DZ_Y_WS,P_AG,P_T_G,WV_T_G)
    if(MAB==2)WV_T_G(1:3,1,ID)=WV_T_G(1:3,1,ID)-ETA_G(ID)*WV(1:3,ID)
    P_T_G(1,ID)=0.0    
    
    if(I_VISCO==1)then
    if(irec==1)then
        P_T_G(1,ID)=-sum(Var_MEMORY_WS(ITX,1:NSLS,1))
    else
        P_T_G(1,ID)=-(sum(Var_MEMORY_WS(ITX,1:NSLS,1))+sum(Var_MEMORY_WS(ITX,1:NSLS,2)))/2.d0
            Var_MEMORY_WS(ITX,1:NSLS,1)=Var_MEMORY_WS(ITX,1:NSLS,2)
    endif
    endif
    
    call WS_T(ID,NPT,ZP,NNX,ITX,NORDX,IDXG,DEXG,NORDZ,IDZG,DEZG,DFSG,&
        IDZ_WS,DEZ_WS,DFS_WS,DX_Y,DY_Y,DZ_Y,DX_Y_WS,DY_Y_WS,DZ_Y_WS,&
        PM,WV,P,V,ST1,ST2,WV_T_G,P_T_G,V_T_G,ST_T_G)
    
endif
!======================================
    !-Add source
    FORCE1=0.D0
    rr=sqrt((xp(itx)-xs)**2+(zp(1,ID)-zs)**2)
    if(rr<=RS)then
      AA=3*pi/(4*(pi*pi-6)*(RS)**3)
      !write(91,*)xp(itx),zp(1,ID),xs,zs,AA*(1+dcos(pi*rr/RS))
      if(MNOG(ID)==3)then
        Force1(1)=SR2*AA*(1+dcos(pi*rr/RS)) !Water
        P_T_G(1,ID)=P_T_G(1,ID)+FORCE1(1)
      endif
      if(MNOG(ID)==1)then
        Force1(NGV)=SR1*AA*(1+dcos(pi*rr/RS)) !Solid
        V_T_G(NGV,1,ID)=V_T_G(NGV,1,ID)+FORCE1(NGV)!*(dcos(45.0*pi/180.0))
      endif
    endif
    
    !-Time stepping, Adams-Bathforth 4th order
    if(MNOG(ID)==1.or.MNOG(ID)==2)then
    
    V_N(1:3,ID)=V(1:3,ID)+dt/12.0*(23.0*V_T_G(1:3,1,ID)-16.0*V_T_G(1:3,2,ID)&
                              +5.0*V_T_G(1:3,3,ID))
    
    ST1_N(1:3,ID)=ST1(1:3,ID)+dt/12.0*(23.0*ST_T_G(1:3,1,ID)-16.0*ST_T_G(1:3,2,ID)&
                              +5.0*ST_T_G(1:3,3,ID))
    ST2_N(1:3,ID)=ST2(1:3,ID)+dt/12.0*(23.0*ST_T_G(4:6,1,ID)-16.0*ST_T_G(4:6,2,ID)&
                              +5.0*ST_T_G(4:6,3,ID))
    endif
    
    if(MNOG(ID)==3.or.MNOG(ID)==2)then
    P_N(ID)=P(ID)+dt/12.0*(23.0*P_T_G(1,ID)-16.0*P_T_G(2,ID)&
                              +5.0*P_T_G(3,ID))
    
    WV_N(1:3,ID)=WV(1:3,ID)+dt/12.0*(23.0*WV_T_G(1:3,1,ID)-16.0*WV_T_G(1:3,2,ID)&
                              +5.0*WV_T_G(1:3,3,ID))
    
    endif
    
    if(NFS==1.and.ITZ==KFS)then
        if(MNOG(ID)==1)V_N(3,ID)=(V_N(3,ID)+V_N(3,ID-1))/2.0
        if(MNOG(ID)==3)P_N(ID)=0.0
    endif
    
    !-Edge of grid = 0
    if(ITX==1.or.ITX==NNX)then
        ST1_N(1:3,ID)=0.0
        ST2_N(1:3,ID)=0.0
        V_N(1:3,ID)=0.0
        WV_N(1:3,ID)=0.0
        P_N(ID)=0.0
    endif
    
    if(ITZ==1.or.(NFS==0.and.ITZ==NNZ))then
        ST1_N(1:3,ID)=0.d0
        ST2_N(1:3,ID)=0.d0
        V_N(1:3,ID)=0.0
        WV_N(1:3,ID)=0.0
        P_N(ID)=0.0
    endif
!-------------  
23 continue 
22 continue 
21 continue 
20 continue

!==================
!--Gaussian Filter on W-S & FS-S
DO I=1,NX-1
DO K=1,NZ-1
  NO=(I-1)*(NORDX-1)*NNZ+(K-1)*(NORDZ-1)+1
      !---
  NORD1=NORDX-1; NORD3=NORDZ-1
  IF(I.EQ.NX-1)NORD1=NORDX
  IF(K.EQ.NZ-1)NORD3=NORDZ

  DO I1=1,NORD1
  DO K1=1,NORD3
  !-------------------
  !-Skip the b.d to set 0
    ITX=(I-1)*(NORDX-1)+I1;ITZ=(K-1)*(NORDZ-1)+K1
    !-------------------------
    ID=NO+(I1-1)*NNZ+(K1-1)    !global index
    IT=(I1-1)*NORDZ+K1         !local  index
   
    
    if(ITX>1.and.ITX<NNX)then
    if(MNOG(ID)==2.or.(MNOG(ID)==1.and.NFS==1.and.ITZ==KFS))then
    !if(I1==1.or.I1==3)then
     
    call ITYPE(ITX,NORDX,NNX,ISTART,IEND,KIX)
    
    !WK(:)=4.2/(sqrt(pi)*DX)*exp(-(4.2/DX*(XP(ITX)-XP(ISTART:IEND)))**2)
    WK(:)=5.2/(sqrt(pi)*DX)*exp(-(5.2/DX*(XP(ITX)-XP(ISTART:IEND)))**2)
    WK(1:NORDX)=WK(1:NORDX)/sum(WK)
    
    !write(91,*)I,I1,ITX,WK(1:NORDX)
    
    GF(1:13,ITX)=0.0
    GF(1,ITX)=dot_product(WK(1:NORDX),WV_N(1,IDXG(ID,1:NORDX)))
    GF(2,ITX)=dot_product(WK(1:NORDX),WV_N(2,IDXG(ID,1:NORDX)))
    GF(3,ITX)=dot_product(WK(1:NORDX),WV_N(3,IDXG(ID,1:NORDX))) !WV_3
    !GF(4,ITX)=dot_product(WK(1:NORDX),P_N(IDXG(ID,1:NORDX)))
    GF(5,ITX)=dot_product(WK(1:NORDX),V_N(1,IDXG(ID,1:NORDX)))
    GF(6,ITX)=dot_product(WK(1:NORDX),V_N(2,IDXG(ID,1:NORDX)))
    GF(7,ITX)=dot_product(WK(1:NORDX),V_N(3,IDXG(ID,1:NORDX)))
    endif
    endif
    enddo
    enddo
enddo
enddo

!-Gaussian filter update
DO I=1,NX-1
DO K=1,NZ-1
  NO=(I-1)*(NORDX-1)*NNZ+(K-1)*(NORDZ-1)+1
      !---
  NORD1=NORDX-1; NORD3=NORDZ-1
  IF(I.EQ.NX-1)NORD1=NORDX
  IF(K.EQ.NZ-1)NORD3=NORDZ
 
!----------------------
  DO I1=1,NORD1
  DO K1=1,NORD3
  !-------------------
  !-Skip the b.d to set 0
    ITX=(I-1)*(NORDX-1)+I1;ITZ=(K-1)*(NORDZ-1)+K1
    !-------------------------
    ID=NO+(I1-1)*NNZ+(K1-1)    !global index
    IT=(I1-1)*NORDZ+K1         !local  index
    
    if(ITX>1.and.ITX<NNX)then
    if(MNOG(ID)==2)then
        WV_N(1,ID)=GF(1,ITX)
        V_N(3,ID)=GF(7,ITX)
    endif
    
    if(MNOG(ID)==1.and.NFS==1.and.ITZ==KFS)then
    !    V_N(1,ID)=GF(5,ITX)
        V_N(3,ID)=GF(7,ITX)
        
    endif
    endif
    
    enddo
    enddo
enddo
enddo

!==========================
!-Each time-step renew parameters
    V_T_G(1:3,2:5,1:NPT)=V_T_G(1:3,1:4,1:NPT)
    V_T_G(1:3,1,1:NPT)=0.0
    
    ST_T_G(1:6,2:5,1:NPT)=ST_T_G(1:6,1:4,1:NPT)
    ST_T_G(1:6,1,1:NPT)=0.0
    
    P_T_G(2:5,1:NPT)=P_T_G(1:4,1:NPT)
    P_T_G(1,1:NPT)=0.0
    
    WV_T_G(1:3,2:5,1:NPT)=WV_T_G(1:3,1:4,1:NPT)
    WV_T_G(1:3,1,1:NPT)=0.0
    
V(1:3,1:NPT)=V_N; ST1(1:3,1:NPT)=ST1_N;ST2(1:3,1:NPT)=ST2_N
P(1:NPT)=P_N; WV(1:3,1:NPT)=WV_N;

do ir=1,NSR
    FTG(ir,1:3,NT) = WV(1:3,IDG(ir))
    FTG(ir,4,NT) = P(IDG(ir))
    FTG(ir,5:7,NT) = V(1:3,IDG(ir))
    FTG(ir,8:10,NT) = ST1(1:3,IDG(ir))
    FTG(ir,11:13,NT) = ST2(1:3,IDG(ir))
enddo

!-Snapshot
if(MODD==0)then
  !write(900,*)real(FTG(1:4,NT))
  if(NT>1000)then
  if(mod(NT,1000)==1)then
  
    write(*,*)ifile,'|',T,'/',(NTD-1)*dt
    do ii=1,NNX
      write(100+ifile,*)real(V_N(1,1+(ii-1)*NNZ:ii*NNZ))
      write(300+ifile,*)real(V_N(3,1+(ii-1)*NNZ:ii*NNZ))
      write(400+ifile,*)real(P_N(1+(ii-1)*NNZ:ii*NNZ))
      write(500+ifile,*)real(WV_N(1,1+(ii-1)*NNZ:ii*NNZ))
      write(700+ifile,*)real(WV_N(3,1+(ii-1)*NNZ:ii*NNZ))
    enddo
    
    close(100+ifile)
    close(200+ifile)
    close(300+ifile)
    close(400+ifile)
    close(500+ifile)
    close(700+ifile)
    
    ifile=ifile+1
  endif
  endif
endif


enddo!T-step
!================
!-Output test for any check
do ii=1,NNX
    !write(91,*)MNOG(1+(ii-1)*NNZ:ii*NNZ)
    !write(91,*)TEST(1,1+(ii-1)*NNZ:ii*NNZ)
    !write(92,*)TEST(2,1+(ii-1)*NNZ:ii*NNZ)
    !write(93,*)TEST(3,1+(ii-1)*NNZ:ii*NNZ)
    !write(93,*)TEST(2,1+(ii-1)*NNZ:ii*NNZ)
enddo


RETURN
END subroutine

!-----------------------------------
subroutine output_SNAPSHOT(FKY,NNX,NNZ,snapshot)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
character(len=10) :: char_i
dimension snapshot(NNX*NNZ)


write(char_i, '(F6.3)') FKY
open(87, file='snapshot_real_'//trim(adjustl(char_i))//'.out')
do ii=1,NNX
      write(87,*)real(snapshot(1+(ii-1)*NNZ:ii*NNZ))
enddo
close(87)

!open(87, file='snapshot_imag_'//trim(adjustl(char_i))//'.out')
!do ii=1,NNX
!      write(87,*)imag(snapshot(1+(ii-1)*NNZ:ii*NNZ))
!enddo
!close(87)

end subroutine