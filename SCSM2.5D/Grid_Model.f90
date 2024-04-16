!=======================================================================
!-------------
Subroutine EXTERNAL_GRID_MODEL2(NNX,NNZ,INDP,ITTI,PT_TEMP,PM)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION PT_TEMP(NNX*NNZ,INDP+1),PM(22,NNX*NNZ)

CHARACTER LAB*80
Doubleprecision,DIMENSION(:,:),allocatable :: rho, ZTA0
Doubleprecision,DIMENSION(:,:,:),allocatable:: V0
Doubleprecision,DIMENSION(:,:,:):: V1(INDP,NNX,NNZ), PM0(22,NNX*NNZ)
DOUBLEPRECISION,DIMENSION(:)::PT11(6)
DOUBLEPRECISION,DIMENSION(:)::PT2(21)
!--------------------------------------------------------
OPEN(1,FILE='rho.inp',STATUS='UNKNOWN')
OPEN(11,FILE='c11.inp',STATUS='UNKNOWN')
OPEN(12,FILE='c13.inp',STATUS='UNKNOWN')
OPEN(13,FILE='c33.inp',STATUS='UNKNOWN')
OPEN(14,FILE='c44.inp',STATUS='UNKNOWN')
OPEN(15,FILE='c66.inp',STATUS='UNKNOWN')

if(ITTI==1)OPEN(16,FILE='theta.inp',STATUS='UNKNOWN')

print *,'Model input start'
PI=4.d0*atan(1.d0)
PM(:,:)=0.d0
PT_TEMP(:,:)=0.0
!-DENSITY

!Check the grid corresponding
    read(1,*)LAB
    read(1,*)NX,NZ
    
if(NX/=NNX)then
print *,'rho GRID X /= EXTERNAL X'
stop
endif

if(NZ/=NNZ)then
print *,'rho GRID Z /= EXTERNAL Z'
stop
endif


allocate(rho(NNX,NNZ))

!read Density
RHO(:,:)=0.d0;
do j=1,NNX
    read(1,*)rho(j,1:NNZ)
enddo

close(1)

!---------------------------------
!-P_velocity
allocate(V0(INDP,NNX,NNZ))
V0(:,:,:)=0.d0
do ii=1,INDP
!print *, INDP,ii
read(ii+10,*)LAB
read(ii+10,*)NX,NZ

!Check the grid corresponding
if(NX/=NNX)then
print *,'v GRID X /= EXTERNAL X'
stop
endif

if(NZ/=NNZ)then
print *,'v GRID Z /= EXTERNAL Z'
stop
endif

!read Density
do j=1,NNX
    read(ii+10,*)V0(ii,j,1:NNZ)
enddo

close(ii+10)
enddo

V1(:,:,:)=0.d0

allocate(zta0(NNX,NNZ))
    ZTA0(:,:)=0.d0; 
if(ITTI==1)then
    read(16,*)LAB
    read(16,*)NX,NZ
!read Density
do j=1,NNX
    read(16,*)ZTA0(j,1:NNZ)
    !write(91,*)rho0(j,1:NX0)
enddo
close(16)
endif

print *, 'External grid input end'


!------------
DO I=1,NNX
DO K=1,NNZ

    !-------------------------
    ID=(I-1)*NNZ+K    !global index
    
    PT11(:)=0.d0; PT2(:)=0.d0
    
    PM(1,ID)=RHO(I,K)
    PT11(1:INDP)=V0(1:INDP,I,K)
   
    call TTI(INDP,PT11,ITTI,tan(pi*ZTA0(I,K)/180.0),PT2)
    PM(2:22,ID)=PT2(1:21)
    
enddo
enddo

do ii=1,NNX
  !write(91,*)PM(1,1+(ii-1)*NNZ:ii*NNZ)
  !write(92,*)PM(2,1+(ii-1)*NNZ:ii*NNZ)
  !write(93,*)PM(4,1+(ii-1)*NNZ:ii*NNZ)
  !write(94,*)PM(13,1+(ii-1)*NNZ:ii*NNZ)
  !write(95,*)PM(17,1+(ii-1)*NNZ:ii*NNZ)
  !write(96,*)PM(22,1+(ii-1)*NNZ:ii*NNZ)
enddo

do ii=1,NZ
  !write(94,*)RHO1(ii,1:NX)
  !write(95,*)VP1(ii,1:NX)
  !write(96,*)VS1(ii,1:NX)
enddo

end subroutine


!-------------
Subroutine EXTERNAL_GRID_Relax_time(NX,NZ,DX,NORDX,NORDZ,NNX,NNZ,NEXTD,NFS,ITTI,PT_TEMP,nani,nsls,tau_ep,tau_sig,tau_TTI)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION PT_TEMP(nnx*nnz,NANI+1)
DIMENSION tau_ep(nnx*nnz,nsls,21),tau_sig(nnx*nnz,nsls,21)
DIMENSION TAU_TTI(nnx*nnz,NSLS,21)

CHARACTER LAB*80
Doubleprecision,DIMENSION(:), allocatable :: RELAX_C, PT11, PT2,PT0
Doubleprecision,DIMENSION(:,:,:,:),allocatable :: tau_ep0,tau_sig0
Doubleprecision,DIMENSION(:,:,:,:):: tau_ep1(NZ,NX,NSLS,21),tau_sig1(NZ,NX,NSLS,21)
Doubleprecision,DIMENSION(:,:,:):: TAU_TTI0(NNX*NNZ,NSLS,21)

!--------------------------------------------------------
!-open the Marmousi relax-time input
OPEN(21,FILE='tau_ep_11.inp',STATUS='UNKNOWN')
OPEN(22,FILE='tau_ep_13.inp',STATUS='UNKNOWN')
OPEN(23,FILE='tau_ep_33.inp',STATUS='UNKNOWN')
OPEN(24,FILE='tau_ep_44.inp',STATUS='UNKNOWN')
OPEN(25,FILE='tau_ep_66.inp',STATUS='UNKNOWN')

OPEN(31,FILE='tau_sig_11.inp',STATUS='UNKNOWN')
OPEN(32,FILE='tau_sig_13.inp',STATUS='UNKNOWN')
OPEN(33,FILE='tau_sig_33.inp',STATUS='UNKNOWN')
OPEN(34,FILE='tau_sig_44.inp',STATUS='UNKNOWN')
OPEN(35,FILE='tau_sig_66.inp',STATUS='UNKNOWN')

tau_ep(:,:,:)=0.d0; tau_sig(:,:,:)=0.d0
PI=4.d0*atan(1.d0); TAU_TTI0(:,:,:)=0.d0

do i=1,Nani
!-P-Tau-ep
    read(i+20,*)LAB
    read(i+20,*)NX0,NZ0,NSLS0

    if(i==1)then
        allocate(tau_ep0(NZ0,NX0,NSLS0,21),tau_sig0(NZ0,NX0,NSLS0,21))
        tau_ep0(1:NZ0,1:NX0,1:NSLS0,1:21)=0.0
        tau_sig0(1:NZ0,1:NX0,1:NSLS0,1:21)=0.0
    endif
    
    !Check the grid corresponding
    if(NX/=NX0+2*NEXTD)then
    print *,'GRID X /= EXTERNAL X'
    stop
    endif

    if(NZ/=NZ0+2*NEXTD .and. NFS==0)then
    !print *, NZ, NZ0, NEXTD
    print *,'rho GRID Z /= EXTERNAL Z'
    stop
    else if(NZ/=NZ0+NEXTD .and. NFS==1)then
    !print *, NZ, NZ0, NEXTD
    print *,'rho GRID Z /= EXTERNAL Z'
    stop
    endif

    if(NSLS/=NSLS0)then
    print *, 'GRID DX /= EXTERNAL DX'
    stop
    endif

    do N=1,NSLS0
    read(i+20,*)LAB
    do j=1,NZ0
    read(i+20,*)tau_ep0(j,1:NX0,N,i)
    !write(91,*)rho0(j,1:NX0)
    enddo
    enddo
    close(i+20)

    !-------------------------
    !-P-Tau-sig
    read(i+30,*)LAB
    read(i+30,*)NX0,NZ0,NSLS0

    !Check the grid corresponding
    if(NX/=NX0+2*NEXTD)then
    print *,'GRID X /= EXTERNAL X'
    stop
    endif

   if(NZ/=NZ0+2*NEXTD .and. NFS==0)then
    print *,'rho GRID Z /= EXTERNAL Z'
    stop
    else if(NZ/=NZ0+NEXTD .and. NFS==1)then
    print *,'rho GRID Z /= EXTERNAL Z'
    stop
    endif

    if(NSLS/=NSLS0)then
    print *, 'GRID DX /= EXTERNAL DX'
    stop
    endif

    do N=1,NSLS0
    read(i+30,*)LAB
    do j=1,NZ0
        read(i+30,*)tau_sig0(j,1:NX0,N,i)
    !write(91,*)rho0(j,1:NX0)
    enddo
    enddo
    close(i+30)
enddo
    

print *, 'External grid input end'
 !stop 
!-------------------------------------
!-Grid Assign

TAU_EP1(1:NZ,1:NX,1:NSLS,1:21)=0.d0

if(NFS==0)then
TAU_EP1(NEXTD+1:NEXTD+NZ0,NEXTD+1:NEXTD+NX0,1:NSLS,1:21)=TAU_EP0
do j=1,NEXTD
do i=1,NEXTD
    TAU_EP1(j,i,:,:)=TAU_EP0(1,1,:,:)
enddo
enddo

do j=NEXTD+NZ0,NZ
do i=1,NEXTD
TAU_EP1(j,i,:,:)=TAU_EP0(NZ0,1,:,:)
enddo
enddo

do j=1,NEXTD
do i=NEXTD+NX0,NX
TAU_EP1(j,i,:,:)=TAU_EP0(1,NX0,:,:)
enddo
enddo

do j=NEXTD+NZ0,NZ
do i=NEXTD+NX0,NX
TAU_EP1(j,i,:,:)=TAU_EP0(NZ0,NX0,:,:)
enddo
enddo

else
TAU_EP1(1:NZ0,NEXTD+1:NEXTD+NX0,1:NSLS,1:21)=TAU_EP0
do j=1,NEXTD
do i=1,NEXTD
    TAU_EP1(1,i,:,:)=TAU_EP0(1,1,:,:)
enddo
enddo

do j=NZ0,NZ
do i=1,NEXTD
TAU_EP1(j,i,:,:)=TAU_EP0(NZ0,1,:,:)
enddo
enddo

do j=1,NEXTD
do i=NEXTD+NX0,NX
TAU_EP1(1,i,:,:)=TAU_EP0(1,NX0,:,:)
enddo
enddo

do j=NZ0,NZ
do i=NEXTD+NX0,NX
TAU_EP1(j,i,:,:)=TAU_EP0(NZ0,NX0,:,:)
enddo
enddo
endif


if(NFS==0)then
do i=1,NEXTD
    TAU_EP1(NEXTD+1:NEXTD+NZ0,i,:,:)=TAU_EP0(:,1,:,:)
    TAU_EP1(NEXTD+1:NEXTD+NZ0,NEXTD+NX0+i,:,:)=TAU_EP0(:,NX0,:,:)
    TAU_EP1(i,NEXTD+1:NEXTD+NX0,:,:)=TAU_EP0(1,:,:,:)
    TAU_EP1(NEXTD+NZ0+i,NEXTD+1:NEXTD+NX0,:,:)=TAU_EP0(NZ0,:,:,:)
enddo
else
do i=1,NEXTD
    TAU_EP1(1:NZ0,i,:,:)=TAU_EP0(:,1,:,:)
    TAU_EP1(1:NZ0,NEXTD+NX0+i,:,:)=TAU_EP0(:,NX0,:,:)
    !TAU_EP1(i,NEXTD+1:NEXTD+NX0,:,:)=TAU_EP0(1,:,:,:)
    TAU_EP1(NZ0+i,NEXTD+1:NEXTD+NX0,:,:)=TAU_EP0(NZ0,:,:,:)
enddo
endif

print *, 'TAU_EP1'


!----------------------------
TAU_SIG1(1:NZ,1:NX,1:NSLS,1:21)=0.d0

if(NFS==0)then
TAU_SIG1(NEXTD+1:NEXTD+NZ0,NEXTD+1:NEXTD+NX0,1:NSLS,1:21)=TAU_SIG0
do j=1,NEXTD
do i=1,NEXTD
    TAU_SIG1(j,i,:,:)=TAU_SIG0(1,1,:,:)
enddo
enddo

do j=NEXTD+NZ0,NZ
do i=1,NEXTD
TAU_SIG1(j,i,:,:)=TAU_SIG0(NZ0,1,:,:)
enddo
enddo

do j=1,NEXTD
do i=NEXTD+NX0,NX
TAU_SIG1(j,i,:,:)=TAU_SIG0(1,NX0,:,:)
enddo
enddo

do j=NEXTD+NZ0,NZ
do i=NEXTD+NX0,NX
TAU_SIG1(j,i,:,:)=TAU_SIG0(NZ0,NX0,:,:)
enddo
enddo

else
TAU_SIG1(1:NZ0,NEXTD+1:NEXTD+NX0,1:NSLS,1:21)=TAU_SIG0
do j=1,NEXTD
do i=1,NEXTD
    TAU_SIG1(1,i,:,:)=TAU_SIG0(1,1,:,:)
enddo
enddo

do j=NZ0,NZ
do i=1,NEXTD
TAU_SIG1(j,i,:,:)=TAU_SIG0(NZ0,1,:,:)
enddo
enddo

do j=1,NEXTD
do i=NEXTD+NX0,NX
TAU_SIG1(1,i,:,:)=TAU_SIG0(1,NX0,:,:)
enddo
enddo

do j=NZ0,NZ
do i=NEXTD+NX0,NX
TAU_SIG1(j,i,:,:)=TAU_SIG0(NZ0,NX0,:,:)
enddo
enddo
endif

if(NFS==0)then
do i=1,NEXTD
    TAU_SIG1(NEXTD+1:NEXTD+NZ0,i,:,:)=TAU_SIG0(:,1,:,:)
    TAU_SIG1(NEXTD+1:NEXTD+NZ0,NEXTD+NX0+i,:,:)=TAU_SIG0(:,NX0,:,:)
    TAU_SIG1(i,NEXTD+1:NEXTD+NX0,:,:)=TAU_SIG0(1,:,:,:)
    TAU_SIG1(NEXTD+NZ0+i,NEXTD+1:NEXTD+NX0,:,:)=TAU_SIG0(NZ0,:,:,:)
enddo
else
do i=1,NEXTD
    TAU_SIG1(1:NZ0,i,:,:)=TAU_SIG0(:,1,:,:)
    TAU_SIG1(1:NZ0,NEXTD+NX0+i,:,:)=TAU_SIG0(:,NX0,:,:)
    !TAU_EP1(i,NEXTD+1:NEXTD+NX0,:,:)=TAU_EP0(1,:,:,:)
    TAU_SIG1(NZ0+i,NEXTD+1:NEXTD+NX0,:,:)=TAU_SIG0(NZ0,:,:,:)
enddo
endif

print *, 'TAU_SIG1'

!-------------------------
DO I=1,NX-1
DO K=1,NZ-1
NO=(I-1)*(NORDX-1)*NNZ+(K-1)*(NORDZ-1)+1
NORD1=NORDX-1; NORD3=NORDZ-1
IF(I.EQ.NX-1)NORD1=NORDX
IF(K.EQ.NZ-1)NORD3=NORDZ

!--------------------------
    DO I1=1,NORD1
    DO K1=1,NORD3
    ITX=(I-1)*(NORDX-1)+I1
    ITZ=(K-1)*(NORDZ-1)+K1
    !-------------------------
    ID=NO+(I1-1)*NNZ+(K1-1)    !global index
    IT=(I1-1)*NORDZ+K1         !local  index
    
    
    IIN=1; IP=I-NEXTD+1; KP=K-NEXTD+1
    If(I<=NEXTD)then
    IP=1;IIN=0
    KP=K-NEXTD+1
    else if(I>=NX-NEXTD) then
    IP=NX-2*NEXTD;IIN=0
    KP=K-NEXTD+1
    end if
    
    If(K<=NEXTD)then
    KP=1; IIN=0
    !IP=I-NEXTD+1
    else if(K>=NZ-NEXTD)then
    KP=NZ-2*NEXTD; IIN=0
    !IP=I-NEXTD+1
    end if
    
    if(IIN==1)then
        IP=I-NEXTD+1; KP=K-NEXTD+1
        !print *, IP,KP
        TAU_EP(ID,1:NSLS,1:21)=TAU_EP1(K,I,1:NSLS,1:21)
        TAU_SIG(ID,1:NSLS,1:21)=TAU_SIG1(K,I,1:NSLS,1:21)
        
    else
        !print *, IP,KP
        TAU_EP(ID,1:NSLS,1:21)=TAU_EP1(K,I,1:NSLS,1:21)
        TAU_SIG(ID,1:NSLS,1:21)=TAU_SIG1(K,I,1:NSLS,1:21)
        
    endif
    TAU_EP(ID,1:NSLS,1:21)=TAU_EP1(K,I,1:NSLS,1:21)
    TAU_SIG(ID,1:NSLS,1:21)=TAU_SIG1(K,I,1:NSLS,1:21)
    
    enddo
    enddo
enddo
enddo

print *, '---------'
!stop
!TTI application
allocate(Relax_C(6),PT11(Nani),PT2(21),PT0(Nani))
tau_TTI(1:NNX*NNZ,1:NSLS,1:21)=0.d0
  do i=1,NNX
  do j=1,NNZ
    !------------------------
    ID=(i-1)*NNZ+j
    
    !print *, PT1(mn,2:6)
    PT0(1:NANI)=PT_TEMP(ID,1:NANI)
    ZTA0=PT_TEMP(ID,NANI+1)
    !C_Relax(1:5)=0.0
    
    do k=1,Nani
    PT11(k)=1.d0+sum(TAU_EP(ID,1:NSLS,k)/TAU_SIG(ID,1:NSLS,k)-1.d0)/dble(NSLS)
    Relax_C(k)=PT0(k)/PT11(k)
    enddo
    
    
    do k=1,NSLS
    PT11=0.d0;  PT2=0.d0
    PT11(1:Nani)=(1.d0-TAU_EP(ID,k,1:Nani)/TAU_SIG(ID,k,1:Nani))/TAU_SIG(ID,k,1:Nani)*RELAX_C(1:Nani)
    
    call TTI(Nani,PT11,ITTI,tan(pi*ZTA0/180.0),PT2)

    tau_TTI(ID,k,1:21)=PT2(1:21)
    !if(ID==45313)write(91,*)tau_TTI(ID,k,1:21)
    enddo
    
    !-----------------------
  enddo
  enddo
  
  TAU_TTI0(1:NNX*NNZ,1:NSLS,1:21)=TAU_TTI(1:NNX*NNZ,1:NSLS,1:21)
  TAU_TTI(:,:,:)=0.0
  
  do i=1,NNX
    do j=1,NNZ
        ID0=(I-1)*NNZ+J
        ID=(I-1)*NNZ+NNZ+1-J
        TAU_TTI(ID,1:NSLS,1:21) = TAU_TTI0(ID0,1:NSLS,1:21)
    enddo
  enddo
  
 
end subroutine

!=======================================================================
!-------------
Subroutine EXTERNAL_GRID_MODEL(NX,NZ,DX,NORDX,NORDZ,NNX,NNZ,NEXTD,NFS,INDP,ITTI,PT_TEMP,PM)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION PT_TEMP(NNX*NNZ,INDP+1),PM(22,NNX*NNZ)

CHARACTER LAB*80
Doubleprecision,DIMENSION(:,:),allocatable :: rho0, ZTA0
Doubleprecision,DIMENSION(:,:):: rho1(NZ,NX), ZTA1(NZ,NX)
Doubleprecision,DIMENSION(:,:,:),allocatable:: V0
Doubleprecision,DIMENSION(:,:,:):: V1(INDP,NZ,NX), PM0(22,NNX*NNZ)
DOUBLEPRECISION,DIMENSION(:)::PT11(6)
DOUBLEPRECISION,DIMENSION(:)::PT2(21)
!--------------------------------------------------------
OPEN(1,FILE='rho.inp',STATUS='UNKNOWN')
OPEN(11,FILE='c11.inp',STATUS='UNKNOWN')
OPEN(12,FILE='c13.inp',STATUS='UNKNOWN')
OPEN(13,FILE='c33.inp',STATUS='UNKNOWN')
OPEN(14,FILE='c44.inp',STATUS='UNKNOWN')
OPEN(15,FILE='c66.inp',STATUS='UNKNOWN')

if(ITTI==1)OPEN(16,FILE='theta.inp',STATUS='UNKNOWN')

print *,'BP2007 model input start'
PI=4.d0*atan(1.d0)
PM(:,:)=0.d0
PT_TEMP(:,:)=0.0
!-DENSITY

!Check the grid corresponding
    read(1,*)LAB
    read(1,*)NX0,NZ0,DX0
    
if(NX/=NX0+2*NEXTD)then
print *,'rho GRID X /= EXTERNAL X'
stop
endif

if(NZ/=NZ0+2*NEXTD .and. NFS==0)then
print *,'rho GRID Z /= EXTERNAL Z'
stop
else if(NZ/=NZ0+NEXTD .and. NFS==1)then
print *,'rho GRID Z /= EXTERNAL Z'
stop
endif

if(DX/=DX0)then
print *, 'rho GRID DX /= EXTERNAL DX'
stop
endif

allocate(rho0(NZ0,NX0))

!read Density
RHO0(:,:)=0.d0; RHO1(:,:)=0.d0
do j=1,NZ0
    
    read(1,*)rho0(j,1:NX0)
    !write(91,*)rho0(j,1:NX0)
enddo

close(1)

!---------------------------------
!-P_velocity
allocate(V0(INDP,NZ0,NX0))
V0(:,:,:)=0.d0
do ii=1,INDP
!print *, INDP,ii
read(ii+10,*)LAB
read(ii+10,*)NX0,NZ0,DX0

!Check the grid corresponding
if(NX/=NX0+2*NEXTD)then
print *,'v GRID X /= EXTERNAL X'
stop
endif

if(NZ/=NZ0+2*NEXTD.and. NFS==0)then
print *,'v GRID Z /= EXTERNAL Z'
stop
else if(NZ/=NZ0+NEXTD .and. NFS==1)then
print *,'rho GRID Z /= EXTERNAL Z'
stop
endif

if(DX/=DX0)then
print *, 'v GRID DX /= EXTERNAL DX'
stop
endif

!read Density
do j=1,NZ0
    read(ii+10,*)V0(ii,j,1:NX0)
    !write(92,*)VP0(j,1:NX0)
enddo

close(ii+10)
enddo
V1(:,:,:)=0.d0

allocate(zta0(NZ0,NX0))
ZTA0(:,:)=0.d0; ZTA1(:,:)=0.d0
if(ITTI==1)then
    read(16,*)LAB
    read(16,*)NX0,NZ0,DX0
!read Density
do j=1,NZ0
    read(16,*)ZTA0(j,1:NX0)
    !write(91,*)rho0(j,1:NX0)
enddo
close(16)
endif

print *, 'External grid input end'
!-------------------------------------
!-Grid Assign

if(NFS==0)then
RHO1(NEXTD+1:NEXTD+NZ0,NEXTD+1:NEXTD+NX0)=RHO0
! Four edge
RHO1(1:NEXTD,1:NEXTD)=RHO0(1,1)
RHO1(NEXTD+NZ0:NZ,1:NEXTD)=RHO0(NZ0,1)
RHO1(1:NEXTD,NEXTD+NX0:NX)=RHO0(1,NX0)
RHO1(NEXTD+NZ0:NZ,NEXTD+NX0:NX)=RHO0(NZ0,NX0)
else
RHO1(1:NZ0,NEXTD+1:NEXTD+NX0)=RHO0

RHO1(1,1:NEXTD)=RHO0(1,1)
RHO1(NZ0:NZ,1:NEXTD)=RHO0(NZ0,1)
RHO1(1,NEXTD+NX0:NX)=RHO0(1,NX0)
RHO1(NZ0:NZ,NEXTD+NX0:NX)=RHO0(NZ0,NX0)
endif

if(NFS==0)then
do i=1,NEXTD
    RHO1(NEXTD+1:NEXTD+NZ0,i)=RHO0(:,1)
    RHO1(NEXTD+1:NEXTD+NZ0,NEXTD+NX0+i)=RHO0(:,NX0)
    RHO1(i,NEXTD+1:NEXTD+NX0)=RHO0(1,:)
    RHO1(NEXTD+NZ0+i,NEXTD+1:NEXTD+NX0)=RHO0(NZ0,:)
enddo
else
do i=1,NEXTD
    RHO1(1:NZ0,i)=RHO0(:,1)
    RHO1(1:NZ0,NEXTD+NX0+i)=RHO0(:,NX0)
    !RHO1(i,NEXTD+1:NEXTD+NX0)=RHO0(1,:)
    RHO1(NZ0+i,NEXTD+1:NEXTD+NX0)=RHO0(NZ0,:)
enddo
endif
print *, 'rho1'

do ii=1,INDP
if(NFS==0)then
V1(ii,NEXTD+1:NEXTD+NZ0,NEXTD+1:NEXTD+NX0)=V0(ii,:,:)
V1(ii,1:NEXTD,1:NEXTD)=V0(ii,1,1)
V1(ii,NEXTD+NZ0:NZ,1:NEXTD)=V0(ii,NZ0,1)
V1(ii,1:NEXTD,NEXTD+NX0:NX)=V0(ii,1,NX0)
V1(ii,NEXTD+NZ0:NZ,NEXTD+NX0:NX)=V0(ii,NZ0,NX0)
else
V1(ii,1:NZ0,NEXTD+1:NEXTD+NX0)=V0(ii,:,:)

V1(ii,1,1:NEXTD)=V0(ii,1,1)
V1(ii,NZ0:NZ,1:NEXTD)=V0(ii,NZ0,1)
V1(ii,1,NEXTD+NX0:NX)=V0(ii,1,NX0)
V1(ii,NZ0:NZ,NEXTD+NX0:NX)=V0(ii,NZ0,NX0)
endif

if(NFS==0)then
do i=1,NEXTD
    V1(ii,NEXTD+1:NEXTD+NZ0,i)=V0(ii,:,1)
    V1(ii,NEXTD+1:NEXTD+NZ0,NEXTD+NX0+i)=V0(ii,:,NX0)
    V1(ii,i,NEXTD+1:NEXTD+NX0)=V0(ii,1,:)
    V1(ii,NEXTD+NZ0+i,NEXTD+1:NEXTD+NX0)=V0(ii,NZ0,:)
enddo
else
do i=1,NEXTD
    V1(ii,1:NZ0,i)=V0(ii,:,1)
    V1(ii,1:NZ0,NEXTD+NX0+i)=V0(ii,:,NX0)
    !V1(ii,i,NEXTD+1:NEXTD+NX0)=V0(ii,1,:)
    V1(ii,NZ0+i,NEXTD+1:NEXTD+NX0)=V0(ii,NZ0,:)
enddo
endif
enddo
print *, 'V'

if(ITTI==1)then
if(NFS==0)then
ZTA1(NEXTD+1:NEXTD+NZ0,NEXTD+1:NEXTD+NX0)=ZTA0
ZTA1(1:NEXTD,1:NEXTD)=ZTA0(1,1)
ZTA1(NEXTD+NZ0:NZ,1:NEXTD)=ZTA0(NZ0,1)
ZTA1(1:NEXTD,NEXTD+NX0:NX)=ZTA0(1,NX0)
ZTA1(NEXTD+NZ0:NZ,NEXTD+NX0:NX)=ZTA0(NZ0,NX0)
else
ZTA1(1:NZ0,NEXTD+1:NEXTD+NX0)=ZTA0
ZTA1(1,1:NEXTD)=ZTA0(1,1)
ZTA1(NZ0:NZ,1:NEXTD)=ZTA0(NZ0,1)
ZTA1(1,NEXTD+NX0:NX)=ZTA0(1,NX0)
ZTA1(NZ0:NZ,NEXTD+NX0:NX)=ZTA0(NZ0,NX0)
endif

if(NFS==0)then
do i=1,NEXTD
    ZTA1(NEXTD+1:NEXTD+NZ0,i)=ZTA0(:,1)
    ZTA1(NEXTD+1:NEXTD+NZ0,NEXTD+NX0+i)=ZTA0(:,NX0)
    ZTA1(i,NEXTD+1:NEXTD+NX0)=ZTA0(1,:)
    ZTA1(NEXTD+NZ0+i,NEXTD+1:NEXTD+NX0)=ZTA0(NZ0,:)
enddo
else
do i=1,NEXTD
    ZTA1(1:NZ0,i)=ZTA0(:,1)
    ZTA1(1:NZ0,NEXTD+NX0+i)=ZTA0(:,NX0)
    !ZTA1(i,NEXTD+1:NEXTD+NX0)=ZTA0(1,:)
    ZTA1(NZ0+i,NEXTD+1:NEXTD+NX0)=ZTA0(NZ0,:)
enddo
endif
print *, 'ZTA1'
endif

!------------
DO I=1,NX-1
DO K=1,NZ-1
NO=(I-1)*(NORDX-1)*NNZ+(K-1)*(NORDZ-1)+1
NORD1=NORDX-1; NORD3=NORDZ-1
IF(I.EQ.NX-1)NORD1=NORDX
IF(K.EQ.NZ-1)NORD3=NORDZ

!--------------------------
    DO I1=1,NORD1
    DO K1=1,NORD3
    ITX=(I-1)*(NORDX-1)+I1
    ITZ=(K-1)*(NORDZ-1)+K1
    !-------------------------
    ID=NO+(I1-1)*NNZ+(K1-1)    !global index
    IT=(I1-1)*NORDZ+K1         !local  index
    
    PT11(:)=0.d0; PT2(:)=0.d0
    
    IIN=1; IP=I-NEXTD+1; KP=K-NEXTD+1
    If(I<=NEXTD)then
    IP=1;IIN=0
    KP=K-NEXTD+1
    else if(I>=NX-NEXTD) then
    IP=NX-2*NEXTD;IIN=0
    KP=K-NEXTD+1
    end if
    
    If(K<=NEXTD)then
    KP=1; IIN=0
    !IP=I-NEXTD+1
    else if(K>=NZ-NEXTD)then ! Above the abs layer
    if(NFS==0)KP=NZ-2*NEXTD;
    if(NFS==1)KP=NZ-NEXTD;
    IIN=0
    !IP=I-NEXTD+1
    end if
    
    if(IIN==1)then
        IP=I-NEXTD+1; 
        KP=K-NEXTD+1
       ! print *, K,IP,KP
        PM(1,ID)=RHO1(K,I)
        PT11(1:INDP)=V1(1:INDP,K,I)
    else
        !print *, IP,KP
        PM(1,ID)=RHO1(K,I)
        PT11(1:INDP)=V1(1:INDP,K,I)
        !print *,PM(1,ID),PT11(1)
    endif
    
    PT_TEMP(ID,1:INDP)=PT11(1:INDP)
    PT_TEMP(ID,INDP+1)=ZTA1(K,I)
    call TTI(INDP,PT11,ITTI,tan(pi*ZTA1(K,I)/180.0),PT2)
    PM(2:22,ID)=PT2(1:21)
    
    enddo
    enddo
enddo
enddo

PM0(1:22,:)=PM(1:22,:);
PM(:,:)=0.d0;

do ii=1,NNX
    do jj=1,NNZ
        ID0=(II-1)*NNZ+jj
        ID=(II-1)*NNZ+NNZ+1-JJ
        PM(1:22,ID)=PM0(1:22,ID0)
    enddo
enddo


do ii=1,NNX
 ! write(91,*)PM(1,1+(ii-1)*NNZ:ii*NNZ)
 ! write(92,*)PM(22,1+(ii-1)*NNZ:ii*NNZ)
 ! write(92,*)PM(3,1+(ii-1)*NNZ:ii*NNZ)
 ! write(93,*)PM(4,1+(ii-1)*NNZ:ii*NNZ)
 ! write(94,*)PM(13,1+(ii-1)*NNZ:ii*NNZ)
 ! write(95,*)PM(17,1+(ii-1)*NNZ:ii*NNZ)
 ! write(96,*)PM(22,1+(ii-1)*NNZ:ii*NNZ)
enddo

do ii=1,NZ
  !write(94,*)RHO1(ii,1:NX)
  !write(95,*)VP1(ii,1:NX)
  !write(96,*)VS1(ii,1:NX)
enddo

end subroutine


!------------------------------------------------------
Subroutine MODEL_CDF(NFS,ITZ,NMD,KT,MT,MNO,K_layer,MFS)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
Dimension KT(NMD),MT(NMD)

Integer,DIMENSION(:) :: KT1(NMD+1)
KT1(1)=1; KT1(2:NMD+1)=KT(1:NMD);
    MNO=MT(1); 
    
    if(NMD>1)then
        do KM=1,NMD
            if(ITZ>=KT1(KM).and.ITZ<=KT1(KM+1))then
            !print *, ITZ, KT(KM), MT(KM)
            MNO=MT(KM)
            endif
            if(KM==NMD.and. ITZ>=KT1(KM+1))MNO=MT(NMD)
        enddo
        
        !-Water-solid layer
        KWS=0; MWS=0
        do KM=1,NMD-1
        if(MT(KM)/=MT(KM+1))then
            MWS=1
            KWS=KT1(KM+1)
        endif
        enddo
        
        if(MWS==1.and.ITZ==KWS)MNO=2
        
    endif
    
end subroutine


!C----------------------------------------------------------------------C
!C                                                                      C
!C     this sburoutine builds up a grid for sub-domain spectral method  C
!C                                                                      C
!C     Entries:                                                         C
!C                                                                      C
!C     (1) IFDM=1 or 0....................create FD or Chbyshev grid;   C
!C     (2) DX,DZ,HIGH,DEPTH............block sizes, air-high $ depth;   C
!C     (3) NEXTD......................block number in Extension zone;   C
!C     (3) NORDX,NORDZ....................number of Chebyshev points;   C
!C     (3) NSF=NSF0+3.......NSF0 is the number of the inter-surfaces    C
!C                          pluse 3 boundaries of the extended zones;   C
!C     (4) NSP(*),XTO(*,*),ZTO(*,*),..............free & interfaces ;   C
!C     (5) NSR,XSR(*),ZSR(*).............sources (S) & recievers (R);   C
!C     (6) ICSR(NSR),VSR(NSR,3,3)........components of S & R vectors;   C
!C                                                                      C
!C                                                                      C
!C     Returns:                                                         C
!C                                                                      C
!C     (1) ASX(NORDX),ASZ(NORDZ).....................Chebyshev points:  C
!C     (2) XM(MX),MZ(MX),ZM(*)..................inner grid for output;  C
!C     (3) NX,X(*),NZ,NZL(*).................sub-domains' coordinates;  C
!C     (5) XP(*),ZP(3,*)........coordinates of total Chebyshev points;  C
!C     (7) XB(NXB),ZB(*),NZB(*)......................the inner points:  C
!C     (8) MSR(NSR),MSR1(NSR,*),FSR(NSR,*)..S&R numbers of the points;  C
!C     (9) IOUT......................=1 for output, -0 dosen't output.  C
!c     (10)IAIR................=1 or 0 for output the air part or not.  C
!C                                                                      C
!C----------------------------------------------------------------------C
SUBROUTINE GRID_2D(DX,DZ,NEXTD,NORDX,NORDZ,NSF,NSP,NFS,XTO,ZTO,XMIN,XMAX,&
              ZMIN,ZMAX,NX,NZ,NZL,XP,ZP,NNX,NNZ)

  IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER LAB*80
      DIMENSION XTO(NSF,NSP),ZTO(NSF,NSP)
      DIMENSION NZL(NSF-1),XP(*),ZP(3,*)

      DIMENSION XTI(:),ZTI(:),XT(:),ZT(:),SF(:,:),X(NX)
      DIMENSION XT1(:),ZT1(:),XT2(:),ZT2(:),A1(:,:),A2(:,:)
      DIMENSION F1(:),F2(:),DHB(:,:),RD(:,:),DMX(:),DMZ(:)

      Dimension DXZ(:,:),DLX_LC(NORDX),IDZ(NORDZ),AASX(NORDX),AASZ(NORDZ)
      dimension ZE_GL(NORDZ),XE_GL(NORDX),IDX(NORDX),DLZ_LC(NORDZ)
      DIMENSION ASX(NORDX),ASZ(NORDZ)

      ALLOCATABLE XTI,ZTI,XT,ZT,XT1,ZT1,XT2,ZT2
      ALLOCATABLE SF,F1,F2,A1,A2,DHB,RD,DMX,DMZ

      allocatable DXZ
!-----------------------------------------
!1) extending the x-direction ----------
!--------------------------------------
      ALLOCATE (DMX(NEXTD))

      DO I=1,NEXTD
      DMX(I)=DX
      ENDDO

      X1=XMIN
      X2=XMAX
      DO I=1,NEXTD
      X1=X1-DMX(I)
      X2=X2+DMX(I)
      ENDDO


      X0=0.5D0*(X1+X2)
      X1=X1-X0
      X2=X2-X0
      XMIN=XMIN-X0
      XMAX=XMAX-X0

      X(1)=X1
      DO I=1,NEXTD
      X(I+1)=X(I)+DMX(NEXTD+1-I)
      ENDDO

      NX0=INT((XMAX-XMIN)/DX)+1
      DO I=1,NX0
      X(I+NEXTD)=X(NEXTD+1)+DFLOAT(I-1)*DX
      ENDDO

      I0=NEXTD+NX0
      DO I=1,NEXTD
      X(I0+I)=X(I0+I-1)+DMX(I)
      ENDDO

      DEALLOCATE (DMX)


!--------------------------------
!2) Extending Z-direcion
!------------------------------
    ALLOCATE(DMZ(NEXTD))
    DMZ(1:NEXTD)=DZ

    ZBOT=ZMIN;ZTOP=ZMAX

    do k=1,NEXTD
      ZBOT=ZBOT-DMZ(k)
      if(NFS==0)ZTOP=ZTOP+DMZ(k)
    enddo
    write(*,*)Zbot, ztop
    !first added interface
    XTO(1,1:NSP)=XTO(2,1:NSP)
    ZTO(1,1:NSP)=ZBOT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    !air top interface
    !XTO(NSF-1,1:NSP)=XTO(NSF-2,1:NSP)
    !ZTO(NSF-1,1:NSP)=ZMAX
    !the end surface in z-direction
    if(NFS==0)then
      XTO(NSF,1:NSP)=XTO(NSF-1,1:NSP)
      ZTO(NSF,1:NSP)=ZTOP
    endif
    !set z0 as the original point

    Z0=ZBOT

    do  i=1,NSF
      ZTO(i,:)=ZTO(i,:)-Z0
    enddo

    ZMIN=ZMIN-Z0
    ZMAX=ZMAX-Z0
    ZBOT=ZBOT-Z0
    ZTOP=ZTOP-Z0

    WRITE(*,1) X(1),X(NX)
    WRITE(*,2) ZBOT,ZTOP
    1  FORMAT('        X-range:',F15.4,1X,F15.4)
    2  FORMAT('        Z-range:',F15.4,1X,F15.4)


!--------------------------------------
!3) NZL
!-------------------------------------
    ALLOCATE (XT1(NSP),ZT1(NSP))
    ALLOCATE (XT2(NSP),ZT2(NSP))

    NZL(  1)=NEXTD  ! 1st layer

    if(NFS==0)NZL(NSF-1)=NEXTD
    NSF1=NSF-1

    DO L=2,NSF1    !mid_interfaces
        XT2(1:NSP)=XTO(L+1,1:NSP)
        ZT2(1:NSP)=ZTO(L+1,1:NSP)
        XT1(1:NSP)=XTO(L,1:NSP)
        ZT1(1:NSP)=ZTO(L,1:NSP)
        H=0.D0
      DO I=1,NX
        XI=X(I)
        H2=ZI(NSP,XT2,ZT2,XI)
        H1=ZI(NSP,XT1,ZT1,XI)
        DH=H2-H1
        H=Max(H,DH)

      ENDDO
      !write(91,*)ZT1(1:NSP)
      !write(92,*)ZT2(1:NSP)
      NZL(L)=MAX(INT(H/DZ),1)
    ENDDO

    NZ=1
    DO L=1,NSF-1
        NZ=NZ+NZL(L)
    ENDDO

!--------------------------
!4) Chebyshev grid
!-------------------------
    NNX=(NX-1)*(NORDX-1)+1
    NNZ=(NZ-1)*(NORDZ-1)+1
    NPT=NNX*NNZ

    PI=3.1415926535897932384626433D0

    FNX=Dfloat(NORDX-1);FNZ=Dfloat(NORDZ-1)

    DO J=1,NORDX      !chebyshev points in x-axis
        FJ=DFLOAT(J-1)
        ASX(J)=DCOS((FNX-FJ)*PI/FNX)
        IF(J==(NORDX+1)/2) ASX(J)=0.d0
    ENDDO

    DO J=1,NORDZ      !chebyshev points in z-axis
      FJ=DFLOAT(J-1)
      ASZ(J)=DCOS((FNZ-FJ)*PI/FNZ)
      IF(J==(NORDX+1)/2) ASZ(J)=0.d0
    ENDDO

  !----------------
  ALLOCATE (F1(NX),F2(NX))
  ALLOCATE (A1(3,NX),A2(3,NX))

!============================
!========================================
  DO  I=1,NX-1   !loop for blocks
    X1=X(I)
    X2=X(I+1)

    F1(1:NX)=0.D0    !top surface of a block
    A1(1:3,1:NX)=0.D0  !spline coeffs b(i) of curves
    K=0
  !--------------
    DO L=1,NSF-1 !loop for layers
      XT2(1:NSP)=XTO(L+1,1:NSP)
      ZT2(1:NSP)=ZTO(L+1,1:NSP)

      XT1(1:NSP)=XTO(L,1:NSP)
      ZT1(1:NSP)=ZTO(L,1:NSP)

    !------------------
      DO L1=1,NZL(L)!loop for blocks
        K=K+1
        NO=(I-1)*(NORDX-1)*NNZ+(K-1)*(NORDZ-1)+1

        !global f2(x) & f1(x)
        !F1(1:NX)=F2(1:NX)
        !A1(1:3,1:NX)=A2(1:3,1:NX)

        !-Choose F2

          DO II=1,NX
            H1=ZI(NSP,XT1,ZT1,X(II))
            H2=ZI(NSP,XT2,ZT2,X(II))
            DM=(H2-H1)/DFLOAT(NZL(L))
            F2(II)=F1(II)+DM
            
          ENDDO
        
        !-Choose A2 from F2
        CALL SPLINE_1D_A(NX,X,F2,A2)
       !----------------------------- 
    !if(NFS==0.and.L==NSF-1)A2=0.d0;A1=0.d0
    
    ! if(I<=NEXTD.or.I>(NX-1)-NEXTD)then
    !    A2=0.d0;A1=0.d0
    ! endif
     NORD1=NORDX-1;NORD3=NORDZ-1
     I2=1;K2=1
     !if(I==(NX-1)-NEXTD.or.I==NX-1)NORD1=NORDX
     !if(I==(NX)-NEXTD)I2=2
     if(K==(NZ-1))NORD3=NORDZ
     if(I==NX-1)NORD1=NORDX
     !-----------------------------
    !==============================
     !---------------------------
      !loop for Ch-points
        DO I1=I2,NORD1
          XI=0.5D0*(X2-X1)*ASX(I1)+0.5D0*(X1+X2)
          II=(I-1)*(NORDX-1)+I1
          XP(II)=XI
            !write(*,*)XP(II)
          !spline for dz/dx
          CALL SPLINE_1D_F(XI,NX,X,F1,A1,Z1,DZX1,DZXX1)
          CALL SPLINE_1D_F(XI,NX,X,F2,A2,Z2,DZX2,DZXX2)
          
          DO K1=K2,NORD3
            ID=NO+(I1-1)*NNZ+K1-1
            !differentiations of zk(x) for every zk
            ZP(1,ID)=0.5D0*(Z2-Z1)*ASZ(K1)+0.5D0*(Z1+Z2)            !zk
            ZP(2,ID)=(DZX2-DZX1)*ASZ(K1)+(DZX1+DZX2)    !dzk/dx
            ZP(3,ID)=(DZXX2-DZXX1)*ASZ(K1)+(DZXX1+DZXX2)    !dzk/dxx
          ENDDO !NORDZ -K1-
        ENDDO !NORDX -I1-
      !-------------------------
    !===============================
      F1(1:NX)=F2(1:NX)
      A1(1:3,1:NX)=A2(1:3,1:NX)
      enddo!NZL -L1-
    enddo!NSF -L-
  enddo!NX -I-
write(*,*)' Grid Generating ends'

!=============================
  write(88,*)XP(1:NNX)
  do k=1,NNX
    write(89,*)ZP(1,1+(k-1)*NNZ:k*NNZ)
  !write(92,*)ZP(2,1+(k-1)*NNZ:k*NNZ)
 ! write(93,*)ZP(3,1+(k-1)*NNZ:k*NNZ)
  enddo
!--------------------------------
DO I=1,NX-1
DO 12 K=1,NZ-1
  NO=(I-1)*(NORDX-1)*NNZ+(K-1)*(NORDZ-1)+1
      !---
  NORD1=NORDX-1; NORD3=NORDZ-1
  IF(I.EQ.NX-1)NORD1=NORDX
  IF(K.EQ.NZ-1)NORD3=NORDZ
  
  DO  I1=1,NORD1     
  DO 10 K1=1,NORD3
  !-------------------
  ITX=(I-1)*(NORDX-1)+I1;ITZ=(K-1)*(NORDZ-1)+K1
    !-------------------------
    ID=NO+(I1-1)*NNZ+(K1-1)    !global index
    IT=(I1-1)*NORDZ+K1         !local  index

!===========================================
  
   call ITYPE(ITX,NORDX,NNX,ISTART,IEND,KIX)
    
    I3=ISTART
    do IX1=1,NORDX
      IDX(IX1)=(I3-1)*NNZ+ITZ
      XE_GL(IX1)=XP(I3)
      I3=I3+1
    enddo
    XM=(XE_GL(NORDX)+XE_GL(1))/2.d0
    XD=(XE_GL(NORDX)-XE_GL(1))/2.d0
    AASX(1:NORDX)=(XE_GL(1:NORDX)-XM)/XD
    DXLC=XD
    
    call CDLI(AASX(KIX),NORDX,AASX,1,DLX_LC)
    !=============================
    DXZ1=0.d0
    do ii=1,NORDX
      DXZ1=DXZ1+DLX_LC(ii)*ZP(1,IDX(ii))
    enddo
    DXZ1=DXZ1/DXLC
    
    DLX_LC=0.d0
    call CDLI(AASX(KIX),NORDX,AASX,2,DLX_LC)
    DXZ2=0.d0
    do ii=1,NORDX
      DXZ2=DXZ2+DLX_LC(ii)*ZP(1,IDX(ii))
    enddo
    
    DXZ2=DXZ2/DXLC/DXLC
    
    ZP(2,ID)=DXZ1;ZP(3,ID)=DXZ2

!-------------------
10  continue
enddo
12  continue 
enddo
   !ZP(2:3,1:NPT)=DXZ(1:2,1:NPT)
   
  ! do k=1,NNX
  !  write(94,*)ZP(2,1+(k-1)*NNZ:k*NNZ)
  !  write(95,*)ZP(3,1+(k-1)*NNZ:k*NNZ)
  !enddo
!---------------------------------------------
    DEALLOCATE (F1,F2,A1,A2)
    DEALLOCATE (XT2,ZT2)
    DEALLOCATE (XT1,ZT1)

      RETURN
      END


!--------------------------
!
!   FN:ZI
!
!-------------------------
      FUNCTION ZI(NSP,XTO,ZTO,XI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XTO(*),ZTO(*),A(:,:)
      ALLOCATABLE A

      IF(XI.LE.XTO(1))THEN
      ZI=ZTO(1)
      GO TO 5
      ENDIF

      IF(XI.GE.XTO(NSP))THEN
      ZI=ZTO(NSP)
      GO TO 5
      ENDIF

      ALLOCATE (A(1:3,NSP))
      CALL SPLINE_1D_A(NSP,XTO,ZTO,A)
      CALL SPLINE_1D_F(XI,NSP,XTO,ZTO,A,FI,FX,FXX)
      ZI=FI

      DEALLOCATE(A)
   5  RETURN
      END

!====================================
!
!
!
!------------------------------------
SUBROUTINE VTI(INDP,P,C)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION P(INDP),C(21),S(3,3)
      dimension AV(6,6),AT(6,6)

      PI=3.1415926535897932384626433D0
      AT(1:6,1:6)=0.d0;
    !---------------------------------

if(INDP==1)then
  AT(1,1)=P(1)
endif

if(INDP==2)then
  AT(1,1)=P(1)+2*P(2);  AT(1,2)=P(1); AT(1,3)=P(1)
                        AT(2,2)=P(1)+2*P(2);  AT(2,3)=P(1);
                                      AT(3,3)=P(1)+2*P(2);
  AT(4,4)=P(2); AT(5,5)=P(2); AT(6,6)=P(2);
endif

!= VTI1, 2, ORT
if(INDP==3)then !VTI2 case
  AT(1,1)=P(1);  AT(1,2)=P(1)-2*P(3); AT(1,3)=P(1)-2*P(2)
                        AT(2,2)=P(1);  AT(2,3)=P(1)-2*P(2);
                                      AT(3,3)=P(1);
  AT(4,4)=P(2); AT(5,5)=P(2); AT(6,6)=P(3);
endif

if(INDP==4)then !VTI1 case
  AT(1,1)=P(1);  AT(1,2)=P(1)-2*P(4); AT(1,3)=-P(3)
                        AT(2,2)=P(1);  AT(2,3)=-P(3);
                                      AT(3,3)=P(2);
  AT(4,4)=P(3); AT(5,5)=P(3); AT(6,6)=P(4);
endif

if(INDP==6)then !ORT case
  AT(1,1)=P(1);  AT(1,2)=-P(6); AT(1,3)=-P(5)
                        AT(2,2)=P(2);  AT(2,3)=-P(4);
                                      AT(3,3)=P(3);
  AT(4,4)=P(4); AT(5,5)=P(5); AT(6,6)=P(6);
endif

28  IP=0
    DO I=1,6
    DO J=I,6
      IP=IP+1
      C(IP)=AT(I,J)
    ENDDO
    ENDDO

RETURN
END

!====================================
!
!
!
!------------------------------------
SUBROUTINE TTI(INDP,P,ITTI,ZTA0,C)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION P(INDP),C(21),S(3,3)
      dimension AV(6,6),AT(6,6)

      PI=3.1415926535897932384626433D0
      AT(1:6,1:6)=0.d0;
    !---------------------------------
AV(1:6,1:6)=0.d0

if(INDP==1)then
  AV(1,1)=P(1)
endif

if(INDP==2)then
  AV(1,1)=P(1)+2*P(2);  AV(1,2)=P(1); AV(1,3)=P(1)
  AV(2,1)=P(1); AV(2,2)=P(1)+2*P(2);  AV(2,3)=P(1);
  AV(3,1)=P(1); AV(3,2)=P(1); AV(3,3)=P(1)+2*P(2);
  AV(4,4)=P(2); AV(5,5)=P(2); AV(6,6)=P(2);
endif

!= VTI1, 2, ORT
if(INDP==3)then !VTI2 case
  AV(1,1)=P(1);  AV(1,2)=P(1)-2*P(3); AV(1,3)=P(1)-2*P(2)
  AV(2,1)=P(1)-2*P(3); AV(2,2)=P(1);  AV(2,3)=P(1)-2*P(2);
  Av(3,1)=P(1)-2*P(2); AV(3,2)=P(1)-2*P(2); AV(3,3)=P(1);
  AV(4,4)=P(2); AV(5,5)=P(2); AV(6,6)=P(3);
endif

if(INDP==4)then !VTI1 case
  AV(1,1)=P(1);  AV(1,2)=P(1)-2*P(4); AV(1,3)=-P(3)
  AV(2,1)=P(1)-2*P(4); AV(2,2)=P(1);  AV(2,3)=-P(3);
  AV(3,1)=-P(3); AV(3,2)=-P(3);     AV(3,3)=P(2);
  AV(4,4)=P(3); AV(5,5)=P(3); AV(6,6)=P(4);
endif

if(INDP==6)then !ORT case
  AV(1,1)=P(1);  AV(1,2)=-P(6); AV(1,3)=-P(5)
  AV(2,1)=-P(6); AV(2,2)=P(2);  AV(2,3)=-P(4);
  AV(3,1)=-P(5); AV(3,2)=-P(4); AV(3,3)=P(3);
  AV(4,4)=P(4); AV(5,5)=P(5); AV(6,6)=P(6);
endif

if(INDP==5)then !GENERAL VTI case
    !AV(1:6,1:6)=0.d0
    AV(1,1)=P(1);AV(1,2)=P(1)-2*P(5);AV(1,3)=P(2);
    AV(2,1)=P(1)-2*P(5);AV(2,2)=P(1);AV(2,3)=P(2);
    AV(3,1)=P(2); AV(3,2)=P(2);AV(3,3)=P(3);
    AV(4,4)=P(4);AV(5,5)=P(4);AV(6,6)=P(5);
endif

if(ITTI==0)then
    AT(:,:)=AV(:,:)
    goto 28
endif
   
    ZTA=atan(ZTA0)
    S(1:3,1:3)=0.d0
    S(1,1)=DCOS(ZTA);  S(1,3)=-DSIN(ZTA)
    S(2,2)=1.d0;
    S(3,1)=DSIN(ZTA); S(3,3)=DCOS(ZTA)
      
     AT(:,:)=0.d0
DO I=1,3
DO J=1,3
  IJ=0
  IF(I.EQ.J)IJ=1
  II=I*IJ+(9-I-J)*(1-IJ)

  DO K=1,3
  DO 30 L=1,3
    KL=0
    IF(K.EQ.L)KL=1
    JJ=K*KL+(9-K-L)*(1-KL)
  if(AT(II,JJ)/=0.d0)goto 30
  !AT(II,JJ)=0.d0

AT(II,JJ)=AV(1,1)*(S(1,I)*S(1,J)*S(1,K)*S(1,L))&
        +AV(1,2)*(S(1,I)*S(1,J)*S(2,K)*S(2,L)+S(2,I)*S(2,J)*S(1,K)*S(1,L))&
        +AV(1,3)*(S(1,I)*S(1,J)*S(3,K)*S(3,L)+S(3,I)*S(3,J)*S(1,K)*S(1,L))&
        +AV(2,2)*(S(2,I)*S(2,J)*S(2,K)*S(2,L))&
        +AV(2,3)*(S(2,I)*S(2,J)*S(3,K)*S(3,L)+S(3,I)*S(3,J)*S(2,K)*S(2,L))&
        +AV(3,3)*(S(3,I)*S(3,J)*S(3,K)*S(3,L))&
        +AV(4,4)*(S(2,I)*S(3,J)*S(2,K)*S(3,L)+S(2,I)*S(3,J)*S(3,K)*S(2,L)&
                +S(3,I)*S(2,J)*S(3,K)*S(2,L)+S(3,I)*S(2,J)*S(2,K)*S(3,L))&
        +AV(5,5)*(S(1,I)*S(3,J)*S(1,K)*S(3,L)+S(1,I)*S(3,J)*S(3,K)*S(1,L)&
                +S(3,I)*S(1,J)*S(3,K)*S(1,L)+S(3,I)*S(1,J)*S(1,K)*S(3,L))&
        +AV(6,6)*(S(2,I)*S(1,J)*S(2,K)*S(1,L)+S(2,I)*S(1,J)*S(1,K)*S(2,L)&
                +S(1,I)*S(2,J)*S(1,K)*S(2,L)+S(1,I)*S(2,J)*S(2,K)*S(1,L))
  

  if(II/=JJ)AT(JJ,II)=AT(II,JJ)
  !----------------------------
  !==============================


30 continue
   ENDDO
!-----------------------------
  ENDDO
  ENDDO


28  IP=0
    DO I=1,6
    DO J=I,6
      IP=IP+1
      C(IP)=AT(I,J)
    ENDDO
    ENDDO

RETURN
END
