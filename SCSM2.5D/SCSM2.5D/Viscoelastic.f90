!-----------------
Subroutine Memory_var_T3(ID,IREC,NPT,NSLS,dt,DX_Y,DY_Y,DZ_Y,NORDX,IDX,NORDZ,&
            IDZ,DFS,tau_sig,tau_ep,AT_memory,V_kl,Var_Memory,Var_Memory_T)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    !-Input
    Dimension DX_Y(13),DY_Y(13),DZ_Y(13)
    dimension IDX(NPT,NORDX),IDZ(NPT,NORDZ),DFS(NPT,2,NORDX,4)
    dimension TAU_SIG(NPT,NSLS,21),TAU_ep(NPT,NSLS,21)
    dimension AT_memory(NPT,NSLS,6,6)
    dimension V_kl(NPT,6)
    dimension Var_Memory(NPT,Nsls,7)
    dimension Var_Memory_T(NPT,Nsls,7)
    
    !-In subroutine
    dimension DF_V(6),TEMP_V(7),ETA(7)
    double precision :: V_kk,V_kk1
    
    !-----------------------
    !-Point location 
    IX=0; IZ=0
    do ii=1,NORDX
        if(ID==IDX(ID,II))IX=ii
    enddo
    
    do jj=1,NORDZ
        if(ID==IDZ(ID,jj))IZ=jj
    enddo
    
    if(IX==0.or.IZ==0)then
        print *, 'Memory_var_T: Point not searched', ID
        print *, IDX(ID,1:NORDX)
        print *, IDZ(ID,1:NORDZ)
        stop
    endif
        
    
    !DxVx
    DF_V(1)=DX_Y(5)
    
    !DyVx
    DF_V(2)=(DY_Y(5)+DX_Y(6))/2.0
    
    !DyVy
    DF_V(3)=DY_Y(6)
    
    !DzVx
    DF_V(4)=(DZ_Y(5)+DX_Y(7))/2.0
    
    !DzVy
    DF_V(5)=(DZ_Y(6)+DY_Y(7))/2.0
    
    !DzVz
    DF_V(6)=DZ_Y(7)
    
    TEMP_V(1:6)=(DF_V(1:6)+V_kl(ID,1:6))/2.d0
    !TEMP_V(2:7)=(DF_V(1:6)+V_kl(ID,1:6))/2.d0
    
    do l=1,NSLS
     !ETA(1:6)=matmul(AT_MEMORY(ID,l,1:6,1:6),TEMP_V(1:6))
     
     ETA(1)=AT_MEMORY(ID,l,1,1)*TEMP_V(1)+ AT_MEMORY(ID,l,1,3)*TEMP_V(3) &
            +AT_MEMORY(ID,l,1,4)*TEMP_V(4)+AT_MEMORY(ID,l,1,6)*TEMP_V(6)
     
     ETA(2)=AT_MEMORY(ID,l,2,2)*TEMP_V(2)+AT_MEMORY(ID,l,2,5)*TEMP_V(5)
     
     ETA(3)=AT_MEMORY(ID,l,3,1)*TEMP_V(1)+ AT_MEMORY(ID,l,3,3)*TEMP_V(3) &
            + AT_MEMORY(ID,l,3,4)*TEMP_V(4)+AT_MEMORY(ID,l,3,6)*TEMP_V(6)
            
     ETA(4)=AT_MEMORY(ID,l,4,1)*TEMP_V(1)+AT_MEMORY(ID,l,4,3)*TEMP_V(3)&
            +AT_MEMORY(ID,l,4,4)*TEMP_V(4)+AT_MEMORY(ID,l,4,6)*TEMP_V(6)
     
     ETA(5)=AT_MEMORY(ID,l,5,2)*TEMP_V(2)+AT_MEMORY(ID,l,5,5)*TEMP_V(5)
     
     ETA(6)=AT_MEMORY(ID,l,6,1)*TEMP_V(1)+ AT_MEMORY(ID,l,6,3)*TEMP_V(3) &
            + AT_MEMORY(ID,l,6,4)*TEMP_V(4)+AT_MEMORY(ID,l,6,6)*TEMP_V(6)
     
     Var_Memory(id,l,1:6)=( (ETA(1:6)-VAR_MEMORY_T(ID,l,1:6)/2.d0/TAU_SIG(ID,L,1))*dt&
        +VAR_MEMORY_T(ID,L,1:6) ) / (1+dt/TAU_sig(ID,L,1)/2.d0)
    enddo
    
    V_kl(ID,1:6)=DF_V(1:6)
    
return
end subroutine

!-----------------
Subroutine Memory_Matrix_VTI(ID,IREC,NPT,Nani,NSLS,PD,tau_sig,tau_ep,tau_TTI,AT_Memory)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    !-Input
    dimension TAU_SIG(NPT,NSLS,21),TAU_ep(NPT,NSLS,21),TAU_TTI(NPT,NSLS,21),PD(22)
    DIMENSION AT_MEMORY(NPT,NSLS,6,6)
    
    AT_MEMORY(ID,1:NSLS,1:6,1:6)=0.d0
    

    do i=1,NSLS
    AT_MEMORY(ID,i,1,1)=TAU_TTI(ID,i,1); AT_MEMORY(ID,i,1,4)=2.0*TAU_TTI(ID,i,5)
    AT_MEMORY(ID,i,1,3)=TAU_TTI(ID,i,2); AT_MEMORY(ID,i,1,6)=TAU_TTI(ID,i,3)
    
    AT_MEMORY(ID,i,2,2)=2.0*TAU_TTI(ID,i,21); AT_MEMORY(ID,i,2,5)=2.0*TAU_TTI(ID,i,18)
    
    AT_MEMORY(ID,i,3,1)=TAU_TTI(ID,i,2); AT_MEMORY(ID,i,3,3)=TAU_TTI(ID,i,7)
    AT_MEMORY(ID,i,3,4)=2.0*TAU_TTI(ID,i,10); AT_MEMORY(ID,i,3,6)=TAU_TTI(ID,i,8)
    
    AT_MEMORY(ID,i,4,1)=TAU_TTI(ID,i,5); AT_MEMORY(ID,i,4,3)=TAU_TTI(ID,i,10)
    AT_MEMORY(ID,i,4,4)=2.0*TAU_TTI(ID,i,19); AT_MEMORY(ID,i,4,6)=TAU_TTI(ID,i,14)
    
    AT_MEMORY(ID,i,5,2)=2.0*TAU_TTI(ID,i,18); AT_MEMORY(ID,i,5,5)=2.0*TAU_TTI(ID,i,16)
    
    AT_MEMORY(ID,i,6,1)=TAU_TTI(ID,i,3); AT_MEMORY(ID,i,6,3)=TAU_TTI(ID,i,8)
    AT_MEMORY(ID,i,6,4)=2.0*TAU_TTI(ID,i,14); AT_MEMORY(ID,i,6,6)=TAU_TTI(ID,i,12)
    enddo

    
    
    AT_MEMORY(ID,1:NSLS,1:6,1:6)=AT_MEMORY(ID,1:NSLS,1:6,1:6)/dble(NSLS)
return
end subroutine

!-----------------
Subroutine Memory_Matrix(ID,IREC,NPT,Nani,NSLS,PD,tau_sig,tau_ep,AT_Memory)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    !-Input
    dimension TAU_SIG(NPT,NSLS,21),TAU_ep(NPT,NSLS,21),PD(22)
    DIMENSION AT_MEMORY(NPT,NSLS,6,6)
    
    AT_MEMORY(ID,1:NSLS,1:6,1:6)=0.d0
    
if(Nani==2)then !Isotropic
    fac_lambda=1.d0+sum(TAU_EP(ID,1:NSLS,1)/TAU_SIG(ID,1:NSLS,1)-1.d0)/dble(NSLS)
    fac_mu=1.d0+sum(TAU_EP(ID,1:NSLS,2)/TAU_SIG(ID,1:NSLS,2)-1.d0)/dble(NSLS)
    Relaxed_Lambda=PD(3)/fac_lambda;
    Relaxed_Mu=PD(22)/fac_mu;
    
    do i=1,NSLS
    fac_l=(1.d0-TAU_EP(ID,i,1)/TAU_SIG(ID,i,1))/TAU_SIG(ID,i,1)*Relaxed_Lambda
    fac_m=(1.d0-TAU_EP(ID,i,2)/TAU_SIG(ID,i,2))/TAU_SIG(ID,i,2)*Relaxed_Mu
    
    !AT_MEMORY(ID,i,1,1)=fac_l
    
    AT_MEMORY(ID,i,1,1)=fac_l+2.d0*fac_m
    AT_MEMORY(ID,i,1,3)=fac_l
    AT_MEMORY(ID,i,1,6)=fac_l
    
    AT_MEMORY(ID,i,2,2)=2.d0*fac_m
    
    
    AT_MEMORY(ID,i,3,1)=fac_l
    AT_MEMORY(ID,i,3,3)=fac_l+2.d0*fac_m
    AT_MEMORY(ID,i,3,6)=fac_l
    
    AT_MEMORY(ID,i,4,4)=2.d0*fac_m
    
    AT_MEMORY(ID,i,5,5)=2.d0*fac_m
    
    AT_MEMORY(ID,i,6,1)=fac_l
    AT_MEMORY(ID,i,6,3)=fac_l
    AT_MEMORY(ID,i,6,6)=fac_l+2.d0*fac_m
    enddo
    
elseif(Nani==3)then !VTI2
    fac_a11=1.d0+sum(TAU_EP(ID,1:NSLS,1)/TAU_SIG(ID,1:NSLS,1)-1.d0)/dble(NSLS)
    fac_a44=1.d0+sum(TAU_EP(ID,1:NSLS,2)/TAU_SIG(ID,1:NSLS,2)-1.d0)/dble(NSLS)
    fac_a66=1.d0+sum(TAU_EP(ID,1:NSLS,3)/TAU_SIG(ID,1:NSLS,3)-1.d0)/dble(NSLS)
    Relaxed_a11=PD(2)/fac_a11;
    Relaxed_a44=PD(20)/fac_a44;
    Relaxed_a66=PD(22)/fac_a66;
    
    do i=1,NSLS
    fac_11=(1.d0-TAU_EP(ID,i,1)/TAU_SIG(ID,i,1))/TAU_SIG(ID,i,1)*Relaxed_a11
    fac_44=(1.d0-TAU_EP(ID,i,2)/TAU_SIG(ID,i,2))/TAU_SIG(ID,i,2)*Relaxed_a44
    fac_66=(1.d0-TAU_EP(ID,i,3)/TAU_SIG(ID,i,3))/TAU_SIG(ID,i,3)*Relaxed_a66
    
    AT_MEMORY(ID,i,1,1)=fac_11
    AT_MEMORY(ID,i,1,3)=fac_11-2.d0*fac_66
    AT_MEMORY(ID,i,1,6)=fac_11-2.d0*fac_44
    
    AT_MEMORY(ID,i,2,2)=2.d0*fac_66
    
    
    AT_MEMORY(ID,i,3,1)=fac_11-2.d0*fac_66
    AT_MEMORY(ID,i,3,3)=fac_11
    AT_MEMORY(ID,i,3,6)=fac_11-2.d0*fac_44
    
    AT_MEMORY(ID,i,4,4)=2.d0*fac_44
    
    AT_MEMORY(ID,i,5,5)=2.d0*fac_44
    
    AT_MEMORY(ID,i,6,1)=fac_11-2.d0*fac_44
    AT_MEMORY(ID,i,6,3)=fac_11-2.d0*fac_44
    AT_MEMORY(ID,i,6,6)=fac_11
    enddo

elseif(Nani==4)then !VTI1
    fac_a11=1.d0+sum(TAU_EP(ID,1:NSLS,1)/TAU_SIG(ID,1:NSLS,1)-1.d0)/dble(NSLS)
    fac_a33=1.d0+sum(TAU_EP(ID,1:NSLS,2)/TAU_SIG(ID,1:NSLS,2)-1.d0)/dble(NSLS)
    fac_a44=1.d0+sum(TAU_EP(ID,1:NSLS,3)/TAU_SIG(ID,1:NSLS,3)-1.d0)/dble(NSLS)
    fac_a66=1.d0+sum(TAU_EP(ID,1:NSLS,4)/TAU_SIG(ID,1:NSLS,4)-1.d0)/dble(NSLS)
    
    Relaxed_a11=PD(2)/fac_a11;
    Relaxed_a33=PD(13)/fac_a33;
    Relaxed_a44=PD(20)/fac_a44;
    Relaxed_a66=PD(22)/fac_a66;
    
    do i=1,NSLS
    fac_11=(1.d0-TAU_EP(ID,i,1)/TAU_SIG(ID,i,1))/TAU_SIG(ID,i,1)*Relaxed_a11
    fac_33=(1.d0-TAU_EP(ID,i,2)/TAU_SIG(ID,i,2))/TAU_SIG(ID,i,2)*Relaxed_a33
    fac_44=(1.d0-TAU_EP(ID,i,3)/TAU_SIG(ID,i,3))/TAU_SIG(ID,i,3)*Relaxed_a44
    fac_66=(1.d0-TAU_EP(ID,i,4)/TAU_SIG(ID,i,4))/TAU_SIG(ID,i,4)*Relaxed_a66
    
    AT_MEMORY(ID,i,1,1)=fac_11
    AT_MEMORY(ID,i,1,3)=fac_11-2.d0*fac_66
    AT_MEMORY(ID,i,1,6)=-fac_44
    
    AT_MEMORY(ID,i,2,2)=2.d0*fac_66
    
    
    AT_MEMORY(ID,i,3,1)=fac_11-2.d0*fac_66
    AT_MEMORY(ID,i,3,3)=fac_11
    AT_MEMORY(ID,i,3,6)=-fac_44
    
    AT_MEMORY(ID,i,4,4)=2.d0*fac_44
    
    AT_MEMORY(ID,i,5,5)=2.d0*fac_44
    
    AT_MEMORY(ID,i,6,1)=-fac_44
    AT_MEMORY(ID,i,6,3)=-fac_44
    AT_MEMORY(ID,i,6,6)=fac_33
    enddo
elseif(Nani==6)then !ORT
    fac_a11=1.d0+sum(TAU_EP(ID,1:NSLS,1)/TAU_SIG(ID,1:NSLS,1)-1.d0)/dble(NSLS)
    fac_a22=1.d0+sum(TAU_EP(ID,1:NSLS,2)/TAU_SIG(ID,1:NSLS,2)-1.d0)/dble(NSLS)
    fac_a33=1.d0+sum(TAU_EP(ID,1:NSLS,3)/TAU_SIG(ID,1:NSLS,3)-1.d0)/dble(NSLS)
    fac_a44=1.d0+sum(TAU_EP(ID,1:NSLS,4)/TAU_SIG(ID,1:NSLS,4)-1.d0)/dble(NSLS)
    fac_a55=1.d0+sum(TAU_EP(ID,1:NSLS,5)/TAU_SIG(ID,1:NSLS,5)-1.d0)/dble(NSLS)
    fac_a66=1.d0+sum(TAU_EP(ID,1:NSLS,6)/TAU_SIG(ID,1:NSLS,6)-1.d0)/dble(NSLS)
    
    Relaxed_a11=PD(2)/fac_a11;
    Relaxed_a22=PD(8)/fac_a22;
    Relaxed_a33=PD(13)/fac_a33;
    Relaxed_a44=PD(17)/fac_a44;
    Relaxed_a55=PD(20)/fac_a55;
    Relaxed_a66=PD(22)/fac_a66;
    
    do i=1,NSLS
    fac_11=(1.d0-TAU_EP(ID,i,1)/TAU_SIG(ID,i,1))/TAU_SIG(ID,i,1)*Relaxed_a11
    fac_22=(1.d0-TAU_EP(ID,i,2)/TAU_SIG(ID,i,2))/TAU_SIG(ID,i,2)*Relaxed_a22
    fac_33=(1.d0-TAU_EP(ID,i,3)/TAU_SIG(ID,i,3))/TAU_SIG(ID,i,3)*Relaxed_a33
    fac_44=(1.d0-TAU_EP(ID,i,4)/TAU_SIG(ID,i,4))/TAU_SIG(ID,i,4)*Relaxed_a44
    fac_55=(1.d0-TAU_EP(ID,i,5)/TAU_SIG(ID,i,5))/TAU_SIG(ID,i,5)*Relaxed_a55
    fac_66=(1.d0-TAU_EP(ID,i,6)/TAU_SIG(ID,i,6))/TAU_SIG(ID,i,6)*Relaxed_a66
    
    AT_MEMORY(ID,i,1,1)=fac_11
    AT_MEMORY(ID,i,1,3)=-fac_66
    AT_MEMORY(ID,i,1,6)=-fac_55
    
    AT_MEMORY(ID,i,2,2)=2.d0*fac_66
    
    AT_MEMORY(ID,i,3,1)=-fac_66
    AT_MEMORY(ID,i,3,3)=fac_22
    AT_MEMORY(ID,i,3,6)=-fac_44
    
    AT_MEMORY(ID,i,4,4)=2.d0*fac_55
    
    AT_MEMORY(ID,i,5,5)=2.d0*fac_44
    
    AT_MEMORY(ID,i,6,1)=-fac_55
    AT_MEMORY(ID,i,6,3)=-fac_44
    AT_MEMORY(ID,i,6,6)=fac_33
    enddo
endif
    
    
    AT_MEMORY(ID,1:NSLS,1:6,1:6)=AT_MEMORY(ID,1:NSLS,1:6,1:6)/dble(NSLS)
return
end subroutine

!----------------------
Subroutine Memory_var_WS_T(ID,ITX,IREC,NNX,NPT,NSLS,wky,dt,ZP,DX_Y,DY_Y,DZ_Y,&
DX_Y_WS,DY_Y_WS,DZ_Y_WS,NORDX,IDX,NORDZ,IDZ,DFS,PD,tau_sig,tau_ep,V_kl,Var_Memory)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    !-Input
    Dimension ZP(3,NPT)
    Dimension DX_Y(13),DY_Y(13),DZ_Y(13),DX_Y_WS(13),DY_Y_WS(13),DZ_Y_WS(13)
    dimension IDX(NPT,NORDX),IDZ(NPT,NORDZ),DFS(NPT,2,NORDX,4)
    dimension TAU_SIG(NPT,NSLS,21), TAU_ep(NPT,NSLS,21),PD(22,NPT)
    dimension V_kl(NNX,6)
    dimension Var_Memory(NNX,Nsls,7)
    
    !-In subroutine
    dimension DF_V(3)
    
    
    fac=1.d0+sum(TAU_EP(ID+1,1:NSLS,1)/TAU_SIG(ID+1,1:NSLS,1)-1.d0)/dble(NSLS)
    Bulk_Me=PD(2,ID+1)/fac;
    
    
    !-----------------------
    !-Point location 
    IX=0; IZ=0
    do ii=1,NORDX
        if(ID==IDX(ID,II))IX=ii
    enddo
    
    do jj=1,NORDZ
        if(ID==IDZ(ID,jj))IZ=jj
    enddo
    
    if(IX==0.or.IZ==0)then
        print *, 'Memory_var_T: Point not searched', ID
        print *, IDX(ID,1:NORDX)
        print *, IDZ(ID,1:NORDZ)
        stop
    endif
        
    
    !DxVx
    DF_V(1)=DX_Y_WS(1)
    !DF_V(1)=DF_V(1)+DFS(ID,2,IZ,3)*DZ_Y_WS(1)
    
    !DyVy
    DF_V(2)=DY_Y_WS(2)
    
    !DzVz
    !DF_V(3)=DFS(ID,2,IZ,4)*DZ_Y_WS(3)
    DF_V(3)=DZ_Y(7)+ZP(2,ID)*(DZ_Y_WS(1)-DZ_Y(5))
    
    !Leap Frog
    do l=1,NSLS
        Var_Memory(ITX,l,2)=((Bulk_Me/dble(NSLS)*(1.d0-TAU_ep(ID+1,l,1)/TAU_sig(ID+1,l,1))&
        *(sum(DF_V)+V_kl(ITX,1))/2.d0-VAR_MEMORY(ITX,l,1)/2.d0)/TAU_SIG(ID+1,l,1)*dt&
        +VAR_MEMORY(ITX,l,1))/(1+dt/TAU_sig(ID+1,l,1)/2.d0)
    enddo
    
    V_kl(ITX,1)=sum(DF_V)
    
return
end subroutine

!----------------------
Subroutine Memory_var_W_T2(ID,IREC,NPT,NSLS,wky,dt,DX_Y,DY_Y,DZ_Y,NORDX,IDX,NORDZ,&
IDZ,DFS,PD,tau_sig,tau_ep,V_kl,Var_Memory,Var_Memory_T)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    !-Input
    Dimension DX_Y(13),DY_Y(13),DZ_Y(13)
    dimension IDX(NPT,NORDX),IDZ(NPT,NORDZ),DFS(NPT,2,NORDX,4)
    dimension TAU_SIG(NPT,NSLS,21), TAU_ep(NPT,NSLS,21),PD(22,NPT)
    dimension V_kl(NPT,6)
    dimension Var_Memory(NPT,Nsls,7)
    dimension Var_Memory_T(NPT,Nsls,7)
    
    !-In subroutine
    dimension DF_V(3)
    
    
    fac=1.d0+sum(TAU_EP(ID,1:NSLS,1)/TAU_SIG(ID,1:NSLS,1)-1.d0)/dble(NSLS)
    Bulk_Me=PD(2,ID)/fac;
    
    
    !-----------------------
    !-Point location 
    IX=0; IZ=0
    do ii=1,NORDX
        if(ID==IDX(ID,II))IX=ii
    enddo
    
    do jj=1,NORDZ
        if(ID==IDZ(ID,jj))IZ=jj
    enddo
    
    if(IX==0.or.IZ==0)then
        print *, 'Memory_var_T: Point not searched', ID
        print *, IDX(ID,1:NORDX)
        print *, IDZ(ID,1:NORDZ)
        stop
    endif
        
    
    !DxVx
    DF_V(1)=DX_Y(1)
    !DF_V(1)=DF_V(1)+DFS(ID,2,IZ,3)*DZ_Y(1)
    
    !DyVy
    DF_V(2)=DY_Y(2)
    
    !DzVz
    DF_V(3)=DZ_Y(3)
    
    !Leap Frog
    do l=1,NSLS
        Var_Memory(id,l,2)=((Bulk_Me/dble(NSLS)*(1.d0-TAU_ep(ID,l,1)/TAU_sig(ID,l,1))&
        *(sum(DF_V)+V_kl(ID,1))/2.d0-VAR_MEMORY(ID,l,1)/2.d0)/TAU_SIG(ID,l,1)*dt&
        +VAR_MEMORY(ID,l,1))/(1+dt/TAU_sig(ID,l,1)/2.d0)
    enddo
    
    V_kl(ID,1)=sum(DF_V)
    
return
end subroutine

