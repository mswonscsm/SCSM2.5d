!-----------------
Subroutine TSRC_2(ID,IREC,NPT,NSLS,dt,DX_Y,DY_Y,DZ_Y,NORDX,IDX,NORDZ,&
            IDZ,DFS,tau_sig,tau_ep,TSRC_ABC,AT_memory,V_kl,Var_Memory,Var_Memory_T)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    !-Input
    Dimension DX_Y(13),DY_Y(13),DZ_Y(13)
    dimension IDX(NPT,NORDX),IDZ(NPT,NORDZ),DFS(NPT,2,NORDX,4)
    dimension TAU_SIG(NPT,NSLS,21),TAU_ep(NPT,NSLS,21)
    dimension AT_memory(NPT,NSLS,6,6),TSRC_ABC(NPT,NSLS,4)
    dimension V_kl(NPT,6)
    dimension Var_Memory(NPT,Nsls,7)
    dimension Var_Memory_T(2,NPT,Nsls,7)
    
    !-In subroutine
    double precision, dimension(:):: DF_V(6),TEMP_V(7),ETA(7)
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
    
    do l=1,NSLS
     !ETA(1:6)=matmul(AT_MEMORY(ID,l,1:6,1:6),TEMP_V(1:6))
     
     ETA(1)=AT_MEMORY(ID,l,1,1)*DF_V(1)+ AT_MEMORY(ID,l,1,3)*DF_V(3) &
            + AT_MEMORY(ID,l,1,4)*DF_V(4)+AT_MEMORY(ID,l,1,6)*DF_V(6)
     
     ETA(2)=AT_MEMORY(ID,l,2,2)*DF_V(2)+AT_MEMORY(ID,l,2,5)*DF_V(5)
     
     ETA(3)=AT_MEMORY(ID,l,3,1)*DF_V(1)+ AT_MEMORY(ID,l,3,3)*DF_V(3) &
            +AT_MEMORY(ID,l,3,4)*DF_V(4)+AT_MEMORY(ID,l,3,6)*DF_V(6)
            
     ETA(4)=AT_MEMORY(ID,l,4,1)*DF_V(1)+AT_MEMORY(ID,l,4,3)*DF_V(3)&
            +AT_MEMORY(ID,l,4,4)*DF_V(4)+AT_MEMORY(ID,l,4,6)*DF_V(6)
     
     ETA(5)=AT_MEMORY(ID,l,5,2)*DF_V(2)+AT_MEMORY(ID,l,5,5)*DF_V(5)
     
     ETA(6)=AT_MEMORY(ID,l,6,1)*DF_V(1)+ AT_MEMORY(ID,l,6,3)*DF_V(3) &
            +AT_MEMORY(ID,l,6,4)*DF_V(4)+AT_MEMORY(ID,l,6,6)*DF_V(6)
     
     Var_Memory(id,l,1:6)=TSRC_ABC(ID,l,1)*Var_Memory(id,l,1:6)+TSRC_ABC(ID,l,2)*Var_Memory_T(1,ID,l,1:6) &
                        +TSRC_ABC(ID,l,3)*(ETA(1:6)-Var_Memory_T(2,ID,l,1:6))&
                        +TSRC_ABC(ID,l,4)*(ETA(1:6)-2.0*Var_Memory_T(1,ID,l,1:6)+Var_Memory_T(2,ID,l,1:6))
     
     Var_Memory_T(2,ID,l,1:6)=Var_Memory_T(1,ID,l,1:6)
     Var_Memory_T(1,ID,l,1:6)=ETA(1:6)                       
    enddo
    
    !V_kl(ID,1:6)=DF_V(1:6)
    
return
end subroutine

!-----------------
Subroutine TSRC_1(ID,IREC,NPT,NSLS,dt,DX_Y,DY_Y,DZ_Y,NORDX,IDX,NORDZ,&
            IDZ,DFS,tau_sig,tau_ep,TSRC_ABC,AT_memory,V_kl,Var_Memory,Var_Memory_T)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    !-Input
    Dimension DX_Y(13),DY_Y(13),DZ_Y(13)
    dimension IDX(NPT,NORDX),IDZ(NPT,NORDZ),DFS(NPT,2,NORDX,4)
    dimension TAU_SIG(NPT,NSLS,21),TAU_ep(NPT,NSLS,21)
    dimension AT_memory(NPT,NSLS,6,6),TSRC_ABC(NPT,NSLS,4)
    dimension V_kl(NPT,6)
    dimension Var_Memory(NPT,Nsls,7)
    dimension Var_Memory_T(NPT,Nsls,7)
    
    !-In subroutine
    double precision, dimension(:):: DF_V(6),TEMP_V(7),ETA(7)
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
    
    do l=1,NSLS
     !ETA(1:6)=matmul(AT_MEMORY(ID,l,1:6,1:6),TEMP_V(1:6))
     
     ETA(1)=AT_MEMORY(ID,l,1,1)*DF_V(1)+ AT_MEMORY(ID,l,1,3)*DF_V(3) &
            + AT_MEMORY(ID,l,1,4)*DF_V(4)+AT_MEMORY(ID,l,1,6)*DF_V(6)
     
     ETA(2)=AT_MEMORY(ID,l,2,2)*DF_V(2)+AT_MEMORY(ID,l,2,5)*DF_V(5)
     
     ETA(3)=AT_MEMORY(ID,l,3,1)*DF_V(1)+ AT_MEMORY(ID,l,3,3)*DF_V(3) &
            +AT_MEMORY(ID,l,3,4)*DF_V(4)+AT_MEMORY(ID,l,3,6)*DF_V(6)
            
     ETA(4)=AT_MEMORY(ID,l,4,1)*DF_V(1)+AT_MEMORY(ID,l,4,3)*DF_V(3)&
            +AT_MEMORY(ID,l,4,4)*DF_V(4)+AT_MEMORY(ID,l,4,6)*DF_V(6)
     
     ETA(5)=AT_MEMORY(ID,l,5,2)*DF_V(2)+AT_MEMORY(ID,l,5,5)*DF_V(5)
     
     ETA(6)=AT_MEMORY(ID,l,6,1)*DF_V(1)+ AT_MEMORY(ID,l,6,3)*DF_V(3) &
            + AT_MEMORY(ID,l,6,4)*DF_V(4)+AT_MEMORY(ID,l,6,6)*DF_V(6)
     
     Var_Memory(id,l,1:6)=TSRC_ABC(ID,l,1)*Var_Memory(id,l,1:6)+TSRC_ABC(ID,l,2)*ETA(1:6) &
                            +TSRC_ABC(ID,l,3)*Var_Memory_T(ID,l,1:6)
     
     Var_Memory_T(ID,l,1:6)=ETA(1:6)                       
    enddo
    
   ! V_kl(ID,1:6)=DF_V(1:6)
    
return
end subroutine

!-----------------
Subroutine TSRC_Matrix(ID,IREC,NPT,NSLS,dt,PD,tau_sig,tau_ep,TSRC_ABC, N_order)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    !-Input
    dimension TAU_SIG(NPT,NSLS,21),TAU_ep(NPT,NSLS,21),PD(22)
    DIMENSION TSRC_ABC(NPT,NSLS,4)
    
    TSRC_ABC(ID,1:NSLS,1:4)=0.0
    
    if(N_order==1)then  !1st order
    TSRC_ABC(ID,1:NSLS,1)=exp(-dt/tau_sig(ID,1:NSLS,1))
    TSRC_ABC(ID,1:NSLS,2)=tau_sig(ID,1:NSLS,1)+tau_sig(ID,1:NSLS,1)**2/dt*&
                            (exp(-dt/tau_sig(ID,1:NSLS,1))-1.0)
    TSRC_ABC(ID,1:NSLS,3)=tau_sig(ID,1:NSLS,1)**2/dt*(1.0-exp(-dt/tau_sig(ID,1:NSLS,1)))&
                            -tau_sig(ID,1:NSLS,1)*exp(-dt/tau_sig(ID,1:NSLS,1))
    elseif(N_order==2)then  !2nd order
    TSRC_ABC(ID,1:NSLS,1)=exp(-dt/tau_sig(ID,1:NSLS,1))
    
    TSRC_ABC(ID,1:NSLS,2)=-tau_sig(ID,1:NSLS,1)*(exp(-dt/tau_sig(ID,1:NSLS,1))-1.0)
    
    TSRC_ABC(ID,1:NSLS,3)=(-tau_sig(ID,1:NSLS,1))**2&
                    *(exp(-dt/tau_sig(ID,1:NSLS,1))-1.0+dt/tau_sig(ID,1:NSLS,1))/(2.0*dt)
                    
    TSRC_ABC(ID,1:NSLS,4)=(-tau_sig(ID,1:NSLS,1))**3&
                    *(exp(-dt/tau_sig(ID,1:NSLS,1))-1.0+dt/tau_sig(ID,1:NSLS,1) &
                    -(dt/tau_sig(ID,1:NSLS,1))**2/2.0)/(dt**2)
                            
    endif
return
end subroutine

!-------------------------------
!-Viscoacoustic
Subroutine Memory_var_TSRC_WS_T(ID,ITX,IREC,NNX,NPT,NSLS,wky,dt,ZP,DX_Y,DY_Y,DZ_Y,&
DX_Y_WS,DY_Y_WS,DZ_Y_WS,NORDX,IDX,NORDZ,IDZ,DFS,PD,tau_sig,tau_ep,TSRC_ABC,V_kl,Var_Memory)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    !-Input
    Dimension ZP(3,NPT)
    Dimension DX_Y(13),DY_Y(13),DZ_Y(13),DX_Y_WS(13),DY_Y_WS(13),DZ_Y_WS(13)
    dimension IDX(NPT,NORDX),IDZ(NPT,NORDZ),DFS(NPT,2,NORDX,4)
    dimension TAU_SIG(NPT,NSLS,21),TAU_ep(NPT,NSLS,21), PD(22,NPT), TSRC_ABC(NPT,NSLS,4)
    dimension V_kl(NNX,6)
    dimension Var_Memory(NNX,Nsls,7)
    
    !-In subroutine
    dimension DF_V(3)
    
    fac=1.d0+sum(TAU_EP(ID+1,1:NSLS,1)/TAU_SIG(ID+1,1:NSLS,1)-1.d0)/dble(NSLS)
    Bulk_Me=PD(2,ID+1)/fac;
    !if(id==57085)print *, sqrt(Bulk_Me/2000.d0)
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
    
    !DyVy
    DF_V(2)=DY_Y_WS(2)
    
    !DzVz
    DF_V(3)=DZ_Y(7)+ZP(2,ID)*(DZ_Y_WS(1)-DZ_Y(5))
    
    !-I
    Var_Memory(ITX,1:NSLS,2)=TSRC_ABC(ID+1,1:NSLS,1)*Var_Memory(ITX,1:NSLS,2) &
        +TSRC_ABC(ID+1,1:NSLS,2)*sum(DF_V)+TSRC_ABC(ID+1,1:NSLS,3)*V_kl(ITX,1)
    !Recursive conv
    Var_Memory(ITX,1:NSLS,1)=Bulk_Me/(dble(NSLS)*tau_sig(ID+1,1:NSLS,1))*(1.d0-tau_ep(ID+1,1:NSLS,1)/tau_sig(ID+1,1:NSLS,1)) &
    *Var_Memory(ITX,1:NSLS,2)
    
    V_kl(ITX,1)=sum(DF_V)
    
return
end subroutine
!-------------------------------
!-Viscoacoustic
Subroutine Memory_var_TSRC_T(ID,IREC,NPT,NSLS,wky,dt,DX_Y,DY_Y,DZ_Y,NORDX,IDX,NORDZ,&
IDZ,DFS,PD,tau_sig,tau_ep,TSRC_ABC,V_kl,Var_Memory)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    !-Input
    Dimension DX_Y(13),DY_Y(13),DZ_Y(13)
    dimension IDX(NPT,NORDX),IDZ(NPT,NORDZ),DFS(NPT,2,NORDX,4)
    dimension TAU_SIG(NPT,NSLS,21),TAU_ep(NPT,NSLS,21), PD(22,NPT), TSRC_ABC(NPT,NSLS,3)
    dimension V_kl(NPT,6)
    dimension Var_Memory(NPT,Nsls,7)
    
    !-In subroutine
    dimension DF_V(3)
    
    fac=1.d0+sum(TAU_EP(ID,1:NSLS,1)/TAU_SIG(ID,1:NSLS,1)-1.d0)/dble(NSLS)
    Bulk_Me=PD(2,ID)/fac;
    !if(id==57085)print *, sqrt(Bulk_Me/2000.d0)
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
    
    
    !-I
    Var_Memory(ID,1:NSLS,2)=TSRC_ABC(ID,1:NSLS,1)*Var_Memory(ID,1:NSLS,2) &
        +TSRC_ABC(ID,1:NSLS,2)*sum(DF_V)+TSRC_ABC(ID,1:NSLS,3)*V_kl(ID,1)
    !Recursive conv
    Var_Memory(ID,1:NSLS,1)=Bulk_Me/(dble(NSLS)*tau_sig(ID,1:NSLS,1))*(1.d0-tau_ep(ID,1:NSLS,1)/tau_sig(ID,1:NSLS,1)) &
    *Var_Memory(ID,1:NSLS,2)
    
    V_kl(ID,1)=sum(DF_V)
    
return
end subroutine

