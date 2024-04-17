!----------------------------
Subroutine FS_T(ID,NPT,ZP,NNX,ITX,NORDX,IDXG,DEXG,NORDZ,IDZG,DEZG,DFSG,&
                    DX_Y,DY_Y,DZ_Y,PM,V,ST1,ST2,V_T,ST_T)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    !-Input
    Dimension ZP(3,NPT), PM(22,NPT)
    Dimension IDXG(NPT,NORDX),IDZG(NPT,NORDZ),DEXG(NPT,NORDX),DEZG(NPT,NORDZ),DFSG(NPT,2,NORDX,4)
    Dimension V(3,NPT),ST1(3,NPT),ST2(3,NPT)
    dimension DX_Y(13),DY_Y(13),DZ_Y(13)
    dimension V_T(3,5,NPT), ST_T(6,5,NPT)
    
    DIMENSION R(6,6), ST_R(6), ST_TEMP(6),R_T(6,6)
    
    IX=0; 
    do ii=1,NORDX
        if(ID==IDXG(ID,II))IX=ii
    enddo
    
    if(IX==0)then
        print *, 'Memory_var_T: Point not searched', ID
        print *, IDXG(ID,1:NORDX)
        stop
    endif
    
    !V_T(1,1,ID)=(DX_Y(8)+ZP(2,ID)*(DZ_Y(8))+DY_Y(9))/PM(1,ID)
    !V_T(3,1,ID)=(ZP(3,ID)*(ST1(1,ID)+ZP(2,ID)*DX_Y(8))&
    !            +ZP(2,ID)*DY_Y(9)+ZP(2,ID)**2*DZ_Y(8))/PM(1,ID)
    !V_T(3,1,ID)=(DX_Y(11)&
    !            +ZP(2,ID)*DY_Y(9)+ZP(2,ID)**2*DZ_Y(8))/PM(1,ID)
    !V_T(3,1,ID)=(ZP(3,ID)*ST1(1,ID)+ZP(2,ID)*DX_Y(8)&
    !            +ZP(2,ID)*DY_Y(9)+DZ_Y(13))/PM(1,ID)
    
    ST_T(4,1,ID)=ZP(2,ID)*(ST_T(1,1,ID))
    ST_T(5,1,ID)=ZP(2,ID)*ST_T(2,1,ID)
    ST_T(6,1,ID)=ZP(2,ID)**2*ST_T(1,1,ID)
    
    !Rotation for BC satisfactory
    ! 1. Normal velocity continuity
    S=sqrt(DFSG(ID,1,IX,3)**2+DFSG(ID,1,IX,4)**2)
    
    !VT=DFSG(ID,1,IX,4)/s*V_T(1,1,ID)+DFSG(ID,1,IX,3)/s*V_T(3,1,ID)
    !VN=-DFSG(ID,1,IX,3)/s*V_T(1,1,ID)+DFSG(ID,1,IX,4)/s*V_T(3,1,ID)
    !VN=0.d0
    
    !V_T(1,1,ID)=DFSG(ID,1,IX,4)/s*VT-DFSG(ID,1,IX,3)/s*VN
    !V_T(3,1,ID)=DFSG(ID,1,IX,3)/s*VT+DFSG(ID,1,IX,4)/s*VN
    
    ! 2. Normal Stress
    SIG_T1T2=(DFSG(ID,1,IX,4)/s)*ST_T(2,1,ID)-(DFSG(ID,1,IX,3)/s)*ST_T(5,1,ID)
    SIG_TT=(DFSG(ID,1,IX,4)/s)**2*ST_T(1,1,ID)-2.0*(DFSG(ID,1,IX,4)/s)*(DFSG(ID,1,IX,3)/s)*ST_T(4,1,ID) &
        +(DFSG(ID,1,IX,3)/s)**2*ST_T(6,1,ID)
    SIG_NN=(DFSG(ID,1,IX,3)/s)**2*ST_T(1,1,ID)+2.0*(DFSG(ID,1,IX,4)/s)*(DFSG(ID,1,IX,3)/s)*ST_T(4,1,ID) &
        +(DFSG(ID,1,IX,4)/s)**2*ST_T(6,1,ID)
    SIG_NT=(DFSG(ID,1,IX,3)/s)*(DFSG(ID,1,IX,4)/s)*ST_T(1,1,ID)&
        +(DFSG(ID,1,IX,4)/s**2-DFSG(ID,1,IX,3)/s**2)*ST_T(4,1,ID) &
        -(DFSG(ID,1,IX,4)/s)*(DFSG(ID,1,IX,3)/s)*ST_T(6,1,ID)
        
    SIG_NN=0.0
    SIG_NT=0.0
    
    !-reverse
    ST_T(1,1,ID)=(DFSG(ID,1,IX,4)/s)**2*SIG_TT&
               +2.0*(DFSG(ID,1,IX,3)/s)*(DFSG(ID,1,IX,4)/s)*SIG_NT&
                +(DFSG(ID,1,IX,3)/s)**2*SIG_NN
    ST_T(2,1,ID)=(DFSG(ID,1,IX,4)/s)*SIG_T1T2
    ST_T(4,1,ID)=-DFSG(ID,1,IX,4)/s*DFSG(ID,1,IX,3)/s*SIG_TT&
                 +(DFSG(ID,1,IX,4)/s**2-DFSG(ID,1,IX,3)/s**2)*SIG_NT&
                 +DFSG(ID,1,IX,4)/s*DFSG(ID,1,IX,3)/s*SIG_NN
    ST_T(5,1,ID)=-(DFSG(ID,1,IX,3)/s)*SIG_T1T2
    ST_T(6,1,ID)=(DFSG(ID,1,IX,3)/s)**2*SIG_TT&
               -2.0*(DFSG(ID,1,IX,3)/s)*(DFSG(ID,1,IX,4)/s)*SIG_NT&
                +(DFSG(ID,1,IX,4)/s)**2*SIG_NN
    
    !ST_T(2,1,ID)=0.0;   
    !ST_T(5,1,ID)=0.0
    
return
end subroutine

!----------------------------
Subroutine WS_T(ID,NPT,ZP,NNX,ITX,NORDX,IDXG,DEXG,NORDZ,IDZG,DEZG,DFSG,&
                    IDZ_WS,DEZ_WS,DFS_WS,&
                    DX_Y,DY_Y,DZ_Y,DX_Y_WS,DY_Y_WS,DZ_Y_WS,&
                    PM,WV,P,V,ST1,ST2,WV_T,P_T,V_T,ST_T)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    !-Input
    Dimension ZP(3,NPT), PM(22,NPT)
    Dimension IDXG(NPT,NORDX),IDZG(NPT,NORDZ),DEXG(NPT,NORDX),DEZG(NPT,NORDZ),DFSG(NPT,2,NORDX,4)
    Dimension IDZ_WS(NNX,NORDZ),DEZ_WS(NNX,NORDZ),DFS_WS(NNX,2,NORDX,4)
    Dimension WV(3,NPT),P(NPT),V(3,NPT),ST1(3,NPT),ST2(3,NPT)
    dimension DX_Y(13),DY_Y(13),DZ_Y(13),DX_Y_WS(13),DY_Y_WS(13),DZ_Y_WS(13)
    dimension  WV_T(3,5,NPT), P_T(5,NPT),V_T(3,5,NPT), ST_T(6,5,NPT)
    
    DIMENSION R(6,6), ST_R(6), ST_TEMP(6),R_T(6,6), PD(22)
    
    IX=0; 
    do ii=1,NORDX
        if(ID==IDXG(ID,II))IX=ii
    enddo
    
    if(IX==0)then
        print *, 'Memory_var_T: Point not searched', ID
        print *, IDXG(ID,1:NORDX)
        stop
    endif
    
    PD(1:22)=PM(1:22,ID+1)
    PD(8)=PM(8,ID+1); PD(13)=PM(13,ID+1)
    PD(1:22)=(PD(1:22)+PM(1:22,ID))/2.0
                
    !ST_T(2,1,ID)=PM(7,ID)*DX_Y(5)+PM(22,ID)*DY_Y(5)+PM(21,ID)*DZ_Y(5)&
    !            +PM(22,ID)*DX_Y(6)+PM(12,ID)*DY_Y(6)+PM(19,ID)*DZ_Y(6)&
    !            +PM(21,ID)*(ZP(2,ID)*(DX_Y(5)-DX_Y_WS(1))+DX_Y_WS(3)+ZP(3,ID)*(V(1,ID)-WV(1,ID)))&
    !            +PM(19,ID)*(ZP(2,ID)*(DY_Y(5)-DY_Y_WS(1))+DY_Y_WS(3))&
    !            +PM(16,ID)*(ZP(2,ID)*(DZ_Y(5)-DZ_Y_WS(1))+DZ_Y_WS(3))
    
    !ST_T(3,1,ID)=PM(3,ID)*DX_Y(5)+PM(12,ID)*DY_Y(5)+PM(11,ID)*DZ_Y(5)&
    !            +PM(12,ID)*DX_Y(6)+PM(8,ID)*DY_Y(6)+PM(10,ID)*DZ_Y(6)&
    !            +PM(11,ID)*(ZP(2,ID)*(DX_Y(5)-DX_Y_WS(1))+DX_Y_WS(3)+ZP(3,ID)*(V(1,ID)-WV(1,ID)))&
    !            +PM(10,ID)*(ZP(2,ID)*(DY_Y(5)-DY_Y_WS(1))+DY_Y_WS(3))&
    !            +PM(9,ID)*(ZP(2,ID)*(DZ_Y(5)-DZ_Y_WS(1))+DZ_Y_WS(3))
    
    !ST_T(2,1,ID)=PD(7)*DX_Y(5)+PD(22)*DY_Y(5)+PD(21)*DZ_Y(5)&
    !            +PD(22)*DX_Y(6)+PD(12)*DY_Y(6)+PD(19)*DZ_Y(6)&
    !            +PD(21)*(ZP(2,ID)*(DX_Y(5)-DX_Y_WS(1))+DX_Y_WS(3)+ZP(3,ID)*(V(1,ID)-WV(1,ID)))&
    !            +PD(19)*(ZP(2,ID)*(DY_Y(5)-DY_Y_WS(1))+DY_Y_WS(3))&
    !            +PD(16)*(ZP(2,ID)*(DZ_Y(5)-DZ_Y_WS(1))+DZ_Y_WS(3))
    
    !ST_T(3,1,ID)=PD(3)*DX_Y(5)+PD(12)*DY_Y(5)+PD(11)*DZ_Y(5)&
    !            +PD(12)*DX_Y(6)+PD(8)*DY_Y(6)+PD(10)*DZ_Y(6)&
    !            +PD(11)*(ZP(2,ID)*(DX_Y(5)-DX_Y_WS(1))+DX_Y_WS(3)+ZP(3,ID)*(V(1,ID)-WV(1,ID)))&
    !            +PD(10)*(ZP(2,ID)*(DY_Y(5)-DY_Y_WS(1))+DY_Y_WS(3))&
    !            +PD(9)*(ZP(2,ID)*(DZ_Y(5)-DZ_Y_WS(1))+DZ_Y_WS(3))
    
    !V_T(1,1,ID)=(DX_Y(8)+ZP(2,ID)*(DZ_Y_WS(4)+DZ_Y(8))+DY_Y(9))/PM(1,ID+1)
    !V_T(2,1,ID)=(DX_Y(9)+DY_Y(10)+ZP(2,ID)*DZ_Y(9))/(PM(1,ID))
    V_T(3,1,ID)=(ZP(2,ID)*(DX_Y_WS(4)+DX_Y(8))+ZP(3,ID)*(P(ID)+ST1(1,ID))+(ZP(2,ID)**2-1.0)*DZ_Y_WS(4)&
                +ZP(2,ID)**2*DZ_Y(8)+ZP(2,ID)*DY_Y(9))/(PM(1,ID+1))
    
    !P_T(1,ID)=P_T(1,ID)-((PM(3,ID)+2.0/3.0*PM(22,ID)+PM(4,ID)+2.0/3.0*PM(20,ID))/2.0)&
    !            *(DX_Y_WS(1)+ZP(2,ID)*(DZ_Y_WS(1)-DZ_Y(5))+DZ_Y(7)+DY_Y_WS(2))
    !P_T(1,ID)=P_T(1,ID)-(PM(2,ID+1))&
    !            *(DX_Y_WS(1)+ZP(2,ID)*(DZ_Y_WS(1)-DZ_Y(5))+DZ_Y(7)+DY_Y_WS(2))
    P_T(1,ID)=P_T(1,ID)-((PM(2,ID)+PM(8,ID)+PM(13,ID))/3.0)&
               *(DX_Y_WS(1)+ZP(2,ID)*(DZ_Y_WS(1)-DZ_Y(5))+DZ_Y(7)+DY_Y_WS(2))            
    
    ST_T(4,1,ID)=ZP(2,ID)*(P_T(1,ID)+ST_T(1,1,ID))
    ST_T(5,1,ID)=ZP(2,ID)*ST_T(2,1,ID)
    ST_T(6,1,ID)=(ZP(2,ID)**2-1.0)*P_T(1,ID)+ZP(2,ID)**2*ST_T(1,1,ID)
    WV_T(3,1,ID)=ZP(2,ID)*(WV_T(1,1,ID)-V_T(1,1,ID))+V_T(3,1,ID)
    
    !Rotation for BC satisfactory
    ! 1. Normal velocity continuity
    S=sqrt(DFSG(ID,1,IX,3)**2+DFSG(ID,1,IX,4)**2)
    vt_w=DFSG(ID,1,IX,4)/s*WV_T(1,1,ID)+DFSG(ID,1,IX,3)/s*WV_T(3,1,ID)
    vn_w=-DFSG(ID,1,IX,3)/s*WV_T(1,1,ID)+DFSG(ID,1,IX,4)/s*WV_T(3,1,ID)
    
    vt_s=DFSG(ID,1,IX,4)/s*V_T(1,1,ID)+DFSG(ID,1,IX,3)/s*V_T(3,1,ID)
    vn_s=-DFSG(ID,1,IX,3)/s*V_T(1,1,ID)+DFSG(ID,1,IX,4)/s*V_T(3,1,ID)
    
    VN=(vn_w+vn_s)/2.0
    WV_T(1,1,ID)=DFSG(ID,1,IX,4)/s*VT_W-DFSG(ID,1,IX,3)/s*VN
    !WV_T(3,1,ID)=DFSG(ID,1,IX,3)/s*VT_W+DFSG(ID,1,IX,4)/s*VN
    
    !VN_S=(vn_w+vn_s)/2.0
    !V_T(1,1,ID)=DFSG(ID,1,IX,4)/s*VT_S-DFSG(ID,1,IX,3)/s*VN
    !V_T(3,1,ID)=DFSG(ID,1,IX,3)/s*VT_S+DFSG(ID,1,IX,4)/s*VN
    
    ! 2. Normal Stress
    SIG_T1T2=(DFSG(ID,1,IX,4)/s)*ST_T(2,1,ID)-(DFSG(ID,1,IX,3)/s)*ST_T(5,1,ID)
        
    SIG_TT=(DFSG(ID,1,IX,4)/s)**2*ST_T(1,1,ID)-2.0*(DFSG(ID,1,IX,4)/s)*(DFSG(ID,1,IX,3)/s)*ST_T(4,1,ID) &
        +(DFSG(ID,1,IX,3)/s)**2*ST_T(6,1,ID)
    SIG_NN=(DFSG(ID,1,IX,3)/s)**2*ST_T(1,1,ID)+2.0*(DFSG(ID,1,IX,4)/s)*(DFSG(ID,1,IX,3)/s)*ST_T(4,1,ID) &
        +(DFSG(ID,1,IX,4)/s)**2*ST_T(6,1,ID)
    SIG_NT=(DFSG(ID,1,IX,3)/s)*(DFSG(ID,1,IX,4)/s)*ST_T(1,1,ID)&
        +(DFSG(ID,1,IX,4)/s**2-DFSG(ID,1,IX,3)/s**2)*ST_T(4,1,ID) &
        -(DFSG(ID,1,IX,4)/s)*(DFSG(ID,1,IX,3)/s)*ST_T(6,1,ID)
        
    SIG_NN=(SIG_NN-P_T(1,ID))/2.0;
    !SIG_TT=(SIG_TT+SIG_NN)/2.0;
    SIG_NT=0.0
    
    !-reverse
    ST_T(1,1,ID)=(DFSG(ID,1,IX,4)/s)**2*SIG_TT&
               +2.0*(DFSG(ID,1,IX,3)/s)*(DFSG(ID,1,IX,4)/s)*SIG_NT&
                +(DFSG(ID,1,IX,3)/s)**2*SIG_NN
    ST_T(2,1,ID)=(DFSG(ID,1,IX,4)/s)*SIG_T1T2
    ST_T(4,1,ID)=-DFSG(ID,1,IX,4)/s*DFSG(ID,1,IX,3)/s*SIG_TT&
                 +(DFSG(ID,1,IX,4)/s**2-DFSG(ID,1,IX,3)/s**2)*SIG_NT&
                 +DFSG(ID,1,IX,4)/s*DFSG(ID,1,IX,3)/s*SIG_NN
    ST_T(5,1,ID)=-(DFSG(ID,1,IX,3)/s)*SIG_T1T2
    
    ST_T(6,1,ID)=(DFSG(ID,1,IX,3)/s)**2*SIG_TT&
               -2.0*(DFSG(ID,1,IX,3)/s)*(DFSG(ID,1,IX,4)/s)*SIG_NT&
                +(DFSG(ID,1,IX,4)/s)**2*SIG_NN
    
return
end subroutine

!----------------------------
Subroutine ELASTIC_T2(ID,NPT,V,ST1,ST2,DX_Y,DY_Y,DZ_Y,V_A,ST_A,V_T,ST_T)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    !-Input
    Dimension V(3,NPT),ST1(3,NPT),ST2(3,NPT)
    dimension DX_Y(13),DY_Y(13),DZ_Y(13)
    dimension V_A(NPT,6,6), ST_A(NPT,3,3,6)
    dimension  V_T(3,5,NPT), ST_T(6,5,NPT)
    
    !-In subroutine
    dimension  ST_Y(6),ST_DX0(6),ST_DZ0(6),DF_V(6)
    
    !-Assign parameters point
    ST_Y(1:3)=DY_Y(8:10)
    ST_Y(4:6)=DY_Y(11:13)
    
    !-Composite DxST, DzST
    ST_DX0(1:3)=DX_Y(8:10)
    ST_DZ0(1:3)=DZ_Y(8:10)
    
    ST_DX0(4:6)=DX_Y(11:13)
    ST_DZ0(4:6)=DZ_Y(11:13)
    
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
    
    !-Matrix multiplication -> V_t, St_t
    !V_T(1:3,1,ID)=matmul(ST_A(ID,1,1:3,1:6),ST_DX0(1:6))+matmul(ST_A(ID,2,1:3,1:6),ST_DZ0(1:6)) &
    !        +matmul(ST_A(ID,3,1:3,1:6),ST_Y(1:6))
    !ST_T(1:6,1,ID)=matmul(V_A(ID,1,1:6,1:3),V_DX0(1:3))+matmul(V_A(ID,2,1:6,1:3),V_DZ0(1:3))&
    !        +matmul(V_A(ID,3,1:6,1:3),V_Y(1:3))
    
    V_T(1,1,ID)=ST_A(ID,1,1,1)*ST_DX0(1)&
            +ST_A(ID,2,1,4)*ST_DZ0(4)+ST_A(ID,3,1,2)*ST_Y(2)
    
    V_T(2,1,ID)=ST_A(ID,1,2,2)*ST_DX0(2)&
            +ST_A(ID,2,2,5)*ST_DZ0(5)+ST_A(ID,3,2,3)*ST_Y(3)
            
    V_T(3,1,ID)=ST_A(ID,1,3,4)*ST_DX0(4)&
            +ST_A(ID,2,3,6)*ST_DZ0(6)+ST_A(ID,3,3,5)*ST_Y(5)
    
    !ST_T(1:6,1,ID)=matmul(V_A(ID,1:6,1:6),DF_V(1:6))
    
    ST_T(1,1,ID)=V_A(ID,1,1)*DF_V(1)+V_A(ID,1,3)*DF_V(3)&
                +V_A(ID,1,4)*DF_V(4)+V_A(ID,1,6)*DF_V(6)
                
    ST_T(2,1,ID)=V_A(ID,2,2)*DF_V(2)+V_A(ID,2,5)*DF_V(5)
    
    ST_T(3,1,ID)=V_A(ID,3,1)*DF_V(1)+V_A(ID,3,3)*DF_V(3)&
                +V_A(ID,3,4)*DF_V(4)+V_A(ID,3,6)*DF_V(6)
                
    ST_T(4,1,ID)=V_A(ID,4,1)*DF_V(1)+V_A(ID,4,3)*DF_V(3)&
                +V_A(ID,4,4)*DF_V(4)+V_A(ID,4,6)*DF_V(6)
    
    ST_T(5,1,ID)=V_A(ID,5,2)*DF_V(2)+V_A(ID,5,5)*DF_V(5)
    
    ST_T(6,1,ID)=V_A(ID,6,1)*DF_V(1)+V_A(ID,6,3)*DF_V(3)&
                +V_A(ID,6,4)*DF_V(4)+V_A(ID,6,6)*DF_V(6)
                
return
end subroutine

!-----------------------------------------------
Subroutine ELASTIC_A(ID,NPT,NORDX,IX,IZ,WKy,DFS,PD,A,B)

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    DIMENSION PD(22),DFS(2,NORDX,4)
    DIMENSION A(NPT,3,3,6),B(NPT,6,6)
        
!-Init
A(ID,:,:,:)=0.d0; B(ID,:,:)=0.d0

!-for Velocity time derivatives
    !A(ID,1,1,1)=DFS(1,IX,2); A(ID,1,2,2)=DFS(1,IX,2); A(ID,1,3,4)=DFS(1,IX,2); 
    
    !A(ID,2,1,1)=DFS(2,IZ,3); A(ID,2,1,4)=DFS(2,IZ,4);
    !A(ID,2,2,2)=DFS(2,IZ,3); A(ID,2,2,5)=DFS(2,IZ,4);
    !A(ID,2,3,4)=DFS(2,IZ,3); A(ID,2,3,6)=DFS(2,IZ,4);
    
    A(ID,1,1,1)=1; A(ID,1,2,2)=1; A(ID,1,3,4)=1; 
    
    A(ID,2,1,4)=1;
    A(ID,2,2,5)=1;
    A(ID,2,3,6)=1;
    
    A(ID,3,1,2)=1.d0;A(ID,3,2,3)=1.d0;A(ID,3,3,5)=1.d0
    
    A(ID,:,:,:)=A(ID,:,:,:)/PD(1)
    
!-For Stress time derivatives
    
    B(ID,1,1)=PD(2);B(ID,1,2)=2.0*PD(7);B(ID,1,3)=PD(3)
    B(ID,1,4)=2.0*PD(6);B(ID,1,5)=2.0*PD(5);B(ID,1,6)=PD(4)
    
    B(ID,2,1)=PD(7);B(ID,2,2)=2.0*PD(22);B(ID,2,3)=PD(12)
    B(ID,2,4)=2.0*PD(21);B(ID,2,5)=2.0*PD(19);B(ID,2,6)=PD(16)
    
    B(ID,3,1)=PD(3);B(ID,3,2)=2.0*PD(12);B(ID,3,3)=PD(8)
    B(ID,3,4)=2.0*PD(11);B(ID,3,5)=2.0*PD(10);B(ID,3,6)=PD(9)
    
    B(ID,4,1)=PD(6);B(ID,4,2)=2.0*PD(21);B(ID,4,3)=PD(11)
    B(ID,4,4)=2.0*PD(20);B(ID,4,5)=2.0*PD(18);B(ID,4,6)=PD(15)
    
    B(ID,5,1)=PD(5);B(ID,5,2)=2.0*PD(19);B(ID,5,3)=PD(10)
    B(ID,5,4)=2.0*PD(18);B(ID,5,5)=2.0*PD(17);B(ID,5,6)=PD(14)
    
    B(ID,6,1)=PD(4);B(ID,6,2)=2.0*PD(16);B(ID,6,3)=PD(9)
    B(ID,6,4)=2.0*PD(15);B(ID,6,5)=2.0*PD(14);B(ID,6,6)=PD(13)
    
!-----------------
return
end subroutine


!------------------------------------------------------------
Subroutine Acoustic_T2(ID,NPT,P,WV,DX_Y,DY_Y,DZ_Y,P_A,P_T,WV_T)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    !-Input
    dimension P(NPT),P_T(5,NPT)
    dimension WV(3,NPT),WV_T(3,5,NPT)
    dimension P_A(NPT,3,4,4)
    dimension DX_Y(13),DY_Y(13),DZ_Y(13)
    
    !-----
    dimension Y(4),YX(4),YZ(4),YY(4)
    
    !-Matrix multiplication -> V_t, St_t
    YX(1:4)=0.d0; YZ(1:4)=0.d0; YY(1:4)=0.d0
    
    YX(1)=DX_Y(4); YX(2:4)=DX_Y(1:3)
    YZ(1)=DZ_Y(4); YZ(2:4)=DZ_Y(1:3)
    YY(1)=DY_Y(4); YY(2:4)=DY_Y(1:3)
    
    
    !Y=matmul(P_A(ID,1,:,:),YX)+matmul(P_A(ID,2,:,:),YZ) &
    !        +matmul(P_A(ID,3,:,:),YY)
    
    P_T(1,ID)=P_A(ID,1,1,2)*YX(2)+P_A(ID,2,1,2)*YZ(2)+P_A(ID,2,1,4)*YZ(4)&
                +P_A(ID,3,1,3)*YY(3)
    
    WV_T(1,1,ID)=P_A(ID,1,2,1)*YX(1)+P_A(ID,2,2,1)*YZ(1)
    WV_T(2,1,ID)=P_A(ID,3,3,1)*YY(1)
    WV_T(3,1,ID)=P_A(ID,2,4,1)*YZ(1)
    
    !P_T(1,ID)=Y(1)
    !WV_T(1:3,1,ID)=Y(2:4)
return
end subroutine


!-----------
Subroutine Acoustic_A(ID,NPT,NORDX,IX,IZ,WKy,DFS,PD,A)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION PD(22),DFS(2,NORDX,4)
    DIMENSION A(NPT,3,4,4)
    
    A(ID,:,:,:)=0.d0
    A(ID,1,1,2)=-PD(2)
    A(ID,1,2,1)=-1.d0/PD(1)
    
    A(ID,2,1,4)=-PD(2)
    A(ID,2,4,1)=-1.d0/PD(1)
    
    
    A(ID,3,1,3)=-PD(2);
    A(ID,3,3,1)=-1.d0/PD(1);
    
return
end subroutine

!----------------------
!-Choose Y
!----------------------
Subroutine Assign_Y_WS(ID,NPT,NNX,ITX,NORDX,IDX,DEX,NORDZ,IDZ,DEZ,DFS,wky,WV,P,V,ST1,ST2,DX_Y,DY_Y,DZ_Y)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    dimension IDX(NPT,NORDX),IDZ(NNX,NORDZ),DEX(NPT,NORDX),DEZ(NNX,NORDZ),DFS(NNX,2,NORDX,4)
    Dimension WV(3,NPT),V(3,NPT),ST1(3,NPT),ST2(3,NPT)
    dimension P(NPT)
    dimension DX_Y(13),DY_Y(13),DZ_Y(13)
    
    !===================
    IX=0; IZ=0
    do ii=1,NORDX
        if(ID==IDX(ID,II))IX=ii
    enddo
    
    do jj=1,NORDZ
        if(ID==IDZ(ITX,jj))IZ=jj
    enddo
    
    if(IX==0.or.IZ==0)then
        print *, 'Memory_var_T: Point not searched', ID
        print *, IDX(ID,1:NORDX)
        print *, IDZ(ITX,1:NORDZ)
        stop
    endif
    
    DX_Y(1:13)=0.d0; DY_Y(1:13)=0.d0; DZ_Y(1:13)=0.d0
    !==========================================
    !-for W
    
    do ii=1,3
    DX_Y(ii)=dot_product(DEX(ID,1:NORDX),WV(ii,IDX(ID,1:NORDX)))
    !DX_Y(ii)=DX_Y(ii)+DFS(2,IZ,3)*dot_product(DEZ(1:NORDZ),WV(ii,IDZ(1:NORDZ)))
    enddo
    
    DY_Y(2)=-wky*WV(2,ID)
    
    do jj=1,3
    DZ_Y(jj)=dot_product(DEZ(ITX,1:NORDZ),WV(jj,IDZ(ITX,1:NORDZ)))
    enddo
    
    DX_Y(4)=dot_product(DEX(ID,1:NORDX),P(IDX(ID,1:NORDX)))
    !DX_Y(4)=DX_Y(4)+DFS(2,IZ,3)*dot_product(DEZ(1:NORDZ),P(IDZ(1:NORDZ)))
    
    DY_Y(4)=wky*P(ID)
    
    DZ_Y(4)=dot_product(DEZ(ITX,1:NORDZ),P(IDZ(ITX,1:NORDZ)))
   
    
    !-------------
    !---Only Solid
   
    !velocity
    do ii=5,7
        DX_Y(ii)=dot_product(DEX(ID,1:NORDX),V(ii-4,IDX(ID,1:NORDX)))
    enddo
    
    DY_Y(5:7)=wky*V(1:3,ID)
    DY_Y(6)=-DY_Y(6)
    
    do jj=5,7
        DZ_Y(jj)=dot_product(DEZ(ITX,1:NORDZ),V(jj-4,IDZ(ITX,1:NORDZ)))
    enddo 
    
    !Stress 1
    do ii=8,10
    DX_Y(ii)=dot_product(DEX(ID,1:NORDX),ST1(ii-7,IDX(ID,1:NORDX)))
    enddo
    
    DY_Y(8:10)=wky*ST1(1:3,ID)
    DY_Y(9)=-DY_Y(9)
    
    do jj=8,10
    DZ_Y(jj)=dot_product(DEZ(ITX,1:NORDZ),ST1(jj-7,IDZ(ITX,1:NORDZ)))
    enddo
    
    !Stress 2
    do ii=11,13
        DX_Y(ii)=dot_product(DEX(ID,1:NORDX),ST2(ii-10,IDX(ID,1:NORDX)))
    enddo
    
    DY_Y(11:13)=wky*ST2(1:3,ID)
    DY_Y(12)=-DY_Y(12)
    
    do jj=11,13
    DZ_Y(jj)=dot_product(DEZ(ITX,1:NORDZ),ST2(jj-10,IDZ(ITX,1:NORDZ)))
    enddo
    
   
    
    DX_Y(1:13)=DFS(ITX,1,IX,2)*DX_Y(1:13)+DFS(ITX,1,IX,3)*DZ_Y(1:13)
    DZ_Y(1:13)=DFS(ITX,2,IZ,4)*DZ_Y(1:13)
    
return
end subroutine
!----------------------
!-Choose Y
!----------------------
Subroutine Assign_Y(ID,MNO,NPT,NORDX,IDX,DEX,NORDZ,IDZ,DEZ,DFS,wky,WV,P,V,ST1,ST2,DX_Y,DY_Y,DZ_Y)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    dimension MNO(NPT),IDX(NPT,NORDX),IDZ(NPT,NORDZ),DEX(NPT,NORDX),DEZ(NPT,NORDZ),DFS(NPT,2,NORDX,4)
    Dimension WV(3,NPT),V(3,NPT),ST1(3,NPT),ST2(3,NPT)
    dimension P(NPT)
    dimension DX_Y(13),DY_Y(13),DZ_Y(13)
    
    !===================
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
    
    DX_Y(1:13)=0.d0; DY_Y(1:13)=0.d0; DZ_Y(1:13)=0.d0
    !==========================================
    !-for W
    if(MNO(ID)==3.or.MNO(ID)==2)then
    do ii=1,3
    DX_Y(ii)=dot_product(DEX(ID,1:NORDX),WV(ii,IDX(ID,1:NORDX)))
    !DX_Y(ii)=DX_Y(ii)+DFS(2,IZ,3)*dot_product(DEZ(1:NORDZ),WV(ii,IDZ(1:NORDZ)))
    enddo
    
    DY_Y(2)=-wky*WV(2,ID)
    
    do jj=1,3
    DZ_Y(jj)=dot_product(DEZ(ID,1:NORDZ),WV(jj,IDZ(ID,1:NORDZ)))
    enddo
    
    DX_Y(4)=dot_product(DEX(ID,1:NORDX),P(IDX(ID,1:NORDX)))
    !DX_Y(4)=DX_Y(4)+DFS(2,IZ,3)*dot_product(DEZ(1:NORDZ),P(IDZ(1:NORDZ)))
    
    DY_Y(4)=wky*P(ID)
    
    DZ_Y(4)=dot_product(DEZ(ID,1:NORDZ),P(IDZ(ID,1:NORDZ)))
    endif
    
    !-------------
    !---Only Solid
    if(MNO(ID)==1.or.MNO(ID)==2)then
    !velocity
    do ii=5,7
        DX_Y(ii)=dot_product(DEX(ID,1:NORDX),V(ii-4,IDX(ID,1:NORDX)))
    enddo
    
    DY_Y(5:7)=wky*V(1:3,ID)
    DY_Y(6)=-DY_Y(6)
    
    do jj=5,7
        DZ_Y(jj)=dot_product(DEZ(ID,1:NORDZ),V(jj-4,IDZ(ID,1:NORDZ)))
    enddo 
    
    !Stress 1
    do ii=8,10
    DX_Y(ii)=dot_product(DEX(ID,1:NORDX),ST1(ii-7,IDX(ID,1:NORDX)))
    enddo
    
    DY_Y(8:10)=wky*ST1(1:3,ID)
    DY_Y(9)=-DY_Y(9)
    
    do jj=8,10
    DZ_Y(jj)=dot_product(DEZ(ID,1:NORDZ),ST1(jj-7,IDZ(ID,1:NORDZ)))
    enddo
    
    !Stress 2
    do ii=11,13
        DX_Y(ii)=dot_product(DEX(ID,1:NORDX),ST2(ii-10,IDX(ID,1:NORDX)))
    enddo
    
    DY_Y(11:13)=wky*ST2(1:3,ID)
    DY_Y(12)=-DY_Y(12)
    
    do jj=11,13
    DZ_Y(jj)=dot_product(DEZ(ID,1:NORDZ),ST2(jj-10,IDZ(ID,1:NORDZ)))
    enddo
    
    endif
    
    DX_Y(1:13)=DFS(ID,1,IX,2)*DX_Y(1:13)+DFS(ID,1,IX,3)*DZ_Y(1:13)
    DZ_Y(1:13)=DFS(ID,2,IZ,4)*DZ_Y(1:13)
    
return
end subroutine

!-----------------------------------------------
subroutine invert3(a,LDA,na,det_l)
 implicit none
 complex*16, intent(inout) :: a (LDA,na)
 integer, intent(in)             :: LDA
 integer, intent(in)             :: na
 complex*16, intent(inout) :: det_l
 complex*16 :: b(4,3)
 integer :: i
 double precision :: f
 det_l = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
        -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
        +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
 do i=1,4
  b(i,1) = a(i,1)
  b(i,2) = a(i,2)
  b(i,3) = a(i,3)
 enddo
 a(1,1) =  b(2,2)*b(3,3) - b(2,3)*b(3,2)
 a(2,1) =  b(2,3)*b(3,1) - b(2,1)*b(3,3)
 a(3,1) =  b(2,1)*b(3,2) - b(2,2)*b(3,1)

 a(1,2) =  b(1,3)*b(3,2) - b(1,2)*b(3,3)
 a(2,2) =  b(1,1)*b(3,3) - b(1,3)*b(3,1)
 a(3,2) =  b(1,2)*b(3,1) - b(1,1)*b(3,2)

 a(1,3) =  b(1,2)*b(2,3) - b(1,3)*b(2,2)
 a(2,3) =  b(1,3)*b(2,1) - b(1,1)*b(2,3)
 a(3,3) =  b(1,1)*b(2,2) - b(1,2)*b(2,1)

end

subroutine invert3_Real(a,LDA,na,det_l)
 implicit none
 doubleprecision, intent(inout) :: a (LDA,na)
 integer, intent(in)             :: LDA
 integer, intent(in)             :: na
 doubleprecision, intent(inout) :: det_l
 doubleprecision :: b(4,3)
 integer :: i
 double precision :: f
 det_l = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
        -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
        +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
 do i=1,4
  b(i,1) = a(i,1)
  b(i,2) = a(i,2)
  b(i,3) = a(i,3)
 enddo
 a(1,1) =  b(2,2)*b(3,3) - b(2,3)*b(3,2)
 a(2,1) =  b(2,3)*b(3,1) - b(2,1)*b(3,3)
 a(3,1) =  b(2,1)*b(3,2) - b(2,2)*b(3,1)

 a(1,2) =  b(1,3)*b(3,2) - b(1,2)*b(3,3)
 a(2,2) =  b(1,1)*b(3,3) - b(1,3)*b(3,1)
 a(3,2) =  b(1,2)*b(3,1) - b(1,1)*b(3,2)

 a(1,3) =  b(1,2)*b(2,3) - b(1,3)*b(2,2)
 a(2,3) =  b(1,3)*b(2,1) - b(1,1)*b(2,3)
 a(3,3) =  b(1,1)*b(2,2) - b(1,2)*b(2,1)

end
    