!************************************************************************
!
! User material subroutine (UMAT) for large-deformation non-linear response of polymers 
! that incorporates relaxation. This UMAT is not for use in plane stress or in any other
!  situation in which there are more strain terms than stress terms. This UMAT does 
!  not support energy outputs.
!
! Pinkesh Malhotra & Alexander Landauer, April 2016
!
!************************************************************************
! Usage:
!************************************************************************
!
! Throughout, we employ a notation in which a quantity with subscript 
!  tau denotes that quantity at the end of the increment, while a 
!  subscript t denotes the quantity at the beginning of the increment.
!
!     Material Properties Vector (*user material, constants = 7)
!     --------------------------------------------------------------
!     Geq       = props(1)  ! Shear modulus for equilibrium response (Neo-Hookean Model used here)
!     Keq       = props(2)  ! Bulk modulus for equilibrium response (Neo-Hookean Model used here)
!     Gneq1     = props(3)  ! Shear modulus of first branch for non-equilibrium response 
!     Gneq2     = props(4)  ! Shear modulus of second branch for non-equilibrium response
!     Gneq3     = props(5)  ! Shear modulus of third branch for non-equilibrium response
!     t1        = props(6)  ! Relaxation time for first branch of non-equilibrium response
!     t2        = props(7)  ! Relaxation time for second branch of non-equilibrium response
!     t3        = props(8)  ! Relaxation time for third branch of non-equilibrium response
!    
!      *NOTE: Only 3 non-equilibrium mechanisms/components have been assumed.
!
!
!                         Rheological Analogue Schematic
!
!                         
!              |      Neo-Hookean equilibrium branch        |
!              |             /\    /\    /\                 |
!              |------------/  \  /  \  /  \  /-------------|
!              |                \/    \/    \/              |
!              |                                            |
!  		       |                                            |
!			   | Linder et al.-based Dissipative Branch (1) | 
!              |        /\    /\               |-----       |
!              |-------/  \  /  \  /-----------|  |---------|
!              |           \/    \/            |-----       |
!      <------ |                                            | ------>
!              | Linder et al.-based Dissipative Branch (2) |
!              |        /\    /\               |-----       |
!              |-------/  \  /  \  /-----------|  |---------|
!              |           \/    \/            |-----       |
!              |                                            |
!              | Linder et al.-based Dissipative Branch (3) |
!              |        /\    /\               |-----       |
!              |-------/  \  /  \  /-----------|  |---------|
!              |           \/    \/            |-----       |
!              |                                            |
!			   
!     State Variables (*depvar 18) :
!     --------------------------------------------------------------
!     statev(1) = A1(1,1) ---- Internal Variable(1,1)
!     statev(2) = A1(2,2) ---- Internal Variable(2,2)
!     statev(3) = A1(3,3) ---- Internal Variable(3,3)
!     statev(4) = A1(2,3) ---- Internal Variable(2,3)
!     statev(5) = A1(1,3) ---- Internal Variable(1,3)
!     statev(6) = A1(1,2) ---- Internal Variable(1,2)
!     statev(7) = A2(1,1) ---- Internal Variable(1,1)
!     statev(8) = A2(2,2) ---- Internal Variable(2,2)
!     statev(9) = A2(3,3) ---- Internal Variable(3,3)
!     statev(10) = A2(2,3) ---- Internal Variable(2,3)
!     statev(11) = A2(1,3) ---- Internal Variable(1,3)
!     statev(12) = A2(1,2) ---- Internal Variable(1,2)
!     statev(13) = A3(1,1) ---- Internal Variable(1,1)
!     statev(14) = A3(2,2) ---- Internal Variable(2,2)
!     statev(15) = A3(3,3) ---- Internal Variable(3,3)
!     statev(16) = A3(2,3) ---- Internal Variable(2,3)
!     statev(17) = A3(1,3) ---- Internal Variable(1,3)
!     statev(18) = A3(1,2) ---- Internal Variable(1,2)
!
!     *NOTE: the matrix A is always symmetric and 6 values need be stored
!
!************************************************************************
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     +     rpl,ddsddt,drplde,drpldt,
     +     stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     +     ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     +     celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
      !
      include 'aba_param.inc'
      !
      dimension stress(ntens),statev(nstatv),
     +     ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     +     stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     +     props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
      !
      character*8 cmname
      !
      ! Variables defined and used in the UMAT
      !
      integer i,j
      !
      real*8 Geq,Keq,Gneq1,Gneq2,Gneq3,t1,t2,t3,stat,
     +  Iden(3,3), A1(3,3),A2(3,3),A3(3,3),C_tau(3,3),B_tau(3,3),
     +  J_tau,F_tau(3,3),Finv_tau(3,3),Cinv_tau(3,3),Jsq_tau,tr_B_tau,
     +  Teq(3,3),TneqT(3,3),Tneq1(3,3),Tneq2(3,3),Tneq3(3,3),Q1(3,3),
     +  Q2(3,3),Q3(3,3),R1,R2,R3,
     +  T(3,3),Ctan_eq(3,3,3,3),Ctan_neq1(3,3,3,3),Ctan_neq2(3,3,3,3),
     +  Ctan_neq3(3,3,3,3),Ctan(3,3,3,3) 

      real*8 zero,one,two,three,fourth,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,fourth=1.d0/4.d0,
     +    third=1.d0/3.d0,half=1.d0/2.d0,twothird=2.d0/3.d0,
     +    twonine=2.d0/9.d0,small=0.0000001)

      ! Do nothing if a "dummy" step
      if(dtime.eq.zero) return
 

      ! Identity matrix
      call onem(Iden)
	 
      ! Material Properties
      Geq       = props(1)  
      Keq       = props(2)  
      Gneq1     = props(3)  
      Gneq2     = props(4)  
      Gneq3     = props(5)  
      t1        = props(6)  
      t2        = props(7)  
      t3        = props(8) 
      !write(*,*) Geq
      !write(*,*) Keq
      !write(*,*) Gneq1
      !write(*,*) Gneq2
      !write(*,*) Gneq3
      !write(*,*) t1
      !write(*,*) t2
      !write(*,*) t3
	  
	! Get the deformation gradient at the end of increment.
      F_tau=0.d0
      F_tau(1,1)=dfgrd1(1,1)
      F_tau(2,2)=dfgrd1(2,2)
      F_tau(3,3)=dfgrd1(3,3)
      if (nshr.ne.0) then
          F_tau(1,2)=dfgrd1(1,2)
          F_tau(2,1)=dfgrd1(2,1)
        if (nshr.ne.1) then
            F_tau(1,3)=dfgrd1(1,3)
            F_tau(3,1)=dfgrd1(3,1)
          if (nshr.ne.2) then
              F_tau(2,3)=dfgrd1(2,3)
              F_tau(3,2)=dfgrd1(3,2)
          endif
        endif
      endif
      !write(*,*) F_tau
    ! Find the inverse of F and Jacobian at the end of the increment.
      Finv_tau=0.d0
      J_tau=0.d0
      call matInv3D(F_tau,Finv_tau,J_tau,stat)

    ! Find the Cauchy-Green Tensors.
      C_tau=0.d0
      B_tau=0.d0
      Cinv_tau=0.d0
      C_tau = matmul(transpose(F_tau),F_tau)
      B_tau = matmul(F_tau,transpose(F_tau))
      call matInv3D(C_tau,Cinv_tau,Jsq_tau,stat)
      !write(*,*) C_tau
      !write(*,*) B_tau
      ! At the start of an Abaqus calculation, the state
      !  variables are passed into UMAT with zero values.
      !  Initialize the state variables. At this point,
      !  the time total_time and step_time both have a value
		!  equal to zero and the step counter, kstep, is 
      !  equal to 1.
      !
      if ((time(1).eq.0.d0).and.(kstep.eq.1)) then
        statev(1)  = 1.d0
        statev(2)  = 1.d0
        statev(3)  = 1.d0 
        statev(4)  = 0.d0
        statev(5)  = 0.d0
        statev(6)  = 0.d0
        statev(7)  = 1.d0
        statev(8)  = 1.d0
        statev(9)  = 1.d0
        statev(10)  = 0.d0
        statev(11)  = 0.d0
        statev(12)  = 0.d0
        statev(13)  = 1.d0
        statev(14)  = 1.d0
        statev(15)  = 1.d0
        statev(16)  = 0.d0
        statev(17)  = 0.d0
        statev(18)  = 0.d0
      end if
      
      ! Store the state variables in the matrices, A1, A2 and A3.
      A1(1,1) = statev(1)
      A1(2,2) = statev(2)
      A1(3,3) = statev(3)
      A1(1,2) = statev(6)
      A1(2,1) = statev(6)
      A1(1,3) = statev(5)
      A1(3,1) = statev(5)
      A1(2,3) = statev(4)
      A1(3,2) = statev(4)
      !
      A2(1,1) = statev(7)
      A2(2,2) = statev(8)
      A2(3,3) = statev(9)
      A2(1,2) = statev(12)
      A2(2,1) = statev(12)
      A2(1,3) = statev(11)
      A2(3,1) = statev(11)
      A2(2,3) = statev(10)
      A2(3,2) = statev(10)
      !
      A3(1,1) = statev(13)
      A3(2,2) = statev(14)
      A3(3,3) = statev(15)
      A3(1,2) = statev(18)
      A3(2,1) = statev(18)
      A3(1,3) = statev(17)
      A3(3,1) = statev(17)
      A3(2,3) = statev(16)
      A3(3,2) = statev(16)
      ! 
      
      ! Update the internal variable Matrix, A, the state variables. 
	  !    NOTE: No suffix 't' or 'tau' is used.
      A1(1,1)=(1.d0/(1.d0+(dtime/t1)))*(A1(1,1)+(dtime/t1)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(1,1))
      A1(2,2)=(1.d0/(1.d0+(dtime/t1)))*(A1(2,2)+(dtime/t1)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(2,2))
      A1(3,3)=(1.d0/(1.d0+(dtime/t1)))*(A1(3,3)+(dtime/t1)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(3,3))
      A1(1,2)=(1.d0/(1.d0+(dtime/t1)))*(A1(1,2)+(dtime/t1)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(1,2))
      A1(2,1)=A1(1,2)
      A1(1,3)=(1.d0/(1.d0+(dtime/t1)))*(A1(1,3)+(dtime/t1)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(1,3))
      A1(3,1)=A1(1,3)
      A1(2,3)=(1.d0/(1.d0+(dtime/t1)))*(A1(2,3)+(dtime/t1)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(2,3))
      A1(3,2)=A1(2,3)
      ! Repeat for the second branch
      A2(1,1)=(1.d0/(1.d0+(dtime/t2)))*(A2(1,1)+(dtime/t2)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(1,1))
      A2(2,2)=(1.d0/(1.d0+(dtime/t2)))*(A2(2,2)+(dtime/t2)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(2,2))
      A2(3,3)=(1.d0/(1.d0+(dtime/t2)))*(A2(3,3)+(dtime/t2)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(3,3))
      A2(1,2)=(1.d0/(1.d0+(dtime/t2)))*(A2(1,2)+(dtime/t2)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(1,2))
      A2(2,1)=A2(1,2)
      A2(1,3)=(1.d0/(1.d0+(dtime/t2)))*(A2(1,3)+(dtime/t2)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(1,3))
      A2(3,1)=A2(1,3)
      A2(2,3)=(1.d0/(1.d0+(dtime/t2)))*(A2(2,3)+(dtime/t2)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(2,3))
      A2(3,2)=A2(2,3)
      !Again for the third branch
      A3(1,1)=(1.d0/(1.d0+(dtime/t3)))*(A3(1,1)+(dtime/t3)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(1,1))
      A3(2,2)=(1.d0/(1.d0+(dtime/t3)))*(A3(2,2)+(dtime/t3)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(2,2))
      A3(3,3)=(1.d0/(1.d0+(dtime/t3)))*(A3(3,3)+(dtime/t3)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(3,3))
      A3(1,2)=(1.d0/(1.d0+(dtime/t3)))*(A3(1,2)+(dtime/t3)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(1,2))
      A3(2,1)=A3(1,2)
      A3(1,3)=(1.d0/(1.d0+(dtime/t3)))*(A3(1,3)+(dtime/t3)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(1,3))
      A3(3,1)=A3(1,3)
      A3(2,3)=(1.d0/(1.d0+(dtime/t3)))*(A3(2,3)+(dtime/t3)*
     + (J_tau**(2.d0/3.d0))*Cinv_tau(2,3))
      A3(3,2)=A3(2,3)

      ! Find the equilibrium Cauchy Stress for Neo-Hookean Model.
      Teq=0.d0
      if(J_tau.lt.small) then
          write(*,*) 'J_tau numerically unstable' 
      end if
      tr_B_tau=B_tau(1,1)+B_tau(2,2)+B_tau(3,3)
      Teq(1,1)=(Geq/(J_tau**(5.d0/3.d0)))*(B_tau(1,1)-
     + (tr_B_tau/3.d0)*Iden(1,1))+Keq*(J_tau-1)*Iden(1,1)
      Teq(2,2)=(Geq/(J_tau**(5.d0/3.d0)))*(B_tau(2,2)-
     + (tr_B_tau/3.d0)*Iden(2,2))+Keq*(J_tau-1)*Iden(2,2)
      Teq(3,3)=(Geq/(J_tau**(5.d0/3.d0)))*(B_tau(3,3)-
     + (tr_B_tau/3.d0)*Iden(3,3))+Keq*(J_tau-1)*Iden(3,3)
      Teq(1,2)=(Geq/(J_tau**(5.d0/3.d0)))*(B_tau(1,2)-
     + (tr_B_tau/3.d0)*Iden(1,2))+Keq*(J_tau-1)*Iden(1,2)
      Teq(1,3)=(Geq/(J_tau**(5.d0/3.d0)))*(B_tau(1,3)-
     + (tr_B_tau/3.d0)*Iden(1,3))+Keq*(J_tau-1)*Iden(1,3)
      Teq(2,3)=(Geq/(J_tau**(5.d0/3.d0)))*(B_tau(2,3)-
     + (tr_B_tau/3.d0)*Iden(2,3))+Keq*(J_tau-1)*Iden(2,3)
	 
      ! Find the non-equilibrium Cauchy Stress
      TneqT=0.d0
      Tneq1=0.d0
      Tneq2=0.d0
      Tneq3=0.d0
      Q1=matmul(F_tau,matmul(A1,transpose(F_tau)))
      Q2=matmul(F_tau,matmul(A2,transpose(F_tau)))
      Q3=matmul(F_tau,matmul(A3,transpose(F_tau)))
      call contraction3D(A1,C_tau,R1)
      call contraction3D(A2,C_tau,R2)
      call contraction3D(A3,C_tau,R3)
      ! First branch 
      Tneq1(1,1)=(Gneq1/(J_tau**(5/3)))*(Q1(1,1)-
     + (R1/3.d0)*Iden(1,1))
      Tneq1(2,2)=(Gneq1/(J_tau**(5/3)))*(Q1(2,2)-
     + (R1/3.d0)*Iden(2,2))
      Tneq1(3,3)=(Gneq1/(J_tau**(5/3)))*(Q1(3,3)-
     + (R1/3.d0)*Iden(3,3))
      Tneq1(1,2)=(Gneq1/(J_tau**(5/3)))*(Q1(1,2)-
     + (R1/3.d0)*Iden(1,2))
      Tneq1(1,3)=(Gneq1/(J_tau**(5/3)))*(Q1(1,3)-
     + (R1/3.d0)*Iden(1,3))
      Tneq1(2,3)=(Gneq1/(J_tau**(5/3)))*(Q1(2,3)-
     + (R1/3.d0)*Iden(2,3))
      ! second branch
      Tneq2(1,1)=(Gneq2/(J_tau**(5/3)))*(Q2(1,1)-
     + ((R2)/3.d0)*Iden(1,1))
      Tneq2(2,2)=(Gneq2/(J_tau**(5/3)))*(Q2(2,2)-
     + ((R2)/3.d0)*Iden(2,2))
      Tneq2(3,3)=(Gneq2/(J_tau**(5/3)))*(Q2(3,3)-
     + ((R2)/3.d0)*Iden(3,3))
      Tneq2(1,2)=(Gneq2/(J_tau**(5/3)))*(Q2(1,2)-
     + ((R2)/3.d0)*Iden(1,2))
      Tneq2(1,3)=(Gneq2/(J_tau**(5/3)))*(Q2(1,3)-
     + ((R2)/3.d0)*Iden(1,3))
      Tneq2(2,3)=(Gneq2/(J_tau**(5/3)))*(Q2(2,3)-
     + ((R2)/3.d0)*Iden(2,3))
      ! third branch
      Tneq3(1,1)=(Gneq3/(J_tau**(5/3)))*(Q3(1,1)-
     + ((R3)/3.d0)*Iden(1,1))
      Tneq3(2,2)=(Gneq3/(J_tau**(5/3)))*(Q3(2,2)-
     + ((R3)/3.d0)*Iden(2,2))
      Tneq3(3,3)=(Gneq3/(J_tau**(5/3)))*(Q3(3,3)-
     + ((R3)/3.d0)*Iden(3,3))
      Tneq3(1,2)=(Gneq3/(J_tau**(5/3)))*(Q3(1,2)-
     + ((R3)/3.d0)*Iden(1,2))
      Tneq3(1,3)=(Gneq3/(J_tau**(5/3)))*(Q3(1,3)-
     + ((R3)/3.d0)*Iden(1,3))
      Tneq3(2,3)=(Gneq3/(J_tau**(5/3)))*(Q3(2,3)-
     + ((R3)/3.d0)*Iden(2,3))     
      
      ! Find total Cauchy Stress
      T=0.d0
      T(1,1)=Teq(1,1)+Tneq1(1,1)+Tneq2(1,1)+Tneq3(1,1)
      T(2,2)=Teq(2,2)+Tneq1(2,2)+Tneq2(2,2)+Tneq3(2,2)
      T(3,3)=Teq(3,3)+Tneq1(3,3)+Tneq2(3,3)+Tneq3(3,3)
      T(1,2)=Teq(1,2)+Tneq1(1,2)+Tneq2(1,2)+Tneq3(1,2)
      T(1,3)=Teq(1,3)+Tneq1(1,3)+Tneq2(1,3)+Tneq3(1,3)
      T(2,3)=Teq(2,3)+Tneq1(2,3)+Tneq2(2,3)+Tneq3(2,3)
      T(2,1)=T(1,2)
      T(3,1)=T(1,3)
      T(3,2)=T(2,3)
      
      ! Convert to 1st P-K stress. NOTE: The same symbol is used.
      T=matmul(T,transpose(Finv_tau))
      T(1,1)=J_tau*T(1,1)
      T(2,2)=J_tau*T(2,2)
      T(3,3)=J_tau*T(3,3)
      T(1,2)=J_tau*T(1,2)
      T(2,1)=J_tau*T(2,1)
      T(1,3)=J_tau*T(1,3)
      T(3,1)=J_tau*T(3,1)
      T(2,3)=J_tau*T(2,3)
      T(3,2)=J_tau*T(3,2)
      write(*,*)T(2,2)
      ! Update the stress	
      stress = 0.d0
      do i=1,ndi
        stress(i) = T(i,i)
      end do
      if (nshr.ne.0) then
        stress(ndi+1) = T(1,2)
        if (nshr.ne.1) then
          stress(ndi+2) = T(1,3)
          if (nshr.ne.2) then
            stress(ndi+3) = T(2,3)
          endif
        endif
      endif
      
      ! Find the tangents
      Ctan_eq=0.d0
      Ctan_neq1=0.d0
      Ctan_neq2=0.d0
      Ctan_neq3=0.d0
      Ctan=0.d0
      do i=1,3
      	do j=1,3
      		do k=1,3
      			do l=1,3
      				Ctan_eq(i,j,k,l)=(Geq*(J_tau**(-5.d0/3.d0)))*
     +                 ((1.d0/2.d0)*(Iden(i,k)*B_tau(j,l)+Iden(j,l)*
     +                 B_tau(i,k)+Iden(i,l)*B_tau(j,k)+
     +                 Iden(j,k)*B_tau(i,l))-
     +                 (2.d0/3.d0)*(Iden(i,j)*B_tau(k,l)+
     +                 Iden(k,l)*B_tau(i,j))+
     +                 (2.d0/9.d0)*(tr_B_tau*Iden(i,j)*Iden(k,l)))+
     +                 Keq*(2.d0*J_tau-1.d0)*Iden(i,j)*Iden(k,l)
                      !
      				Ctan_neq1(i,j,k,l)=(Gneq1*(J_tau**(-5.d0/3.d0)))*
     +                 (-(2.d0/3.d0)*(Q1(i,j)*
     +                 Iden(k,l)-(R1/3.d0)*Iden(i,j)*Iden(k,l))+
     +                 (1.d0/2.d0)*(Q1(l,j)*Iden(i,k)+Q1(i,l)*Iden(j,k)+
     +                 Q1(k,j)*Iden(i,l)+Q1(i,k)*Iden(j,l))-
     +                 (1.d0/3.d0)*(Q1(l,k)+Q1(k,l))*Iden(i,j))
                      !
      				Ctan_neq2(i,j,k,l)=(Gneq2*(J_tau**(-5.d0/3.d0)))*
     +                 (-(2.d0/3.d0)*(Q2(i,j)*
     +                 Iden(k,l)-(R2/3.d0)*Iden(i,j)*Iden(k,l))+
     +                 (1.d0/2.d0)*(Q2(l,j)*Iden(i,k)+Q2(i,l)*Iden(j,k)+
     +                 Q2(k,j)*Iden(i,l)+Q2(i,k)*Iden(j,l))-
     +                 (1.d0/3.d0)*(Q2(l,k)+Q2(k,l))*Iden(i,j))
                      !
      				Ctan_neq3(i,j,k,l)=(Gneq3*(J_tau**(-5.d0/3.d0)))*
     +                 (-(2.d0/3.d0)*(Q3(i,j)*
     +                 Iden(k,l)-(R3/3.d0)*Iden(i,j)*Iden(k,l))+
     +                 (1.d0/2.d0)*(Q3(l,j)*Iden(i,k)+Q3(i,l)*Iden(j,k)+
     +                 Q3(k,j)*Iden(i,l)+Q3(i,k)*Iden(j,l))-
     +                 (1.d0/3.d0)*(Q3(l,k)+Q3(k,l))*Iden(i,j))
                      !
      				Ctan(i,j,k,l)=Ctan(i,j,k,l)+Ctan_eq(i,j,k,l)+
     +                  Ctan_neq1(i,j,k,l)+Ctan_neq2(i,j,k,l)+
     +                 Ctan_neq3(i,j,k,l)
      			end do
      		end do
      	end do
      end do

      ! Update the tangents
      ddsdde = 0.d0
      !
      do i=1,ndi
        do j=1,ndi
          ddsdde(i,j) = Ctan(i,i,j,j)
        end do
      end do
      !
      if (nshr.ne.0) then
        do i=1,ndi
          ddsdde(i,ndi+1) = Ctan(i,i,1,2)
          ddsdde(ndi+1,i) = Ctan(1,2,i,i)
        end do
        ddsdde(ndi+1,ndi+1) = Ctan(1,2,1,2)
        if (nshr.ne.1) then
          do i=1,ndi
            ddsdde(i,ndi+2) = Ctan(i,i,1,3)
            ddsdde(ndi+2,i) = Ctan(1,3,i,i)
          end do
          ddsdde(ndi+2,ndi+2) = Ctan(1,3,1,3)
          ddsdde(ndi+1,ndi+2) = Ctan(1,2,1,3)
          ddsdde(ndi+2,ndi+1) = Ctan(1,3,1,2)
          if (nshr.ne.2) then
            do i=1,ndi
              ddsdde(i,ndi+3) = Ctan(i,i,2,3)
              ddsdde(ndi+3,i) = Ctan(2,3,i,i)
            end do
            ddsdde(ndi+3,ndi+3) = Ctan(2,3,2,3)
            ddsdde(ndi+1,ndi+3) = Ctan(1,2,2,3)
            ddsdde(ndi+3,ndi+1) = Ctan(2,3,1,2)
            ddsdde(ndi+2,ndi+3) = Ctan(1,3,2,3)
            ddsdde(ndi+3,ndi+2) = Ctan(2,3,1,3)
          endif
        endif
      endif

      !Update the state variables
      statev(1) = A1(1,1) 
      statev(2) = A1(2,2) 
      statev(3) = A1(3,3) 
      statev(4) = A1(2,3) 
      statev(5) = A1(1,3) 
      statev(6) = A1(1,2)
      statev(7) = A2(1,1) 
      statev(8) = A2(2,2) 
      statev(9) = A2(3,3) 
      statev(10) = A2(2,3) 
      statev(11) = A2(1,3) 
      statev(12) = A2(1,2)
      statev(13) = A3(1,1) 
      statev(14) = A3(2,2) 
      statev(15) = A3(3,3) 
      statev(16) = A3(2,3) 
      statev(17) = A3(1,3) 
      statev(18) = A3(1,2)

      return
      end subroutine umat

!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv3D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D

!****************************************************************************

      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv

      
      istat = 1
      
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = 1.d0/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
      end subroutine matInv2D

!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet

!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
         do J=1,3
            if (i .eq. j) then
              A(i,j) = 1.d0
            else
              A(i,j) = 0.d0
            end if
         end do
      end do


      return
      end subroutine onem

!****************************************************************************

        subroutine multiply3D(A,B,C)
		!
		! Multiplies 3x3 matrices A and B to give B
		!
        implicit none
		!
        real*8 A(3,3), B(3,3), C(3,3)
		!
        C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
        C(1,2)=A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
        C(1,3)=A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)
        C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
        C(2,2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
        C(2,3)=A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
        C(3,1)=A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
        C(3,2)=A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
        C(3,3)=A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)
		
      return
      end subroutine multiply3D
		
!*****************************************************************************

      subroutine contraction3D(A,B,C)
		!
		! C=A:B
		! 
        implicit none
		!
        real*8 A(3,3),B(3,3),C
		!
        C=A(1,1)*B(1,1)+A(1,2)*B(1,2)+A(1,3)*B(1,3)+
     +	   A(2,1)*B(2,1)+A(2,2)*B(2,2)+A(2,3)*B(2,3)+
     +     A(3,1)*B(3,1)+A(3,2)*B(3,2)+A(3,3)*B(3,3)
	 
      return
      end subroutine contraction3D