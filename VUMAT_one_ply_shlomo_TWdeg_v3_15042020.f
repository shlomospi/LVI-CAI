		  subroutine vumat(
c Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
c Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
c
c*************************************************************************
c*    Tsai-Wu progressive degradation Model for a Composite layer        *                                                                        
c*                                                                       *
c*                         By Shlomo Spitzer                             *
c*                      Last Update    18/11/2019                        *
c*   Abaqus Notation: for stress/strain (11,22,33,12,23,31)              *                                                  
c*************************************************************************
      include 'vaba_param.inc'
c************************************ 
c*   Global Arrays/Variables
c*   All arrays dimensioned by (*) are not used in this algorithm
c************************************
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(*), strainInc(nblock,ndir+nshr),
     2  relSpinInc(*), tempOld(*),
     3  stretchOld(*),
     4  defgradOld(*),
     5  fieldOld(*), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(*),
     6  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     7  enerInelasOld(*), tempNew(*),
     8  stretchNew(*),
     8  defgradNew(*),
     9  fieldNew(*),
     2  enerInternNew(*), enerInelasNew(*)
c
      character*80 cmname

c***************************************************************************** 
c   Local Arrays/scalars Variables    
c***************************************************************************** 
       integer i1,i,j
      real*8 E11,E22,E33,G23,G13,G12,nu23,nu13,nu12               ! Effective of one layer, 1 is the fiber direction. ti is interface layer's thickness & t is total thickness
	real*8 X_t,X_c,Y_t,Y_c,S_i,S_o                                  ! Maximum stresses for a layer, by direction ( for the Tsai Wu criterion):
                                                                      ! Fiber tension, Fiber compression, Matrix tension, Matrix compression, In plain shear, Out of Plain shear
	real*8 f1,f2,f11,f22,f66,f12                                       ! Tsai Wu parameters, calculated from: X_t,X_c,Y_t,Y_c,S_i,S_o
	real*8 CL_A_Layer(6,6)
	real*8 EG(6),SG(6)                                              ! golbal strain & stress vectors. abaqus notation                                       
	real*8 EG_layer(6),SG_layer(6)                                    ! holds the effective OOP stress in the interface layer
	real*8 TW(6),Max_TW,Phi,i_Max_TW ! TW parameters
      real*8 D11,D22,D33
      real*8 beta1t,beta1c,beta2t,beta2c,beta12
      real*8 gamma1t,gamma1c,gamma2t,gamma2c,gamma12
      real*8 E11t,E11c,E22t,E22c,E33t,E33c
c***************************************************************************** 
c     The state variables 
c***************************************************************************** 
c SDV(1)  == fiber damage counter
C StateNew(i,1) == E11
C StateNew(i,2) == E22
C StateNew(i,3) == E33
C StateNew(i,4) == E12
C StateNew(i,5) == E23
C StateNew(i,6) == E13
C StateNew(i,7) == Flag for damage initiation and reduction to zero of non diagonal components
C StateNew(i,8) == TW value
C StateNew(i,9) == TW dominant value index (1-fiber 2-matrix 3-shear)
C StateNew(i,10) == Flag for deletion (1=element exists; 0=element failed)
C StateNew(i,11) == D1t
C StateNew(i,12) == D1c
C StateNew(i,13) == D2t
C StateNew(i,14) == D2c
C StateNew(i,15) == D12
C StateNew(i,16) == E11_I_1t
C StateNew(i,17) == E11_I_1c
C StateNew(i,18) == E22_I_2t
C StateNew(i,19) == E22_I_2c
C StateNew(i,20) == E12_I_12
C StateNew(i,21) == TW fiber value
C StateNew(i,22) == TW matrix value
C StateNew(i,23) == TW shear value    
C******************************************************************************
C     Define Material
C******************************************************************************
C     Make sure that stresses are represented in units of [Pa]. This is the case
C     when the lengths are in meters in the analysis.
C     
C
      !Density
C       Density should be defined in the material outside the VUMAT, for this material it is:
C
C       density = 1584.9 [kg/m^3]
C
      ! Orthotropic layer's properties: IM7/977-3 !MI7 !(T300/NY9200Z) 
                       !  E11=163810d6 !137e9 ! [Pa]
                       !  E22=8790d6   !9.3e9 ! [Pa]
                       !  E33=8820d6   !9.3e9 ! [Pa]
        G12 = 5677.4d6 ! 4160d6   !5.3e9 ! [Pa]
        G23 = 3141.8d6 ! 3110d6   !3.8e9 ! [Pa]
        G13 = 5215.7d6 ! 4160d6   !5.3e9 ! [Pa]
        nu12=0.32      !.31
        nu23=0.461     !.3
        nu13=0.329     !.31 
C
        E11t = 158d9
        E22t = 8.97d9
        E33t = 8.89d9
        E11c = 130d9
        E22c = 8.69d9
        E33c = 8.63d9
c
      ! critical stress value: IM7/977-3 ! IM7 !(T300/NY9200Z)
        X_t = 2704d6   !2764.3d6  !1747e6    ! [Pa]
        X_c = 1764.3d6 !2150.01d6 !1357e6    ! [Pa]
        Y_t = 95.7d6   !79.75d6   !67e6 !78  ! [Pa]
        Y_c = 236.5d6  !249.08d6  !170e6     ! [Pa]
        S_i = 117.9d6  ! 260d6               ! [Pa]
        S_o = 96.5d6   !98.19d6   !124e6     ! [Pa]
c		
      ! the material parameters of modified "Tsai-Wu" failure criteria, 1 is the fiber direction:
	  f1=1/X_t-1/X_c
	  f2=1/Y_t-1/Y_c
        f11=1/(X_t*X_c)
	  f22=1/(Y_t*Y_c)
	  f66=1/(S_i)**2
	  f12=-0.5*dsqrt(DABS(f11*f22))
        !write(*,*) "f1,f2,f12", f1,f2,f12
C
        ! progressive damage parameters. # corresponds to damage mode#. t for tension and c for compression. 
        ! the value represents how many times the maximal strain is larger than the strain at damage initiation, assuming the element is not deleted on the way.  
        beta1t = 2.5
        beta1c = 2.5
        beta2t = 3
        beta2c = 3
        beta12 = 3
        ! Gamma parameters represent the maximal damage progression allowed in each mode
        gamma1t = 0.8
        gamma1c = 0.8
        gamma2t = 0.8
        gamma2c = 0.8
        gamma12 = 0.8
C
c****************************************************************************
C
      do 100 i= 1,nblock 
c  
c******************************************************************************

c calculating the strain vector EG(6) at the material point i
c according to Abaqus VUMAT notation 
c
      EG(1) = stateold(i,1) + strainInc(i,1)
      EG(2) = stateold(i,2) + strainInc(i,2)
	EG(3) = stateold(i,3) + strainInc(i,3)
	EG(4) = stateold(i,4) + strainInc(i,4)
	EG(5) = stateold(i,5) + strainInc(i,5)
	EG(6) = stateold(i,6) + strainInc(i,6)
C
      stateNew(i,1)=EG(1)
      stateNew(i,2)=EG(2)
	stateNew(i,3)=EG(3)
	stateNew(i,4)=EG(4)
	stateNew(i,5)=EG(5)
	stateNew(i,6)=EG(6)
C
C*****************************************************************************
c calculating the Local stiffness matrix [CL_A] for the N+1 orthotropic layers
c of the sublaminate model according to Abaqus VUMAT notation [CL_A]	
c
      do i1 = 1,6
          do j = 1,6
              CL_A_layer(i1,j) = 0
          end do
      end do
C
      if (EG(1).lt.0) then
          E11 = E11c
      else
          E11 = E11t
      end if
C
      if (EG(2).lt.0) then
          E22 = E22c
      else
          E22 = E22t
      end if
C
      if (EG(3).lt.0) then
          E33 = E33c
      else
          E33 = E33t
      end if
C
      if (stateOld(i,7).gt.0.9) then ! damage initiated
C
          
          if (EG(1).gt.0) then ! fiber tension
              D11 = StateOld(i,11)
          else ! fiber compression
              D11 = StateOld(i,12)
          end if
C
          if (EG(2).gt.0) then ! Matrix tension
              D22 = StateOld(i,13)
          else ! Matrix compression
              D22 = StateOld(i,14)
          end if 
C
          D33 = StateOld(i,15) ! one mode for shear
C
 
          call get_C_mat_OfLayer_local_Aba_damaged(E11*(1-D11),
     &    E22*(1-D22),E33,G12*(1-D33),G23,G13
     &    ,CL_A_Layer)
C
      else 
          call get_C_mat_OfLayer_local_Aba(E11,E22
     &    ,E33,G12,G23,G13,nu12,nu23,nu13,CL_A_Layer)  !builds the stiffness matrix in local coord' of a layer (VUMAT notation)
         ! write(*,*) "CL_A",CL_A_layer
      end if
C      

c****************************************************************************     
c calculating the stress vector stressNew(6) at the material point i
c according to Abaqus VUMAT notation
c
       call Get_globalS_aba_no_rotation(CL_A_Layer,EG,SG)
c
	stressNew(i,1)=SG(1)
	stressNew(i,2)=SG(2)
	stressNew(i,3)=SG(3)
	stressNew(i,4)=SG(4)
	stressNew(i,5)=SG(5)
	stressNew(i,6)=SG(6)
   !    write(*,*)'stress=',sg
C
C***************************************************************************
C                           Failure calculations                           *
C***************************************************************************
C
C***************************************************************************        
C     Tsai Wu Calculations
C*************************************************************************** 
      phi=0
      I_max_TW=0
      Max_TW=0
      do j=1,6
          TW(J)=0
      end do
C
       call Get_TW(SG,TW,Phi,Max_TW,I_MAX_TW,f1,f2,f11,f22,f66,f12)
          statenew(i,8) = phi
          statenew(i,9) = I_MAX_TW
 !         statenew(i,19) = Max_TW
 !     statenew(i,20) = phi
          statenew(i,21) = TW(1)+TW(2)
          statenew(i,22) = TW(3)+TW(4)
          statenew(i,23) = TW(6)
 !         statenew(i,24) = TW(4)
 !         statenew(i,25) = TW(5)
 !         statenew(i,26) = TW(6)
   !   write(*,*) "TW right after exit of subroutine"
   !   write(*,*) Phi
C***************************************************************************
C     In plane tw based damage models 
C***************************************************************************
 
      StateNew(i,7) = StateOld(i,7)
 !     write(*,*) "stresses=", SG
 !     write(*,*) "Strains=", EG
 !     write(*,*)'phi=',phi
 !     write(*,*) "SDV17=", StateNew(i,17)
      if (phi.ge.1.05) then ! tw reached critical value. only the uncoupled parts of the stiffnes matrix stay and are degraded acording to the dominant mode.
          stateNew(i,7) = 1 ! degradation initiation flag
 
      end if
C         
      if ((phi.ge.1.05).and.(StateOld(i,7).gt.0.95)) then
          ! mode 1t initiation
          if ((I_max_TW.eq.1).and.
     &        (EG(1).gt.0).and.
     &        (StateOld(i,16).eq.0)) then
              StateNew(i,16) = EG(1)
              StateNew(i,17) = StateOld(i,17)
              StateNew(i,18) = StateOld(i,18)
              StateNew(i,19) = StateOld(i,19)
              StateNew(i,20) = StateOld(i,20)
C
          ! mode 1c initiation
          else if ((I_max_TW.eq.1).and.
     &        (EG(1).lt.0).and.
     &        (StateOld(i,17).eq.0)) then
              StateNew(i,16) = StateOld(i,16)
              StateNew(i,17) = EG(1)
              StateNew(i,18) = StateOld(i,18)
              StateNew(i,19) = StateOld(i,19)
              StateNew(i,20) = StateOld(i,20)
C
          ! mode 2t initiation
          else if ((I_max_TW.eq.2).and.
     &        (EG(2).gt.0).and.
     &        (StateOld(i,18).eq.0)) then
              StateNew(i,16) = StateOld(i,16)
              StateNew(i,17) = StateOld(i,17)
              StateNew(i,18) = EG(2)
              StateNew(i,19) = StateOld(i,19)
              StateNew(i,20) = StateOld(i,20)
C
          ! mode 2c initiation
          else if ((I_max_TW.eq.2).and.
     &        (EG(2).lt.0).and.
     &        (StateOld(i,19).eq.0)) then
              StateNew(i,16) = StateOld(i,16)
              StateNew(i,17) = StateOld(i,17)
              StateNew(i,18) = StateOld(i,18)
              StateNew(i,19) = EG(2)
              StateNew(i,20) = StateOld(i,20)
C
          ! mode 12 initiation
          else if ((I_max_TW.eq.3).and.
     &        (StateOld(i,20).eq.0)) then
              StateNew(i,16) = StateOld(i,16)
              StateNew(i,17) = StateOld(i,17)
              StateNew(i,18) = StateOld(i,18)
              StateNew(i,19) = StateOld(i,19)
              StateNew(i,20) = ABS(EG(4))
          else ! TW is gt 1, but no damage initiation occurred
               do i1 = 16,20
                  StateNew(i,i1) = StateOld(i,i1) 
               end do
          end if
      else    ! Tw is less than 1
          do i1 = 16,20
              StateNew(i,i1) = StateOld(i,i1) ! by this time, all damage initiation parameters were updated
          end do
      end if
C ************************************************************
C Updating damage parameters (D)
      ! fiber tension
      if ((EG(1).gt.0).and.(StateNew(i,16).gt.0)) then 
          StateNew(i,11) = 
     &     DMAX1(beta1t*(EG(1)-StateNew(i,16))/(EG(1)*(beta1t-1))
     &          ,StateOld(i,11))
          if (StateNew(i,11).GE.gamma1t) then 
              StateNew(i,11) = gamma1t !no deletion
              !StateNew(i,10) = 0 !deletion
          end if
      end if
C     ! fiber compression
      if ((EG(1).lt.0).and.(StateNew(i,17).lt.0)) then 
          StateNew(i,12) = 
     &     DMAX1(beta1c*(EG(1)-StateNew(i,17))/(EG(1)*(beta1c-1))
     &          ,StateOld(i,12))
          if (StateNew(i,12).GE.gamma1c) then 
              StateNew(i,12) = gamma1c !no deletion
              !StateNew(i,10) = 0 ! deletion
          end if
      end if
C     ! Matrix tension
      if ((EG(2).gt.0).and.(StateNew(i,18).gt.0)) then 
          StateNew(i,13) = 
     &     DMAX1(beta2t*(EG(2)-StateNew(i,18))/(EG(2)*(beta2t-1))
     &          ,StateOld(i,13))
          if (StateNew(i,13).GE.gamma2t) then 
              StateNew(i,13) = gamma2t !no deletion

              !StateNew(i,10) = 0 ! deletion
          end if
      end if
C     ! Matrix compression
      if ((EG(2).lt.0).and.(StateNew(i,19).lt.0)) then ! 
          StateNew(i,14) = 
     &     DMAX1(beta2c*(EG(2)-StateNew(i,19))/(EG(2)*(beta2c-1))
     &          ,StateOld(i,14))
          if (StateNew(i,14).GE.gamma2c) then 
              StateNew(i,14) = gamma2c !no deletion
              !StateNew(i,10) = 0 ! deletion
          end if
      end if
C     ! shear
      if ((StateNew(i,2).NE.0).and.(StateNew(i,20).GT.0)) then
          StateNew(i,15) = 
     &  DMAX1(beta12*(ABS(EG(4))-StateNew(i,20))/(ABS(EG(4))*(beta12-1))
     &          ,StateOld(i,15))
          if (StateNew(i,15).GE.gamma12) then 
              StateNew(i,15) = gamma12 !no deletion
              ! StateNew(i,10) = 0  ! deletion
          end if
      end if
C
      !update other Ds
      do j= 11,15
          if (stateNew(i,j)==0) then
              StateNew(i,j) = StateOld(i,j)
          end if
      end do
      
      
      
CCCCCCCCCCC CCCCCCCCCCCCCCCCCCCCCC CCCCCCCCCCCCCCCCCC CCCCCCCC CCCCCCCCCCCCCC         
!          if (I_max_Tw.eq.1) then ! fiber damage
!              if (StateOld(i,5).eq.0) then ! Max-stress strain needs to be defined
!                  StateNew(i,5) = EG(1)
!              else
!                  StateNew(i,5) = StateOld(i,5)
!              end if
!              StateNew(i,6) = StateOld(i,6)
!              StateNew(i,7) = StateOld(i,7)
!              if (EG(1).ge.0) then ! tension mode (probelm might rise if the mode changes form t to c)
!                  StateNew(i,1) = DMAX1(beta1t*(EG(1)-StateNew(i,5))/
!     &         (EG(1)*(beta1t-1)),StateOld(i,1))
!                  if (StateNew(i,1).GE.gamma1t) then 
!                      StateNew(i,10) = 0
!                  end if
!              else ! compression mode
!                  StateNew(i,1) = DMAX1(beta1c*(EG(1)-StateNew(i,5))/
c     &  (EG(1)*(beta1c-1)),StateOld(i,1))
c                  if (StateNew(i,1).GE.gamma1c) then
c                      StateNew(i,10) = 0
c                  end if
c              end if
cC
c !                write(*,*)'stat1=',StateNew(i,1)
c          else if (I_max_Tw.eq.2) then ! matrix damage
c              if (StateOld(i,6).eq.0) then ! Max-stress strain needs to be defined
c                  StateNew(i,6) = EG(2)
c              else
c                  StateNew(i,6) = StateOld(i,6)
c              end if
c              StateNew(i,5) = StateOld(i,5)
c              StateNew(i,7) = StateOld(i,7)
c              if (EG(2).ge.0) then ! tension mode (probelm might rise if the mode changes form t to c)
c                  StateNew(i,2) = DMAX1(beta2t*(EG(2)-StateNew(i,6))/
c     &         (EG(2)*(beta2t-1)),StateOld(i,2))
c                  if (StateNew(i,2).GE.gamma2t) then 
c                      StateNew(i,10) = 0
c                  end if
c              else ! compression mode
c                  StateNew(i,2) = DMAX1(beta2c*(EG(2)-StateNew(i,6))/
c     &  (EG(2)*(beta2c-1)),StateOld(i,2))
c                  if (StateNew(i,2).GE.gamma2c) then
c                      StateNew(i,10) = 0
c                  end if
c              end if
cC
c         else if (I_max_Tw.eq.3) then ! Sheer damage
c              if (StateOld(i,7).eq.0) then ! Max-stress strain needs to be defined
c                  StateNew(i,7) = DABS(EG(4))
c              else
c                  StateNew(i,7) = StateOld(i,7)
c              end if
c              StateNew(i,5) = StateOld(i,5)
c              StateNew(i,6) = StateOld(i,6)
c              StateNew(i,3) = DMAX1(beta12*(DABS(EG(4))-StateNew(i,7))/
c     &         (EG(4)*(beta12-1)),StateOld(i,3))
c                  if (StateNew(i,3).GE.gamma12) then 
c                      StateNew(i,10) = 0
c                  end if
c              end if
c          else
c  !       write(*,*) "error: TW value greater than 1 but no damage
  !   & mode was found"
  !            write(*,*) phi
  !            write(*,*) I_MAX_TW
c          end if

C***************************************************************************
  100 continue
                  
      return
      end
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
C
C                                 SUBROUTINES
C
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************      
c***************************************************************************
C     get_C_mat_OfLayer_local_Aba
c***************************************************************************
c  the following subroutine calculates the local stiffness matrix of 
c  an orthotropic material by givining the nine elastic properties according to the: 
c   Abaqus VUMAT notation [CL_A]
c
c
c         [  C11       C12       C13      0          0         0   ]
c         [  C21       C22       C23      0          0         0   ]
c [C]=    [  C31       C32       C33      0          0         0   ]
c         [   0         0         0       C44        0         0   ]
c         [   0         0         0       0         C55        0   ]
c         [   0         0         0       0          0         C66 ]
c  
c
c         [  CL_A11       CL_A12       CL_A13      0            0            0   ]
c         [  CL_A21       CL_A22       CL_A23      0            0            0   ]
c [CL_A]= [  CL_A31       CL_A32       CL_A33      0            0            0   ]
c         [   0            0            0         CL_A44        0            0   ]
c         [   0            0            0          0           CL_A55        0   ]
c         [   0            0            0          0            0         CL_A66 ]
c  get_C_mat_OfLayer_local_Aba(E11,E22,E33,G12,G23,G13,nu12,nu23,nu13,CL_A_Layer)
c
c  nu_ij/E_ii=nu_ji/E_jj
c
      subroutine get_C_mat_OfLayer_local_Aba(E11,E22,E33,G12,G23,G13,
     & nu12,nu23,nu13,CL_A)
	implicit real*8 (a-h,o-z)
	integer i,j
      real*8 C(6,6),CL_A(6,6) !  C: stiffness matrix of an orthotropic material
	real*8 const_1,nu21,nu32,nu31,G31,E11,E22,E33,G12,G23,G13,nu12,
     &  nu23,nu13
	do i=1,6
	   do j=1,6
	      C(i,j)=0
	      CL_A(i,j)=0
	   end do
	end do
C
	nu21=nu12*E22/E11
	nu32=nu23*E33/E22
	nu31=nu13*E33/E11
	G31=G13
	const_1=(1-nu12*nu21-nu23*nu32-nu31*nu13-2*nu12*nu23*nu31)
     &    /(E11*E22*E33)
C         
	CL_A(1,1)=(1-nu23*nu32)/(E22*E33*const_1)
	CL_A(1,2)=(nu21+nu31*nu23)/(E22*E33*const_1)
	CL_A(1,3)=(nu31+nu21*nu32)/(E22*E33*const_1)
	CL_A(2,1)=(nu12+nu13*nu32)/(E33*E11*const_1)
	CL_A(2,2)=(1-nu31*nu13)/(E33*E11*const_1)
	CL_A(2,3)=(nu32+nu31*nu12)/(E33*E11*const_1)
	CL_A(3,1)=(nu13+nu12*nu23)/(E11*E22*const_1)
	CL_A(3,2)=(nu23+nu13*nu21)/(E11*E22*const_1)
	CL_A(3,3)=(1-nu12*nu21)/(E11*E22*const_1)
	CL_A(4,4)=2*G12
	CL_A(5,5)=2*G23
	CL_A(6,6)=2*G31
	return
      end
c***************************************************************************
C     get_C_mat_OfLayer_local_Aba_damaged
c***************************************************************************
c  the following subroutine calculates the local stiffness matrix of 
c  an orthotropic material by givining the nine elastic properties according to the: 
c   Abaqus VUMAT notation [CL_A]
c
c
c         [  C11        0         0       0          0         0   ]
c         [  0         C22        0       0          0         0   ]
c [C]=    [  0          0         C33     0          0         0   ]
c         [   0         0         0       C44        0         0   ]
c         [   0         0         0       0         C55        0   ]
c         [   0         0         0       0          0         C66 ]
c  
c
c         [  CL_A11       CL_A12       CL_A13      0            0            0   ]
c         [  CL_A21       CL_A22       CL_A23      0            0            0   ]
c [CL_A]= [  CL_A31       CL_A32       CL_A33      0            0            0   ]
c         [   0            0            0         CL_A44        0            0   ]
c         [   0            0            0          0           CL_A55        0   ]
c         [   0            0            0          0            0         CL_A66 ]
c  get_C_mat_OfLayer_local_Aba(E11,E22,E33,G12,G23,G13,nu12,nu23,nu13,CL_A_Layer)
c
c  nu_ij/E_ii=nu_ji/E_jj
c
      subroutine get_C_mat_OfLayer_local_Aba_damaged(E11,E22,E33,G12,G23
     & ,G13,CL_A)
	implicit real*8 (a-h,o-z)
	integer i,j
      real*8 CL_A(6,6) !  C: stiffness matrix of an orthotropic material
	real*8 E11,E22,E33,G12,G23,G13
	do i=1,6
	   do j=1,6
	      CL_A(i,j)=0
	   end do
	end do
C         
	CL_A(1,1)=E11
	CL_A(1,2)=0
	CL_A(1,3)=0
	CL_A(2,1)=0
	CL_A(2,2)=E22
	CL_A(2,3)=0
	CL_A(3,1)=0
	CL_A(3,2)=0
	CL_A(3,3)=E33
	CL_A(4,4)=2*G12
	CL_A(5,5)=2*G23
	CL_A(6,6)=2*G13
	return
      end
c*************************************************************************
c the following subroutine calculates the stress vector SG(6) at 
C the material point according to Abaqus VUMAT notation
c
      subroutine Get_globalS_aba_no_rotation(eff_CG_A,EG,SG)
      implicit real*8 (a-h,o-z)
	real*8 eff_CG_A(6,6),SG(6),EG(6)
      SG(1)=eff_CG_A(1,1)*EG(1)+eff_CG_A(1,2)*EG(2)+eff_CG_A(1,3)*EG(3)

      SG(2)=eff_CG_A(2,1)*EG(1)+eff_CG_A(2,2)*EG(2)+eff_CG_A(2,3)*EG(3)

	SG(3)=eff_CG_A(3,1)*EG(1)+eff_CG_A(3,2)*EG(2)+eff_CG_A(3,3)*EG(3)

	SG(4)=eff_CG_A(4,4)*EG(4)

	SG(5)=eff_CG_A(5,5)*EG(5)
      SG(6)=eff_CG_A(6,6)*EG(6)
      return
      end
c*****************************************************************************
C     Tsai-Wu
C*****************************************************************************
c  the following subroutine calculates the "Modified Tsai-Wu" Value for 
C  a  given vector  of stresses in  abaqus  notation.  used  per  layer.
C
c "Modified Tsai-Wu" failure criterion for anisotropic materials
c "Modified Tsai-Wu" polinomial is:
c
c   according to Abaqus VUMAT notation:
c
c Phi=(f1*S1+f2*S2+f3*S3)+(f11*S1^2+f22*S2^2+f33*S3^2)+(f44*S31^2+f55*S12^2+
c    f66*S23^2)+(2*f12*S1*S2+2*f13*S1*S3+2*f23*S2*S3) 
c where: tau4=S31, tau5=S12, tau6=S23
c
      subroutine Get_TW(layer_local_stress_aba,TW_calculated_components,
     & Phi_TotalTW,max_component,max_component_index,f1,f2,f11,f22,
     & f66,f12)
      implicit real*8 (a-h,o-z)
      real*8 layer_local_stress_aba(6),TW_calculated_components(6)
      real*8 Phi_TotalTW,max_component,max_component_index
      real*8 f1,f2,f11,f22,f66,f12
      integer i1
C
      TW_calculated_components(1)=f1*layer_local_stress_aba(1)
      TW_calculated_components(2)=f11*layer_local_stress_aba(1)**2
      TW_calculated_components(3)=f2*layer_local_stress_aba(2)
      TW_calculated_components(4)=f22*layer_local_stress_aba(2)**2
      TW_calculated_components(5)=2*f12*layer_local_stress_aba(1)*
     & layer_local_stress_aba(2)
      TW_calculated_components(6)=f66*layer_local_stress_aba(4)**2
C   
	max_component=0 !max_component stores the maximal component of "Modified Tsai-Wu" polynomial
	Phi_TotalTW=0
      do i1=1,6
	 Phi_TotalTW=Phi_TotalTW+TW_calculated_components(i1)
  !     Write(*,*) "after adding one part of TW"
  !     write(*,*) Phi_TOTALTW
      end do
  !    write(*,*) "total TW inside subroutine"
      if ((TW_calculated_components(3)+TW_calculated_components(4)).gt.
     & (TW_calculated_components(1)+TW_calculated_components(2))) then ! 2>1
          if (TW_calculated_components(6).gt.
     & (TW_calculated_components(3)+TW_calculated_components(4))) then ! 3>2>1
              max_component_index = 3
              max_component = 
     &         TW_calculated_components(6)
          else ! 2>1>3 or 2>3>1
              max_component_index = 2
                           max_component = 
     &         TW_calculated_components(3)+TW_calculated_components(4)
          end if
      else if ((TW_calculated_components(1)+TW_calculated_components(2)! 1>2
     & ).gt.(TW_calculated_components(6))) then ! 1>2>3 or 1>3>2
          max_component_index = 1
                       max_component = 
     &         TW_calculated_components(1)+TW_calculated_components(2)
      else ! 3>1>2
          max_component_index = 3
           max_component = TW_calculated_components(6)
      endif
C
	return
      end