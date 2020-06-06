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
c
c*************************************************************************
c*  	Sub-laminate Model for Multi-layered Composites                  *                                
c*     N-layers model where the N+1 layer is cohesive interface          *                                             
c*                                                                       *
C* 			By Shlomo Spitzer, based on Yarden's code					 *
C*					'Version 3' 25.04.2020				   			     *	
c*   Abaqus Notation: for stress/strain (11,22,33,12,23,31)              *
c*   Sub-laminate notation (i,o) - (11,22,12)-(33,13,23)                 *                                                 
c*************************************************************************
!Version updates:

C "1" !Fixed the In Plane strain transformation subroutine
C "2" !Fixed the matrix compresion damage criterion
C "3" !Fixed the Layer C transformation subroutine
C "4" !Added a freeze to the In Plane damage when delamination starts
C "5" !fixed a mismatch in beta12 and gamma12 variable names
C "6" !fixed a bug where the "Update_layer_E_I" subroutine was called with the wrong parameters

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
     7  enerInelasOld(*), tempNew(*),
     8  stretchNew(*),
     8  defgradNew(*),
     9  fieldNew(*),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(*), enerInelasNew(*)
c
      character*80 cmname
c
c   Local Arrays/scalars Variables

      integer i1,i,k,j                                                ! indexing
      parameter N_layer=4 !<-------------------------------------------------------------FILL HERE
      ! Number of layers
		real*8 tk(N_layer),theta(N_layer),ti,t,N_layerXX				! Vector of thickness, Vector of orientation.
		real*8 E11,E22,E33,G23,G13,G12,nu23,nu13,nu12             ! Effective of one layer, 1 is the fiber direction. ti is interface layer's thickness & t is total thickness
		real*8 E11t,E11c,E22t,E22c,E33t,E33c
		real*8 X_t,X_c,Y_t,Y_c,S_i,S_o                                  ! Maximum stresses for a layer, by direction ( for the Tsai Wu criterion):
                                                                      ! Fiber tension, Fiber compression, Matrix tension, Matrix compression, In plain shear, Out of Plain shear
		real*8 f1,f2,f11,f22,f66,f12                                       ! Tsai Wu parameters, calculated from: X_t,X_c,Y_t,Y_c,S_i,S_o
		real*8 So_critical												! parameters that define interface degredation and delamination
		real*8 C_L_A_damaged(6,6),CL_A_Layer(6,6),CL_A_interface(6,6)	! holds the out plain multiplication value for the interface properties & layer properties
		real*8 CG_A_TEMP(6,6),CG_Ab(6*(N_layer+1),6)					! Global stiffness matrix
		real*8 CG_S_interface(6,6),CG_S_TEMP(6,6),CG_SL(6*(N_layer+1),6)! Global stiffness matrix
	                                                                ! according to Sublaminate model notation [CG_S] 
		real*8  eff_G_ABED(6,6),!G_ABED1(6,6),G_ABED2(6,6),G_ABED3(6,6), ! effective ABED matrix. according to Sublaminate model notation 
     &  G_ABED_TEMP(6,6),G_ABED_interface(6,6),G_ABED_SET(6*(N_layer+1),6)
		
		real*8 eff_CG_S(6,6),eff_CG_A(6,6)                              ! EFFECTIVE tangent stiffness matrix: Sublaminate notation, abaqus notation.
		real*8 EG_SET(6,(N_layer+1)),SG_SET(6,(N_layer+1))
		real*8 EL_SET(6,(N_layer+1)),SL_SET(6,(N_layer+1))				 ! strain & stress vectors SETS for each layer. G for global L for local.  abaqus notation
		real*8 EG(6),SG(6)                                              ! golbal strain & stress vectors. abaqus notation
		real*8 eff_G_EiSo(6)                                            ! global vector of IP strain and OOP stress
		real*8 EG_layer(6),SG_layer(6),EG_interface(6),SG_interface(6)  ! temp vectors 
		real*8 EL_layer(6),SL_layer(6),EL_interface(6),SL_interface(6)	! temp vectors 
		real*8 TW_value,i_Max_TW,Max_TW									! TW parameters

		real*8 E_IP_G(3),E_IP_L(3)
		real*8 alpha!, maximal OOP stress
		real*8 D11,D22,D12,temp_E_I,temp_E_atm,D_old,beta,gamma,D_new
		real*8 beta1t,beta1c,beta2t,beta2c,beta12,beta_dela,gamma_dela
		real*8 gamma1t,gamma1c,gamma2t,gamma2c,gamma12
      
		real*8 eff_OOP_stress,eff_OOP_strain,eff_OOP_strain_I
		real*8 local_layer_S_SL(6),local_layer_E_SL(6)
		real*8 D1t_new,D1c_new,D2t_new,D2c_new,D12_new,temp_D_new
		real*8 D1t_old,D1c_old,D2t_old,D2c_old,D12_old  
		real*8 E_I1t_old,E_I1c_old,E_I2t_old,E_I2c_old,E_I12_old
		real*8 E_I1t_new,E_I1c_new,E_I2t_new,E_I2c_new,E_I12_new
c*****************************************************************************     
c     The state variables 
c***************************************************************************** 
c StateNew(i,1)  == strain_aba 11
c StateNew(i,2)  == strain_aba 22
c StateNew(i,3)  == strain_aba 33
c StateNew(i,4)  == strain_aba 12
c StateNew(i,5)  == strain_aba 23
c StateNew(i,6)  == strain_aba 13
c StateNew(i,7)  == flag for degrading of material properties in the interface layer. holds 0 for no damage or the effective initialization strain eff_E_I
c StateNew(i,8)  == D of interface layer in tension
c StateNew(i,9)  == -empty-
c StateNew(i,10) == material point or elemnt deletion flag: 1=element exists; 0=element failed
C StateNew(i,11*k) == flag for initialization of progressive degradation of layer k
c StateNew(i,1+11*k) == D1t of layer k 
C StateNew(i,2+11*k) == D1c of layer k
C StateNew(i,3+11*k) == D2t of layer k
C StateNew(i,4+11*k) == D2c of layer k
C StateNew(i,5+11*k) == D12 of layer k
C StateNew(i,6+11*k) == E11_I_1t of layer k  !17 28 39 50
C StateNew(i,7+11*k) == E11_I_1c of layer k  !18 29 40 51
C StateNew(i,8+11*k) == E22_I_2t of layer k  !19 30 41 52
C StateNew(i,9+11*k) == E22_I_2c of layer k  !20 31 42 53
C StateNew(i,10+11*k) == E12_I_12 of layer k !21 32 43 54
c*****************************************************************************
C******************************************************************************
C******************************************************************************
C     Define Geometry
C     
C     t(k)     is the k layer's relative thickness and is being automatically computed
C     ti       is the interface's relative thickness and can be left as a small constant
C     theta(k) should be defined by user in the next lines. In the future as VUMAT parameters.
C******************************************************************************

      N_layerXX = real(N_layer)
 


    ! Orientation of each layer in the model (Order is not important) !<--------------------------------------Update orientations
      theta(1) = 90 
      tk(1) = 1
C
      theta(2) = 0
      tk(2) = 3
C
      theta(3) = 45
      tk(3) = 1
C
      theta(4) = -45
      tk(4) = 1
C
      t = 0
      do K=1,N_layer
          t = t+tk(k)
      end do
		ti= t*1.d-10 !(E-10)
      t = t + ti
C******************************************************************************
C******************************************************************************
c     Define Material
C******************************************************************************
C     Make sure that stresses are represented in units of [Pa]. This is the case
C     when the lengths are in meters in the analysis.
C
      ! Orthotropic layer's properties: MI7 !(T300/NY9200Z) <------------------------------------------------------------------------------FILL HERE
        G12 = 5677.4d6 ! 4160d6   !5.3e9 ! [Pa]
        G23 = 3141.8d6 ! 3110d6   !3.8e9 ! [Pa]
        G13 = 5215.7d6 ! 4160d6   !5.3e9 ! [Pa]
C
        nu12=0.32      !.31
        nu23=0.461     !.3
        nu13=0.329     !.31
C
        E11t = 158d9
        E22t = 8.97d9
        E33t = 8.89d9
C
        E11c = 158d9!130d9
        E22c = 8.97d9!8.69d9
        E33c = 8.89d9!8.63d9
c
      ! critical stress value: IM7/977-3 ! IM7 !(T300/NY9200Z)
        X_t = 2704d6   !2764.3d6  !1747e6    ! [Pa]
        X_c = 1764.3d6 !2150.01d6 !1357e6    ! [Pa]
        Y_t = 95.7d6   !79.75d6   !67e6 !78  ! [Pa]
        Y_c = 236.5d6  !249.08d6  !170e6     ! [Pa]
        S_i = 117.9d6  ! 260d6               ! [Pa]
        S_o = 96.5d6   !98.19d6   !124e6     ! [Pa]
c		
        So_critical= 0.8*64d6!1*Y_t ! compared to the effective stress OOP in the Interface layer. Yarden used 1.3*Y_t
c
      ! the material parameters of modified "Tsai-Wu" failure criteria, 1 is the fiber direction:
		f1=1/X_t-1/X_c
		f2=1/Y_t-1/Y_c
        f11=1/(X_t*X_c)
		f22=1/(Y_t*Y_c)
		f66=1/(S_i)**2
		f12=-0.5*dsqrt(DABS(f11*f22))
	    ! progressive damage parameters. # corresponds to damage mode#. t for tension and c for compression. 
        ! the value represents how many times the maximal strain is larger than the strain at damage initiation, assuming the element is not deleted on the way.  
        beta1t = 1.5
        beta1c = 1.5
        beta2t = 2
        beta2c = 2.
        beta12 = 2.
        beta_dela = 2
        ! Gamma parameters represent the maximal damage progression allowed in each mode
        gamma1t = 0.8
        gamma1c = 0.8
        gamma2t = 0.8
        gamma2c = 0.8
        gamma12 = 0.8
        gamma_dela = 1.-1.D-10
C
        alpha = 0.28!2 ! determines how dominant are the shear strains in comparion to the tensile ones for delamination determination
c****************************************************************************
C
        do 100 i= 1,nblock 
c  
*****************************************************************************
C
c******************************************************************************

c calculating the strain vector EG(6) at the material point i
c according to Abaqus VUMAT notation 
c
       EG(1) = stateOld(i,1) + strainInc(i,1)
       EG(2) = stateOld(i,2) + strainInc(i,2)
	   EG(3) = stateOld(i,3) + strainInc(i,3)
	   EG(4) = stateOld(i,4) + strainInc(i,4)
	   EG(5) = stateOld(i,5) + strainInc(i,5)
	   EG(6) = stateOld(i,6) + strainInc(i,6)
C
      stateNew(i,1)=EG(1)
      stateNew(i,2)=EG(2)
	  stateNew(i,3)=EG(3)
	  stateNew(i,4)=EG(4)
	  stateNew(i,5)=EG(5)
	  stateNew(i,6)=EG(6)
C
C ! building the 2D IP strain vector for determining fiber and matrix mode
      E_IP_G(1) = EG(1)
      E_IP_G(2) = EG(2)
      E_IP_G(3) = EG(4)

C
c******************************************************************************

c****************************************************************************
c calculating the Global stiffness matrix [CG] of the orthotropic layers
c of the sublaminate model:according to Abaqus VUMAT notation [CG_A]
c
		do k=1,6*(N_layer+1) !initializing CG_ab
			do j=1,6
			CG_Ab(k,j)=0
			end do
		end do
c	 
		do k=1,N_layer ! put all CG_A matrices in a set [6*N_layer,6]
			call Get_local_IP_strain(theta(k),E_IP_G,E_IP_L) ! check local IP strain for tension or compression mode, IP global strain is equal in every layeer

			if (E_IP_L(1).lt.0) then
                  E11 = E11c
			else
					  E11 = E11t
			end if
C
			if (E_IP_L(2).lt.0) then
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

          if (StateOld(i,11*k).ge.0.9) then ! This layer started degradation 
                  !Determine relevant damage parameter
              if (E_IP_L(1).gt.0) then ! fiber tension
                      D11 = StateOld(i,1+11*k) 
              else					 ! fiber compression
                      D11 = StateOld(i,2+11*k) 
              end if
C
              if (E_IP_L(2).gt.0) then ! Matrix tension
                      D22 = StateOld(i,3+11*k)
              else					 ! Matrix compression
                      D22 = StateOld(i,4+11*k) 
              end if 
C					! one mode for shear
					D12 = StateOld(i,5+11*k)
C! Create damaged local layer C matrix, aba notation

              call Get_local_aba_C_matrix_prog_damage(E11,E22,E33,G12,
     &                                G23,G13,D11,D22,D12,C_L_A_damaged)

C! rotate the C matrix  to global layer C matrix, aba notation
              call get_C_mat_ofLayer_global_Aba(theta(k),
     &                            C_L_A_damaged,CG_A_TEMP)
              CG_Ab(k*6-5:k*6,1:6)=CG_A_TEMP
C
          else ! No degradation
C
			call get_C_mat_OfLayer_local_Aba(E11,E22,E33,G12,G23,G13
     &		,nu12,nu23,nu13,CL_A_Layer)  !builds the stiffness matrix in local coord' of a layer (VUMAT notation)
		    call get_C_mat_ofLayer_global_Aba(theta(k),  ! if no damage happend, take the general layer C and rotate
     &        CL_A_Layer,CG_A_TEMP)
		    CG_Ab(k*6-5:k*6,1:6)=CG_A_TEMP   

          end if
      end do
c 
      if (stateOld(i,7).gt.0) then ! damage initialized in the OOP     
          if (EG(3).gt.0) then ! delamination
              temp_D_new = stateOld(i,8)
          else                 ! Crushing
              temp_D_new = 0!temp_D_new = StateOld(i,9)
          end if
          call get_interface_layer_C(E33t,temp_D_new,CL_A_interface)

      else ! no damage to interface
          temp_D_new = 0
          call get_interface_layer_C(E33t,temp_D_new,CL_A_interface)
      end if
      CG_Ab((N_layer+1)*6-5:(N_layer+1)*6,1:6)= CL_A_interface ! add to the set of global C matrices abaqus notation

c****************************************************************************
c calculating the EFFECTIVE G matrix (eff_G_ABED) of the Global [G_ABED) matrix for 
c the orthotropic layers of the sublaminate model according to Sublaminate 
c model notation [CG_S]
c
c
	do k=1,6*(N_layer+1) !initializing CG_SL
			do j=1,6
			CG_SL(k,j)=0
			end do
      end do
C		
	do k=1,N_layer ! convert the abaqus notation to SL notation and put all CG_SL matrices in a set [6*N_layer,6]
		call get_CG_S(CG_ab(k*6-5:k*6,1:6),CG_S_TEMP)
		CG_SL(k*6-5:k*6,1:6)=CG_S_TEMP
	end do
C		
		call get_CG_S(CG_ab((N_layer+1)*6-5:(N_layer+1)*6,1:6),
     & CG_S_interface) !	put CG_S_interface
		CG_SL((N_layer+1)*6-5:(N_layer+1)*6,1:6)=CG_S_interface
C		
	do k=1,6*(N_layer+1) !initializing G_ABED_SET
			do j=1,6
			G_ABED_SET(k,j)=0
			end do
      end do
C		
	do k=1,N_layer ! put all G_ABED matrices in a set [6*N_layer,6]
		call get_G_ABED(CG_SL(k*6-5:k*6,1:6),G_ABED_TEMP)
		G_ABED_SET(k*6-5:k*6,1:6)=G_ABED_TEMP
	end do
C		
      call get_G_ABED(CG_SL((N_layer+1)*6-5:(N_layer+1)*6,1:6),
     &   G_ABED_interface) !	put G_ABED_interface
		G_ABED_SET((N_layer+1)*6-5:(N_layer+1)*6,1:6)=G_ABED_interface
C****************************************************************************		
C    Combine the Global ABED matrices to the effective one of the SL
	do k=1,6 !initializing eff_G_ABED
			do j=1,6
			eff_G_ABED(k,j)=0
			end do
      end do
c 
	do k=1,N_layer
		eff_G_ABED=eff_G_ABED+(tk(k)/t)*G_ABED_SET(k*6-5:k*6,1:6)
	end do
c
       eff_G_ABED=eff_G_ABED+(ti/t)*G_ABED_SET((N_layer+1)*6-5:
     & (N_layer+1)*6,1:6)
c
c****************************************************************************
c calculating the EFFECTIVE tangent stiffness matrix [eff_CG_S] of 
c the orthotropic layers of the sublaminate model according to Sublaminate 
c model notation [CG_S]
c (partial invertion back to global C matrix in SL notation)
      call get_eff_CG_S(eff_G_ABED,eff_CG_S)	  
c	
c****************************************************************************
c calculating the EFFECTIVE tangent stiffness matrix [eff_CG_A] of 
c the orthotropic layers of the sublaminate model according to 
c Abaqus VUMAT notation [CG_A]
c (convert the global C matrix to abaqus notation
      call get_eff_CG_A(eff_CG_S,eff_CG_A)
c     
c****************************************************************************     
c calculating the stress vector stressNew(6) at the material point i
c according to Abaqus VUMAT notation
c
      call get_G_stress_A(eff_CG_A,EG,SG)
c
	stressNew(i,1)=SG(1)
	stressNew(i,2)=SG(2)
	stressNew(i,3)=SG(3)
	stressNew(i,4)=SG(4)
	stressNew(i,5)=SG(5)
	stressNew(i,6)=SG(6)
C****************************************************************************
C delamination determination - should result in updated sdv 7 & 8
      
C     
      if (EG(3).lt.0) then ! If OOP in compression, no degradation is allowed
          StateNew(i,7) = StateOld(i,7)
          StateNew(i,8) = StateOld(i,8)
      else ! else it's in tension, check if should initialize or progress
C
          eff_OOP_strain_I = StateOld(i,7) ! delamination initiation strain
C
          if (eff_OOP_strain_I.eq.0) then !if no delamination yet happened, check if needed to be started
              eff_OOP_stress = DSQRT(SG(3)**2+alpha*(SG(5)**2+SG(6)**2))
              if (eff_OOP_stress.gt.So_critical) then ! if delamination criterion was met, initiate delamination
                  eff_OOP_strain_I = 
     &                        DSQRT(EG(3)**2+alpha*(EG(5)**2+EG(6)**2))
                  StateNew(i,7) = eff_OOP_strain_I
                  StateNew(i,8) = StateOld(i,8)
              else ! delamination didn't start and shouldn't yet, keep SDVs
                  StateNew(i,7) = StateOld(i,7)
                  StateNew(i,8) = StateOld(i,8)
              end if   
C	
          else ! else, delamination was already initiated and eff_OOP_strain_I.ne.0, progress delamination
C
			  StateNew(i,7) = Stateold(i,7) ! keep old initation strain
              eff_OOP_strain=DSQRT(EG(3)**2+alpha*(EG(5)**2+EG(6)**2)) !Current effective Strain
C
				if (eff_OOP_strain.gt.eff_OOP_strain_I) then ! should progress the damage if strain reached a new maximum
					D_old = StateOld(i,8)   !previous delamination D
					call get_updated_D(eff_OOP_strain_I,eff_OOP_strain,
     &                               D_old,beta_dela,gamma_dela,D_new) ! Update D
					StateNew(i,8) = D_new
				else ! delamination initiated but strain is lower then the initialization one
					StateNew(i,8) = StateOld(i,8)
				end if
C
          end if
      end if
!statenew(i,9)  = DSQRT(EG(3)**2+alpha*(EG(5)**2+EG(6)**2)) ! temp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
c**************************************************************************** 
c calculating the stress & strain EG, SG vectors for each 
c layer of the sublaminate model according to Sublaminate model notation and GLOBAL
c coordinate system
c
c    {G_SiEo(k)}=[G_ABED(k)]*{eff_G_EiSo}
c where: {eff_G_EiSo}={EG(1),EG(2),EG(4),SG(3),SG(6),SG(5)} (8.37)
c
		eff_G_EiSo(1)=EG(1)
		eff_G_EiSo(2)=EG(2) 
		eff_G_EiSo(3)=EG(4) 
		eff_G_EiSo(4)=SG(3)
		eff_G_EiSo(5)=SG(5)
		eff_G_EiSo(6)=SG(6)
c	
		do k=1,6 !initializing Global SG&EG sets 
          do j=1,(N_layer+1)
			EG_SET(k,j)=0
			SG_SET(k,j)=0
          end do
		end do
c		
	do k=1,N_layer ! put all G_E_S vectors in EG,SG sets [6,N_layer+1]

        call get_G_E_S_vectors(G_ABED_SET(k*6-5:k*6,1:6),eff_G_EiSo,
     &            EG_layer,SG_layer)
		EG_SET(1:6,k)=EG_layer
		SG_SET(1:6,k)=SG_layer
	end do
c		
		call get_G_E_S_vectors(G_ABED_SET((N_layer+1)*6-5:(N_layer+1)*6,1:6),
     &   eff_G_EiSo,EG_interface,SG_interface) !	put SG_interface, EG_interface
		EG_SET(1:6,(N_layer+1))=EG_interface
		SG_SET(1:6,(N_layer+1))=SG_interface
c      
c****************************************************************************
c calculating the stress & strain EL, SL vectors for each
c layer of the sublaminate model according to Sublaminate model notation and LOCAL
c coordinate system
c 
      do k=1,6 !initializing Global SL&EL sets 
          do j=1,(N_layer+1)
			EL_SET(k,j)=0
			SL_SET(k,j)=0
          end do
      end do
c		
	do k=1,N_layer ! put all G_E_L vectors in EL,SL sets [6,N_layer+1] abaqus notation
		call get_L_E_S_vectors(theta(k),EG_SET(1:6,k),SG_SET(1:6,k),EL_layer,
     & SL_layer)
		EL_SET(1:6,k)=EL_layer
		SL_SET(1:6,k)=SL_layer
 
	end do
c		
		call get_L_E_S_vectors(0,EG_SET(1:6,(N_layer+1)),SG_SET(1:6,
     &  (N_layer+1)),EL_interface,SL_interface) !	put SL_interface, EL_interface
		EL_SET(1:6,(N_layer+1))=EL_interface
		SL_SET(1:6,(N_layer+1))=SL_interface
c		
C***************************************************************************
C                       End of Stess calculations                          *
C***************************************************************************
C
! TEMP
 !     Write(*,*) " "
 !     write(*,*) "CL_A"
!      write(*,*) CL_A_LAYER
 !     write(*,*) "CG_A_LAYER1"
!      write(*,*) CG_Ab(1:6,1:6)
 !     write(*,*) "CG_A_LAYER2"
!      write(*,*) CG_Ab(7:12,1:6)
 !     write(*,*) "CG_A_LAYER3"
 !     write(*,*) CG_Ab(13:18,1:6)
 !     write(*,*) "CG_A_LAYER4"
!      write(*,*) CG_Ab(19:24,1:6)
!1      write(*,*) "CG_SL_LAYER1"
!      write(*,*) CG_SL(1:6,1:6)
!      write(*,*) "CG_SL_LAYER2"
!      write(*,*) CG_sl(7:12,1:6)
!      write(*,*) "CG_SL_LAYER3"
!      write(*,*) CG_sl(13:18,1:6)
!      write(*,*) "CG_SL_LAYER4"
!      write(*,*) CG_sl(19:24,1:6)
!      write(*,*) "G_ABED_SET_layer1"
!      write(*,*) G_ABED_SET(1:6,1:6)
!      write(*,*) "G_ABED_SET_layer2"
!      write(*,*) G_ABED_SET(7:12,1:6)!
!      write(*,*) "G_ABED_SET_layer3"
!      write(*,*) G_ABED_SET(13:18,1:6)
!      write(*,*) "G_ABED_SET_layer4"
!      write(*,*) G_ABED_SET(19:24,1:6)
!      write(*,*) "eff_G_ABED"
!      write(*,*) eff_G_ABED
!      write(*,*) "eff_CG_S"
!     write(*,*) eff_CG_S
!      write(*,*) "effective C Global abaqus"
!      write(*,*) eff_CG_A    
!      write(*,*) "global strain"
!      write(*,*) EG    !
!	  write(*,*) "global stress"
!      write(*,*) SG    
!	  write(*,*) " "
!      write(*,*) " "  
C***************************************************************************
C***************************************************************************
C                          Layer Failure calculations                      *
C***************************************************************************
C***************************************************************************
c
C for each layer in the sublaminate, calculate TW_polynomial and update damage
C 
C If delamaination was initiated, the other damage modes are frozen, 
C this is controled by SDV7, the flag for delamination
C
	if (StateNew(i,7).gt.0.000001) then ! if delamination started, freeze In plane damage
		do k=1,N_layer ! per layer k
			do j=1,6 ! per damage parameter
				StateNew(i,11*k+j-1)=Stateold(i,11*k+j-1)
			end do
		end do
	else ! There was no delamination - progress In plane damage
	
      do k=1,N_layer ! per layer k
          TW_value=0
          I_max_TW=0
          Max_TW=0
          local_layer_S_SL = SL_SET(1:6,k) !create the stress vector of the layer, SL, local
          local_layer_e_SL = EL_SET(1:6,k) !create the strain vector of the layer, SL, local
          call Get_TW(local_layer_S_SL,TW_value,Max_TW,I_MAX_TW,
     &         f1,f2,f11,f22,f66,f12)
C keep old degradation flag
          StateNew(i,11*k) = StateOld(i,11*k)
          if (TW_value.ge.1.05) then ! tw reached critical value. only the uncoupled parts of the stiffnes matrix stay and are degraded acording to the dominant mode.
              StateNew(i,11*k) = 1 ! degradation initiation flag
 !             write(*,*) "degradation initiation flag"
          end if
C prepare parameters for subroutine
          if (StateOld(i,11*k).gt.0.95) then ! degradation started already, check for E_I updates
!			 write(*,*) "SDV's","",stateOld(i,1+11*k),stateOld(i,2+11*k)
!			 write(*,*) stateOld(i,3+11*k),stateOld(i,4+11*k),""
c
              E_I1t_old = StateOld(i,6+11*k)
              E_I1c_old = StateOld(i,7+11*k)
              E_I2t_old = StateOld(i,8+11*k)
              E_I2c_old = StateOld(i,9+11*k)
              E_I12_old = StateOld(i,10+11*k)
C
              call Update_layer_E_I(local_layer_e_SL,I_max_TW,
     &			E_I1t_old,E_I1c_old,E_I2t_old,E_I2c_old,E_I12_old,
     &			E_I1t_new,E_I1c_new,E_I2t_new,E_I2c_new,E_I12_new)  
C get new E_I's
			StateNew(i,6+11*k) = E_I1t_new
			StateNew(i,7+11*k) = E_I1c_new
			StateNew(i,8+11*k) = E_I2t_new
			StateNew(i,9+11*k) = E_I2c_new
			StateNew(i,10+11*k) = E_I12_new
          else    ! degradation did not start
              do i1 = 6,10
                  StateNew(i,i1+11*k) = StateOld(i,i1+11*k) ! just keep old values of E_I's
              end do
          end if
C ! update D's
          if (StateOld(i,11*k).lt.0.95) then ! no degradation started
              do i1 = 1,5
                  StateNew(i,i1+11*k) = StateOld(i,i1+11*k)
              end do
          else ! degradation started, check each mode and update damage
		    D1t_old = StateOld(i,1+11*k)
		    D1c_old = StateOld(i,2+11*k)
		    D2t_old = StateOld(i,3+11*k)
		    D2c_old = StateOld(i,4+11*k)
		    D12_old = StateOld(i,5+11*k)
      if ((local_layer_e_SL(1).gt.0).and.(E_I1t_new.gt.1E-10)) then 
                  E_atm = local_layer_e_SL(1)
		        call get_updated_D(E_I1t_new,E_atm,D1t_old,beta1t,gamma1t,
     &             D1t_new)
                  StateNew(i,1+11*k) = D1t_new
              end if
C
      if ((local_layer_e_SL(1).lt.0).and.(E_I1c_new.lt.-1E-10)) then 
                  E_atm = local_layer_e_SL(1)
		        call get_updated_D(E_I1c_new,E_atm,D1c_old,beta1c,gamma1c,
     &             D1c_new)
                  StateNew(i,2+11*k) = D1c_new
              end if
C
      if ((local_layer_e_SL(2).gt.0).and.(E_I2t_new.gt.1E-10)) then 
                  E_atm = local_layer_e_SL(2)
		        call get_updated_D(E_I2t_new,E_atm,D2t_old,beta2t,gamma2t,
     &             D2t_new)
                  StateNew(i,3+11*k) = D2t_new
              end if
C
      if ((local_layer_e_SL(2).lt.0).and.(E_I2c_new.lt.-1E-10)) then 
                  E_atm = local_layer_e_SL(2)
		        call get_updated_D(E_I2c_new,E_atm,D2c_old,beta2c,gamma2c,
     &             D2c_new)
                  StateNew(i,4+11*k) = D2c_new
              end if
C
      if (E_I12_new.gt.1E-10) then 
	  	!			  write(*,*) "E_I12_new"
				  !write(*,*) E_I12_new
                  E_atm = DABS(local_layer_e_SL(3))
		        call get_updated_D(E_I12_new,E_atm,D12_old,beta12,gamma12,
     &              D12_new)
                  StateNew(i,5+11*k) = D12_new 
		!		  write(*,*) "D12_new"
		!		  write(*,*) D12_new
              end if
C update other D's which did not progress
              do j= 1,5
                  if (stateNew(i,j+11*k)==0) then
                      StateNew(i,j+11*k) = StateOld(i,j+11*k)
                  end if
              end do
          end if    
      end do
	  end if
          
C***************************************************************************

  100 continue
 !     write(*,*)"exit VUMAT"
	  
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
C
c*****************************************************************************
C     Tsai-Wu
C*****************************************************************************
c  the following subroutine calculates the "Tsai-Wu" Value for 
C  a  given vector of local stresses in Sublaminate  notation. used per layer.

      subroutine Get_TW(layer_local_stress_SL,
     & Phi_TotalTW,max_component,max_component_index,f1,f2,f11,f22,
     & f66,f12)
      implicit real*8 (a-h,o-z)
      real*8 layer_local_stress_SL(6),TW_calculated_components(6)
      real*8 Phi_TotalTW,max_component,max_component_index
      real*8 f1,f2,f11,f22,f66,f12
      integer i1
C
      TW_calculated_components(1)=f1*layer_local_stress_SL(1)
      TW_calculated_components(2)=f11*layer_local_stress_SL(1)**2
      TW_calculated_components(3)=f2*layer_local_stress_SL(2)
      TW_calculated_components(4)=f22*layer_local_stress_SL(2)**2
      TW_calculated_components(5)=2*f12*layer_local_stress_SL(1)*
     &                            layer_local_stress_SL(2)
      TW_calculated_components(6)=f66*layer_local_stress_SL(3)**2
C   
	max_component = 0 !max_component stores the maximal component of "Tsai-Wu" polynomial, unused unless needed for debugging
      max_component_index = 0
	Phi_TotalTW=0
      do i1=1,6
          Phi_TotalTW = Phi_TotalTW+TW_calculated_components(i1)
      end do
c Determine the dominant mode of the TW polynomiyal
      if ((TW_calculated_components(3)+TW_calculated_components(4)).gt.
     & (TW_calculated_components(1)+TW_calculated_components(2))) then ! 2>1
          if (TW_calculated_components(6).gt.
     & (TW_calculated_components(3)+TW_calculated_components(4))) then ! 3>2>1
              max_component_index = 3
              max_component = TW_calculated_components(6)
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
C     get_C_mat_ofLayer_global_Aba
c*************************************************************************** 
c  the following subroutine calculates the Global stiffness matrix [GC] 
c  for an orthotropic layer by giving the local stiffness matrix and its angle
c  refered to the global axes according to Abaqus VUMAT notation [CG_A]
c
c         [  CG11      CG12      CG13      CG14       0         0   ]
c         [  CG21      CG22      CG23      CG24       0         0   ]
c [CG_A]= [  CG31      CG32      CG33      CG34       0         0   ]
c         [  CG41      CG42      CG43      CG44       0         0   ]
c         [   0         0         0         0        CG55      CG56 ]
c         [   0         0         0         0        CG65      CG66 ]
      subroutine get_C_mat_ofLayer_global_Aba(theta,CL_A,CG_A)
	implicit real*8 (a-h,o-z)
	real*8 theta,CG_A(6,6),CL_A(6,6),T(6,6),TI(6,6),theta_radian
	real*8 c,s,pi,a,b
	
      if (theta.eq.0) then 
          do i=1,6
	        do j=1,6
	            CG_A(i,j)=CL_A(i,j)
	        end do
          end do
          goto 99
      elseif (theta.eq.90) then
          c=0
          s=1
      else
          pi=3.1415926535897932
! Since the transformation is back to global coord. The MINUS is added here. 
	    theta_radian=  -theta*pi/180
	    c=dcos(theta_radian)
          s=dsin(theta_radian)
      end if
      do i=1,6
	 do j=1,6
	  CG_A(i,j)=0
	 end do
      end do
	  
	  ! The following lines calculates - CG_A=matmul(T,matmul(CL_A,transpose(T)))
	
		 CG_A(1,1)=c**2*(CL_A(1,1)*c**2 + CL_A(2,1)*s**2 + 2*CL_A(4,1)*c*s) +
     &   s**2*(CL_A(1,2)*c**2 + CL_A(2,2)*s**2 + 2*CL_A(4,2)*c*s) + 
     &  2*c*s*(CL_A(1,4)*c**2 + CL_A(2,4)*s**2 + 2*CL_A(4,4)*c*s) 

		 CG_A(1,2)=c**2*(CL_A(1,2)*c**2 + CL_A(2,2)*s**2 + 2*CL_A(4,2)*c*s) +
     &   s**2*(CL_A(1,1)*c**2 + CL_A(2,1)*s**2 + 2*CL_A(4,1)*c*s) -
     & 2*c*s*(CL_A(1,4)*c**2 + CL_A(2,4)*s**2 + 2*CL_A(4,4)*c*s)

		 CG_A(1,3)=CL_A(1,3)*c**2 + CL_A(2,3)*s**2 + 2*CL_A(4,3)*c*s

		 CG_A(1,4)=(c**2 - s**2)*(CL_A(1,4)*c**2 + CL_A(2,4)*s**2 +
     & 2*CL_A(4,4)*c*s) - c*s*(CL_A(1,1)*c**2 + CL_A(2,1)*s**2 + 
     & 2*CL_A(4,1)*c*s)+c*s*(CL_A(1,2)*c**2+CL_A(2,2)*s**2+2*CL_A(4,2)
     & *c*s)

	 CG_A(2,1)=c**2*(CL_A(2,1)*c**2 + CL_A(1,1)*s**2 - 2*CL_A(4,1)*c*s) + 
     & s**2*(CL_A(2,2)*c**2 + CL_A(1,2)*s**2 - 2*CL_A(4,2)*c*s) +
     & 2*c*s*(CL_A(2,4)*c**2 + CL_A(1,4)*s**2 - 2*CL_A(4,4)*c*s)

	 CG_A(2,2)=c**2*(CL_A(2,2)*c**2 + CL_A(1,2)*s**2 - 2*CL_A(4,2)*c*s) +
     &  s**2*(CL_A(2,1)*c**2 + CL_A(1,1)*s**2 - 2*CL_A(4,1)*c*s) -
     & 2*c*s*(CL_A(2,4)*c**2 + CL_A(1,4)*s**2 - 2*CL_A(4,4)*c*s)

		 CG_A(2,3)=CL_A(2,3)*c**2 + CL_A(1,3)*s**2 - 2*CL_A(4,3)*c*s
		 
		 CG_A(2,4)=(c**2 - s**2)*(CL_A(2,4)*c**2 + CL_A(1,4)*s**2 -
     & 2*CL_A(4,4)*c*s) - c*s*(CL_A(2,1)*c**2 + CL_A(1,1)*s**2 -
     & 2*CL_A(4,1)*c*s) + c*s*(CL_A(2,2)*c**2 + CL_A(1,2)*s**2 -
     & 2*CL_A(4,2)*c*s)
		 
		 CG_A(3,1)=CL_A(3,1)*c**2 + CL_A(3,2)*s**2 + 2*CL_A(3,4)*c*s
		 
		 CG_A(3,2)=CL_A(3,2)*c**2 + CL_A(3,1)*s**2 - 2*CL_A(3,4)*c*s
		 
		 CG_A(3,3)=CL_A(3,3)
		 
		 CG_A(3,4)=CL_A(3,4)*(c**2 - s**2) - CL_A(3,1)*c*s + CL_A(3,2)*c*s
		 
		 CG_A(4,1)=c**2*(CL_A(4,1)*(c**2 - s**2) - CL_A(1,1)*c*s +
     & CL_A(2,1)*c*s) + s**2*(CL_A(4,2)*(c**2 - s**2) - CL_A(1,2)*c*s +
     & CL_A(2,2)*c*s) + 2*c*s*(CL_A(4,4)*(c**2 - s**2) - CL_A(1,4)*c*s
     & + CL_A(2,4)*c*s)

		 CG_A(4,2)=c**2*(CL_A(4,2)*(c**2 - s**2) - CL_A(1,2)*c*s +
     & CL_A(2,2)*c*s) + s**2*(CL_A(4,1)*(c**2 - s**2) - CL_A(1,1)*c*s +
     & CL_A(2,1)*c*s) - 2*c*s*(CL_A(4,4)*(c**2 - s**2) - CL_A(1,4)*c*s
     & + CL_A(2,4)*c*s)
		 
		 CG_A(4,3)=CL_A(4,3)*(c**2 - s**2) - CL_A(1,3)*c*s + CL_A(2,3)*c*s

		 CG_A(4,4)=(c**2 - s**2)*(CL_A(4,4)*(c**2 - s**2) - CL_A(1,4)*c*s +
     & CL_A(2,4)*c*s) - c*s*(CL_A(4,1)*(c**2 - s**2) - CL_A(1,1)*c*s +
     & CL_A(2,1)*c*s) + c*s*(CL_A(4,2)*(c**2 - s**2) - CL_A(1,2)*c*s +
     & CL_A(2,2)*c*s)

		 CG_A(5,5)=c*(CL_A(5,5)*c - CL_A(6,5)*s)-s*(CL_A(5,6)*c-CL_A(6,6)*s)

		 CG_A(5,6)=c*(CL_A(5,6)*c - CL_A(6,6)*s)+s*(CL_A(5,5)*c-CL_A(6,5)*s)
		 
		 CG_A(6,5)=c*(CL_A(6,5)*c + CL_A(5,5)*s)-s*(CL_A(6,6)*c+ CL_A(5,6)*s)

		 CG_A(6,6)=c*(CL_A(6,6)*c + CL_A(5,6)*s)+ s*(CL_A(6,5)*c+CL_A(5,5)*s)
99    return
	end

c***************************************************************************
C     probabaly not needed anymore
c  the following subroutine calculates the inverse matrix of 
c  the stress transformation matrix [T] of (6*6)
c*
c*  c=cos(theta); s=sin(theta)
c*
c*
c*       [  T11=+c^2     T12=+s^2    0      T14=+2*s*c       0          0      ]
c*       [  T21=+s^2     T22=+c^2    0      T24=-2*s*c       0          0      ]
c*  [T]= [   0            0          1       0               0          0      ]
c*       [  T41=-s*c     T42=+s*c    0      T44=c^2-s^2      0          0      ]
c*       [   0            0          0       0              T55=+c     T56=-s  ] 
c*       [   0            0          0       0              T65=+s     T66=+c  ]
c*
c
      subroutine inverse_matrix(theta,T,TI)
	implicit real*8 (a-h,o-z)
	real*8 theta,T(6,6),TI(6,6),theta_radian
	real*8 c,s,pi,a,b
	pi=3.1415926535897932
	theta_radian=theta*pi/180
	c=dcos(theta_radian)
      s=dsin(theta_radian)
	if (theta.eq.90) then
		s = 1
          c = 0
	else if (theta.eq.0) then
		s = 0
          c = 1
	end if
	do i=1,6
	 do j=1,6
	  T(i,j)=0
	  TI(i,j)=0
	 end do
      end do
	T(1,1)=c**2
	T(1,2)=s**2 
	T(1,4)=2*s*c
	T(2,1)=s**2
	T(2,2)=c**2 
	T(2,4)=-2*s*c
	T(3,3)=1
	T(4,1)=-1*s*c
	T(4,2)=s*c
	T(4,4)=c**2-s**2
	T(5,5)=c
	T(5,6)=-1*s
	T(6,5)=s
	T(6,6)=c

	a=T(1,1)*T(2,2)*T(4,4)-T(1,1)*T(2,4)*T(4,2)-T(1,2)*T(2,1)*T(4,4)  
     &  +T(1,2)*T(2,4)*T(4,1)+T(1,4)*T(2,1)*T(4,2)-T(1,4)*T(2,2)*T(4,1)
	b=T(5,5)*T(6,6)-T(5,6)*T(6,5)
	TI(1,1)=(T(2,2)*T(4,4)-T(2,4)*T(4,2))/a	
	TI(1,2)=-(T(1,2)*T(4,4)-T(1,4)*T(4,2))/a	
	TI(1,4)=(T(1,2)*T(2,4)-T(1,4)*T(2,2))/a	
	TI(2,1)=-(T(2,1)*T(4,4)-T(2,4)*T(4,1))/a	
	TI(2,2)=(T(1,1)*T(4,4)-T(1,4)*T(4,1))/a	
	TI(2,4)=-(T(1,1)*T(2,4)-T(1,4)*T(2,1))/a
	TI(3,3)=1
	TI(4,1)=(T(2,1)*T(4,2)-T(2,2)*T(4,1))/a
	TI(4,2)=-(T(1,1)*T(4,2)-T(1,2)*T(4,1))/a
	TI(4,4)=(T(1,1)*T(2,2)-T(1,2)*T(2,1))/a
	TI(5,5)=T(6,6)/b
	TI(5,6)=-1*T(5,6)/b
	TI(6,5)=-1*T(6,5)/b
	TI(6,6)=T(5,5)/b 
	    
      return
	end
c*************************************************************************** 
c  the following subroutine calculates the Global stiffness matrix [CG_S] 
c  Sublaminate model notation
c
c         [  CG11      CG12      CG14      CG13       0         0   ]
c         [  CG21      CG22      CG24      CG23       0         0   ]
c [CG_S]= [  CG41      CG42      CG44      CG43       0         0   ]
c         [  CG31      CG32      CG34      CG33       0         0   ]
c         [   0         0         0         0        CG66      CG65 ]
c         [   0         0         0         0        CG56      CG55 ]
c
      subroutine get_CG_S(CG_A,CG_S)
      implicit real*8 (a-h,o-z)
	real*8 CG_A(6,6),CG_S(6,6)
      do i=1,6
	 do j=1,6
	  CG_S(i,j)=0
	 end do
      end do
	CG_S(1,1)=CG_A(1,1)
	CG_S(1,2)=CG_A(1,2)
	CG_S(1,3)=CG_A(1,4)
	CG_S(1,4)=CG_A(1,3)
      CG_S(2,1)=CG_A(2,1)
	CG_S(2,2)=CG_A(2,2)
	CG_S(2,3)=CG_A(2,4) 
	CG_S(2,4)=CG_A(2,3)
	CG_S(3,1)=CG_A(4,1)
	CG_S(3,2)=CG_A(4,2)
	CG_S(3,3)=CG_A(4,4)
	CG_S(3,4)=CG_A(4,3)
	CG_S(4,1)=CG_A(3,1)
	CG_S(4,2)=CG_A(3,2)
	CG_S(4,3)=CG_A(3,4)
	CG_S(4,4)=CG_A(3,3)
	CG_S(5,5)=CG_A(6,6)
	CG_S(5,6)=CG_A(6,5)
	CG_S(6,5)=CG_A(5,6)
	CG_S(6,6)=CG_A(5,5)
      return
	end

c*******************************************************************************
c  the following subroutine calculates the GLOBAL [G_ABED] matrix
c  for an orthotropic layer according to the Sublaminate model notation [CG_S]
c
c  AG=CG_II-CG_IO*inv(CG_OO)*CG_OI
c  BG=CG_IO*inv(CG_OO)
c  EG=-1*inv(CG_OO)*CG_OI
c  DG=inv(CG_OO)
c
c  where:
c
c         [  CG11     CG12     CG13 ]           [  CG14     0     0   ]   
c  CG_II= [  CG21     CG22     CG23 ] ;  CL_IO= [  CG24     0     0   ] 
c         [  CG41     CG42     CG44 ]           [  CG43     0     0   ] 
c
c         [  CG31    CG32    CG34]           [  CG33     0     0   ]   
c  CG_OI= [  0       0       0   ] ;  CL_OO= [  0       CG66  CG65 ] 
c         [  0       0       0   ]           [  0       CG56  CG55 ]
c 
c           [  AG11       AG12      AG13  |    BG11        0         0   ]
c           [  AG21       AG22      AG23  |    BG21        0         0   ]
c [G_ABED]= [  AG31       AG32      AG33  |    BG31        0         0   ]
c           [ -----------------------------------------------------------]
c           [  EG11       EG12      EG13  |    DG11        0         0   ]
c           [   0          0          0   |    0         DG22       DG23 ]
c           [   0          0          0   |    0         DG32       DG33 ]
c
      subroutine get_G_ABED(CG_S,G_ABED)
      implicit real*8 (a-h,o-z)
	real*8 CG_S(6,6),s(6,6),G_ABED(6,6)
	real*8 AG(3,3),BG(3,3),EG(3,3),DG(3,3)
	do i=1,3
	 do j=1,3
	  AG(i,j)=0
        BG(i,j)=0
	  EG(i,j)=0
	  DG(i,j)=0
	 end do
      end do
	do i=1,6
	 do j=1,6
	  s(i,j)=CG_S(i,j)
	  G_ABED(i,j)=0
       end do
      end do

	
c	 AG
	AG(1,1)=s(1,1) - (s(1,4)*s(4,1))/s(4,4)
	AG(1,2)=s(1,2) - (s(1,4)*s(4,2))/s(4,4)
	AG(1,3)=s(1,3) - (s(1,4)*s(4,3))/s(4,4)
	AG(2,1)=s(2,1) - (s(2,4)*s(4,1))/s(4,4)
	AG(2,2)=s(2,2) - (s(2,4)*s(4,2))/s(4,4)
	AG(2,3)=s(2,3) - (s(2,4)*s(4,3))/s(4,4)
	AG(3,1)=s(3,1) - (s(3,4)*s(4,1))/s(4,4)
	AG(3,2)=s(3,2) - (s(3,4)*s(4,2))/s(4,4)
	AG(3,3)=s(3,3) - (s(3,4)*s(4,3))/s(4,4)
c
c      BG
      BG(1,1)=s(1,4)/s(4,4)
	BG(2,1)=s(2,4)/s(4,4)
	BG(3,1)=s(3,4)/s(4,4) 
c
c      EG
      EG(1,1)=-s(1,4)/s(4,4)
	EG(1,2)=-s(2,4)/s(4,4)
	EG(1,3)=-s(3,4)/s(4,4)
c
c      DG
      DG(1,1)=1/s(4,4)
	DG(2,2)=s(6,6)/(s(5,5)*s(6,6) - s(5,6)*s(6,5)) 
	DG(2,3)=-s(5,6)/(s(5,5)*s(6,6) - s(5,6)*s(6,5))
	DG(3,2)=-s(6,5)/(s(5,5)*s(6,6) - s(5,6)*s(6,5))
	DG(3,3)=s(5,5)/(s(5,5)*s(6,6) - s(5,6)*s(6,5)) 
c
c      G_ABED
      G_ABED(1,1)=AG(1,1)
	G_ABED(1,2)=AG(1,2)
	G_ABED(1,3)=AG(1,3)
	G_ABED(2,1)=AG(2,1)
	G_ABED(2,2)=AG(2,2)
	G_ABED(2,3)=AG(2,3)
	G_ABED(3,1)=AG(3,1)
	G_ABED(3,2)=AG(3,2)
	G_ABED(3,3)=AG(3,3)
	G_ABED(1,4)=BG(1,1)
	G_ABED(2,4)=BG(2,1)
	G_ABED(3,4)=BG(3,1)
	G_ABED(4,1)=EG(1,1)
	G_ABED(4,2)=EG(1,2)
	G_ABED(4,3)=EG(1,3)
	G_ABED(4,4)=DG(1,1)
	G_ABED(5,5)=DG(2,2)
	G_ABED(5,6)=DG(2,3)
	G_ABED(6,5)=DG(3,2)
	G_ABED(6,6)=DG(3,3)
	return
	end
c*************************************************************************
c the following subroutine calculates the EFFECTIVE tangent stiffness matrix
c [eff_CG_S] of the N+1 orthotropic layers of the sublaminate model according
c to Sublaminate model notation [CG_S]
c
c  eff_CG_OO=inv(eff_DG)
c  eff_CG_IO=inv(eff_DG)*eff_BG
c  eff_CG_OI=-1*inv(eff_DG)*eff_EG
c  eff_CG_II=eff_AG+eff_CG_IO*inv(eff_CG_OO)*eff_CG_OI

      subroutine get_eff_CG_S(eff_G_ABED,eff_CG_S)
      implicit real*8 (a-h,o-z)
	real*8 eff_CG_S(6,6),eff_G_ABED(6,6),f(6,6),A_bar(3,3),B_bar(3,3),
     & D_bar(3,3),D_bar_inv(3,3),C_ii_bar(3,3),C_io_bar(3,3),
     & C_oi_bar(3,3),C_oo_bar(3,3),C_bar(6,6)
	do i=1,6
	 do j=1,6
	  eff_CG_S(i,j)=0
        f(i,j)=eff_G_ABED(i,j)
       end do
      end do

      eff_CG_S(1,1)=f(1,4)**2/f(4,4) + f(1,1)

	eff_CG_S(1,2)=f(1,2) + (f(1,4)*f(2,4))/f(4,4)

	eff_CG_S(1,3)=f(1,3) + (f(1,4)*f(3,4))/f(4,4)

	eff_CG_S(1,4)=f(1,4)/f(4,4)

      eff_CG_S(2,1)=f(2,1) + (f(1,4)*f(2,4))/f(4,4)

	eff_CG_S(2,2)=f(2,4)**2/f(4,4) + f(2,2)

	eff_CG_S(2,3)=f(2,3) + (f(2,4)*f(3,4))/f(4,4)
 
	eff_CG_S(2,4)=f(2,4)/f(4,4)

	eff_CG_S(3,1)=f(3,1) + (f(1,4)*f(3,4))/f(4,4)

	eff_CG_S(3,2)=f(3,2) + (f(2,4)*f(3,4))/f(4,4)

	eff_CG_S(3,3)=f(3,4)**2/f(4,4) + f(3,3)

	eff_CG_S(3,4)=f(3,4)/f(4,4)

	eff_CG_S(4,1)=f(1,4)/f(4,4)

	eff_CG_S(4,2)=f(2,4)/f(4,4)
	eff_CG_S(4,3)=f(3,4)/f(4,4)

      eff_CG_S(4,4)=1/f(4,4)
    

	eff_CG_S(5,5)=f(6,6)/(f(5,5)*f(6,6) - f(5,6)*f(6,5))

	eff_CG_S(5,6)=-f(5,6)/(f(5,5)*f(6,6) - f(5,6)*f(6,5))

	eff_CG_S(6,5)=-f(6,5)/(f(5,5)*f(6,6) - f(5,6)*f(6,5))

	eff_CG_S(6,6)=f(5,5)/(f(5,5)*f(6,6) - f(5,6)*f(6,5))
	
	return
	end
c*************************************************************************
c the following subroutine calculates the EFFECTIVE tangent stiffness matrix [eff_CG_A]
c of the orthotropic layers of the sublaminate model according to 
c Abaqus VUMAT notation [CG_A]
c
      subroutine get_eff_CG_A(eff_CG_S,eff_CG_A)
      implicit real*8 (a-h,o-z)
	real*8 eff_CG_S(6,6),eff_CG_A(6,6)
      do i=1,6
	 do j=1,6
	  eff_CG_A(i,j)=0
       end do
      end do
c
      eff_CG_A(1,1)=eff_CG_S(1,1)
	eff_CG_A(1,2)=eff_CG_S(1,2)
	eff_CG_A(1,3)=eff_CG_S(1,4)
	eff_CG_A(1,4)=eff_CG_S(1,3)
      eff_CG_A(2,1)=eff_CG_S(2,1)
	eff_CG_A(2,2)=eff_CG_S(2,2)
	eff_CG_A(2,3)=eff_CG_S(2,4) 
	eff_CG_A(2,4)=eff_CG_S(2,3)
	eff_CG_A(3,1)=eff_CG_S(4,1)
	eff_CG_A(3,2)=eff_CG_S(4,2)
	eff_CG_A(3,3)=eff_CG_S(4,4)
	eff_CG_A(3,4)=eff_CG_S(4,3)
	eff_CG_A(4,1)=eff_CG_S(3,1)
	eff_CG_A(4,2)=eff_CG_S(3,2)
	eff_CG_A(4,3)=eff_CG_S(3,4)
	eff_CG_A(4,4)=eff_CG_S(3,3)
	eff_CG_A(5,5)=eff_CG_S(6,6)
	eff_CG_A(5,6)=eff_CG_S(6,5)
	eff_CG_A(6,5)=eff_CG_S(5,6)
	eff_CG_A(6,6)=eff_CG_S(5,5)
      return
	end
c*************************************************************************
c the following subroutine calculates the stress vector SG(6) at 
C the material point according to Abaqus VUMAT notation
c
      subroutine get_G_stress_A(eff_CG_A,EG,SG)
      implicit real*8 (a-h,o-z)
	real*8 eff_CG_A(6,6),SG(6),EG(6)
      SG(1)=eff_CG_A(1,1)*EG(1)+eff_CG_A(1,2)*EG(2)+eff_CG_A(1,3)*EG(3)+
     & eff_CG_A(1,4)*EG(4)

      SG(2)=eff_CG_A(2,1)*EG(1)+eff_CG_A(2,2)*EG(2)+eff_CG_A(2,3)*EG(3)+
     & eff_CG_A(2,4)*EG(4)

	SG(3)=eff_CG_A(3,1)*EG(1)+eff_CG_A(3,2)*EG(2)+eff_CG_A(3,3)*EG(3)+
     & eff_CG_A(3,4)*EG(4)

	SG(4)=eff_CG_A(4,1)*EG(1)+eff_CG_A(4,2)*EG(2)+eff_CG_A(4,3)*EG(3)+
     & eff_CG_A(4,4)*EG(4)

	SG(5)=eff_CG_A(5,5)*EG(5)+eff_CG_A(5,6)*EG(6)

      SG(6)=eff_CG_A(6,5)*EG(5)+eff_CG_A(6,6)*EG(6)
      return
	end
c**************************************************************************
c the following subroutine calculates the stress & strain vectors EG, SG for an
c orthotropic layer of the sublaminate model according to 
c Sublaminate model notation and Global coordinate system
c
c    {G_SiEo(k)}=[G_ABED(k)]*{eff_G_EiSo}
c where: {eff_G_EiSo}={EG(1),EG(2),EG(4),SG(3),SG(6),SG(5)}
c
      subroutine get_G_E_S_vectors(G_ABED,eff_G_EiSo,EG,SG)
      implicit real*8 (a-h,o-z)
	real*8 G_ABED(6,6),eff_G_EiSo(6),G_SiEo(6),EG(6),SG(6)
	
		 !G_SiEo=matmul(G_ABED,eff_G_EiSo) - (8.42)
	G_SiEo(1)=G_ABED(1,1)*eff_G_EiSo(1)+G_ABED(1,2)*eff_G_EiSo(2)+
     & G_ABED(1,3)*eff_G_EiSo(3)+G_ABED(1,4)*eff_G_EiSo(4)

	G_SiEo(2)=G_ABED(2,1)*eff_G_EiSo(1)+G_ABED(2,2)*eff_G_EiSo(2)+
     & G_ABED(2,3)*eff_G_EiSo(3)+G_ABED(2,4)*eff_G_EiSo(4)

      G_SiEo(3)=G_ABED(3,1)*eff_G_EiSo(1)+G_ABED(3,2)*eff_G_EiSo(2)+
     & G_ABED(3,3)*eff_G_EiSo(3)+G_ABED(3,4)*eff_G_EiSo(4)
      
	G_SiEo(4)=G_ABED(4,1)*eff_G_EiSo(1)+G_ABED(4,2)*eff_G_EiSo(2)+
     & G_ABED(4,3)*eff_G_EiSo(3)+G_ABED(4,4)*eff_G_EiSo(4)

	G_SiEo(5)=G_ABED(5,5)*eff_G_EiSo(5)+G_ABED(5,6)*eff_G_EiSo(6)

	G_SiEo(6)=G_ABED(6,5)*eff_G_EiSo(5)+G_ABED(6,6)*eff_G_EiSo(6)
c EG
      EG(1)=eff_G_EiSo(1)
      EG(2)=eff_G_EiSo(2)
      EG(3)=eff_G_EiSo(3)
	EG(4)=G_SiEo(4)
      EG(5)=G_SiEo(5)
	EG(6)=G_SiEo(6)
c SG
      SG(1)=G_SiEo(1)
      SG(2)=G_SiEo(2)
      SG(3)=G_SiEo(3)
	SG(4)=eff_G_EiSo(4)
      SG(5)=eff_G_EiSo(5)
	SG(6)=eff_G_EiSo(6)
	
	return
	end
c*****************************************************************************
c the following subroutine calculates the stress & strain vectors EL, SL for an
c orthotropic layer of the sublaminate model according to 
c Sublaminate model notation and Local coordinate system     
 
      subroutine get_L_E_S_vectors(theta,EG,SG,EL,SL)
      implicit real*8 (a-h,o-z)
	real*8 theta,EG(6),SG(6),EL(6),SL(6),T(6,6)
      real*8 c,s,pi
	pi=3.1415926535897932
	theta_radian=theta*pi/180
	c=dcos(theta_radian)
      s=dsin(theta_radian)    
      do i=1,6
	 do j=1,6
	  T(i,j)=0
	 end do
      end do
	  ! This is the Rotation Matrix of E according to the sublaminate notation
        ! (NOTE THAT VUMAT USE epsilon shear strains) 
	T(1,1)=c**2
	T(1,2)=s**2 
	T(1,3)=2*s*c
	T(2,1)=s**2
	T(2,2)=c**2 
	T(2,3)=-2*s*c
	T(3,1)=-2*s*c
	T(3,2)=2*s*c
	T(3,3)=c**2-s**2
	T(4,4)=1
	T(5,5)=2*c
	T(5,6)=2*s
	T(6,5)=-2*s
	T(6,6)=2*c

      EL(1)=EG(1)*T(1,1)+EG(2)*T(1,2)+EG(3)*T(1,3)
      EL(2)=EG(1)*T(2,1)+EG(2)*T(2,2)+EG(3)*T(2,3)
      EL(3)=EG(1)*T(3,1)+EG(2)*T(3,2)+EG(3)*T(3,3)
	EL(4)=EG(4)
      EL(5)=EG(5)*T(5,5)+EG(6)*T(5,6)
	EL(6)=EG(5)*T(6,5)+EG(6)*T(6,6)

      SL(1)=SG(1)*T(1,1)+SG(2)*T(1,2)+SG(3)*T(1,3)
      SL(2)=SG(1)*T(2,1)+SG(2)*T(2,2)+SG(3)*T(2,3)
      SL(3)=SG(1)*T(3,1)+SG(2)*T(3,2)+SG(3)*T(3,3)
	SL(4)=SG(4)
      SL(5)=SG(5)*0.5*T(5,5)+SG(6)*0.5*T(5,6)
	SL(6)=SG(5)*0.5*T(6,5)+SG(6)*0.5*T(6,6)
	return
	end
! ***************************************************************************
! The followng subroutine calculate the local abaqus notation C matrix. given the initial props and the five damage parameters.
! It is required that the D's and E's will take into consideration the mode (tension or compression) before calling this subroutine.
      subroutine Get_local_aba_C_matrix_prog_damage(E11,E22,E33,G12,G23,
     & G13,D11,D22,D12,CL_A_damaged)
! 	a layer reaching this subroutine has initialized damage progression. Hence, only the diagonal parameters remain.
	real*8 E11,E22,E33,G12,G23,G13,D11,D22,D12
	real*8 CL_A_damaged(6,6)
	do i=1,6
	   do j=1,6
	      CL_A_damaged(i,j)=0
	   end do
	end do
C         
	CL_A_damaged(1,1)=E11*(1-D11)
	CL_A_damaged(2,2)=E22*(1-D22)
	CL_A_damaged(3,3)=E33
	CL_A_damaged(4,4)=2*G12*(1-D12)
	CL_A_damaged(5,5)=2*G23
	CL_A_damaged(6,6)=2*G13
	return
      end
C***************************************************************************
C The following subroutine transforms the global IP strain vector to the local one of the layer.
C This is done for determination of tension or compression modes in fiber and matrix. Sublaminate notation
! This is a 2D version of the Get_L_E_S_vectors

      subroutine Get_local_IP_strain(theta,Global_E,local_E)
      real*8 theta,Global_E(3),local_E(3),costheta,sintheta,pi
      real*8 theta_radian
      pi=3.1415926535897932
      theta_radian=  theta*pi/180
      costheta=dcos(theta_radian)
      sintheta=dsin(theta_radian)
      if (theta.eq.90) then
		sintheta = 1
          costheta = 0
	else if (theta.eq.0) then
		sintheta = 0
          costheta = 1
	end if
C
      local_E(1) = global_E(1)*costheta**2 + global_E(2)*sintheta**2 
     &					+ 2*global_E(3)*sintheta*costheta
      local_E(2) = global_E(1)*sintheta**2 + global_E(2)*costheta**2
     &					- 2*global_E(3)*sintheta*costheta
      local_E(3) = (global_E(2)-global_E(1))*sintheta*costheta 
     &					+ global_E(3)*(costheta**2-sintheta**2)
		! This is the shear strain in the SL notation
      return
      end

C***************************************************************************
C This Subroutine creates the interface layer C matrix
      subroutine get_interface_layer_C(E33,damage,CL_A_interface)
      real*8 E33,damage,CL_A_interface(6,6)
      real*8 adj_dmg
C
      do k=1,6
          do j=1,6
              CL_A_interface(k,j)=0
          end do
      end do
C
      CL_A_interface(1,1) = E33*1.0d-10
      CL_A_interface(2,2) = E33*1.0d-10
      CL_A_interface(4,4) = E33*1.0d-10
C
      if (damage.eq.0) then     
          CL_A_interface(3,3) = E33
          CL_A_interface(5,5) = E33
          CL_A_interface(6,6) = E33
      else
		adj_dmg = (-11.555555-20.444444*damage)

          CL_A_interface(3,3) = E33*DEXP(adj_dmg)
          CL_A_interface(5,5) = E33!*DEXP(adj_dmg)
          CL_A_interface(6,6) = E33!*DEXP(adj_dmg)
		  
! Important notice: since the OOP stiffness is calculated as a sum of the tk/E33k, the delamination damage will effect the material only when it's close to 1.
! given ti =tk*1E-10. significant effect, the effective damage needs to be at least 1-1E-10 for a 50% reduction etc
      end if
C
      return
      end
C***************************************************************************
C This Subroutine calculates the new damage parameter given beta gamma strain_I and current strain 
C E_I is strain at initialization of damage, D is the damage fraction, beta is many time maximal strain is larger than E_I, gamma is maximal allowed damage
      subroutine get_updated_D(E_I,E_atm,D_old,beta,gamma,D_new)
      real*8 E_I,E_atm,D_old,beta,gamma,D_new
C
      D_new = DMAX1(beta*(E_atm-E_I)/(E_atm*(beta-1)),D_old)
      if (D_new.gt.gamma) then 
          D_new = gamma
      end if
C    
      return
      end
C*************************************************************************** 
C This subroutine updates initial 5 E's and 5 D's of a layer. according to the TW polynomial
C
      subroutine Update_layer_E_I(local_E_SL,I_max_TW,
     &   E_I1t_old,E_I1c_old,E_I2t_old,E_I2c_old,E_I12_old,
     &   E_I1t_new,E_I1c_new,E_I2t_new,E_I2c_new,E_I12_new)
C
      real*8 local_E_SL(6),I_max_TW
      real*8 E_I1t_old,E_I1c_old,E_I2t_old,E_I2c_old,E_I12_old
      real*8 E_I1t_new,E_I1c_new,E_I2t_new,E_I2c_new,E_I12_new
C         
          ! mode 1t initiation
          if ((I_max_TW.eq.1).and.
     &        (local_E_SL(1).gt.0).and.
     &        (E_I1t_Old.eq.0)) then
              E_I1t_new = local_E_SL(1)
              E_I1c_new = E_I1c_Old
              E_I2t_new = E_I2t_Old
              E_I2c_new = E_I2c_Old
              E_I12_new = E_I12_Old
C
          ! mode 1c initiation
          else if ((I_max_TW.eq.1).and.
     &        (local_E_SL(1).lt.0).and.
     &        (E_I1c_Old.eq.0)) then
              E_I1t_new = E_I1t_Old
              E_I1c_new = local_E_SL(1)
              E_I2t_new = E_I2t_Old
              E_I2c_new = E_I2c_Old
              E_I12_new = E_I12_Old
C
          ! mode 2t initiation
          else if ((I_max_TW.eq.2).and.
     &        (local_E_SL(2).gt.0).and.
     &        (E_I2t_Old.eq.0)) then
              E_I1t_new = E_I1t_Old
              E_I1c_new = E_I1c_Old
              E_I2t_new = local_E_SL(2)
              E_I2c_new = E_I2c_Old
              E_I12_new = E_I12_Old
C
          ! mode 2c initiation
          else if ((I_max_TW.eq.2).and.
     &        (local_E_SL(2).lt.0).and.
     &        (E_I2c_Old.eq.0)) then
              E_I1t_new = E_I1t_Old
              E_I1c_new = E_I1c_Old
              E_I2t_new = E_I2t_Old
              E_I2c_new = local_E_SL(2)
              E_I12_new = E_I12_Old
C
          ! mode 12 initiation
          else if ((I_max_TW.eq.3).and.
     &        (E_I12_Old.eq.0)) then
              E_I1t_new = E_I1t_Old
              E_I1c_new = E_I1c_Old
              E_I2t_new = E_I2t_Old
              E_I2c_new = E_I2c_Old
              E_I12_new = ABS(local_E_SL(3))
          else ! TW is gt 1, but no damage initiation occurred
			E_I1t_new = E_I1t_Old
              E_I1c_new = E_I1c_Old
              E_I2t_new = E_I2t_Old
              E_I2c_new = E_I2c_Old
              E_I12_new = E_I12_Old
          end if
      return
      end
C ************************************************************