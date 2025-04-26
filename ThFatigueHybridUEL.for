      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),PROPS(NPROPS),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
      
      PARAMETER (ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,THREE=3.D0,
     1 TOLER=1.0D-8,FOUR=4.D0,RP25=0.25D0,HALF=0.5D0,SIX=6.D0,
     2 N_ELEM=310707,AlphaEle=1821,NSTVTO=6,NSTVTT=10,NSTV=18)  

      INTEGER I,J,L,K,K1,K2,K3,K4,ODD,EVEN


      REAL*8 WT(NNODE),LC(4,2),XI(2),dNdxi(NNODE,2),
     1 VJACOB(2,2),dNdx(NNODE,2),N(NNODE,1),BPh(2,NDOFEL),DP(2),
     2 SDV(NSTV),BB(3,NDOFEL),CMAT(3,3),STRAIN(3),STRESS(3),
     3 ELSTRAIN(3) 
      
      REAL*8 DETJACOBIAN,INVJACOBIAN(2,2),EMOD,ENU,THCK,kVAL,LPAR,GC,
     1 HIST,ENG,EIVL(2),TensileStr(2),CompressiveStr(2),ALPHAI(2),
     2 alphaTEMP,DELTEMP,alph,alphn,alphBn,STRAINTEMP,TMIN
      
      COMMON/KUSER/USRVAR(N_ELEM,NSTV,4),iter

	   
      IF (JTYPE.EQ.ONE) THEN
           PreviIncm = USRVAR(JELEM,18,1)
        IF (KINC.NE.PreviIncm) THEN
C	    write(6,*) KINC,PreviIncm
            USRVAR(JELEM,18,1) = KINC
	    iter = ZERO
	    USRVAR(JELEM,17,1) = zero
	    else
	    USRVAR(JELEM,17,1) = USRVAR(JELEM,17,1) + one
	  ENDIF
	    iter = USRVAR(JELEM,17,1)
	  
C        write(6,*) KINC,JELEM,iter
C     ==================================================================
C     Material parameters
C     ==================================================================
       LPAR = PROPS(1)
       GC = PROPS(2)
       THCK = PROPS(3)
       kVAL = PROPS(4)
       kFlagF=props(5)
        
       if (kFlagF.eq.1) then
        alphT=GC/(2.d0*6.d0*LPAR)
      else
        alphT=1.d10
      endif
      
	  
C     Initialize the Matrices********************************************         
      DO K1 = 1, NDOFEL                      
        DO KRHS = 1, NRHS
         RHS(K1,KRHS) = ZERO
        END DO
        DO K2 = 1, NDOFEL
         AMATRX(K2,K1) = ZERO

        END DO
      END DO 
      
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================      
      
       LC(1,1) = -ONE/THREE**HALF
       LC(1,2) = -ONE/THREE**HALF
       LC(2,1) = ONE/THREE**HALF
       LC(2,2) = -ONE/THREE**HALF
       LC(3,1) = ONE/THREE**HALF
       LC(3,2) = ONE/THREE**HALF
       LC(4,1) = -ONE/THREE**HALF
       LC(4,2) = ONE/THREE**HALF
       DO I=1,4
        WT(I) = ONE
       END DO  
C      *****************************************************************       
C      Calculate at the integration points 
C      *****************************************************************
       DO INPT=1,NNODE
    
C     Local coordinates of the integration point
        XI(1) = LC(INPT,1)
        XI(2) = LC(INPT,2)      
C     Call the function to obtain the shape functions and its derivatives        
      call SHAPEFUN(N,dNdxi,XI) 
      
C     Lets find the jacobian
     
      DO I = 1,2
         DO J = 1,2
          VJACOB(I,J) = ZERO
          DO K = 1,NNODE
           VJACOB(I,J) = VJACOB(I,J) + COORDS(J,K)*dNdxi(K,I)
          END DO
         END DO
       END DO
      
      
      DETJACOBIAN = ZERO
      
      DETJACOBIAN = VJACOB(1,1)*VJACOB(2,2) - VJACOB(2,1)*VJACOB(1,2) 
      IF (DETJACOBIAN.LT.ZERO) THEN
         WRITE(7,*) 'Negative Jacobian',DETJACOBIAN
         CALL XIT	
        END IF
C     *****************************************************************
C     Inverse of Jacobian
C     ***************************************************************** 
           
      INVJACOBIAN(1,1) = VJACOB(2,2)/DETJACOBIAN
      INVJACOBIAN(1,2) = -VJACOB(1,2)/DETJACOBIAN
      INVJACOBIAN(2,1) = -VJACOB(2,1)/DETJACOBIAN
      INVJACOBIAN(2,2) = VJACOB(1,1)/DETJACOBIAN

C     *****************************************************************      
C     Calculate the B Matrix
C     *****************************************************************
      DO K = 1,NNODE
         DO I = 1,2
          dNdx(K,I) = ZERO
          DO J = 1,2
           dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*INVJACOBIAN(I,J)
          END DO
         END DO
      END DO      
      
C    ***********************************************************
    
      DO I=1,NNODE
          BPh(1,I)  = dNdx(I,1)
          BPh(2,I) = dNdx(I,2)
      END DO
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
        PHASE=ZERO
        DPHASE=ZERO
      DO I=1,4
         PHASE=PHASE+N(I,1)*U(I)
      END DO

C	  IF (PHASE.GT.1.d0) THEN
C	     PHASE=1.d0
C	  ENDIF  

       SDV(1) = PHASE
C     **************************************************************      
C     Gradient of phi
C     **************************************************************      
        DO I=1,2
         DP(I)=ZERO
        END DO
        DO I=1,2
         DO J=1,4
          DP(I)=DP(I)+BPh(I,J)*U(J)
         END DO
      END DO
		
C     ******************************************************************      
C     Calculate elastic energy      
C     ******************************************************************

        IF (iter.EQ.ZERO) THEN

         ENGN = USRVAR(JELEM,8,INPT)
         
        ELSE

         ENGN = USRVAR(JELEM,2,INPT)
         
        ENDIF
        
        HISTN = USRVAR(JELEM,2,INPT)
        IF (ENGN.GT.HISTN) THEN
         HIST=ENGN
        ELSE
         HIST = HISTN
        ENDIF
        SDV(2) = HIST

!     =================================================================
!     Fatigue Fracture
!     ================================================================= 
       alph = USRVAR(JELEM,8,INPT)*((1-PHASE)**TWO+kVAL) 
       IF (iter.EQ.ZERO) THEN
          alphn = USRVAR(JELEM,3,INPT)
		  alphBn = USRVAR(JELEM,5,INPT)
		  SDV(4) = alphn
		  SDV(6) = alphBn
      else
          alphn = USRVAR(JELEM,4,INPT)
		  alphBn = USRVAR(JELEM,6,INPT)
		  SDV(4) = alphn
		  SDV(6) = alphBn
      endif   
      
C       alphBn = USRVAR(JELEM,5,INPT)

!     Update fatigue history variable
       if ((alph.ge.alphn).and.(dtime.gt.0.d0)) then
        alphB = alphBn+abs(alph-alphn)
       else
        alphB=alphBn
       endif      
           
       if (alphB.lt.alphT) then  
        Fdeg= 1.d0
       else
        Fdeg=(2.d0*alphT/(alphB+alphT))**2.d0
      endif  
      
C      if (JELEM.eq.AlphaEle) then
C        cln(2) = alphB
C      endif
      
       SDV(3) = alph
       SDV(5) = alphB
C     ******************************************************************      
C     Calculating element stiffness matrix
C     ==================================================================
        vl=DETJACOBIAN*THCK*WT(INPT)
        AMATRX(1:4,1:4)=AMATRX(1:4,1:4)+vl*(matmul(N,transpose(N))*
     1    (GC/LPAR*Fdeg+TWO*HIST)+(matmul(transpose(BPh),BPh)*
     2   GC*LPAR*Fdeg))
            
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================     
       DO I=1,NDOFEL
         
         RHS(I,1) = RHS(I,1) - (-TWO*(1 - PHASE)*N(I,1)*HIST+Fdeg*GC*
     1      ((ONE/LPAR)*PHASE*N(I,1)))*WT(INPT)*DETJACOBIAN*THCK   
         
         DO J=1,2
          RHS(I,1)=RHS(I,1)- Fdeg*GC*(LPAR*BPh(J,I)*DP(J))*WT(INPT)*
     1      DETJACOBIAN*THCK
         END DO
       END DO

C     ==================================================================
C     Uploading solution dep. variables
C     ==================================================================
        DO I=1,NSTVTO
         SVARS(NSTVTO*(INPT-1)+I) = SDV(I)
         USRVAR(JELEM,I,INPT) = SDV(I)

        END DO
       END DO       
      
      ELSEIF (JTYPE.EQ.TWO) THEN
               
       EMOD = PROPS(1)
       ENU = PROPS(2)
       THCK = PROPS(3)
       kVAL = PROPS(4)  
       alphaTEMP = PROPS(5)	
       DELTEMP = PROPS(6)
!	   For a temperature cycle of (T - To) = 100 K
       TDELTEMP = 100.D0 !T
	   TMIN = ZERO       !To
	   
       ELAMEL=EMOD*ENU/((ONE+ENU)*(ONE-TWO*ENU))
       ELAMEG=EMOD/(TWO*(ONE+ENU))
      
       DO K1 = 1, NDOFEL                      
        DO KRHS = 1, NRHS
         RHS(K1,KRHS) = ZERO
        END DO
        DO K2 = 1, NDOFEL
         AMATRX(K2,K1) = ZERO
        END DO
       END DO

C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
       LC(1,1) = -ONE/THREE**HALF
       LC(1,2) = -ONE/THREE**HALF
       LC(2,1) = ONE/THREE**HALF
       LC(2,2) = -ONE/THREE**HALF
       LC(3,1) = ONE/THREE**HALF
       LC(3,2) = ONE/THREE**HALF
       LC(4,1) = -ONE/THREE**HALF
       LC(4,2) = ONE/THREE**HALF
       DO I=1,FOUR
        WT(I) = ONE
       END DO  
C      *****************************************************************       
C      Calculate at the integration points 
C      *****************************************************************
       DO INPT=1,FOUR
    
C     Local coordinates of the integration point
        XI(1) = LC(INPT,1)
        XI(2) = LC(INPT,2)  
                
C     Call the function to obtain the shape functions and its derivatives        
      call SHAPEFUN(N,dNdxi,XI)      
      
C     Lets find the jacobian
      DO I = 1,2
         DO J = 1,2
          VJACOB(I,J) = ZERO
          DO K = 1,NNODE
           VJACOB(I,J) = VJACOB(I,J) + COORDS(J,K)*dNdxi(K,I)
          END DO
         END DO
       END DO
      
      DETJACOBIAN = ZERO
      
      DETJACOBIAN = VJACOB(1,1)*VJACOB(2,2)- VJACOB(2,1)*VJACOB(1,2) 
	  IF (DETJACOBIAN.LT.ZERO) THEN
         WRITE(7,*) 'Negative Jacobian',DETJACOBIAN
         CALL XIT	
       ENDIF
C     *****************************************************************
C     Inverse of Jacobian
C     *****************************************************************

           
      INVJACOBIAN(1,1) = VJACOB(2,2)/DETJACOBIAN
      INVJACOBIAN(1,2) = -VJACOB(1,2)/DETJACOBIAN
      INVJACOBIAN(2,1) = -VJACOB(2,1)/DETJACOBIAN
      INVJACOBIAN(2,2) = VJACOB(1,1)/DETJACOBIAN
C     *****************************************************************      
C     Calculate the B Matrix      
C     *****************************************************************
      DO K = 1,NNODE
         DO I = 1,2
          dNdx(K,I) = ZERO
          DO J = 1,2
           dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*INVJACOBIAN(I,J)
          END DO
         END DO
      END DO

      
      EVEN=0;
      DO I=1,NNODE
          ODD = EVEN+1
          EVEN = ODD+1
          BB(1,ODD)  = dNdx(I,1)
          BB(1,EVEN) = ZERO
          BB(2,ODD)  = ZERO
          BB(2,EVEN) = dNdx(I,2)
          BB(3,ODD)  = dNdx(I,2)
          BB(3,EVEN) = dNdx(I,1)
          
      END DO
      
C     ==================================================================
C     Calculating materials stiffness matrix (plane strain)
C     ==================================================================
        DO I=1,3
         DO J=1,3
          CMAT(I,J)=ZERO
         END DO
        END DO
C        CMAT(1,1)=EMOD/((ONE-ENU**TWO))
C        CMAT(2,2)=EMOD/((ONE-ENU**TWO))
C        CMAT(3,3)=EMOD/(TWO*(ONE-ENU**TWO))*(ONE-ENU)
C        CMAT(1,2)=EMOD/((ONE-ENU**TWO))*ENU
C        CMAT(2,1)=EMOD/((ONE-ENU**TWO))*ENU

        CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        

		       
C     *******************************************************************
C     Calculating Strain
C     *******************************************************************
        DO J=1,3
         STRAIN(J)=ZERO
        END DO
        DO I=1,3
         DO J=1,NDOFEL
          STRAIN(I)=STRAIN(I)+BB(I,J)*U(J)    
         END DO
        END DO
        
C        SDV(4) = STRAIN(1)
C        SDV(5) = STRAIN(2)
		
C	  ------------------------------------------------------------------	
C     Elastic strain		
C	  ------------------------------------------------------------------
            TTEMPT = USRVAR(JELEM-N_ELEM,14,INPT)
             flag  = USRVAR(JELEM-N_ELEM,16,INPT)
          if (flag.le.zero.and.iter.eq.zero) then
             TTEMPT = TTEMPT + DELTEMP
             diffT = TDELTEMP - TTEMPT
            if (TTEMPT.lt.TDELTEMP.and.diffT.gt.(DELTEMP/two)) then

               STRAINTEMP = TTEMPT*alphaTEMP 
C	       STRAIN(1) = STRAIN(1)-STRAINTEMP
C              STRAIN(2) = STRAIN(2)-STRAINTEMP
C            elseif (TTEMPT.ge.TDELTEMP) then
             else
              flag=one

              STRAINTEMP = TTEMPT*alphaTEMP 
            endif
C            write(6,*) TTEMPT,STRAINTEMP
          elseif (flag.eq.one.and.iter.eq.zero) then 
            
             TTEMPT = TTEMPT - DELTEMP
             diffT = TTEMPT - TMIN
             if (TTEMPT.gt.TMIN.and.diffT.gt.(DELTEMP/two)) then
C             STRAINTEMP = -ABS(TTEMPT-TDELTEMP)*alphaTEMP 
              STRAINTEMP = TTEMPT*alphaTEMP
C	     STRAIN(1) = STRAIN(1)+STRAINTEMP
C            STRAIN(2) = STRAIN(2)+STRAINTEMP
C             elseif (TTEMPT.le.zero) then
             else
              flag=zero
C              STRAINTEMP = -ABS(TTEMPT-TDELTEMP)*alphaTEMP 
               STRAINTEMP = TTEMPT*alphaTEMP
            endif
C            write(6,*) TTEMPT,STRAINTEMP
          else

          STRAINTEMP = USRVAR(JELEM-N_ELEM,15,INPT)  
          	
	  endif
	  
	  
	  STRAIN(1) = STRAIN(1)-STRAINTEMP
      STRAIN(2) = STRAIN(2)-STRAINTEMP	
	  SDV(8) = TTEMPT
	  SDV(9) = STRAINTEMP
	  SDV(10) = flag
C        DO J=1,3
C         SDV(7)=STRAIN(2)
C       END DO
C     *****************************************************************
C     Eigen Decomposition     
C     *****************************************************************
       EIVL(1)=ZERO
       EIVL(2)=ZERO
       EIVL(1)=(STRAIN(1)+STRAIN(2)+SQRT((STRAIN(1)+
     1 STRAIN(2))**TWO-FOUR*(STRAIN(1)*STRAIN(2)-
     2 STRAIN(3)**TWO/FOUR)))/TWO
      
       EIVL(2)=(STRAIN(1)+STRAIN(2)-SQRT((STRAIN(1)+
     1 STRAIN(2))**TWO-FOUR*(STRAIN(1)*STRAIN(2)-
     2 STRAIN(3)**TWO/FOUR)))/TWO
C     ****************************************************************
C      Tensile and compressive strain*********************************

      TensileStr(1) = ZERO
      TensileStr(2) = ZERO
	  CompressiveStr(1) = ZERO
      CompressiveStr(2) = ZERO
      TensileStr(1) = (EIVL(1)+ABS(EIVL(1)))/TWO
      TensileStr(2) = (EIVL(2)+ABS(EIVL(2)))/TWO
      CompressiveStr(1) = (EIVL(1)-ABS(EIVL(1)))/TWO
      CompressiveStr(2) = (EIVL(2)-ABS(EIVL(2)))/TWO
	  EIVLSUM = EIVL(1)+EIVL(2)
	  Trace = STRAIN(1)+STRAIN(2)
	  ALPHA=ZERO
	  BETA=ZERO

       IF (Trace.GT.ZERO) THEN
        ALPHA=ONE
		ELSE
		BETA=ONE
       ENDIF
     
C     *****************************************************************      
C     Calculating stress      
C     *****************************************************************
       DO K1=1,3
         STRESS(k1) = ZERO
        END DO
        DO K1=1,3
         DO K2=1,3
          STRESS(K1) = STRESS(K1)+CMAT(K1,K2)*STRAIN(K2)
         END DO
       END DO
	   
C	   STRESS(1) = STRESS(1) - alphaTEM*DELTEMP*EMOD/(1-TWO*ENU)
C	   STRESS(2) = STRESS(2) - alphaTEM*DELTEMP*EMOD/(1-TWO*ENU)
	   
C        DO J=1,3
C         SDV(J+3) = STRESS(J)
C      END DO
      SDV(4) = STRESS(1)
      SDV(5) = STRESS(2)
      
C     ----------------------------------------------------------------
C     Equivalent stress      
C     ----------------------------------------------------------------
       SDV(6)=SQRT(((STRESS(1)-STRESS(2))**TWO+STRESS(2)**TWO+
     1  STRESS(1)**TWO+SIX*STRESS(3)**TWO)/TWO)
     
C       DO J=1,3
C         SDV(J+6) = STRESS(J)*((ONE-PHASE)**TWO+kVAL)
C      END DO
C     ==================================================================
C     Calculating elastic ENERGY
C     ==================================================================
        ENG=ZERO
        DO K2=1,3
         ENG=ENG+STRESS(K2)*STRAIN(K2)*HALF
        END DO
		
C        SDV(1)=ENG*((ONE-PHASE)**TWO+kVAL)
C        SDV(2)=ENG

	    SDV(2)=((ELAMEL*(ALPHA*Trace)**TWO)/TWO+ELAMEG*
     1   ((TensileStr(1))**TWO+(TensileStr(2))**TWO))

	    SDV(3)=((ELAMEL*(BETA*Trace)**TWO)/TWO+ELAMEG*
     1   ((CompressiveStr(1))**TWO+(CompressiveStr(2))**TWO))

        ENERGY(2)=ENG   
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
        IF (iter.EQ.ZERO) THEN
	      IF (SDV(2).lt.SDV(3)) THEN
		     PHASE=USRVAR(JELEM-N_ELEM,7,INPT)
	         
	      else
		    PHASE = USRVAR(JELEM-N_ELEM,1,INPT)
	      endif
         
        ELSE
		
         PHASE=USRVAR(JELEM-N_ELEM,7,INPT)
		 
        ENDIF

C         PHASE = USRVAR(JELEM-N_ELEM,1,INPT)

        SDV(1)=PHASE 
C     *****************************************************************      
C     Calculating element stiffness matrix      
C     *****************************************************************
      

         AMATRX(1:8,1:8)=AMATRX(1:8,1:8)+
     1   DETJACOBIAN*THCK*WT(INPT)*(((ONE-PHASE)**TWO+kVAL)*
     2   matmul(matmul(transpose(BB),CMAT),BB))
       
C     ******************************************************************
C      Calculate the Residual Vector
C     ******************************************************************
       DO K1=1,NDOFEL
         DO K4=1,3
           RHS(K1,1) = RHS(K1,1) - WT(INPT)*BB(K4,K1)*STRESS(K4)*
     1    DETJACOBIAN*THCK*((ONE-PHASE)**TWO+kVAL)
         END DO
        END DO
C     ******************************************************************
C     Solution dependent variables
C     ******************************************************************
      DO I=1,NSTVTT
        SVARS(NSTVTT*(INPT-1)+I) = SDV(I)
	USRVAR(JELEM-N_ELEM,I+NSTVTO,INPT) = SDV(I)
      END DO      
      END DO
C       elem=(JELEM-N_ELEM)     
C	write(6,*) KINC,iter,elem
      
      ENDIF
      
       RETURN
      END
      
      SUBROUTINE SHAPEFUN(N,dNdxi,xi)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 N(4,1),dNdxi(4,2)
      Real*8 XI(2)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)
      
      
      N(1,1) = ONE/FOUR*(ONE-XI(1))*(ONE-XI(2))
      N(2,1) = ONE/FOUR*(ONE+XI(1))*(ONE-XI(2))
      N(3,1) = ONE/FOUR*(ONE+XI(1))*(ONE+XI(2))
      N(4,1) = ONE/FOUR*(ONE-XI(1))*(ONE+XI(2))
C
C     Derivatives of shape functions respect to local coordinates
      DO I=1,4
        DO J=1,2
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  -ONE/FOUR*(ONE-XI(2))
      dNdxi(1,2) =  -ONE/FOUR*(ONE-XI(1))
      dNdxi(2,1) =  ONE/FOUR*(ONE-XI(2))
      dNdxi(2,2) =  -ONE/FOUR*(ONE+XI(1))
      dNdxi(3,1) =  ONE/FOUR*(ONE+XI(2))
      dNdxi(3,2) =  ONE/FOUR*(ONE+XI(1))
      dNdxi(4,1) =  -ONE/FOUR*(ONE+XI(2))
      dNdxi(4,2) =  ONE/FOUR*(ONE-XI(1))
      RETURN
      END
      
C     UMAT Subroutine      
       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
       CHARACTER*80 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

       PARAMETER (ONE=1.0,TWO=2.0,THREE=3.0,SIX=6.0, HALF=0.5,
     1 N_ELEM=310707,NSTV=18) 
C       DATA NEWTON,TOLER/40,1.D-6/
       
       COMMON/KUSER/USRVAR(N_ELEM,NSTV,4)
       
      
      NELEMAN=NOEL-TWO*N_ELEM
C       DO I=1,1
        STATEV(1) = USRVAR(NELEMAN,7,NPT) !Phase field
C        STATEV(2) = USRVAR(NELEMAN,10,NPT)
C        STATEV(3) = USRVAR(NELEMAN,11,NPT) 
C        STATEV(4) = USRVAR(NELEMAN,12,NPT)  
	    STATEV(2) = USRVAR(NELEMAN,14,NPT)  ! Temperature cycle
C        STATEV(5) = USRVAR(NELEMAN,15,NPT)
C       END DO

      RETURN
      END      
