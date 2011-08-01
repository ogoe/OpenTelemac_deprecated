C                       ************************
                        SUBROUTINE COST_FUNCTION
C                       ************************
C
     *(JCOUT,OPTION,MODE)
C     
C***********************************************************************
C PROGICIEL : TELEMAC 2D        25/06/93  E. BARROS
C             UPGRADE TO 5.2    02/10/00  A. LEOPARDI (UNINA)                                                                       
C***********************************************************************
C
C     FUNCTION : PARTIAL COMPUTATION (ONE TIME STEP) OF COST FUNCTION
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |    JCOUT       |<-->| COST FUNCTION
C |    U,V,H       | -->| VELOCITY AND DEPTH
C |    UD,VD,HD    | -->| MEASURES
C |    NPOIN       | -->| NUMBER OF POINTS IN THE MESH
C |    ALPHA       | -->| WEIGHT FUNCTIONS
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C     APPELE PAR : HOMERE_ADJ_T2D
C     
C     SOUS-PROGRAMME APPELE : RIEN
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
      USE INTERFACE_TELEMAC2D, EX_COST_FUNCTION => COST_FUNCTION
C
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION , INTENT(INOUT) :: JCOUT
      INTEGER , INTENT(IN)             :: OPTION
      CHARACTER(LEN=3) , INTENT(IN)    :: MODE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      COMMON/INFO/LNG,LU
      INTEGER LNG,LU
C
      INTEGER I,J
      DOUBLE PRECISION C
C
C=======================================================================
C
      IF(MODE.EQ.'FCT') THEN
C
C       HERE U,V AND H ARE GIVEN BY A CALL TO LITENR
C
        IF(OPTION.EQ.1) THEN
C
          DO I=1,NPOIN
C
            JCOUT = JCOUT + ALPHA1%R(I) * (H%R(I)-HD%R(I))**2
     *                    + ALPHA2%R(I) * (U%R(I)-UD%R(I))**2
     *                    + ALPHA3%R(I) * (V%R(I)-VD%R(I))**2
C
          ENDDO
C
        ELSEIF(OPTION.EQ.2) THEN
C
          DO I=1,NPOIN
C
            JCOUT = JCOUT + ALPHA1%R(I) * GRAV * 
     *              (SQRT(MAX(H%R(I),0.D0))-SQRT(MAX(HD%R(I),0.D0)))**2
     *                    + ALPHA2%R(I) * (U%R(I)-UD%R(I))**2
     *                    + ALPHA3%R(I) * (V%R(I)-VD%R(I))**2
C
          ENDDO 
C
        ELSE
C
          IF(LNG.EQ.1) THEN 
            WRITE(LU,*) 'COST_FUNCTION : OPTION NON PREVUE : ',MODE
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'COST_FUNCTION: UNEXPECTED OPTION : ',MODE
          ENDIF
          CALL PLANTE(1)
          STOP
C
C       TEST ON OPTION
        ENDIF
C
C=======================================================================
C
      ELSEIF(MODE.EQ.'GRD') THEN
C
      IF(    INCLU2(ESTIME,'FROTTEMENT')
     *   .OR.INCLU2(ESTIME,'FRICTION'  )  ) THEN
C
C     IDENTIFICATION DU FROTTEMENT
C                                                                                                                                                         
      IF (KFROT.EQ.3.OR.KFROT.EQ.2) THEN
C
      CALL SLOPES(TE3,ZF,MESH)
      CALL VECTOR(T1,'=','MASBAS          ',U%ELM,1.D0,
     *            T3,T3,T3,T3,T3,T3,MESH,.TRUE.,TE3)
C                                                                       
      CALL FRICTI(T3,T4,T2,T2,UN,VN,HN,CF,MESH,T2,T5,VERTIC,UNSV2D,MSK,
     *            MASKEL,HFROT)
C                                                                       
      CALL OS( 'X=XY    ' , T3 , T1 , T1 , C )             
      CALL OS( 'X=XY    ' , T4 , T1 , T1 , C )    
C                                                                      
      CALL OS( 'X=CXYZ  ' , T3 , QQ , UU , -2.D0 )               
      CALL OS( 'X=CXYZ  ' , T4 , RR , VV , -2.D0 )  
C                                                                       
      CALL OS( 'X=X+Y   ' , T3 , T4 , T4 , C )
      CALL OS( 'X=Y/Z   ' , T3 , T3 , CHESTR , C )               
C
      ELSE
C
       IF(LNG.EQ.1) THEN
         WRITE(LU,*) 'COST_FUNCTION : FROTTEMENT NON TRAITE : ',KFROT
       ENDIF
       IF(LNG.EQ.2) THEN
         WRITE(LU,*) 'COST_FUNCTION: UNEXPECTED FRICTION LAW: ',KFROT
       ENDIF
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
      ELSE
C
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'COST_FUNCTION : PARAMETRE NON PREVU :'
          WRITE(LU,*) ESTIME
          WRITE(LU,*) 'VERIFIER LE MOT CLEF : ESTIMATION DE PARAMETRES'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'COST_FUNCTION: UNEXPECTED PARAMETER :'
          WRITE(LU,*) ESTIME
          WRITE(LU,*) 'CHECK THE KEY-WORD: PARAMETER ESTIMATION'
        ENDIF
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C COMPUTATION OF A GRADIENT FOR EVERY ZONE, AFTER BUILDING
C A GRADIENT VECTOR FOR EVERY POINT (IN T3)
C
      IF(NZONE.GT.0) THEN
C       SI ON IDENTIFIE UN SEUL PARAMETRE S POUR UN ENSEMBLE DE POINTS
C       LE GRADIENT DJ/DS EST LA SOMME DES GRADIENTS DE CHACUN DES
C       POINTS DE L'ENSEMBLE
        DO J=1,NZONE
          DO I=1,NPOIN
            IF(ZONE%I(I).EQ.J) GRADJ%R(J)=GRADJ%R(J)+T3%R(I)             
          ENDDO
        ENDDO
      ELSE
C       NOTE JMH : ICI ON SUPPOSE QUE NPARAM = NPOIN
        CALL OS('X=X+Y   ',GRADJ,T3,T3,0.D0)
      ENDIF
C
C=======================================================================
C
      ELSEIF(MODE.EQ.'RHS') THEN
C
C           IT    IT    IT
C  TERMS 2 W   ( X   - M   ) OR EQUIVALENT DEPENDING ON THE OPTION
C           IP    IP    IP
C
      IF(OPTION.EQ.1) THEN
C
        CALL OS( 'X=Y-Z   ', CV1 , HH , HD , C )
        CALL OS( 'X=Y-Z   ', CV2 , UU , UD , C )
        CALL OS( 'X=Y-Z   ', CV3 , VV , VD , C )
C
        CALL OS( 'X=CXY   ', CV1 , ALPHA1 , ALPHA1 , 2.D0 )
        CALL OS( 'X=CXY   ', CV2 , ALPHA2 , ALPHA2 , 2.D0 )
        CALL OS( 'X=CXY   ', CV3 , ALPHA3 , ALPHA3 , 2.D0 )
C
      ELSEIF(OPTION.EQ.2) THEN
C
C       HERE COST FUNCTION COMPUTED WITH CELERITY INSTEAD OF DEPTH
        CALL OS( 'X=SQR(Y)', T1  , HH , T1 , C   )
        CALL OS( 'X=SQR(Y)', T2  , HD , T2 , C   )
        CALL OS( 'X=Y-Z   ', T3  , T1 , T2 , C )
        CALL OS( 'X=Y/Z   ', CV1 , T3 , T1 , C )
        CALL OS( 'X=Y-Z   ', CV2 , UU , UD , C )
        CALL OS( 'X=Y-Z   ', CV3 , VV , VD , C )
C
        CALL OS( 'X=CXY   ', CV1 , ALPHA1 , ALPHA1 , GRAV )
        CALL OS( 'X=CXY   ', CV2 , ALPHA2 , ALPHA2 , 2.D0 )
        CALL OS( 'X=CXY   ', CV3 , ALPHA3 , ALPHA3 , 2.D0 )
C
      ENDIF
C
C=======================================================================
C
      ELSE
C
C=======================================================================
C
       IF(LNG.EQ.1) WRITE(LU,*) 'COST_FUNCTION : MODE NON PREVU : ',MODE
       IF(LNG.EQ.2) WRITE(LU,*) 'COST_FUNCTION: UNEXPECTED MODE : ',MODE
       CALL PLANTE(1)
       STOP
C
C=======================================================================
C
C     TEST ON MODE
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
