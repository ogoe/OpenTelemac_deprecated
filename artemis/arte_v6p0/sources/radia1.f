C                       *****************
                        SUBROUTINE RADIA1
C                       *****************
C
     *(LISHHO)
C                                                                       
C***********************************************************************
C
C  ARTEMIS VERSION 5.1   04/06/99    D. AELBRECHT (LNH) 01 30 87 74 12
C                                    F. BECQ      (LNH) 
C
C  LINKED TO BIEF VERS. 5.0          J-M HERVOUET (LNH) 01 30 87 80 18
C    
C***********************************************************************
C
C     FONCTION  : CALCUL DES CONTRIANTES DE RADIATION ET
C     --------    DES FORCES MOTRICES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    OMEGA       | -->|  PULSATION DE LA HOULE                       |
C |    LISHHO      |<-->|  NOMBRE DE LISSAGES POUR LA HAUTEUR DE HOULE |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR            : ARTEMIS
C
C  SOUS-PROGRAMME APPELE : OS, VECTOR, FILTER
C
C
C  REFERENCES : M.W. DINGEMANS, A.C. RADDER and H.J. DE VRIEND
C  ----------   COMPUTATION OF THE DRIVING FORCES OF WACE-INDUCED
C               CURRENTS. Coastal Engineering, 11 (1987) pp 539-563.
C
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_ARTEMIS
      USE INTERFACE_ARTEMIS , EX_RADIA1 => RADIA1
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER I
      INTEGER LISHHO
C
      DOUBLE PRECISION BID 
C
C VARIABLES INTERNES A RADIA1
C
      DOUBLE PRECISION COE , COCO, COSI, SISI , Z(1)
      DOUBLE PRECISION OMEG2
      INTEGER          IRADIA , LISRAD
C
      LOGICAL MAS
C
      INTRINSIC SQRT, ATAN2, DMOD, ABS, COS, SIN
C
C-----------------------------------------------------------------------
C
C     ON FIXE EN DUR LA METHODE DE CALCUL DES CONTRAINTES DE
C     RADIATION APRES LES TEST EFFCTUES : METHODE 2
C
      IRADIA = 2
C
C-----------------------------------------------------------------------
C
C ON DONNE POUR MEMOIRE LA PROGRAMMATION DE LA METHODE 1
C CI-DESSOUS MEME SI ELLE N'EST PAS UTILISEE
C
C=======================================================================
C     CONTRAINTES DE RADIATION........METHODE 1
C=======================================================================
C
      IF(IRADIA.EQ.1) THEN
C
      CALL OS('X=YZ    ' , T3 , PHII , PHII , BID )
      CALL OS('X=X+YZ  ' , T3 , PHIR , PHIR , BID )
      CALL VECTOR(T2 , '=' , 'GRADF          X' , IELM ,
     *            1.D0 , T3 , T1 , T1 , T1 , T1 , T1 ,
     *            MESH , MSK , MASKEL)
      CALL VECTOR(T2 , '=' , 'GRADF          X' , IELM ,
     *            1.D0 , T2 , T1 , T1 , T1 , T1 , T1 ,
     *            MESH , MSK , MASKEL)
      CALL VECTOR(T4 , '=' , 'GRADF          Y' , IELM ,
     *            1.D0 , T3 , T1 , T1 , T1 , T1 , T1 ,
     *            MESH , MSK , MASKEL)
      CALL VECTOR(T4 , '=' , 'GRADF          Y' , IELM ,
     *            1.D0 , T4 , T1 , T1 , T1 , T1 , T1 ,
     *            MESH , MSK , MASKEL)
C
      CALL VECTOR(T1 , '=' , 'MASBAS          ' , IELM ,
     *            1.D0 , T3 , T3 , T3 , T3 , T3 , T3 ,
     *            MESH , MSK , MASKEL)
C
      CALL OS('X=Y+Z   ' , T4 , T2 , T4 , BID   )
      CALL OS('X=Y/Z   ' , T4 , T4 , T1 , BID   )
      CALL OS('X=Y/Z   ' , T4 , T4 , T1 , BID   )
C
      CALL OS('X=CY/Z  ' , T2 , CG , C  , 2.D0  )
      CALL OS('X=X+C   ' , T2 , T3 , T3 , -1.D0 )
      OMEG2 = OMEGA*OMEGA
      CALL OS('X=CX    ' , T2 , T3 , T3 , OMEG2 )
C
      COE = 1.D0/(8.D0*GRAV)
C
C NECESSITE DE RECALCULER U1 ET V1 CAR CES VARIABLES N'EXISTENT 
C PLUS !!
C
C      DO I = 1,NPOIN
C         SXX(I) = COE*(2.D0*C(I)*CG(I)*(U0(I)*U0(I) + U1(I)*U1(I))
C     *             + T2(I)*T3(I) - (GRAV*ZF(I)+ C(I)*CG(I))*T4(I))
C         SXY(I) = COE*(2.D0*C(I)*CG(I)*(V0(I)*U0(I) + V1(I)*U1(I)))
C         SYY(I) = COE*(2.D0*C(I)*CG(I)*(V0(I)*V0(I) + V1(I)*V1(I))
C     *             + T2(I)*T3(I) - (GRAV*ZF(I)+ C(I)*CG(I))*T4(I))
C      END DO
C
C
C=======================================================================
C     CONTRAINTES DE RADIATION........METHODE 2 (IDENTIQUE A TOMAWAC)
C=======================================================================
C
      ELSE
C
      CALL OS('X=Y     ',T3,HHO,SBID,BID)
C
C -------------------------------------------------------------
C LISSAGE DE LA HAUTEUR DE HOULE POUR ELIMINER LES OSCILLATIONS
C PARASITES
C -------------------------------------------------------------
C
      IF(LISHHO.GT.0) THEN
         MAS = .TRUE.
         CALL FILTER(T3,MAS,T1,T2,AM1,'MATMAS          ',
     *               1.D0,SBID,SBID,SBID,SBID,SBID,SBID,
     *               MESH,MSK,MASKEL,LISHHO)
      ENDIF
C   
      CALL OS('X=Y     ',HHO,T3,SBID,BID)
C
C -------------------------------------------------------------
C CALCUL DES CONTRAINTES SXX, SXY ET SYY
C -------------------------------------------------------------
C
      CALL OS('X=Y/Z   ' , T1 , CG , C  , BID )
      DO I=1,NPOIN
         COCO=COS(INCI%R(I))*COS(INCI%R(I))
         COSI=COS(INCI%R(I))*SIN(INCI%R(I))
         SISI=SIN(INCI%R(I))*SIN(INCI%R(I))
         COE=GRAV*HHO%R(I)*HHO%R(I)/8.D0
C
C ON A LE COEFF. 1/8 CI-DESSUS CAR HHO REPRESENTE LA HAUTEUR 
C DE HOULE ENERGETIQUE EN HOULE REGULIERE
C
         SXX%R(I)= SXX%R(I) + (T1%R(I)*(1.D0+COCO)-0.5D0)*COE
         SXY%R(I)= SXY%R(I) + (T1%R(I)*COSI)*COE
         SYY%R(I)= SYY%R(I) + (T1%R(I)*(1.D0+SISI)-0.5D0)*COE
      END DO
      END IF
C
C
C=======================================================================
C GRADIENTS SPATIAUX DES CONTRAINTES DE RADIATION.
C=======================================================================
C
C  -----------------------------------------------
C  LISSAGES EVENTUELS DES CONTRAINTES DE RADIATION
C  -----------------------------------------------
C
      LISRAD = 3
C
      CALL OS('X=Y     ',T3,SXX,SBID,BID)
      IF(LISRAD.GT.0) THEN
         MAS = .TRUE.
         CALL FILTER(T3,MAS,T1,T2,AM1,'MATMAS          ',
     *               1.D0,SBID,SBID,SBID,SBID,SBID,SBID,
     *               MESH,MSK,MASKEL,LISRAD)
      ENDIF
      CALL OS('X=Y     ',SXX,T3,SBID,BID)
C
      CALL OS('X=Y     ',T3,SXY,SBID,BID)
      IF(LISRAD.GT.0) THEN
         MAS = .TRUE.
         CALL FILTER(T3,MAS,T1,T2,AM1,'MATMAS          ',
     *               1.D0,SBID,SBID,SBID,SBID,SBID,SBID,
     *               MESH,MSK,MASKEL,LISRAD)
      ENDIF
      CALL OS('X=Y     ',SXY,T3,SBID,BID)
C
      CALL OS('X=Y     ',T3,SYY,SBID,BID)
      IF(LISRAD.GT.0) THEN
         MAS = .TRUE.
         CALL FILTER(T3,MAS,T1,T2,AM1,'MATMAS          ',
     *               1.D0,SBID,SBID,SBID,SBID,SBID,SBID,
     *               MESH,MSK,MASKEL,LISRAD)
      ENDIF
      CALL OS('X=Y     ',SYY,T3,SBID,BID)
C
C FIN DU LISSAGE DES CONTRAINTES DE RADIATION
C -------------------------------------------------------
C
C=======================================================================
C FORCES MOTRICES FX ET FY POUR LES COURANTS DE HOULE
C=======================================================================
C
      CALL VECTOR(T1 , '=' , 'MASBAS          ' , IELM ,
     *            1.D0 , T3 , T3 , T3 , T3 , T3 , T3 ,
     *            MESH , MSK , MASKEL )
C
      CALL VECTOR(T2 , '=' , 'GRADF          X' , IELM ,
     *            1.D0 , SXX, T4 , T4 , T4 , T4 , T4 ,
     *            MESH , MSK , MASKEL )
      CALL OS('X=Y/Z   ',T2,T2,T1,BID)
C
      CALL VECTOR
     * (T3,'=','GRADF          Y',IELM,1.D0,SXY,T4,T4,T4,T4,T4,
     *  MESH , MSK , MASKEL )
      CALL OS('X=Y/Z   ',T3,T3,T1,BID)
C     ----------------------------------
C     FORCE FX = - (dSXX/dx + dSXY/dy) / h
C     ----------------------------------
      CALL OS('X=Y+Z   ',FX,T2,T3,BID)
      CALL OS('X=CY/Z  ',FX,FX,H,-1.D0)
C
C     ----------------------------------
C
      CALL VECTOR(T1 , '=' , 'MASBAS          ' , IELM ,
     *            1.D0 , T3 , T3 , T3 , T3 , T3 , T3 ,
     *            MESH , MSK , MASKEL )
C
      CALL VECTOR
     * (T2,'=','GRADF          X',IELM,1.D0,SXY,T4,T4,T4,T4,T4,
     *  MESH , MSK , MASKEL )
      CALL OS('X=Y/Z   ',T2,T2,T1,BID)
C
      CALL VECTOR
     * (T3,'=','GRADF          Y',IELM,1.D0,SYY,T4,T4,T4,T4,T4,
     *  MESH , MSK , MASKEL )
      CALL OS('X=Y/Z   ',T3,T3,T1,BID)
C
C     ----------------------------------
C     FORCE FY = - (dSXY/dx + dSYY/dy) / h
C     ----------------------------------
      CALL OS('X=Y+Z   ',FY,T2,T3,BID)
      CALL OS('X=CY/Z  ',FY,FY,H,-1.D0)
C
C=======================================================================
C
      RETURN
      END
