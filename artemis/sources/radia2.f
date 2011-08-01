C                       *****************
                        SUBROUTINE RADIA2
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
C |    LISHHO      |<-->|  lissage de la hauteur de houle              |
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
      USE INTERFACE_ARTEMIS, EX_RADIA2=> RADIA2 
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
C VARIABLES INTERNES A RADIA2
C
      DOUBLE PRECISION COE , COCO, COSI, SISI , Z(1)
      INTEGER          LISRAD
C
      LOGICAL MAS 
C
      INTRINSIC SQRT, ATAN2, DMOD, ABS, COS, SIN
C
C
C=======================================================================
C CONTRAINTES DE RADIATION........METHODE 2 (IDENTIQUE A TOMAWAC)
C=======================================================================
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
         COE=GRAV*HALE%R(I)*HALE%R(I)/16.D0
C
C ON A LE COEFF. 1/16 CI-DESSUS CAR HHO REPRESENTE LA HAUTEUR 
C DE HOULE SIGNIFICATIVE ET NON ENERGETIQUE EN HOULE ALEATOIRE
C
         SXX%R(I)= (T1%R(I)*(1.D0+COCO)-0.5D0)*COE
         SXY%R(I)= (T1%R(I)*COSI)*COE
         SYY%R(I)= (T1%R(I)*(1.D0+SISI)-0.5D0)*COE
      END DO
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
      CALL VECTOR
     * (T2,'=','GRADF          X',IELM,1.D0,SXX,T4,T4,T4,T4,T4,
     *  MESH , MSK , MASKEL )
      CALL OS('X=Y/Z   ',T2,T2,T1,BID)
C
      CALL VECTOR
     * (T3,'=','GRADF          Y',IELM,1.D0,SXY,T4,T4,T4,T4,T4,
     *  MESH , MSK , MASKEL )
      CALL OS('X=Y/Z   ',T3,T3,T1,BID)
C     ------------------------------------
C     FORCE FX = - (dSXX/dx + dSXY/dy) / h
C     ------------------------------------
      CALL OS('X=Y+Z   ',FX,T2,T3,BID)
      CALL OS('X=CY/Z  ',FX,FX,H,-1.D0)
C
C     ------------------------------------
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
C     ------------------------------------
C     FORCE FY = - (dSXY/dx + dSYY/dy) / h
C     ------------------------------------
      CALL OS('X=Y+Z   ',FY,T2,T3,BID)
      CALL OS('X=CY/Z  ',FY,FY,H,-1.D0)
C
C=======================================================================
C
      RETURN
      END
