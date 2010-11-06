C                       *****************
                        SUBROUTINE CALRES
C                       *****************
C
C***********************************************************************
C
C  ARTEMIS VERSION 5.1   04/06/99    D. AELBRECHT (LNH) 01 30 87 74 12
C
C  LINKED TO BIEF VERS. 5.0          J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  : CALCULE LA HAUTEUR ET LA PHASE DE LA HOULE,
C                 LES VITESSES ET LA COTE DE LA SURFACE LIBRE
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : ARTEMIS
C
C  SOUS-PROGRAMME APPELE : OS, VECTOR
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_ARTEMIS
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER I
C
      DOUBLE PRECISION PI,RADDEG
      DOUBLE PRECISION ZERO, BID
      DOUBLE PRECISION A1, A2 ,ALPHA0, D1, D2, PHI, PHI1, PHI2 ,MODPHI
      DOUBLE PRECISION TETA01, XU1 ,XU2, XV1, XV2, WT0
C
      INTRINSIC SQRT, ATAN2, DMOD, ABS, COS, SIN
C
C-----------------------------------------------------------------------
C
C
      PARAMETER (ZERO = 1.D-10)
      PARAMETER (PI = 3.1415926535897932384626433D0)
      PARAMETER (RADDEG = 57.29577951D0)
C
C=======================================================================
C HAUTEUR DE HOULE
C=======================================================================
C
      CALL OS( 'X=N(Y,Z)', T1, PHIR, PHII , BID             )
      CALL OS( 'X=CY    ', HHO  ,T1, SBID ,(2.D0*OMEGA/GRAV))
C
C=======================================================================
C PHASE DU POTENTIEL  (EN RADIAN)
C=======================================================================
C
      DO 10 I=1,NPOIN
         IF (T1%R(I).LT.ZERO) THEN
            PHAS%R(I) = 0.D0
         ELSE
            PHAS%R(I) = ATAN2( PHII%R(I),PHIR%R(I) )
         ENDIF
10    CONTINUE
C
C=======================================================================
C COTE DE LA SURFACE LIBRE
C=======================================================================
C
      DO 20 I=1,NPOIN
         S%R(I) = -OMEGA/GRAV*PHII%R(I) + H%R(I) + ZF%R(I)
20    CONTINUE
C
C=======================================================================
C VITESSES EN SURFACE  (A T=0 ET A T=OMEGA/4)
C=======================================================================
C
C CALCUL DES GRADIENTS ( DE PHIR ET PHII)
C
C
      CALL VECTOR(U0 , '=' , 'GRADF          X' , IELM ,
     *            1.D0 , PHIR , SBID, SBID , SBID , SBID , SBID ,
     *            MESH , MSK , MASKEL )
C
      CALL VECTOR(V0 , '=' , 'GRADF          Y' , IELM ,
     *            1.D0 , PHIR , SBID , SBID , SBID , SBID , SBID ,
     *            MESH , MSK , MASKEL )
C
C     ON MET L'ANCIENNE VARIABLE U1 DANS T3
C     CAR ELLE SERT POUR CALCULER INCI
C
      CALL VECTOR(T3 , '=' , 'GRADF          X' , IELM ,
     *            1.D0 , PHII , SBID , SBID , SBID , SBID , SBID ,
     *            MESH , MSK , MASKEL )
C
C     ON MET L'ANCIENNE VARIABLE V1 DANS T4
C     CAR ELLE SERT POUR CALCULER INCI
C
      CALL VECTOR(T4 , '=' , 'GRADF          Y' , IELM ,
     *            1.D0 , PHII , SBID , SBID , SBID , SBID , SBID ,
     *            MESH , MSK , MASKEL )
C
      CALL VECTOR(T1 , '=' , 'MASBAS          ' , IELM ,
     *            1.D0 , SBID , SBID , SBID , SBID , SBID , SBID ,
     *            MESH , MSK , MASKEL )
C
      CALL OS( 'X=Y/Z   ' , U0    , U0    , T1 , BID )
      CALL OS( 'X=Y/Z   ' , V0    , V0    , T1 , BID )
      CALL OS( 'X=Y/Z   ' , T3    , T3    , T1 , BID )
      CALL OS( 'X=Y/Z   ' , T4    , T4    , T1 , BID )
C
C=======================================================================
C CALCUL DE L'INCIDENCE DE LA HOULE
C=======================================================================
C
C ON A : U0 (D(PHIR)/DX) : A      U1 (D(PHII)/DX): B
C        V0 (D(PHIR)/DY) : C      V1 (D(PHII)/DY): D
C PASSAGE DE U= A COS WT + B SIN WT  A : U = A1 COS ( WT - PHI1)
C            V= C COS WT + D SIN WT      V = A2 COS ( WT - PHI2)
C
      DO 30 I=1,NPOIN
        A1 = SQRT ( U0%R(I)*U0%R(I) + T3%R(I)*T3%R(I) )
        PHI1 = ATAN2( T3%R(I),U0%R(I) )
        A2 = SQRT ( V0%R(I)*V0%R(I) + T4%R(I)*T4%R(I) )
        PHI2 = ATAN2( T4%R(I),V0%R(I) )
C
C ECRITURE SOUS FORME : U = A1 COS ( (WT - PHI1))
C                       V = A2 COS ( (WT - PHI1) - PHI )
C ON CHOISIT PHI ENTRE 0 ET 2*PI
C
        PHI = PHI2 - PHI1
        IF (PHI.LT.0.D0)   PHI = PHI+2.D0*PI
C
C RECHERCHE DE LA DIRECTION ET DU (WT0) OU LE GRAND AXE DE L'ELLIPSE
C EST ATTEINT.
C TRAITEMENT DE CAS PARTICULIERS (POLARISATION RECTILIGNE)
C
        MODPHI = DMOD( PHI, PI )
        IF ( (MODPHI.LT.ZERO).OR.((PI-MODPHI).LT.ZERO) ) THEN
          WT0 = PHI1
          IF ( (PHI.LT.2D0*ZERO).OR.((2.D0*PI-PHI).LT.2D0*ZERO) )THEN
            ALPHA0 = ATAN2( A2,A1 )
          ELSE
C           ON A : (ABS(PHI-PI).LT.2D0*ZERO)
            ALPHA0 = 2.D0*PI - ATAN2( A2,A1 )
          ENDIF
        ELSE
C CAS GENERAL : POLARISATION ELLIPTIQUE
C ON A : TAN(2*(WT0 - PHI1)) = A2**2*SIN(2*PHI)/(A1**2+A2**2*COS(2*PHI))
          TETA01 = ATAN2( (A2*A2*SIN(2*PHI)) ,
     *                    (A1*A1 + A2*A2*COS(2*PHI)) ) / 2.D0
          XU1 = A1 * COS ( TETA01)
          XV1 = A2 * COS ( TETA01 - PHI )
          XU2 = -A1 * SIN ( TETA01)
          XV2 = -A2 * SIN ( TETA01 - PHI )
          D1 = XU1*XU1 + XV1*XV1
          D2 = XU2*XU2 + XV2*XV2
          IF (D2.GT.D1) THEN
             TETA01 = TETA01 + PI/2.D0
             XU1    = XU2
             XV1    = XV2
          ENDIF
          WT0    = TETA01 + PHI1
          ALPHA0 = ATAN2( XV1,XU1 )
        ENDIF
        INCI%R(I)  = ALPHA0
        T2%R(I) = WT0
 30   CONTINUE
C
C SURFACE LIBRE EN PHASE AVEC ALPHA0
C ON DONNE UN SIGNE POSTIF A L'INCIDENCE SI LA SURFACE LIBRE EST
C POSITIVE.
C
      DO 40 I=1,NPOIN
         A1 = -(PHII%R(I)*COS(T2%R(I))-PHIR%R(I)*SIN(T2%R(I)))
         IF (A1.LT.0.D0) THEN
           IF (INCI%R(I).GE.0.D0) THEN
             INCI%R(I) = INCI%R(I) - PI
           ELSE
             INCI%R(I) = INCI%R(I) + PI
           ENDIF
         ENDIF
 40   CONTINUE
C
C=======================================================================
C
      RETURN
      END
