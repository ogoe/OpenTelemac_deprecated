C                       ***************
                        FUNCTION QBBJ78
C                       ***************
C
     *( B     , IQBBJ )
C
C**********************************************************************
C  TOMAWAC - V1.1         F.BECQ & M. BENOIT (EDF/DER/LNH)  -  14/02/96
C**********************************************************************
C
C  FONCTION : CALCUL DE LA PROBABILITE DE DEFERLEMENT : QB EST UTILISEE
C  ********** DANS LE MODELE DE BATTJES ET JANSSEN (1978).
C
C  ARGUMENTS :
C  ***********
C +-------------+----+--------------------------------------------+
C ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C +-------------+----+--------------------------------------------+
C ! B           ! -->! (HE/HM)                                    !
C ! IQBBJ       ! -->! INDICE DE LA METHODE UTILISEE              !
C +-------------+----+--------------------------------------------+
C ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  QBREK1
C  ********    - PROGRAMME(S) APPELE(S) :    -
C
C  REFERENCES : - BATTJES ET JANSSEN (1978) : ENERGY LOSS AND SET-UP
C  ************   DUE TO BREAKING OF RANDOM WAVES. (ICCE'78)
C
C               - DINGEMANS (1983) : VERIFICATION OF NUMERICAL WAVE
C                 PROPAGATION MODELS WITH FIELD MEASUREMENTS.
C                 CREDIZ VERIFICATION HARINGVLIET.
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      DOUBLE PRECISION QBBJ78, B
      INTEGER  IQBBJ
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      DOUBLE PRECISION F     , CB    , EPS   , QMAX  , QMIN  , Q0
      DOUBLE PRECISION B2    , EXPO
C
C
      EPS = 1.D-7
C
      IF (B.GE.1.D0) THEN
        QBBJ78 = 1.D0
        RETURN
      ENDIF
C
      IF(IQBBJ.EQ.0) THEN
C       =========================
C       RESOLUTION PAR DICHOTOMIE
C       =========================
        QMIN  = 0.D0
        QMAX  = B
   10   CONTINUE
        QBBJ78 = (QMIN+QMAX)/2.D0
        F      = 1.D0 - QBBJ78 + B*B*DLOG(QBBJ78)
        IF (ABS(F).LT.EPS) RETURN
        IF (F.GT.0.D0) THEN
           QMAX = QBBJ78
        ELSE
           QMIN = QBBJ78
        ENDIF
        GOTO 10
C
      ELSEIF(IQBBJ.EQ.1) THEN
C       ======================================================
C       FORMULATION EXPLICITE 1 (INSPIREE DE CREDIZ VERSION-1)
C       ======================================================
      CB = 0.5D0
        IF (B.GE.CB) THEN
          QBBJ78 = ((B-CB)/(1.D0-CB))**2
        ELSE
          QBBJ78 = 0.D0
        ENDIF
C
      ELSEIF(IQBBJ.EQ.2) THEN
C       ======================================================
C       FORMULATION EXPLICITE 2 (INSPIREE DE CREDIZ VERSION-2)
C       ======================================================
        CB = 0.3D0
        IF (B.LT.CB) THEN
          QBBJ78 = 0.D0
        ELSEIF (B.LT.0.5D0) THEN
          B2     = B**2
          EXPO   = DEXP(-1.D0/B2)
          QBBJ78 = B2*EXPO/(B2-EXPO)
        ELSEIF (B.LT.0.9D0) THEN
          Q0     = (2.D0*B-1.D0)**2
          B2     = B**2
          EXPO   = DEXP((Q0-1.D0)/B2)
          QBBJ78 = Q0 - B2*(Q0-EXPO)/(B2-EXPO)
        ELSE
          QBBJ78 = (2.D0*B-1.D0)**2
        ENDIF
C
      ELSEIF(IQBBJ.EQ.3) THEN
C       ======================================================
C       FORMULATION EXPLICITE 3 (INSPIREE DE CREDIZ VERSION-3)
C       ======================================================
        QBBJ78 = 2.4D0*B**7
C
      ENDIF
C
      RETURN
      END
