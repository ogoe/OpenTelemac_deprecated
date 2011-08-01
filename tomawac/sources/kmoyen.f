C                       *****************
                        SUBROUTINE KMOYEN
C                       *****************
C
     *( XKMOY , XK    , F     , FREQ  , DFREQ , TAILF , NF    , NPLAN ,
     *  NPOIN2, AUX1  , AUX2  , AUX3  )
C
C**********************************************************************
C  TOMAWAC - V1.0    P. THELLIER & M. BENOIT (EDF/DER/LNH)  -  04/04/95
C**********************************************************************
C
C  FONCTION : CALCUL DU NOMBRE D'ONDE MOYEN EN TOUS LES POINTS DU 
C  ********   MAILLAGE SPATIAL 2D.
C
C  ARGUMENTS :
C  ***********
C  +-------------+----+--------------------------------------------+
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C  +-------------+----+--------------------------------------------+
C  ! XKMOY(-)    !<-- ! TABLEAU DES NOMBRES D'ONDE MOYEN           !
C  ! XK(-,-)     ! -->! TABLEAU DES NOMBRES D'ONDE                 !
C  ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL                       !
C  ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !
C  ! DFREQ(-)    ! -->! TABLEAU DES PAS DE FREQUENCE               !
C  ! TAILF       ! -->! FACTEUR DE QUEUE (TAILF = 4 OU 5)          !
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C  ! AUX1(-)     !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C  ! AUX2(-)     !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C  ! AUX3(-)     !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C  +-------------+----+--------------------------------------------+
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C  +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  SEMIMP, PRE2D
C  ********    - PROGRAMME(S) APPELE(S) :    -
C
C  REMARQUES :
C  ***********
C  - LA PARTIE HAUTES-FREQUENCES DU SPECTRE N'EST PRISE EN COMPTE QUE
C    SI LE FACTEUR DE QUEUE (TAILF) EST STRICTEMENT SUPERIEUR A 1.
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NF    , NPLAN , NPOIN2
      DOUBLE PRECISION TAILF
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF), XK(NPOIN2,NF)
      DOUBLE PRECISION FREQ(NF)  , DFREQ(NF) , XKMOY(NPOIN2)
      DOUBLE PRECISION AUX1(NPOIN2) , AUX2(NPOIN2) , AUX3(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  IPLAN , JF    , IP
      DOUBLE PRECISION COEFF , PI    , SEUIL , CTE1  , CTE2  , AUX4
C
C
      PI = 3.141592654D0
      SEUIL = 1.D-20         
      COEFF = SQRT(9.806D0)/(2.D0*PI)
      DO 30 IP = 1,NPOIN2
        AUX1(IP) = 0.D0
        AUX2(IP) = 0.D0
   30 CONTINUE
C
C.....SOMMATIONS SUR LA PARTIE DISCRETISEE DU SPECTRE.
C     """"""""""""""""""""""""""""""""""""""""""""""""
      DO 20 JF = 1,NF
        AUX4=DFREQ(JF)
C
        DO 15 IP=1,NPOIN2
          AUX3(IP) = 0.D0
   15   CONTINUE
        DO 10 IPLAN = 1,NPLAN
          DO 5 IP=1,NPOIN2
            AUX3(IP) = AUX3(IP) + F(IP,IPLAN,JF)
    5     CONTINUE
   10   CONTINUE
C
        DO 25 IP = 1,NPOIN2
          AUX1(IP)=AUX1(IP)+AUX3(IP)*AUX4
          AUX2(IP)=AUX2(IP)+AUX3(IP)/SQRT(XK(IP,JF))*AUX4
   25   CONTINUE
C
   20 CONTINUE
C
C.....PRISE EN COMPTE EVENTUELLE DE LA PARTIE HAUTES FREQUENCES.
C     """"""""""""""""""""""""""""""""""""""""""""""""""""""""""
      IF (TAILF.GT.1.D0) THEN
        CTE1=FREQ(NF)/(TAILF-1.D0)
        CTE2=COEFF/TAILF
      ELSE
        CTE1=0.D0
        CTE2=0.D0
      ENDIF
      DO 45 IP=1,NPOIN2
        AUX1(IP) = AUX1(IP) + AUX3(IP)*CTE1
        AUX2(IP) = AUX2(IP) + AUX3(IP)*CTE2
   45 CONTINUE
C
C.....CALCUL DU NOMBRE D'ONDE MOYEN. 
C     """"""""""""""""""""""""""""""
      DO 50 IP=1,NPOIN2
        IF (AUX2(IP).LT.SEUIL) THEN
          XKMOY(IP) = 1.D0
        ELSE
          XKMOY(IP) = (AUX1(IP)/AUX2(IP))**2
        ENDIF
   50 CONTINUE
C
      RETURN
      END
