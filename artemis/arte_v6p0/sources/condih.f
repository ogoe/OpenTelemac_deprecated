C                       *****************
                        SUBROUTINE CONDIH
C                       *****************
C
C***********************************************************************
C
C  ARTEMIS VERSION 5.1   02/06/99    D. AELBRECHT (LNH) 01 30 87 74 12
C
C  LINKED TO BIEF VERS. 5.0          J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  : INITIALISATION DES TABLEAUX DES GRANDEURS PHYSIQUES
C
C-----------------------------------------------------------------------
C
C APPELE PAR : ARTEMIS
C
C SOUS-PROGRAMME APPELE : OS
C
C FONCTION APPELEE : INCLUS,MAJUS
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
      DOUBLE PRECISION COTE
      DOUBLE PRECISION PI,BID,DHTEST
      PARAMETER( PI = 3.1415926535897932384626433D0 )
C
      INTRINSIC SINH, SQRT
C
C-----------------------------------------------------------------------
C
      CALL MAJUS(CDTINI)
C
C-----------------------------------------------------------------------
C
C   INITIALISATION DE H , LA HAUTEUR D'EAU
C
      IF(INCLUS(CDTINI,'COTE NULLE').OR.
     *   INCLUS(CDTINI,'ZERO ELEVATION') ) THEN
        COTE = 0.D0
        CALL OS( 'X=C     ' , H , SBID , SBID , COTE )
        CALL OS( 'X=X-Y   ' , H , ZF  , SBID , BID  )
      ELSEIF(INCLUS(CDTINI,'COTE CONSTANTE').OR.
     *       INCLUS(CDTINI,'CONSTANT ELEVATION') ) THEN
        COTE = COTINI
        CALL OS( 'X=C     ' , H , SBID , SBID , COTE )
        CALL OS( 'X=X-Y   ' , H , ZF  , SBID , BID  )
      ELSEIF(INCLUS(CDTINI,'HAUTEUR NULLE').OR.
     *       INCLUS(CDTINI,'ZERO DEPTH') ) THEN
        CALL OS( 'X=C     ' , H , SBID , SBID , 0.D0 )
      ELSEIF(INCLUS(CDTINI,'HAUTEUR CONSTANTE').OR.
     *       INCLUS(CDTINI,'CONSTANT DEPTH') ) THEN
        CALL OS( 'X=C     ' , H , SBID , SBID , HAUTIN )
      ELSEIF(INCLUS(CDTINI,'PARTICULIERES').OR.
     *       INCLUS(CDTINI,'SPECIAL')        ) THEN
C  ZONE A MODIFIER
        IF(LNG.EQ.1) WRITE(LU,10)
        IF(LNG.EQ.2) WRITE(LU,11)
10      FORMAT(1X,'CONDIH : AVEC DES CONDITIONS INITIALES PARTICULIERES'
     *         ,/,'         VOUS DEVEZ MODIFIER CONDIH')
11      FORMAT(1X,'CONDIH : WITH SPECIAL INITIAL CONDITIONS'
     *         ,/,'         YOU HAVE TO MODIFY CONDIH')
        CALL PLANTE(0)
        STOP
C  FIN DE LA ZONE A MODIFIER
      ELSE
        IF(LNG.EQ.1) WRITE(LU,20) CDTINI
        IF(LNG.EQ.2) WRITE(LU,21) CDTINI
20      FORMAT(1X,'CONDIH : CONDITION INITIALE INCONNUE :',/,A72)
21      FORMAT(1X,'CONDIH : UNKNOWN INITIAL CONDITION :',/,A72)
        CALL PLANTE(0)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C  CLIPPING DE H (PAS DE VALEURS INFERIEURES A 1.D-2)
C
      CALL CLIP(H,1.D-2,.TRUE.,1.D6,.FALSE.,NPOIN)
C
C-----------------------------------------------------------------------
C
C   CALCUL DU NOMBRE D'ONDE K
C   ON UTILISE UNE FORMULE EXPLICITE (CF L'EXCELLENT RAPPORT EDF DE
C   F. DHELLEMMES 'PRECIS SUR LES VAGUES' )
C


      OMEGA = 2.D0*PI/PER
      CALL OS('X=CY    ', T1 , H , SBID , OMEGA**2/GRAV )
C
C     INITIALISATION MAXIMISANTE DE DHTEST
C 
      DHTEST = 1.D6
C
      DO 100 I=1,NPOIN
         T2%R(I) = 1.D0 + T1%R(I) *( 0.6522D0 +
     *                    T1%R(I) *( 0.4622D0 +
     *                    T1%R(I) *
     *                    T1%R(I) *( 0.0864D0 +
     *                    T1%R(I) *( 0.0675D0 ) )))
         T2%R(I) = SQRT( T1%R(I)*(T1%R(I) + 1.D0/T2%R(I)) )
         K%R(I)  = T2%R(I)/H%R(I)
         DHTEST  = MIN( DHTEST , H%R(I) )
100   CONTINUE
C
C   TEST POUR DETECTER S'IL Y A EU CLIPPING SUR H
C
      IF (DHTEST.LE.1.01D-2) THEN
         IF(LNG.EQ.1) WRITE(LU,120)
         IF(LNG.EQ.2) WRITE(LU,121)
120      FORMAT(1X,'CONDIH : ATTENTION !! VOUS AVEZ ATTEINT LE SEUIL '
     *          ,/,'         MINI DE HAUTEUR D''EAU (1 cm).'
     *          ,/,'         VERIFIEZ BATHY OU CONDITIONS INITIALES')
121      FORMAT(1X,'CONDIH : WARNING !! YOU REACHED MINIMUM THRESHOLD'
     *          ,/,'         FOR WATER DEPTH (1 cm). CHECK THE'
     *          ,/,'         BATHYMETRY OR INITIAL CONDITIONS')
      ENDIF
C
C-----------------------------------------------------------------------
C
C   CALCUL DE LA VITESSE DE PHASE
C
      CALL OS('X=CY    ', T1    , K     , SBID , 1.D0/OMEGA )
      CALL OS('X=1/Y   ', C     , T1    , SBID , BID        )
C
C-----------------------------------------------------------------------
C
C   CALCUL DE LA VITESSE DE GROUPE
C
      DO 200 I=1,NPOIN
         CG%R(I) = C%R(I)/2.D0 *
     *             (1.D0 + 2.D0*K%R(I)*H%R(I)/SINH(2.D0*K%R(I)*H%R(I)))
200   CONTINUE
C     
    

C-----------------------------------------------------------------------
C
      RETURN
      END
