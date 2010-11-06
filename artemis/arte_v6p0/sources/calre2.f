C                       *****************
                        SUBROUTINE CALRE2
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
C     FONCTION  : CALCULE DES GRANDEURS CARACTERISTIQUES MOYENNES
C                 DU SPECTRE DE HOULE EN HOULE ALEATOIRE : 
C                    K : NOMBRE D'ONDE MOYEN
C                    C : CELERITE DE PHASE MOYENNE
C                    CG : CELERITE DE GROUPE MOYENNE  
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |                |    |                                              |
C | T01            |--> | PERIODE MOYENNE ISSUE DU MOMENT D'ORDRE 1    |
C | T02            |--> | PERIODE MOYENNE ISSUE DU MOMENT D'ORDRE 2    |
C | TM             |--> | PERIODE MOYENNE ISSUE DU MOMENT D'ORDRE 1    |
C |                |    | EN PERIODE                                   |
C | T1,T2          |--> | TABLEAUX DE TRAVAIL                          |
C | H              |--> | HAUTEUR D'EAU AU REPOS                       |
C | K              |<-- | NOMBRE D'ONDE MOYEN                          |
C | C              |<-- | CELERITE DE PHASE MOYENNE                    |
C | CG             |<-- | CELERITE DE GROUPE MOYENNE                   |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : ARTEMI
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
      DOUBLE PRECISION PI,DEUXPI, DHTEST
C
      DOUBLE PRECISION BID
C
      INTRINSIC SQRT, ATAN2, DMOD, ABS, COS, SIN
C
C-----------------------------------------------------------------------
C
      PI = 3.1415926535897932384626433D0 
      DEUXPI = 2.D0*PI
      GRAV = 9.81D0
C-----------------------------------------------------------------------
C
C   CALCUL DU NOMBRE D'ONDE K
C   ON UTILISE UNE FORMULE EXPLICITE (CF L'EXCELLENT RAPPORT EDF DE
C   F. DHELLEMMES 'PRECIS SUR LES VAGUES' )
C
C-----------------------------------------------------------------------
C
C ON MET OMEGA MOYEN DANS T1
C
      CALL OS( 'X=1/Y   ', T1 , T01  , SBID  , BID )
      CALL OS( 'X=CX    ', T1 , SBID  , SBID  , DEUXPI )
C
C ON CALCULE OMEGA**2 * H / GRAV
C
      CALL OS( 'X=Y**C  ', T2 , T1 , SBID , 2.D0 )
      CALL OS( 'X=CXY   ', T2 , H  , SBID , 1.D0/GRAV )
C
C     INITIALISATION MAXIMISANTE DE DHTEST
C
      DHTEST = 1.D6
C
      DO 100 I=1,NPOIN
         T1%R(I) = 1.D0 + T2%R(I) *( 0.6522D0 +
     *                  T2%R(I) *( 0.4622D0 +
     *                  T2%R(I) *
     *                  T2%R(I) *( 0.0864D0 +
     *                  T2%R(I) *( 0.0675D0 ) )))
         T1%R(I) = SQRT( T2%R(I)*(T2%R(I) + 1.D0/T1%R(I)) )
         K%R(I)  = T1%R(I)/H%R(I)
         DHTEST  = MIN( DHTEST , H%R(I) )
100   CONTINUE
C
C
C=======================================================================
C CALCUL DE C MOYEN
C=======================================================================
C
      CALL OS( 'X=1/Y   ', T1 , T01  , SBID  , BID )
      CALL OS( 'X=CX    ', T1 , SBID  , SBID  , DEUXPI )
      CALL OS( 'X=Y/Z   ', C  , T1   , K    , BID )
C
C
C=======================================================================
C CALCUL DE CG MOYEN
C=======================================================================
C
      DO 200 I=1,NPOIN
         CG%R(I) = C%R(I)/2.D0 *
     *             (1.D0 + 2.D0*K%R(I)*H%R(I)/SINH(2.D0*K%R(I)*H%R(I)))
200   CONTINUE
C                                                                       
      RETURN                                                            
      END                                                               
