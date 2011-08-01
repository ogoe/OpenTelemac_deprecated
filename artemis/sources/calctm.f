C                       *****************
                        SUBROUTINE CALCTM
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
C     FONCTION  : CALCULE LES DIFFERENTES VALEURS DE PERIODES
C                 MOYENNES DU SPECTRE : 
C                    T01 = m0/m1
C                    T02 = SQRT(m0/m2)
C                    TM   
C                 (VOIR DEFINITIONS DANS LA LISTE
C                 DES PARAMETRES D'ETATS DE MER ETABLIS PAR L'AIRH)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     T01        |<-->|  MOMENT D'ORDRE 1 PUIS VALEUR DE T01         |
C |     T02        |<-->|  MOMENT D'ORDRE 2 PUIS VALEUR DE T02         |
C |     TM         |<-->|  MOMENT D'ORDRE 1 EN PERIODE PUIS VALEUR TM  |
C |     T1         | -->|  TABLEAU DE TRAVAIL                          |
C |     T2         | -->|  TABLEAU DE TRAVAIL                          |
C |     T3         | -->|  TABLEAU DE TRAVAIL                          |
C |     HHO        | -->|  HAUTEURS DE HOULE SIGNIFICATIVES            |
C |     HALE       | -->|  ICI, ENERGIE DE HOULE SANS LA DERN. COMPOS. | 
C |     NPALE      | -->|  NOMBRE DE PERIODES DE DISCRETISATION        |
C |                |    |  DU SPECTRE DE HOULE                         |
C |     NDALE      | -->|  NOMBRE DE DIRECTIONS DE DISCRETISATION      |
C |                |    |  DU SPECTRE DE HOULE                         |
C |     PER        | -->|  PERIODE DE HOULE EN COURS DE CALCUL         |
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
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_ARTEMIS
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      DOUBLE PRECISION PONDER,RADDEG
C
      DOUBLE PRECISION BID
C
      INTRINSIC SQRT, ATAN2, DMOD, ABS, COS, SIN
C
C-----------------------------------------------------------------------
C
C STRUCTURES
C
C
C-----------------------------------------------------------------------
C
      PONDER = 1.D0/DBLE(NDALE*NPALE)
      RADDEG = 180.D0/3.141592654D0
C
C=======================================================================
C CALCUL DU MOMENT m0 MIS DANS T2
C=======================================================================
C
      CALL OS( 'X=Y**C  ', T1 , HHO  , SBID  , 2.D0 )
      CALL OS( 'X=CX    ', T1 , SBID  , SBID  , PONDER )
      CALL OS( 'X=Y+Z   ', T2 , HALE , T1   , BID )
C
C=======================================================================
C PERIODE T01 = m0 / m1 
C=======================================================================
C
      CALL OS( 'X=Y     ', T3 , T01 , SBID  , BID )
      CALL OS( 'X=Y/Z   ', T01, T2  , T3   , BID )
C
C=======================================================================
C PERIODE T02 = SQRT( m0 / m2 )
C=======================================================================
C
      CALL OS( 'X=Y     ', T3 , T02 , SBID , BID )
      CALL OS( 'X=Y/Z   ', T1 , T2  , T3  , BID )
      CALL OS( 'X=SQR(Y)', T02, T1  , SBID , BID )
C
C=======================================================================
C PERIODE TM =  mt1 / m0
C=======================================================================
C
      CALL OS( 'X=Y     ', T3 , TM , SBID , BID )
      CALL OS( 'X=Y/Z   ', TM , T3 , T2  , BID )
C
C
C=======================================================================
C DIRECTION MOYENNE INCI 
C=======================================================================
C
      CALL OS( 'X=A(Y,Z)',INCI , MSIN, MCOS , BID )
C
      RETURN
      END
