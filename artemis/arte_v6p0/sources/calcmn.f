C                       ******************
                        SUBROUTINE CALCMN
C                       ******************
C
C***********************************************************************
C
C  ARTEMIS VERSION 5.1   04/06/99    D. AELBRECHT (LNH) 01 30 87 74 12
C
C  LINKED TO BIEF VERS. 5.0          J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  : CALCULE LES VALEURS APPROCHEES DES MOMENTS m0, m1, ET
C                 m2 DU SPECTRE DE HOULE POUR CALCULER LES PERIODES
C                 MOYENNES DU SPECTRE ET LA DIRECTION MOYENNE
C                 (VOIR DEFINITIONS DANS LA LISTE DES PARAMETRES
C                  ETABLIS PAR L'AIRH)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
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
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_ARTEMIS
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      DOUBLE PRECISION FREQ, FREQ2, PONDER
C
      DOUBLE PRECISION BID
C
      INTRINSIC SQRT, ATAN2, DMOD, ABS, COS, SIN
C
C-----------------------------------------------------------------------
C
C STRUCTURES
C
C-----------------------------------------------------------------------
C
      FREQ = 1.D0/PER
      FREQ2 = FREQ * FREQ
      PONDER = 1.D0/DBLE(NDALE*NPALE)
C
C=======================================================================
C MOMENT m1 = INTEGRALE ( f * S(f) * df )
C=======================================================================
C
      CALL OS( 'X=Y**C  ', T1, HHO , SBID , 2.D0 )
      CALL OS( 'X=CY    ', T2, T1  , SBID , FREQ )
      CALL OS( 'X=CX    ', T2, SBID , SBID , PONDER )
      CALL OS( 'X=X+Y   ', T01, T2 , SBID , 1.D0 )
C
C=======================================================================
C MOMENT m2 = INTEGRALE ( f**2 * S(f) * df )
C=======================================================================
C
      CALL OS( 'X=CY    ', T2, T1  , SBID , FREQ2 )
      CALL OS( 'X=CX    ', T2, SBID , SBID , PONDER )
      CALL OS( 'X=X+Y   ', T02, T2 , SBID , 1.D0 )
C
C=======================================================================
C MOMENT mt1 = INTEGRALE ( T * S(f) * df ) 
C=======================================================================
C
      CALL OS( 'X=CY    ', T2 , T1  , SBID , PER )
      CALL OS( 'X=CX    ', T2 , SBID , SBID , PONDER )
      CALL OS( 'X=X+Y   ', TM , T2  , SBID , BID )
C
C=======================================================================
C MOMENT MCOS = INTEGRALE ( COS(INCI) * S(f) * df ) 
C=======================================================================
C
      CALL OS( 'X=COS(Y)',  T2 , INCI , SBID , BID )
      CALL OS( 'X=CXY   ',  T2 , T1   , SBID , PONDER )
      CALL OS( 'X=X+Y   ', MCOS, T2   , SBID , BID )
C 
C=======================================================================
C MOMENT MSIN = INTEGRALE ( SIN(INCI) * S(f) * df ) 
C=======================================================================
C
      CALL OS( 'X=SIN(Y)',  T2 , INCI , SBID , BID )
      CALL OS( 'X=CXY   ',  T2 , T1   , SBID , PONDER )
      CALL OS( 'X=X+Y   ', MSIN, T2   , SBID , BID )
C 
      RETURN
      END
