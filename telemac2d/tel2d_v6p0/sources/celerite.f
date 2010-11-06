C                       *******************
                        SUBROUTINE CELERITE
C                       *******************
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  :  CALCUL DU CHAMP DE CELERITE DE REFERENCE (C0)
C                  QUI SERT DANS LE CAS D'ONDE INCIDENTE.
C
C                  ICI C0 = RACINE(G H INITIAL)
C                  MAIS C0 EST UN TABLEAU DE BORD.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      C0        |<-- | CELERITE INITIALE                            |
C |      H         | -->| HAUTEUR INITIALE
C |      GRAV      | -->| PESANTEUR                                    |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS PROGRAMME APPELE : OV
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C-----------------------------------------------------------------------
C
      CALL OSBD( 'X=CY    ' , C0 , H  , H , GRAV , MESH )
      CALL OS  ( 'X=+(Y,C)' , C0 , C0 , H , 0.D0        )
      CALL OS  ( 'X=SQR(Y)' , C0 , C0 , H , 0.D0        )
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
