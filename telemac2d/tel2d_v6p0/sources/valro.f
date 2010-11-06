C                       ****************
                        SUBROUTINE VALRO
C                       ****************
C
     *(RO,S,ROEAU)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    01/09/94    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C      FONCTION: CALCUL DE LA MASSE VOLUMIQUE EN FONCTION DE LA
C                SALINITE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    RO          |<-- |  TABLEAU DE LA MASSE VOLUMIQUE.
C |    S           | -->|  BLOC DES TRACEURS.
C |    ROEAU       | -->|  MASSE VOLUMIQUE DE L'EAU A LA TEMPERATURE.
C |                |    |  MOYENNE, QUAND LA SALINITE EST NULLE.
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : PREDON
C
C  SOUS-PROGRAMME APPELE : OS
C
C**********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(IN)    :: S
      TYPE(BIEF_OBJ), INTENT(INOUT) :: RO
      DOUBLE PRECISION, INTENT(IN)  :: ROEAU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     BEWARE: HERE IT IS ASSUMED THAT SALINITY IS THE FIRST TRACER
C
      CALL OS( 'X=CY    ' , X=RO , Y=S%ADR(1)%P , C=0.749979D0 )
      CALL OS( 'X=X+C   ' , X=RO , C=ROEAU      )
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
