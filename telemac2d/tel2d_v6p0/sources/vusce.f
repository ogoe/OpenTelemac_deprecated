C                       *******************************
                        DOUBLE PRECISION FUNCTION VUSCE
C                       *******************************
C
C
     *( TIME , I )
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C FONCTION  : DONNE LA VALEUR DE LA VITESSE DES SOURCES SUIVANT X
C
C             PERMET DE PROGRAMMER DES VITESSES VARIABLES EN FONCTION
C             DU TEMPS ET DE LA PROFONDEUR.
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C FUNCTION  : GIVES THE VALUE OF VELOCITY ALONG X AT SOURCES.
C
C             PERMET DE PROGRAMMER DES VITESSES VARIABLES EN FONCTION
C             DU TEMPS ET DE LA PROFONDEUR.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   TIME         | -->| TIME
C |   I            | -->| NUMBER OF THE SOURCE
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : BORD
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
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(IN) :: TIME
      INTEGER         , INTENT(IN) :: I
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     USCE STEMS FROM THE PARAMETER FILE, BUT VUSCE MAY BE MODIFIED HERE
C     (READ IN A FILE, ETC.)
C 
      VUSCE = USCE(I)
C
C-----------------------------------------------------------------------
C
      RETURN
      END       
