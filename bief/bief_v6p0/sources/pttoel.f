C                       *****************
                        SUBROUTINE PTTOEL
C                       *****************
C
     *(XEL,X,MESH)
C
C***********************************************************************
C BIEF VERSION 5.6           10/01/95    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C FONCTION : PASSAGE D'UN VECTEUR PAR POINTS A UN VECTEUR PAR
C            ELEMENTS.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      XEL       |<-- |  VECTEUR SUR LES ELEMENTS
C |      X         | -->|  VECTEUR PAR POINTS
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : PTEL11
C
C**********************************************************************
C
      USE BIEF, EX_PTTOEL => PTTOEL
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C  STRUCTURES DE VECTEURS
C
      TYPE(BIEF_OBJ), INTENT(IN)    :: X
      TYPE(BIEF_OBJ), INTENT(INOUT) :: XEL
C
C-----------------------------------------------------------------------
C
C  STRUCTURE DE MAILLAGE
C
      TYPE(BIEF_MESH), INTENT(IN) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      IF(X%ELM.EQ.11) THEN
C
        CALL PTEL11(XEL%R,X%R,MESH%IKLE%I,MESH%NELMAX,MESH%NELEM)
C
      ELSE
       IF (LNG.EQ.1) WRITE(LU,100) X%ELM
       IF (LNG.EQ.2) WRITE(LU,101) X%ELM
100    FORMAT(1X,'PTTOEL (BIEF) : IELM = ',1I6,' ELEMENT NON PREVU')
101    FORMAT(1X,'PTTOEL (BIEF) : IELM = ',1I6,' ELEMENT NOT AVAILABLE')
       CALL PLANTE(1)
       STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
