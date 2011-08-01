C                       *****************
                        SUBROUTINE SLOPES
C                       *****************
C
     *(COEF,Z,MESH)
C
C***********************************************************************
C  BIEF VERSION 5.5          17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C FONCTION : CALCUL DU COEFFICIENT 1 / COS(ALFA)
C
C            OU ALFA EST LA PENTE D'UN ELEMENT TRIANGULAIRE.
C
C            CE COEFFICIENT EST UTILISE POUR LES TERMES DE FROTTEMENT
C            SUR LE FOND.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      COEF      |<-- |  RESULTAT
C |      Z         | -->|  COTE DU FOND
C |      MESH      | -->|  MAILLAGE
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : PROPAG
C
C  SOUS-PROGRAMME APPELE : NEANT
C
C-----------------------------------------------------------------------
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_SLOPES => SLOPES
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  STRUCTURES DE VECTEUR
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: COEF
      TYPE(BIEF_OBJ), INTENT(IN)    :: Z
C
C-----------------------------------------------------------------------
C
C  STRUCTURE DE MAILLAGE
C
      TYPE(BIEF_MESH), INTENT(IN)   :: MESH
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER NELEM,NELMAX,IELM
C
C-----------------------------------------------------------------------
C
      IELM   = MESH%TYPELM
      NELEM  = MESH%NELEM
      NELMAX = MESH%NELMAX
C
C-----------------------------------------------------------------------
C
      COEF%ELM = IELM
      COEF%DIM1= NBPTS(IELM)
C
C-----------------------------------------------------------------------
C
      IF(IELM.EQ.10) THEN
        CALL SLOP10(COEF%R,MESH%XEL%R,MESH%YEL%R,
     *              Z%R,MESH%IKLE%I,NELEM,NELMAX)
      ELSE
        IF(LNG.EQ.1) WRITE(LU,300) MESH%TYPELM
        IF(LNG.EQ.2) WRITE(LU,301) MESH%TYPELM
300     FORMAT(1X,'SLOPES (BIEF) : ELEMENT NON PREVU : ',1I6)
301     FORMAT(1X,'SLOPES (BIEF) : UNKNOWN ELEMENT:',1I6)
        CALL PLANTE(0)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
