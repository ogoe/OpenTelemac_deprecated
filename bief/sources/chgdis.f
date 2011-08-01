C                       *****************
                        SUBROUTINE CHGDIS
C                       *****************
C
     *(X,OLDELT,NEWELT,MESH)
C
C***********************************************************************
C BIEF VERSION 5.9        13/02/08    J-M HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : CHANGEMENT DE DISCRETISATION POUR UN VECTEUR
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      X         |<-- | VECTEUR A MODIFIER
C |      OLDELT    | -->| ANCIEN TYPE
C |      NEWELT    | -->| NOUVEAU TYPE
C |      MESH      | -->| STRUCTURE DE MAILLAGE
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR :
C
C  SOUS-PROGRAMME APPELE : CG1112
C
C**********************************************************************
C
      USE BIEF, EX_CHGDIS => CHGDIS
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X
      INTEGER, INTENT(IN)           :: OLDELT,NEWELT
      TYPE(BIEF_MESH), INTENT(IN)   :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      IF(NBPTS(NEWELT).GT.X%MAXDIM1) THEN
        IF(LNG.EQ.1) WRITE(LU,200) X%NAME
        IF(LNG.EQ.2) WRITE(LU,201) X%NAME
200     FORMAT(1X,'CHGDIS (BIEF) : EXTENSION IMPOSSIBLE POUR ',A6)
201     FORMAT(1X,'CHGDIS (BIEF) : EXTENSION IMPOSSIBLE FOR ',A6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      X%DIM1 = NBPTS(NEWELT)

      IF(OLDELT.EQ.11.AND.NEWELT.EQ.12) THEN
C
        CALL CG1112(X%R,X%DIM1,X%DIM2,
     *              MESH%IKLE%I,MESH%NELEM,MESH%NELMAX)
C
      ELSEIF(OLDELT.EQ.11.AND.NEWELT.EQ.13) THEN
C
        CALL CG1113(X%R,X%DIM1,X%DIM2,
     *              MESH%IKLE%I,MESH%NELEM,MESH%NELMAX)
C
      ELSEIF((OLDELT.EQ.12.OR.OLDELT.EQ.13).AND.NEWELT.EQ.11) THEN
C
C       ON NE FAIT RIEN
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,10) OLDELT,NEWELT
        IF(LNG.EQ.2) WRITE(LU,11) OLDELT,NEWELT
10      FORMAT(1X,'CHGDIS : CAS NON PREVU :'    ,I6,' ',I6)
11      FORMAT(1X,'CHGDIS: CASE NOT IMPLEMENTED:',I6,' ',I6)
        WRITE(LU,*) 'STRUCTURE X = ',X%NAME
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      X%ELM=NEWELT
      X%DIM1=NBPTS(NEWELT)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
