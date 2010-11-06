C                       *****************
                        SUBROUTINE PREBDT
C                       *****************
C
     *(X,A,B,D,MESH,PREXSM,DIADON,S)
C
C***********************************************************************
C BIEF VERSION 5.1           23/12/94  J.M. HERVOUET (LNH)  30 87 80 18
C***********************************************************************
C
C    FONCTION:  PRECONDITIONNEMENT DIAGONAL-BLOC D'UN SYSTEME A X = B
C               QUI PEUT ETRE COMPOSE DE BLOCS DE 4 OU 9 MATRICES.
C
C    EXEMPLE D'UN BLOC DE 4 :
C
C         (   A11   A12  ) ( X1 ) = ( B1 )
C         (              ) (    )   (    )
C         (              ) (    )   (    )
C         (   A21   A22  ) ( X2 ) = ( B2 )
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   X            |<-->|  BLOC DES VECTEURS INCONNUS
C |   A            | -->|  BLOC DE MATRICES
C |   B            | -->|  BLOC DES SECONDS MEMBRES DU SYSTEME.
C |   D            |<-- |  BLOC DE DIAGONALES
C |   MESH         | -->|  BLOC DES TABLEAUX D'ENTIERS DU MAILLAGE.
C |   PREXSM       | -->|  .TRUE. : ON PRECONDITIONNE X,X2,X3 ET SM
C |   DIADON       | -->|  .TRUE. : LES DIAGONALES SONT DONNEES.
C |   S            | -->|  0 : SYSTEME NORMAL       (INTERDIT ICI)
C |                |    |  1 : BLOC A UNE MATRICE   (INTERDIT ICI)
C |                |    |  2 : BLOC A 4   MATRICES
C |                |    |  3 : BLOC A 9   MATRICES
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : OS , OM
C
C**********************************************************************
C
      USE BIEF, EX_PREBDT => PREBDT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: S
C
      LOGICAL, INTENT(IN) :: PREXSM,DIADON
C
C-----------------------------------------------------------------------
C
C  STRUCTURES DE VECTEURS OU DE BLOCS DE VECTEURS
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X,B,D
C
C-----------------------------------------------------------------------
C
C  STRUCTURE DE MATRICE OU DE BLOC DE MATRICES
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: A
C
C-----------------------------------------------------------------------
C
C  STRUCTURE DE MAILLAGE
C
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C  BLOCS DE 4 MATRICES :
C
      IF(S.EQ.2) THEN
C
C  BLOCS DE 4 MATRICES :
C
        CALL PREBD4(X%ADR(1)%P,X%ADR(2)%P,
     *              A%ADR(1)%P,A%ADR(2)%P,A%ADR(3)%P,A%ADR(4)%P,
     *              B%ADR(1)%P,B%ADR(2)%P,
     *              D%ADR(1)%P,D%ADR(2)%P,D%ADR(3)%P,D%ADR(4)%P,
     *              MESH,PREXSM,DIADON)
C
      ELSEIF(S.EQ.3) THEN
C
C  BLOCS DE 9 MATRICES :
C
        CALL PREBD9(X%ADR(1)%P,X%ADR(2)%P,X%ADR(3)%P,
     *              A%ADR(1)%P,A%ADR(2)%P,A%ADR(3)%P,
     *              A%ADR(4)%P,A%ADR(5)%P,A%ADR(6)%P,
     *              A%ADR(7)%P,A%ADR(8)%P,A%ADR(9)%P,
     *              B%ADR(1)%P,B%ADR(2)%P,B%ADR(3)%P,
     *              D%ADR(1)%P,D%ADR(2)%P,D%ADR(3)%P,
     *              D%ADR(4)%P,D%ADR(5)%P,D%ADR(6)%P,
     *              D%ADR(7)%P,D%ADR(8)%P,D%ADR(9)%P,
     *              MESH,PREXSM,DIADON)
C
      ELSE
C
C-----------------------------------------------------------------------
C
C  ERREUR
C
        IF(LNG.EQ.1) WRITE(LU,100) S
        IF(LNG.EQ.2) WRITE(LU,200) S
100     FORMAT(1X,'PREBDT (BIEF) : S NON PREVU :',1I6)
200     FORMAT(1X,'PREBDT (BIEF) : UNEXPECTED S :',1I6)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C                                                            -1
      RETURN
      END 
 
 
