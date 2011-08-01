C                       *****************
                        SUBROUTINE PRECDT
C                       *****************
C
     *(X,A,B,D,MESH,PRECON,PREXSM,DIADON,S)
C
C***********************************************************************
C BIEF VERSION 5.1           24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C   FONCTION:  PRECONDITIONNEMENT DIAGONAL D'UN SYSTEME A X = B
C              QUI PEUT ETRE COMPOSE D'UNE MATRICE OU DE
C              BLOCS DE 4 OU 9 MATRICES.
C
C   SI DIADON=.TRUE. LES DIAGONALES DE PRECONDITIONNEMENT
C   SONT DONNEES PAR L'UTILISATEUR.
C
C   SI DIADON=.FALSE. :
C
C   SI PRECON EST DIVISIBLE PAR 2 : PRECONDITIONNEMENT DIAGONAL AVEC
C                                   LA DIAGONALE DE A.
C   SI PRECON EST DIVISIBLE PAR 3 : PRECONDITIONNEMENT BLOC-DIAGONAL
C                                   CE PRECONDITIONNEMENT EST COMMENCE
C                                   PAR LE SOUS-PROGRAMME PREBDT, MAIS
C                                   SE TERMINE PAR UN PRECONDITIONNEMENT
C                                   DIAGONAL CLASSIQUE QUE L'ON FAIT ICI
C   SI PRECON EST DIVISIBLE PAR 5 : PRECONDITIONNEMENT DIAGONAL AVEC
C                                   LA VALEUR ABSOLUE DE LA DIAGONALE DE
C                                   A.
C
C-----------------------------------------------------------------------
C
C   EXEMPLE D'UN BLOC DE 4 :
C
C        (  A11    A12  ) ( X1 ) = ( B1 )
C        (              ) (    )   (    )
C        (              ) (    )   (    )
C        (  A21    A22  ) ( X2 ) = ( B2 )
C
C   LES MATRICES DIAGONALES DE PRECONDITIONNEMENT SONT D1 ET D2
C
C   CES DIAGONALES SONT DONNEES PAR L'UTILISATEUR: DIADON=.TRUE.)
C   OU BIEN CE SONT LES DIAGONALES DE A11,A22 : PRECON DIVISIBLE PAR 2
C
C                                                       -1
C    ON EFFECTUE LE CHANGEMENT DE VARIABLE X1PRIME =  D1    X1
C                                                       -1
C                                          X2PRIME =  D2    X2
C
C    PUIS TOUT LE SYSTEME EST MULTIPLIE A GAUCHE PAR D
C
C    ON EFFECTUE AINSI LE PRODUIT:
C
C  ( D1   0  )       (  A11   A12  )     ( D1   0  )
C  (         )       (             )     (         )
C  (         )   X   (             )     (         )
C  ( 0   D2  )       (  A21   A22  )  X  ( 0   D2  )
C
C
C   CE QUI DONNE :
C
C           (  D1  A11 D1       D1  A12 D2  )
C           (                               )
C           (                               )
C           (  D2  A21 D1       D2  A22 D2  )
C
C
C   ATTENTION : TOUTES LES MATRICES AIJ SONT REMPLACEES
C               PAR LEUR VALEUR APRES PRECONDITIONNEMENT.
C
C   TRAITEMENT DES SECONDS MEMBRES :
C
C   ( D1   0 )       ( B1 )
C   (        )  X    (    )
C   ( 0   D2 )       ( B2 )
C
C   NOTE : NE PAS OUBLIER APRES RESOLUTION D'INVERSER LE CHANGEMENT
C          DE VARIABLE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   X            |<-->|  VALEURS A L' ETAPE N+1.
C |   A            | -->| MATRICE OU BLOC DE MATRICES
C |   B            | -->|  SECONDS MEMBRES DU SYSTEME.
C |   D            |<-- |  STOCKAGE DE MATRICES DIAGONALES
C |                |    |  (TOUJOURS DANS UN BLOC)
C |   MESH         | -->|  MAILLAGE.
C |   PRECON       | -->|  VARIANTE DE PRECONDITIONNEMENT
C |                |    |  ICI PRECON DOIT ETRE MULTIPLE DE 2,3 OU 5
C |   PREXSM       | -->|  .TRUE. : ON PRECONDITIONNE X1,X2,X3 ET SM
C |   DIADON       | -->|  .TRUE. : LES DIAGONALES SONT DONNEES.
C |   S            | -->|  0 : SYSTEME NORMAL
C |                |    |  1 : BLOC A UNE MATRICE
C |                |    |  2 : BLOC A 4   MATRICES
C |                |    |  3 : BLOC A 9   MATRICES
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : PRECD1, PRECD4, PRECD9.
C
C**********************************************************************
C
      USE BIEF, EX_PRECDT => PRECDT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: PRECON,S
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
C  CAS CLASSIQUE
C
      IF(S.EQ.0) THEN
C
        CALL PRECD1(X,A,B,D%ADR(1)%P,MESH,PRECON,PREXSM,DIADON)
C
C-----------------------------------------------------------------------
C
C  CAS OU A X,B ET D SONT DES BLOCS DE 1 OBJET
C
      ELSEIF(S.EQ.1) THEN
C
        CALL PRECD1(X%ADR(1)%P,A%ADR(1)%P,B%ADR(1)%P,D%ADR(1)%P,MESH,
     *              PRECON,PREXSM,DIADON)
C
      ELSEIF(S.EQ.2) THEN
C
        CALL PRECD4(X%ADR(1)%P,X%ADR(2)%P,
     *              A%ADR(1)%P,A%ADR(2)%P,A%ADR(3)%P,A%ADR(4)%P,
     *              B%ADR(1)%P,B%ADR(2)%P,
     *              D%ADR(1)%P,D%ADR(2)%P,
     *              MESH,PRECON,PREXSM,DIADON)
C
      ELSEIF(S.EQ.3) THEN
C
        CALL PRECD9(X%ADR(1)%P,X%ADR(2)%P,X%ADR(3)%P,
     *              A%ADR(1)%P,A%ADR(2)%P,A%ADR(3)%P,A%ADR(4)%P,
     *              A%ADR(5)%P,A%ADR(6)%P,A%ADR(7)%P,A%ADR(8)%P,
     *              A%ADR(9)%P,
     *              B%ADR(1)%P,B%ADR(2)%P,B%ADR(3)%P,
     *              D%ADR(1)%P,D%ADR(2)%P,D%ADR(3)%P,
     *              MESH,PRECON,PREXSM,DIADON)
C
      ELSE
C
C-----------------------------------------------------------------------
C
C  ERREUR
C
        IF(LNG.EQ.1) WRITE(LU,100) S
        IF(LNG.EQ.2) WRITE(LU,200) S
100     FORMAT(1X,'PRECDT (BIEF) : S NON PREVU :',1I6)
200     FORMAT(1X,'PRECDT (BIEF) : UNEXPECTED S :',1I6)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
