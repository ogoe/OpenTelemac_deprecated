C                       ***************
                        SUBROUTINE UM1X
C                       ***************
C
     *(X,D,S)
C
C***********************************************************************
C BIEF VERSION 5.1           23/12/94  J.M. HERVOUET (LNH)  30 87 80 18
C***********************************************************************
C
C    FONCTION : FIN DES OPERATIONS DU PRECONDITIONNEMENT BLOC-DIAGONAL
C               OU TOUTE OPERATION SEMBLABLE
C
C                              -1   PRIME
C               OPERATION X = U    X
C
C    EXEMPLE D'UN BLOC DE 4 :
C
C              (   I     D12  )
C              (              )
C         U =  (              )
C              (              )
C              (   0      I   )
C                                   PRIME          PRIME         PRIME
C              (   I    -D12  ) ( X1      )    ( X1    -   D12 X2 )
C  -1  PRIME   (              ) (         )    (                  )
C U   X     =  (              ) (         )  = (                  )
C              (              ) (   PRIME )    (                  )
C              (   0      I   ) ( X2      )    ( X2               )
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   X            |<-->| X ET X' (TRANSFORMATION SUR PLACE)
C |   D            |<-- |  BLOC DE MATRICES DIAGONALES
C |   S            | -->|  2 : BLOC A 4   MATRICES
C |                |    |  3 : BLOC A 9   MATRICES
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : OS , OM
C
C**********************************************************************
C
      USE BIEF, EX_UM1X => UM1X
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
C     STRUCTURES DE VECTEURS OU DE BLOCS DE VECTEURS
C
      INTEGER, INTENT(IN)           :: S
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X
      TYPE(BIEF_OBJ), INTENT(IN)    :: D
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
C     BLOCS DE 4 MATRICES :
C
      IF(S.EQ.2) THEN
C
C     BLOCS DE 4 MATRICES :
C
        CALL UM1X04(X%ADR(1)%P,X%ADR(2)%P,D%ADR(3)%P)
C
      ELSEIF(S.EQ.3) THEN
C
C     BLOCS DE 9 MATRICES :
C
        CALL UM1X09(X%ADR(1)%P,X%ADR(2)%P,X%ADR(3)%P,
     *              D%ADR(4)%P,D%ADR(5)%P,D%ADR(7)%P)
C
      ELSE
C
C-----------------------------------------------------------------------
C
C  ERREUR
C
        IF(LNG.EQ.1) WRITE(LU,100) S
        IF(LNG.EQ.2) WRITE(LU,200) S
100     FORMAT(1X,'UM1X (BIEF) : S NON PREVU :',1I6)
200     FORMAT(1X,'UM1X (BIEF) : UNEXPECTED S :',1I6)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C                                                            -1
      RETURN
      END  
 
 
