C                       *****************
                        SUBROUTINE PREBD9
C                       *****************
C
     *(X1,X2,X3,A11,A12,A13,A21,A22,A23,A31,A32,A33,
     * B1,B2,B3,D11,D12,D13,D21,D22,D23,D31,D32,D33,
     * MESH,PREXSM,DIADON)
C
C***********************************************************************
C BIEF VERSION 5.1           23/12/94  J.M. HERVOUET (LNH)  30 87 80 18
C***********************************************************************
C
C    FONCTION:  PRECONDITIONNEMENT BLOC-DIAGONAL D'UN SYSTEME A X = B
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   X1,2,3       |<-->|  VALEURS A L' ETAPE N+1.
C |   A11,...      | -->|  MATRICES DU BLOC DE 9
C |   B1,2,3       | -->|  SECONDS MEMBRES DU SYSTEME.
C |   D11,...      |<-- |  STOCKAGE DE MATRICES DIAGONALES
C |   MESH         | -->|  BLOC DES ENTIERS DU MAILLAGE.
C |   PREXSM       | -->|  .TRUE. : ON PRECONDITIONNE X,X2,X3 ET SM
C |   DIADON       | -->|  .TRUE. : LES DIAGONALES SONT DONNEES.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : OS , OM
C
C**********************************************************************
C
      USE BIEF, EX_PREBD9 => PREBD9
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+     
C
      LOGICAL, INTENT(IN) :: PREXSM,DIADON
C
C-----------------------------------------------------------------------
C
C  STRUCTURES DE VECTEURS
C
      TYPE(BIEF_OBJ), INTENT(IN)    :: X3,B1
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X1,X2,B2,B3
      TYPE(BIEF_OBJ), INTENT(INOUT) :: D11,D12,D13,D21,D22,D23
      TYPE(BIEF_OBJ), INTENT(INOUT) :: D31,D32,D33
C
C-----------------------------------------------------------------------
C
C  STRUCTURES DE MATRICE
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: A11,A12,A13,A21,A22
      TYPE(BIEF_OBJ), INTENT(INOUT) :: A23,A31,A32,A33
C
C-----------------------------------------------------------------------
C
C  STRUCTURE DE MAILLAGE
C
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+     
C
      INTEGER I,NPOIN1,NPOIN2,NPOIN3
      DOUBLE PRECISION C      
C
C-----------------------------------------------------------------------
C
      NPOIN1 = X1%DIM1
      NPOIN2 = X2%DIM1
      NPOIN3 = X3%DIM1
C
      IF(NPOIN2.NE.NPOIN1.AND.NPOIN3.NE.NPOIN1) THEN
        IF(LNG.EQ.1) WRITE(LU,100)
        IF(LNG.EQ.2) WRITE(LU,200)
100     FORMAT(1X,'PREBD9 (BIEF) : MATRICES RECTANGULAIRES',/,1X,
     *  'PRECONDITIONNEMENT BLOC-DIAGONAL IMPOSSIBLE DANS CE CAS')
200     FORMAT(1X,'PREBD9 (BIEF) : RECTANGULAR MATRICES',/,1X,
     *  'BLOCK-DIAGONAL PRECONDITIONING IMPOSSIBLE IN THIS CASE')
        CALL PLANTE(0)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C  PREPARATION DES DIAGONALES :
C
      IF(.NOT.DIADON) THEN
C
        CALL OS( 'X=Y     ' , D11 , A11%D ,D11,C)
        CALL OS( 'X=Y     ' , D12 , A12%D ,D12,C)
        CALL OS( 'X=Y     ' , D13 , A13%D ,D13,C)
        CALL OS( 'X=Y     ' , D21 , A21%D ,D21,C)
        CALL OS( 'X=Y     ' , D22 , A22%D ,D22,C)
        CALL OS( 'X=Y     ' , D23 , A23%D ,D23,C)
        CALL OS( 'X=Y     ' , D31 , A31%D ,D31,C)
        CALL OS( 'X=Y     ' , D32 , A32%D ,D32,C)
        CALL OS( 'X=Y     ' , D33 , A33%D ,D33,C)
C
C  TEST POUR SE RAMENER AU PRECONDITIONNEMENT DIAGONAL
C
C       CALL OS( 'X=Y     ' , D11 , A11%D ,Z,C   )
C       CALL OS( 'X=C     ' , D12 , A12%D ,Z,0.D0)
C       CALL OS( 'X=C     ' , D13 , A13%D ,Z,0.D0)
C       CALL OS( 'X=C     ' , D21 , A21%D ,Z,0.D0)
C       CALL OS( 'X=Y     ' , D22 , A22%D ,Z,C   )
C       CALL OS( 'X=C     ' , D23 , A23%D ,Z,0.D0)
C       CALL OS( 'X=C     ' , D31 , A31%D ,Z,0.D0)
C       CALL OS( 'X=C     ' , D32 , A32%D ,Z,0.D0)
C       CALL OS( 'X=Y     ' , D33 , A33%D ,Z,C   )
C
C  FIN DU TEST POUR SE RAMENER AU PRECONDITIONNEMENT DIAGONAL
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  DECOMPOSITION L D U DU BLOC DES DIAGONALES :
C
C     D11 NE SERT DESORMAIS QUE PAR SON INVERSE
      CALL OS( 'X=1/Y   ' , D11 , D11 , D11 , C )
C
      DO 10 I = 1,NPOIN1
C
        D21%R(I) =  D21%R(I) * D11%R(I)
        D31%R(I) =  D31%R(I) * D11%R(I)
        D22%R(I) =  D22%R(I) - D21%R(I) * D12%R(I)
C
10    CONTINUE
C
C     D22 NE SERT QUE PAR SON INVERSE
      CALL OS( 'X=1/Y   ' , D22 , D22 , D22 , C )
C
      DO 11 I = 1,NPOIN1
C
        D32%R(I) = (D32%R(I) - D31%R(I) * D12%R(I)) * D22%R(I)
        D23%R(I) =  D23%R(I) - D21%R(I) * D13%R(I)
        D33%R(I) =  D33%R(I)
     *             -D31%R(I)*D13%R(I)-D32%R(I)*D23%R(I)
        D12%R(I) =  D12%R(I) * D11%R(I)
        D13%R(I) =  D13%R(I) * D11%R(I)
        D23%R(I) =  D23%R(I) * D22%R(I)
C
11    CONTINUE
C
C-----------------------------------------------------------------------
C
C CHANGEMENT DE VARIABLES :
C
      IF(PREXSM) THEN
C
        CALL OS( 'X=X+YZ  ' , X1 , X2 , D12 , C )
        CALL OS( 'X=X+YZ  ' , X1 , X3 , D13 , C )
        CALL OS( 'X=X+YZ  ' , X2 , X3 , D23 , C )
C
      ENDIF
C
C  ON CALCULE LA RACINE ET
C  ON INVERSE D11,D22,D33. DESORMAIS ILS NE SERVENT PLUS QUE SOUS
C  CETTE FORME
C
C     INVERSE DE D11 : DEJA FAIT PLUS HAUT.
C     INVERSE DE D22 : DEJA FAIT PLUS HAUT.
      CALL OS( 'X=1/Y   ' , D33 , D33 , D33 , C )
      CALL OS( 'X=SQR(Y)' , D11 , D11 , D11 , C )
      CALL OS( 'X=SQR(Y)' , D22 , D22 , D22 , C )
      CALL OS( 'X=SQR(Y)' , D33 , D33 , D33 , C )
C
C=======================================================================
C MULTIPLICATION DE A A GAUCHE PAR L'INVERSE DE L
C=======================================================================
C
C  A21 :
      CALL OM( 'M=M-DN  ' , A21 , A11 , D21 , C , MESH)
C  A22 :
      CALL OM( 'M=M-DN  ' , A22 , A12 , D21 , C , MESH)
C  A23 :
      CALL OM( 'M=M-DN  ' , A23 , A13 , D21 , C , MESH)
C  A31 :
      CALL OM( 'M=M-DN  ' , A31 , A11 , D31 , C , MESH)
      CALL OM( 'M=M-DN  ' , A31 , A21 , D32 , C , MESH)
C  A32 :
      CALL OM( 'M=M-DN  ' , A32 , A12 , D31 , C , MESH)
      CALL OM( 'M=M-DN  ' , A32 , A22 , D32 , C , MESH)
C  A33 :
      CALL OM( 'M=M-DN  ' , A33 , A13 , D31 , C , MESH)
      CALL OM( 'M=M-DN  ' , A33 , A23 , D32 , C , MESH)
C
C=======================================================================
C MULTIPLICATION DE A A DROITE PAR L'INVERSE DE U
C=======================================================================
C
C  A12 :
      CALL OM( 'M=M-ND  ' , A12 , A11 , D12 , C , MESH)
C  A22 :
      CALL OM( 'M=M-ND  ' , A22 , A21 , D12 , C , MESH)
C  A32 :
      CALL OM( 'M=M-ND  ' , A32 , A31 , D12 , C , MESH)
C  A13 :
      CALL OM( 'M=M-ND  ' , A13 , A11 , D13 , C , MESH)
      CALL OM( 'M=M-ND  ' , A13 , A12 , D23 , C , MESH)
C  A23 :
      CALL OM( 'M=M-ND  ' , A23 , A21 , D13 , C , MESH)
      CALL OM( 'M=M-ND  ' , A23 , A22 , D23 , C , MESH)
C  A33 :
      CALL OM( 'M=M-ND  ' , A33 , A31 , D13 , C , MESH)
      CALL OM( 'M=M-ND  ' , A33 , A32 , D23 , C , MESH)
C
C-----------------------------------------------------------------------
C
C NOUVEAU SECOND MEMBRE
C
      IF(PREXSM) THEN
C
      DO 21 I = 1,NPOIN1
        B2%R(I) = B2%R(I)-D21%R(I)*B1%R(I)
        B3%R(I) = B3%R(I)-D31%R(I)*B1%R(I)-D32%R(I)*B2%R(I)
21    CONTINUE
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
