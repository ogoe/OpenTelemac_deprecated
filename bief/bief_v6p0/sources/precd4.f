C                       *****************
                        SUBROUTINE PRECD4
C                       *****************
C
     *(X1,X2,A11,A12,A21,A22,
     * B1,B2,D1,D2,MESH,PRECON,PREXSM,DIADON)
C
C***********************************************************************
C BIEF VERSION 6.0     06/07/2009    J-M HERVOUET (LNHE)  01 30 87 80 18
C***********************************************************************
C
C    FONCTION:  PRECONDITIONNEMENT DIAGONAL D'UN SYSTEME A X = B
C               (VOIR EXPLICATIONS DANS PRECDT).
C
C    ICI A EST UN BLOC DE 4 MATRICES.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   X1,2         |<-->|  VALEURS A L' ETAPE N+1.
C |   A11,12,21,22 | -->|  MATRICES COMPOSANT LE BLOC A
C |   B1,2         | -->|  SECONDS MEMBRES DU SYSTEME.
C |   D1,2         |<-- |  STOCKAGE DE MATRICES DIAGONALES
C |   MESH         | -->|  MAILLAGE.
C |   PRECON       | -->|  VARIANTE DE PRECONDITIONNEMENT
C |   PREXSM       | -->|  .TRUE. : ON PRECONDITIONNE X1,X2 ET SM
C |   DIADON       | -->|  .TRUE. : LES DIAGONALES SONT DONNEES.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : OS , OM
C
C**********************************************************************
C
      USE BIEF, EX_PRECD4 => PRECD4
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: PRECON
C
      LOGICAL, INTENT(IN) :: PREXSM,DIADON
C
C-----------------------------------------------------------------------
C
C  STRUCTURES DE VECTEURS
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X1,X2,B1,B2,D1,D2
C
C-----------------------------------------------------------------------
C
C  STRUCTURES DE MATRICE
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: A11,A12,A21,A22
C
C-----------------------------------------------------------------------
C
C  STRUCTURE DE MAILLAGE
C
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION C
C
C-----------------------------------------------------------------------
C
C  PREPARATION DES DIAGONALES :
C
      IF(.NOT.DIADON) THEN
C
C CALCUL DES RACINES CARREES DES VALEURS ABSOLUES
C
        IF(PRECON.EQ.5) THEN
          CALL OS( 'X=ABS(Y)' , D1 , A11%D , D1 , C )
          CALL OS( 'X=ABS(Y)' , D2 , A22%D , D2 , C )
        ELSE
          CALL OS( 'X=Y     ' , D1 , A11%D , D1 , C )
          CALL OS( 'X=Y     ' , D2 , A22%D , D2 , C )
        ENDIF
C
C PARALLELISME : DIAGONALE COMPLETE AVANT DE FAIRE LA RACINE
C
        IF(NCSIZE.GT.1) THEN
          CALL PARCOM(D1,2,MESH)
          CALL PARCOM(D2,2,MESH)
        ENDIF
C
        CALL OS( 'X=SQR(Y)' , D1 , D1 , D1 , C )
        CALL OS( 'X=SQR(Y)' , D2 , D2 , D2 , C )
C
C-----------------------------------------------------------------------
C                                                            -1
C  CHANGEMENT DE VARIABLE (DANS D,D2 ET D3 IL Y A EN FAIT D ,...
C
        IF(PREXSM) THEN
          CALL OS( 'X=XY    ' , X1 , D1 , D1 , C )
          CALL OS( 'X=XY    ' , X2 , D2 , D2 , C )
        ENDIF
C
C-----------------------------------------------------------------------
C
C CALCUL DES INVERSES DES RACINES CARREES DES DIAGONALES
C ON AURA AINSI LES VRAIES D1,D2 ET NON PLUS LEURS INVERSES
C
        CALL OS( 'X=1/Y   ' , D1 , D1 , D1 , C , 2 , 1.D0 , 1.D-10)
        CALL OS( 'X=1/Y   ' , D2 , D2 , D2 , C , 2 , 1.D0 , 1.D-10)
C
      ELSE
C
C  CAS OU D EST DONNE, CHANGEMENT DE VARIABLES
C  CHANGEMENT DE VARIABLE (DANS D1,D2   IL Y A VRAIMENT D1,D2)
C
        IF(PREXSM) THEN
          CALL OS( 'X=Y/Z   ' , X1 , X1 , D1 , C )
          CALL OS( 'X=Y/Z   ' , X2 , X2 , D2 , C )
        ENDIF
C
      ENDIF
C
C=======================================================================
C PRECONDITIONNEMENT DE A11 :
C=======================================================================
C
      CALL OM( 'M=DMD   ' , A11,A11 ,D1,C,MESH)
C
C=======================================================================
C PRECONDITIONNEMENT DE A12 :
C=======================================================================
C
      CALL OM( 'M=DM    ' , A12,A12 ,D1,C,MESH)
      CALL OM( 'M=MD    ' , A12,A12 ,D2,C,MESH)
C
C=======================================================================
C PRECONDITIONNEMENT DE A21 :
C=======================================================================
C
      CALL OM( 'M=DM    ' , A21,A21 ,D2,C,MESH)
      CALL OM( 'M=MD    ' , A21,A21 ,D1,C,MESH)
C
C=======================================================================
C PRECONDITIONNEMENT DE A22 :
C=======================================================================
C
      CALL OM( 'M=DMD   ' , A22,A22 ,D2,C,MESH)
C
C=======================================================================
C
C     DIFFERENTS CAS OU LES DIAGONALES SONT CONNUES
C     (VALABLE SEULEMENT AVEC UN SEUL DOMAINE)
C
      IF(NCSIZE.LE.1.OR.NPTIR.EQ.0) THEN
C
C       SI PRECON = 2 OU 3
        IF(2*(PRECON/2).EQ.PRECON.AND..NOT.DIADON) THEN
          A11%TYPDIA='I'
          A22%TYPDIA='I'
        ELSEIF(3*(PRECON/3).EQ.PRECON.AND..NOT.DIADON) THEN
          A11%TYPDIA='I'
          A22%TYPDIA='I'
          A12%TYPDIA='0'
          A21%TYPDIA='0'
        ENDIF
C
      ENDIF
C
C=======================================================================
C
C PRECONDITIONNEMENT DU SECOND MEMBRE
C
      IF(PREXSM) THEN
        CALL OS( 'X=XY    ' , B1 , D1 , D1 , C )
        CALL OS( 'X=XY    ' , B2 , D2 , D2 , C )
      ENDIF
C
C=======================================================================
C
      RETURN
      END 
 
 
