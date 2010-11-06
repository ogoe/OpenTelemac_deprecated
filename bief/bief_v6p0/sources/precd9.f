C                       *****************
                        SUBROUTINE PRECD9
C                       *****************
C
     *(X1,X2,X3,A11,A12,A13,A21,A22,A23,A31,A32,A33,
     * B1,B2,B3,D1,D2,D3,MESH,PRECON,PREXSM,DIADON)
C
C***********************************************************************
C BIEF VERSION 5.1           24/04/97    J-M HERVOUET (LNH)  30 87 80 18
C***********************************************************************
C
C    FONCTION:  PRECONDITIONNEMENT DIAGONAL D'UN SYSTEME A X = B
C               (VOIR EXPLICATIONS DANS PRECDT).
C
C    ICI A EST UN BLOC DE 9 MATRICES.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   X1,2,3       |<-->|  VALEURS A L' ETAPE N+1.
C |   A11,...      | -->|  MATRICES COMPOSANT LE BLOC A
C |   B1,2,3       | -->|  SECONDS MEMBRES DU SYSTEME.
C |   D1,2,3       |<-- |  STOCKAGE DE MATRICES DIAGONALES
C |   MESH         | -->|  MAILLAGE.
C |   PRECON       | -->|  VARIANTE DE PRECONDITIONNEMENT
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
      USE BIEF, EX_PRECD9 => PRECD9
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
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X1,X2,X3,B1,B2,B3,D1,D2,D3
C
C-----------------------------------------------------------------------
C
C  STRUCTURES DE MATRICE
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: A11,A12,A13,A21
      TYPE(BIEF_OBJ), INTENT(INOUT) :: A22,A23,A31,A32,A33
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
C  PROBLEME EN PARALLELISME, IL FAUDRAIT FAIRE ABS APRES PARCOM
          CALL OS( 'X=ABS(Y)' , D1 , A11%D , D1 , C )
          CALL OS( 'X=ABS(Y)' , D2 , A22%D , D2 , C )
          CALL OS( 'X=ABS(Y)' , D3 , A33%D , D3 , C )
        ELSE
          CALL OS( 'X=Y     ' , D1 , A11%D , D1 , C )
          CALL OS( 'X=Y     ' , D2 , A22%D , D2 , C )
          CALL OS( 'X=Y     ' , D3 , A33%D , D3 , C )
        ENDIF
C
C PARALLELISME : DIAGONALES COMPLETES AVANT DE FAIRE LA RACINE
C
        IF(NCSIZE.GT.1) THEN
          CALL PARCOM(D1,2,MESH)
          CALL PARCOM(D2,2,MESH)
          CALL PARCOM(D3,2,MESH)
        ENDIF
C
        CALL OS( 'X=SQR(Y)' , D1 , D1 , D1 , C )
        CALL OS( 'X=SQR(Y)' , D2 , D2 , D2 , C )
        CALL OS( 'X=SQR(Y)' , D3 , D3 , D3 , C )
C
C-----------------------------------------------------------------------
C                                                            -1
C  CHANGEMENT DE VARIABLE (DANS D1,D2 ET D3 IL Y A EN FAIT D ,...
C
        IF(PREXSM) THEN
          CALL OS( 'X=XY    ' , X1 , D1 , D1 , C )
          CALL OS( 'X=XY    ' , X2 , D2 , D2 , C )
          CALL OS( 'X=XY    ' , X3 , D3 , D3 , C )
        ENDIF
C
C-----------------------------------------------------------------------
C
C CALCUL DES INVERSES DES RACINES CARREES DES DIAGONALES
C ON AURA AINSI LES VRAIES D1,D2,D3 ET NON PLUS LEURS INVERSES
C
        CALL OS( 'X=1/Y   ' , D1 , D1 , D1 , C , 2 , 1.D0 , 1.D-10)
        CALL OS( 'X=1/Y   ' , D2 , D2 , D2 , C , 2 , 1.D0 , 1.D-10)
        CALL OS( 'X=1/Y   ' , D3 , D3 , D3 , C , 2 , 1.D0 , 1.D-10)
C
      ELSE
C
C  CAS OU D1,D2,D3 SONT DONNEES, CHANGEMENT DE VARIABLES
C  CHANGEMENT DE VARIABLE (DANS D1,D2,D3 IL Y A VRAIMENT D1,D2,D3 )
C
        IF(PREXSM) THEN
          CALL OS( 'X=Y/Z   ' , X1 , X1 , D1 , C )
          CALL OS( 'X=Y/Z   ' , X2 , X2 , D2 , C )
          CALL OS( 'X=Y/Z   ' , X3 , X3 , D3 , C )
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
C PRECONDITIONNEMENT DE A13 :
C=======================================================================
C
      CALL OM( 'M=DM    ' , A13,A13 ,D1,C,MESH)
      CALL OM( 'M=MD    ' , A13,A13 ,D3,C,MESH)
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
C PRECONDITIONNEMENT DE A23 :
C=======================================================================
C
      CALL OM( 'M=DM    ' , A23,A23 ,D2,C,MESH)
      CALL OM( 'M=MD    ' , A23,A23 ,D3,C,MESH)
C
C=======================================================================
C PRECONDITIONNEMENT DE A31 :
C=======================================================================
C
      CALL OM( 'M=DM    ' , A31,A31 ,D3,C,MESH)
      CALL OM( 'M=MD    ' , A31,A31 ,D1,C,MESH)
C
C=======================================================================
C PRECONDITIONNEMENT DE A32 :
C=======================================================================
C
      CALL OM( 'M=DM    ' , A32,A32 ,D3,C,MESH)
      CALL OM( 'M=MD    ' , A32,A32 ,D2,C,MESH)
C
C=======================================================================
C PRECONDITIONNEMENT DE A33 :
C=======================================================================
C
      CALL OM( 'M=DMD   ' , A33,A33 ,D3,C,MESH)
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
          A33%TYPDIA='I'
        ELSEIF(3*(PRECON/3).EQ.PRECON.AND..NOT.DIADON) THEN
          A11%TYPDIA='I'
          A22%TYPDIA='I'
          A33%TYPDIA='I'
          A12%TYPDIA='0'
          A13%TYPDIA='0'
          A21%TYPDIA='0'
          A23%TYPDIA='0'
          A31%TYPDIA='0'
          A32%TYPDIA='0'
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
        CALL OS( 'X=XY    ' , B3 , D3 , D3 , C )
      ENDIF
C
C=======================================================================
C
      RETURN
      END 
 
 
