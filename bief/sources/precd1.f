C                       *****************
                        SUBROUTINE PRECD1
C                       *****************
C
     *(X,A,B,D,MESH,PRECON,PREXSM,DIADON)
C
C***********************************************************************
C BIEF VERSION 6.0     06/07/2009    J-M HERVOUET (LNHE)  01 30 87 80 18
C***********************************************************************
C
C    FONCTION:  PRECONDITIONNEMENT DIAGONAL D'UN SYSTEME A X = B
C               (VOIR EXPLICATIONS DANS PRECDT).
C
C    ICI A EST UNE MATRICE SIMPLE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   X            |<-->|  VALEURS A L' ETAPE N+1.
C |   A            | -->|  MATRICE
C |   B            | -->|  SECONDS MEMBRES DU SYSTEME.
C |   D            |<-- |  STOCKAGE DE MATRICES DIAGONALES
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
      USE BIEF, EX_PRECD1 => PRECD1
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
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X,B,D
C
C-----------------------------------------------------------------------
C
C  STRUCTURE DE MATRICE
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
      DOUBLE PRECISION C
C
C-----------------------------------------------------------------------
C
C  PREPARATION DES DIAGONALES :
C
      IF(.NOT.DIADON) THEN
C
C CALCUL DES RACINES CARREES DES VALEURS ABSOLUES OU DES VALEURS
C
        IF(PRECON.EQ.5) THEN
          CALL OS( 'X=ABS(Y)' , X=D , Y=A%D )
        ELSE
          CALL OS( 'X=Y     ' , X=D , Y=A%D )
        ENDIF
C
C PARALLELISME : DIAGONALE COMPLETE AVANT DE FAIRE LA RACINE
C
        IF(NCSIZE.GT.1) THEN
          CALL PARCOM(D,2,MESH)
        ENDIF
C
        CALL OS( 'X=SQR(Y)' , X=D , Y=D )
C
C-----------------------------------------------------------------------
C                                                 -1
C  CHANGEMENT DE VARIABLE (DANS D IL Y A EN FAIT D
C
        IF(PREXSM) CALL OS( 'X=XY    ' , X , D , D , C )
C
C-----------------------------------------------------------------------
C
C CALCUL DES INVERSES DES RACINES CARREES DES DIAGONALES
C ON AURA AINSI LA VRAIE D ET NON PLUS SON INVERSE
C
        CALL OS( 'X=1/Y   ' , D , D , D , C , 2 , 1.D0 , 1.D-10 )
C
      ELSE
C
C  CAS OU D EST DONNE, CHANGEMENT DE VARIABLES
C  CHANGEMENT DE VARIABLE (DANS D IL Y A VRAIMENT D )
C
        IF(PREXSM) THEN
          CALL OS( 'X=Y/Z   ' , X=X , Y=X , Z=D )
        ENDIF
C
      ENDIF
C
C=======================================================================
C PRECONDITIONNEMENT DE A :
C=======================================================================
C
      CALL OM( 'M=DMD   ' , A , A , D , C , MESH )
C     SI PRECON = 2 OU 3
      IF((2*(PRECON/2).EQ.PRECON.OR.3*(PRECON/3).EQ.PRECON).AND.
     *                                                 .NOT.DIADON) THEN
C       VALABLE SEULEMENT AVEC UN SEUL DOMAINE
        IF(NCSIZE.LE.1.OR.NPTIR.EQ.0) A%TYPDIA='I'
      ENDIF
C
C=======================================================================
C
C PRECONDITIONNEMENT DU SECOND MEMBRE
C
      IF(PREXSM) CALL OS( 'X=XY    ' , X=B , Y=D )
C
C=======================================================================
C
      RETURN
      END
