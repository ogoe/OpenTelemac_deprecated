C                       ***************
                        SUBROUTINE LUMP
C                       ***************
C
     *(DIAG,A,MESH,XMUL)
C
C***********************************************************************
C BIEF VERSION 5.5           08/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C     FONCTION  : SOMME DES TERMES PAR LIGNE D'UNE MATRICE A
C                 LE RESULTAT EST MULTIPLIE PAR XMUL.
C
C                 ON FAIT SIMPLEMENT DIAG = A X (VECTEUR EGAL A XMUL)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |    DIAG        |<-- | VECTEUR RESULTAT.
C |    A           | -->| MATRICE
C |    MESH        | -->| BLOC DES TABLEAUX ENTIERS DU MAILLAGE.
C |    XMUL        | -->| COEFFICIENT MULTIPLICATEUR
C |    MSK         | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |    MASKEL      | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
C SOUS-PROGRAMMES APPELES : ELMMAT , CPSTVC
C
C***********************************************************************
C
      USE BIEF, EX_LUMP => LUMP
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DOUBLE PRECISION, INTENT(IN)   :: XMUL
      TYPE(BIEF_OBJ), INTENT(IN)     :: A
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: DIAG
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DOUBLE PRECISION C
C
C-----------------------------------------------------------------------
C
      IF(A%ELMLIN.NE.A%ELMCOL) THEN
        IF (LNG.EQ.1) WRITE(LU,50)
        IF (LNG.EQ.2) WRITE(LU,51)
50      FORMAT(1X,'LUMP (BIEF) : A N''EST PAS UNE MATRICE CARREE')
51      FORMAT(1X,'LUMP (BIEF) : A IS NOT A SQUARE MATRIX')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      CALL CPSTVC(A%D,DIAG)
C
C-----------------------------------------------------------------------
C
C  CONSTRUCTION D'UN VECTEUR QUI VAUT PARTOUT XMUL
C
      CALL OS( 'X=C     ', X=DIAG , C=XMUL )
C
C  DIAG EST LE PRODUIT DE A PAR CE VECTEUR
C  DIAG JOUE ICI LE ROLE DE X ET Y (ON PEUT LE FAIRE).
C
      CALL MATVEC('X=AY    ',DIAG,A,DIAG,C,MESH,.TRUE.)
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
