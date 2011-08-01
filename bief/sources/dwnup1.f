C                       *****************
                        SUBROUTINE DWNUP1
C                       *****************
C
     *(X, A,B ,DITR,MESH)
C
C***********************************************************************
C BIEF VERSION 5.5         26/02/04    J-M HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : RESOLUTION DU SYSTEME A X = B
C
C            ICI LA MATRICE A EST LE RESULTAT D'UNE DECOMPOSITION
C            EFFECTUEE PAR LE SOUS-PROGRAMME DECLDU.
C
C            CHAQUE MATRICE ELEMENTAIRE A ETE DECOMPOSEE SOUS LA FORME :
C
C            LE X DE X UE
C
C            LE : TRIANGULAIRE INFERIEURE AVEC DES 1 SUR LA DIAGONALE.
C            DE : DIAGONALE
C            UE : TRIANGULAIRE SUPERIEURE AVEC DES 1 SUR LA DIAGONALE.
C
C                                                T
C            SI LA MATRICE EST SYMETRIQUE : LE =  UE
C
C            LES MATRICES "DE" SONT CONSIDEREES COMME DES DIAGONALES
C            DE TAILLE NPOIN X NPOIN QU'IL FAUT IMAGINER COMPLETEES
C            AVEC DES 1 POUR LES POINTS QUI N'APPARTIENNENT PAS A
C            L'ELEMENT CONSIDERE.
C
C            ON A EFFECTUE ENSUITE LE PRODUIT DE TOUTES CES DIAGONALES
C            CE QUI A DONNE LA DIAGONALE DB.
C
C !!!!!!!!!  ENFIN : DB A ETE INVERSEE CAR C'EST SOUS CETTE FORME
C                    QU'ELLE EST UTILISEE ICI.
C
C            LA MATRICE A EST ICI :
C
C            LE PRODUIT DE 1 A NELEM DE TOUTES LES MATRICES LE
C
C            MULTIPLIE PAR :
C
C            LA DIAGONALE DB
C
C            MULTIPLIE PAR :
C
C            LE PRODUIT DE NELEM A 1 DE TOUTES LES MATRICES UE.
C
C-----------------------------------------------------------------------
C  SIGNIFICATION DE IELM :
C
C  TYPE D'ELEMENT      NOMBRE DE POINTS          PROGRAMME ICI
C
C  11 : TRIANGLE P1            3                       OUI
C  12 : TRIANGLE QUASI-BULLE   4                       OUI
C  21 : QUADRILATERE Q1        4                       OUI
C  41 : PRISMES TELEMAC-3D     6                       OUI
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      X         |<-- |  SOLUTION DU SYSTEME AX = B
C |      A         |<-- |  MATRICE A SOUS FORME LDU
C |      B         |<-- |  SECOND MEMBRE DU SYSTEME A RESOUDRE.
C |      DITR      | -->|  CARACTERE  'D' : ON CALCULE AVEC A
C |                |    |             'T' : ON CALCULE AVEC A TRANSPOSEE
C |      MESH      | -->|  BLOC DES TABLEAUX D'ENTIERS DU MAILLAGE.
C |      COPY      | -->|  SI .TRUE. B EST RECOPIE SUR X.
C |                |    |  AU PREALABLE.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : DESCEN , REMONT , PLANTE
C
C**********************************************************************
C
      USE BIEF, EX_DWNUP1 => DWNUP1
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X
      TYPE(BIEF_OBJ), INTENT(IN)    :: B
      TYPE(BIEF_OBJ), INTENT(IN)    :: A
      TYPE(BIEF_MESH), INTENT(IN)   :: MESH
      CHARACTER(LEN=1), INTENT(IN)  :: DITR
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELM,NPOIN,NELEM,NELMAX
C
      DOUBLE PRECISION C
C
      CHARACTER*1 TYPD,TYPX
C
C-----------------------------------------------------------------------
C
      TYPD  = A%TYPDIA
      TYPX  = A%TYPEXT
      NPOIN = A%D%DIM1
      IELM  = A%ELMLIN
      NELEM = MESH%NELEM
      NELMAX= MESH%NELMAX
      CALL CPSTVC(B,X)
C
C-----------------------------------------------------------------------
C
C 1) DESCENTE AVEC RECOPIE DE B DANS X
C
      IF(A%STO.EQ.1) THEN
        CALL DESCEN(X%R, A%X%R,TYPX,B%R,
     *       MESH%IKLE%I,NELEM,NELMAX,NPOIN,IELM,DITR,.TRUE.,MESH%LV)
      ELSEIF(A%STO.EQ.3) THEN
        CALL DESSEG(X%R, A%X%R,TYPX,B%R,
     *              MESH%GLOSEG%I,MESH%NSEG,NPOIN,DITR,.TRUE.)
      ENDIF
C
C-----------------------------------------------------------------------
C
C 2) SUITE D'INVERSIONS DES MATRICES DIAGONALES
C
      IF(TYPD(1:1).NE.'I') THEN
        CALL OV( 'X=XY    ' , X%R , A%D%R , X%R , C , NPOIN )
      ENDIF
C
C-----------------------------------------------------------------------
C
C 3) REMONTEE SANS RECOPIE PREALABLE DE B DANS X
C
      IF(A%STO.EQ.1) THEN
        CALL REMONT(X%R, A%X%R,TYPX,B%R,
     *       MESH%IKLE%I,NELEM,NELMAX,NPOIN,IELM,DITR,.FALSE.,MESH%LV)
      ELSE
        CALL REMSEG(X%R, A%X%R,TYPX,B%R,
     *              MESH%GLOSEG%I,MESH%NSEG,NPOIN,DITR,.FALSE.)
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
