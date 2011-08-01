C                       ****************
                        SUBROUTINE MER41
C                       ****************
C
     *(X, XA1 ,XA2 ,XA3 ,XA4 ,XA5 ,
     *    XA6 ,XA7 ,XA8 ,XA9 ,XA10,
     *    XA11,XA12,XA13,XA14,XA15,
     *    IKLE1,IKLE2,IKLE3,IKLE4,IKLE5,IKLE6,
     *    NELEM,NELMAX,NPOIN,LV)
C
C***********************************************************************
C BIEF VERSION 5.1           05/02/91    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : PRODUIT X = U B (ATTENTION : ELEMENT PAR ELEMENT)
C
C            ICI ELEMENT PRISME P1 OU ELEMENT A SIX POINTS
C
C            OPERATION INVERSE DU SOUS-PROGRAMME REMONT, D'OU LE NOM
C
C            ICI LA MATRICE U EST LE RESULTAT D'UNE DECOMPOSITION
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
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      XA1,..15  | -->|  TERMES EXTRADIAGONAUX DE LA MATRICE A
C |                |    |  CORRESPONDANT A LA PARTIE INFERIEURE
C |      IKLE1,.,6 | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      NPOIN     | -->|  DIMENSION DES TABLEAUX
C |      LV        | -->|  LONGUEUR DU VECTEUR POUR LA VECTORISATION
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES :
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPOIN,NELEM,NELMAX,LV
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX),IKLE3(NELMAX)
      INTEGER, INTENT(IN) :: IKLE4(NELMAX),IKLE5(NELMAX),IKLE6(NELMAX)
C
      DOUBLE PRECISION, INTENT(INOUT) :: X(NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: XA1(NELMAX),XA2(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XA3(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XA4(NELMAX),XA5(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XA6(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XA7(NELMAX),XA8(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XA9(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XA10(NELMAX),XA11(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XA12(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XA13(NELMAX),XA14(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XA15(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IB
C
      INTRINSIC MIN
C
C-----------------------------------------------------------------------
C
C SUITE D'INVERSIONS DES MATRICES TRIANGULAIRES INFERIEURES
C
C-----------------------------------------------------------------------
C BOUCLE EN MODE SCALAIRE (LV=1) OU AVEC VECTORISATION FORCEE
C-----------------------------------------------------------------------
C
      IF(LV.EQ.1) THEN
C
C  MODE SCALAIRE
C
      DO 10 IELEM = 1 , NELEM
          X(IKLE1(IELEM))=X(IKLE1(IELEM))+XA5 (IELEM)*X(IKLE6(IELEM))
     *                                   +XA4 (IELEM)*X(IKLE5(IELEM))
     *                                   +XA3 (IELEM)*X(IKLE4(IELEM))
     *                                   +XA2 (IELEM)*X(IKLE3(IELEM))
     *                                   +XA1 (IELEM)*X(IKLE2(IELEM))
          X(IKLE2(IELEM))=X(IKLE2(IELEM))+XA9 (IELEM)*X(IKLE6(IELEM))
     *                                   +XA8 (IELEM)*X(IKLE5(IELEM))
     *                                   +XA7 (IELEM)*X(IKLE4(IELEM))
     *                                   +XA6 (IELEM)*X(IKLE3(IELEM))
          X(IKLE3(IELEM))=X(IKLE3(IELEM))+XA12(IELEM)*X(IKLE6(IELEM))
     *                                   +XA11(IELEM)*X(IKLE5(IELEM))
     *                                   +XA10(IELEM)*X(IKLE4(IELEM))
          X(IKLE4(IELEM))=X(IKLE4(IELEM))+XA14(IELEM)*X(IKLE6(IELEM))
     *                                   +XA13(IELEM)*X(IKLE5(IELEM))
          X(IKLE5(IELEM))=X(IKLE5(IELEM))+XA15(IELEM)*X(IKLE6(IELEM))
10    CONTINUE
C
      ELSE
C
C  MODE VECTORIEL
C
      DO 20 IB = 1,(NELEM+LV-1)/LV
*VOCL LOOP,NOVREC
CDIR$ IVDEP
      DO 30 IELEM = 1+(IB-1)*LV , MIN(NELEM,IB*LV)
          X(IKLE1(IELEM))=X(IKLE1(IELEM))+XA5 (IELEM)*X(IKLE6(IELEM))
     *                                   +XA4 (IELEM)*X(IKLE5(IELEM))
     *                                   +XA3 (IELEM)*X(IKLE4(IELEM))
     *                                   +XA2 (IELEM)*X(IKLE3(IELEM))
     *                                   +XA1 (IELEM)*X(IKLE2(IELEM))
          X(IKLE2(IELEM))=X(IKLE2(IELEM))+XA9 (IELEM)*X(IKLE6(IELEM))
     *                                   +XA8 (IELEM)*X(IKLE5(IELEM))
     *                                   +XA7 (IELEM)*X(IKLE4(IELEM))
     *                                   +XA6 (IELEM)*X(IKLE3(IELEM))
          X(IKLE3(IELEM))=X(IKLE3(IELEM))+XA12(IELEM)*X(IKLE6(IELEM))
     *                                   +XA11(IELEM)*X(IKLE5(IELEM))
     *                                   +XA10(IELEM)*X(IKLE4(IELEM))
          X(IKLE4(IELEM))=X(IKLE4(IELEM))+XA14(IELEM)*X(IKLE6(IELEM))
     *                                   +XA13(IELEM)*X(IKLE5(IELEM))
          X(IKLE5(IELEM))=X(IKLE5(IELEM))+XA15(IELEM)*X(IKLE6(IELEM))
30    CONTINUE
20    CONTINUE
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
