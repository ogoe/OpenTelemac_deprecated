C                       ****************
                        SUBROUTINE REM21
C                       ****************
C
     *(X, XA1,XA2,XA3,XA4,XA5,XA6 , IKLE1,IKLE2,IKLE3,IKLE4,
     * NELEM,NELMAX,NPOIN,LV)
C
C***********************************************************************
C BIEF VERSION 5.1           05/02/91    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : PRODUIT X = U B (ATTENTION : ELEMENT PAR ELEMENT)
C
C            ICI ELEMENT Q1 OU ELEMENT A QUATRE POINTS
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
C |      X         |<-->|  VECTEUR DONNEE ET RESULTAT.
C |      XA1,..6   | -->|  TERMES EXTRADIAGONAUX DE LA MATRICE A
C |                |    |  CORRESPONDANT A LA PARTIE INFERIEURE
C |      IKLE1,2,34| -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
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
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX)
      INTEGER, INTENT(IN) :: IKLE3(NELMAX),IKLE4(NELMAX)
C
      DOUBLE PRECISION, INTENT(INOUT) :: X(NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: XA1(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XA2(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XA3(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XA4(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XA5(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XA6(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IB
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
      DO 10 IELEM = NELEM , 1 , -1
        X(IKLE3(IELEM))=X(IKLE3(IELEM))-XA6(IELEM)*X(IKLE4(IELEM))
        X(IKLE2(IELEM))=X(IKLE2(IELEM))-XA5(IELEM)*X(IKLE4(IELEM))
     *                                 -XA4(IELEM)*X(IKLE3(IELEM))
        X(IKLE1(IELEM))=X(IKLE1(IELEM))-XA3(IELEM)*X(IKLE4(IELEM))
     *                                 -XA2(IELEM)*X(IKLE3(IELEM))
     *                                 -XA1(IELEM)*X(IKLE2(IELEM))
10    CONTINUE
C
      ELSE
C
C  MODE VECTORIEL
C
      DO 20 IB = (NELEM+LV-1)/LV , 1 , -1
*VOCL LOOP,NOVREC
CDIR$ IVDEP
      DO 30 IELEM = MIN(NELEM,IB*LV) , 1+(IB-1)*LV , -1
        X(IKLE3(IELEM))=X(IKLE3(IELEM))-XA6(IELEM)*X(IKLE4(IELEM))
        X(IKLE2(IELEM))=X(IKLE2(IELEM))-XA5(IELEM)*X(IKLE4(IELEM))
     *                                 -XA4(IELEM)*X(IKLE3(IELEM))
        X(IKLE1(IELEM))=X(IKLE1(IELEM))-XA3(IELEM)*X(IKLE4(IELEM))
     *                                 -XA2(IELEM)*X(IKLE3(IELEM))
     *                                 -XA1(IELEM)*X(IKLE2(IELEM))
30    CONTINUE
20    CONTINUE
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
