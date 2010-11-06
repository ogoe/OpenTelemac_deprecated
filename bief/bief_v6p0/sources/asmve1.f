C                       *****************
                        SUBROUTINE ASMVE1
C                       *****************
C
     *(X, IKLE,W, NPOIN,NELEM,NELMAX,LV)
C
C***********************************************************************
C BIEF VERSION 5.1           17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : BOUCLE D'ASSEMBLAGE MULTIPLICATIF DE VECTEUR.
C                                -------------
C
C CETTE FORME D'ASSEMBLAGE EST UTILISEE AVEC DES PRECONDITIONNEMENTS DE
C TYPE 'CROUT' ELEMENTAIRE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  X             |<-->| VECTEUR ASSEMBLE
C |  IKLE          | -->| CORRESPONDANCE NUMEROTATION LOCALE-GLOBALE
C |  W             | -->| TABLEAUX DE TRAVAIL CONTENANT LE VECTEUR SOUS
C |                |    | FORME NON ASSEMBLEE
C |  NPOIN         | -->| DIMENSION DU TABLEAU X
C |  NELEM         | -->| NOMBRE D'ELEMENTS DANS LE MAILLAGE.
C |  NELMAX        | -->| PREMIERE DIMENSION DE IKLE ET W.
C |                |    | (CAS D'UN MAILLAGE ADAPTATIF)
C |  LV            | -->| LONGUEUR DU VECTEUR POUR LA VECTORISATION
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELANT : ASMVEC
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)    :: NELEM,NELMAX,NPOIN,LV
      DOUBLE PRECISION, INTENT(INOUT) :: X(NPOIN)
      INTEGER         , INTENT(IN)    :: IKLE(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: W(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IB
C
      INTRINSIC MIN
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
        X(IKLE(IELEM)) = X(IKLE(IELEM)) * W(IELEM)
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
        X(IKLE(IELEM)) = X(IKLE(IELEM)) * W(IELEM)
30    CONTINUE
20    CONTINUE
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
