C                       *****************
                        SUBROUTINE ASSVE1
C                       *****************
C
     *(X, IKLE,W, NELEM,NELMAX,LV,MSK,MASKEL)
C
C***********************************************************************
C BIEF VERSION 5.1           18/08/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : BOUCLE D'ASSEMBLAGE DE VECTEUR.
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
C |  NNNNN         | -->| ANCIENNE DIMENSION DU TABLEAU X
C |                |    | (NE SERT PLUS, CONSERVE POUR COMPATIBILITE
C |                |    |  AVEC LES ANCIENNES VERSIONS).
C |  NELEM         | -->| NOMBRE D'ELEMENTS DANS LE MAILLAGE.
C |  NELMAX        | -->| PREMIERE DIMENSION DE IKLE ET W.
C |                |    | (CAS D'UN MAILLAGE ADAPTATIF)
C |  LV            | -->| LONGUEUR DU VECTEUR POUR LA VECTORISATION
C |  MSK           | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |  MASKEL        | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELANT : ASSVEC
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(INOUT) :: X(*)
      INTEGER         , INTENT(IN)    :: NELEM,NELMAX,LV
      INTEGER         , INTENT(IN)    :: IKLE(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: W(NELMAX),MASKEL(NELMAX)
      LOGICAL         , INTENT(IN)    :: MSK
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
C  AVEC MASQUAGE
C
      IF(MSK) THEN
C
      IF(LV.EQ.1) THEN
C
C  MODE SCALAIRE
C
      DO 10 IELEM = 1 , NELEM
        X(IKLE(IELEM)) = X(IKLE(IELEM)) + W(IELEM) * MASKEL(IELEM)
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
        X(IKLE(IELEM)) = X(IKLE(IELEM)) + W(IELEM) * MASKEL(IELEM)
30    CONTINUE
20    CONTINUE
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  SANS MASQUAGE
C
      ELSE
C
      IF(LV.EQ.1) THEN
C
C  MODE SCALAIRE
C
      DO 40 IELEM = 1 , NELEM
        X(IKLE(IELEM)) = X(IKLE(IELEM)) + W(IELEM)
40    CONTINUE
C
      ELSE
C
C  MODE VECTORIEL
C
      DO 60 IB = 1,(NELEM+LV-1)/LV
*VOCL LOOP,NOVREC
CDIR$ IVDEP
      DO 50 IELEM = 1+(IB-1)*LV , MIN(NELEM,IB*LV)
        X(IKLE(IELEM)) = X(IKLE(IELEM)) + W(IELEM)
50    CONTINUE
60    CONTINUE
C
      ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
