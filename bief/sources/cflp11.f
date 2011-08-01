C                       *****************
                        SUBROUTINE CFLP11
C                       *****************
C
     *(U,V,X,Y,IKLE,NELEM,NELMAX,W1)
C
C***********************************************************************
C BIEF VERSION 5.1           17/08/94      C MOULIN   (LNH) 30 87 83 81
C                                          + MODIFS JMH LE 17/08/94
C***********************************************************************
C
C  FONCTION  : CALCULE LE NOMBRE DE COURANT EN CHAQUE POINT DU MAILLAGE
C              POUR CHAQUE PAS DE TEMPS.
C
C              LE CRITERE DE STABILITE DU SCHEMA DISTRIBUTIF N
C              EST ICI UTILISE COMME EVALUATION DU NOMBRE DE COURANT.
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      U         | -->| VITESSE SUIVANT X.                           |
C |      V         | -->| VITESSE SUIVANT Y.                           |
C |      X         | -->| ABSCISSES DES POINTS DU MAILLAGE PAR ELEMENTS|
C |      Y         | -->| ORDONNEES DES POINTS DU MAILLAGE PAR ELEMENTS|
C |      IKLE      | -->| NUMEROS DES NOEUDS DE CHAQUE ELEMENT.        |
C |      NELEM     | -->| NOMBRE D'ELEMENTS DU MAILLAGE.               |
C |      NELMAX    | -->| NOMBRE D'ELEMENTS MAXIMUM DU MAILLAGE        |
C |                |    | (CAS D'UN MAILLAGE ADAPTATIF).               |
C |      W1        | -->| RESULTAT PARTIEL                             |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : TELMAC
C
C SOUS-PROGRAMME APPELE : LISSAG
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)  :: NELEM,NELMAX
      DOUBLE PRECISION, INTENT(IN)  :: U(*),V(*)
      DOUBLE PRECISION, INTENT(IN)  :: X(NELMAX,*),Y(NELMAX,*)
      INTEGER         , INTENT(IN)  :: IKLE(NELMAX,*)
      DOUBLE PRECISION, INTENT(OUT) :: W1(NELMAX,*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
C
      DOUBLE PRECISION U1,U2,U3,V1,V2,V3,USUR2,VSUR2
      DOUBLE PRECISION SUR6,K1,K2,K3,L12,L13,L21,L23,L31,L32
      DOUBLE PRECISION X2,X3,Y2,Y3
C
      INTRINSIC MAX,MIN
C
C-----------------------------------------------------------------------
C
C
      SUR6 = 1.D0 / 6.D0
C
C BOUCLE SUR LES ELEMENTS
C
        DO 10 IELEM = 1, NELEM
C
          X2 = X(IELEM,2)
          X3 = X(IELEM,3)
          Y2 = Y(IELEM,2)
          Y3 = Y(IELEM,3)
C
          U1 = U(IKLE(IELEM,1))
          U2 = U(IKLE(IELEM,2))
          U3 = U(IKLE(IELEM,3))
          V1 = V(IKLE(IELEM,1))
          V2 = V(IKLE(IELEM,2))
          V3 = V(IKLE(IELEM,3))
C
          USUR2 = (U1+U2+U3)*SUR6
          VSUR2 = (V1+V2+V3)*SUR6
C
          K1 = USUR2 * (Y2-Y3) - VSUR2 * (X2-X3)
          K2 = USUR2 * (Y3   ) - VSUR2 * (X3   )
          K3 = USUR2 * (  -Y2) - VSUR2 * (  -X2)
C
          L12 = MAX(  MIN(K1,-K2) , 0.D0 )
          L13 = MAX(  MIN(K1,-K3) , 0.D0 )
          L21 = MAX(  MIN(K2,-K1) , 0.D0 )
          L23 = MAX(  MIN(K2,-K3) , 0.D0 )
          L31 = MAX(  MIN(K3,-K1) , 0.D0 )
          L32 = MAX(  MIN(K3,-K2) , 0.D0 )
C
          W1(IELEM,1) = L12 + L13
          W1(IELEM,2) = L21 + L23
          W1(IELEM,3) = L31 + L32
C
10      CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
