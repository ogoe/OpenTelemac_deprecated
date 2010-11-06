C                       *****************
                        SUBROUTINE CFLP12
C                       *****************
C
     *(U,V,X,Y,IKLE,NELEM,NELMAX,W1)
C
C***********************************************************************
C BIEF VERSION 5.6           29/12/05      C MOULIN   (LNH) 30 87 83 81
C                                          + MODIFS JMH LE 29/12/05
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
      DOUBLE PRECISION, INTENT(IN)  :: X(NELMAX*3),Y(NELMAX*3)
      INTEGER         , INTENT(IN)  :: IKLE(NELMAX*4)
      DOUBLE PRECISION, INTENT(OUT) :: W1(NELMAX*4)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IT,IAD1,IAD2,IAD3,IG1,IG2,IG3
C
      DOUBLE PRECISION USUR2,VSUR2
      DOUBLE PRECISION SUR6,K1,K2,K3,L12,L13,L21,L23,L31,L32
      DOUBLE PRECISION X1,X2,X3,Y1,Y2,Y3,TIERS
C
      INTRINSIC MAX,MIN
C
C-----------------------------------------------------------------------
C
C     POUR UN TRIANGLE QUASI-BULLE : NUMEROS DES SOMMETS DES
C     SOUS-TRIANGLES DANS LE TRIANGLE INITIAL.
C     IL(NUMERO DU SOUS-TRIANGLE,NUMERO LOCAL DANS LE SOUS-TRIANGLE)
C
      INTEGER IL(3,3)
      DATA IL /1,2,3,2,3,1,4,4,4/
C
C-----------------------------------------------------------------------
C
      TIERS= 1.D0 / 3.D0
      SUR6 = 1.D0 / 6.D0
C
C     INITIALISATIONS DES W
C
      DO 32 IELEM = 1 , 4*NELMAX
        W1(IELEM) = 0.D0
32    CONTINUE
C
C     CALCUL SELON LE SCHEMA PSI, BOUCLE SUR LES TROIS SOUS-TRIANGLES
C     ET PREASSEMBLAGE.
C
      DO 10 IT=1,3
CDIR$ IVDEP
      DO 33 IELEM = 1 , NELEM
C
C       ADRESSES DANS UN TABLEAU (NELMAX,*)
        IAD1= IELEM + (IL(IT,1)-1)*NELMAX
        IAD2= IELEM + (IL(IT,2)-1)*NELMAX
        IAD3= IELEM + (IL(IT,3)-1)*NELMAX
C       NUMEROS GLOBAUX DANS LE TRIANGLE INITIAL
        IG1 = IKLE(IAD1)
        IG2 = IKLE(IAD2)
        IG3 = IKLE(IAD3)
C       COORDONNEES DES SOMMETS DU SOUS-TRIANGLE
        X1 = X(IAD1)
        X2 = X(IAD2) - X1
        Y1 = Y(IAD1)
        Y2 = Y(IAD2) - Y1
C       LE POINT 3 EST TOUJOURS LE CENTRE DU TRIANGLE INITIAL
        X3=TIERS*(X(IELEM)+X(IELEM+NELMAX)+X(IELEM+2*NELMAX))-X1
        Y3=TIERS*(Y(IELEM)+Y(IELEM+NELMAX)+Y(IELEM+2*NELMAX))-Y1
C
        USUR2 = (U(IG1)+U(IG2)+U(IG3))*SUR6
        VSUR2 = (V(IG1)+V(IG2)+V(IG3))*SUR6
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
        W1(IAD1) = W1(IAD1) + L12 + L13
        W1(IAD2) = W1(IAD2) + L21 + L23
        W1(IAD3) = W1(IAD3) + L31 + L32
C
33    CONTINUE
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
