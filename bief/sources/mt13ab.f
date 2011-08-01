C                       *****************
                        SUBROUTINE MT13AB
C                       *****************
C
     *(  A11 , A12 , A13 , A14 ,
     *   A21 , A22 , A23 , A24 ,
     *   A31 , A32 , A33 , A34 ,
     *   XMUL,XEL,YEL,NELEM,NELMAX,ICOORD)
C
C***********************************************************************
C BIEF VERSION 5.1           06/02/95    J-M HERVOUET (LNH) 30 87 80 18
C                                          C MOULIN   (LNH) 30 87 83 81
C***********************************************************************
C
C  FONCTION : CALCUL DES COEFFICIENTS DE LA MATRICE SUIVANTE
C
C  EXEMPLE AVEC ICOORD = 1
C
C
C                  /           D
C A(I,J)= XMUL *  /  PSI1(I) * --( PSI2(J) ) D(OMEGA)
C                /OMEGA        DX
C
C ATTENTION : L'EXPRESSION DE CETTE MATRICE EST LA TRANSPOSEE
C             DE CELLE DE MT13BB
C
C  AVEC ICOORD=2 ON AURAIT UNE DERIVEE SUIVANT Y
C
C  PSI1 : BASES DE TYPE TRIANGLE LINEAIRE
C  PSI2 : BASES DE TYPE TRIANGLE QUASI-BULLE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     A11,A12... |<-- |  ELEMENTS DE LA MATRICE
C |     XMUL       | -->|  FACTEUR MULTIPLICATIF
C |     SF,SG,SH   | -->|  STRUCTURES DE F,G ET H.
C |     SU,SV,SW   | -->|  STRUCTURES DE U,V ET W.
C |     F,G,H      | -->|  FONCTIONS INTERVENANT DANS LE CALCUL DE LA
C |                |    |  MATRICE.
C |     U,V,W      | -->|  COMPOSANTES D'UN VECTEUR INTERVENANT DANS LE
C |                |    |  CALCUL DE LA MATRICE.
C |     XEL,YEL,ZEL| -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |     SURFAC     | -->|  SURFACE DES TRIANGLES.
C |     IKLE1..3   | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |     NELEM      | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |     NELMAX     | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |     ICOORD     | -->|  1: DERIVEE SUIVANT X, 2:SUIVANT Y
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : ASSVEC , OV
C
C-----------------------------------------------------------------------
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_MT13AB => MT13AB
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX,ICOORD
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*),A14(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A21(*),A22(*),A23(*),A24(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A31(*),A32(*),A33(*),A34(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
      DOUBLE PRECISION X2,X3,Y2,Y3
      DOUBLE PRECISION XSUR6,XSUR18
C
C-----------------------------------------------------------------------
C
      XSUR6  = XMUL/6.D0
      XSUR18 = XMUL/18.D0
C
C================================
C  CAS DE LA DERIVEE SUIVANT X  =
C================================
C
        IF(ICOORD.EQ.1) THEN
C
C   BOUCLE SUR LES ELEMENTS
C
        DO 1 IELEM = 1 , NELEM
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM)=(    Y2 + 2*Y3 ) * XSUR18
        A13(IELEM)=( -2*Y2 -   Y3 ) * XSUR18
        A14(IELEM)=( -  Y2 +   Y3 ) * XSUR6
        A21(IELEM)=(  3*Y2 - 2*Y3 ) * XSUR18
        A23(IELEM)=( -3*Y2 +   Y3 ) * XSUR18
        A24(IELEM)=(       -   Y3 ) * XSUR6
        A31(IELEM)=(  2*Y2 - 3*Y3 ) * XSUR18
        A32(IELEM)=( -  Y2 + 3*Y3 ) * XSUR18
        A34(IELEM)=     Y2          * XSUR6
C
C   TERMES DIAGONAUX
C   LA SOMME DE CHAQUE LIGNE EST NULLE
C
        A11(IELEM) = - A12(IELEM) - A13(IELEM) - A14(IELEM)
        A22(IELEM) = - A21(IELEM) - A23(IELEM) - A24(IELEM)
        A33(IELEM) = - A31(IELEM) - A32(IELEM) - A34(IELEM)
C
1     CONTINUE
C
        ELSEIF(ICOORD.EQ.2) THEN
C
C================================
C  CAS DE LA DERIVEE SUIVANT Y  =
C================================
C
        DO 2 IELEM = 1 , NELEM
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
        X2  =  XEL(IELEM,2)
        X3  =  XEL(IELEM,3)
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM)=( -  X2 - 2*X3 ) * XSUR18
        A13(IELEM)=(  2*X2 +   X3 ) * XSUR18
        A14(IELEM)=(    X2 -   X3 ) * XSUR6
        A21(IELEM)=( -3*X2 + 2*X3 ) * XSUR18
        A23(IELEM)=(  3*X2 -   X3 ) * XSUR18
        A24(IELEM)=            X3   * XSUR6
        A31(IELEM)=( -2*X2 + 3*X3 ) * XSUR18
        A32(IELEM)=(    X2 - 3*X3 ) * XSUR18
        A34(IELEM)=( -  X2        ) * XSUR6
C
C   TERMES DIAGONAUX
C   LA SOMME DE CHAQUE LIGNE DE LA MATRICE EST NULLE
C
        A11(IELEM) = - A12(IELEM) - A13(IELEM) - A14(IELEM)
        A22(IELEM) = - A21(IELEM) - A23(IELEM) - A24(IELEM)
        A33(IELEM) = - A31(IELEM) - A32(IELEM) - A34(IELEM)
C
2       CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
          CALL PLANTE(0)
C
        ENDIF
C
200       FORMAT(1X,'MT13AB (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'MT13AB (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
