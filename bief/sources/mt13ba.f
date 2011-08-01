C                       *****************
                        SUBROUTINE MT13BA
C                       *****************
C
     *(  A11 , A12 , A13 ,
     *   A21 , A22 , A23 ,
     *   A31 , A32 , A33 ,
     *   A41 , A42 , A43 ,
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
C  PSI1 : BASES DE TYPE TRIANGLE QUASI-BULLE
C  PSI2 : BASES DE TYPE IELM2
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
      USE BIEF, EX_MT13BA => MT13BA
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX,ICOORD
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A21(*),A22(*),A23(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A31(*),A32(*),A33(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A41(*),A42(*),A43(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
      DOUBLE PRECISION X2,X3,Y2,Y3
      DOUBLE PRECISION XSUR9,XSUR6
C
C-----------------------------------------------------------------------
C
      XSUR6 = XMUL/6.D0
      XSUR9 = XMUL/9.D0
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
        A12(IELEM) =  Y3*XSUR9
        A13(IELEM) = -Y2*XSUR9
        A21(IELEM) = -(Y3-Y2)*XSUR9
        A23(IELEM) = -Y2*XSUR9
        A31(IELEM) = -(Y3-Y2)*XSUR9
        A32(IELEM) =  Y3*XSUR9
        A41(IELEM) = -(Y3-Y2)*XSUR6
        A42(IELEM) =  Y3*XSUR6
        A43(IELEM) = -Y2*XSUR6
C
C   TERMES DIAGONAUX
C   LA SOMME DES COLONNES DE LA MATRICE EST LE VECTEUR NUL
C
        A11(IELEM) = - A12(IELEM) - A13(IELEM)
        A22(IELEM) = - A21(IELEM) - A23(IELEM)
        A33(IELEM) = - A31(IELEM) - A32(IELEM)
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
        A12(IELEM) = -X3*XSUR9
        A13(IELEM) =  X2*XSUR9
        A21(IELEM) = -(X2-X3)*XSUR9
        A23(IELEM) =  X2*XSUR9
        A31(IELEM) = -(X2-X3)*XSUR9
        A32(IELEM) = -X3*XSUR9
        A41(IELEM) = -(X2-X3)*XSUR6
        A42(IELEM) = -X3*XSUR6
        A43(IELEM) =  X2*XSUR6
C
C   TERMES DIAGONAUX
C   LA SOMME DES COLONNES DE LA MATRICE EST LE VECTEUR NUL
C
        A11(IELEM) = - A12(IELEM) - A13(IELEM)
        A22(IELEM) = - A21(IELEM) - A23(IELEM)
        A33(IELEM) = - A31(IELEM) - A32(IELEM)
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
200       FORMAT(1X,'MT13BA (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'MT13BA (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
