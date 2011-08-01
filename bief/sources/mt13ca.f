C                       *****************
                        SUBROUTINE MT13CA
C                       *****************
C
     *(  A11 , A12 , A13 ,
     *   A21 , A22 , A23 ,
     *   A31 , A32 , A33 ,
     *   A41 , A42 , A43 ,
     *   A51 , A52 , A53 ,
     *   A61 , A62 , A63 ,
     *   XMUL,XEL,YEL,NELEM,NELMAX,ICOORD)
C
C***********************************************************************
C BIEF VERSION 5.9           09/07/2008    ALGIANE FROEHLY
C                                          
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
C 
C
C  AVEC ICOORD=2 ON AURAIT UNE DERIVEE SUIVANT Y
C
C  PSI1 : BASES DE TYPE TRIANGLE P2
C  PSI2 : BASES DE TYPE TRIANGLE P1
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
      USE BIEF!, EX_MT13CA => MT13CA
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
      DOUBLE PRECISION, INTENT(INOUT) :: A51(*),A52(*),A53(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A61(*),A62(*),A63(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
      DOUBLE PRECISION X2,X3,Y2,Y3
      DOUBLE PRECISION XSUR6
C
C-----------------------------------------------------------------------
C
      XSUR6 = XMUL/6.D0
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
        A12(IELEM) = 0.D0
        A13(IELEM) = 0.D0
        A21(IELEM) = 0.D0
        A23(IELEM) = 0.D0
        A31(IELEM) = 0.D0
        A32(IELEM) = 0.D0
        A42(IELEM) =   Y3*XSUR6
        A43(IELEM) = - Y2*XSUR6
        A41(IELEM) = - A42(IELEM) - A43(IELEM)
        A51(IELEM) =   A41(IELEM)
        A52(IELEM) =   A42(IELEM)
        A53(IELEM) =   A43(IELEM)
        A61(IELEM) =   A41(IELEM)
        A62(IELEM) =   A42(IELEM)
        A63(IELEM) =   A43(IELEM)
C
C   TERMES DIAGONAUX
C
        A11(IELEM) = 0.D0
        A22(IELEM) = 0.D0
        A33(IELEM) = 0.D0
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
        A12(IELEM) = 0.D0
        A13(IELEM) = 0.D0
        A21(IELEM) = 0.D0
        A23(IELEM) = 0.D0
        A31(IELEM) = 0.D0
        A32(IELEM) = 0.D0
        A42(IELEM) = - X3*XSUR6
        A43(IELEM) =   X2*XSUR6
        A41(IELEM) = - A42(IELEM) - A43(IELEM)
        A51(IELEM) =   A41(IELEM)
        A52(IELEM) =   A42(IELEM)
        A53(IELEM) =   A43(IELEM)
        A61(IELEM) =   A41(IELEM)
        A62(IELEM) =   A42(IELEM)
        A63(IELEM) =   A43(IELEM)
C
C   TERMES DIAGONAUX
C
        A11(IELEM) = 0.D0
        A22(IELEM) = 0.D0
        A33(IELEM) = 0.D0
C
2       CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
          CALL PLANTE(1)
C
        ENDIF
C
200       FORMAT(1X,'MT13CA (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'MT13CA (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
C
C-----------------------------------------------------------------------
C
      RETURN
      END
