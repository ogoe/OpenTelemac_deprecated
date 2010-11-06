C                       *****************
                        SUBROUTINE MT06PP
C                       *****************
C
     *( T,XM,XMUL,SF,F,Z,SURFAC,IKLE,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.1          28/11/94    J-M HERVOUET (LNH) 30 87 80 18
C Arnaud Desitter - University of Bristol - April 1998
C***********************************************************************
C
C FONCTION : CONSTRUCTION DE LA MATRICE DE MASSE POUR DES PRISMES P1
C
C-----------------------------------------------------------------------
C
C      FONCTION:
C      ========:
C
C    CE SOUS-PROGRAMME CALCULE LES COEFFICIENTS DE LA MATRICE SUIVANTE:
C
C                              /
C                    A    =   /  F * (P *P )*J(X,Y) DXDY
C                     I J    /S        I  J
C
C    PAR MAILLE ELEMENTAIRE.
C
C     J(X,Y) : JACOBIEN DE LA TRANSFORMATION ISOPARAMETRIQUE
C
C     L'ELEMENT EST LE PRISME P1
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
C |     X,Y,Z      | -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |     SURFAC     | -->|  SURFACE DES TRIANGLES.
C |     IKLE1..6   | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |     NELEM      | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |     NELMAX     | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : AUCUN
C
C**********************************************************************
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_MT06PP => MT06PP
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE(NELMAX,6)
C
      DOUBLE PRECISION, INTENT(INOUT) :: T(NELMAX,6), XM(NELMAX,30)
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
      DOUBLE PRECISION, INTENT(IN) :: F(*)
C
C     STRUCTURE DE F
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF
C
      DOUBLE PRECISION, INTENT(IN) :: Z(*)
      DOUBLE PRECISION, INTENT(IN) :: SURFAC(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      DOUBLE PRECISION SUR2160,COEF,H1,H2,H3,F1,F2,F3,F4,F5,F6
C
      INTEGER IELEM,IELMF
C
C-----------------------------------------------------------------------
C
      SUR2160 = XMUL / 2160.D0
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
      IF(IELMF.NE.41) THEN
        IF (LNG.EQ.1) WRITE(LU,100) IELMF
        IF (LNG.EQ.2) WRITE(LU,101) IELMF
100     FORMAT(1X,'MT06PP (BIEF) :',/,
     *         1X,'DISCRETISATION DE F : ',1I6,' NON PREVUE')
101     FORMAT(1X,'MT06PP (BIEF) :',/,
     *         1X,'DISCRETIZATION OF F : ',1I6,' NOT AVAILABLE')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C BOUCLE SUR LES ELEMENTS
C
      DO IELEM = 1,NELEM
C
         COEF = SURFAC(IELEM)* SUR2160
C
         H1 = (Z(IKLE(IELEM,4)) - Z(IKLE(IELEM,1))) * COEF
         H2 = (Z(IKLE(IELEM,5)) - Z(IKLE(IELEM,2))) * COEF
         H3 = (Z(IKLE(IELEM,6)) - Z(IKLE(IELEM,3))) * COEF
C
         F1 = F(IKLE(IELEM,1))
         F2 = F(IKLE(IELEM,2))
         F3 = F(IKLE(IELEM,3))
         F4 = F(IKLE(IELEM,4))
         F5 = F(IKLE(IELEM,5))
         F6 = F(IKLE(IELEM,6))
C
C-----------------------------------------------------------------------
C
C  TERMES EXTRA-DIAGONAUX
C
      XM(IELEM,01) =
     &     (9*F1+6*F2+3*F3+3*F4+2*F5+F6)*h1+
     &     (6*F1+9*F2+3*F3+2*F4+3*F5+F6)*h2+
     &     (3*F1+3*F2+3*F3+F4+F5+F6)*h3
      XM(IELEM,02) =
     &     (9*F1+3*F2+6*F3+3*F4+F5+2*F6)*h1+
     &     (3*F1+3*F2+3*F3+F4+F5+F6)*h2+
     &     (6*F1+3*F2+9*F3+2*F4+F5+3*F6)*h3
      XM(IELEM,03) =
     &     (12*F1+3*F2+3*F3+12*F4+3*F5+3*F6)*h1+
     &     (3*F1+2*F2+F3+3*F4+2*F5+F6)*h2+
     &     (3*F1+F2+2*F3+3*F4+F5+2*F6)*h3
      XM(IELEM,04) =
     &     (3*F1+2*F2+F3+3*F4+2*F5+F6)*h1+
     &     (2*F1+3*F2+F3+2*F4+3*F5+F6)*h2+
     &     (F1+F2+F3+F4+F5+F6)*h3
      XM(IELEM,05) =
     &     (3*F1+F2+2*F3+3*F4+F5+2*F6)*h1+
     &     (F1+F2+F3+F4+F5+F6)*h2+
     &     (2*F1+F2+3*F3+2*F4+F5+3*F6)*h3
      XM(IELEM,06) =
     &     (3*F1+3*F2+3*F3+F4+F5+F6)*h1+
     &     (3*F1+9*F2+6*F3+F4+3*F5+2*F6)*h2+
     &     (3*F1+6*F2+9*F3+F4+2*F5+3*F6)*h3
      XM(IELEM,07) =
     &     (3*F1+2*F2+F3+3*F4+2*F5+F6)*h1+
     &     (2*F1+3*F2+F3+2*F4+3*F5+F6)*h2+
     &     (F1+F2+F3+F4+F5+F6)*h3
      XM(IELEM,08) =
     &     (2*F1+3*F2+F3+2*F4+3*F5+F6)*h1+
     &     (3*F1+12*F2+3*F3+3*F4+12*F5+3*F6)*h2+
     &     (F1+3*F2+2*F3+F4+3*F5+2*F6)*h3
      XM(IELEM,09) =
     &     (F1+F2+F3+F4+F5+F6)*h1+
     &     (F1+3*F2+2*F3+F4+3*F5+2*F6)*h2+
     &     (F1+2*F2+3*F3+F4+2*F5+3*F6)*h3
      XM(IELEM,10) =
     &     (3*F1+F2+2*F3+3*F4+F5+2*F6)*h1+
     &     (F1+F2+F3+F4+F5+F6)*h2+
     &     (2*F1+F2+3*F3+2*F4+F5+3*F6)*h3
      XM(IELEM,11) =
     &     (F1+F2+F3+F4+F5+F6)*h1+
     &     (F1+3*F2+2*F3+F4+3*F5+2*F6)*h2+
     &     (F1+2*F2+3*F3+F4+2*F5+3*F6)*h3
      XM(IELEM,12) =
     &     (2*F1+F2+3*F3+2*F4+F5+3*F6)*h1+
     &     (F1+2*F2+3*F3+F4+2*F5+3*F6)*h2+
     &     (3*F1+3*F2+12*F3+3*F4+3*F5+12*F6)*h3
      XM(IELEM,13) =
     &     (3*F1+2*F2+F3+9*F4+6*F5+3*F6)*h1+
     &     (2*F1+3*F2+F3+6*F4+9*F5+3*F6)*h2+
     &     (F1+F2+F3+3*F4+3*F5+3*F6)*h3
      XM(IELEM,14) =
     &     (3*F1+F2+2*F3+9*F4+3*F5+6*F6)*h1+
     &     (F1+F2+F3+3*F4+3*F5+3*F6)*h2+
     &     (2*F1+F2+3*F3+6*F4+3*F5+9*F6)*h3
      XM(IELEM,15) =
     &     (F1+F2+F3+3*F4+3*F5+3*F6)*h1+
     &     (F1+3*F2+2*F3+3*F4+9*F5+6*F6)*h2+
     &     (F1+2*F2+3*F3+3*F4+6*F5+9*F6)*h3
C
C  TERMES DIAGONAUX
C
      T(IELEM,1) =
     &     (36*F1+9*F2+9*F3+12*F4+3*F5+3*F6)*h1+
     &     (9*F1+6*F2+3*F3+3*F4+2*F5+F6)*h2+
     &     (9*F1+3*F2+6*F3+3*F4+F5+2*F6)*h3
      T(IELEM,2) =
     &     (6*F1+9*F2+3*F3+2*F4+3*F5+F6)*h1+
     &     (9*F1+36*F2+9*F3+3*F4+12*F5+3*F6)*h2+
     &     (3*F1+9*F2+6*F3+F4+3*F5+2*F6)*h3
      T(IELEM,3) =
     &     (6*F1+3*F2+9*F3+2*F4+F5+3*F6)*h1+
     &     (3*F1+6*F2+9*F3+F4+2*F5+3*F6)*h2+
     &     (9*F1+9*F2+36*F3+3*F4+3*F5+12*F6)*h3
      T(IELEM,4) =
     &     (12*F1+3*F2+3*F3+36*F4+9*F5+9*F6)*h1+
     &     (3*F1+2*F2+F3+9*F4+6*F5+3*F6)*h2+
     &     (3*F1+F2+2*F3+9*F4+3*F5+6*F6)*h3
      T(IELEM,5) =
     &     (2*F1+3*F2+F3+6*F4+9*F5+3*F6)*h1+
     &     (3*F1+12*F2+3*F3+9*F4+36*F5+9*F6)*h2+
     &     (F1+3*F2+2*F3+3*F4+9*F5+6*F6)*h3
      T(IELEM,6) =
     &     (2*F1+F2+3*F3+6*F4+3*F5+9*F6)*h1+
     &     (F1+2*F2+3*F3+3*F4+6*F5+9*F6)*h2+
     &     (3*F1+3*F2+12*F3+9*F4+9*F5+36*F6)*h3
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
