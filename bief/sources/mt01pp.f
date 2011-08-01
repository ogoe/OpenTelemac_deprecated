C                       *****************
                        SUBROUTINE MT01PP
C                       *****************
C
     *( T,XM,XMUL,Z,SURFAC,IKLE,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 6.0      21/05/2010    J-M HERVOUET (LNHE) 01 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
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
C                                 /
C                    A    = XMUL /  (P *P )*J(X,Y) DXDY
C                     I J       /S    I  J
C
C    PAR MAILLE ELEMENTAIRE.
C
C     J(X,Y) : JACOBIEN DE LA TRANSFORMATION ISOPARAMETRIQUE
C
C     L'ELEMENT EST LE TRIANGLE P1
C
C    MAILLAGE REEL
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     A11,A12... |<-- |  ELEMENTS DE LA MATRICE
C |     XMUL       | -->|  FACTEUR MULTIPLICATIF
C |         Z      | -->|  COORDONNEES DES POINTS DANS L'ELEMENT
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
      USE BIEF, EX_MT01PP => MT01PP
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NELEM,NELMAX
      INTEGER, INTENT(IN)             :: IKLE(NELMAX,6)
      DOUBLE PRECISION, INTENT(INOUT) :: T(NELMAX,6),XM(NELMAX,30)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL,Z(*),SURFAC(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      INTEGER IELEM
      DOUBLE PRECISION SUR360,COEF,H1,H2,H3,HT,Z41,Z52,Z63
C
      DOUBLE PRECISION EPSILON
      DATA EPSILON/1.D-3/
C
C-----------------------------------------------------------------------
C
      SUR360 = XMUL / 360.D0
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 1 IELEM = 1,NELEM
C
         COEF = SURFAC(IELEM) * SUR360
C
C        ICI TRAITEMENT DES ZONES SECHES
C
         H1 = MAX(Z(IKLE(IELEM,4)) - Z(IKLE(IELEM,1)),EPSILON) * COEF
         H2 = MAX(Z(IKLE(IELEM,5)) - Z(IKLE(IELEM,2)),EPSILON) * COEF
         H3 = MAX(Z(IKLE(IELEM,6)) - Z(IKLE(IELEM,3)),EPSILON) * COEF
         HT = H1 + H2 + H3
C
C-----------------------------------------------------------------------
C
C  TERMES EXTRA-DIAGONAUX
C
         XM(IELEM,4)  = H1 + H2 + HT
         XM(IELEM,5)  = H1 + H3 + HT
         XM(IELEM,9)  = H2 + H3 + HT
         XM(IELEM,7)  = XM(IELEM,4)
         XM(IELEM,10) = XM(IELEM,5)
         XM(IELEM,11) = XM(IELEM,9)
C
         XM(IELEM,3)  =  4*H1 + HT + HT
         XM(IELEM,8)  =  4*H2 + HT + HT
         XM(IELEM,12) =  4*H3 + HT + HT
C
         XM(IELEM,1)  = XM(IELEM,4) + XM(IELEM,4)
         XM(IELEM,2)  = XM(IELEM,5) + XM(IELEM,5)
         XM(IELEM,6)  = XM(IELEM,9) + XM(IELEM,9)
         XM(IELEM,13) = XM(IELEM,1)
         XM(IELEM,14) = XM(IELEM,2)
         XM(IELEM,15) = XM(IELEM,6)
C
C  TERMES DIAGONAUX
C
         T(IELEM,1) = XM(IELEM,3)  + XM(IELEM,3)
         T(IELEM,2) = XM(IELEM,8)  + XM(IELEM,8)
         T(IELEM,3) = XM(IELEM,12) + XM(IELEM,12)
         T(IELEM,4) = T(IELEM,1)
         T(IELEM,5) = T(IELEM,2)
         T(IELEM,6) = T(IELEM,3)
C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
1     CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
