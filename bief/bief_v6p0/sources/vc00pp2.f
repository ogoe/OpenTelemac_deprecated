C                       ******************
                        SUBROUTINE VC00PP2
C                       ******************
C
     *( XMUL,Z,SURFAC,
     *  IKLE1,IKLE2,IKLE3,IKLE4,IKLE5,IKLE6,NELEM,NELMAX,
     *  W1,W2,W3,W4,W5,W6)
C
C***********************************************************************
C BIEF VERSION 5.5         11/10/04    J-M HERVOUET (LNH) 01 30 87 80 18
C                                      (IDEA AND THEORY: ASTRID DECOENE)
C***********************************************************************
C
C    FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C                    /
C    VEC(I) = XMUL  /    PSI(I)  D(OMEGA)
C                  /OMEGA
C
C    IN A FORM COMPATIBLE WITH MASS-LUMPING IN 2D   !!!!!!!!!!!!!!!!!!!
C    (COMPARE WITH VC00PP)
C
C    PSI(I) EST UNE BASE DE TYPE PRISME P1
C
C
C    ATTENTION : LE RESULTAT EST DANS W  SOUS FORME NON ASSEMBLEE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      XMUL      | -->|  COEFFICIENT MULTIPLICATEUR.
C |      Z         | -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |      SURFAC    | -->|  SURFACE DES ELEMENTS.
C |      IKLE1,... | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE.
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE.
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE.
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      W1,2,3    |<-- |  VECTEUR RESULTAT SOUS FORME NON ASSEMBLEE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C  PROGRAMMES APPELES : ASSVEC
C
C-----------------------------------------------------------------------
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX),IKLE3(NELMAX)
      INTEGER, INTENT(IN) :: IKLE4(NELMAX),IKLE5(NELMAX),IKLE6(NELMAX)
C
      DOUBLE PRECISION, INTENT(IN) :: Z(*),SURFAC(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT)::W1(NELMAX),W2(NELMAX),W3(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT)::W4(NELMAX),W5(NELMAX),W6(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER IELEM
      DOUBLE PRECISION XSUR6,H1,H2,H3,COEF
C
C-----------------------------------------------------------------------
C
      XSUR6  = XMUL/6.D0
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 3 IELEM = 1 , NELEM
C
         H1 = Z(IKLE4(IELEM)) - Z(IKLE1(IELEM))
         H2 = Z(IKLE5(IELEM)) - Z(IKLE2(IELEM))
         H3 = Z(IKLE6(IELEM)) - Z(IKLE3(IELEM))
C
         COEF = XSUR6 * SURFAC(IELEM)
C
         W1(IELEM) = COEF * H1
         W2(IELEM) = COEF * H2
         W3(IELEM) = COEF * H3
         W4(IELEM) = W1(IELEM)
         W5(IELEM) = W2(IELEM)
         W6(IELEM) = W3(IELEM)
C
3     CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
