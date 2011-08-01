C                       *****************
                        SUBROUTINE VC00AA
C                       *****************
C
     *( XMUL,SURFAC,NELEM,NELMAX,W1,W2,W3 )
C
C***********************************************************************
C BIEF VERSION 5.1           05/02/91    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C                    /
C    VEC(I) = XMUL  /    PSI(I)  D(OMEGA)
C                  /OMEGA
C
C    PSI(I) EST UNE BASE DE TYPE TRIANGLE P1
C
C    ATTENTION : LE RESULTAT EST DANS W  SOUS FORME NON ASSEMBLEE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      XMUL      | -->|  COEFFICIENT MULTIPLICATEUR.
C |      SURFAC    | -->|  SURFACE DES ELEMENTS.
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE.
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE.
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      W1,2,3    | -->|  VECTEUR RESULTAT SOUS FORME NON ASSEMBLEE.
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
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
C
      DOUBLE PRECISION, INTENT(INOUT) :: W1(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W2(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W3(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: SURFAC(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER IELEM
C
      DOUBLE PRECISION XSUR3,COEF
C
C-----------------------------------------------------------------------
C
      XSUR3 = XMUL / 3.D0
C
C-----------------------------------------------------------------------
C
      DO 3 IELEM = 1 , NELEM
C
        COEF = XSUR3 * SURFAC(IELEM)
C
        W1(IELEM) = COEF
        W2(IELEM) = COEF
        W3(IELEM) = COEF
C
3     CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
