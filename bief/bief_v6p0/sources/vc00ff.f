C                       *****************
                        SUBROUTINE VC00FF
C                       *****************
C
     *( XMUL,X,Y,Z,
     *  IKLE1,IKLE2,IKLE3,IKLE4,NBOR,NELEM,NELMAX,W1,W2,W3,W4 )
C
C***********************************************************************
C BIEF VERSION 5.4           05/02/91    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C                    /
C    VEC(I) = XMUL  /    PSI(I)  D(OMEGA)
C                  /OMEGA
C
C
C    PSI(I) EST UNE BASE DE TYPE SEGMENT P1
C
C    F EST UN VECTEUR DE TYPE IELMF
C
C    ATTENTION : LE RESULTAT EST DANS W SOUS FORME NON ASSEMBLEE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      XMUL      | -->|  COEFFICIENT MULTIPLICATEUR.
C |      XEL,YEL,..| -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |      SURFAC    | -->|  SURFACE DES ELEMENTS.
C |      IKLE1,... | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE.
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
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: NBOR(*)
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX)
      INTEGER, INTENT(IN) :: IKLE3(NELMAX),IKLE4(NELMAX)
C
      DOUBLE PRECISION, INTENT(IN)    :: X(*),Y(*),Z(*)
      DOUBLE PRECISION, INTENT(INOUT) :: W1(NELMAX),W2(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W3(NELMAX),W4(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER IELEM,I1,I2,I3,I4
      DOUBLE PRECISION XSUR24,H1,H2,AL
C
      INTRINSIC SQRT
C
C***********************************************************************
C
C NOTE ON PARALLELISM : WITH PARALLELISM NELEM MAY BE WRONG AS BOUNDARY
C                       ELEMENTS MAY BE IN ANOTHER SUBDOMAIN. IN THIS
C                       CASE WE HAVE I1=I2 AND THUS AL=0.D0, SO WRONG
C                       ELEMENTS DO NOT CONTRIBUTE.
C
C***********************************************************************
C
      XSUR24 = XMUL/24.D0
C
C   BOUCLE SUR LES FACES DE BORD
C
      DO 1 IELEM = 1,NELEM
C
C  NUMEROTATION GLOBALE DES SOMMETS DE LA FACE
C
         I1 = IKLE1(IELEM)
         I2 = IKLE2(IELEM)
         I3 = IKLE3(IELEM)
         I4 = IKLE4(IELEM)
C
         AL = SQRT((X(NBOR(I2))-X(NBOR(I1)))**2
     *            +(Y(NBOR(I2))-Y(NBOR(I1)))**2) * XSUR24
C
         H1 = Z(NBOR(I4)) - Z(NBOR(I1))
         H2 = Z(NBOR(I3)) - Z(NBOR(I2))
C
         W1(IELEM) = (3.D0*H1+H2)*AL
         W2(IELEM) = (3.D0*H2+H1)*AL
         W3(IELEM) = (3.D0*H2+H1)*AL
         W4(IELEM) = (3.D0*H1+H2)*AL
C
1     CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
