C                       *****************
                        SUBROUTINE VC00FT
C                       *****************
C
     *( XMUL,X,Y,Z,
     *  IKLE1,IKLE2,IKLE3,NBOR,NELEM,NELMAX,W1,W2,W3)
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
C    PSI(I) EST UNE BASE DE TYPE SEGMENT P1 MAIS SUR UN MAILLAGE
C    DE TRIANGLES VERTICAUX DANS L'ESPACE X,Y,Z
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
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX),IKLE3(NELMAX)
C
      DOUBLE PRECISION, INTENT(IN)    :: X(*),Y(*),Z(*)
      DOUBLE PRECISION, INTENT(INOUT) :: W1(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W2(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W3(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER IELEM,I1,I2,I3
      DOUBLE PRECISION XSUR3,COEF,X1,X2,X3,Y1,Y2
      DOUBLE PRECISION Y3,Z1,Z2,Z3,S
C
      INTRINSIC SQRT
C
C***********************************************************************
C
         XSUR3 = XMUL/3.D0
C
C   BOUCLE SUR LES FACES DE BORD
C
         DO 1 IELEM = 1,NELEM
C
C           NUMEROTATION GLOBALE DES SOMMETS DE LA FACE
C
            I1 = NBOR(IKLE1(IELEM))
            I2 = NBOR(IKLE2(IELEM))
            I3 = NBOR(IKLE3(IELEM))
C
            X1 = X(I1)
            Y1 = Y(I1)
            Z1 = Z(I1)
C
            X2 = X(I2)-X1
            X3 = X(I3)-X1            
            Y2 = Y(I2)-Y1
            Y3 = Y(I3)-Y1            
            Z2 = Z(I2)-Z1
            Z3 = Z(I3)-Z1
C
C           CALCUL DE LA SURFACE DU TRIANGLE (PAR PRODUIT VECTORIEL)
C
            S=0.5D0*SQRT(  (Y2*Z3-Y3*Z2)**2
     *                    +(X3*Z2-X2*Z3)**2  )
C    *                    +(X2*Y3-X3*Y2)**2  )  CE TERME EST NUL
C
            COEF=XSUR3*S
C
            W1(IELEM) = COEF
            W2(IELEM) = COEF
            W3(IELEM) = COEF
C
1        CONTINUE
C
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C
C     NOTE : POUR UN TRIANGLE PLAN (VC00AA):
C
C     XSUR3 = XMUL / 3.D0
C
C     DO 3 IELEM = 1 , NELEM
C
C
C       COEF = XSUR3 * SURFAC(IELEM)
C
C       W1(IELEM) = COEF 
C       W2(IELEM) = COEF 
C       W3(IELEM) = COEF 
C
C3     CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
