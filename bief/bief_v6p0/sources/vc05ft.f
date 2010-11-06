C                       *****************
                        SUBROUTINE VC05FT
C                       *****************
C
     *( XMUL,SU,SV,U,V,X,Y,Z,
     *  IKLE1,IKLE2,IKLE3,NBOR,NELEM,NELMAX,W1,W2,W3)
C
C***********************************************************************
C BIEF VERSION 6.0           05/02/91    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C                    /         ->
C    VEC(I) = XMUL  /    (U,V).N  PSI(I) D(GAMMA)
C                  /GAMMA
C
C
C    PSI(I) EST UNE BASE DE TYPE TRIANGLE P1 SUR UN MAILLAGE VERTICAL
C    LATERAL D'UN MAILLAGE DE PRISMES DECOUPES EN TETRAEDRES
C
C    ATTENTION : LE RESULTAT EST DANS W1,2,3 SOUS FORME NON ASSEMBLEE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      XMUL      | -->|  COEFFICIENT MULTIPLICATEUR.
C |      SF,SG,SH  | -->|  STRUCTURES DES FONCTIONS F,G ET H
C |      SU,SV,SW  | -->|  STRUCTURES DES FONCTIONS U,V ET W
C |      F,G,H     | -->|  FONCTIONS INTERVENANT DANS LA FORMULE.
C |      U,V,W     | -->|  COMPOSANTES D'UN VECTEUR
C |                |    |  INTERVENANT DANS LA FORMULE.
C |      X,Y,Z     | -->|  COORDONNEES DES POINTS
C |      SURFAC    | -->|  SURFACE DES ELEMENTS.
C |      IKLE1,... | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE.
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE.
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE.
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      W1,2,3,4  |<-- |  VECTEUR RESULTAT SOUS FORME NON ASSEMBLEE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C  PROGRAMMES APPELES :
C
C-----------------------------------------------------------------------
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_VC05FT => VC05FT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NBOR(*),NELEM,NELMAX 
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX),IKLE3(NELMAX)
C
      DOUBLE PRECISION, INTENT(IN)   :: X(*),Y(*),Z(*)
      DOUBLE PRECISION, INTENT(INOUT):: W1(NELMAX),W2(NELMAX),W3(NELMAX)
      DOUBLE PRECISION, INTENT(IN)   :: XMUL
C
C-----------------------------------------------------------------------
C
C     STRUCTURES DE U,V ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN)   :: SU,SV
      DOUBLE PRECISION, INTENT(IN) :: U(*),V(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION U1,U2,U3,V1,V2,V3,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3
      DOUBLE PRECISION XSUR24
      INTEGER IELMU,IELMV,IELEM,N1,N2,N3,I1,I2,I3
C
C-----------------------------------------------------------------------
C
      IELMU=SU%ELM
      IELMV=SV%ELM
C
C-----------------------------------------------------------------------
C
      XSUR24 = XMUL/24.D0
C
C     U LINEAIRE PAR TETRAEDRES
C
C-----------------------------------------------------------------------
C
      IF(IELMU.EQ.61.AND.IELMV.EQ.61) THEN
C
C-----------------------------------------------------------------------
C
C   BOUCLE SUR LES FACES DE BORD
C
         DO IELEM = 1,NELEM
C
C           NUMEROTATION DE BORD DES SOMMETS DE LA FACE
C
            I1 = IKLE1(IELEM)
            I2 = IKLE2(IELEM)
            I3 = IKLE3(IELEM)
C
C           NUMEROTATION GLOBALE DES SOMMETS DE LA FACE
C
            N1 = NBOR(I1)
            N2 = NBOR(I2)
            N3 = NBOR(I3)
C
            U1 = U(I1)
            U2 = U(I2)
            U3 = U(I3)
            V1 = V(I1)
            V2 = V(I2)
            V3 = V(I3)
C
            X1 = X(N1)
            X2 = X(N2)-X1
            X3 = X(N3)-X1
            Y1 = Y(N1)
            Y2 = Y(N2)-Y1
            Y3 = Y(N3)-Y1
            Z1 = Z(N1)
            Z2 = Z(N2)-Z1
            Z3 = Z(N3)-Z1
C
C           CALCUL DE LA SURFACE DU TRIANGLE (PAR PRODUIT VECTORIEL)
C
C           REMARQUE : VECTEUR NORMAL AU TRIANGLE,
C                      DONT LA NORME EST LA SURFACE :
C
C                      0.5  (Y2*Z3-Y3*Z2)
C                      0.5  (X3*Z2-X2*Z3)
C                      0.5  (X2*Y3-X3*Y2)  : CE TERME EST NUL
C
C           ON S'INSPIRE ENSUITE DE MASVEC SUR DES TRIANGLES
C
C
            W1(IELEM) = XSUR24* (  (Y2*Z3-Y3*Z2)*(2*U1+  U2+  U3)
     *                            +(X3*Z2-X2*Z3)*(2*V1+  V2+  V3) )
            W2(IELEM) = XSUR24* (  (Y2*Z3-Y3*Z2)*(  U1+2*U2+  U3)
     *                            +(X3*Z2-X2*Z3)*(  V1+2*V2+  V3) )
            W3(IELEM) = XSUR24* (  (Y2*Z3-Y3*Z2)*(  U1+  U2+2*U3)
     *                            +(X3*Z2-X2*Z3)*(  V1+  V2+2*V3) )
C
         ENDDO
C
C-----------------------------------------------------------------------
C
      ELSEIF(IELMU.EQ.51.AND.IELMV.EQ.51) THEN
C
C-----------------------------------------------------------------------
C
C   BOUCLE SUR LES FACES DE BORD
C
         DO IELEM = 1,NELEM
C
C           NUMEROTATION GLOBALE DES SOMMETS DE LA FACE
C
            N1 = NBOR(IKLE1(IELEM))
            N2 = NBOR(IKLE2(IELEM))
            N3 = NBOR(IKLE3(IELEM))
C
            U1 = U(N1)
            U2 = U(N2)
            U3 = U(N3)
            V1 = V(N1)
            V2 = V(N2)
            V3 = V(N3)
            X1 = X(N1)
            X2 = X(N2)-X1
            X3 = X(N3)-X1
            Y1 = Y(N1)
            Y2 = Y(N2)-Y1
            Y3 = Y(N3)-Y1
            Z1 = Z(N1)
            Z2 = Z(N2)-Z1
            Z3 = Z(N3)-Z1
C
C           CALCUL DE LA SURFACE DU TRIANGLE (PAR PRODUIT VECTORIEL)
C
C           REMARQUE : VECTEUR NORMAL AU TRIANGLE,
C                      DONT LA NORME EST LA SURFACE :
C
C                      0.5  (Y2*Z3-Y3*Z2)
C                      0.5  (X3*Z2-X2*Z3)
C                      0.5  (X2*Y3-X3*Y2)  : CE TERME EST NUL
C
C           ON S'INSPIRE ENSUITE DE MASVEC SUR DES TRIANGLES
C
C
            W1(IELEM) = XSUR24* (  (Y2*Z3-Y3*Z2)*(2*U1+  U2+  U3)
     *                            +(X3*Z2-X2*Z3)*(2*V1+  V2+  V3) )
            W2(IELEM) = XSUR24* (  (Y2*Z3-Y3*Z2)*(  U1+2*U2+  U3)
     *                            +(X3*Z2-X2*Z3)*(  V1+2*V2+  V3) )
            W3(IELEM) = XSUR24* (  (Y2*Z3-Y3*Z2)*(  U1+  U2+2*U3)
     *                            +(X3*Z2-X2*Z3)*(  V1+  V2+2*V3) )
C
         ENDDO
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
         IF (LNG.EQ.1) WRITE(LU,100) IELMU,SU%NAME
         IF (LNG.EQ.2) WRITE(LU,101) IELMU,SU%NAME
100      FORMAT(1X,'VC05FT (BIEF) :',/,
     *          1X,'DISCRETISATION DE U NON PREVUE : ',1I6,
     *          1X,'NOM REEL : ',A6)
101      FORMAT(1X,'VC05FT (BIEF) :',/,
     *          1X,'DISCRETIZATION OF U NOT AVAILABLE:',1I6,
     *          1X,'REAL NAME: ',A6)
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
