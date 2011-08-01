C                       *****************
                        SUBROUTINE MT06FT2
C                       *****************
C
     *( A11 , A12 , A13 ,
     *        A22 , A23 , 
     *              A33 , 
     *  XMUL,SF,F,SG,G,X,Y,Z,IKLE1,IKLE2,IKLE3,NBOR,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.5         26/04/04    J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
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
C     PAR MAILLE ELEMENTAIRE .
C
C     J(X,Y) : JACOBIEN DE LA TRANSFORMATION ISOPARAMETRIQUE
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     L'ELEMENT EST LE TRIANGLE P1, MAIS DANS UN MAILLAGE DE PRISMES
C     DECOUPE EN TETRAEDRES
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_MT06FT2 => MT06FT2
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
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*)
      DOUBLE PRECISION, INTENT(INOUT) ::        A22(*),A23(*)
      DOUBLE PRECISION, INTENT(INOUT) ::               A33(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: F(*),G(*)
C
C     STRUCTURE DE F
      TYPE(BIEF_OBJ), INTENT(IN) :: SF,SG
C
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*),Z(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTRINSIC SQRT
C
C-----------------------------------------------------------------------
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      INTEGER IELMF,IELMG,I1,I2,I3,IELEM
C
      DOUBLE PRECISION SUR60,S,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,F1,F2,F3,F123
      DOUBLE PRECISION DET1,DET2
C
C**********************************************************************
C
      IELMF=SF%ELM
      IELMG=SG%ELM
C
C-----------------------------------------------------------------------
C
C     F LINEAIRE PAR FACE DE BORD
C
      IF( (IELMF.EQ.61.OR.IELMF.EQ.81) .AND. IELMG.EQ.10 ) THEN
C
         SUR60  = XMUL/60.D0
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
            F1 = F(IKLE1(IELEM)) * G(IELEM)
            F2 = F(IKLE2(IELEM)) * G(IELEM)
            F3 = F(IKLE3(IELEM)) * G(IELEM)
            F123  = F1 + F2 + F3
C
C           CALCUL DE LA SURFACE DU TRIANGLE (PAR PRODUIT VECTORIEL)
C
            S=0.5D0*SQRT(  (Y2*Z3-Y3*Z2)**2 
     *                    +(X3*Z2-X2*Z3)**2  
     *                    +(X2*Y3-X3*Y2)**2  )
C 
            DET1 = S * SUR60
            DET2 = DET1 + DET1
C
C***********************************************************************
C
C          ELEMENTS EXTERIEURS A LA DIAGONALE
C
           A12(IELEM) = DET1 * (F123+F123-F3)
           A13(IELEM) = DET1 * (F123+F123-F2)
           A23(IELEM) = DET1 * (F123+F123-F1)
C
C          TERMES DIAGONAUX
C
           A11(IELEM) = DET2 * (F123+F1+F1)
           A22(IELEM) = DET2 * (F123+F2+F2)
           A33(IELEM) = DET2 * (F123+F3+F3)
C
1        CONTINUE
C
C-----------------------------------------------------------------------
C
C     AUTRES TYPES DE DISCRETISATION DE F
C
      ELSE
C
         IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME,SG%NAME
         IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME,SG%NAME
100      FORMAT(1X,'MT06FT2 (BIEF) :',/,
     *          1X,'DISCRETISATION DE F NON PREVUE : ',1I6,
     *          1X,'NOM REELS : ',A6,' ET ',A6)
101      FORMAT(1X,'MT06FT2 (BIEF) :',/,
     *          1X,'DISCRETIZATION OF F NOT AVAILABLE:',1I6,
     *          1X,'REAL NAME: ',A6,' AND ',A6)
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C     NOTE : SUR UN MAILLAGE DE TRIANGLE EN PLAN (X,Y)
C
C     DO 1 IELEM = 1 , NELEM
C
C     F1 = F(IKLE1(IELEM))
C     F2 = F(IKLE2(IELEM))
C     F3 = F(IKLE3(IELEM))
C
C     F123 = F1 + F2 + F3
C
C     DET1 = SURFAC(IELEM) * SUR60
C     DET2 = DET1 + DET1
C
C***********************************************************************
C
C  ELEMENTS EXTERIEURS A LA DIAGONALE
C
C     A12(IELEM) = DET1 * (F123+F123-F3)
C     A13(IELEM) = DET1 * (F123+F123-F2)
C     A23(IELEM) = DET1 * (F123+F123-F1)
C
C  TERMES DIAGONAUX
C
C     A11(IELEM) = DET2 * (F123+F1+F1)
C     A22(IELEM) = DET2 * (F123+F2+F2)
C     A33(IELEM) = DET2 * (F123+F3+F3)
C
C1     CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
