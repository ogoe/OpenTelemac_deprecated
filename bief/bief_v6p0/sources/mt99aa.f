C                       *****************
                        SUBROUTINE MT99AA
C                       *****************
C
     *( A11 , A12 , A13 ,
     *  A21 , A22 , A23 ,
     *  A31 , A32 , A33 ,
     *  XMUL,SF,F,XEL,YEL,
     *  SURFAC,IKLE1,IKLE2,IKLE3,NELEM,NELMAX,FORMUL,TDIA,TEXT)
C
C***********************************************************************
C BIEF VERSION 5.1           28/11/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : 
C
C-----------------------------------------------------------------------
C
C      FONCTION:
C      ========:
C
C    CE SOUS-PROGRAMME CALCULE LES COEFFICIENTS DE LA MATRICE SUIVANTE:
C
C                                 /     DF             D
C                    A    = XMUL /  F * -- * PSI2(J) * --( PSI1(I) ) D(O
C                     I J       /S      DX             DX
C
C    PAR MAILLE ELEMENTAIRE.
C
C     J(X,Y) : JACOBIEN DE LA TRANSFORMATION ISOPARAMETRIQUE
C
C     L'ELEMENT EST LE TRIANGLE P1
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
C PROGRAMMES APPELES : AUCUN
C
C**********************************************************************
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_MT99AA => MT99AA
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX),IKLE3(NELMAX)
C
      CHARACTER(LEN=1), INTENT(INOUT) :: TDIA,TEXT
      CHARACTER(LEN=16), INTENT(IN) :: FORMUL
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A21(*),A22(*),A23(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A31(*),A32(*),A33(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: F(*)
C
C     STRUCTURE DE F
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF
C
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
      DOUBLE PRECISION, INTENT(IN) :: SURFAC(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      INTEGER IELEM
      DOUBLE PRECISION SUR48,DET,F1,F2,F3,X2,X3,Y2,Y3,F123,SUR24
C
C=======================================================================
C
      SUR24 = XMUL/24.D0
      SUR48 = XMUL/48.D0
C
C-----------------------------------------------------------------------
C
      IF(SF%ELM.NE.11) THEN
C
        IF (LNG.EQ.1) WRITE(LU,2000) SF%ELM
        IF (LNG.EQ.2) WRITE(LU,2001) SF%ELM
2000    FORMAT(1X,'MT99AA (BIEF) : TYPE DE F : ',I6,' NON PREVU')
2001    FORMAT(1X,'MT99AA (BIEF) : TYPE OF F:',I6,' NOT IMPLEMENTED')
        CALL PLANTE(0)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C        
      IF(FORMUL(8:16).EQ.'     0XX0') THEN
C
      TDIA='Q'
      TEXT='Q'
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 1 IELEM = 1 , NELEM
C
      DET = SUR48 / SURFAC(IELEM)
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      Y2 = YEL(IELEM,2)
      Y3 = YEL(IELEM,3)
C
      F1 = F(IKLE1(IELEM))
      F2 = F(IKLE2(IELEM))
      F3 = F(IKLE3(IELEM))
C
      A11(IELEM)=-(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(Y2-Y3)*(2*F1+F2+F3)*DET
      A12(IELEM)=-(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*Y3*(2*F1+F2+F3)*DET
      A13(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*Y2*(2*F1+F2+F3)*DET
      A21(IELEM)=-(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(Y2-Y3)*(F1+2*F2+F3)*DET
      A22(IELEM)=-(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*Y3*(F1+2*F2+F3)*DET
      A23(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*Y2*(F1+2*F2+F3)*DET
      A31(IELEM)=-(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(Y2-Y3)*(F1+F2+2*F3)*DET
      A32(IELEM)=-(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*Y3*(F1+F2+2*F3)*DET
      A33(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*Y2*(F1+F2+2*F3)*DET
C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
1     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSEIF(FORMUL(8:16).EQ.'     0YY0') THEN
C
      TDIA='Q'
      TEXT='Q'
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 2 IELEM = 1 , NELEM
C
      DET = SUR48 / SURFAC(IELEM)
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      X2 = XEL(IELEM,2)
      X3 = XEL(IELEM,3)
C
      F1 = F(IKLE1(IELEM))
      F2 = F(IKLE2(IELEM))
      F3 = F(IKLE3(IELEM))
C
      A11(IELEM)=-(X2*F1-X3*F1+X3*F2-X2*F3)*(-X2+X3)*(2*F1+F2+F3)*DET
      A12(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)*X3*(2*F1+F2+F3)*DET
      A13(IELEM)=-(X2*F1-X3*F1+X3*F2-X2*F3)*X2*(2*F1+F2+F3)*DET
      A21(IELEM)=-(X2*F1-X3*F1+X3*F2-X2*F3)*(-X2+X3)*(F1+2*F2+F3)*DET
      A22(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)*X3*(F1+2*F2+F3)*DET
      A23(IELEM)=-(X2*F1-X3*F1+X3*F2-X2*F3)*X2*(F1+2*F2+F3)*DET
      A31(IELEM)=-(X2*F1-X3*F1+X3*F2-X2*F3)*(-X2+X3)*(F1+F2+2*F3)*DET
      A32(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)*X3*(F1+F2+2*F3)*DET
      A33(IELEM)=-(X2*F1-X3*F1+X3*F2-X2*F3)*X2*(F1+F2+2*F3)*DET
C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
2     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSEIF(FORMUL(8:16).EQ.'     XX00') THEN
C
      TDIA='Q'
      TEXT='Q'
C     SYMETRIE NON PRISE EN COMPTE
C     TEXT='S'
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 3 IELEM = 1 , NELEM
C
      DET = SUR48 / SURFAC(IELEM)
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      Y2 = YEL(IELEM,2)
      Y3 = YEL(IELEM,3)
C
      F1 = F(IKLE1(IELEM))
      F2 = F(IKLE2(IELEM))
      F3 = F(IKLE3(IELEM))
      A11(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)**2*2*DET
      A12(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)**2*DET
      A13(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)**2*DET
      A21(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)**2*DET
      A22(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)**2*2*DET
      A23(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)**2*DET
      A31(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)**2*DET
      A32(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)**2*DET
      A33(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)**2*2*DET

C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
3     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSEIF(FORMUL(8:16).EQ.'     0X0Y') THEN
C
      TDIA='Q'
      TEXT='Q'
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 5 IELEM = 1 , NELEM
C
      DET = SUR48 / SURFAC(IELEM)
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      Y2 = YEL(IELEM,2)
      Y3 = YEL(IELEM,3)
      X2 = XEL(IELEM,2)
      X3 = XEL(IELEM,3)
C
      F1 = F(IKLE1(IELEM))
      F2 = F(IKLE2(IELEM))
      F3 = F(IKLE3(IELEM))
C
      A11(IELEM)=-(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(-X2+X3)*(2*F1+F2+F3)*DET
      A12(IELEM)=-(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(-X2+X3)*(F1+2*F2+F3)*DET
      A13(IELEM)=-(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(-X2+X3)*(F1+F2+2*F3)*DET
      A21(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*X3*(2*F1+F2+F3)*DET
      A22(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*X3*(F1+2*F2+F3)*DET
      A23(IELEM)=(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*X3*(F1+F2+2*F3)*DET
      A31(IELEM)=-(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*X2*(2*F1+F2+F3)*DET
      A32(IELEM)=-(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*X2*(F1+2*F2+F3)*DET
      A33(IELEM)=-(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*X2*(F1+F2+2*F3)*DET
C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
5     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSEIF(FORMUL(8:16).EQ.'     XY00') THEN
C
      TDIA='Q'
      TEXT='Q'
C     SYMETRIE NON PRISE EN COMPTE
C     TEXT='S'
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 6 IELEM = 1 , NELEM
C
      DET = SUR48 / SURFAC(IELEM)
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      Y2 = YEL(IELEM,2)
      Y3 = YEL(IELEM,3)
      X2 = XEL(IELEM,2)
      X3 = XEL(IELEM,3)
C
      F1 = F(IKLE1(IELEM))
      F2 = F(IKLE2(IELEM))
      F3 = F(IKLE3(IELEM))
C
      A11(IELEM)=
     *(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(X2*F1-X3*F1+X3*F2-X2*F3)*2*DET
      A12(IELEM)=
     *(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(X2*F1-X3*F1+X3*F2-X2*F3)*DET
      A13(IELEM)=
     *(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(X2*F1-X3*F1+X3*F2-X2*F3)*DET
      A21(IELEM)=
     *(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(X2*F1-X3*F1+X3*F2-X2*F3)*DET
      A22(IELEM)=
     *(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(X2*F1-X3*F1+X3*F2-X2*F3)*2*DET
      A23(IELEM)=
     *(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(X2*F1-X3*F1+X3*F2-X2*F3)*DET
      A31(IELEM)=
     *(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(X2*F1-X3*F1+X3*F2-X2*F3)*DET
      A32(IELEM)=
     *(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(X2*F1-X3*F1+X3*F2-X2*F3)*DET
      A33(IELEM)=
     *(-Y2*F1+Y3*F1-Y3*F2+Y2*F3)*(X2*F1-X3*F1+X3*F2-X2*F3)*2*DET
C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
6     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSEIF(FORMUL(8:16).EQ.'     YY00') THEN
C
      TDIA='Q'
      TEXT='Q'
C     SYMETRIE NON PRISE EN COMPTE
C     TEXT='S'
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 7 IELEM = 1 , NELEM
C
      DET = SUR48 / SURFAC(IELEM)
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      X2 = XEL(IELEM,2)
      X3 = XEL(IELEM,3)
C
      F1 = F(IKLE1(IELEM))
      F2 = F(IKLE2(IELEM))
      F3 = F(IKLE3(IELEM))
C
      A11(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)**2*2*DET
      A12(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)**2*DET
      A13(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)**2*DET
      A21(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)**2*DET
      A22(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)**2*2*DET
      A23(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)**2*DET
      A31(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)**2*DET
      A32(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)**2*DET
      A33(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)**2*2*DET
C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
7     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSEIF(FORMUL(8:16).EQ.'     0Y0X') THEN
C
      TDIA='Q'
      TEXT='Q'
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 9 IELEM = 1 , NELEM
C
      DET = SUR48 / SURFAC(IELEM)
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      Y2 = YEL(IELEM,2)
      Y3 = YEL(IELEM,3)
      X2 = XEL(IELEM,2)
      X3 = XEL(IELEM,3)
C
      F1 = F(IKLE1(IELEM))
      F2 = F(IKLE2(IELEM))
      F3 = F(IKLE3(IELEM))
C
      A11(IELEM)=-(X2*F1-X3*F1+X3*F2-X2*F3)*(Y2-Y3)*(2*F1+F2+F3)*DET
      A12(IELEM)=-(X2*F1-X3*F1+X3*F2-X2*F3)*(Y2-Y3)*(F1+2*F2+F3)*DET
      A13(IELEM)=-(X2*F1-X3*F1+X3*F2-X2*F3)*(Y2-Y3)*(F1+F2+2*F3)*DET
      A21(IELEM)=-(X2*F1-X3*F1+X3*F2-X2*F3)*Y3*(2*F1+F2+F3)*DET
      A22(IELEM)=-(X2*F1-X3*F1+X3*F2-X2*F3)*Y3*(F1+2*F2+F3)*DET
      A23(IELEM)=-(X2*F1-X3*F1+X3*F2-X2*F3)*Y3*(F1+F2+2*F3)*DET
      A31(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)*Y2*(2*F1+F2+F3)*DET
      A32(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)*Y2*(F1+2*F2+F3)*DET
      A33(IELEM)=(X2*F1-X3*F1+X3*F2-X2*F3)*Y2*(F1+F2+2*F3)*DET
C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
9     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSEIF(FORMUL(8:16).EQ.'00XX+00YY') THEN
C
      TDIA='Q'
      TEXT='Q'
C     SYMETRIE NON PRISE EN COMPTE
C     TEXT='S'
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 10 IELEM = 1 , NELEM
C
      DET = SUR24 / SURFAC(IELEM)
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      X2 = XEL(IELEM,2)
      X3 = XEL(IELEM,3)
      Y2 = YEL(IELEM,2)
      Y3 = YEL(IELEM,3)
C
      F1 = F(IKLE1(IELEM))
      F2 = F(IKLE2(IELEM))
      F3 = F(IKLE3(IELEM))
      F123 = F2**2+F1*F2+F2*F3+F3**2+F1**2+F1*F3
C
C  ELEMENTS EXTERIEURS A LA DIAGONALE
C
      A12(IELEM) = (Y2*Y3-Y3**2+X2*X3-X3**2)*F123*DET
      A13(IELEM) = -(Y2**2-Y2*Y3+X2**2-X2*X3)*F123*DET
      A23(IELEM) = -(Y2*Y3+X2*X3)*F123*DET
C
C  TERMES DIAGONAUX
C
      A11(IELEM) = (Y2**2-2*Y2*Y3+Y3**2+X2**2-2*X2*X3+X3**2)*F123*DET
      A22(IELEM) = (Y3**2+X3**2)*F123*DET
      A33(IELEM) = (Y2**2+X2**2)*F123*DET
C
C  SYMETRIES
C
      A21(IELEM) = A12(IELEM)
      A31(IELEM) = A13(IELEM)
      A32(IELEM) = A23(IELEM)
C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSEIF(FORMUL(8:16).EQ.'     00XX') THEN
C
      TDIA='Q'
      TEXT='Q'
C     SYMETRIE NON PRISE EN COMPTE
C     TEXT='S'
C
C   BOUCLE SUR LES ELEMENTS
C
      DO IELEM = 1 , NELEM
C
      DET = SUR24 / SURFAC(IELEM)
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      Y2 = YEL(IELEM,2)
      Y3 = YEL(IELEM,3)
C
      F1 = F(IKLE1(IELEM))
      F2 = F(IKLE2(IELEM))
      F3 = F(IKLE3(IELEM))
      F123 = F2**2+F1*F2+F2*F3+F3**2+F1**2+F1*F3
C
C  ELEMENTS EXTERIEURS A LA DIAGONALE
C
      A12(IELEM) = (Y2*Y3-Y3**2)*F123*DET
      A13(IELEM) = -(Y2**2-Y2*Y3)*F123*DET
      A23(IELEM) = -(Y2*Y3)*F123*DET
C
C  TERMES DIAGONAUX
C
      A11(IELEM) = (Y2**2-2*Y2*Y3+Y3**2)*F123*DET
      A22(IELEM) = (Y3**2)*F123*DET
      A33(IELEM) = (Y2**2)*F123*DET
C
C  SYMETRIES
C
      A21(IELEM) = A12(IELEM)
      A31(IELEM) = A13(IELEM)
      A32(IELEM) = A23(IELEM)
C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      ELSEIF(FORMUL(8:16).EQ.'     00YY') THEN
C
      TDIA='Q'
      TEXT='Q'
C     SYMETRIE NON PRISE EN COMPTE
C     TEXT='S'
C
C   BOUCLE SUR LES ELEMENTS
C
      DO IELEM = 1 , NELEM
C
      DET = SUR24 / SURFAC(IELEM)
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      X2 = XEL(IELEM,2)
      X3 = XEL(IELEM,3)
C
      F1 = F(IKLE1(IELEM))
      F2 = F(IKLE2(IELEM))
      F3 = F(IKLE3(IELEM))
      F123 = F2**2+F1*F2+F2*F3+F3**2+F1**2+F1*F3
C
C  ELEMENTS EXTERIEURS A LA DIAGONALE
C
      A12(IELEM) = (X2*X3-X3**2)*F123*DET
      A13(IELEM) = -(X2**2-X2*X3)*F123*DET
      A23(IELEM) = -(X2*X3)*F123*DET
C
C  TERMES DIAGONAUX
C
      A11(IELEM) = (X2**2-2*X2*X3+X3**2)*F123*DET
      A22(IELEM) = (X3**2)*F123*DET
      A33(IELEM) = (X2**2)*F123*DET
C
C  SYMETRIES
C
      A21(IELEM) = A12(IELEM)
      A31(IELEM) = A13(IELEM)
      A32(IELEM) = A23(IELEM)
C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      ELSEIF(FORMUL(8:16).EQ.'     00XY') THEN
C
      TDIA='Q'
      TEXT='Q'
C
C   BOUCLE SUR LES ELEMENTS
C
      DO IELEM = 1 , NELEM
C
      DET = SUR24 / SURFAC(IELEM)
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      X2 = XEL(IELEM,2)
      X3 = XEL(IELEM,3)
      Y2 = YEL(IELEM,2)
      Y3 = YEL(IELEM,3)
C
      F1 = F(IKLE1(IELEM))
      F2 = F(IKLE2(IELEM))
      F3 = F(IKLE3(IELEM))
      F123 = F2**2+F1*F2+F2*F3+F3**2+F1**2+F1*F3
C
      A11(IELEM) =  (Y2-Y3)*(-X2+X3)*F123*DET
      A12(IELEM) =      Y3 *(-X2+X3)*F123*DET
      A13(IELEM) =  -Y2    *(-X2+X3)*F123*DET
      A21(IELEM) = -(Y2-Y3)*     X3 *F123*DET
      A22(IELEM) =     -Y3 *     X3 *F123*DET
      A23(IELEM) =   Y2    *     X3 *F123*DET
      A31(IELEM) =  (Y2-Y3)*  X2    *F123*DET
      A32(IELEM) =      Y3 *  X2    *F123*DET
      A33(IELEM) =  -Y2    *  X2    *F123*DET
C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
      ENDDO
C
C   CAS NON PREVU
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
        IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
1000    FORMAT(1X,'MT99AA (BIEF) : MATRICE NON PREVUE : ',A16)
1001    FORMAT(1X,'MT99AA (BIEF) : MATRIX NOT IMPLEMENTED:',A16)
        CALL PLANTE(0)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
