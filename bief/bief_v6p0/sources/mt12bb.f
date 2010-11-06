C                       *****************
                        SUBROUTINE MT12BB
C                       *****************
C
     *(  A11 , A12 , A13 , A14 ,
     *   A21 , A22 , A23 , A24 ,
     *   A31 , A32 , A33 , A34 ,
     *   A41 , A42 , A43 , A44 ,
     *   XMUL,SF,SU,SV,F,U,V,
     *   XEL,YEL,
     *   IKLE1,IKLE2,IKLE3,IKLE4,
     *   NELEM,NELMAX,ICOORD)
C
C***********************************************************************
C BIEF VERSION 5.1           09/12/94    J-M HERVOUET (LNH) 30 87 80 18
C                                          C MOULIN   (LNH) 30 87 83 81
C***********************************************************************
C
C  FONCTION : CALCUL DES COEFFICIENTS DE LA MATRICE SUIVANTE
C
C  EXEMPLE AVEC ICOORD = 1
C
C
C                /           D       -> -->
C A(I,J)= XMUL  /  PSI2(J) * --(F) * U.GRAD(PSI1(I))  D(OMEGA)
C              /OMEGA        DX
C
C
C  AVEC ICOORD=2 ON AURAIT UNE DERIVEE SUIVANT Y
C
C  PSI1 : BASES DE TYPE TRIANGLE P1
C  PSI2 : BASES DE TYPE IELM2
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
      USE BIEF, EX_MT12BB => MT12BB
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX,ICOORD
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX)
      INTEGER, INTENT(IN) :: IKLE3(NELMAX),IKLE4(NELMAX)
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*),A14(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A21(*),A22(*),A23(*),A24(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A31(*),A32(*),A33(*),A34(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A41(*),A42(*),A43(*),A44(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: F(*),U(*),V(*)
C
C     STRUCTURES DE F,U,V
      TYPE(BIEF_OBJ), INTENT(IN) :: SF,SU,SV
C
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF,IELMU,IELMV
      DOUBLE PRECISION X2,X3,Y2,Y3,F1,F2,F3,F4
      DOUBLE PRECISION U1,U2,U3,U4,V1,V2,V3,V4,UX,UY
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
      IELMU=SU%ELM
      IELMV=SV%ELM
C
C-----------------------------------------------------------------------
C  CAS OU F EST DE DISCRETISATION P1
C-----------------------------------------------------------------------
C
      IF(IELMF.EQ.11) THEN
C
      IF(IELMU.EQ.10.AND.IELMV.EQ.10) THEN
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
        X2 = XEL(IELEM,2)
        X3 = XEL(IELEM,3)
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM)) - F1
        F3  =  F(IKLE3(IELEM)) - F1
C
        UX  =  U(IELEM)
        UY  =  V(IELEM)
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM) = ((UX*Y3-2*UX*Y2+2*UY*X2-UY*X3)*(F3*Y2-F2*Y3))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A13(IELEM) = ((2*UX*Y3-UX*Y2+UY*X2-2*UY*X3)*(F3*Y2-F2*Y3))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A14(IELEM) = ((UX*Y3-UX*Y2+UY*X2-UY*X3)*(F3*Y2-F2*Y3))
     * *XMUL/(6*(X2* Y3-X3*Y2))
C
      A21(IELEM) = (-(UX*Y3+UX*Y2-UY*X2-UY*X3)*(F3*Y2-F2*Y3))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A23(IELEM) = (-(2*UX*Y3-UX*Y2+UY*X2-2*UY*X3)*(F3*Y2-F2*Y3))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A24(IELEM) = (-(UX*Y3-UY*X3)*(F3*Y2-F2*Y3))
     * *XMUL/(6*(X2*Y3-X3*Y2))
C
      A31(IELEM) = ((UX*Y3+UX*Y2-UY*X2-UY*X3)*(F3*Y2-F2*Y3))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A32(IELEM) = (-(UX*Y3-2*UX*Y2+2*UY*X2-UY*X3)*(F3*Y2-F2*Y3))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A34(IELEM) = ((UX*Y2-UY*X2)*(F3*Y2-F2*Y3))
     * *XMUL/(6*(X2*Y3-X3*Y2))
C
      A41(IELEM) = (-(UX*Y3-UX*Y2+UY*X2-UY*X3)*(F3*Y2-F2*Y3))
     * *XMUL/(6*(X2*Y3-X3*Y2))
C
      A42(IELEM) = ((UX*Y3-UY*X3)*(F3*Y2-F2*Y3))
     * *XMUL/(6*(X2*Y3-X3*Y2))
C
      A43(IELEM) = (-(UX*Y2-UY*X2)*(F3*Y2-F2*Y3))
     * *XMUL/(6*(X2*Y3-X3*Y2))
C
C   TERMES DIAGONAUX
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
C
        A11(IELEM) = - A21(IELEM)  - A31(IELEM) - A41(IELEM)
        A22(IELEM) = - A12(IELEM)  - A32(IELEM) - A42(IELEM)
        A33(IELEM) = - A13(IELEM)  - A23(IELEM) - A43(IELEM)
        A44(IELEM) = - A14(IELEM)  - A24(IELEM) - A34(IELEM)
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
        Y2  =  YEL(IELEM,2)
        Y3  =  YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM)) - F1
        F3  =  F(IKLE3(IELEM)) - F1
C
        UX  =  U(IELEM)
        UY  =  V(IELEM)
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM) = (-(UX*Y3-2*UX*Y2+2*UY*X2-UY*X3)*(X2*F3-X3*F2))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A13(IELEM) = (-(2*UX*Y3-UX*Y2+UY*X2-2*UY*X3)*(X2*F3-X3*F2))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A14(IELEM) = (-(UX*Y3-UX*Y2+UY*X2-UY*X3)*(X2*F3-X3*F2))
     * *XMUL/(6*(X2*Y3-X3*Y2))
C
      A21(IELEM) = ((UX*Y3+UX*Y2-UY*X2-UY*X3)*(X2*F3-X3*F2))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A23(IELEM) = ((2*UX*Y3-UX*Y2+UY*X2-2*UY*X3)*(X2*F3-X3*F2))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A24(IELEM) = ((UX*Y3-UY*X3)*(X2*F3-X3*F2))
     * *XMUL/(6*(X2*Y3-X3*Y2))
C
      A31(IELEM) = (-(UX*Y3+UX*Y2-UY*X2-UY*X3)*(X2*F3-X3*F2))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A32(IELEM) = ((UX*Y3-2*UX*Y2+2*UY*X2-UY*X3)*(X2*F3-X3*F2))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A34(IELEM) = (-(UX*Y2-UY*X2)*(X2*F3-X3*F2))
     * *XMUL/(6*(X2*Y3-X3*Y2))
C
      A41(IELEM) = ((UX*Y3-UX*Y2+UY*X2-UY*X3)*(X2*F3-X3*F2))
     * *XMUL/(6*(X2*Y3-X3*Y2))
C
      A42(IELEM) = (-(UX*Y3-UY*X3)*(X2*F3-X3*F2))
     * *XMUL/(6*(X2*Y3-X3*Y2))
C
      A43(IELEM) = ((UX*Y2-UY*X2)*(X2*F3-X3*F2))
     * *XMUL/(6*(X2*Y3-X3*Y2))
C
C   TERMES DIAGONAUX
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
C
        A11(IELEM) = - A21(IELEM)  - A31(IELEM) - A41(IELEM)
        A22(IELEM) = - A12(IELEM)  - A32(IELEM) - A42(IELEM)
        A33(IELEM) = - A13(IELEM)  - A23(IELEM) - A43(IELEM)
        A44(IELEM) = - A14(IELEM)  - A24(IELEM) - A34(IELEM)
C
2       CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
          CALL PLANTE(0)
          STOP
        ENDIF
C
C
      ELSEIF(IELMU.EQ.12) THEN
C
C================================
C  CAS DE LA DERIVEE SUIVANT X  =
C================================
C
        IF(ICOORD.EQ.1) THEN
C
C   BOUCLE SUR LES ELEMENTS
C
        DO 3 IELEM = 1 , NELEM
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
        X2 = XEL(IELEM,2)
        X3 = XEL(IELEM,3)
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM)) - F1
        F3  =  F(IKLE3(IELEM)) - F1
C
        U1  =  U(IKLE1(IELEM))
        U2  =  U(IKLE2(IELEM))
        U3  =  U(IKLE3(IELEM))
        U4  =  U(IKLE4(IELEM))
        V1  =  V(IKLE1(IELEM))
        V2  =  V(IKLE2(IELEM))
        V3  =  V(IKLE3(IELEM))
        V4  =  V(IKLE4(IELEM))
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM) = (-(2*X2*V4+4*X2*V2+2*X2*V1-X3*V4-2*X3*V2-X3
     * *V1+U4*Y3-2*U4*Y2+2*U2*Y3-4*U2*Y2+U1*Y3-2*U1*Y2)*(F3+
     * F2)*(Y3+Y2))*XMUL/(216*(X2*Y3-X3*Y2))
C
      A13(IELEM) = ((2*X2*V3+X2*V4+X2*V1-4*X3*V3-2*X3*V4-2*X3*
     * V1+4*U3*Y3-2*U3*Y2+2*U4*Y3-U4*Y2+2*U1*Y3-U1*Y2)*(F3*
     * Y2-F2*Y3))*XMUL/(72*(X2*Y3-X3*Y2))
C
      A14(IELEM) = (3*(X2*V3+2*X2*V4+X2*V1-2*X3*V3-4*X3*V4-2*
     * X3*V1+2*U3*Y3-U3*Y2+4*U4*Y3-2*U4*Y2+2*U1*Y3-U1*Y2)*(
     * F3*Y2-F2*Y3)-(4*X2*V4+2*X2*V2+2*X2*V1-2*X3*V4-X3*V2-
     * X3*V1+2*U4*Y3-4*U4*Y2+U2*Y3-2*U2*Y2+U1*Y3-2*U1*Y2)*(
     * F3+F2)*(Y3+Y2))*XMUL/(216*(X2*Y3-X3*Y2))
C
      A21(IELEM) = (-(X2*V4+X2*V2+2*X2*V1+X3*V4+X3*V2+2*X3*V1-U4
     * *Y3-U4*Y2-U2*Y3-U2*Y2-2*U1*Y3-2*U1*Y2)*(F3+F2)*(Y3+Y2))
     * *XMUL/(216*(X2*Y3-X3*Y2))
C
      A23(IELEM) = (-(2*X2*V3+X2*V4+X2*V2-4*X3*V3-2*X3*V4-2*X3
     * *V2+4*U3*Y3-2*U3*Y2+2*U4*Y3-U4*Y2+2*U2*Y3-U2*Y2)*(F3*
     * Y2-F2*Y3))*XMUL/(72*(X2*Y3-X3*Y2))
C
      A24(IELEM) = (-(3*(X2*V3+2*X2*V4+X2*V2-2*X3*V3-4*X3*V4-
     * 2*X3*V2+2*U3*Y3-U3*Y2+4*U4*Y3-2*U4*Y2+2*U2*Y3-U2*Y2)*
     * (F3*Y2-F2*Y3)+(2*X2*V4+X2*V2+X2*V1+2*X3*V4+X3*V2+X3*V1-
     * 2*U4*Y3-2*U4*Y2-U2*Y3-U2*Y2-U1*Y3-U1*Y2)*(F3+F2)*(Y3+Y2)
     * ))*XMUL/(216*(X2*Y3-X3*Y2))
C
      A31(IELEM) = (-(X2*V3+X2*V4+2*X2*V1+X3*V3+X3*V4+2*X3*V1-U3
     * *Y3-U3*Y2-U4*Y3-U4*Y2-2*U1*Y3-2*U1*Y2)*(F3*Y2-F2*Y3))*XMUL/(
     * 72*(X2*Y3-X3*Y2))
C
      A32(IELEM) = (-(2*X2*V3+2*X2*V4+4*X2*V2-X3*V3-X3*V4-2*X3
     * *V2+U3*Y3-2*U3*Y2+U4*Y3-2*U4*Y2+2*U2*Y3-4*U2*Y2)*(F3*
     * Y2-F2*Y3))*XMUL/(72*(X2*Y3-X3*Y2))
C
      A34(IELEM) = (-(3*X2*V3+6*X2*V4+2*X2*V2+X2*V1-X3*V2+X3*V1
     * -3*U3*Y2-6*U4*Y2+U2*Y3-2*U2*Y2-U1*Y3-U1*Y2)*(F3*Y2-F2*
     * Y3))*XMUL/(72*(X2*Y3-X3*Y2))
C
      A41(IELEM) = ((X2*V4+X2*V2+2*X2*V1-U4*Y2-U2*Y2-2*U1*Y2)*(
     * F3+F2)*(Y3+Y2)+3*(X3*V3+X3*V4+2*X3*V1-U3*Y3-U4*Y3-2*U1
     * *Y3)*(F3*Y2-F2*Y3))*XMUL/(72*(X2*Y3-X3*Y2))
C
      A42(IELEM) = (3*(X2*V3+X2*V4+2*X2*V2-X3*V3-X3*V4-2*X3*V2+
     * U3*Y3-U3*Y2+U4*Y3-U4*Y2+2*U2*Y3-2*U2*Y2)*(F3*Y2-F2*Y3)+
     * (X2*V4+2*X2*V2+X2*V1-U4*Y2-2*U2*Y2-U1*Y2)*(F3+F2)*(Y3+
     * Y2))*XMUL/(72*(X2*Y3-X3*Y2))
C
      A43(IELEM) = ((2*X2*V3+X2*V4+X2*V2-X3*V2+X3*V1-2*U3*Y2-U4*
     * Y2+U2*Y3-U2*Y2-U1*Y3)*(F3*Y2-F2*Y3))*XMUL/(24*(X2*Y3-X3*Y2))
C
C   TERMES DIAGONAUX
C   LA SOMME DE CHAQUE COLONNE EST NULLE.
C
        A11(IELEM) = - A21(IELEM)  - A31(IELEM) - A41(IELEM)
        A22(IELEM) = - A12(IELEM)  - A32(IELEM) - A42(IELEM)
        A33(IELEM) = - A13(IELEM)  - A23(IELEM) - A43(IELEM)
        A44(IELEM) = - A14(IELEM)  - A24(IELEM) - A34(IELEM)
C
3     CONTINUE
C
        ELSEIF(ICOORD.EQ.2) THEN
C
C================================
C  CAS DE LA DERIVEE SUIVANT Y  =
C================================
C
        DO 4 IELEM = 1 , NELEM
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
        X2  =  XEL(IELEM,2)
        X3  =  XEL(IELEM,3)
        Y2  =  YEL(IELEM,2)
        Y3  =  YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM)) - F1
        F3  =  F(IKLE3(IELEM)) - F1
C
        U1  =  U(IKLE1(IELEM))
        U2  =  U(IKLE2(IELEM))
        U3  =  U(IKLE3(IELEM))
        U4  =  U(IKLE4(IELEM))
        V1  =  V(IKLE1(IELEM))
        V2  =  V(IKLE2(IELEM))
        V3  =  V(IKLE3(IELEM))
        V4  =  V(IKLE4(IELEM))
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM) = (2*X2*V4+4*X2*V2+2*X2*V1-X3*V4-2*X3*V2-X3*
     * V1+U4*Y3-2*U4*Y2+2*U2*Y3-4*U2*Y2+U1*Y3-2*U1*Y2)*(X2+
     * X3)*(F3+F2)*XMUL/(216*(X2*Y3-X3*Y2))
C
      A13(IELEM) = -(2*X2*V3+X2*V4+X2*V1-4*X3*V3-2*X3*V4-2*X3
     * *V1+4*U3*Y3-2*U3*Y2+2*U4*Y3-U4*Y2+2*U1*Y3-U1*Y2)*(X2*
     * F3-X3*F2)*XMUL/(72*(X2*Y3-X3*Y2))
C
      A14(IELEM) = -(3*(X2*V3+2*X2*V4+X2*V1-2*X3*V3-4*X3*V4-
     * 2*X3*V1+2*U3*Y3-U3*Y2+4*U4*Y3-2*U4*Y2+2*U1*Y3-U1*Y2)*
     * (X2*F3-X3*F2)-(4*X2*V4+2*X2*V2+2*X2*V1-2*X3*V4-X3*V2-
     * X3*V1+2*U4*Y3-4*U4*Y2+U2*Y3-2*U2*Y2+U1*Y3-2*U1*Y2)*(
     * X2+X3)*(F3+F2))*XMUL/(216*(X2*Y3-X3*Y2))
C
      A21(IELEM) = (X2*V4+X2*V2+2*X2*V1+X3*V4+X3*V2+2*X3*V1-U4*
     * Y3-U4*Y2-U2*Y3-U2*Y2-2*U1*Y3-2*U1*Y2)*(X2+X3)*(F3+F2)*XMUL/
     * (216*(X2*Y3-X3*Y2))
C
      A23(IELEM) = (2*X2*V3+X2*V4+X2*V2-4*X3*V3-2*X3*V4-2*X3*
     * V2+4*U3*Y3-2*U3*Y2+2*U4*Y3-U4*Y2+2*U2*Y3-U2*Y2)*(X2*
     * F3-X3*F2)*XMUL/(72*(X2*Y3-X3*Y2))
C
      A24(IELEM) = (3*(X2*V3+2*X2*V4+X2*V2-2*X3*V3-4*X3*V4-2*
     * X3*V2+2*U3*Y3-U3*Y2+4*U4*Y3-2*U4*Y2+2*U2*Y3-U2*Y2)*(
     * X2*F3-X3*F2)+(2*X2*V4+X2*V2+X2*V1+2*X3*V4+X3*V2+X3*V1-
     * 2*U4*Y3-2*U4*Y2-U2*Y3-U2*Y2-U1*Y3-U1*Y2)*(X2+X3)*(F3+F2)
     * )*XMUL/(216*(X2*Y3-X3*Y2))
C
      A31(IELEM) = (X2*V3+X2*V4+2*X2*V1+X3*V3+X3*V4+2*X3*V1-U3*
     * Y3-U3*Y2-U4*Y3-U4*Y2-2*U1*Y3-2*U1*Y2)*(X2*F3-X3*F2)*XMUL/(
     * 72*(X2*Y3-X3*Y2))
C
      A32(IELEM) = (2*X2*V3+2*X2*V4+4*X2*V2-X3*V3-X3*V4-2*X3*
     * V2+U3*Y3-2*U3*Y2+U4*Y3-2*U4*Y2+2*U2*Y3-4*U2*Y2)*(X2*
     * F3-X3*F2)*XMUL/(72*(X2*Y3-X3*Y2))
C
      A34(IELEM) = (3*X2*V3+6*X2*V4+2*X2*V2+X2*V1-X3*V2+X3*V1-
     * 3*U3*Y2-6*U4*Y2+U2*Y3-2*U2*Y2-U1*Y3-U1*Y2)*(X2*F3-X3*F2
     * )*XMUL/(72*(X2*Y3-X3*Y2))
C
      A41(IELEM) = -((X2*V4+X2*V2+2*X2*V1-U4*Y2-U2*Y2-2*U1*Y2)*
     * (X2+X3)*(F3+F2)+3*(X2*F3-X3*F2)*(X3*V3+X3*V4+2*X3*V1-U3
     * *Y3-U4*Y3-2*U1*Y3))*XMUL/(72*(X2*Y3-X3*Y2))
C
      A42(IELEM) = -(3*(X2*V3+X2*V4+2*X2*V2-X3*V3-X3*V4-2*X3*
     * V2+U3*Y3-U3*Y2+U4*Y3-U4*Y2+2*U2*Y3-2*U2*Y2)*(X2*F3-X3*
     * F2)+(X2*V4+2*X2*V2+X2*V1-U4*Y2-2*U2*Y2-U1*Y2)*(X2+X3)*(
     * F3+F2))*XMUL/(72*(X2*Y3-X3*Y2))
C
      A43(IELEM) = -(2*X2*V3+X2*V4+X2*V2-X3*V2+X3*V1-2*U3*Y2-U4
     * *Y2+U2*Y3-U2*Y2-U1*Y3)*(X2*F3-X3*F2)*XMUL/(24*(X2*Y3-X3*Y2))
C
C   TERMES DIAGONAUX
C   LA SOMME DE CHAQUE COLONNE EST NULLE.
C
        A11(IELEM) = - A21(IELEM)  - A31(IELEM) - A41(IELEM)
        A22(IELEM) = - A12(IELEM)  - A32(IELEM) - A42(IELEM)
        A33(IELEM) = - A13(IELEM)  - A23(IELEM) - A43(IELEM)
        A44(IELEM) = - A14(IELEM)  - A24(IELEM) - A34(IELEM)
C
4       CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
          CALL PLANTE(0)
          STOP
        ENDIF
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,300) IELMU
        IF (LNG.EQ.2) WRITE(LU,301) IELMU
        CALL PLANTE(0)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C  CAS OU F EST DE DISCRETISATION QUASI-BULLE
C-----------------------------------------------------------------------
C
      ELSEIF(IELMF.EQ.12) THEN
C
      IF (IELMU.EQ.10) THEN
C
C================================
C  CAS DE LA DERIVEE SUIVANT X  =
C================================
C
        IF(ICOORD.EQ.1) THEN
C
C   BOUCLE SUR LES ELEMENTS
C
        DO 5 IELEM = 1 , NELEM
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
        X2 = XEL(IELEM,2)
        X3 = XEL(IELEM,3)
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM)) - F1
        F3  =  F(IKLE3(IELEM)) - F1
        F4  =  F(IKLE4(IELEM)) - F1
C
        UX  =  U(IELEM)
        UY  =  V(IELEM)
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM) = ((UX*Y3-2*UX*Y2+2*UY*X2-UY*X3)*(3*F4*Y2-F2*Y3-F2
     * *Y2))*XMUL/(18*(X2*Y3-X3*Y2))
C
      A13(IELEM) = ((2*UX*Y3-UX*Y2+UY*X2-2*UY*X3)*(F3*Y3+F3*Y2-3*F4
     * *Y3))*XMUL/(18*(X2*Y3-X3*Y2))
C
      A14(IELEM) = ((2*UX*Y3-UX*Y2+UY*X2-2*UY*X3)*(F3*Y3+F3*Y2-3*F4
     * *Y3)+(UX*Y3-2*UX*Y2+2*UY*X2-UY*X3)*(3*F4*Y2-F2*Y3-F2*Y2))*XMUL/
     * (18*(X2*Y3-X3*Y2))
C
      A21(IELEM) = (-(UX*Y3+UX*Y2-UY*X2-UY*X3)*(3*F4*Y2-F2*Y3-F2*Y2))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A23(IELEM) = ((2*UX*Y3-UX*Y2+UY*X2-2*UY*X3)*(F3*Y3-2*F3*Y2-3
     * *F4*Y3+3*F4*Y2+2*F2*Y3-F2*Y2))*XMUL/(18*(X2*Y3-X3*Y2))
C
      A24(IELEM) = ((2*UX*Y3-UX*Y2+UY*X2-2*UY*X3)*(F3*Y3-2*F3*Y2-3
     * *F4*Y3+3*F4*Y2+2*F2*Y3-F2*Y2)-(UX*Y3+UX*Y2-UY*X2-UY*X3)*(3
     * *F4*Y2-F2*Y3-F2*Y2))*XMUL/(18*(X2*Y3-X3*Y2))
C
      A31(IELEM) = ((UX*Y3+UX*Y2-UY*X2-UY*X3)*(F3*Y3+F3*Y2-3*F4*Y3))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A32(IELEM) = ((UX*Y3-2*UX*Y2+2*UY*X2-UY*X3)*(F3*Y3-2*F3*Y2-3
     * *F4*Y3+3*F4*Y2+2*F2*Y3-F2*Y2))*XMUL/(18*(X2*Y3-X3*Y2))
C
      A34(IELEM) = ((UX*Y3+UX*Y2-UY*X2-UY*X3)*(F3*Y3+F3*Y2-3*F4*Y3)+(
     * UX*Y3-2*UX*Y2+2*UY*X2-UY*X3)*(F3*Y3-2*F3*Y2-3*F4*Y3+3*F4
     * *Y2+2*F2*Y3-F2*Y2))*XMUL/(18*(X2*Y3-X3*Y2))
C
      A41(IELEM) = (-((UX*Y3-UY*X3)*(F3*Y3+F3*Y2-3*F4*Y3)-(UX*Y2-UY*
     * X2)*(3*F4*Y2-F2*Y3-F2*Y2)))*XMUL/(6*(X2*Y3-X3*Y2))
C
      A42(IELEM) = (-((UX*Y3-UX*Y2+UY*X2-UY*X3)*(F3*Y3-2*F3*Y2-3*F4*
     * Y3+3*F4*Y2+2*F2*Y3-F2*Y2)-(UX*Y2-UY*X2)*(3*F4*Y2-F2*Y3-
     * F2*Y2)))*XMUL/(6*(X2*Y3-X3*Y2))
C
      A43(IELEM) = (-((UX*Y3-UX*Y2+UY*X2-UY*X3)*(F3*Y3-2*F3*Y2-3*F4*
     * Y3+3*F4*Y2+2*F2*Y3-F2*Y2)+(UX*Y3-UY*X3)*(F3*Y3+F3*Y2-3*
     * F4*Y3)))*XMUL/(6*(X2*Y3-X3*Y2))

C
C
C   TERMES DIAGONAUX
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
C
        A11(IELEM) = - A21(IELEM)  - A31(IELEM) - A41(IELEM)
        A22(IELEM) = - A12(IELEM)  - A32(IELEM) - A42(IELEM)
        A33(IELEM) = - A13(IELEM)  - A23(IELEM) - A43(IELEM)
        A44(IELEM) = - A14(IELEM)  - A24(IELEM) - A34(IELEM)
C
5       CONTINUE
C
        ELSEIF(ICOORD.EQ.2) THEN
C
C================================
C  CAS DE LA DERIVEE SUIVANT Y  =
C================================
C
        DO 6 IELEM = 1 , NELEM
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
        X2  =  XEL(IELEM,2)
        X3  =  XEL(IELEM,3)
        Y2  =  YEL(IELEM,2)
        Y3  =  YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM)) - F1
        F3  =  F(IKLE3(IELEM)) - F1
        F4  =  F(IKLE4(IELEM)) - F1
C
        UX  =  U(IELEM)
        UY  =  V(IELEM)
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM) = (-(UX*Y3-2*UX*Y2+2*UY*X2-UY*X3)*(3*X2*F4-X2*F2-
     * X3*F2))*XMUL/(18*(X2*Y3-X3*Y2))
C
      A13(IELEM) = (-(2*UX*Y3-UX*Y2+UY*X2-2*UY*X3)*(X2*F3+X3*F3-3*
     * X3*F4))*XMUL/(18*(X2*Y3-X3*Y2))
C
      A14(IELEM) = (-((2*UX*Y3-UX*Y2+UY*X2-2*UY*X3)*(X2*F3+X3*F3-3*
     * X3*F4)+(UX*Y3-2*UX*Y2+2*UY*X2-UY*X3)*(3*X2*F4-X2*F2-X3*F2)
     * ))*XMUL/(18*(X2*Y3-X3*Y2))
C
      A21(IELEM) = ((UX*Y3+UX*Y2-UY*X2-UY*X3)*(3*X2*F4-X2*F2-X3*F2))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A23(IELEM) = ((2*UX*Y3-UX*Y2+UY*X2-2*UY*X3)*(2*X2*F3-3*X2*F4
     * +X2*F2-X3*F3+3*X3*F4-2*X3*F2))*XMUL/(18*(X2*Y3-X3*Y2))
C
      A24(IELEM) = ((2*UX*Y3-UX*Y2+UY*X2-2*UY*X3)*(2*X2*F3-3*X2*F4
     * +X2*F2-X3*F3+3*X3*F4-2*X3*F2)+(UX*Y3+UX*Y2-UY*X2-UY*X3)*(3
     * *X2*F4-X2*F2-X3*F2))*XMUL/(18*(X2*Y3-X3*Y2))
C
      A31(IELEM) = (-(UX*Y3+UX*Y2-UY*X2-UY*X3)*(X2*F3+X3*F3-3*X3*F4))
     * *XMUL/(18*(X2*Y3-X3*Y2))
C
      A32(IELEM) = ((UX*Y3-2*UX*Y2+2*UY*X2-UY*X3)*(2*X2*F3-3*X2*F4
     * +X2*F2-X3*F3+3*X3*F4-2*X3*F2))*XMUL/(18*(X2*Y3-X3*Y2))
C
      A34(IELEM) = (-((UX*Y3+UX*Y2-UY*X2-UY*X3)*(X2*F3+X3*F3-3*X3*F4)
     * -(UX*Y3-2*UX*Y2+2*UY*X2-UY*X3)*(2*X2*F3-3*X2*F4+X2*F2-X3*
     * F3+3*X3*F4-2*X3*F2)))*XMUL/(18*(X2*Y3-X3*Y2))
C
      A41(IELEM) = ((UX*Y3-UY*X3)*(X2*F3+X3*F3-3*X3*F4)-(UX*Y2-UY*X2)
     * *(3*X2*F4-X2*F2-X3*F2))*XMUL/(6*(X2*Y3-X3*Y2))
C
      A42(IELEM) = (-((UX*Y3-UX*Y2+UY*X2-UY*X3)*(2*X2*F3-3*X2*F4+X2*
     * F2-X3*F3+3*X3*F4-2*X3*F2)+(UX*Y2-UY*X2)*(3*X2*F4-X2*F2-
     * X3*F2)))*XMUL/(6*(X2*Y3-X3*Y2))
C
      A43(IELEM) = (-((UX*Y3-UX*Y2+UY*X2-UY*X3)*(2*X2*F3-3*X2*F4+X2*
     * F2-X3*F3+3*X3*F4-2*X3*F2)-(UX*Y3-UY*X3)*(X2*F3+X3*F3-3*
     * X3*F4)))*XMUL/(6*(X2*Y3-X3*Y2))
C
C
C   TERMES DIAGONAUX
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
C
        A11(IELEM) = - A21(IELEM)  - A31(IELEM) - A41(IELEM)
        A22(IELEM) = - A12(IELEM)  - A32(IELEM) - A42(IELEM)
        A33(IELEM) = - A13(IELEM)  - A23(IELEM) - A43(IELEM)
        A44(IELEM) = - A14(IELEM)  - A24(IELEM) - A34(IELEM)
C
6       CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,200) ICOORD
          CALL PLANTE(0)
          STOP
C
        ENDIF
C
C
      ELSEIF(IELMU.EQ.12) THEN
C
C================================
C  CAS DE LA DERIVEE SUIVANT X  =
C================================
C
        IF(ICOORD.EQ.1) THEN
C
C   BOUCLE SUR LES ELEMENTS
C
        DO 7 IELEM = 1 , NELEM
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
        X2 = XEL(IELEM,2)
        X3 = XEL(IELEM,3)
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM)) - F1
        F3  =  F(IKLE3(IELEM)) - F1
        F4  =  F(IKLE4(IELEM)) - F1
C
        U1  =  U(IKLE1(IELEM))
        U2  =  U(IKLE2(IELEM))
        U3  =  U(IKLE3(IELEM))
        U4  =  U(IKLE4(IELEM))
        V1  =  V(IKLE1(IELEM))
        V2  =  V(IKLE2(IELEM))
        V3  =  V(IKLE3(IELEM))
        V4  =  V(IKLE4(IELEM))
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM) = (-(2*X2*V4+4*X2*V2+2*X2*V1-X3*V4-2*X3*V2-X3
     * *V1+U4*Y3-2*U4*Y2+2*U2*Y3-4*U2*Y2+U1*Y3-2*U1*Y2)*(Y3+
     * Y2)*F4)*XMUL/(72*(X2*Y3-X3*Y2))
C
      A13(IELEM) = ((2*X2*V3+X2*V4+X2*V1-4*X3*V3-2*X3*V4-2*X3*
     * V1+4*U3*Y3-2*U3*Y2+2*U4*Y3-U4*Y2+2*U1*Y3-U1*Y2)*(F3*
     * Y3+F3*Y2-3*F4*Y3))*XMUL/(72*(X2*Y3-X3*Y2))
C
      A14(IELEM) = ((X2*V3+2*X2*V4+X2*V1-2*X3*V3-4*X3*V4-2*X3*
     * V1+2*U3*Y3-U3*Y2+4*U4*Y3-2*U4*Y2+2*U1*Y3-U1*Y2)*(F3*
     * Y3+F3*Y2-3*F4*Y3)-(4*X2*V4+2*X2*V2+2*X2*V1-2*X3*V4-
     * X3*V2-X3*V1+2*U4*Y3-4*U4*Y2+U2*Y3-2*U2*Y2+U1*Y3-2*U1*
     * Y2)*(Y3+Y2)*F4)*XMUL/(72*(X2*Y3-X3*Y2))
C
      A21(IELEM) = (-(X2*V4+X2*V2+2*X2*V1+X3*V4+X3*V2+2*X3*V1-U4
     * *Y3-U4*Y2-U2*Y3-U2*Y2-2*U1*Y3-2*U1*Y2)*(Y3+Y2)*F4)*XMUL/(72
     * *(X2*Y3-X3*Y2))
C
      A23(IELEM) = ((2*X2*V3+X2*V4+X2*V2-4*X3*V3-2*X3*V4-2*X3*
     * V2+4*U3*Y3-2*U3*Y2+2*U4*Y3-U4*Y2+2*U2*Y3-U2*Y2)*(F3*
     * Y3-2*F3*Y2-3*F4*Y3+3*F4*Y2+2*F2*Y3-F2*Y2))*XMUL/(72*(X2*
     * Y3-X3*Y2))
C
      A24(IELEM) = ((X2*V3+2*X2*V4+X2*V2-2*X3*V3-4*X3*V4-2*X3*
     * V2+2*U3*Y3-U3*Y2+4*U4*Y3-2*U4*Y2+2*U2*Y3-U2*Y2)*(F3*
     * Y3-2*F3*Y2-3*F4*Y3+3*F4*Y2+2*F2*Y3-F2*Y2)-(2*X2*V4+
     * X2*V2+X2*V1+2*X3*V4+X3*V2+X3*V1-2*U4*Y3-2*U4*Y2-U2*Y3-
     * U2*Y2-U1*Y3-U1*Y2)*(Y3+Y2)*F4)*XMUL/(72*(X2*Y3-X3*Y2))
C
      A31(IELEM) = (-(X2*V3+X2*V4+2*X2*V1+X3*V3+X3*V4+2*X3*V1-U3
     * *Y3-U3*Y2-U4*Y3-U4*Y2-2*U1*Y3-2*U1*Y2)*(F3*Y3+F3*Y2-3*
     * F4*Y3))*XMUL/(72*(X2*Y3-X3*Y2))
C
      A32(IELEM) = ((2*X2*V3+2*X2*V4+4*X2*V2-X3*V3-X3*V4-2*X3*
     * V2+U3*Y3-2*U3*Y2+U4*Y3-2*U4*Y2+2*U2*Y3-4*U2*Y2)*(F3*
     * Y3-2*F3*Y2-3*F4*Y3+3*F4*Y2+2*F2*Y3-F2*Y2))*XMUL/(72*(X2*
     * Y3-X3*Y2))
C
      A34(IELEM) = ((2*X2*V3+4*X2*V4+2*X2*V2-X3*V3-2*X3*V4-X3*
     * V2+U3*Y3-2*U3*Y2+2*U4*Y3-4*U4*Y2+U2*Y3-2*U2*Y2)*(F3*
     * Y3-2*F3*Y2-3*F4*Y3+3*F4*Y2+2*F2*Y3-F2*Y2)-(X2*V3+2*
     * X2*V4+X2*V1+X3*V3+2*X3*V4+X3*V1-U3*Y3-U3*Y2-2*U4*Y3-2*
     * U4*Y2-U1*Y3-U1*Y2)*(F3*Y3+F3*Y2-3*F4*Y3))*XMUL/(72*(X2*Y3-X3
     * *Y2))
C
      A41(IELEM) = ((X2*V4+X2*V2+2*X2*V1-U4*Y2-U2*Y2-2*U1*Y2)*(
     * Y3+Y2)*F4+(X3*V3+X3*V4+2*X3*V1-U3*Y3-U4*Y3-2*U1*Y3)*(F3
     * *Y3+F3*Y2-3*F4*Y3))*XMUL/(24*(X2*Y3-X3*Y2))
C
      A42(IELEM) = (-((X2*V3+X2*V4+2*X2*V2-X3*V3-X3*V4-2*X3*V2+
     * U3*Y3-U3*Y2+U4*Y3-U4*Y2+2*U2*Y3-2*U2*Y2)*(F3*Y3-2*F3*
     * Y2-3*F4*Y3+3*F4*Y2+2*F2*Y3-F2*Y2)-(X2*V4+2*X2*V2+X2*
     * V1-U4*Y2-2*U2*Y2-U1*Y2)*(Y3+Y2)*F4))*XMUL/(24*(X2*Y3-X3*Y2))
C
      A43(IELEM) = (-((2*X2*V3+X2*V4+X2*V2-2*X3*V3-X3*V4-X3*V2+
     * 2*U3*Y3-2*U3*Y2+U4*Y3-U4*Y2+U2*Y3-U2*Y2)*(F3*Y3-2*F3*Y2
     * -3*F4*Y3+3*F4*Y2+2*F2*Y3-F2*Y2)-(2*X3*V3+X3*V4+X3*V1-
     * 2*U3*Y3-U4*Y3-U1*Y3)*(F3*Y3+F3*Y2-3*F4*Y3)))*XMUL/(24*(X2*Y3
     * -X3*Y2))
C
C
C   TERMES DIAGONAUX
C   LA SOMME DE CHAQUE COLONNE EST NULLE
C
        A11(IELEM) = - A21(IELEM)  - A31(IELEM) - A41(IELEM)
        A22(IELEM) = - A12(IELEM)  - A32(IELEM) - A42(IELEM)
        A33(IELEM) = - A13(IELEM)  - A23(IELEM) - A43(IELEM)
        A44(IELEM) = - A14(IELEM)  - A24(IELEM) - A34(IELEM)
C
7       CONTINUE
C
        ELSEIF(ICOORD.EQ.2) THEN
C
C================================
C  CAS DE LA DERIVEE SUIVANT Y  =
C================================
C
        DO 8 IELEM = 1 , NELEM
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
        X2  =  XEL(IELEM,2)
        X3  =  XEL(IELEM,3)
        Y2  =  YEL(IELEM,2)
        Y3  =  YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM)) - F1
        F3  =  F(IKLE3(IELEM)) - F1
        F4  =  F(IKLE4(IELEM)) - F1
C
        U1  =  U(IKLE1(IELEM))
        U2  =  U(IKLE2(IELEM))
        U3  =  U(IKLE3(IELEM))
        U4  =  U(IKLE4(IELEM))
        V1  =  V(IKLE1(IELEM))
        V2  =  V(IKLE2(IELEM))
        V3  =  V(IKLE3(IELEM))
        V4  =  V(IKLE4(IELEM))
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM) = ((2*X2*V4+4*X2*V2+2*X2*V1-X3*V4-2*X3*V2-X3*
     * V1+U4*Y3-2*U4*Y2+2*U2*Y3-4*U2*Y2+U1*Y3-2*U1*Y2)*(X2+
     * X3)*F4)*XMUL/(72*(X2*Y3-X3*Y2))
C
      A13(IELEM) = (-(X2*F3+X3*F3-3*X3*F4)*(2*X2*V3+X2*V4+X2*V1-
     * 4*X3*V3-2*X3*V4-2*X3*V1+4*U3*Y3-2*U3*Y2+2*U4*Y3-U4*
     * Y2+2*U1*Y3-U1*Y2))*XMUL/(72*(X2*Y3-X3*Y2))
C
      A14(IELEM) = (-((X2*F3+X3*F3-3*X3*F4)*(X2*V3+2*X2*V4+X2*V1
     * -2*X3*V3-4*X3*V4-2*X3*V1+2*U3*Y3-U3*Y2+4*U4*Y3-2*U4
     * *Y2+2*U1*Y3-U1*Y2)-(4*X2*V4+2*X2*V2+2*X2*V1-2*X3*V4-
     * X3*V2-X3*V1+2*U4*Y3-4*U4*Y2+U2*Y3-2*U2*Y2+U1*Y3-2*U1*
     * Y2)*(X2+X3)*F4))*XMUL/(72*(X2*Y3-X3*Y2))
C
      A21(IELEM) = ((X2*V4+X2*V2+2*X2*V1+X3*V4+X3*V2+2*X3*V1-U4*
     * Y3-U4*Y2-U2*Y3-U2*Y2-2*U1*Y3-2*U1*Y2)*(X2+X3)*F4)*XMUL/(72*
     * (X2*Y3-X3*Y2))
C
      A23(IELEM) = ((2*X2*F3-3*X2*F4+X2*F2-X3*F3+3*X3*F4-2*X3*
     * F2)*(2*X2*V3+X2*V4+X2*V2-4*X3*V3-2*X3*V4-2*X3*V2+4*
     * U3*Y3-2*U3*Y2+2*U4*Y3-U4*Y2+2*U2*Y3-U2*Y2))*XMUL/(72*(X2*
     * Y3-X3*Y2))
C
      A24(IELEM) = ((2*X2*F3-3*X2*F4+X2*F2-X3*F3+3*X3*F4-2*X3*
     * F2)*(X2*V3+2*X2*V4+X2*V2-2*X3*V3-4*X3*V4-2*X3*V2+2*
     * U3*Y3-U3*Y2+4*U4*Y3-2*U4*Y2+2*U2*Y3-U2*Y2)+(2*X2*V4+
     * X2*V2+X2*V1+2*X3*V4+X3*V2+X3*V1-2*U4*Y3-2*U4*Y2-U2*Y3-
     * U2*Y2-U1*Y3-U1*Y2)*(X2+X3)*F4)*XMUL/(72*(X2*Y3-X3*Y2))
C
      A31(IELEM) = ((X2*F3+X3*F3-3*X3*F4)*(X2*V3+X2*V4+2*X2*V1+
     * X3*V3+X3*V4+2*X3*V1-U3*Y3-U3*Y2-U4*Y3-U4*Y2-2*U1*Y3-2*
     * U1*Y2))*XMUL/(72*(X2*Y3-X3*Y2))
C
      A32(IELEM) = ((2*X2*F3-3*X2*F4+X2*F2-X3*F3+3*X3*F4-2*X3*
     * F2)*(2*X2*V3+2*X2*V4+4*X2*V2-X3*V3-X3*V4-2*X3*V2+U3*
     * Y3-2*U3*Y2+U4*Y3-2*U4*Y2+2*U2*Y3-4*U2*Y2))*XMUL/(72*(X2*
     * Y3-X3*Y2))
C
      A34(IELEM) = ((2*X2*F3-3*X2*F4+X2*F2-X3*F3+3*X3*F4-2*X3*
     * F2)*(2*X2*V3+4*X2*V4+2*X2*V2-X3*V3-2*X3*V4-X3*V2+U3*
     * Y3-2*U3*Y2+2*U4*Y3-4*U4*Y2+U2*Y3-2*U2*Y2)+(X2*F3+X3*
     * F3-3*X3*F4)*(X2*V3+2*X2*V4+X2*V1+X3*V3+2*X3*V4+X3*V1-
     * U3*Y3-U3*Y2-2*U4*Y3-2*U4*Y2-U1*Y3-U1*Y2))*XMUL/(72*(X2*Y3-
     * X3*Y2))
C
      A41(IELEM) = (-((X2*F3+X3*F3-3*X3*F4)*(X3*V3+X3*V4+2*X3*V1
     * -U3*Y3-U4*Y3-2*U1*Y3)+(X2*V4+X2*V2+2*X2*V1-U4*Y2-U2*Y2-
     * 2*U1*Y2)*(X2+X3)*F4))*XMUL/(24*(X2*Y3-X3*Y2))
C
      A42(IELEM) = (-((2*X2*F3-3*X2*F4+X2*F2-X3*F3+3*X3*F4-2*
     * X3*F2)*(X2*V3+X2*V4+2*X2*V2-X3*V3-X3*V4-2*X3*V2+U3*Y3-
     * U3*Y2+U4*Y3-U4*Y2+2*U2*Y3-2*U2*Y2)+(X2*V4+2*X2*V2+X2*
     * V1-U4*Y2-2*U2*Y2-U1*Y2)*(X2+X3)*F4))*XMUL/(24*(X2*Y3-X3*Y2))
C
      A43(IELEM) = (-((2*X2*F3-3*X2*F4+X2*F2-X3*F3+3*X3*F4-2*
     * X3*F2)*(2*X2*V3+X2*V4+X2*V2-2*X3*V3-X3*V4-X3*V2+2*U3*
     * Y3-2*U3*Y2+U4*Y3-U4*Y2+U2*Y3-U2*Y2)+(X2*F3+X3*F3-3*X3*
     * F4)*(2*X3*V3+X3*V4+X3*V1-2*U3*Y3-U4*Y3-U1*Y3)))*XMUL/(24*(
     * X2*Y3-X3*Y2))
C
C
C   TERMES DIAGONAUX
C   LA SOMME DE CHAQUE COLONNE EST NULLE.
C
        A11(IELEM) = - A21(IELEM)  - A31(IELEM) - A41(IELEM)
        A22(IELEM) = - A12(IELEM)  - A32(IELEM) - A42(IELEM)
        A33(IELEM) = - A13(IELEM)  - A23(IELEM) - A43(IELEM)
        A44(IELEM) = - A14(IELEM)  - A24(IELEM) - A34(IELEM)
C
8       CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
          CALL PLANTE(0)
          STOP
        ENDIF
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,300) IELMU
        IF (LNG.EQ.2) WRITE(LU,301) IELMU
        CALL PLANTE(0)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
       IF (LNG.EQ.1) WRITE(LU,100) IELMF
       IF (LNG.EQ.2) WRITE(LU,101) IELMF
100    FORMAT(1X,'MT12BB (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,' NON PREVUE')
101    FORMAT(1X,'MT12BB (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F : ',1I6,' NOT AVAILABLE')
       CALL PLANTE(0)
       STOP
      ENDIF
C
200       FORMAT(1X,'MT12BB (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'MT12BB (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
C
300    FORMAT(1X,'MT12BB (BIEF) :',/,
     *        1X,'DISCRETISATION DE U : ',1I6,' NON PREVUE')
301    FORMAT(1X,'MT12BB (BIEF) :',/,
     *        1X,'DISCRETIZATION OF U : ',1I6,' NOT AVAILABLE')
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
