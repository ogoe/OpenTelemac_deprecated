C                       *****************
                        SUBROUTINE MT12AB
C                       *****************
C
     *(  A11 , A12 , A13 , A14 ,
     *   A21 , A22 , A23 , A24 ,
     *   A31 , A32 , A33 , A34 ,
     *   XMUL,SF,SU,SV,F,U,V,
     *   XEL,YEL,SURFAC,
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
C  PSI2 : BASES DE TYPE QUASI-BULLE
C  F    : FONCTION DE TYPE TRIANGLE P1
C  U    : VECTEUR DE TYPE P0 OU QUASI-BULLE
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
      USE BIEF, EX_MT12AB => MT12AB
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
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: F(*),U(*),V(*)
C
C     STRUCTURES DE F,U,V
      TYPE(BIEF_OBJ), INTENT(IN) :: SF,SU,SV
C
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
      DOUBLE PRECISION, INTENT(IN) :: SURFAC(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF,IELMU,IELMV
      DOUBLE PRECISION X2,X3,Y2,Y3,F1,F2,F3
      DOUBLE PRECISION U1,U2,U3,U4,V1,V2,V3,V4,UX,UY
      DOUBLE PRECISION XSUR18,XSUR12,XSUR72,XSU144
      DOUBLE PRECISION AUX18,AUX12,AUX144,AUX72,UNSURF
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
      IELMU=SU%ELM
      IELMV=SV%ELM
C
      XSUR18 = XMUL /  18.D0
      XSUR12 = XMUL /  12.D0
      XSUR72 = XMUL /  72.D0
      XSU144 = XMUL / 144.D0
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
        UNSURF= 1.D0 / SURFAC(IELEM)
        AUX12 = XSUR12 * UNSURF
        AUX18 = XSUR18 * UNSURF
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM)=-( (Y3-Y2)*UX+X2*UY-X3*UY )*(Y3*F2-Y2*F3)*AUX18
        A14(IELEM)=-( (Y3-Y2)*UX+X2*UY-X3*UY )*(Y3*F2-Y2*F3)*AUX12
        A21(IELEM)=-(X3*UY-UX*Y3)*(Y3*F2-Y2*F3)*AUX18
        A24(IELEM)=-(X3*UY-UX*Y3)*(Y3*F2-Y2*F3)*AUX12
        A31(IELEM)= (X2*UY-UX*Y2)*(Y3*F2-Y2*F3)*AUX18
        A34(IELEM)= (X2*UY-UX*Y2)*(Y3*F2-Y2*F3)*AUX12
        A13(IELEM)= A12(IELEM)
        A23(IELEM)= A21(IELEM)
        A32(IELEM)= A31(IELEM)
C
C   TERMES DIAGONAUX
C   LA SOMME DE CHAQUE COLONNE EST NULLE
C
        A11(IELEM) = - A21(IELEM) - A31(IELEM)
        A22(IELEM) = - A12(IELEM) - A32(IELEM)
        A33(IELEM) = - A13(IELEM) - A23(IELEM)
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
        UNSURF= 1.D0 / SURFAC(IELEM)
        AUX12 = XSUR12 * UNSURF
        AUX18 = XSUR18 * UNSURF
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM)=-(X2*UY-X3*UY+UX*Y3-UX*Y2)*(X2*F3-X3*F2)*AUX18
        A14(IELEM)=-(X2*UY-X3*UY+UX*Y3-UX*Y2)*(X2*F3-X3*F2)*AUX12
        A21(IELEM)=-(X2*F3-X3*F2)*(X3*UY-UX*Y3)*AUX18
        A24(IELEM)=-(X2*F3-X3*F2)*(X3*UY-UX*Y3)*AUX12
        A31(IELEM)= (X2*UY-UX*Y2)*(X2*F3-X3*F2)*AUX18
        A34(IELEM)= (X2*UY-UX*Y2)*(X2*F3-X3*F2)*AUX12
        A13(IELEM)= A12(IELEM)
        A23(IELEM)= A21(IELEM)
        A32(IELEM)= A31(IELEM)
C
C   TERMES DIAGONAUX
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
C
        A11(IELEM) = - A21(IELEM) - A31(IELEM)
        A22(IELEM) = - A12(IELEM) - A32(IELEM)
        A33(IELEM) = - A13(IELEM) - A23(IELEM)
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
        UNSURF= 1.D0 / SURFAC(IELEM)
        AUX144= XSU144 * UNSURF
        AUX72 = XSUR72 * UNSURF
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM)=(-((V3+V4+2*V2)*X2-(V3+V4+2*V2)*X3+(V4+2*V2+
     *   V1)*X2-(V4+2*V2+V1)*X3+(Y3-Y2)*U3+2*(Y3-Y2)*U4+4*(Y3-
     *   Y2)*U2+(Y3-Y2)*U1)*(Y3*F2-Y2*F3))*AUX144
        A13(IELEM)=(-((2*V3+V4+V2)*X2-(2*V3+V4+V2)*X3+(2*V3+V4+
     *   V1)*X2-(2*V3+V4+V1)*X3+4*(Y3-Y2)*U3+2*(Y3-Y2)*U4+(Y3-
     *   Y2)*U2+(Y3-Y2)*U1)*(Y3*F2-Y2*F3))*AUX144
        A14(IELEM)=(-(X2*V3+3*X2*V4+X2*V2+X2*V1-X3*V3-3*X3*V4-X3
     *   *V2-X3*V1+U3*Y3-U3*Y2+3*U4*Y3-3*U4*Y2+U2*Y3-U2*Y2+U1*Y3
     *   -U1*Y2)*(Y3*F2-Y2*F3))*AUX72
        A21(IELEM)=(-(X3*V3+2*X3*V4+X3*V2+4*X3*V1-U3*Y3-2*U4*Y3
     *   -U2*Y3-4*U1*Y3)*(Y3*F2-Y2*F3))*AUX144
        A23(IELEM)=(-(4*X3*V3+2*X3*V4+X3*V2+X3*V1-4*U3*Y3-2*U4
     *   *Y3-U2*Y3-U1*Y3)*(Y3*F2-Y2*F3))*AUX144
        A24(IELEM)=(-(X3*V3+3*X3*V4+X3*V2+X3*V1-U3*Y3-3*U4*Y3-U2
     *   *Y3-U1*Y3)*(Y3*F2-Y2*F3))*AUX72
        A31(IELEM)=((X2*V3+2*X2*V4+X2*V2+4*X2*V1-U3*Y2-2*U4*Y2-
     *   U2*Y2-4*U1*Y2)*(Y3*F2-Y2*F3))*AUX144
        A32(IELEM)=((X2*V3+2*X2*V4+4*X2*V2+X2*V1-U3*Y2-2*U4*Y2-
     *   4*U2*Y2-U1*Y2)*(Y3*F2-Y2*F3))*AUX144
        A34(IELEM)=((X2*V3+3*X2*V4+X2*V2+X2*V1-U3*Y2-3*U4*Y2-U2*
     *   Y2-U1*Y2)*(Y3*F2-Y2*F3))*AUX72
C
C   TERMES DIAGONAUX
C   LA SOMME DE CHAQUE COLONNE EST NULLE
C
        A11(IELEM) = - A21(IELEM) - A31(IELEM)
        A22(IELEM) = - A12(IELEM) - A32(IELEM)
        A33(IELEM) = - A13(IELEM) - A23(IELEM)
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
        UNSURF= 1.D0 / SURFAC(IELEM)
        AUX144= XSU144 * UNSURF
        AUX72= XSUR72 * UNSURF
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM)=(-(X2*V3+2*X2*V4+4*X2*V2+X2*V1-X3*V3-2*X3*V4
     *   -4*X3*V2-X3*V1+U3*Y3-U3*Y2+2*U4*Y3-2*U4*Y2+4*U2*Y3-4
     *   *U2*Y2+U1*Y3-U1*Y2)*(X2*F3-X3*F2))*AUX144
        A13(IELEM)=(-(4*X2*V3+2*X2*V4+X2*V2+X2*V1-4*X3*V3-2*X3
     *   *V4-X3*V2-X3*V1+4*U3*Y3-4*U3*Y2+2*U4*Y3-2*U4*Y2+U2*Y3
     *   -U2*Y2+U1*Y3-U1*Y2)*(X2*F3-X3*F2))*AUX144
        A14(IELEM)=(-(X2*V3+3*X2*V4+X2*V2+X2*V1-X3*V3-3*X3*V4-X3
     *   *V2-X3*V1+U3*Y3-U3*Y2+3*U4*Y3-3*U4*Y2+U2*Y3-U2*Y2+U1*Y3
     *   -U1*Y2)*(X2*F3-X3*F2))*AUX72
        A21(IELEM)=(-(X2*F3-X3*F2)*(X3*V3+2*X3*V4+X3*V2+4*X3*V1-
     *   U3*Y3-2*U4*Y3-U2*Y3-4*U1*Y3))*AUX144
        A23(IELEM)=(-(X2*F3-X3*F2)*(4*X3*V3+2*X3*V4+X3*V2+X3*V1-
     *   4*U3*Y3-2*U4*Y3-U2*Y3-U1*Y3))*AUX144
        A24(IELEM)=(-(X2*F3-X3*F2)*(X3*V3+3*X3*V4+X3*V2+X3*V1-U3*
     *   Y3-3*U4*Y3-U2*Y3-U1*Y3))*AUX72
        A31(IELEM)=((X2*V3+2*X2*V4+X2*V2+4*X2*V1-U3*Y2-2*U4*Y2-
     *   U2*Y2-4*U1*Y2)*(X2*F3-X3*F2))*AUX144
        A32(IELEM)=((X2*V3+2*X2*V4+4*X2*V2+X2*V1-U3*Y2-2*U4*Y2-
     *   4*U2*Y2-U1*Y2)*(X2*F3-X3*F2))*AUX144
        A34(IELEM)=((X2*V3+3*X2*V4+X2*V2+X2*V1-U3*Y2-3*U4*Y2-U2*
     *   Y2-U1*Y2)*(X2*F3-X3*F2))*AUX72
C
C   TERMES DIAGONAUX
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
C
        A11(IELEM) = - A21(IELEM) - A31(IELEM)
        A22(IELEM) = - A12(IELEM) - A32(IELEM)
        A33(IELEM) = - A13(IELEM) - A23(IELEM)
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
C
      ELSE
       IF (LNG.EQ.1) WRITE(LU,100) IELMF
       IF (LNG.EQ.2) WRITE(LU,101) IELMF
100    FORMAT(1X,'MT12AB (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,' NON PREVUE')
101    FORMAT(1X,'MT12AB (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F : ',1I6,' NOT AVAILABLE')
       CALL PLANTE(0)
       STOP
      ENDIF
C
200       FORMAT(1X,'MT12AB (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'MT12AB (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
C
300    FORMAT(1X,'MT12AB (BIEF) :',/,
     *        1X,'DISCRETISATION DE U : ',1I6,' NON PREVUE')
301    FORMAT(1X,'MT12AB (BIEF) :',/,
     *        1X,'DISCRETIZATION OF U : ',1I6,' NOT AVAILABLE')
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
