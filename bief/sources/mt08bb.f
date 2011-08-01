C                       *****************
                        SUBROUTINE MT08BB
C                       *****************
C
     *(  A11 , A12 , A13 , A14 ,
     *   A21 , A22 , A23 , A24 ,
     *   A31 , A32 , A33 , A34 ,
     *   A41 , A42 , A43 , A44 ,
     *   XMUL,SF,F,XEL,YEL,IKLE1,IKLE2,IKLE3,IKLE4,NELEM,NELMAX,ICOORD)
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
C                 /                     D
C A(I,J)=-XMUL   /  PSI2(J) *    F    * --( PSI1(I) ) D(OMEGA)
C               /OMEGA                  DX
C
C  ATTENTION AU SIGNE MOINS ||
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
      USE BIEF, EX_MT08BB => MT08BB
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
      DOUBLE PRECISION, INTENT(IN) :: F(*)
C
C     STRUCTURE DE F
      TYPE(BIEF_OBJ), INTENT(IN) :: SF
C
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF
      DOUBLE PRECISION X2,X3,Y2,Y3,F1,F2,F3,F4
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
C
C-----------------------------------------------------------------------
C  CAS OU F EST DE DISCRETISATION P1
C-----------------------------------------------------------------------
C
      IF(IELMF.EQ.11) THEN
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
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM))
        F3  =  F(IKLE3(IELEM))
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM) = (2*Y2*(-F3-7*F2-4*F1)+Y3*(F3+7*F2+4*F1))*XMUL/216
        A13(IELEM) = (Y2*(-7*F3-F2-4*F1)+2*Y3*(7*F3+F2+4*F1))*XMUL/216
        A14(IELEM) = (Y2*(-3*F3-4*F2-5*F1)+Y3*(4*F3+3*F2+5*F1))*XMUL/72
        A21(IELEM) = (Y2*(-F3-4*F2-7*F1)+Y3*(-F3-4*F2-7*F1))*XMUL/216
        A23(IELEM) = (Y2*(7*F3+4*F2+F1)+2*Y3*(-7*F3-4*F2-F1))*XMUL/216
        A24(IELEM) = (Y2*(F3-F1)+Y3*(-4*F3-5*F2-3*F1))*XMUL/72
        A31(IELEM) = (Y2*(4*F3+F2+7*F1)+Y3*(4*F3+F2+7*F1))*XMUL/216
        A32(IELEM) = (2*Y2*(4*F3+7*F2+F1)+Y3*(-4*F3-7*F2-F1))*XMUL/216
        A34(IELEM) = (Y2*(5*F3+4*F2+3*F1)+Y3*(-F2+F1))*XMUL/72
        A41(IELEM) = (Y2*(F3+4*F2+7*F1)+Y3*(-4*F3-F2-7*F1))*XMUL/72
        A42(IELEM) = (3*Y2*(-F3+F1)+Y3*(4*F3+7*F2+F1))*XMUL/72
        A43(IELEM) = (Y2*(-7*F3-4*F2-F1)+3*Y3*(F2-F1))*XMUL/72
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
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM))
        F3  =  F(IKLE3(IELEM))
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM) = (2*X2*(F3+7*F2+4*F1)+X3*(-F3-7*F2-4*F1))*XMUL/216
        A13(IELEM) = (X2*(7*F3+F2+4*F1)+2*X3*(-7*F3-F2-4*F1))*XMUL/216
        A14(IELEM) = (X2*(3*F3+4*F2+5*F1)+X3*(-4*F3-3*F2-5*F1))*XMUL/72
        A21(IELEM) = (X2*(F3+4*F2+7*F1)+X3*(F3+4*F2+7*F1))*XMUL/216
        A23(IELEM) = (X2*(-7*F3-4*F2-F1)+2*X3*(7*F3+4*F2+F1))*XMUL/216
        A24(IELEM) = (X2*(-F3+F1)+X3*(4*F3+5*F2+3*F1))*XMUL/72
        A31(IELEM) = (X2*(-4*F3-F2-7*F1)+X3*(-4*F3-F2-7*F1))*XMUL/216
        A32(IELEM) = (2*X2*(-4*F3-7*F2-F1)+X3*(4*F3+7*F2+F1))*XMUL/216
        A34(IELEM) = (X2*(-5*F3-4*F2-3*F1)+X3*(F2-F1))*XMUL/72
        A41(IELEM) = (X2*(-F3-4*F2-7*F1)+X3*(4*F3+F2+7*F1))*XMUL/72
        A42(IELEM) = (3*X2*(F3-F1)+X3*(-4*F3-7*F2-F1))*XMUL/72
        A43(IELEM) = (X2*(7*F3+4*F2+F1)+3*X3*(-F2+F1))*XMUL/72
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
C-----------------------------------------------------------------------
C  CAS OU F EST DE DISCRETISATION QUASI-BULLE
C-----------------------------------------------------------------------
C
      ELSEIF(IELMF.EQ.12) THEN
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
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM))
        F3  =  F(IKLE3(IELEM))
        F4  =  F(IKLE4(IELEM))
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM) = (2*Y2*(-F4-2*F2-F1)+Y3*(F4+2*F2+F1))*XMUL/72
        A13(IELEM) = (Y2*(-2*F3-F4-F1)+2*Y3*(2*F3+F4+F1))*XMUL/72
        A14(IELEM) = (Y2*(-F3-6*F4-2*F2-3*F1)+Y3*(2*F3+6*F4+F2+
     *                3*F1))*XMUL/72
        A21(IELEM) = (Y2*(-F4-F2-2*F1)+Y3*(-F4-F2-2*F1))*XMUL/72
        A23(IELEM) = (Y2*(2*F3+F4+F2)+2*Y3*(-2*F3-F4-F2))*XMUL/72
        A24(IELEM) = (Y2*(F3-F1)+Y3*(-2*F3-6*F4-3*F2-F1))*XMUL/72
        A31(IELEM) = (Y2*(F3+F4+2*F1)+Y3*(F3+F4+2*F1))*XMUL/72
        A32(IELEM) = (2*Y2*(F3+F4+2*F2)+Y3*(-F3-F4-2*F2))*XMUL/72
        A34(IELEM) = (Y2*(3*F3+6*F4+2*F2+F1)+Y3*(-F2+F1))*XMUL/72
        A41(IELEM) = (Y2*(F4+F2+2*F1)+Y3*(-F3-F4-2*F1))*XMUL/24
        A42(IELEM) = (Y2*(-F3+F1)+Y3*(F3+F4+2*F2))*XMUL/24
        A43(IELEM) = (Y2*(-2*F3-F4-F2)+Y3*(F2-F1))*XMUL/24
C
C   TERMES DIAGONAUX
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
C
        A11(IELEM) = - A21(IELEM)  - A31(IELEM) - A41(IELEM)
        A22(IELEM) = - A12(IELEM)  - A32(IELEM) - A42(IELEM)
        A33(IELEM) = - A13(IELEM)  - A23(IELEM) - A43(IELEM)
        A44(IELEM) = - A14(IELEM)  - A24(IELEM) - A34(IELEM)
C
3       CONTINUE
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
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM))
        F3  =  F(IKLE3(IELEM))
        F4  =  F(IKLE4(IELEM))
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM) = (2*X2*(F4+2*F2+F1)+X3*(-F4-2*F2-F1))*XMUL/72
        A13(IELEM) = (X2*(2*F3+F4+F1)+2*X3*(-2*F3-F4-F1))*XMUL/72
        A14(IELEM) = (X2*(F3+6*F4+2*F2+3*F1)+X3*(-2*F3-6*F4-F2-
     *                3*F1))*XMUL/72
        A21(IELEM) =  (X2*(F4+F2+2*F1)+X3*(F4+F2+2*F1))*XMUL/72
        A23(IELEM) =  (X2*(-2*F3-F4-F2)+2*X3*(2*F3+F4+F2))*XMUL/72
        A24(IELEM) =  (X2*(-F3+F1)+X3*(2*F3+6*F4+3*F2+F1))*XMUL/72
        A31(IELEM) =  (X2*(-F3-F4-2*F1)+X3*(-F3-F4-2*F1))*XMUL/72
        A32(IELEM) =  (2*X2*(-F3-F4-2*F2)+X3*(F3+F4+2*F2))*XMUL/72
        A34(IELEM) =  (X2*(-3*F3-6*F4-2*F2-F1)+X3*(F2-F1))*XMUL/72
        A41(IELEM) =  (X2*(-F4-F2-2*F1)+X3*(F3+F4+2*F1))*XMUL/24
        A42(IELEM) =  (X2*(F3-F1)+X3*(-F3-F4-2*F2))*XMUL/24
        A43(IELEM) =  (X2*(2*F3+F4+F2)+X3*(-F2+F1))*XMUL/24
C
C   TERMES DIAGONAUX
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
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
C-----------------------------------------------------------------------
C
      ELSE
       IF (LNG.EQ.1) WRITE(LU,100) IELMF
       IF (LNG.EQ.2) WRITE(LU,101) IELMF
100    FORMAT(1X,'MT08BB (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,' NON PREVUE')
101    FORMAT(1X,'MT08BB (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F : ',1I6,' NOT AVAILABLE')
       CALL PLANTE(0)
       STOP
      ENDIF
C
200   FORMAT(1X,'MT08BB (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *       1I6,' VERIFIER ICOORD')
201   FORMAT(1X,'MT08BB (BIEF) : IMPOSSIBLE COMPONENT ',
     *       1I6,' CHECK ICOORD')
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
