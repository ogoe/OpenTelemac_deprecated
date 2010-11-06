C                       *****************
                        SUBROUTINE MT12AA
C                       *****************
C
     *(  A11 , A12 , A13 ,
     *   A21 , A22 , A23 ,
     *   A31 , A32 , A33 ,
     *   XMUL,SF,SU,SV,F,U,V,
     *   XEL,YEL,SURFAC,
     *   IKLE1,IKLE2,IKLE3,
     *   NELEM,NELMAX,ICOORD)
C
C***********************************************************************
C BIEF VERSION 5.1           09/12/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        C  MOULIN    (LNH) 30 87 83 81
C***********************************************************************
C
C  FONCTION : CALCUL DES COEFFICIENTS DE LA MATRICE SUIVANTE
C
C
C                 /                -->  --->
C A(I,J) = XMUL  /   PSIJ GRAD(F)   U .GRAD(PSII) D(OMEGA)
C               /OMEGA
C
C                                -->
C  F VECTEUR                      U  VECTEUR DE COMPOSANTES U,V,W
C
C  ATTENTION AU SIGNE MOINS ||
C
C  AVEC ICOORD=2 ON AURAIT UNE DERIVEE SUIVANT Y
C  AVEC ICOORD=3 ON AURAIT UNE DERIVEE SUIVANT Z
C
C  PSI1 : LINEAIRES
C  PSI2 : LINEAIRES
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
      USE BIEF, EX_MT12AA => MT12AA
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX,ICOORD
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX),IKLE3(NELMAX)
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A21(*),A22(*),A23(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A31(*),A32(*),A33(*)
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
      DOUBLE PRECISION XSUR12,XSUR48,X2,X3,Y2,Y3,F1,F2,F3,DEN
      DOUBLE PRECISION U1,U2,U3,V1,V2,V3,U123,V123
C
C-----------------------------------------------------------------------
C
      XSUR12 = XMUL/12.D0
      XSUR48 = XMUL/48.D0
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
      IELMU=SU%ELM
      IELMV=SV%ELM
C
C  CAS OU F EST DE DISCRETISATION P1 ET U P0
C
      IF(IELMF.EQ.11.AND.IELMU.EQ.10.AND.IELMV.EQ.10) THEN
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
      DEN = (F3*Y2 - F2*Y3) * XSUR12 / SURFAC(IELEM)
C
C   TERMES EXTRADIAGONAUX
C
      A23(IELEM) = (  X3 *V(IELEM) -  Y3    *U(IELEM) )*DEN
      A31(IELEM) =-(  X2 *V(IELEM) -     Y2 *U(IELEM) )*DEN
C
      A12(IELEM) = - A23(IELEM) - A31(IELEM)
      A13(IELEM) =   A12(IELEM)
      A21(IELEM) =   A23(IELEM)
      A32(IELEM) =   A31(IELEM)
C
C TERMES DIAGONAUX (LA SOMME DE CHAQUE COLONNE EST NULLE)
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
      DEN = (F3*X2 - F2*X3) * XSUR12 / SURFAC(IELEM)
C
C   TERMES EXTRADIAGONAUX
C
      A23(IELEM) = -( X3*V(IELEM) - Y3*U(IELEM) ) * DEN
      A31(IELEM) =  ( X2*V(IELEM) - Y2*U(IELEM) ) * DEN
C
      A12(IELEM) = - A23(IELEM) - A31(IELEM)
      A13(IELEM) =   A12(IELEM)
      A21(IELEM) =   A23(IELEM)
      A32(IELEM) =   A31(IELEM)
C
C TERMES DIAGONAUX (LA SOMME DE CHAQUE COLONNE EST NULLE)
C
      A11(IELEM) = - A21(IELEM) - A31(IELEM)
      A22(IELEM) = - A12(IELEM) - A32(IELEM)
      A33(IELEM) = - A13(IELEM) - A23(IELEM)
C
2     CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
200       FORMAT(1X,'MT12AA (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'MT12AA (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
          CALL PLANTE(0)
          STOP
C
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(IELMF.EQ.11.AND.IELMU.EQ.11) THEN
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
C
      V1  =  V(IKLE1(IELEM))
      V2  =  V(IKLE2(IELEM))
      V3  =  V(IKLE3(IELEM))
C
      U123 = U1 + U2 + U3
      V123 = V1 + V2 + V3
C
      DEN = (F3*Y2 - F2*Y3) * XSUR48 / SURFAC(IELEM)
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM) = ( (X2-X3)*(V123+V2) + (Y3-Y2)*(U123+U2) )*DEN
C
      A13(IELEM) = ( (X2-X3)*(V123+V3) + (Y3-Y2)*(U123+U3) )*DEN
C
      A23(IELEM) = ( X3*(V123+V3) - Y3*(U123+U3) )*DEN
C
      A21(IELEM) = ( X3*(V123+V1) - Y3*(U123+U1) )*DEN
C
      A31(IELEM) =-( X2*(V123+V1) - Y2*(U123+U1) )*DEN
C
      A32(IELEM) =-( X2*(V123+V2) - Y2*(U123+U2) )*DEN
C
C TERMES DIAGONAUX (LA SOMME DE CHAQUE COLONNE EST NULLE)
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
C
      V1  =  V(IKLE1(IELEM))
      V2  =  V(IKLE2(IELEM))
      V3  =  V(IKLE3(IELEM))
C
      U123 = U1 + U2 + U3
      V123 = V1 + V2 + V3
C
      DEN = (F3*X2 - F2*X3) * XSUR48 / SURFAC(IELEM)
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM) =-( (X2-X3)*(V123+V2) + (Y3-Y2)*(U123+U2) )*DEN
C
      A13(IELEM) =-( (X2-X3)*(V123+V3) + (Y3-Y2)*(U123+U3) )*DEN
C
      A23(IELEM) =-( X3*(V123+V3) - Y3*(U123+U3) )*DEN
C
      A21(IELEM) =-( X3*(V123+V1) - Y3*(U123+U1) )*DEN
C
      A31(IELEM) = ( X2*(V123+V1) - Y2*(U123+U1) )*DEN
C
      A32(IELEM) = ( X2*(V123+V2) - Y2*(U123+U2) )*DEN
C
C TERMES DIAGONAUX (LA SOMME DE CHAQUE COLONNE EST NULLE)
C
      A11(IELEM) = - A21(IELEM) - A31(IELEM)
      A22(IELEM) = - A12(IELEM) - A32(IELEM)
      A33(IELEM) = - A13(IELEM) - A23(IELEM)
C
4     CONTINUE
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
          CALL PLANTE(0)
          STOP
        ENDIF
C
C     AUTRES TYPES DE FONCTIONS F
C
C-----------------------------------------------------------------------
C
      ELSE
       IF (LNG.EQ.1) WRITE(LU,100) IELMF,IELMU
       IF (LNG.EQ.2) WRITE(LU,101) IELMF,IELMU
100    FORMAT(1X,'MT12AA (BIEF) :',/,
     *        1X,'COMBINAISON DE F ET U: ',1I6,2X,1I6,' NON PREVUE')
101    FORMAT(1X,'MT12AA (BIEF) :',/,
     *        1X,'COMBINATION OF F AND U: ',1I6,2X,1I6,' NOT AVAILABLE')
       CALL PLANTE(0)
       STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
