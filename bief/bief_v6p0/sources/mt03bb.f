C                       *****************
                        SUBROUTINE MT03BB
C                       *****************
C
     *( A11 , A12 , A13 , A14 ,
     *  A21 , A22 , A23 , A24 ,
     *  A31 , A32 , A33 , A34 ,
     *  A41 , A42 , A43 , A44 ,
     *  XMUL,SF,SG,SU,SV,F,G,U,V,
     *  XEL,YEL,IKLE1,IKLE2,IKLE3,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.1           12/04/93    J-M HERVOUET (LNH) 30 87 80 18
C                                            C MOULIN (LNH) 30 87 83 81
C***********************************************************************
C
C  FONCTION : CALCUL DES COEFFICIENTS DE LA MATRICE SUIVANTE
C
C
C                 /  -->   - -->          -->  --->
C A(I,J) = XMUL  /   KEL . GRAD(PSI1(I)) * U . GRAD(PSI2(J)) D(OMEGA)
C               /OMEGA
C         -->
C         KEL VECTEUR CONSTANT SUR L'ELEMENT, DE COMPOSANTES F ET G
C
C         PSI1 DE DISCRETISATION P1
C         PSI2 DE DISCRETISATION P2
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
C |     IKLE1..4   | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |     NELEM      | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |     NELMAX     | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : NEANT
C
C-----------------------------------------------------------------------
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_MT03BB => MT03BB
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
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*),A14(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A21(*),A22(*),A23(*),A24(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A31(*),A32(*),A33(*),A34(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A41(*),A42(*),A43(*),A44(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: F(*),G(*),U(*),V(*)
C
C     STRUCTURES DE      F, G, U, V
      TYPE(BIEF_OBJ), INTENT(IN) :: SF,SG,SU,SV
C
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      INTEGER IELEM,IELMF,IELMG,IELMU,IELMV
C
      DOUBLE PRECISION X2,X3,Y2,Y3,U1,U2,U3,U4,V1,V2,V3,V4
      DOUBLE PRECISION KXEL,KYEL
C
C-----------------------------------------------------------------------
C
      IELMF = SF%ELM
      IELMG = SG%ELM
      IELMU = SU%ELM
      IELMV = SV%ELM
C
C-----------------------------------------------------------------------
C CAS OU U EST DE DISCRETISATION P1
C-----------------------------------------------------------------------
C
      IF(IELMF.EQ.10.AND.IELMG.EQ.10.AND.
     *   IELMU.EQ.11.AND.IELMV.EQ.11) THEN
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 1 IELEM = 1 , NELEM
C
      X2  =   XEL(IELEM,2)
      X3  =   XEL(IELEM,3)
      Y2  =   YEL(IELEM,2)
      Y3  =   YEL(IELEM,3)
C
      U1   =  U(IKLE1(IELEM))
      U2   =  U(IKLE2(IELEM))
      U3   =  U(IKLE3(IELEM))
      V1   =  V(IKLE1(IELEM))
      V2   =  V(IKLE2(IELEM))
      V3   =  V(IKLE3(IELEM))
C
      KXEL =  F(IELEM)
      KYEL =  G(IELEM)
C
C ON CALCULE 9 DES 16 TERMES
C
      A11(IELEM)=
     *  (X2**2*KYEL*(8*V3+17*V2+20*V1)+4*X2*X3*KYEL*(-
     *  5*V3-5*V2-8*V1)+X2*Y2*KXEL*(-8*V3-17*V2-20*V1)+X2*Y2*
     *  KYEL*(-8*U3-17*U2-20*U1)+2*X2*Y3*KXEL*(5*V3+5*V2+8*V1)
     *  +2*X2*Y3*KYEL*(5*U3+5*U2+8*U1)+X3**2*KYEL*(17*V3+8*V2+
     *  20*V1)+2*X3*Y2*KXEL*(5*V3+5*V2+8*V1)+2*X3*Y2*KYEL*(5*U3
     *  +5*U2+8*U1)+X3*Y3*KXEL*(-17*V3-8*V2-20*V1)+X3*Y3*KYEL*(-
     *  17*U3-8*U2-20*U1)+Y2**2*KXEL*(8*U3+17*U2+20*U1)+4*Y2*
     *  Y3*KXEL*(-5*U3-5*U2-8*U1)+Y3**2*KXEL*(17*U3+8*U2+20*U1)
     *  )*XMUL/(54*X2*Y3-54*X3*Y2)
C
      A12(IELEM)=
     * (2*X2**2*KYEL*(V3+4*V2+4*V1)+X2*X3*KYEL*(V3+4*
     *  V2+4*V1)+2*X2*Y2*KXEL*(-V3-4*V2-4*V1)+2*X2*Y2*KYEL*(-U3-
     *  4*U2-4*U1)+X2*Y3*KXEL*(V3+4*V2+4*V1)+2*X2*Y3*KYEL*(-U3-4
     *  *U2-4*U1)+X3**2*KYEL*(-V3-4*V2-4*V1)+2*X3*Y2*KXEL*(-V3-4
     *  *V2-4*V1)+X3*Y2*KYEL*(U3+4*U2+4*U1)+X3*Y3*KXEL*(V3+4*V2+
     *  4*V1)+X3*Y3*KYEL*(U3+4*U2+4*U1)+2*Y2**2*KXEL*(U3+4*U2+4*
     *  U1)+Y2*Y3*KXEL*(U3+4*U2+4*U1)+Y3**2*KXEL*(-U3-4*U2-4*U1))
     *  *XMUL/(54*X2*Y3-54*X3*Y2)
C
      A13(IELEM)=
     *  (X2**2*KYEL*(-4*V3-V2-4*V1)+X2*X3*KYEL*(4*V3+V2+
     *  4*V1)+X2*Y2*KXEL*(4*V3+V2+4*V1)+X2*Y2*KYEL*(4*U3+U2+4*U1)
     *  +2*X2*Y3*KXEL*(-4*V3-V2-4*V1)+X2*Y3*KYEL*(4*U3+U2+4*U1)+
     *  2*X3**2*KYEL*(4*V3+V2+4*V1)+X3*Y2*KXEL*(4*V3+V2+4*V1)+2*
     *  X3*Y2*KYEL*(-4*U3-U2-4*U1)+2*X3*Y3*KXEL*(-4*V3-V2-4*V1)+
     *  2*X3*Y3*KYEL*(-4*U3-U2-4*U1)+Y2**2*KXEL*(-4*U3-U2-4*U1)+
     *  Y2*Y3*KXEL*(4*U3+U2+4*U1)+2*Y3**2*KXEL*(4*U3+U2+4*U1))
     *  *XMUL/(54*X2*Y3-54*X3*Y2)
C
      A21(IELEM)=
     *  (2*X2**2*KYEL*(V3+4*V2+4*V1)+X2*X3*KYEL*(V3+4*
     *  V2+4*V1)+2*X2*Y2*KXEL*(-V3-4*V2-4*V1)+2*X2*Y2*KYEL*(-U3-
     *  4*U2-4*U1)+2*X2*Y3*KXEL*(-V3-4*V2-4*V1)+X2*Y3*KYEL*(U3+4
     *  *U2+4*U1)+X3**2*KYEL*(-V3-4*V2-4*V1)+X3*Y2*KXEL*(V3+4*V2+
     *  4*V1)+2*X3*Y2*KYEL*(-U3-4*U2-4*U1)+X3*Y3*KXEL*(V3+4*V2+4
     *  *V1)+X3*Y3*KYEL*(U3+4*U2+4*U1)+2*Y2**2*KXEL*(U3+4*U2+4*
     *  U1)+Y2*Y3*KXEL*(U3+4*U2+4*U1)+Y3**2*KXEL*(-U3-4*U2-4*U1))
     *  *XMUL/(54*X2*Y3-54*X3*Y2)
C
      A23(IELEM)=
     *  (2*X2**2*KYEL*(4*V3+4*V2+V1)+5*X2*X3*KYEL*(-4*
     *  V3-4*V2-V1)+2*X2*Y2*KXEL*(-4*V3-4*V2-V1)+2*X2*Y2*KYEL*(-
     *  4*U3-4*U2-U1)+4*X2*Y3*KXEL*(4*V3+4*V2+V1)+X2*Y3*KYEL*(4*
     *  U3+4*U2+U1)+2*X3**2*KYEL*(4*V3+4*V2+V1)+X3*Y2*KXEL*(4*V3
     *  +4*V2+V1)+4*X3*Y2*KYEL*(4*U3+4*U2+U1)+2*X3*Y3*KXEL*(-4*
     *  V3-4*V2-V1)+2*X3*Y3*KYEL*(-4*U3-4*U2-U1)+2*Y2**2*KXEL*(
     *  4*U3+4*U2+U1)+5*Y2*Y3*KXEL*(-4*U3-4*U2-U1)+2*Y3**2*KXEL*
     *  (4*U3+4*U2+U1))*XMUL/(54*X2*Y3-54*X3*Y2)
C
      A24(IELEM)=
     *  (X2**2*KYEL*(-5*V3-8*V2-5*V1)+X2*X3*KYEL*(11*V3
     *  +8*V2-V1)+X2*Y2*KXEL*(5*V3+8*V2+5*V1)+X2*Y2*KYEL*(5*U3+
     *  8*U2+5*U1)+X2*Y3*KXEL*(-7*V3-4*V2+2*V1)+X2*Y3*KYEL*(-4*
     *  U3-4*U2-U1)+2*X3**2*KYEL*(-4*V3-4*V2-V1)+X3*Y2*KXEL*(-4*
     *  V3-4*V2-V1)+X3*Y2*KYEL*(-7*U3-4*U2+2*U1)+2*X3*Y3*KXEL*(
     *  4*V3+4*V2+V1)+2*X3*Y3*KYEL*(4*U3+4*U2+U1)+Y2**2*KXEL*(-5
     *  *U3-8*U2-5*U1)+Y2*Y3*KXEL*(11*U3+8*U2-U1)+2*Y3**2*KXEL*(
     *  -4*U3-4*U2-U1))*XMUL/(18*X2*Y3-18*X3*Y2)
C
      A31(IELEM)=
     *  (X2**2*KYEL*(-4*V3-V2-4*V1)+X2*X3*KYEL*(4*V3+V2+
     *  4*V1)+X2*Y2*KXEL*(4*V3+V2+4*V1)+X2*Y2*KYEL*(4*U3+U2+4*U1)
     *  +X2*Y3*KXEL*(4*V3+V2+4*V1)+2*X2*Y3*KYEL*(-4*U3-U2-4*U1)+
     *  2*X3**2*KYEL*(4*V3+V2+4*V1)+2*X3*Y2*KXEL*(-4*V3-V2-4*V1)
     *  +X3*Y2*KYEL*(4*U3+U2+4*U1)+2*X3*Y3*KXEL*(-4*V3-V2-4*V1)+
     *  2*X3*Y3*KYEL*(-4*U3-U2-4*U1)+Y2**2*KXEL*(-4*U3-U2-4*U1)+
     *  Y2*Y3*KXEL*(4*U3+U2+4*U1)+2*Y3**2*KXEL*(4*U3+U2+4*U1))
     *  *XMUL/(54*X2*Y3-54*X3*Y2)
C
      A32(IELEM)=
     *  (2*X2**2*KYEL*(4*V3+4*V2+V1)+5*X2*X3*KYEL*(-4*
     *  V3-4*V2-V1)+2*X2*Y2*KXEL*(-4*V3-4*V2-V1)+2*X2*Y2*KYEL*(-
     *  4*U3-4*U2-U1)+X2*Y3*KXEL*(4*V3+4*V2+V1)+4*X2*Y3*KYEL*(4*
     *  U3+4*U2+U1)+2*X3**2*KYEL*(4*V3+4*V2+V1)+4*X3*Y2*KXEL*(4
     *  *V3+4*V2+V1)+X3*Y2*KYEL*(4*U3+4*U2+U1)+2*X3*Y3*KXEL*(-4*
     *  V3-4*V2-V1)+2*X3*Y3*KYEL*(-4*U3-4*U2-U1)+2*Y2**2*KXEL*(
     *  4*U3+4*U2+U1)+5*Y2*Y3*KXEL*(-4*U3-4*U2-U1)+2*Y3**2*KXEL*
     *  (4*U3+4*U2+U1))*XMUL/(54*X2*Y3-54*X3*Y2)
C
      A34(IELEM)=
     *  (2*X2**2*KYEL*(-4*V3-4*V2-V1)+X2*X3*KYEL*(8*V3+
     *  11*V2-V1)+2*X2*Y2*KXEL*(4*V3+4*V2+V1)+2*X2*Y2*KYEL*(4*U3
     *  +4*U2+U1)+X2*Y3*KXEL*(-4*V3-4*V2-V1)+X2*Y3*KYEL*(-4*U3-7
     *  *U2+2*U1)+X3**2*KYEL*(-8*V3-5*V2-5*V1)+X3*Y2*KXEL*(-4*V3
     *  -7*V2+2*V1)+X3*Y2*KYEL*(-4*U3-4*U2-U1)+X3*Y3*KXEL*(8*V3+
     *  5*V2+5*V1)+X3*Y3*KYEL*(8*U3+5*U2+5*U1)+2*Y2**2*KXEL*(-4
     *  *U3-4*U2-U1)+Y2*Y3*KXEL*(8*U3+11*U2-U1)+Y3**2*KXEL*(-8*U3
     *  -5*U2-5*U1))*XMUL/(18*X2*Y3-18*X3*Y2)
C
C  ON UTILISE ICI LES PROPRIETES DE "CARRE MAGIQUE" D'UNE MATRICE
C  DE TYPE DIFFUSION (SOMME DE CHAQUE LIGNE ET CHAQUE COLONNE NULLE)
C
      A14(IELEM) = - A11(IELEM) - A12(IELEM) - A13(IELEM)
      A22(IELEM) = - A21(IELEM) - A23(IELEM) - A24(IELEM)
      A33(IELEM) = - A31(IELEM) - A32(IELEM) - A34(IELEM)
      A41(IELEM) = - A11(IELEM) - A21(IELEM) - A31(IELEM)
      A42(IELEM) = - A12(IELEM) - A22(IELEM) - A32(IELEM)
      A43(IELEM) = - A13(IELEM) - A23(IELEM) - A33(IELEM)
      A44(IELEM) = - A14(IELEM) - A24(IELEM) - A34(IELEM)
C
1     CONTINUE
C
C-----------------------------------------------------------------------
C CAS OU U EST DE DISCRETISATION QUASI-BULLE
C-----------------------------------------------------------------------
C
      ELSEIF(IELMF.EQ.10.AND.IELMG.EQ.10.AND.
     *       IELMU.EQ.12.AND.IELMV.EQ.12) THEN
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 2 IELEM = 1 , NELEM
C
      X2  =   XEL(IELEM,2)
      X3  =   XEL(IELEM,3)
      Y2  =   YEL(IELEM,2)
      Y3  =   YEL(IELEM,3)
C
      U1   =  U(IKLE1(IELEM))
      U2   =  U(IKLE2(IELEM))
      U3   =  U(IKLE3(IELEM))
      U4   =  U(IKLE3(IELEM))
      V1   =  V(IKLE1(IELEM))
      V2   =  V(IKLE2(IELEM))
      V3   =  V(IKLE3(IELEM))
      V4   =  V(IKLE3(IELEM))
C
      KXEL =  F(IELEM)
      KYEL =  G(IELEM)
C
C ON CALCULE 9 DES 16 TERMES
C
      A11(IELEM)=
     *  (X2**2*KYEL*(V3+5*V4+4*V2+5*V1)+4*X2*X3*KYEL*(-
     *  V3-2*V4-V2-2*V1)+X2*Y2*KXEL*(-V3-5*V4-4*V2-5*V1)+X2*Y2
     *  *KYEL*(-U3-5*U4-4*U2-5*U1)+2*X2*Y3*KXEL*(V3+2*V4+V2+2*
     *  V1)+2*X2*Y3*KYEL*(U3+2*U4+U2+2*U1)+X3**2*KYEL*(4*V3+5*V4
     *  +V2+5*V1)+2*X3*Y2*KXEL*(V3+2*V4+V2+2*V1)+2*X3*Y2*KYEL*(
     *  U3+2*U4+U2+2*U1)+X3*Y3*KXEL*(-4*V3-5*V4-V2-5*V1)+X3*Y3
     *  *KYEL*(-4*U3-5*U4-U2-5*U1)+Y2**2*KXEL*(U3+5*U4+4*U2+5*
     *  U1)+4*Y2*Y3*KXEL*(-U3-2*U4-U2-2*U1)+Y3**2*KXEL*(4*U3+5*
     *  U4+U2+5*U1))*XMUL/(18*X2*Y3-18*X3*Y2)
C
      A12(IELEM)=
     *  (2*X2**2*KYEL*(V4+V2+V1)+X2*X3*KYEL*(V4+V2+V1)-(2
     *  *X2*Y2*KXEL)*(V4+V2+V1)-(2*X2*Y2*KYEL)*(U4+U2+U1)+X2*Y3*KXEL*(
     *  V4+V2+V1)-(2*X2*Y3*KYEL)*(U4+U2+U1)-(X3**2*KYEL)*(V4+V2+V1)-
     *  (2*X3*Y2*KXEL)*(V4+V2+V1)+X3*Y2*KYEL*(U4+U2+U1)+X3*Y3*KXEL*(V4
     *  +V2+V1)+X3*Y3*KYEL*(U4+U2+U1)+2*Y2**2*KXEL*(U4+U2+U1)+Y2*Y3*
     *  KXEL*(U4+U2+U1)-(Y3**2*KXEL)*(U4+U2+U1))
     *  *XMUL/(18*X2*Y3-18*X3*Y2)
C
      A13(IELEM)=
     *  (-(X2**2*KYEL)*(V3+V4+V1)+X2*X3*KYEL*(V3+V4+V1)+X2*
     *  Y2*KXEL*(V3+V4+V1)+X2*Y2*KYEL*(U3+U4+U1)-(2*X2*Y3*KXEL)*(V3+V4
     *  +V1)+X2*Y3*KYEL*(U3+U4+U1)+2*X3**2*KYEL*(V3+V4+V1)+X3*Y2*KXEL*
     *  (V3+V4+V1)-(2*X3*Y2*KYEL)*(U3+U4+U1)-(2*X3*Y3*KXEL)*(V3+V4+
     *  V1)-(2*X3*Y3*KYEL)*(U3+U4+U1)-(Y2**2*KXEL)*(U3+U4+U1)+Y2*Y3*
     *  KXEL*(U3+U4+U1)+2*Y3**2*KXEL*(U3+U4+U1))
     *  *XMUL/(18*X2*Y3-18*X3*Y2)
C
      A21(IELEM)=
     *  (2*X2**2*KYEL*(V4+V2+V1)+X2*X3*KYEL*(V4+V2+V1)-(2
     *  *X2*Y2*KXEL)*(V4+V2+V1)-(2*X2*Y2*KYEL)*(U4+U2+U1)-(2*X2*Y3*
     *  KXEL)*(V4+V2+V1)+X2*Y3*KYEL*(U4+U2+U1)-(X3**2*KYEL)*(V4+V2+V1)+
     *  X3*Y2*KXEL*(V4+V2+V1)-(2*X3*Y2*KYEL)*(U4+U2+U1)+X3*Y3*KXEL*(V4
     *  +V2+V1)+X3*Y3*KYEL*(U4+U2+U1)+2*Y2**2*KXEL*(U4+U2+U1)+Y2*Y3*
     *  KXEL*(U4+U2+U1)-(Y3**2*KXEL)*(U4+U2+U1))
     *  *XMUL/(18*X2*Y3-18*X3*Y2)
C
      A23(IELEM)=
     *  (2*X2**2*KYEL*(V3+V4+V2)-(5*X2*X3*KYEL)*(V3+V4+V2
     *  )-(2*X2*Y2*KXEL)*(V3+V4+V2)-(2*X2*Y2*KYEL)*(U3+U4+U2)+4*X2
     *  *Y3*KXEL*(V3+V4+V2)+X2*Y3*KYEL*(U3+U4+U2)+2*X3**2*KYEL*(V3+V4+
     *  V2)+X3*Y2*KXEL*(V3+V4+V2)+4*X3*Y2*KYEL*(U3+U4+U2)-(2*X3*Y3*
     *  KXEL)*(V3+V4+V2)-(2*X3*Y3*KYEL)*(U3+U4+U2)+2*Y2**2*KXEL*(U3+
     *  U4+U2)-(5*Y2*Y3*KXEL)*(U3+U4+U2)+2*Y3**2*KXEL*(U3+U4+U2))
     *  *XMUL/(18*X2*Y3-18*X3*Y2)
C
      A24(IELEM)=
     *  (X2**2*KYEL*(-V3-2*V4-2*V2-V1)+X2*X3*KYEL*(3*V3+
     *  2*V4+2*V2-V1)+X2*Y2*KXEL*(V3+2*V4+2*V2+V1)+X2*Y2*KYEL*(U3+
     *  2*U4+2*U2+U1)+X2*Y3*KXEL*(-2*V3-V4-V2+V1)-(X2*Y3*KYEL)*(U3+
     *  U4+U2)-(2*X3**2*KYEL)*(V3+V4+V2)-(X3*Y2*KXEL)*(V3+V4+V2)+X3*
     *  Y2*KYEL*(-2*U3-U4-U2+U1)+2*X3*Y3*KXEL*(V3+V4+V2)+2*X3*Y3*
     *  KYEL*(U3+U4+U2)+Y2**2*KXEL*(-U3-2*U4-2*U2-U1)+Y2*Y3*KXEL*(3*
     *  U3+2*U4+2*U2-U1)-(2*Y3**2*KXEL)*(U3+U4+U2))
     *  *XMUL/(6*X2*Y3-6*X3*Y2)
C
      A31(IELEM)=
     *  (-(X2**2*KYEL)*(V3+V4+V1)+X2*X3*KYEL*(V3+V4+V1)+X2*
     *  Y2*KXEL*(V3+V4+V1)+X2*Y2*KYEL*(U3+U4+U1)+X2*Y3*KXEL*(V3+V4+V1)-
     *  (2*X2*Y3*KYEL)*(U3+U4+U1)+2*X3**2*KYEL*(V3+V4+V1)-(2*X3*Y2
     *  *KXEL)*(V3+V4+V1)+X3*Y2*KYEL*(U3+U4+U1)-(2*X3*Y3*KXEL)*(V3+V4+
     *  V1)-(2*X3*Y3*KYEL)*(U3+U4+U1)-(Y2**2*KXEL)*(U3+U4+U1)+Y2*Y3*
     *  KXEL*(U3+U4+U1)+2*Y3**2*KXEL*(U3+U4+U1))
     *  *XMUL/(18*X2*Y3-18*X3*Y2)
C
      A32(IELEM)=
     *  (2*X2**2*KYEL*(V3+V4+V2)-(5*X2*X3*KYEL)*(V3+V4+V2
     *  )-(2*X2*Y2*KXEL)*(V3+V4+V2)-(2*X2*Y2*KYEL)*(U3+U4+U2)+X2*Y3
     *  *KXEL*(V3+V4+V2)+4*X2*Y3*KYEL*(U3+U4+U2)+2*X3**2*KYEL*(V3+V4+
     *  V2)+4*X3*Y2*KXEL*(V3+V4+V2)+X3*Y2*KYEL*(U3+U4+U2)-(2*X3*Y3*
     *  KXEL)*(V3+V4+V2)-(2*X3*Y3*KYEL)*(U3+U4+U2)+2*Y2**2*KXEL*(U3+
     *  U4+U2)-(5*Y2*Y3*KXEL)*(U3+U4+U2)+2*Y3**2*KXEL*(U3+U4+U2))
     *  *XMUL/(18*X2*Y3-18*X3*Y2)
C
      A34(IELEM)=
     *  (-(2*X2**2*KYEL)*(V3+V4+V2)+X2*X3*KYEL*(2*V3+2*
     *  V4+3*V2-V1)+2*X2*Y2*KXEL*(V3+V4+V2)+2*X2*Y2*KYEL*(U3+U4+U2
     *  )-(X2*Y3*KXEL)*(V3+V4+V2)+X2*Y3*KYEL*(-U3-U4-2*U2+U1)+X3**2*
     *  KYEL*(-2*V3-2*V4-V2-V1)+X3*Y2*KXEL*(-V3-V4-2*V2+V1)-(X3*Y2
     *  *KYEL)*(U3+U4+U2)+X3*Y3*KXEL*(2*V3+2*V4+V2+V1)+X3*Y3*KYEL*(2
     *  *U3+2*U4+U2+U1)-(2*Y2**2*KXEL)*(U3+U4+U2)+Y2*Y3*KXEL*(2*U3
     *  +2*U4+3*U2-U1)+Y3**2*KXEL*(-2*U3-2*U4-U2-U1))
     *  *XMUL/(6*X2*Y3-6*X3*Y2)
C
C
C  ON UTILISE ICI LES PROPRIETES DE "CARRE MAGIQUE" D'UNE MATRICE
C  DE TYPE DIFFUSION (SOMME DE CHAQUE LIGNE ET CHAQUE COLONNE NULLE)
C
      A14(IELEM) = - A11(IELEM) - A12(IELEM) - A13(IELEM)
      A22(IELEM) = - A21(IELEM) - A23(IELEM) - A24(IELEM)
      A33(IELEM) = - A31(IELEM) - A32(IELEM) - A34(IELEM)
      A41(IELEM) = - A11(IELEM) - A21(IELEM) - A31(IELEM)
      A42(IELEM) = - A12(IELEM) - A22(IELEM) - A32(IELEM)
      A43(IELEM) = - A13(IELEM) - A23(IELEM) - A33(IELEM)
      A44(IELEM) = - A14(IELEM) - A24(IELEM) - A34(IELEM)
C
2     CONTINUE
C
C     AUTRES TYPES DE FONCTIONS F ET G
C
C-----------------------------------------------------------------------
C
      ELSE
C
       IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
       IF (LNG.EQ.1) WRITE(LU,110) IELMG,SG%NAME
       IF (LNG.EQ.1) WRITE(LU,200) IELMU,SU%NAME
       IF (LNG.EQ.1) WRITE(LU,300)
       IF (LNG.EQ.1) WRITE(LU,101) IELMF,SF%NAME
       IF (LNG.EQ.1) WRITE(LU,111) IELMG,SG%NAME
       IF (LNG.EQ.1) WRITE(LU,201) IELMU,SU%NAME
       IF (LNG.EQ.1) WRITE(LU,301)
100    FORMAT(1X,'MT03BB (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,
     *        1X,'NOM REEL : ',A6)
110    FORMAT(1X,'DISCRETISATION DE G : ',1I6,
     *        1X,'NOM REEL : ',A6)
200    FORMAT(1X,'DISCRETISATION DE U : ',1I6,
     *        1X,'NOM REEL : ',A6)
300    FORMAT(1X,'CAS NON PREVU')
101    FORMAT(1X,'MT03BB (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F:',1I6,
     *        1X,'REAL NAME: ',A6)
111    FORMAT(1X,'DISCRETIZATION OF G:',1I6,
     *        1X,'REAL NAME: ',A6)
201    FORMAT(1X,'DISCRETIZATION OF U:',1I6,
     *        1X,'REAL NAME: ',A6)
301    FORMAT(1X,'CASE NOT IMPLEMENTED')
       CALL PLANTE(0)
       STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
