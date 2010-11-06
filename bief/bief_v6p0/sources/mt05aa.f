C                       *****************
                        SUBROUTINE MT05AA
C                       *****************
C
     *( A11 , A12 , A13 ,
     *  A21 , A22 , A23 ,
     *  A31 , A32 , A33 ,
     *  XMUL,SU,SV,U,V,
     *  XEL,YEL,IKLE,NELEM,NELMAX,FORMUL)
C
C***********************************************************************
C BIEF VERSION 5.1           05/02/91    J-M HERVOUET (LNH) 30 87 80 18
C              5.9           01/07/08     A FROEHLY (MATMECA)
C***********************************************************************
C
C  FONCTION : CALCUL DES COEFFICIENTS DE LA MATRICE SUIVANTE
C
C
C                       /           ->  --->
C       A(I,J) = XMUL  /  PSI1(I) * U . GRAD(PSI2(J)) D(OMEGA)
C                     /OMEGA
C
C
C
C  PSI1 : BASES DE TYPE TRIANGLE P1
C  PSI2 : BASES DE TYPE TRIANGLE P1
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
C-----------------------------------------------------------------------
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_MT05AA => MT05AA
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE(NELMAX,*)
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A21(*),A22(*),A23(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A31(*),A32(*),A33(*)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL,U(*),V(*)
      TYPE(BIEF_OBJ), INTENT(IN)      :: SU,SV
      CHARACTER(LEN=16) :: FORMUL
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      INTEGER IELMU,IELMV,IELEM
C
      DOUBLE PRECISION SUR24,SUR120,SUR216,SUR12
      DOUBLE PRECISION X2,X3,Y2,Y3
      DOUBLE PRECISION U1,U2,U3,U4,U5,U6,V1,V2,V3,V4,V5,V6
      DOUBLE PRECISION U123,V123,COU1,COV1,COU2,COV2,COU3,COV3
      DOUBLE PRECISION QUATRU,QUATRV,SUR6,USUR2,VSUR2
      DOUBLE PRECISION K1,K2,K3,L12,L13,L21,L23,L31,L32
C
      INTRINSIC MAX,MIN
C
C-----------------------------------------------------------------------
C
      SUR12  = XMUL/12.D0
      SUR6   = XMUL/  6.D0
      SUR24  = XMUL/ 24.D0
      SUR120 = XMUL/120.D0
      SUR216 = XMUL/216.D0
C
C-----------------------------------------------------------------------
C
      IELMU = SU%ELM
      IELMV = SV%ELM
C
C-----------------------------------------------------------------------
C
C  CAS OU U ET V SONT CONSTANTS PAR ELEMENT.
C
      IF(IELMU.EQ.10.AND.IELMV.EQ.10) THEN
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 1 IELEM = 1 , NELEM
C
      X2  =   XEL(IELEM,2) * SUR24
      X3  =   XEL(IELEM,3) * SUR24
      Y2  =   YEL(IELEM,2) * SUR24
      Y3  =   YEL(IELEM,3) * SUR24
C
      QUATRU = 4 * U(IELEM)
      QUATRV = 4 * V(IELEM)
C
C   TERMES DIAGONAUX
C
      A11(IELEM)    = (X3-X2) * QUATRV + (Y2-Y3) * QUATRU
      A22(IELEM)    = -X3     * QUATRV      +Y3  * QUATRU
      A33(IELEM)    =     X2  * QUATRV -  Y2     * QUATRU
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM)  = -X3     * QUATRV +     Y3  * QUATRU
      A13(IELEM)  =     X2  * QUATRV - Y2      * QUATRU
      A23(IELEM)  =     X2  * QUATRV - Y2      * QUATRU
      A21(IELEM)  = (X3-X2) * QUATRV + (Y2-Y3) * QUATRU
      A31(IELEM)  = (X3-X2) * QUATRV + (Y2-Y3) * QUATRU
      A32(IELEM)  = -X3     * QUATRV +     Y3  * QUATRU
C
1     CONTINUE
C
C-----------------------------------------------------------------------
C
C  CAS OU U ET V SONT LINEAIRES, QUASI-BULLE OU QUADRATIQUES ET SCHEMA N
C
      ELSEIF(FORMUL(16:16).EQ.'N'   .AND.
     * ( (IELMU.EQ.11.AND.IELMV.EQ.11).OR.
     *   (IELMU.EQ.12.AND.IELMV.EQ.12).OR.
     *   (IELMU.EQ.13.AND.IELMV.EQ.13)      )  ) THEN
C
C     SCHEMA N : U ET V TRAITES COMME LINEAIRES
C
      DO 33 IELEM = 1 , NELEM
C
         X2 = XEL(IELEM,2)
         X3 = XEL(IELEM,3)
         Y2 = YEL(IELEM,2)
         Y3 = YEL(IELEM,3)
C
         U1 = U(IKLE(IELEM,1))
         U2 = U(IKLE(IELEM,2))
         U3 = U(IKLE(IELEM,3))
         V1 = V(IKLE(IELEM,1))
         V2 = V(IKLE(IELEM,2))
         V3 = V(IKLE(IELEM,3))
C
         USUR2 = (U1+U2+U3)*SUR6
         VSUR2 = (V1+V2+V3)*SUR6
C
         K1 = USUR2 * (Y2-Y3) - VSUR2 * (X2-X3)
         K2 = USUR2 * (Y3   ) - VSUR2 * (X3   )
         K3 = USUR2 * (  -Y2) - VSUR2 * (  -X2)
C
         L12 = MAX(  MIN(K1,-K2) , 0.D0 )
         L13 = MAX(  MIN(K1,-K3) , 0.D0 )
         L21 = MAX(  MIN(K2,-K1) , 0.D0 )
         L23 = MAX(  MIN(K2,-K3) , 0.D0 )
         L31 = MAX(  MIN(K3,-K1) , 0.D0 )
         L32 = MAX(  MIN(K3,-K2) , 0.D0 )
C
C   TERMES DIAGONAUX
C
         A11(IELEM)  = L12 + L13
         A22(IELEM)  = L21 + L23
         A33(IELEM)  = L31 + L32
C
C   TERMES EXTRADIAGONAUX
C
         A12(IELEM)  = - L12
         A13(IELEM)  = - L13
         A23(IELEM)  = - L23
         A21(IELEM)  = - L21
         A31(IELEM)  = - L31
         A32(IELEM)  = - L32
C
33    CONTINUE
C
      ELSEIF(IELMU.EQ.11.AND.IELMV.EQ.11) THEN
C
C   SCHEMA CLASSIQUE, U ET V LINEAIRES
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 2 IELEM = 1 , NELEM
C
      X2  =   XEL(IELEM,2) * SUR24
      X3  =   XEL(IELEM,3) * SUR24
      Y2  =   YEL(IELEM,2) * SUR24
      Y3  =   YEL(IELEM,3) * SUR24
C
      U1   =  U(IKLE(IELEM,1))
      U2   =  U(IKLE(IELEM,2))
      U3   =  U(IKLE(IELEM,3))
      V1   =  V(IKLE(IELEM,1))
      V2   =  V(IKLE(IELEM,2))
      V3   =  V(IKLE(IELEM,3))
C
      U123 =  U1 + U2 + U3
      V123 =  V1 + V2 + V3
C
C   TERMES DIAGONAUX
C
      A11(IELEM)    = (X3-X2) * (V123+V1) + (Y2-Y3) * (U123+U1)
      A22(IELEM)    = -X3     * (V123+V2)      +Y3  * (U123+U2)
      A33(IELEM)    =     X2  * (V123+V3) -  Y2     * (U123+U3)
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM)  = -X3     * (V123+V1) +     Y3  * (U123+U1)
      A13(IELEM)  =     X2  * (V123+V1) -  Y2     * (U123+U1)
      A23(IELEM)  =     X2  * (V123+V2) -  Y2     * (U123+U2)
      A21(IELEM)  = (X3-X2) * (V123+V2) + (Y2-Y3) * (U123+U2)
      A31(IELEM)  = (X3-X2) * (V123+V3) + (Y2-Y3) * (U123+U3)
      A32(IELEM)  = -X3     * (V123+V3) +     Y3  * (U123+U3)
C
2     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSEIF(IELMU.EQ.12.AND.IELMV.EQ.12) THEN
C
C   SCHEMA CLASSIQUE, U ET V QUASI-BULLE
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 3 IELEM = 1 , NELEM
C
      X2  =   XEL(IELEM,2)
      X3  =   XEL(IELEM,3)
      Y2  =   YEL(IELEM,2)
      Y3  =   YEL(IELEM,3)
C
      U1   =  U(IKLE(IELEM,1))
      U2   =  U(IKLE(IELEM,2))
      U3   =  U(IKLE(IELEM,3))
      U4   =  U(IKLE(IELEM,4))
      V1   =  V(IKLE(IELEM,1))
      V2   =  V(IKLE(IELEM,2))
      V3   =  V(IKLE(IELEM,3))
      V4   =  V(IKLE(IELEM,4))
C
      COV1 =  5*V3+12*V4+5*V2+14*V1
      COU1 =  5*U3+12*U4+5*U2+14*U1
      COV2 =  5*V3+12*V4+14*V2+5*V1
      COU2 =  5*U3+12*U4+14*U2+5*U1
      COV3 =  14*V3+12*V4+5*V2+5*V1
      COU3 =  14*U3+12*U4+5*U2+5*U1
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM)  = ( -X3*COV1 + Y3*COU1 )*SUR216
      A13(IELEM)  = (  X2*COV1 - Y2*COU1 )*SUR216
      A21(IELEM)  = ( (X3-X2)*COV2 + (Y2-Y3)*COU2 )*SUR216
      A23(IELEM)  = (  X2*COV2 - Y2*COU2 )*SUR216
      A31(IELEM)  = ( (X3-X2)*COV3 + (Y2-Y3)*COU3 )*SUR216
      A32(IELEM)  = ( -X3*COV3 + Y3*COU3 )*SUR216
C
C   TERMES DIAGONAUX (SOMME DES COLONNES = VECTEUR NUL)
C
      A11(IELEM) = - A12(IELEM) - A13(IELEM)
      A22(IELEM) = - A21(IELEM) - A23(IELEM)
      A33(IELEM) = - A31(IELEM) - A32(IELEM)
C
3     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSEIF(IELMU.EQ.13.AND.IELMV.EQ.13) THEN
C
C   SCHEMA CLASSIQUE, U ET V P2
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 4 IELEM = 1 , NELEM
C
      X2  =   XEL(IELEM,2)
      X3  =   XEL(IELEM,3)
      Y2  =   YEL(IELEM,2)
      Y3  =   YEL(IELEM,3)
C
      U1   =  U(IKLE(IELEM,1))
      U2   =  U(IKLE(IELEM,2))
      U3   =  U(IKLE(IELEM,3))
      U4   =  U(IKLE(IELEM,4))
      U5   =  U(IKLE(IELEM,5))
      U6   =  U(IKLE(IELEM,6))
      V1   =  V(IKLE(IELEM,1))
      V2   =  V(IKLE(IELEM,2))
      V3   =  V(IKLE(IELEM,3))
      V4   =  V(IKLE(IELEM,4))
      V5   =  V(IKLE(IELEM,5))
      V6   =  V(IKLE(IELEM,6))
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM) = ((V3+V2-2.D0*V1-4.D0*V5-8.D0*(V4+V6)) * X3 
     *           -  (U2-4.D0*U5-8.D0*(U4+U6)-2.D0*U1+U3) * Y3) * SUR120
      A13(IELEM) = ((2.D0*V1-V3-V2+4.D0*V5+8.D0*(V4+V6)) * X2
     *           +  (U2-4.D0*U5-8.D0*(U4+U6)-2.D0*U1+U3) * Y2) * SUR120
      A21(IELEM) = ((V1-8.D0*(V4+V5)-4.D0*V6+V3-2.D0*V2) * (X2-X3)
     *           +  (2.D0*U2-U3-U1+8.D0*(U5+U4)+4.D0*U6) * (Y2-Y3))
     *           *   SUR120
      A23(IELEM) = ((8.D0*(V4+V5)-V1+4.D0*V6-V3+2.D0*V2) * X2
     *           -  (2.D0*U2-U3-U1+8.D0*(U5+U4)+4.D0*U6) * Y2) * SUR120
      A31(IELEM) = ((V1+V2-4.D0*V4-2.D0*V3-8.D0*(V6+V5)) * (X2-X3)
     *           +  (4.D0*U4-U2+8.D0*(U6+U5)+2.D0*U3-U1) * (Y2-Y3))
     *           *   SUR120
      A32(IELEM) = ((V1+V2-4.D0*V4-2.D0*V3-8.D0*(V6+V5)) * X3
     *           +  (4.D0*U4-U2+8.D0*(U6+U5)+2.D0*U3-U1) * Y3) * SUR120
C
C   TERMES DIAGONAUX (SOMME DES COLONNES = VECTEUR NUL)
C
      A11(IELEM) = - A12(IELEM) - A13(IELEM)
      A22(IELEM) = - A21(IELEM) - A23(IELEM)
      A33(IELEM) = - A31(IELEM) - A32(IELEM)
C
4     CONTINUE
C
C     AUTRES TYPES DE DISCRETISATION DE U ET V
C
C-----------------------------------------------------------------------
C
      ELSE
C
       IF(IELMU.EQ.IELMV) THEN
       IF (LNG.EQ.1) WRITE(LU,100) IELMU
       IF (LNG.EQ.2) WRITE(LU,101) IELMU
100    FORMAT(1X,'MT05AA (BIEF) :',/,
     *        1X,'DISCRETISATION DE U ET V : ',1I6,' NON PREVUE')
101    FORMAT(1X,'MT05AA (BIEF) :',/,
     *        1X,'DISCRETIZATION OF U AND V : ',1I6,' NOT AVAILABLE')
       ELSE
       IF (LNG.EQ.1) WRITE(LU,200) IELMU,IELMV
       IF (LNG.EQ.2) WRITE(LU,201) IELMU,IELMV
200    FORMAT(1X,'MT05AA (BIEF) :',/,
     *        1X,'U ET V DE DISCRETISATIONS DIFFERENTES :',1I6,3X,1I6)
201    FORMAT(1X,'MT05AA (BIEF) :',/,
     *        1X,'U AND V OF A DIFFERENT DISCRETISATION:',1I6,3X,1I6)
       ENDIF
C
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
