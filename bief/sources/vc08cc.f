C                       *****************
                        SUBROUTINE VC08CC
C                       *****************
C
     *( XMUL,SF,SU,SV,F,U,V,XEL,YEL,
     *  IKLE1,IKLE2,IKLE3,IKLE4,IKLE5,IKLE6,NELEM,NELMAX,
     *  W1,W2,W3,W4, W5,W6,FORMUL )
C
C***********************************************************************
C BIEF VERSION 5.9        01/07/08    A FROEHLY
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C
C                    /                  DF      DF
C      V  =  XMUL   /       PSII  * ( U --  + V -- )   D(OMEGA)
C       I          /OMEGA               DX      DY
C
C
C    PSI(I) EST UNE BASE DE TYPE P2
C
C    ATTENTION : LE RESULTAT EST DANS W  SOUS FORME NON ASSEMBLEE.
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
C**********************************************************************
C
      USE BIEF, EX_VC08CC => VC08CC
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX*3),YEL(NELMAX*3),XMUL
C     W1 IS ALSO USED AS 1-DIMENSIONAL FOR ALL W
      DOUBLE PRECISION, INTENT(INOUT) :: W1(6*NELMAX),W2(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W3(NELMAX),W4(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W5(NELMAX),W6(NELMAX)
C     IKLE1 IS ALSO USED AS A 1-DIMENSIONAL IKLE
      INTEGER, INTENT(IN) :: IKLE1(6*NELMAX),IKLE2(NELMAX)
      INTEGER, INTENT(IN) :: IKLE3(NELMAX),IKLE4(NELMAX)
      INTEGER, INTENT(IN) :: IKLE5(NELMAX),IKLE6(NELMAX)
      CHARACTER(LEN=16), INTENT(IN) :: FORMUL
C
C     STRUCTURES DE F,G,H,U,V,W ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF,SU,SV
      DOUBLE PRECISION, INTENT(IN) :: F(*),U(*),V(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF,IELMU,IELMV
C
      DOUBLE PRECISION K1,K2,K3
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION X2,Y2,X3,Y3,F1,F2,F3,F4,F5,F6
      DOUBLE PRECISION U1,U2,U3,U4,U5,U6,V1,V2,V3,V4,V5,V6
      DOUBLE PRECISION ANS1,SUR6
      DOUBLE PRECISION PHIT,USUR2,VSUR2,XSU90,XSU360,XSU630,XSU2520
      DOUBLE PRECISION L12,L13,L21,L23,L31,L32,BETAN1,BETAN2,BETAN3
      INTRINSIC MAX,MIN
C
C-----------------------------------------------------------------------
C
      SUR6   = 1.D0 / 6.D0
      XSU90  = XMUL/90.D0
      XSU360 = XMUL/360.D0
      XSU630 = XMUL/630.D0
      XSU2520= XMUL/2520.D0
C
      IELMF=SF%ELM
      IELMU=SU%ELM
      IELMV=SV%ELM
C
C-----------------------------------------------------------------------
C
C     FONCTION F ET VECTEUR U P2
C
      IF(IELMF.EQ.13.AND.IELMU.EQ.13.AND.IELMV.EQ.13) THEN
C
      IF(FORMUL(14:16).EQ.'PSI') THEN
C
C     SCHEMA PSI P1 PUIS INTERPOLATION LINEAIRE
C
      DO 33 IELEM = 1 , NELEM
C
         X2 = XEL(IELEM+NELMAX)
         X3 = XEL(IELEM+2*NELMAX)
         Y2 = YEL(IELEM+NELMAX)
         Y3 = YEL(IELEM+2*NELMAX)
C
         F1 = F(IKLE1(IELEM))
         F2 = F(IKLE2(IELEM))
         F3 = F(IKLE3(IELEM))
C
         U1 = U(IKLE1(IELEM))
         U2 = U(IKLE2(IELEM))
         U3 = U(IKLE3(IELEM))
         V1 = V(IKLE1(IELEM))
         V2 = V(IKLE2(IELEM))
         V3 = V(IKLE3(IELEM))
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
         BETAN1 = L12*(F1-F2) + L13*(F1-F3)
         BETAN2 = L21*(F2-F1) + L23*(F2-F3)
         BETAN3 = L31*(F3-F1) + L32*(F3-F2)
C
         PHIT = BETAN1 + BETAN2 + BETAN3
C
         IF(PHIT.GT.0.D0) THEN
           W1(IELEM) =   XMUL * MAX( MIN( BETAN1, PHIT),0.D0 )
           W2(IELEM) =   XMUL * MAX( MIN( BETAN2, PHIT),0.D0 )
           W3(IELEM) =   XMUL * MAX( MIN( BETAN3, PHIT),0.D0 )
         ELSE
           W1(IELEM) = - XMUL * MAX( MIN(-BETAN1,-PHIT),0.D0 )
           W2(IELEM) = - XMUL * MAX( MIN(-BETAN2,-PHIT),0.D0 )
           W3(IELEM) = - XMUL * MAX( MIN(-BETAN3,-PHIT),0.D0 )
         ENDIF
         W4(IELEM) =   (W1(IELEM)+ W2(IELEM))/2.D0
         W5(IELEM) =   (W2(IELEM)+ W3(IELEM))/2.D0
         W6(IELEM) =   (W3(IELEM)+ W1(IELEM))/2.D0
C
33    CONTINUE
C
      ELSE
C
C     CALCUL CLASSIQUE
C
      DO 3 IELEM = 1 , NELEM
C
         X2 = XEL(IELEM+NELMAX)
         X3 = XEL(IELEM+2*NELMAX)
         Y2 = YEL(IELEM+NELMAX)
         Y3 = YEL(IELEM+2*NELMAX)
C
         U1 = U(IKLE1(IELEM))
         U2 = U(IKLE2(IELEM))
         U3 = U(IKLE3(IELEM))
         U4 = U(IKLE4(IELEM))
         U5 = U(IKLE5(IELEM))
         U6 = U(IKLE6(IELEM))
C
         V1 = V(IKLE1(IELEM))
         V2 = V(IKLE2(IELEM))
         V3 = V(IKLE3(IELEM))
         V4 = V(IKLE4(IELEM))
         V5 = V(IKLE5(IELEM))
         V6 = V(IKLE6(IELEM))
C
         F1 = F(IKLE1(IELEM))
         F2 = F(IKLE2(IELEM)) - F1
         F3 = F(IKLE3(IELEM)) - F1
         F4 = F(IKLE4(IELEM)) - F1
         F5 = F(IKLE5(IELEM)) - F1
         F6 = F(IKLE6(IELEM)) - F1
C
      ANS1 =-20.D0*Y2*F3*U5-64.D0*X2*F6*V4+16.D0*Y2*F4*U6+
     *       16.D0*X3*F6*V4+4.D0*Y2*F5*U3-32.D0*Y3*F4*U5-16.D0*X2*F4*V6-
     *       16.D0*Y2*F3*U4+32.D0*X3*F4*V5-32.D0*X2*F4*V4-
     *       24.D0*Y2*F4*U1+24.D0*Y2*F5*U1-80.D0*Y3*F4*U4-
     *       24.D0*X2*F5*V1+16.D0*X2*F3*V4-9.D0*X3*F2*V2+
     *       80.D0*X3*F4*V4-16.D0*Y3*F6*U4+9.D0*X2*F3*V3+
     *       4.D0*Y3*F6*U2-48.D0*X2*F4*V5-32.D0*Y3*F6*U6-
     *       32.D0*Y2*F5*U4+32.D0*Y3*F2*U4-16.D0*Y2*F5*U6+
     *       11.D0*X3*F2*V3+24.D0*Y3*F6*U1+16.D0*X2*F5*V6-
     *       64.D0*Y3*F4*U6-48.D0*Y3*F6*U5+4.D0*X2*F4*V3-
     *       48.D0*X3*F5*V5-18.D0*X3*F2*V1+16.D0*Y3*F2*U6+
     *       32.D0*Y2*F4*U4+20.D0*Y3*F4*U3+11.D0*Y2*F3*U2+
     *       20.D0*Y3*F2*U5-4.D0*X3*F6*V2+32.D0*X2*F5*V4-11.D0*X2*F3*V2-
     *       32.D0*X3*F2*V4-9.D0*Y2*F3*U3+96.D0*Y2*F6*U1+9.D0*Y3*F2*U2-
     *       24.D0*X3*F6*V1+32.D0*Y3*F5*U6-24.D0*Y3*F5*U1+
     *       48.D0*Y2*F4*U5+96.D0*X3*F4*V1-48.D0*Y2*F5*U5+
     *       48.D0*X2*F5*V5-32.D0*X3*F5*V6-20.D0*X3*F2*V5+4.D0*X3*F5*V2+
     *       64.D0*X3*F4*V6-16.D0*X2*F4*V2+18.D0*X2*F3*V1+
     *       24.D0*X3*F5*V1-96.D0*X2*F6*V1+32.D0*X3*F6*V6-
     *       11.D0*Y3*F2*U3-4.D0*X2*F5*V3+20.D0*X2*F3*V5+16.D0*Y3*F5*U3
      W1(IELEM) =( -16.D0*Y3*F6*U3-20.D0*X3*F4*V3-16.D0*X3*F5*V3-
     *              16.D0*Y2*F5*U2+16.D0*Y2*F4*U2+16.D0*X2*F5*V2-
     *              20.D0*Y2*F6*U2+32.D0*Y2*F6*U5+20.D0*X2*F6*V2+
     *              80.D0*Y2*F6*U6-96.D0*Y3*F4*U1-80.D0*X2*F6*V6-
     *              32.D0*X2*F6*V5-18.D0*Y2*F3*U1-4.D0*Y2*F4*U3+
     *              32.D0*X2*F3*V6+16.D0*X3*F6*V3-32.D0*Y2*F3*U6-
     *              16.D0*X3*F2*V6-4.D0*Y3*F5*U2+24.D0*X2*F4*V1+
     *              64.D0*Y2*F6*U4+16.D0*Y3*F5*U4+48.D0*Y3*F5*U5+
     *              48.D0*X3*F6*V5+18.D0*Y3*F2*U1-16.D0*X3*F5*V4+
     *              ANS1) * (-XSU2520)
C
      ANS1 = 32.D0*Y2*F3*U5-16.D0*X2*F6*V4-16.D0*Y2*F4*U6-16.D0*X3*F6*V4
     *+16.D0*Y2*F5*U3-64.D0*Y3*F4*U5+16.D0*X2*F4*V6+16.D0*Y2*F3*U4
     *+64.D0*X3*F4*V5-48.D0*X2*F4*V4-16.D0*Y2*F4*U1+16.D0*Y2*F5*U1
     *-80.D0*Y3*F4*U4-96.D0*Y3*F4*U2-16.D0*X2*F5*V1-16.D0*X2*F3*V4
     *-78.D0*X3*F2*V2+80.D0*X3*F4*V4+16.D0*Y3*F6*U4-9.D0*X2*F3*V3
     *-24.D0*Y3*F6*U2-48.D0*X2*F4*V5+48.D0*Y3*F6*U6-48.D0*Y2*F5*U4
     *+48.D0*Y3*F2*U4+16.D0*Y2*F5*U6+9.D0*X3*F2*V3-4.D0*Y3*F6*U1
     *-16.D0*X2*F5*V6-32.D0*Y3*F4*U6+32.D0*Y3*F6*U5+16.D0*X2*F4*V3
     *+32.D0*X3*F5*V5+9.D0*X3*F2*V1+12.D0*Y3*F2*U6-20.D0*Y2*F6*U3
     *+48.D0*Y2*F4*U4+20.D0*Y3*F4*U3+18.D0*Y2*F3*U2+48.D0*Y3*F2*U5
     *+24.D0*X3*F6*V2+48.D0*X2*F5*V4-18.D0*X2*F3*V2-48.D0*X3*F2*V4
     *+96.D0*X3*F4*V2+9.D0*Y2*F3*U3+20.D0*Y2*F6*U1+78.D0*Y3*F2*U2
     *+4.D0*X3*F6*V1-48.D0*Y3*F5*U6+4.D0*Y3*F5*U1+48.D0*Y2*F4*U5
     *-48.D0*Y2*F5*U5+48.D0*X2*F5*V5+48.D0*X3*F5*V6-48.D0*X3*F2*V5
     *-24.D0*X3*F5*V2+32.D0*X3*F4*V6-120.D0*X2*F4*V2+11.D0*X2*F3*V1
     *-4.D0*X3*F5*V1-20.D0*X2*F6*V1-48.D0*X3*F6*V6-9.D0*Y3*F2*U3
     *-16.D0*X2*F5*V3-32.D0*X2*F3*V5-16.D0*Y3*F5*U3+16.D0*Y3*F6*U3
     *-20.D0*X3*F4*V3+16.D0*X3*F5*V3-120.D0*Y2*F5*U2+120.D0*Y2*F4*U2
     *+120.D0*X2*F5*V2-16.D0*Y2*F6*U5+16.D0*X2*F6*V5
      W2(IELEM) = (-11.D0*Y2*F3*U1-16.D0*Y2*F4*U3-20.D0*X2*F3*V6-
     *              16.D0*X3*F6*V3+20.D0*Y2*F3*U6-12.D0*X3*F2*V6+
     *              24.D0*Y3*F5*U2+16.D0*X2*F4*V1+16.D0*Y2*F6*U4-
     *              16.D0*Y3*F5*U4-32.D0*Y3*F5*U5-32.D0*X3*F6*V5+
     *              20.D0*X2*F6*V3-9.D0*Y3*F2*U1+16.D0*X3*F5*V4+
     *              ANS1) * XSU2520
C
      ANS1 = -48.D0*Y2*F3*U5-32.D0*X2*F6*V4-16.D0*Y2*F4*U6-
     *        16.D0*X3*F6*V4-24.D0*Y2*F5*U3+16.D0*Y3*F4*U5+
     *        16.D0*X2*F4*V6-12.D0*Y2*F3*U4-16.D0*X3*F4*V5+
     *        48.D0*X2*F4*V4+4.D0*Y2*F4*U1-4.D0*Y2*F5*U1+20.D0*Y3*F4*U2+
     *        4.D0*X2*F5*V1+12.D0*X2*F3*V4+9.D0*X3*F2*V2+16.D0*Y3*F6*U4+
     *        78.D0*X2*F3*V3+16.D0*Y3*F6*U2+32.D0*X2*F4*V5-
     *        48.D0*Y3*F6*U6+48.D0*Y2*F5*U4-20.D0*Y3*F2*U4+
     *        16.D0*Y2*F5*U6+18.D0*X3*F2*V3+16.D0*Y3*F6*U1-
     *        16.D0*X2*F5*V6-16.D0*Y3*F4*U6-48.D0*Y3*F6*U5-
     *        24.D0*X2*F4*V3-48.D0*X3*F5*V5-11.D0*X3*F2*V1-
     *        16.D0*Y3*F2*U6+96.D0*Y2*F6*U3-48.D0*Y2*F4*U4+
     *        9.D0*Y2*F3*U2-32.D0*Y3*F2*U5-16.D0*X3*F6*V2-
     *        48.D0*X2*F5*V4-9.D0*X2*F3*V2+20.D0*X3*F2*V4-
     *        20.D0*X3*F4*V2-78.D0*Y2*F3*U3-9.D0*Y3*F2*U2-
     *        16.D0*X3*F6*V1+48.D0*Y3*F5*U6-16.D0*Y3*F5*U1-
     *        32.D0*Y2*F4*U5+20.D0*X3*F4*V1+32.D0*Y2*F5*U5-
     *        32.D0*X2*F5*V5-48.D0*X3*F5*V6+32.D0*X3*F2*V5+
     *        16.D0*X3*F5*V2+16.D0*X3*F4*V6+16.D0*X2*F4*V2-
     *        9.D0*X2*F3*V1+16.D0*X3*F5*V1+48.D0*X3*F6*V6-
     *        18.D0*Y3*F2*U3+24.D0*X2*F5*V3+48.D0*X2*F3*V5
      W3(IELEM) =(120.D0*Y3*F5*U3-120.D0*Y3*F6*U3-120.D0*X3*F5*V3+
     *            16.D0*Y2*F5*U2-16.D0*Y2*F4*U2-16.D0*X2*F5*V2-
     *            20.D0*Y2*F6*U2+64.D0*Y2*F6*U5+20.D0*X2*F6*V2+
     *            80.D0*Y2*F6*U6-20.D0*Y3*F4*U1-80.D0*X2*F6*V6-
     *            64.D0*X2*F6*V5+9.D0*Y2*F3*U1+24.D0*Y2*F4*U3+
     *            48.D0*X2*F3*V6+120.D0*X3*F6*V3-48.D0*Y2*F3*U6+
     *            16.D0*X3*F2*V6-16.D0*Y3*F5*U2-4.D0*X2*F4*V1+
     *            32.D0*Y2*F6*U4-16.D0*Y3*F5*U4+48.D0*Y3*F5*U5+
     *            48.D0*X3*F6*V5-96.D0*X2*F6*V3+11.D0*Y3*F2*U1+
     *            16.D0*X3*F5*V4+ANS1) * XSU2520
C
      ANS1 = 4.D0*Y2*F3*U5-64.D0*X2*F6*V4-32.D0*Y2*F4*U6-
     *       32.D0*X3*F6*V4-12.D0*Y2*F5*U3+16.D0*Y3*F4*U5+
     *       32.D0*X2*F4*V6-24.D0*Y2*F3*U4-16.D0*X3*F4*V5+
     *       96.D0*X2*F4*V4+8.D0*Y2*F4*U1-8.D0*Y2*F5*U1+
     *       20.D0*Y3*F4*U2+8.D0*X2*F5*V1+24.D0*X2*F3*V4+
     *       12.D0*X3*F2*V2+32.D0*Y3*F6*U4-3.D0*X2*F3*V3-
     *       4.D0*Y3*F6*U2+48.D0*X2*F4*V5+32.D0*Y3*F6*U6+
     *       96.D0*Y2*F5*U4-40.D0*Y3*F2*U4+32.D0*Y2*F5*U6-
     *       5.D0*X3*F2*V3-4.D0*Y3*F6*U1-32.D0*X2*F5*V6-16.D0*Y3*F4*U6+
     *       32.D0*Y3*F6*U5-12.D0*X2*F4*V3+32.D0*X3*F5*V5-
     *       8.D0*X3*F2*V1-4.D0*Y3*F2*U6-8.D0*Y2*F6*U3-96.D0*Y2*F4*U4-
     *       4.D0*Y2*F3*U2-20.D0*Y3*F2*U5+4.D0*X3*F6*V2-96.D0*X2*F5*V4+
     *       4.D0*X2*F3*V2+40.D0*X3*F2*V4-20.D0*X3*F4*V2+3.D0*Y2*F3*U3+
     *       16.D0*Y2*F6*U1-12.D0*Y3*F2*U2+4.D0*X3*F6*V1-32.D0*Y3*F5*U6+
     *       4.D0*Y3*F5*U1-48.D0*Y2*F4*U5+20.D0*X3*F4*V1+48.D0*Y2*F5*U5-
     *       48.D0*X2*F5*V5+32.D0*X3*F5*V6+20.D0*X3*F2*V5-4.D0*X3*F5*V2+
     *       16.D0*X3*F4*V6+12.D0*X2*F4*V2+4.D0*X2*F3*V1-4.D0*X3*F5*V1-
     *       16.D0*X2*F6*V1-32.D0*X3*F6*V6+5.D0*Y3*F2*U3+12.D0*X2*F5*V3-
     *       4.D0*X2*F3*V5+4.D0*Y3*F5*U3-4.D0*Y3*F6*U3-4.D0*X3*F5*V3+
     *       12.D0*Y2*F5*U2-12.D0*Y2*F4*U2-12.D0*X2*F5*V2-4.D0*Y2*F6*U2
      W4(IELEM) = (4.D0*X2*F6*V2+16.D0*Y2*F6*U6-20.D0*Y3*F4*U1-
     *            16.D0*X2*F6*V6-4.D0*Y2*F3*U1+12.D0*Y2*F4*U3-
     *            4.D0*X2*F3*V6+4.D0*X3*F6*V3+4.D0*Y2*F3*U6+
     *            4.D0*X3*F2*V6+4.D0*Y3*F5*U2-8.D0*X2*F4*V1+
     *            64.D0*Y2*F6*U4-32.D0*Y3*F5*U4-32.D0*Y3*F5*U5-
     *            32.D0*X3*F6*V5+8.D0*X2*F6*V3+8.D0*Y3*F2*U1+
     *            32.D0*X3*F5*V4+ ANS1) * (-XSU630)
C
      ANS1 = -40.D0*Y2*F3*U5+32.D0*Y2*F4*U6+32.D0*X3*F6*V4+
     *        8.D0*Y2*F5*U3-64.D0*Y3*F4*U5-32.D0*X2*F4*V6-4.D0*Y2*F3*U4+
     *        64.D0*X3*F4*V5-48.D0*X2*F4*V4-12.D0*Y2*F4*U1+
     *        12.D0*Y2*F5*U1-16.D0*Y3*F4*U4-16.D0*Y3*F4*U2-
     *        12.D0*X2*F5*V1+4.D0*X2*F3*V4-12.D0*X3*F2*V2+
     *        16.D0*X3*F4*V4-32.D0*Y3*F6*U4+12.D0*X2*F3*V3+
     *        8.D0*Y3*F6*U2-96.D0*X2*F4*V5-48.D0*Y3*F6*U6-
     *        48.D0*Y2*F5*U4+20.D0*Y3*F2*U4-32.D0*Y2*F5*U6+
     *        8.D0*X3*F2*V3+12.D0*Y3*F6*U1+32.D0*X2*F5*V6-
     *        96.D0*Y3*F6*U5+8.D0*X2*F4*V3-96.D0*X3*F5*V5+
     *        5.D0*X3*F2*V1+4.D0*Y3*F2*U6+16.D0*Y2*F6*U3+48.D0*Y2*F4*U4+
     *        4.D0*Y3*F4*U3+8.D0*Y2*F3*U2+40.D0*Y3*F2*U5-
     *        8.D0*X3*F6*V2+48.D0*X2*F5*V4-8.D0*X2*F3*V2-
     *        20.D0*X3*F2*V4+16.D0*X3*F4*V2-12.D0*Y2*F3*U3-
     *        8.D0*Y2*F6*U1+12.D0*Y3*F2*U2-12.D0*X3*F6*V1+
     *        48.D0*Y3*F5*U6-12.D0*Y3*F5*U1+96.D0*Y2*F4*U5-
     *        8.D0*X3*F4*V1-96.D0*Y2*F5*U5+96.D0*X2*F5*V5-
     *        48.D0*X3*F5*V6-40.D0*X3*F2*V5+8.D0*X3*F5*V2-
     *        12.D0*X2*F4*V2-5.D0*X2*F3*V1+12.D0*X3*F5*V1+
     *        8.D0*X2*F6*V1+48.D0*X3*F6*V6-8.D0*Y3*F2*U3-8.D0*X2*F5*V3
      W5(IELEM) = (40.D0*X2*F3*V5+12.D0*Y3*F5*U3-12.D0*Y3*F6*U3-
     *            4.D0*X3*F4*V3-12.D0*X3*F5*V3-12.D0*Y2*F5*U2+
     *            12.D0*Y2*F4*U2+12.D0*X2*F5*V2-4.D0*Y2*F6*U2+
     *            64.D0*Y2*F6*U5+4.D0*X2*F6*V2+16.D0*Y2*F6*U6+
     *            8.D0*Y3*F4*U1-16.D0*X2*F6*V6-64.D0*X2*F6*V5+
     *            5.D0*Y2*F3*U1-8.D0*Y2*F4*U3+20.D0*X2*F3*V6+
     *            12.D0*X3*F6*V3-20.D0*Y2*F3*U6-4.D0*X3*F2*V6-
     *            8.D0*Y3*F5*U2+12.D0*X2*F4*V1+32.D0*Y3*F5*U4+
     *            96.D0*Y3*F5*U5+96.D0*X3*F6*V5-16.D0*X2*F6*V3-
     *            5.D0*Y3*F2*U1-32.D0*X3*F5*V4+ ANS1) * XSU630
C
      ANS1 = 20.D0*Y2*F3*U5-16.D0*X2*F6*V4-32.D0*Y2*F4*U6-
     *       32.D0*X3*F6*V4-4.D0*Y2*F5*U3+32.D0*X2*F4*V6+4.D0*Y2*F3*U4+
     *       32.D0*X2*F4*V4+4.D0*Y2*F4*U1-4.D0*Y2*F5*U1-16.D0*Y3*F4*U4+
     *       8.D0*Y3*F4*U2+4.D0*X2*F5*V1-4.D0*X2*F3*V4+3.D0*X3*F2*V2+
     *       16.D0*X3*F4*V4+32.D0*Y3*F6*U4-12.D0*X2*F3*V3-
     *       12.D0*Y3*F6*U2+32.D0*X2*F4*V5+96.D0*Y3*F6*U6+
     *       32.D0*Y2*F5*U4-4.D0*Y3*F2*U4+32.D0*Y2*F5*U6-
     *       4.D0*X3*F2*V3-8.D0*Y3*F6*U1-32.D0*X2*F5*V6-64.D0*Y3*F4*U6+
     *       48.D0*Y3*F6*U5-4.D0*X2*F4*V3+48.D0*X3*F5*V5-4.D0*X3*F2*V1+
     *       24.D0*Y3*F2*U6-20.D0*Y2*F6*U3-32.D0*Y2*F4*U4+4.D0*Y3*F4*U3-
     *       5.D0*Y2*F3*U2-4.D0*Y3*F2*U5+12.D0*X3*F6*V2-32.D0*X2*F5*V4+
     *       5.D0*X2*F3*V2+4.D0*X3*F2*V4-8.D0*X3*F4*V2+12.D0*Y2*F3*U3+
     *       20.D0*Y2*F6*U1-3.D0*Y3*F2*U2+8.D0*X3*F6*V1-96.D0*Y3*F5*U6+
     *       8.D0*Y3*F5*U1-32.D0*Y2*F4*U5+16.D0*X3*F4*V1+32.D0*Y2*F5*U5-
     *       32.D0*X2*F5*V5+96.D0*X3*F5*V6+4.D0*X3*F2*V5-12.D0*X3*F5*V2+
     *       64.D0*X3*F4*V6-4.D0*X2*F4*V2+8.D0*X2*F3*V1-8.D0*X3*F5*V1-
     *       20.D0*X2*F6*V1-96.D0*X3*F6*V6+4.D0*Y3*F2*U3+4.D0*X2*F5*V3-
     *       20.D0*X2*F3*V5-12.D0*Y3*F5*U3+12.D0*Y3*F6*U3-4.D0*X3*F4*V3+
     *       12.D0*X3*F5*V3-4.D0*Y2*F5*U2+4.D0*Y2*F4*U2+4.D0*X2*F5*V2-
     *       16.D0*Y2*F6*U5-16.D0*Y3*F4*U1+16.D0*X2*F6*V5-8.D0*Y2*F3*U1
      W6(IELEM) = (4.D0*Y2*F4*U3-40.D0*X2*F3*V6-12.D0*X3*F6*V3+
     *            40.D0*Y2*F3*U6-24.D0*X3*F2*V6+12.D0*Y3*F5*U2-
     *            4.D0*X2*F4*V1+16.D0*Y2*F6*U4-32.D0*Y3*F5*U4-
     *            48.D0*Y3*F5*U5-48.D0*X3*F6*V5+20.D0*X2*F6*V3+
     *            4.D0*Y3*F2*U1+32.D0*X3*F5*V4+ANS1) * (-XSU630)
C
3     CONTINUE
C
      ENDIF
C
C     FONCTION F P2 ET VECTEUR U LINEAIRE
C
      ELSEIF(IELMF.EQ.13.AND.IELMU.EQ.11.AND.IELMV.EQ.11) THEN
C
      IF(FORMUL(14:16).EQ.'PSI') THEN
C
C     SCHEMA PSI P1 PUIS INTERPOLATION LINEAIRE
C
      DO 44 IELEM = 1 , NELEM
C
         X2 = XEL(IELEM+NELMAX)
         X3 = XEL(IELEM+2*NELMAX)
         Y2 = YEL(IELEM+NELMAX)
         Y3 = YEL(IELEM+2*NELMAX)
C
         F1 = F(IKLE1(IELEM))
         F2 = F(IKLE2(IELEM))
         F3 = F(IKLE3(IELEM))
C
         U1 = U(IKLE1(IELEM))
         U2 = U(IKLE2(IELEM))
         U3 = U(IKLE3(IELEM))
         V1 = V(IKLE1(IELEM))
         V2 = V(IKLE2(IELEM))
         V3 = V(IKLE3(IELEM))
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
         BETAN1 = L12*(F1-F2) + L13*(F1-F3)
         BETAN2 = L21*(F2-F1) + L23*(F2-F3)
         BETAN3 = L31*(F3-F1) + L32*(F3-F2)
C
         PHIT = BETAN1 + BETAN2 + BETAN3
C
         IF(PHIT.GT.0.D0) THEN
           W1(IELEM) =   XMUL * MAX( MIN( BETAN1, PHIT),0.D0 )
           W2(IELEM) =   XMUL * MAX( MIN( BETAN2, PHIT),0.D0 )
           W3(IELEM) =   XMUL * MAX( MIN( BETAN3, PHIT),0.D0 )
         ELSE
           W1(IELEM) = - XMUL * MAX( MIN(-BETAN1,-PHIT),0.D0 )
           W2(IELEM) = - XMUL * MAX( MIN(-BETAN2,-PHIT),0.D0 )
           W3(IELEM) = - XMUL * MAX( MIN(-BETAN3,-PHIT),0.D0 )
         ENDIF
         W4(IELEM) =   (W1(IELEM)+ W2(IELEM))/2.D0
         W5(IELEM) =   (W2(IELEM)+ W3(IELEM))/2.D0
         W6(IELEM) =   (W3(IELEM)+ W1(IELEM))/2.D0
C
44    CONTINUE
C
      ELSE
C
C     CALCUL CLASSIQUE
C
      DO IELEM = 1 , NELEM
C
         X2 = XEL(IELEM+NELMAX)
         X3 = XEL(IELEM+2*NELMAX)
         Y2 = YEL(IELEM+NELMAX)
         Y3 = YEL(IELEM+2*NELMAX)
C
         U1 = U(IKLE1(IELEM))
         U2 = U(IKLE2(IELEM))
         U3 = U(IKLE3(IELEM))
         V1 = V(IKLE1(IELEM))
         V2 = V(IKLE2(IELEM))
         V3 = V(IKLE3(IELEM))
C
         F1 = F(IKLE1(IELEM))
         F2 = F(IKLE2(IELEM)) - F1
         F3 = F(IKLE3(IELEM)) - F1
         F4 = F(IKLE4(IELEM)) - F1
         F5 = F(IKLE5(IELEM)) - F1
         F6 = F(IKLE6(IELEM)) - F1
C
      W1(IELEM) = (-4.D0*X3*F6*V2-4.D0*Y3*F5*U2-4.D0*Y2*F6*U2-
     *             8.D0*Y2*F4*U2-4.D0*X2*F5*V3-4.D0*Y2*F4*U3+
     *             8.D0*Y3*F6*U3-X2*F3*V2+4.D0*X3*F5*V2+4.D0*Y2*F5*U3+
     *             4.D0*Y3*F6*U2+8.D0*X3*F5*V3+4.D0*X2*F4*V3-
     *             6.D0*Y3*F2*U1+24.D0*X2*F6*V1+Y2*F3*U2+
     *             6.D0*Y2*F3*U1+8.D0*X2*F4*V2-8.D0*X2*F5*V2+
     *             8.D0*Y3*F4*U2-24.D0*X3*F4*V1+8.D0*Y2*F5*U2-
     *             5.D0*Y3*F2*U2-8.D0*X3*F4*V2-24.D0*Y2*F6*U1+
     *             6.D0*X3*F2*V1-6.D0*X2*F3*V1-Y3*F2*U3+5.D0*Y2*F3*U3-
     *             8.D0*Y2*F6*U3+8.D0*X2*F6*V3+4.D0*X2*F6*V2-
     *             4.D0*X3*F4*V3-8.D0*X3*F6*V3-5.D0*X2*F3*V3+
     *             24.D0*Y3*F4*U1+4.D0*Y3*F4*U3+5.D0*X3*F2*V2-
     *             8.D0*Y3*F5*U3+X3*F2*V3) * XSU360
C
      W2(IELEM)= (-24.D0*Y2*F4*U2-8.D0*Y3*F6*U3+6.D0*X2*F3*V2-
     *            4.D0*Y3*F6*U1-8.D0*X3*F5*V3+4.D0*X3*F6*V1-
     *            3.D0*Y3*F2*U1+4.D0*X2*F6*V1-6.D0*Y2*F3*U2-Y2*F3*U1+
     *            4.D0*Y3*F5*U1+24.D0*X2*F4*V2-24.D0*X2*F5*V2+
     *            24.D0*Y3*F4*U2-8.D0*X3*F4*V1+24.D0*Y2*F5*U2-
     *            18.D0*Y3*F2*U2-24.D0*X3*F4*V2-4.D0*Y2*F6*U1+
     *            3.D0*X3*F2*V1+X2*F3*V1-3.D0*Y3*F2*U3-5.D0*Y2*F3*U3-
     *            4.D0*X3*F5*V1+4.D0*Y2*F6*U3-4.D0*X2*F6*V3-
     *            4.D0*X3*F4*V3+8.D0*X3*F6*V3+5.D0*X2*F3*V3+
     *            8.D0*Y3*F4*U1+4.D0*Y3*F4*U3+18.D0*X3*F2*V2+
     *            8.D0*Y3*F5*U3+3.D0*X3*F2*V3) * (-XSU360)
C
      W3(IELEM)= (4.D0*X2*F5*V1-4.D0*Y2*F6*U2-4.D0*Y2*F5*U1+
     *            8.D0*Y2*F4*U2+24.D0*Y3*F6*U3-3.D0*X2*F3*V2+
     *            24.D0*X3*F5*V3+Y3*F2*U1+8.D0*X2*F6*V1+4.D0*Y2*F4*U1+
     *            3.D0*Y2*F3*U2+3.D0*Y2*F3*U1-4.D0*X2*F4*V1-
     *            8.D0*X2*F4*V2+8.D0*X2*F5*V2-4.D0*Y3*F4*U2-
     *            4.D0*X3*F4*V1-8.D0*Y2*F5*U2+5.D0*Y3*F2*U2+
     *            4.D0*X3*F4*V2-8.D0*Y2*F6*U1-X3*F2*V1-3.D0*X2*F3*V1+
     *            6.D0*Y3*F2*U3+18.D0*Y2*F3*U3-24.D0*Y2*F6*U3+
     *            24.D0*X2*F6*V3+4.D0*X2*F6*V2-24.D0*X3*F6*V3-
     *            18.D0*X2*F3*V3+4.D0*Y3*F4*U1-5.D0*X3*F2*V2-
     *            24.D0*Y3*F5*U3-6.D0*X3*F2*V3)*(-XSU360)
C      
      W4(IELEM) = (4.D0*X3*F6*V2+8.D0*X2*F5*V1+4.D0*Y3*F5*U2-
     *             4.D0*Y2*F6*U2-8.D0*Y2*F5*U1+12.D0*Y2*F4*U2+
     *             4.D0*X2*F5*V3+4.D0*Y2*F4*U3-4.D0*Y3*F6*U3-
     *             2.D0*X2*F3*V2-4.D0*X3*F5*V2-4.D0*Y3*F6*U1-
     *             4.D0*Y2*F5*U3-4.D0*Y3*F6*U2-4.D0*X3*F5*V3+
     *             4.D0*X3*F6*V1-4.D0*X2*F4*V3+2.D0*Y3*F2*U1+
     *             8.D0*X2*F6*V1+8.D0*Y2*F4*U1+2.D0*Y2*F3*U2+
     *             2.D0*Y2*F3*U1+4.D0*Y3*F5*U1-8.D0*X2*F4*V1-
     *             12.D0*X2*F4*V2+12.D0*X2*F5*V2-4.D0*Y3*F4*U2-
     *             4.D0*X3*F4*V1-12.D0*Y2*F5*U2+6.D0*Y3*F2*U2+
     *             4.D0*X3*F4*V2-8.D0*Y2*F6*U1-2.D0*X3*F2*V1-
     *             2.D0*X2*F3*V1+Y3*F2*U3-Y2*F3*U3-4.D0*X3*F5*V1+
     *             4.D0*X2*F6*V2+4.D0*X3*F6*V3+X2*F3*V3+4.D0*Y3*F4*U1-
     *             6.D0*X3*F2*V2+4.D0*Y3*F5*U3-X3*F2*V3)*XSU90
C
      W5(IELEM) = (8.D0*X3*F6*V2+4.D0*X2*F5*V1+8.D0*Y3*F5*U2+
     *             4.D0*Y2*F6*U2-4.D0*Y2*F5*U1+12.D0*Y2*F4*U2+
     *             8.D0*X2*F5*V3+8.D0*Y2*F4*U3-12.D0*Y3*F6*U3+
     *             2.D0*X2*F3*V2-8.D0*X3*F5*V2-4.D0*Y3*F6*U1-
     *             8.D0*Y2*F5*U3-8.D0*Y3*F6*U2-12.D0*X3*F5*V3+
     *             4.D0*X3*F6*V1-8.D0*X2*F4*V3+Y3*F2*U1+4.D0*Y2*F4*U1-
     *             2.D0*Y2*F3*U2-Y2*F3*U1+4.D0*Y3*F5*U1-4.D0*X2*F4*V1-
     *             12.D0*X2*F4*V2+12.D0*X2*F5*V2-8.D0*Y3*F4*U2-
     *             12.D0*Y2*F5*U2+6.D0*Y3*F2*U2+8.D0*X3*F4*V2-X3*F2*V1+
     *             X2*F3*V1+2.D0*Y3*F2*U3-6.D0*Y2*F3*U3-4.D0*X3*F5*V1+
     *             8.D0*Y2*F6*U3-8.D0*X2*F6*V3-4.D0*X2*F6*V2+
     *             4.D0*X3*F4*V3+12.D0*X3*F6*V3+6.D0*X2*F3*V3-
     *             4.D0*Y3*F4*U3-6.D0*X3*F2*V2+12.D0*Y3*F5*U3-
     *             2.D0*X3*F2*V3) * XSU90
      W6(IELEM) = (4.D0*X3*F6*V2+4.D0*X2*F5*V1+4.D0*Y3*F5*U2-
     *             4.D0*Y2*F5*U1+4.D0*Y2*F4*U2+4.D0*X2*F5*V3+
     *             4.D0*Y2*F4*U3-12.D0*Y3*F6*U3+X2*F3*V2-
     *             4.D0*X3*F5*V2-8.D0*Y3*F6*U1-4.D0*Y2*F5*U3-
     *             4.D0*Y3*F6*U2-12.D0*X3*F5*V3+8.D0*X3*F6*V1-
     *             4.D0*X2*F4*V3-2.D0*Y3*F2*U1+4.D0*X2*F6*V1+
     *             4.D0*Y2*F4*U1-Y2*F3*U2-2.D0*Y2*F3*U1+8.D0*Y3*F5*U1-
     *             4.D0*X2*F4*V1-4.D0*X2*F4*V2+4.D0*X2*F5*V2-
     *             8.D0*X3*F4*V1-4.D0*Y2*F5*U2+Y3*F2*U2-
     *             4.D0*Y2*F6*U1+2.D0*X3*F2*V1+2.D0*X2*F3*V1-
     *             2.D0*Y3*F2*U3-6.D0*Y2*F3*U3-8.D0*X3*F5*V1+
     *             4.D0*Y2*F6*U3-4.D0*X2*F6*V3-4.D0*X3*F4*V3+
     *             12.D0*X3*F6*V3+6.D0*X2*F3*V3+8.D0*Y3*F4*U1+
     *             4.D0*Y3*F4*U3-X3*F2*V2+12.D0*Y3*F5*U3+2.D0*X3*F2*V3)
     *             * XSU90
C
      ENDDO
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
       IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
       IF (LNG.EQ.1) WRITE(LU,200) IELMU,SU%NAME
       IF (LNG.EQ.1) WRITE(LU,300)
       IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
       IF (LNG.EQ.1) WRITE(LU,201) IELMU,SU%NAME
       IF (LNG.EQ.1) WRITE(LU,301)
100    FORMAT(1X,'VC08CC (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,
     *        1X,'NOM REEL : ',A6)
200    FORMAT(1X,'DISCRETISATION DE U : ',1I6,
     *        1X,'NOM REEL : ',A6)
300    FORMAT(1X,'CAS NON PREVU')
101    FORMAT(1X,'VC08CC (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F:',1I6,
     *        1X,'REAL NAME: ',A6)
201    FORMAT(1X,'DISCRETIZATION OF U:',1I6,
     *        1X,'REAL NAME: ',A6)
301    FORMAT(1X,'CASE NOT IMPLEMENTED')
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
