C                       *****************
                        SUBROUTINE VC04PP
C                       *****************
C
     *( XMUL,SU,SV,SW,U,V,W,SF,SG,SH,F,G,H,X,Y,Z,
     *  IKLE1,IKLE2,IKLE3,IKLE4,IKLE5,IKLE6,NELEM,NELMAX,
     *  W1,W2,W3,W4,W5,W6,SPECAD,FORMUL)
C
C***********************************************************************
C BIEF VERSION 6.0        21/06/06    J-M HERVOUET (LNHE) 01 30 87 80 18
C                                       
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C 
C
C     SI FORMUL SE TERMINE PAR 'HOR'
C     DANS LE MAILLAGE TRANSFORME
C
C                 /            D(PSII*)           D(PSII*)
C     V  = XMUL  /    H * U * -------- +  H * V * --------   D(OMEGA*)
C      I        /OMEGA*          DX                 DY
C
C     ATTENTION : ICI MAILLAGE TRANSFORME
C
C
C     SI FORMUL SE TERMINE PAR 'TOT'
C     DANS LE MAILLAGE REEL
C
C                 /     ->  --->           
C     V  = XMUL  /      U * GRAD(PSI)   D(OMEGA*)
C      I        /OMEGA           
C
C     ATTENTION : ICI MAILLAGE REEL
C
C     LE TYPE DE DISCRETISATION DE PSII EST LE PRISME P1
C
C     ATTENTION : LE RESULTAT EST DANS W SOUS FORME NON ASSEMBLEE.
C
C     ATTENTION : IF SPECAD=.TRUE., THE ADVECTING FIELD IS NOT ONLY
C                 U AND V BUT U+DM1*GRAD(ZCONV)
C
C
C
C
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
C |      W1,2,3    |<-- |  VECTEUR RESULTAT SOUS FORME NON ASSEMBLEE.
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
      USE BIEF, EX_VC04PP => VC04PP
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX),IKLE3(NELMAX)
      INTEGER, INTENT(IN) :: IKLE4(NELMAX),IKLE5(NELMAX),IKLE6(NELMAX)
C
      DOUBLE PRECISION, INTENT(IN)   ::X(*),Y(*),Z(*)
      DOUBLE PRECISION, INTENT(INOUT)::W1(NELMAX),W2(NELMAX),W3(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT)::W4(NELMAX),W5(NELMAX),W6(NELMAX)
      DOUBLE PRECISION, INTENT(IN)   ::XMUL
C
      LOGICAL, INTENT(IN) :: SPECAD
      CHARACTER(LEN=16), INTENT(IN) :: FORMUL
C
C     STRUCTURES DE U,V ET DONNEES REELLES
C
      TYPE(BIEF_OBJ),   INTENT(IN) :: SU,SV,SW,SF,SG,SH
      DOUBLE PRECISION, INTENT(IN) :: U(*),V(*),W(*),F(*),G(*),H(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION SUR144,X1,X2,X3,Y1,Y2,Y3,H1,H2,H3,SHT
      DOUBLE PRECISION HU1,HU2,HUINF,HUSUP,HV1,HV2,HVINF,HVSUP
      DOUBLE PRECISION U1,U2,U3,U4,U5,U6,V1,V2,V3,V4,V5,V6
      DOUBLE PRECISION Q1,Q2,Q3,Q4,Q5,Q6,Z2,Z3,Z4,Z5,Z6,G2,G3
      DOUBLE PRECISION GRADGX,GRADGY,DET
      INTEGER I1,I2,I3,I4,I5,I6,IELEM2,NELEM2
      INTEGER IELEM,IELMU,IELMV,IELMW
C
C**********************************************************************
C
      SUR144 = XMUL / 144.D0
C
C-----------------------------------------------------------------------
C
      IELMU=SU%ELM
      IELMV=SV%ELM
      IELMW=SW%ELM
      NELEM2 = NBPTS(10)
C
C-----------------------------------------------------------------------
C
C   TERMES HORIZONTAUX
C
      IF(FORMUL(14:16).EQ.'HOR') THEN
C
C   BOUCLE SUR LES ELEMENTS
C
      IF(IELMU.EQ.41.AND.IELMV.EQ.41) THEN
C
C-----------------------------------------------------------------------
C
C  U ET V DISCRETISEES EN PRISME P1 :
C
      IF(.NOT.SPECAD) THEN
C
C     STANDARD CASE
C
      DO IELEM = 1 , NELEM
C
         I1 = IKLE1(IELEM)
         I2 = IKLE2(IELEM)
         I3 = IKLE3(IELEM)
         I4 = IKLE4(IELEM)
         I5 = IKLE5(IELEM)
         I6 = IKLE6(IELEM)
C
         H1 = Z(I4) - Z(I1)
         H2 = Z(I5) - Z(I2)
         H3 = Z(I6) - Z(I3)
         SHT = H1 + H2 + H3
C
         HU1 = (U(I1)+U(I2)+U(I3))*SHT + H1*U(I1) + H2*U(I2) + H3*U(I3)
         HU2 = (U(I4)+U(I5)+U(I6))*SHT + H1*U(I4) + H2*U(I5) + H3*U(I6)
         HUINF = (HU1+HU1+HU2) * SUR144
         HUSUP = (HU1+HU2+HU2) * SUR144
C
         HV1 = (V(I1)+V(I2)+V(I3))*SHT + H1*V(I1) + H2*V(I2) + H3*V(I3)
         HV2 = (V(I4)+V(I5)+V(I6))*SHT + H1*V(I4) + H2*V(I5) + H3*V(I6)
         HVINF = (HV1+HV1+HV2) * SUR144
         HVSUP = (HV1+HV2+HV2) * SUR144
C
         Y1 = Y(I2) - Y(I3)
         Y2 = Y(I3) - Y(I1)
         Y3 = Y(I1) - Y(I2)
C
         X1 = X(I3) - X(I2)
         X2 = X(I1) - X(I3)
         X3 = X(I2) - X(I1)
C
         W1(IELEM) = Y1*HUINF + X1*HVINF
         W2(IELEM) = Y2*HUINF + X2*HVINF
         W3(IELEM) = Y3*HUINF + X3*HVINF
         W4(IELEM) = Y1*HUSUP + X1*HVSUP
         W5(IELEM) = Y2*HUSUP + X2*HVSUP
         W6(IELEM) = Y3*HUSUP + X3*HVSUP
C
      ENDDO
C
      ELSE
C
C     CASE WITH SPECIFIC ADVECTING FIELD
C
      DO IELEM = 1 , NELEM
C
         I1 = IKLE1(IELEM)
         I2 = IKLE2(IELEM)
         I3 = IKLE3(IELEM)
         I4 = IKLE4(IELEM)
         I5 = IKLE5(IELEM)
         I6 = IKLE6(IELEM)
C
         H1 = Z(I4) - Z(I1)
         H2 = Z(I5) - Z(I2)
         H3 = Z(I6) - Z(I3)
         SHT = H1 + H2 + H3
C
         Y2 = Y(I2) - Y(I1)
         Y3 = Y(I3) - Y(I1)
         X2 = X(I2) - X(I1)
         X3 = X(I3) - X(I1)
         DET= X2*Y3-X3*Y2
         IELEM2 = MOD(IELEM-1,NELEM2) + 1
C        G IS PIECE-WISE LINEAR
C        IT IS ZCONV IN TELEMAC-3D
         G2 = G(IELEM2+NELEM2)-G(IELEM2)
         G3 = G(IELEM2+2*NELEM2)-G(IELEM2)
         GRADGX=(G2*Y3-G3*Y2)/DET
         GRADGY=(X2*G3-X3*G2)/DET
C
         U1=U(I1)+F(I1)*GRADGX
         U2=U(I2)+F(I2)*GRADGX
         U3=U(I3)+F(I3)*GRADGX
         U4=U(I4)+F(I4)*GRADGX
         U5=U(I5)+F(I5)*GRADGX
         U6=U(I6)+F(I6)*GRADGX
         V1=V(I1)+F(I1)*GRADGY
         V2=V(I2)+F(I2)*GRADGY
         V3=V(I3)+F(I3)*GRADGY
         V4=V(I4)+F(I4)*GRADGY
         V5=V(I5)+F(I5)*GRADGY
         V6=V(I6)+F(I6)*GRADGY
C
         HU1 = (U1+U2+U3)*SHT + H1*U1 + H2*U2 + H3*U3
         HU2 = (U4+U5+U6)*SHT + H1*U4 + H2*U5 + H3*U6
         HUINF = (HU1+HU1+HU2) * SUR144
         HUSUP = (HU1+HU2+HU2) * SUR144
C
         HV1 = (V1+V2+V3)*SHT + H1*V1 + H2*V2 + H3*V3
         HV2 = (V4+V5+V6)*SHT + H1*V4 + H2*V5 + H3*V6
         HVINF = (HV1+HV1+HV2) * SUR144
         HVSUP = (HV1+HV2+HV2) * SUR144
C
         Y1 = Y(I2) - Y(I3)
         Y2 = Y(I3) - Y(I1)
         Y3 = Y(I1) - Y(I2)
C
         X1 = X(I3) - X(I2)
         X2 = X(I1) - X(I3)
         X3 = X(I2) - X(I1)
C
         W1(IELEM) = Y1*HUINF + X1*HVINF
         W2(IELEM) = Y2*HUINF + X2*HVINF
         W3(IELEM) = Y3*HUINF + X3*HVINF
         W4(IELEM) = Y1*HUSUP + X1*HVSUP
         W5(IELEM) = Y2*HUSUP + X2*HVSUP
         W6(IELEM) = Y3*HUSUP + X3*HVSUP
C
      ENDDO
C
      ENDIF
C
C     ELSEIF(IELMU.EQ.  ) THEN
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
       IF (LNG.EQ.1) WRITE(LU,101) IELMU,SU%NAME
       IF (LNG.EQ.2) WRITE(LU,102) IELMU,SU%NAME
101    FORMAT(1X,'VC04PP (BIEF) :',/,
     *        1X,'DISCRETISATION DE U ET V : ',1I6,' CAS NON PREVU',/,
     *        1X,'NOM REEL DE U : ',A6)
102    FORMAT(1X,'VC04PP (BIEF) :',/,
     *        1X,'DISCRETISATION OF U ET V : ',1I6,' NOT IMPLEMENTED',/,
     *        1X,'REAL NAME OF U : ',A6)
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(FORMUL(14:16).EQ.'TOT') THEN
C
C     TERMES HORIZONTAUX
C
C     BOUCLE SUR LES ELEMENTS
C
C     HERE SIGMAG NOT TAKEN INTO ACCOUNT (IT WOULD BE TO IMPLEMENT IF A
C     COMPATIBILITY WITH 2D CONTINUITY EQUATION REQUIRED)
C
      IF(IELMU.EQ.41.AND.
     *   IELMV.EQ.41.AND.
     *   IELMW.EQ.41) THEN
C
C-----------------------------------------------------------------------
C
C     VITESSE DISCRETISEE EN PRISME P1 :
C
      DO IELEM = 1 , NELEM
C
         I1 = IKLE1(IELEM)
         I2 = IKLE2(IELEM)
         I3 = IKLE3(IELEM)
         I4 = IKLE4(IELEM)
         I5 = IKLE5(IELEM)
         I6 = IKLE6(IELEM)
C
         X2=X(I2)-X(I1)
         X3=X(I3)-X(I1)
         Y2=Y(I2)-Y(I1)
         Y3=Y(I3)-Y(I1)
C
         Z2=Z(I2)-Z(I1)
         Z3=Z(I3)-Z(I1)
         Z4=Z(I4)-Z(I1)
         Z5=Z(I5)-Z(I1)
         Z6=Z(I6)-Z(I1)
C
         U1=U(I1)
         U2=U(I2)
         U3=U(I3)
         U4=U(I4)
         U5=U(I5)
         U6=U(I6)
C
         V1=V(I1)
         V2=V(I2)
         V3=V(I3)
         V4=V(I4)
         V5=V(I5)
         V6=V(I6)
C
         Q1=W(I1)
         Q2=W(I2)
         Q3=W(I3)
         Q4=W(I4)
         Q5=W(I5)
         Q6=W(I6)
!
!        DERIVATIVE IN Z
!
         W1(IELEM) = (-3*Q6-3*Q5-3*Q2-3*Q3-6*Q4-6*Q1)*(Y3*X2-X3*Y2)
         W2(IELEM) = (-3*Q6-6*Q5-6*Q2-3*Q3-3*Q4-3*Q1)*(Y3*X2-X3*Y2)
         W3(IELEM) = (-6*Q6-3*Q1-3*Q4-6*Q3-3*Q2-3*Q5)*(Y3*X2-X3*Y2)
         W4(IELEM) = - W1(IELEM)
         W5(IELEM) = - W2(IELEM)
         W6(IELEM) = - W3(IELEM)
!
!        DERIVATIVE IN X
!
         W1(IELEM) = W1(IELEM) +
     #         (-2*Z2*U1+3*Z6*U3+Z6*U2-4*Z2*U2-Z6*U5+Z5*U4+3*Z4*U5+4*Z5*
     #U2-2*Z2*U5-2*Z3*U5+6*Z4*U4-3*Z6*U4-3*Z3*U6-6*Z3*U1-3*Z3*U4+2*Z5*U3
     #+2*Z5*U5+6*Z4*U1+3*Z4*U2+Z5*U6-2*Z2*U3-4*Z3*U2-Z2*U4-Z2*U6+3*Z4*U3
     #+3*Z4*U6+2*Z5*U1-6*Z3*U3)*Y2+(6*Z2*U1-4*Z6*U3-2*Z6*U6-2*Z6*U2+6*Z2
     #*U2-Z6*U5+3*Z5*U4-3*Z4*U5-3*Z5*U2+3*Z2*U5+Z3*U5-6*Z4*U4-Z6*U4+2*Z3
     #*U6+2*Z3*U1+Z3*U4-Z5*U3-6*Z4*U1-3*Z4*U2+Z5*U6+4*Z2*U3+2*Z3*U2+3*Z2
     #*U4+2*Z2*U6-3*Z4*U3-3*Z4*U6+4*Z3*U3-2*Z6*U1)*Y3
         W2(IELEM) = W2(IELEM) +
     #         (Z4*U3-4*Z3*U2+2*Z4*U2-2*Z6*U4+2*Z4*U6-2*Z3*U3-Z6*U3-2*Z6
     #*U6-2*Z3*U5+2*Z4*U4-2*Z3*U1+Z4*U1-Z3*U6-Z6*U1-2*Z6*U2-4*Z6*U5+4*Z4
     #*U5-Z3*U4)*Y2+(Z4*U3-Z4*U6-2*Z3*U2+6*Z5*U5+4*Z6*U3+2*Z6*U2-2*Z3*U6
     #+3*Z5*U3+3*Z5*U4-Z3*U4+6*Z5*U2+3*Z4*U1+Z6*U4-Z3*U5+2*Z6*U6-4*Z3*U3
     #+2*Z6*U1+3*Z5*U1+3*Z5*U6+Z6*U5-2*Z3*U1-3*Z4*U5)*Y3
         W3(IELEM) = W3(IELEM) +
     #         (-6*Z6*U6-2*Z5*U1-3*Z6*U1-Z5*U6+2*Z2*U5+4*Z2*U2+Z2*U4-2*Z
     #5*U5+2*Z2*U1-2*Z5*U3-3*Z6*U2+3*Z4*U6+Z2*U6-6*Z6*U3-4*Z5*U2+2*Z2*U3
     #-Z5*U4+Z4*U5-Z4*U2-3*Z4*U1-3*Z6*U5-3*Z6*U4)*Y2+(-2*Z4*U3-Z4*U2-4*Z
     #4*U6+2*Z2*U2+2*Z5*U3+4*Z2*U3-2*Z4*U4+2*Z5*U4+2*Z2*U6+Z5*U1+2*Z2*U1
     #-Z4*U1+Z5*U2+Z2*U4+Z2*U5+2*Z5*U5+4*Z5*U6-2*Z4*U5)*Y3
         W4(IELEM) = W4(IELEM) +
     #         (Z3*U2-2*Z2*U6+4*Z5*U5-2*Z2*U2+6*Z6*U6+3*Z6*U1+2*Z6*U2-Z3
     #*U5+2*Z5*U6-2*Z2*U4+2*Z5*U2-3*Z3*U6-Z2*U3+3*Z6*U3-Z2*U1+2*Z5*U4+6*
     #Z6*U4+3*Z3*U1+Z5*U3+4*Z6*U5-4*Z2*U5+Z5*U1)*Y2+(Z3*U2-4*Z5*U6-2*Z6*
     #U3-2*Z5*U3-6*Z5*U5-6*Z5*U4-Z6*U1+4*Z3*U6-Z2*U3+2*Z3*U3-2*Z6*U4+2*Z
     #3*U4-4*Z6*U6+2*Z3*U5+3*Z2*U5-2*Z6*U5+Z2*U6-3*Z5*U2-Z6*U2-3*Z2*U1+Z
     #3*U1-3*Z5*U1)*Y3
         W5(IELEM) = W5(IELEM) +
     #         (2*Z6*U6+Z6*U1-Z4*U1-2*Z4*U6-2*Z4*U4+Z3*U4+2*Z3*U3+2*Z3*U
     #5+4*Z3*U2+Z3*U6-2*Z4*U2+Z6*U3+2*Z6*U4-4*Z4*U5+2*Z6*U2+2*Z3*U1+4*Z6
     #*U5-Z4*U3)*Y2+(2*Z4*U3-Z3*U1-Z3*U2-6*Z2*U2+2*Z6*U3+4*Z4*U6+Z6*U1+3
     #*Z4*U2+6*Z4*U4-4*Z3*U6-3*Z2*U3+3*Z4*U1+2*Z6*U4-2*Z3*U4+6*Z4*U5-2*Z
     #3*U3+Z6*U2-3*Z2*U4-6*Z2*U5+2*Z6*U5-2*Z3*U5+4*Z6*U6-3*Z2*U6-3*Z2*U1
     #)*Y3
         W6(IELEM) = W6(IELEM) +
     #         (3*Z3*U5+2*Z2*U4-2*Z5*U2+3*Z3*U2-Z5*U1-2*Z5*U6+Z2*U1-2*Z4
     #*U2+2*Z2*U2+4*Z2*U5-Z5*U3-6*Z4*U6+2*Z2*U6+3*Z3*U4+6*Z3*U3+6*Z3*U6+
     #Z2*U3-6*Z4*U4+3*Z3*U1-2*Z5*U4-4*Z4*U5-3*Z4*U1-4*Z5*U5-3*Z4*U3)*Y2+
     #(2*Z4*U3-Z5*U1+4*Z4*U6-2*Z2*U2-2*Z5*U3+Z4*U2+Z4*U1-2*Z5*U4-4*Z2*U3
     #-Z5*U2-4*Z5*U6-2*Z5*U5-2*Z2*U6-2*Z2*U1-Z2*U5-Z2*U4+2*Z4*U5+2*Z4*U4
     #)*Y3
!
!        DERIVATIVE IN Y
!
         W1(IELEM) = W1(IELEM) +
     #         (-4*Z5*V2-3*Z4*V5-Z5*V6+Z6*V5+2*Z3*V5+4*Z2*V2+3*Z3*V6+Z2*
     #V4+4*Z3*V2-2*Z5*V3-Z6*V2+Z2*V6+6*Z3*V1-6*Z4*V1-2*Z5*V5+6*Z3*V3+2*Z
     #2*V3-3*Z4*V3-2*Z5*V1+3*Z3*V4-6*Z4*V4-3*Z4*V6+3*Z6*V4-3*Z6*V3-Z5*V4
     #+2*Z2*V1+2*Z2*V5-3*Z4*V2)*X2+(3*Z5*V2+3*Z4*V5-Z5*V6+Z6*V5-Z3*V5-6*
     #Z2*V2-2*Z3*V6-3*Z2*V4-2*Z3*V2+Z5*V3+2*Z6*V2-2*Z2*V6-2*Z3*V1+6*Z4*V
     #1-4*Z3*V3-4*Z2*V3+3*Z4*V3+2*Z6*V6-Z3*V4+6*Z4*V4+3*Z4*V6+Z6*V4+4*Z6
     #*V3-3*Z5*V4-6*Z2*V1+2*Z6*V1-3*Z2*V5+3*Z4*V2)*X3
         W2(IELEM) = W2(IELEM) +
     #         (-2*Z4*V6+2*Z6*V6+Z6*V3+2*Z6*V4-Z4*V3-Z4*V1+2*Z3*V1+Z3*V4
     #-4*Z4*V5+4*Z6*V5+2*Z6*V2+2*Z3*V3+Z3*V6+4*Z3*V2-2*Z4*V4+2*Z3*V5-2*Z
     #4*V2+Z6*V1)*X2+(-3*Z5*V1-6*Z5*V2+2*Z3*V1-3*Z5*V4+2*Z3*V6+4*Z3*V3+Z
     #3*V4+Z3*V5-6*Z5*V5-3*Z4*V1+3*Z4*V5-2*Z6*V2-Z6*V5-Z4*V3-2*Z6*V1-4*Z
     #6*V3-3*Z5*V6-Z6*V4+Z4*V6+2*Z3*V2-3*Z5*V3-2*Z6*V6)*X3
         W3(IELEM) = W3(IELEM) +
     #         (-2*Z2*V5+Z5*V4+6*Z6*V3+6*Z6*V6+3*Z6*V5-2*Z2*V1+3*Z6*V4+2
     #*Z5*V5+Z5*V6+3*Z4*V1-2*Z2*V3-Z4*V5+3*Z6*V1-Z2*V6-3*Z4*V6+2*Z5*V3-4
     #*Z2*V2-Z2*V4+3*Z6*V2+4*Z5*V2+Z4*V2+2*Z5*V1)*X2+(-4*Z5*V6-2*Z2*V2-2
     #*Z2*V1-Z2*V4-Z5*V1+Z4*V1+Z4*V2+2*Z4*V4-2*Z5*V5+4*Z4*V6-Z2*V5-2*Z2*
     #V6+2*Z4*V5-2*Z5*V4-4*Z2*V3-Z5*V2+2*Z4*V3-2*Z5*V3)*X3
         W4(IELEM) = W4(IELEM) +
     #         (-2*Z5*V4-6*Z6*V6-2*Z5*V6+4*Z2*V5-4*Z5*V5-4*Z6*V5+2*Z2*V2
     #+Z2*V3-6*Z6*V4+2*Z2*V6+3*Z3*V6-Z5*V1-2*Z6*V2-Z5*V3-3*Z6*V3+Z2*V1+2
     #*Z2*V4-3*Z6*V1-Z3*V2-2*Z5*V2-3*Z3*V1+Z3*V5)*X2+(3*Z5*V1+2*Z6*V5+3*
     #Z2*V1-Z3*V2-Z3*V1+2*Z6*V4-4*Z3*V6-2*Z3*V3+Z6*V2-2*Z3*V4+2*Z6*V3+6*
     #Z5*V5+Z2*V3-Z2*V6+2*Z5*V3+Z6*V1+4*Z5*V6+6*Z5*V4-2*Z3*V5+3*Z5*V2+4*
     #Z6*V6-3*Z2*V5)*X3
         W5(IELEM) = W5(IELEM) +
     #         (Z4*V3-Z6*V1-2*Z6*V6+4*Z4*V5-2*Z3*V3-Z6*V3+2*Z4*V4-Z3*V4+
     #2*Z4*V6-2*Z6*V4-4*Z6*V5-2*Z3*V1-Z3*V6-4*Z3*V2-2*Z6*V2-2*Z3*V5+2*Z4
     #*V2+Z4*V1)*X2+(-2*Z6*V5+3*Z2*V1-4*Z4*V6+Z3*V1-2*Z4*V3+4*Z3*V6+2*Z3
     #*V3+2*Z3*V4-6*Z4*V5-3*Z4*V1-3*Z4*V2-6*Z4*V4+3*Z2*V4-2*Z6*V3+6*Z2*V
     #2+3*Z2*V3-4*Z6*V6+3*Z2*V6-Z6*V1-Z6*V2+6*Z2*V5-2*Z6*V4+2*Z3*V5+Z3*V
     #2)*X3
         W6(IELEM) = W6(IELEM) +
     #         (6*Z4*V6-2*Z2*V2+2*Z5*V4-4*Z2*V5-6*Z3*V3+2*Z5*V6-Z2*V1+6*
     #Z4*V4+3*Z4*V1+3*Z4*V3+Z5*V1-3*Z3*V1-3*Z3*V4-2*Z2*V6-6*Z3*V6+Z5*V3-
     #Z2*V3-2*Z2*V4+4*Z5*V5-3*Z3*V2+2*Z5*V2-3*Z3*V5+2*Z4*V2+4*Z4*V5)*X2+
     #(2*Z5*V4+2*Z2*V1+Z5*V1-Z4*V2-2*Z4*V4+2*Z5*V5-Z4*V1+4*Z5*V6+2*Z2*V2
     #-2*Z4*V5+2*Z5*V3+Z5*V2+4*Z2*V3+Z2*V5+Z2*V4+2*Z2*V6-4*Z4*V6-2*Z4*V3
     #)*X3
!
         W1(IELEM) = W1(IELEM)*SUR144   
         W2(IELEM) = W2(IELEM)*SUR144
         W3(IELEM) = W3(IELEM)*SUR144
         W4(IELEM) = W4(IELEM)*SUR144
         W5(IELEM) = W5(IELEM)*SUR144
         W6(IELEM) = W6(IELEM)*SUR144
!
      ENDDO
C
C-----------------------------------------------------------------------
C
C     ELSEIF(IELMU.EQ.  ) THEN
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
       IF (LNG.EQ.1) WRITE(LU,301) IELMU,SU%NAME
       IF (LNG.EQ.2) WRITE(LU,302) IELMU,SU%NAME
301    FORMAT(1X,'VC04PP (BIEF) :',/,
     *        1X,'DISCRETISATION DE U,V,W : ',1I6,' CAS NON PREVU',/,
     *        1X,'NOM REEL DE U : ',A6)
302    FORMAT(1X,'VC04PP (BIEF) :',/,
     *        1X,'DISCRETISATION OF U,V,W : ',1I6,' NOT IMPLEMENTED',/,
     *        1X,'REAL NAME OF U : ',A6)
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
       IF (LNG.EQ.1) WRITE(LU,201) FORMUL
       IF (LNG.EQ.2) WRITE(LU,202) FORMUL
201    FORMAT(1X,'VC04PP (BIEF) :',/,
     *        1X,'IL MANQUE HOR OU TOT EN FIN DE FORMULE : ',A16)
202    FORMAT(1X,'VC04PP (BIEF) :',/,
     *        1X,'HOR OR TOT LACKING AT THE END OF THE FORMULA : ',A16)
       CALL PLANTE(1)
       STOP
C
C-----------------------------------------------------------------------
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
