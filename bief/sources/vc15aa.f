C                       *****************
                        SUBROUTINE VC15AA
C                       *****************
C
     *( XMUL,SF,SU,SV,F,U,V,
     *  XEL,YEL,
     *  IKLE1,IKLE2,IKLE3,IKLE4,NELEM,NELMAX,
     *  W1,W2,W3)
C
C***********************************************************************
C BIEF VERSION 5.1           09/12/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C
C                    /                D(FU)    D(FV)
C      V  =  XMUL   /       PSII  * (   --  +    -- )   D(OMEGA)
C       I          /OMEGA               DX       DY
C
C
C    PSI(I) EST UNE BASE DE TYPE TRIANGLE P1
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
C |      W1,2,3    |<-- |  VECTEUR RESULTAT SOUS FORME NON ASSEMBLEE.
C |      FORMUL    | -->|  FORMULE DE CALCUL (ARGUMENT DE VECTOR).
C |                |    |  (NON UTILISE, SERVIRA A CHOISIR ENTRE
C |                |    |   DIFFERENTS SCHEMAS DE CALCUL).
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
      USE BIEF !, EX_VC15AA => VC15AA
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX)
      INTEGER, INTENT(IN) :: IKLE3(NELMAX),IKLE4(NELMAX)
C
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,*),YEL(NELMAX,*)
      DOUBLE PRECISION, INTENT(INOUT) ::W1(NELMAX),W2(NELMAX),W3(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
C     STRUCTURES DE F,G,H,U,V,W ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF,SU,SV
      DOUBLE PRECISION, INTENT(IN) :: F(*),U(*),V(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF,IELMU,IELMV
      DOUBLE PRECISION X2,Y2,X3,Y3,F1,F2,F3,U1,U2,U3,U4,V1,V2,V3,V4
      DOUBLE PRECISION XSUR24,XSU216
C
C-----------------------------------------------------------------------
C
      XSUR24 = XMUL/24.D0
      XSU216 = XMUL/216.D0
C
      IELMF=SF%ELM
      IELMU=SU%ELM
      IELMV=SV%ELM
C
C-----------------------------------------------------------------------
C
C     FONCTION F ET VECTEUR U LINEAIRES
C
      IF(IELMF.EQ.11.AND.IELMU.EQ.11.AND.IELMV.EQ.11) THEN
C
      DO 3 IELEM = 1 , NELEM
C
         X2 = XEL(IELEM,2)
         X3 = XEL(IELEM,3)
         Y2 = YEL(IELEM,2)
         Y3 = YEL(IELEM,3)
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
         W1(IELEM)=(((2*F3+F2+F1)*V3+(F3-F2-4*F1)*V1+(F3-F1)*V2)*X2-((
     *          F3+2*F2+F1)*V2-(F3-F2+4*F1)*V1+(F2-F1)*V3)*X3-((Y3+Y2)*
     *          F3-(Y3+Y2)*F2+4*(Y3-Y2)*F1)*U1+((Y3+Y2)*F1+(Y3-Y2)*F3+2
     *          *F2*Y3)*U2-((Y3+Y2)*F1-(Y3-Y2)*F2+2*F3*Y2)*U3)*XSUR24
         W2(IELEM)=(-(((F3+4*F2-F1)*V2-(F3+F2+2*F1)*V1+(F2-F1)*V3)*X3-
     *          2*((F3+F2)*V3+(F3-F1)*V2-(F2+F1)*V1)*X2+(2*(Y3-Y2)*F1+(
     *          Y3-2*Y2)*F2+F3*Y3)*U1-((Y3-2*Y2)*F3-(Y3-2*Y2)*F1+4*F2
     *          *Y3)*U2-((Y3-2*Y2)*F2-2*F3*Y2-F1*Y3)*U3))*XSUR24
         W3(IELEM)=(((4*F3+F2-F1)*V3-(F3+F2+2*F1)*V1+(F3-F1)*V2)*X2-2
     *          *((F3+F2)*V2-(F3+F1)*V1+(F2-F1)*V3)*X3-((2*Y3-Y2)*F3+2*
     *           (Y3-Y2)*F1-F2*Y2)*U1+((2*Y3-Y2)*F3+2*F2*Y3+F1*Y2)*U2+((
     *             2*Y3-Y2)*F2-(2*Y3-Y2)*F1-4*F3*Y2)*U3)*XSUR24
C
3     CONTINUE
C
C-----------------------------------------------------------------------
C
C     FONCTION F LINEAIRE ET VECTEUR U QUASI-BULLE
C
      ELSEIF(IELMF.EQ.11.AND.IELMU.EQ.12.AND.IELMU.EQ.12) THEN
C
C
      DO 44 IELEM = 1 , NELEM
C
         X2 = XEL(IELEM,2)
         X3 = XEL(IELEM,3)
         Y2 = YEL(IELEM,2)
         Y3 = YEL(IELEM,3)
C
         F1 = F(IKLE1(IELEM))
         F2 = F(IKLE2(IELEM))
         F3 = F(IKLE3(IELEM))
C
         U1 = U(IKLE1(IELEM))
         U2 = U(IKLE2(IELEM))
         U3 = U(IKLE3(IELEM))
         U4 = U(IKLE4(IELEM))
         V1 = V(IKLE1(IELEM))
         V2 = V(IKLE2(IELEM))
         V3 = V(IKLE3(IELEM))
         V4 = V(IKLE4(IELEM))
C
         W1(IELEM)=((14*X2*V3+12*X2*V4+5*X2*V2+5*X2*V1+4*X3*V3-
     *     12*X3*V4-5*X3*V2+13*X3*V1-4*U3*Y3-14*U3*Y2+12*U4*Y3-
     *     12*U4*Y2+5*U2*Y3-5*U2*Y2-13*U1*Y3-5*U1*Y2)*F3+(5*X2*
     *     V3+12*X2*V4-4*X2*V2-13*X2*V1-5*X3*V3-12*X3*V4-14*X3
     *     *V2-5*X3*V1+5*U3*Y3-5*U3*Y2+12*U4*Y3-12*U4*Y2+14*U2
     *     *Y3+4*U2*Y2+5*U1*Y3+13*U1*Y2)*F2+(5*X2*V3+12*X2*V4-
     *     13*X2*V2-40*X2*V1+13*X3*V3-12*X3*V4-5*X3*V2+40*X3*V1
     *     -13*U3*Y3-5*U3*Y2+12*U4*Y3-12*U4*Y2+5*U2*Y3+13*U2*
     *     Y2-40*U1*Y3+40*U1*Y2)*F1)*XSU216
         W2(IELEM)=(18*X2*V3*F3+18*X2*V3*F2+18*X2*V2*F3-18*X2*V2*
     *     F1-18*X2*V1*F2-18*X2*V1*F1-4*X3*V3*F3-13*X3*V3*F2+5*
     *     X3*V3*F1+12*X3*V4*F3+12*X3*V4*F2+12*X3*V4*F1-13*X3*V2
     *     *F3-40*X3*V2*F2+5*X3*V2*F1+5*X3*V1*F3+5*X3*V1*F2+14*
     *     X3*V1*F1+4*U3*Y3*F3+13*U3*Y3*F2-5*U3*Y3*F1-18*U3*Y2*
     *     F3-18*U3*Y2*F2-12*U4*Y3*F3-12*U4*Y3*F2-12*U4*Y3*F1+
     *     13*U2*Y3*F3+40*U2*Y3*F2-5*U2*Y3*F1-18*U2*Y2*F3+18*U2*
     *     Y2*F1-5*U1*Y3*F3-5*U1*Y3*F2-14*U1*Y3*F1+18*U1*Y2*F2+
     *     18*U1*Y2*F1)*XSU216
          W3(IELEM)=(40*X2*V3*F3+13*X2*V3*F2-5*X2*V3*F1-12*X2*V4*F3
     *     -12*X2*V4*F2-12*X2*V4*F1+13*X2*V2*F3+4*X2*V2*F2-5*X2
     *     *V2*F1-5*X2*V1*F3-5*X2*V1*F2-14*X2*V1*F1-18*X3*V3*F2+
     *     18*X3*V3*F1-18*X3*V2*F3-18*X3*V2*F2+18*X3*V1*F3+18*X3
     *     *V1*F1+18*U3*Y3*F2-18*U3*Y3*F1-40*U3*Y2*F3-13*U3*Y2*
     *     F2+5*U3*Y2*F1+12*U4*Y2*F3+12*U4*Y2*F2+12*U4*Y2*F1+18
     *     *U2*Y3*F3+18*U2*Y3*F2-13*U2*Y2*F3-4*U2*Y2*F2+5*U2*Y2*
     *     F1-18*U1*Y3*F3-18*U1*Y3*F1+5*U1*Y2*F3+5*U1*Y2*F2+14*
     *     U1*Y2*F1)*XSU216
C
44    CONTINUE
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
       IF (LNG.EQ.2) WRITE(LU,201) IELMU,SU%NAME
       IF (LNG.EQ.2) WRITE(LU,301)
100    FORMAT(1X,'VC15AA (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,
     *        1X,'NOM REEL : ',A6)
200    FORMAT(1X,'DISCRETISATION DE U : ',1I6,
     *        1X,'NOM REEL : ',A6)
300    FORMAT(1X,'CAS NON PREVU')
101    FORMAT(1X,'VC15AA (BIEF) :',/,
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
 
 
