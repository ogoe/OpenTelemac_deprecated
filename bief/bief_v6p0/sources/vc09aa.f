C                       *****************
                        SUBROUTINE VC09AA
C                       *****************
C
     *( XMUL,SF,SG,SU,SV,F,G,U,V,
     *  XEL,YEL,IKLE1,IKLE2,IKLE3,NELEM,NELMAX,W1,W2,W3 )
C
C***********************************************************************
C BIEF VERSION 5.1           09/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C
C                    /                    DF        DF
C      V  =  XMUL   /       PSII  * ( G U --  + G V -- )   D(OMEGA)
C       I          /OMEGA                 DX        DY
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
      USE BIEF, EX_VC09AA => VC09AA
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
      DOUBLE PRECISION, INTENT(IN)    :: XEL(NELMAX,*),YEL(NELMAX,*)
      DOUBLE PRECISION, INTENT(INOUT) ::W1(NELMAX),W2(NELMAX),W3(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL
C
C     STRUCTURES DE F,G,U,V ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF,SG,SU,SV
      DOUBLE PRECISION, INTENT(IN) :: F(*),G(*),U(*),V(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF,IELMU,IELMG,IELMV
C
      DOUBLE PRECISION X2,Y2,X3,Y3,F1,F2,F3,G1,G2,G3,U1,U2,U3,V1,V2,V3
      DOUBLE PRECISION XS120,T11,T12,T13,T22,T23,T33,FTX,FTY
      DOUBLE PRECISION U123,V123
C
C-----------------------------------------------------------------------
C
      XS120 = XMUL / 120.D0
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
      IELMG=SG%ELM
      IELMU=SU%ELM
      IELMV=SV%ELM
C
C-----------------------------------------------------------------------
C
C     FONCTION F ET G ET VECTEUR U LINEAIRES
C
      IF(      IELMF.EQ.11.AND.IELMG.EQ.11
     *    .AND.IELMU.EQ.11.AND.IELMV.EQ.11  ) THEN
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
         G1 = G(IKLE1(IELEM))
         G2 = G(IKLE2(IELEM))
         G3 = G(IKLE3(IELEM))
C
         U1 = U(IKLE1(IELEM))
         U2 = U(IKLE2(IELEM))
         U3 = U(IKLE3(IELEM))
C
         V1 = V(IKLE1(IELEM))
         V2 = V(IKLE2(IELEM))
         V3 = V(IKLE3(IELEM))
C
         FTX = F1 * (X3-X2) - F2 * X3 + F3*X2
         FTY = F1 * (Y2-Y3) + F2 * Y3 - F3*Y2
C
         U123 = U1 + U2 + U3
         V123 = V1 + V2 + V3
C
         T11 = FTX * ( V123 + V1 + V1 ) + FTY * ( U123 + U1 + U1 )
         T22 = FTX * ( V123 + V2 + V2 ) + FTY * ( U123 + U2 + U2 )
         T33 = FTX * ( V123 + V3 + V3 ) + FTY * ( U123 + U3 + U3 )
C
         T12 = FTX * ( V123 + V123 - V3 ) + FTY * ( U123 + U123 - U3 )
         T13 = FTX * ( V123 + V123 - V2 ) + FTY * ( U123 + U123 - U2 )
         T23 = FTX * ( V123 + V123 - V1 ) + FTY * ( U123 + U123 - U1 )
C
         W1(IELEM) = ( 2 * G1*T11 +     G2*T12 +     G3*T13 ) * XS120
         W2(IELEM) = (     G1*T12 + 2 * G2*T22 +     G3*T23 ) * XS120
         W3(IELEM) = (     G1*T13 +     G2*T23 + 2 * G3*T33 ) * XS120
C
3     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
       IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
       IF (LNG.EQ.1) WRITE(LU,110) IELMG,SG%NAME
       IF (LNG.EQ.1) WRITE(LU,200) IELMU,SU%NAME
       IF (LNG.EQ.1) WRITE(LU,300)
       IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,111) IELMG,SG%NAME
       IF (LNG.EQ.1) WRITE(LU,201) IELMU,SU%NAME
       IF (LNG.EQ.1) WRITE(LU,301)
100    FORMAT(1X,'VC09AA (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,
     *        1X,'NOM REEL : ',A6)
110    FORMAT(1X,'DISCRETISATION DE G : ',1I6,
     *        1X,'NOM REEL : ',A6)
200    FORMAT(1X,'DISCRETISATION DE U : ',1I6,
     *        1X,'NOM REEL : ',A6)
300    FORMAT(1X,'CAS NON PREVU')
101    FORMAT(1X,'VC09AA (BIEF) :',/,
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
 
 
