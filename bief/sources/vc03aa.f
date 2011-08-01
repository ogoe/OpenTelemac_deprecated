C                       *****************
                        SUBROUTINE VC03AA
C                       *****************
C
     *( XMUL,SF,SG,SH,SU,SV,F,G,H,U,V,
     *  XEL,YEL,SURFAC,
     *  IKLE1,IKLE2,IKLE3,NELEM,NELMAX,
     *  W1,W2,W3 )
C
C***********************************************************************
C BIEF VERSION 5.1           09/12/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C
C                    /                     DF      DF
C      V  =  XMUL   /   K GRAD(PSII) * ( U --  + V -- )   D(OMEGA)
C       I          /OMEGA                  DX      DY
C
C
C    PSI(I) EST UNE BASE DE TYPE TRIANGLE P1
C
C    F, U ET V SONT DES VECTEURS
C    K EST UN VECTEUR DE COMPOSANTES G ET H
C
C
C    ATTENTION : LE RESULTAT EST DANS W SOUS FORME NON ASSEMBLEE.
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
      USE BIEF, EX_VC03AA => VC03AA
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
      DOUBLE PRECISION, INTENT(IN)   :: XEL(NELMAX,*),YEL(NELMAX,*)
      DOUBLE PRECISION, INTENT(INOUT):: W1(NELMAX),W2(NELMAX),W3(NELMAX)
      DOUBLE PRECISION, INTENT(IN)   :: SURFAC(NELMAX)
      DOUBLE PRECISION, INTENT(IN)   :: XMUL
C
C     STRUCTURES DE F,G,H,U,V ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF,SG,SH,SU,SV
      DOUBLE PRECISION, INTENT(IN) :: F(*),G(*),H(*),U(*),V(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER IELEM,IELMF,IELMG,IELMU,IELMV,IELMH
      DOUBLE PRECISION X2,Y2,X3,Y3,F1,F2,F3,U123,V123
      DOUBLE PRECISION WX1,WX2,WX3,WY1,WY2,WY3,XSUR12,COEF
C
C-----------------------------------------------------------------------
C
      XSUR12 = XMUL / 12.D0
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
      IELMG=SG%ELM
      IELMU=SU%ELM
      IELMV=SV%ELM
      IELMH=SH%ELM
C
C-----------------------------------------------------------------------
C
C     F LINEAIRE, G ET H P0 ET U ,V LINEAIRES
C
      IF(       IELMF.EQ.11
     *     .AND.IELMG.EQ.10
     *     .AND.IELMH.EQ.10
     *     .AND.IELMU.EQ.11
     *     .AND.IELMV.EQ.11  ) THEN
C
      DO 3 IELEM = 1 , NELEM
C
      X2 = XEL(IELEM,2)
      X3 = XEL(IELEM,3)
      Y2 = YEL(IELEM,2)
      Y3 = YEL(IELEM,3)
C
      F1 = F(IKLE1(IELEM))
      F2 = F(IKLE2(IELEM)) - F1
      F3 = F(IKLE3(IELEM)) - F1
C
C     F1 DESORMAIS NUL (SEUL LE GRADIENT DE F INTERVIENT)
C
      U123 = U(IKLE1(IELEM)) + U(IKLE2(IELEM)) + U(IKLE3(IELEM))
      V123 = V(IKLE1(IELEM)) + V(IKLE2(IELEM)) + V(IKLE3(IELEM))
C
      WX1 = ( - F2*X3*Y2 + F2*X3*Y3 + F3*X2*Y2 - F3*X2*Y3 ) * V123
     *    + ( + F2*Y2*Y3 - F2*Y3*Y3 - F3*Y2*Y2 + F3*Y2*Y3 ) * U123
C
      WY1 = (   F2*X2*X3 - F3*X2*X2 - F2*X3*X3 + F3*X2*X3 ) * V123
     *    + ( - F2*X2*Y3 + F2*X3*Y3 + F3*X2*Y2 - F3*X3*Y2 ) * U123
C
      WX2 = Y3 * ( (F3*X2-F2*X3) * V123 + (F2*Y3-F3*Y2) * U123 )
C
      WY2 = X3 * ( (F2*X3-F3*X2) * V123 + (F3*Y2-F2*Y3) * U123 )
C
      WX3 = Y2 * ( (F2*X3-F3*X2) * V123 + (F3*Y2-F2*Y3) * U123 )
C
      WY3 = X2 * ( (F3*X2-F2*X3) * V123 + (F2*Y3-F3*Y2) * U123 )
C
      COEF = XSUR12 / SURFAC(IELEM)
C
      W1(IELEM) = ( WX1*G(IELEM) + WY1*H(IELEM) ) * COEF
      W2(IELEM) = ( WX2*G(IELEM) + WY2*H(IELEM) ) * COEF
      W3(IELEM) = ( WX3*G(IELEM) + WY3*H(IELEM) ) * COEF
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
100    FORMAT(1X,'VC03AA (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,
     *        1X,'NOM REEL : ',A6)
110    FORMAT(1X,'DISCRETISATION DE G : ',1I6,
     *        1X,'NOM REEL : ',A6)
200    FORMAT(1X,'DISCRETISATION DE U : ',1I6,
     *        1X,'NOM REEL : ',A6)
300    FORMAT(1X,'CAS NON PREVU')
101    FORMAT(1X,'VC03AA (BIEF) :',/,
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
 
 
