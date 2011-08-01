C                       *****************
                        SUBROUTINE VC11BB
C                       *****************
C
     *( XMUL,SF,SG,F,G,XEL,YEL,
     *  IKLE1,IKLE2,IKLE3,IKLE4,NELEM,NELMAX,W1,W2,W3,W4 , ICOORD )
C
C***********************************************************************
C BIEF VERSION 5.1           09/12/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CALCUL DES TERMES SUIVANTS
C
C            (EXEMPLE DE LA COMPOSANTE X QUI CORRESPOND A ICOORD = 1)
C
C                       /            DF
C    VEC(I)  =  XMUL   /  ( G  P  *( --  )) D(OMEGA)
C                     /OMEGA    I    DX
C
C
C
C    P   EST UNE BASE QUASI-BULLE
C     I
C
C    F EST UN VECTEUR DE DISCRETISATION P1 OU AUTRE
C    G EST UN VECTEUR DE DISCRETISATION P1 OU AUTRE
C
C    NOTE IMPORTANTE : SI F EST DE TYPE P0, LE RESULTAT EST NUL
C                      ICI, SI F EST P0, CELA SIGNIFIE QUE F EST
C                      P1, MAIS DONNEE PAR ELEMENTS.
C                      LE DIMENSIONNEMENT DE F DOIT ETRE ALORS :
C                      F(NELMAX,3)
C
C    ATTENTION : LE RESULTAT EST DANS W SOUS FORME NON ASSEMBLEE.
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
C |      ICOORD    | -->|  NUMERO DE LA COORDONNEE POUR LA DERIVEE.
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
      USE BIEF   !, EX_VC11BB => VC11BB
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
      DOUBLE PRECISION, INTENT(IN)    :: XEL(NELMAX,*),YEL(NELMAX,*)
      DOUBLE PRECISION, INTENT(INOUT) :: W1(NELMAX),W2(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W3(NELMAX),W4(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL
C
C     STRUCTURES DE F,G ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN)   :: SF,SG
      DOUBLE PRECISION, INTENT(IN) :: F(*),G(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF,IELMG
C
      DOUBLE PRECISION F1,F2,F3,F4,G1,G2,G3,G4,X2,X3,Y2,Y3,AUX
      DOUBLE PRECISION XSUR72 ,XSU216,XSUR18
C
C-----------------------------------------------------------------------
C
      XSUR72= XMUL / 72.D0
      XSU216= XMUL /216.D0
      XSUR18= XMUL / 18.D0
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
      IELMG=SG%ELM
C
C-----------------------------------------------------------------------
C
C     F ET G LINEAIRES
C
      IF(IELMF.EQ.11.AND.IELMG.EQ.11) THEN
C
C  COORDONNEE X
C
      IF(ICOORD.EQ.1) THEN
C
      DO 1 IELEM = 1 , NELEM
C
        F1 = F(IKLE1(IELEM))
        F2 = F(IKLE2(IELEM)) - F1
        F3 = F(IKLE3(IELEM)) - F1
        G1 = G(IKLE1(IELEM))
        G2 = G(IKLE2(IELEM))
        G3 = G(IKLE3(IELEM))
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
        AUX = F2*Y3 - F3*Y2
C
        W1(IELEM)=(5*G3+5*G2+14*G1)*AUX*XSU216
        W2(IELEM)=(5*G3+14*G2+5*G1)*AUX*XSU216
        W3(IELEM)=(14*G3+5*G2+5*G1)*AUX*XSU216
        W4(IELEM)=(G3+G2+G1)       *AUX*XSUR18
C
1     CONTINUE
C
      ELSEIF(ICOORD.EQ.2) THEN
C
C  COORDONNEE Y
C
      DO 2 IELEM = 1 , NELEM
C
        F1 = F(IKLE1(IELEM))
        F2 = F(IKLE2(IELEM)) - F1
        F3 = F(IKLE3(IELEM)) - F1
        G1 = G(IKLE1(IELEM))
        G2 = G(IKLE2(IELEM))
        G3 = G(IKLE3(IELEM))
        X2 = XEL(IELEM,2)
        X3 = XEL(IELEM,3)
        AUX = X2*F3 - X3*F2
C
        W1(IELEM)=AUX*(5*G3+5*G2+14*G1)*XSU216
        W2(IELEM)=AUX*(5*G3+14*G2+5*G1)*XSU216
        W3(IELEM)=AUX*(14*G3+5*G2+5*G1)*XSU216
        W4(IELEM)=AUX*(G3+G2+G1)       *XSUR18
C
2     CONTINUE
C
      ELSE
C
          IF (LNG.EQ.1) WRITE(LU,20) ICOORD
          IF (LNG.EQ.2) WRITE(LU,21) ICOORD
20        FORMAT(1X,'VC11BB (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
21        FORMAT(1X,'VC11BB (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
          CALL PLANTE(0)
          STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C     F ET G QUASI-BULLE
C
      ELSEIF(IELMF.EQ.12.AND.IELMG.EQ.12) THEN
C
C  COORDONNEE X
C
      IF(ICOORD.EQ.1) THEN
C
      DO 3 IELEM = 1 , NELEM
C
        F1 = F(IKLE1(IELEM))
        F2 = F(IKLE2(IELEM)) - F1
        F3 = F(IKLE3(IELEM)) - F1
        F4 = F(IKLE4(IELEM)) - F1
        G1 = G(IKLE1(IELEM))
        G2 = G(IKLE2(IELEM))
        G3 = G(IKLE3(IELEM))
        G4 = G(IKLE4(IELEM))
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        W1(IELEM)=(3*((Y3-Y2)*G4+2*(Y3-Y2)*G1+G3*Y3-G2*Y2)*F4-(G3+
     *   G4+2*G1)*(Y3+Y2)*F3+(G4+G2+2*G1)*(Y3+Y2)*F2)*XSUR72
        W2(IELEM)=(((2*Y3-Y2)*G3+(Y3+Y2)*G1+3*G4*Y3+6*G2*Y3)*F2-3
     *   *((Y3-Y2)*G3+G4*Y3+2*G2*Y3+G1*Y2)*F4+(G3+G4+2*G2)*(Y3-
     *   2*Y2)*F3)*XSUR72
        W3(IELEM)=(-(((Y3+Y2)*G1-(Y3-2*Y2)*G2+6*G3*Y2+3*G4*Y2)*F3+
     *   3*((Y3-Y2)*G2-2*G3*Y2-G4*Y2-G1*Y3)*F4-(2*G3+G4+G2)*(2*
     *   Y3-Y2)*F2))*XSUR72
        W4(IELEM)=(((2*Y3-Y2)*G3+(Y3+Y2)*G1+6*G4*Y3+3*G2*Y3)*F2-((
     *   Y3+Y2)*G1-(Y3-2*Y2)*G2+3*G3*Y2+6*G4*Y2)*F3+3*((Y3-Y2)
     *   *G1+G3*Y2-G2*Y3)*F4)*XSUR72
C
3     CONTINUE
C
      ELSEIF(ICOORD.EQ.2) THEN
C
C  COORDONNEE Y
C
      DO 4 IELEM = 1 , NELEM
C
        F1 = F(IKLE1(IELEM))
        F2 = F(IKLE2(IELEM)) - F1
        F3 = F(IKLE3(IELEM)) - F1
        F4 = F(IKLE4(IELEM)) - F1
        G1 = G(IKLE1(IELEM))
        G2 = G(IKLE2(IELEM))
        G3 = G(IKLE3(IELEM))
        G4 = G(IKLE4(IELEM))
        X2 = XEL(IELEM,2)
        X3 = XEL(IELEM,3)
C
        W1(IELEM)=(-(3*((G3+G4+2*G1)*X3-(G4+G2+2*G1)*X2)*F4-(X2+X3
     *   )*(G3+G4+2*G1)*F3+(X2+X3)*(G4+G2+2*G1)*F2))*XSUR72
        W2(IELEM)=(-(((2*G3+3*G4+6*G2+G1)*X3-(G3-G1)*X2)*F2-3*((
     *   G3+G4+2*G2)*X3-(G3-G1)*X2)*F4-(2*X2-X3)*(G3+G4+2*G2)*
     *   F3))*XSUR72
        W3(IELEM)=(((6*G3+3*G4+2*G2+G1)*X2-(G2-G1)*X3)*F3-3*((2*
     *   G3+G4+G2)*X2-(G2-G1)*X3)*F4+(X2-2*X3)*(2*G3+G4+G2)*F2)*XSUR72
        W4(IELEM)=(((3*G3+6*G4+2*G2+G1)*X2-(G2-G1)*X3)*F3-((2*G3+
     *   6*G4+3*G2+G1)*X3-(G3-G1)*X2)*F2-3*((G3-G1)*X2-(G2-G1)*
     *   X3)*F4)*XSUR72
C
4     CONTINUE
C
      ELSE
C
          IF (LNG.EQ.1) WRITE(LU,20) ICOORD
          IF (LNG.EQ.2) WRITE(LU,21) ICOORD
          CALL PLANTE(0)
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
       IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
       IF (LNG.EQ.1) WRITE(LU,200) IELMG,SG%NAME
       IF (LNG.EQ.1) WRITE(LU,300)
       IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,201) IELMG,SG%NAME
       IF (LNG.EQ.2) WRITE(LU,301)
100    FORMAT(1X,'VC11BB (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,
     *        1X,'NOM REEL : ',A6)
200    FORMAT(1X,'DISCRETISATION DE G : ',1I6,
     *        1X,'NOM REEL : ',A6)
300    FORMAT(1X,'CAS NON PREVU')
101    FORMAT(1X,'VC11BB (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F:',1I6,
     *        1X,'REAL NAME: ',A6)
201    FORMAT(1X,'DISCRETIZATION OF G:',1I6,
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
 
 
