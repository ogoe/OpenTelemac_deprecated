C                       ******************
                        SUBROUTINE VC11AA2
C                       ******************
C
     *( XMUL,SF,SG,SH,F,G,H,XEL,YEL,
     *  IKLE1,IKLE2,IKLE3,NELEM,NELMAX,W1,W2,W3 , ICOORD )
C
C***********************************************************************
C BIEF VERSION 5.1           09/12/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CALCUL DES TERMES SUIVANTS
C
C            ( EXEMPLE DE LA COMPOSANTE X QUI CORRESPOND A ICOORD =1 )
C
C                       /                  DF
C    VEC(I)  =  XMUL   /  ( G . H .  P  *( --  )) D(OMEGA)
C                     /OMEGA          I    DX
C
C
C
C    P   EST UNE BASE LINEAIRE
C     I
C
C    F EST UN VECTEUR DE DISCRETISATION P1 
C    G EST UN VECTEUR DE DISCRETISATION P1 DISCONTINU
C    H EST UN VECTEUR DE DISCRETISATION P0
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
C |      ICOORD    | -->|  NUMERO DE LA COORDONNEE.
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
      USE BIEF, EX_VC11AA2 => VC11AA2
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
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,*),YEL(NELMAX,*)
      DOUBLE PRECISION,INTENT(INOUT)::W1(NELMAX),W2(NELMAX),W3(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
C     STRUCTURES DE F,G ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF,SG,SH
      DOUBLE PRECISION, INTENT(IN) :: F(*),G(*),H(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF,IELMG,IELMH
      DOUBLE PRECISION XSUR24 ,F1,F2,F3,G1,G2,G3,X2,X3,Y2,Y3,KSAT
C
C-----------------------------------------------------------------------
C
      XSUR24= XMUL / 24.D0
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
      IELMG=SG%ELM
      IELMH=SH%ELM
C
C-----------------------------------------------------------------------
C
C     F P1, G P1 DISCONTINU, H P0
C
      IF (IELMG.EQ.10.AND.SG%DIMDISC.EQ.11.AND.SG%DIM2.EQ.3
     &   .AND.IELMF.EQ.11.AND.
     &   IELMH.EQ.10) THEN
C
C  COORDONNEE X
C
      IF(ICOORD.EQ.1) THEN
C
      DO 1 IELEM = 1 , NELEM
C
        KSAT=H(IELEM)
        G1 = G(IELEM)
        G2 = G(IELEM+NELEM)
        G3 = G(IELEM+2*NELEM)
        F1 = F(IKLE1(IELEM))
        F2 = F(IKLE2(IELEM))
        F3 = F(IKLE3(IELEM))
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        W1(IELEM)=(Y2*(-G3*F3+G3*F1-G2*F3+G2*F1-2*G1*F3+2*G1*F1)+Y3*(
     *             G3*F2-G3*F1+G2*F2-G2*F1+2*G1*F2-2*G1*F1))* XSUR24
        W2(IELEM)=(Y2*(-G3*F3+G3*F1-2*G2*F3+2*G2*F1-G1*F3+G1*F1)+Y3*(
     *             G3*F2-G3*F1+2*G2*F2-2*G2*F1+G1*F2-G1*F1))* XSUR24
        W3(IELEM)=(Y2*(-2*G3*F3+2*G3*F1-G2*F3+G2*F1-G1*F3+G1*F1)+Y3*(
     *             2*G3*F2-2*G3*F1+G2*F2-G2*F1+G1*F2-G1*F1))* XSUR24
C
        W1(IELEM)=KSAT*W1(IELEM)
        W2(IELEM)=KSAT*W2(IELEM)
        W3(IELEM)=KSAT*W3(IELEM)
C
1     CONTINUE
C
      ELSEIF(ICOORD.EQ.2) THEN
C
C  COORDONNEE Y
C
      DO 2 IELEM = 1 , NELEM
C
        KSAT=H(IELEM)
        G1 = G(IELEM)
        G2 = G(IELEM+NELEM)
        G3 = G(IELEM+2*NELEM)
        F1 = F(IKLE1(IELEM))
        F2 = F(IKLE2(IELEM))
        F3 = F(IKLE3(IELEM))
        X2 = XEL(IELEM,2)
        X3 = XEL(IELEM,3)
C
        W1(IELEM)=(X2*(G3*F3-G3*F1+G2*F3-G2*F1+2*G1*F3-2*G1*F1)+X3*(-
     *             G3*F2+G3*F1-G2*F2+G2*F1-2*G1*F2+2*G1*F1)) * XSUR24
        W2(IELEM)=(X2*(G3*F3-G3*F1+2*G2*F3-2*G2*F1+G1*F3-G1*F1)+X3*(-
     *             G3*F2+G3*F1-2*G2*F2+2*G2*F1-G1*F2+G1*F1)) * XSUR24
        W3(IELEM)=(X2*(2*G3*F3-2*G3*F1+G2*F3-G2*F1+G1*F3-G1*F1)+X3*(-
     *             2*G3*F2+2*G3*F1-G2*F2+G2*F1-G1*F2+G1*F1)) * XSUR24
C
        W1(IELEM)=KSAT*W1(IELEM)
        W2(IELEM)=KSAT*W2(IELEM)
        W3(IELEM)=KSAT*W3(IELEM)
C
2     CONTINUE
C
      ELSE
C
          IF (LNG.EQ.1) WRITE(LU,20) ICOORD
          IF (LNG.EQ.2) WRITE(LU,21) ICOORD
20        FORMAT(1X,'VC11AA2 (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
21        FORMAT(1X,'VC11AA2 (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
          CALL PLANTE(0)
          STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
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
100    FORMAT(1X,'VC11AA2 (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,
     *        1X,'NOM REEL : ',A6)
200    FORMAT(1X,'DISCRETISATION DE G : ',1I6,
     *        1X,'NOM REEL : ',A6)
300    FORMAT(1X,'CAS NON PREVU')
101    FORMAT(1X,'VC11AA2 (BIEF) :',/,
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
