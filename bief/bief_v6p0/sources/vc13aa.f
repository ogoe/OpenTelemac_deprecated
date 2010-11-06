C                       *****************
                        SUBROUTINE VC13AA
C                       *****************
C
     *( XMUL,SF,F,XEL,YEL,
     *  IKLE1,IKLE2,IKLE3,IKLE4,NELEM,NELMAX,W1,W2,W3 , ICOORD )
C
C***********************************************************************
C BIEF VERSION 5.1           09/12/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C            (EXEMPLE DE LA COMPOSANTE X QUI CORRESPOND A ICOORD = 1)
C
C
C                       /            DF
C    VEC(I)  =  XMUL   /     ( P  *( --  )) D(OMEGA)
C                     /OMEGA    I    DX
C
C    P   EST UNE BASE LINEAIRE
C     I
C
C    F EST UN VECTEUR DE DISCRETISATION P1 OU AUTRE
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
C |      W1,2,3    |<-- |  VECTEUR RESULTAT SOUS FORME NON ASSEMBLEE.
C |      ICOORD    | -->|  COORDONNEE SUIVANT LAQUELLE ON DERIVE.
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
      USE BIEF  !, EX_VC13AA => VC13AA
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
      DOUBLE PRECISION, INTENT(INOUT) ::W1(NELMAX),W2(NELMAX),W3(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL
C
C     STRUCTURE DE F ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF
      DOUBLE PRECISION, INTENT(IN) :: F(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF,DISCF
      DOUBLE PRECISION XSUR6,XSUR18,F1,F2,F3,F4,X2,X3,Y2,Y3
C
C-----------------------------------------------------------------------
C
      XSUR6 = XMUL / 6.D0
      XSUR18= XMUL /18.D0
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
      DISCF = SF%DIMDISC
C
C-----------------------------------------------------------------------
C
C     F LINEAIRE
C
      IF(IELMF.EQ.11) THEN
C
C  COORDONNEE X
C
      IF(ICOORD.EQ.1) THEN
C
      DO 1 IELEM = 1 , NELEM
C
        W1(IELEM) = ( YEL(IELEM,2) *
     *                (F(IKLE1(IELEM))-F(IKLE3(IELEM)))
     *              + YEL(IELEM,3) *
     *                (F(IKLE2(IELEM))-F(IKLE1(IELEM))) ) * XSUR6
        W2(IELEM) = W1(IELEM)
        W3(IELEM) = W1(IELEM)
C
1     CONTINUE
C
      ELSEIF(ICOORD.EQ.2) THEN
C
C  COORDONNEE Y
C
      DO 2 IELEM = 1 , NELEM
C
         W1(IELEM) = ( XEL(IELEM,2) *
     *                 (F(IKLE3(IELEM))-F(IKLE1(IELEM)))
     *               + XEL(IELEM,3) *
     *                 (F(IKLE1(IELEM))-F(IKLE2(IELEM))) ) * XSUR6
         W2(IELEM)  =  W1(IELEM)
         W3(IELEM)  =  W1(IELEM)
C
2     CONTINUE
C
      ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
200       FORMAT(1X,'VC13AA (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'VC13AA (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
          CALL PLANTE(0)
          STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C     F QUASI-BULLE
C
      ELSEIF(IELMF.EQ.12) THEN
C
C  COORDONNEE X
C
      IF(ICOORD.EQ.1) THEN
C
      DO 5 IELEM = 1 , NELEM
C
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        F1 = F(IKLE1(IELEM))
        F2 = F(IKLE2(IELEM)) - F1
        F3 = F(IKLE3(IELEM)) - F1
        F4 = F(IKLE4(IELEM)) - F1
C
        W1(IELEM)=(Y2*(-3*F4-2*F3+F2)+Y3*(3*F4-F3+2*F2)) * XSUR18
        W2(IELEM)=(-3*Y2*F3+Y3*(-3*F4+F3+4*F2)) * XSUR18
        W3(IELEM)=(Y2*(3*F4-4*F3-F2)+3*Y3*F2) * XSUR18
C
5     CONTINUE
C
      ELSEIF(ICOORD.EQ.2) THEN
C
C  COORDONNEE Y
C
      DO 6 IELEM = 1 , NELEM
C
        X2 = XEL(IELEM,2)
        X3 = XEL(IELEM,3)
C
        F1 = F(IKLE1(IELEM))
        F2 = F(IKLE2(IELEM)) - F1
        F3 = F(IKLE3(IELEM)) - F1
        F4 = F(IKLE4(IELEM)) - F1
C
        W1(IELEM)=(X2*(3*F4+2*F3-F2)+X3*(-3*F4+F3-2*F2)) * XSUR18
        W2(IELEM)=(3*X2*F3+X3*(3*F4-F3-4*F2)) * XSUR18
        W3(IELEM)=(X2*(-3*F4+4*F3+F2)-3*X3*F2) * XSUR18
C
6     CONTINUE
C
      ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
          CALL PLANTE(0)
          STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C     ATTENTION : ICI F LINEAIRE MAIS DISCONTINUE ENTRE LES ELEMENTS
C
      ELSEIF(IELMF.EQ.10.AND.DISCF.EQ.12) THEN
C
C  COORDONNEE X
C
      IF(ICOORD.EQ.1) THEN
C
      DO 7 IELEM = 1 , NELEM
C
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        F1 = F(IELEM)
        F2 = F(IELEM+  NELMAX)-F1
        F3 = F(IELEM+2*NELMAX)-F1
        F4 = F(IELEM+3*NELMAX)-F1
C
        W1(IELEM)=(Y2*(-3*F4-2*F3+F2)+Y3*(3*F4-F3+2*F2)) * XSUR18
        W2(IELEM)=(-3*Y2*F3+Y3*(-3*F4+F3+4*F2)) * XSUR18
        W3(IELEM)=(Y2*(3*F4-4*F3-F2)+3*Y3*F2) * XSUR18
C
7     CONTINUE
C
      ELSEIF(ICOORD.EQ.2) THEN
C
C  COORDONNEE Y
C
      DO 8 IELEM = 1 , NELEM
C
        X2 = XEL(IELEM,2)
        X3 = XEL(IELEM,3)
C
        F1 = F(IELEM)
        F2 = F(IELEM+  NELMAX)-F1
        F3 = F(IELEM+2*NELMAX)-F1
        F4 = F(IELEM+3*NELMAX)-F1
C
        W1(IELEM)=(X2*(3*F4+2*F3-F2)+X3*(-3*F4+F3-2*F2)) * XSUR18
        W2(IELEM)=(3*X2*F3+X3*(3*F4-F3-4*F2)) * XSUR18
        W3(IELEM)=(X2*(-3*F4+4*F3+F2)-3*X3*F2) * XSUR18
C
8     CONTINUE
C
      ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,200) ICOORD
          CALL PLANTE(0)
          STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(IELMF.EQ.10.AND.DISCF.EQ.11) THEN
C
C  COORDONNEE X
C
      IF(ICOORD.EQ.1) THEN
C
      DO 3 IELEM = 1 , NELEM
C
        W1(IELEM) = ( YEL(IELEM,2) * (F(IELEM)-F(IELEM+2*NELMAX))
     *              + YEL(IELEM,3) * (F(IELEM+NELMAX)-F(IELEM)))* XSUR6
        W2(IELEM) = W1(IELEM)
        W3(IELEM) = W1(IELEM)
C
3     CONTINUE
C
      ELSEIF(ICOORD.EQ.2) THEN
C
C  COORDONNEE Y
C
      DO 4 IELEM = 1 , NELEM
C
         W1(IELEM) = ( XEL(IELEM,2) * (F(IELEM+2*NELMAX)-F(IELEM))
     *               + XEL(IELEM,3) * (F(IELEM)-F(IELEM+NELMAX)))*XSUR6
         W2(IELEM)  =  W1(IELEM)
         W3(IELEM)  =  W1(IELEM)
C
4     CONTINUE
C
      ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,200) ICOORD
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
       IF (LNG.EQ.1) WRITE(LU,101) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,102) IELMF,SF%NAME
101    FORMAT(1X,'VC13AA (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,' CAS NON PREVU',/,
     *        1X,'NOM REEL DE F : ',A6)
102    FORMAT(1X,'VC13AA (BIEF) :',/,
     *        1X,'DISCRETISATION OF F : ',1I6,' NOT IMPLEMENTED',/,
     *        1X,'REAL NAME OF F: ',A6)
       CALL PLANTE(0)
       STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
