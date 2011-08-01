C                       *****************
                        SUBROUTINE VC11TT
C                       *****************
C
     *( XMUL,SF,SG,F,G,X,Y,Z,IKLE1,IKLE2,IKLE3,IKLE4,NELEM,NELMAX,
     *  W1,W2,W3,W4,ICOORD )
C
C***********************************************************************
C BIEF VERSION 5.3         25/03/02   J-M HERVOUET (LNH) 01 30 87 80 18
C                                        
C ARNAUD DESITTER - UNIVERSITY OF BRISTOL - APRIL 1998
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C            (EXEMPLE DE LA COMPOSANTE X QUI CORRESPOND A ICOORD = 1)
C
C
C                       /            DF
C    VEC(I)  =  XMUL   /  ( G  P  *( --  )) D(OMEGA)
C                     /OMEGA    I    DX
C
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
C                      F(NELMAX,4)
C
C    ATTENTION : LE RESULTAT EST DANS W SOUS FORME NON ASSEMBLEE.
C
C                MAILLAGE REEL
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
C |      W1,2,3,4  |<-- |  VECTEUR RESULTAT SOUS FORME NON ASSEMBLEE.
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
C-----------------------------------------------------------------------
C
      USE BIEF
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
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*),Z(*),XMUL
      DOUBLE PRECISION, INTENT(INOUT) :: W1(NELMAX),W2(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W3(NELMAX),W4(NELMAX)
C
C     STRUCTURES DE F,G,H,U,V,W ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF,SG
      DOUBLE PRECISION, INTENT(IN) :: F(*),G(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C VARIABLES LOCALES
C
      INTEGER IELEM,IELMF,IELMG
      DOUBLE PRECISION F1,F2,F3,F4
      DOUBLE PRECISION G1,G2,G3,G4,X2,X3,X4,Y2,Y3,Y4,Z2,Z3,Z4
      INTEGER I1,I2,I3,I4
C
      DOUBLE PRECISION XSUR120,F2MF1,F3MF1,F4MF1,G2MG1,G3MG1,G4MG1
C
C-----------------------------------------------------------------------
C INITIALISATIONS
C
      XSUR120 = XMUL/120.D0
C
      IELMF = SF%ELM
      IELMG = SG%ELM
C
C-----------------------------------------------------------------------
C     F ET G LINEAIRES
C
      IF ((IELMF.EQ.31.AND.IELMG.EQ.31).OR.
     *    (IELMF.EQ.51.AND.IELMG.EQ.51)     ) THEN
C
      IF (ICOORD.EQ.1) THEN
C
C-----------------------------------------------------------------------
C  DERIVEE SUIVANT X
C
      DO  IELEM = 1 , NELEM
C
         I1 = IKLE1(IELEM)
         I2 = IKLE2(IELEM)
         I3 = IKLE3(IELEM)
         I4 = IKLE4(IELEM)
C
         F1 = F(I1)
         F2 = F(I2)
         F3 = F(I3)
         F4 = F(I4)
C
         G1 = G(I1)
         G2 = G(I2)
         G3 = G(I3)
         G4 = G(I4)
C
         F2MF1 = F2-F1
         F3MF1 = F3-F1
         F4MF1 = F4-F1
         G2MG1 = G2-G1
         G3MG1 = G3-G1
         G4MG1 = G4-G1
C
C  COORDONNEES REELLES DES POINTS DE L'ELEMENT ( ORIGINE EN 1 )
C
         Y2  =  Y(I2) - Y(I1)
         Y3  =  Y(I3) - Y(I1)
         Y4  =  Y(I4) - Y(I1)
         Z2  =  Z(I2) - Z(I1)
         Z3  =  Z(I3) - Z(I1)
         Z4  =  Z(I4) - Z(I1)
C
         W1(IELEM) = (
     # (5*F2MF1*G1+F2MF1*G2MG1+F2MF1*G3MG1+F2MF1*G4MG1)*(Y3*Z4-Y4*Z3)
     #+(5*F3MF1*G1+F3MF1*G2MG1+F3MF1*G3MG1+F3MF1*G4MG1)*(Z2*Y4-Y2*Z4)
     #+(5*F4MF1*G1+F4MF1*G2MG1+F4MF1*G3MG1+F4MF1*G4MG1)*(Y2*Z3-Z2*Y3)
     #               ) * XSUR120
C
         W2(IELEM) = (
     #-F4MF1*Z2*Y3*G4MG1+F4MF1*Y2*Z3*G4MG1+F3MF1*Z2*Y4*G4MG1
     #-F3MF1*Y2*Z4*G4MG1+F2MF1*Y3*Z4*G4MG1-F2MF1*Y4*Z3*G4MG1
     #+5*F2MF1*Y3*Z4*G1+2*F2MF1*Y3*Z4*G2MG1+F2MF1*Y3*Z4*G3MG1
     #-5*F2MF1*Y4*Z3*G1-2*F2MF1*Y4*Z3*G2MG1-F2MF1*Y4*Z3*G3MG1
     #-5*F3MF1*Y2*Z4*G1-2*F3MF1*Y2*Z4*G2MG1-F3MF1*Y2*Z4*G3MG1
     #+5*F3MF1*Z2*Y4*G1+2*F3MF1*Z2*Y4*G2MG1+F3MF1*Z2*Y4*G3MG1
     #+5*F4MF1*Y2*Z3*G1+2*F4MF1*Y2*Z3*G2MG1+F4MF1*Y2*Z3*G3MG1
     #-5*F4MF1*Z2*Y3*G1-2*F4MF1*Z2*Y3*G2MG1-F4MF1*Z2*Y3*G3MG1
     #               ) * XSUR120
         W3(IELEM) = (
     #-(-F2MF1*Y3*Z4+F2MF1*Y4*Z3+F3MF1*Y2*Z4
     #  -F3MF1*Z2*Y4-F4MF1*Y2*Z3+F4MF1*Z2*Y3)
     #               *(2*G3MG1+G2MG1+G4MG1+5*G1)
     #               ) * XSUR120
         W4(IELEM) = (
     #-2*F4MF1*Z2*Y3*G4MG1
     #+2*F4MF1*Y2*Z3*G4MG1
     #+2*F3MF1*Z2*Y4*G4MG1
     #-2*F3MF1*Y2*Z4*G4MG1
     #+2*F2MF1*Y3*Z4*G4MG1
     #-2*F2MF1*Y4*Z3*G4MG1
     #+5*F2MF1*Y3*Z4*G1
     #+F2MF1*Y3*Z4*G2MG1
     #+F2MF1*Y3*Z4*G3MG1-5*F2MF1*Y4*Z3*G1-F2MF1*Y4*Z3*G2MG1
     #-F2MF1*Y4*Z3*G3MG1-5*F3MF1*Y2*Z4*G1-F3MF1*Y2*Z4*G2MG1
     #-F3MF1*Y2*Z4*G3MG1+5*F3MF1*Z2*Y4*G1+F3MF1*Z2*Y4*G2MG1
     #+F3MF1*Z2*Y4*G3MG1+5*F4MF1*Y2*Z3*G1+F4MF1*Y2*Z3*G2MG1
     #+F4MF1*Y2*Z3*G3MG1-5*F4MF1*Z2*Y3*G1-F4MF1*Z2*Y3*G2MG1
     #-F4MF1*Z2*Y3*G3MG1
     #               ) * XSUR120
C
      ENDDO
C
      ELSE IF (ICOORD.EQ.2) THEN
C
C-----------------------------------------------------------------------
C  DERIVEE SUIVANT Y
C
      DO   IELEM = 1 , NELEM
C
         I1 = IKLE1(IELEM)
         I2 = IKLE2(IELEM)
         I3 = IKLE3(IELEM)
         I4 = IKLE4(IELEM)
C
         F1 = F(I1)
         F2 = F(I2)
         F3 = F(I3)
         F4 = F(I4)
C
         G1 = G(I1)
         G2 = G(I2)
         G3 = G(I3)
         G4 = G(I4)
C
         F2MF1 = F2-F1
         F3MF1 = F3-F1
         F4MF1 = F4-F1
         G2MG1 = G2-G1
         G3MG1 = G3-G1
         G4MG1 = G4-G1
C
C  COORDONNEES REELLES DES POINTS DE L'ELEMENT ( ORIGINE EN 1 )
C
         X2  =  X(I2) - X(I1)
         X3  =  X(I3) - X(I1)
         X4  =  X(I4) - X(I1)
         Z2  =  Z(I2) - Z(I1)
         Z3  =  Z(I3) - Z(I1)
         Z4  =  Z(I4) - Z(I1)
C
         W1(IELEM) = (
     #-F2MF1*X3*Z4*G2MG1+F3MF1*X2*Z4*G3MG1+5*F3MF1*X2*Z4*G1
     #-F2MF1*X3*Z4*G3MG1-5*F2MF1*X3*Z4*G1
     #+F2MF1*X4*Z3*G2MG1+F2MF1*X4*Z3*G3MG1+5*F2MF1*X4*Z3*G1
     #+F3MF1*X2*Z4*G2MG1+F4MF1*Z2*X3*G2MG1+F4MF1*Z2*X3*G3MG1
     #+5*F4MF1*Z2*X3*G1-F4MF1*X2*Z3*G2MG1-F4MF1*X2*Z3*G3MG1
     #-5*F4MF1*X2*Z3*G1-F3MF1*Z2*X4*G2MG1-F3MF1*Z2*X4*G3MG1
     #-5*F3MF1*Z2*X4*G1+F2MF1*X4*Z3*G4MG1+F3MF1*X2*Z4*G4MG1
     #-F4MF1*X2*Z3*G4MG1-F3MF1*Z2*X4*G4MG1-F2MF1*X3*Z4*G4MG1
     #+F4MF1*Z2*X3*G4MG1 ) * XSUR120
         W2(IELEM) = (
     #         -2*F2MF1*X3*Z4*G2MG1+F3MF1*X2*Z4*G3MG1+5*F3MF1*X2*Z4*G1-F
     #2MF1*X3*Z4*G3MG1-5*F2MF1*X3*Z4*G1+2*F2MF1*X4*Z3*G2MG1+F2MF1*X4*Z3*
     #G3MG1+5*F2MF1*X4*Z3*G1+2*F3MF1*X2*Z4*G2MG1+2*F4MF1*Z2*X3*G2MG1+F4M
     #F1*Z2*X3*G3MG1+5*F4MF1*Z2*X3*G1-2*F4MF1*X2*Z3*G2MG1-F4MF1*X2*Z3*G3
     #MG1-5*F4MF1*X2*Z3*G1-2*F3MF1*Z2*X4*G2MG1-F3MF1*Z2*X4*G3MG1-5*F3MF1
     #*Z2*X4*G1+F2MF1*X4*Z3*G4MG1+F3MF1*X2*Z4*G4MG1-F4MF1*X2*Z3*G4MG1-F3
     #MF1*Z2*X4*G4MG1-F2MF1*X3*Z4*G4MG1+F4MF1*Z2*X3*G4MG1 ) * XSUR120
         W3(IELEM) = (
     #         -(F2MF1*X3*Z4-F2MF1*X4*Z3-F3MF1*X2*Z4+F3MF1*Z2*X4+F4MF1*X
     #2*Z3-F4MF1*Z2*X3)*(2*G3MG1+G2MG1+G4MG1+5*G1) ) * XSUR120
         W4(IELEM) = (
     #         -F2MF1*X3*Z4*G2MG1+F3MF1*X2*Z4*G3MG1+5*F3MF1*X2*Z4*G1-F2M
     #F1*X3*Z4*G3MG1-5*F2MF1*X3*Z4*G1+F2MF1*X4*Z3*G2MG1+F2MF1*X4*Z3*G3MG
     #1+5*F2MF1*X4*Z3*G1+F3MF1*X2*Z4*G2MG1+F4MF1*Z2*X3*G2MG1+F4MF1*Z2*X3
     #*G3MG1+5*F4MF1*Z2*X3*G1-F4MF1*X2*Z3*G2MG1-F4MF1*X2*Z3*G3MG1-5*F4MF
     #1*X2*Z3*G1-F3MF1*Z2*X4*G2MG1-F3MF1*Z2*X4*G3MG1-5*F3MF1*Z2*X4*G1+2*
     #F2MF1*X4*Z3*G4MG1+2*F3MF1*X2*Z4*G4MG1-2*F4MF1*X2*Z3*G4MG1-2*F3MF1*
     #Z2*X4*G4MG1-2*F2MF1*X3*Z4*G4MG1+2*F4MF1*Z2*X3*G4MG1 ) * XSUR120
C
      ENDDO
C
      ELSE IF (ICOORD.EQ.3) THEN
C-----------------------------------------------------------------------
C  DERIVEE SUIVANT Z
C
      DO   IELEM = 1 , NELEM
C
         I1 = IKLE1(IELEM)
         I2 = IKLE2(IELEM)
         I3 = IKLE3(IELEM)
         I4 = IKLE4(IELEM)
C
         F1 = F(I1)
         F2 = F(I2)
         F3 = F(I3)
         F4 = F(I4)
C
         G1 = G(I1)
         G2 = G(I2)
         G3 = G(I3)
         G4 = G(I4)
C
         F2MF1 = F2-F1
         F3MF1 = F3-F1
         F4MF1 = F4-F1
         G2MG1 = G2-G1
         G3MG1 = G3-G1
         G4MG1 = G4-G1
C
C  COORDONNEES REELLES DES POINTS DE L'ELEMENT
C
         X2  =  X(I2) - X(I1)
         X3  =  X(I3) - X(I1)
         X4  =  X(I4) - X(I1)
         Y2  =  Y(I2) - Y(I1)
         Y3  =  Y(I3) - Y(I1)
         Y4  =  Y(I4) - Y(I1)
C
         W1(IELEM) = (
     #         5*F2MF1*X3*Y4*G1+F2MF1*X3*Y4*G2MG1+F2MF1*X3*Y4*G3MG1-5*F2
     #MF1*X4*Y3*G1-F2MF1*X4*Y3*G2MG1-F2MF1*X4*Y3*G3MG1-5*F3MF1*X2*Y4*G1-
     #F3MF1*X2*Y4*G2MG1-F3MF1*X2*Y4*G3MG1+5*F3MF1*Y2*X4*G1+F3MF1*Y2*X4*G
     #2MG1+F3MF1*Y2*X4*G3MG1+5*F4MF1*X2*Y3*G1+F4MF1*X2*Y3*G2MG1-5*F4MF1*
     #Y2*X3*G1-F4MF1*Y2*X3*G2MG1-F4MF1*Y2*X3*G3MG1+F4MF1*X2*Y3*G3MG1-F4M
     #F1*Y2*X3*G4MG1-F3MF1*X2*Y4*G4MG1+F4MF1*X2*Y3*G4MG1+F2MF1*X3*Y4*G4M
     #G1+F3MF1*Y2*X4*G4MG1-F2MF1*X4*Y3*G4MG1 ) * XSUR120
         W2(IELEM) = (
     #         5*F2MF1*X3*Y4*G1+2*F2MF1*X3*Y4*G2MG1+F2MF1*X3*Y4*G3MG1-5*
     #F2MF1*X4*Y3*G1-2*F2MF1*X4*Y3*G2MG1-F2MF1*X4*Y3*G3MG1-5*F3MF1*X2*Y4
     #*G1-2*F3MF1*X2*Y4*G2MG1-F3MF1*X2*Y4*G3MG1+5*F3MF1*Y2*X4*G1+2*F3MF1
     #*Y2*X4*G2MG1+F3MF1*Y2*X4*G3MG1+5*F4MF1*X2*Y3*G1+2*F4MF1*X2*Y3*G2MG
     #1-5*F4MF1*Y2*X3*G1-2*F4MF1*Y2*X3*G2MG1-F4MF1*Y2*X3*G3MG1+F4MF1*X2*
     #Y3*G3MG1-F4MF1*Y2*X3*G4MG1-F3MF1*X2*Y4*G4MG1+F4MF1*X2*Y3*G4MG1+F2M
     #F1*X3*Y4*G4MG1+F3MF1*Y2*X4*G4MG1-F2MF1*X4*Y3*G4MG1 ) * XSUR120 
         W3(IELEM) = (
     #         -(-F2MF1*X3*Y4+F2MF1*X4*Y3+F3MF1*X2*Y4-F3MF1*Y2*X4-F4MF1*
     #X2*Y3+F4MF1*Y2*X3)*(2*G3MG1+G2MG1+G4MG1+5*G1) ) * XSUR120
         W4(IELEM) = (
     #         5*F2MF1*X3*Y4*G1+F2MF1*X3*Y4*G2MG1+F2MF1*X3*Y4*G3MG1-5*F2
     #MF1*X4*Y3*G1-F2MF1*X4*Y3*G2MG1-F2MF1*X4*Y3*G3MG1-5*F3MF1*X2*Y4*G1-
     #F3MF1*X2*Y4*G2MG1-F3MF1*X2*Y4*G3MG1+5*F3MF1*Y2*X4*G1+F3MF1*Y2*X4*G
     #2MG1+F3MF1*Y2*X4*G3MG1+5*F4MF1*X2*Y3*G1+F4MF1*X2*Y3*G2MG1-5*F4MF1*
     #Y2*X3*G1-F4MF1*Y2*X3*G2MG1-F4MF1*Y2*X3*G3MG1+F4MF1*X2*Y3*G3MG1-2*F
     #4MF1*Y2*X3*G4MG1-2*F3MF1*X2*Y4*G4MG1+2*F4MF1*X2*Y3*G4MG1+2*F2MF1*X
     #3*Y4*G4MG1+2*F3MF1*Y2*X4*G4MG1-2*F2MF1*X4*Y3*G4MG1 ) * XSUR120
C
      ENDDO
C
      ELSE
C
C-----------------------------------------------------------------------
C
         IF (LNG.EQ.1) WRITE(LU,200) ICOORD
         IF (LNG.EQ.2) WRITE(LU,201) ICOORD
 200     FORMAT(1X,'VC11TT (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *        1I6,' VERIFIER ICOORD')
 201     FORMAT(1X,'VC11TT (BIEF) : IMPOSSIBLE COMPONENT ',
     *        1I6,' CHECK ICOORD')
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C ERREUR
C
      ELSE
C
C-----------------------------------------------------------------------
C
         IF (LNG.EQ.1) WRITE(LU,1100) IELMF,SF%NAME
         IF (LNG.EQ.1) WRITE(LU,1200) IELMG,SG%NAME
         IF (LNG.EQ.1) WRITE(LU,1300)
         IF (LNG.EQ.2) WRITE(LU,1101) IELMF,SF%NAME
         IF (LNG.EQ.2) WRITE(LU,1201) IELMG,SG%NAME
         IF (LNG.EQ.2) WRITE(LU,1301)
         CALL PLANTE(1)
         STOP
 1100    FORMAT(1X,'VC11TT (BIEF) :',/,
     *          1X,'DISCRETISATION DE F : ',1I6,
     *          1X,'NOM REEL : ',A6)
 1200    FORMAT(1X,'DISCRETISATION DE G : ',1I6,
     *          1X,'NOM REEL : ',A6)
 1300    FORMAT(1X,'CAS NON PREVU')
 1101    FORMAT(1X,'VC11TT (BIEF) :',/,
     *          1X,'DISCRETIZATION OF F:',1I6,
     *          1X,'REAL NAME: ',A6)
 1201    FORMAT(1X,'DISCRETIZATION OF G:',1I6,
     *          1X,'REAL NAME: ',A6)
 1301    FORMAT(1X,'CASE NOT IMPLEMENTED')
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
