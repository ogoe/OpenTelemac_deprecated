C                       *****************
                        SUBROUTINE VC13TT
C                       *****************
C
     *( XMUL,SF,F,X,Y,Z,
     *  IKLE1,IKLE2,IKLE3,IKLE4,NELEM,NELMAX,
     *  W1,W2,W3,W4,ICOORD , FORMUL )
C
C***********************************************************************
C BIEF VERSION 5.3           25/03/02  J-M HERVOUET (LNH) 01 30 87 80 18
C                                         
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
C                MAILLAGE REEL.
C
C    ICI L'ELEMENT EST UN TETRAEDRE LINEAIRE
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
      USE BIEF !, EX_VC13TT => VC13TT
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
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*),Z(*)
      DOUBLE PRECISION, INTENT(INOUT) :: W1(NELMAX),W2(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W3(NELMAX),W4(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
C     STRUCTURE DE F ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF
      DOUBLE PRECISION, INTENT(IN) :: F(*)
      CHARACTER(LEN=16), INTENT(IN) :: FORMUL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF
      DOUBLE PRECISION X2,X3,X4,Y2,Y3,Y4,Z2,Z3,Z4,XSUR24
      INTEGER I1,I2,I3,I4
C
C-----------------------------------------------------------------------
C
      XSUR24 = XMUL/24.D0
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
C
C=======================================================================
C
C     F LINEAIRE
C
      IF(IELMF.EQ.31.OR.IELMF.EQ.51) THEN
C
      IF(ICOORD.EQ.1) THEN
C
C-----------------------------------------------------------------------
C
C  DERIVEE SUIVANT X
C
      DO 3 IELEM = 1 , NELEM
C
         I1 = IKLE1(IELEM)
         I2 = IKLE2(IELEM)
         I3 = IKLE3(IELEM)
         I4 = IKLE4(IELEM)
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
         W1(IELEM)=(  (F(I2)-F(I1))*(Y3*Z4-Y4*Z3)
     #               +(F(I3)-F(I1))*(Z2*Y4-Y2*Z4)
     #               +(F(I4)-F(I1))*(Y2*Z3-Z2*Y3)  )*XSUR24
C
         W2(IELEM)=W1(IELEM)
         W3(IELEM)=W1(IELEM)
         W4(IELEM)=W1(IELEM)
C
3     CONTINUE
C
      ELSEIF(ICOORD.EQ.2) THEN
C
C-----------------------------------------------------------------------
C
C  DERIVEE SUIVANT Y
C
      DO IELEM = 1 , NELEM
C
         I1 = IKLE1(IELEM)
         I2 = IKLE2(IELEM)
         I3 = IKLE3(IELEM)
         I4 = IKLE4(IELEM)
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
         W1(IELEM)=(  (F(I2)-F(I1))*(X4*Z3-X3*Z4)
     #               +(F(I3)-F(I1))*(X2*Z4-Z2*X4)
     #               +(F(I4)-F(I1))*(Z2*X3-X2*Z3)  )*XSUR24
C
         W2(IELEM)=W1(IELEM)
         W3(IELEM)=W1(IELEM)
         W4(IELEM)=W1(IELEM)
C
      ENDDO
C
      ELSEIF(ICOORD.EQ.3) THEN
C
C-----------------------------------------------------------------------
C
C  DERIVEE SUIVANT Z
C
      DO 5 IELEM = 1 , NELEM
C
         I1 = IKLE1(IELEM)
         I2 = IKLE2(IELEM)
         I3 = IKLE3(IELEM)
         I4 = IKLE4(IELEM)
C
C  COORDONNEES REELLES DES POINTS DE L'ELEMENT ( ORIGINE EN 1 )
C
         X2  =  X(I2) - X(I1)
         X3  =  X(I3) - X(I1)
         X4  =  X(I4) - X(I1)
         Y2  =  Y(I2) - Y(I1)
         Y3  =  Y(I3) - Y(I1)
         Y4  =  Y(I4) - Y(I1)
C
         W1(IELEM)=(  (F(I2)-F(I1))*(X3*Y4-X4*Y3)
     #               +(F(I3)-F(I1))*(Y2*X4-X2*Y4)
     #               +(F(I4)-F(I1))*(X2*Y3-Y2*X3)  )*XSUR24
C
         W2(IELEM)=W1(IELEM)
         W3(IELEM)=W1(IELEM)
         W4(IELEM)=W1(IELEM)
C
5     CONTINUE
C
      ELSE
C
C-----------------------------------------------------------------------
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
200       FORMAT(1X,'VC13TT (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'VC13TT (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
          CALL PLANTE(1)
          STOP
C
      ENDIF
C
C=======================================================================
C
      ELSE
C
C=======================================================================
C
       IF (LNG.EQ.1) WRITE(LU,101) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,102) IELMF,SF%NAME
101    FORMAT(1X,'VC13TT (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,' CAS NON PREVU',/,
     *        1X,'NOM REEL DE F : ',A6)
102    FORMAT(1X,'VC13TT (BIEF) :',/,
     *        1X,'DISCRETISATION OF F : ',1I6,' NOT IMPLEMENTED',/,
     *        1X,'REAL NAME OF F: ',A6)
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
C=======================================================================
C
C  INCONSISTANCES HYDROSTATIQUES
C
C     IF(FORMUL(6:6).EQ.'2') THEN
C
C     DO IELEM = 1 , NELEM
C
C        I1 = IKLE1(IELEM)
C        I2 = IKLE2(IELEM)
C        I3 = IKLE3(IELEM)
C        I4 = IKLE4(IELEM)
C        I5 = IKLE5(IELEM)
C        I6 = IKLE6(IELEM)
C
C        IF(MAX(Z(I1),Z(I2),Z(I3)).GT.MIN(Z(I4),Z(I5),Z(I6))) THEN
C          W1(IELEM)=0.D0
C          W2(IELEM)=0.D0
C          W3(IELEM)=0.D0
C          W4(IELEM)=0.D0
C          W5(IELEM)=0.D0
C          W6(IELEM)=0.D0
C        ENDIF
C
C     ENDDO
C
C     ENDIF
C
C=======================================================================
C
      RETURN
      END
