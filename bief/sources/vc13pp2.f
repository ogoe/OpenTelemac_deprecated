C                       ******************
                        SUBROUTINE VC13PP2
C                       ******************
C
     *( XMUL,SF,F,X,Y,Z,
     *  IKLE1,IKLE2,IKLE3,IKLE4,IKLE5,IKLE6,NELEM,NELMAX,
     *  W1,W2,W3,W4,W5,W6,ICOORD)
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
C    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C    ICI F EST UNE FONCTION DE X ET Y (PAS DE Z)
C
C    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      USE BIEF !, EX_VC13PP2 => VC13PP2
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX,ICOORD
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX),IKLE3(NELMAX)
      INTEGER, INTENT(IN) :: IKLE4(NELMAX),IKLE5(NELMAX),IKLE6(NELMAX)
C
      DOUBLE PRECISION, INTENT(IN)    :: X(*),Y(*),Z(*)
      DOUBLE PRECISION, INTENT(INOUT) ::W1(NELMAX),W2(NELMAX),W3(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) ::W4(NELMAX),W5(NELMAX),W6(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL
C
C     STRUCTURE DE F ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF
      DOUBLE PRECISION, INTENT(IN) :: F(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION XSUR48,F1,F2,F3
      DOUBLE PRECISION X2,X3,Y2,Y3,H1,H2,H3,H123,FX6,FY6
      INTEGER I1,I2,I3,I4,I5,I6,IELEM,IELMF
C
C-----------------------------------------------------------------------
C
      XSUR48  = XMUL/48.D0
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
C
C=======================================================================
C
C     F LINEAIRE
C
      IF(IELMF.EQ.41) THEN
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
         I5 = IKLE5(IELEM)
         I6 = IKLE6(IELEM)
C
         F1 = F(I1)
         F2 = F(I2)
         F3 = F(I3)
C
C  COORDONNEES REELLES DES POINTS DE L'ELEMENT ( ORIGINE EN 1 )
C
         Y2  =  Y(I2) - Y(I1)
         Y3  =  Y(I3) - Y(I1)
C
         FX6 = Y2*(F1-F3) + Y3*(F2-F1) 
C         
         H1 = Z(I4) - Z(I1)
         H2 = Z(I5) - Z(I2)
         H3 = Z(I6) - Z(I3)
         H123 = H1 + H2 + H3
C
         W1(IELEM) = XSUR48 * (H123+H1) * FX6
         W2(IELEM) = XSUR48 * (H123+H2) * FX6
         W3(IELEM) = XSUR48 * (H123+H3) * FX6
C
         W4(IELEM) = W1(IELEM)
         W5(IELEM) = W2(IELEM)
         W6(IELEM) = W3(IELEM)
C
3     CONTINUE
C
      ELSEIF(ICOORD.EQ.2) THEN
C
C-----------------------------------------------------------------------
C
C  DERIVEE SUIVANT Y
C
      DO 4 IELEM = 1 , NELEM
C
         I1 = IKLE1(IELEM)
         I2 = IKLE2(IELEM)
         I3 = IKLE3(IELEM)
         I4 = IKLE4(IELEM)
         I5 = IKLE5(IELEM)
         I6 = IKLE6(IELEM)
C
         F1 = F(I1)
         F2 = F(I2)
         F3 = F(I3)
C
C  COORDONNEES REELLES DES POINTS DE L'ELEMENT ( ORIGINE EN 1 )
C
         X2  =  X(I2) - X(I1)
         X3  =  X(I3) - X(I1)
C
         FY6 = X2*(F3-F1) + X3*(F1-F2)  
C         
         H1 = Z(I4) - Z(I1)
         H2 = Z(I5) - Z(I2)
         H3 = Z(I6) - Z(I3)
         H123 = H1 + H2 + H3
C
         W1(IELEM) = XSUR48 * (H123+H1) * FY6
         W2(IELEM) = XSUR48 * (H123+H2) * FY6
         W3(IELEM) = XSUR48 * (H123+H3) * FY6
C
         W4(IELEM) = W1(IELEM)
         W5(IELEM) = W2(IELEM)
         W6(IELEM) = W3(IELEM)
C
4        CONTINUE
C
      ELSE
C
C-----------------------------------------------------------------------
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
200       FORMAT(1X,'VC13PP2 (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'VC13PP2 (BIEF) : IMPOSSIBLE COMPONENT ',
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
101    FORMAT(1X,'VC13PP2 (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,' CAS NON PREVU',/,
     *        1X,'NOM REEL DE F : ',A6)
102    FORMAT(1X,'VC13PP2 (BIEF) :',/,
     *        1X,'DISCRETISATION OF F : ',1I6,' NOT IMPLEMENTED',/,
     *        1X,'REAL NAME OF F: ',A6)
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
C=======================================================================
C
      RETURN
      END
