C                       *****************
                        SUBROUTINE VC13CC
C                       *****************
C
     *( XMUL,SF,F,XEL,YEL,
     *  IKLE1,IKLE2,IKLE3,IKLE4,IKLE5,IKLE6,
     *  NELEM,NELMAX,
     *  W1,W2,W3,W4,W5,W6,ICOORD )
C
C***********************************************************************
C BIEF VERSION 5.9         01/07/08   A FROEHLY (MATMECA) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C                    /            --->
C    VEC(I) = XMUL  /    PSI(I) * GRAD(F)  D(OMEGA)
C                  /OMEGA
C
C
C    PSI(I) EST UNE BASE DE TYPE P2
C
C    NOTE IMPORTANTE : SI F EST DE TYPE P0, LE RESULTAT EST NUL
C                      ICI, SI F EST P0, CELA SIGNIFIE QUE F EST
C                      P1, MAIS DONNEE PAR ELEMENTS.
C                      LE DIMENSIONNEMENT DE F DOIT ETRE ALORS :
C                      F(NELMAX,6)
C
C
C    ATTENTION : LE RESULTAT EST DANS W  SOUS FORME NON ASSEMBLEE.
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
C |      W1,2,3... |<-- |  VECTEUR RESULTAT SOUS FORME NON ASSEMBLEE.
C |      ICOORD    | -->|  COORDONNEE SUIVANT LAQUELLE ON DERIVE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF !, EX_VC13CC => VC13CC
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
      DOUBLE PRECISION, INTENT(IN)    :: XEL(NELMAX,*),YEL(NELMAX,*)
      DOUBLE PRECISION, INTENT(INOUT) :: W1(NELMAX),W2(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W3(NELMAX),W4(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W5(NELMAX),W6(NELMAX)
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
      DOUBLE PRECISION F1,F2,F3,F4,F5,F6
      DOUBLE PRECISION X2,X3,Y2,Y3
      DOUBLE PRECISION XSUR6,XSUR30
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
      DISCF = SF%DIMDISC
      XSUR6 = XMUL / 6.D0 
      XSUR30= XMUL / 30.D0                 
C
C-----------------------------------------------------------------------
C     F DE DISCRETISATION P1
C-----------------------------------------------------------------------
C
      IF(IELMF.EQ.11) THEN
C
C================================
C  CAS DE LA DERIVEE SUIVANT X  =
C================================
C
        IF(ICOORD.EQ.1) THEN
C
        DO 1 IELEM = 1 , NELEM
C
        F1  = F(IKLE1(IELEM))
        F2  = F(IKLE2(IELEM))
        F3  = F(IKLE3(IELEM))
C
        Y2  =  YEL(IELEM,2)
        Y3  =  YEL(IELEM,3)
C
        W1(IELEM) = 0.D0
        W2(IELEM) = 0.D0
        W3(IELEM) = 0.D0
        W4(IELEM) = (Y2*(F1-F3)+Y3*(F2-F1))*XSUR6
        W5(IELEM) = W4(IELEM) 
        W6(IELEM) = W4(IELEM)
C
1       CONTINUE
C
C================================
C  CAS DE LA DERIVEE SUIVANT Y  =
C================================
C
      ELSEIF(ICOORD.EQ.2) THEN
C
        DO 2 IELEM = 1 , NELEM
C
        F1  = F(IKLE1(IELEM))
        F2  = F(IKLE2(IELEM))
        F3  = F(IKLE3(IELEM))
C
        X2  =  XEL(IELEM,2)
        X3  =  XEL(IELEM,3)
C
        W1(IELEM) = 0.D0
        W2(IELEM) = 0.D0
        W3(IELEM) = 0.D0
        W4(IELEM) = (X2*(F3-F1)+X3*(F1-F2))*XSUR6
        W5(IELEM) = W4(IELEM) 
        W6(IELEM) = W4(IELEM)
C
2     CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
          CALL PLANTE(1)
          STOP
        ENDIF
C
C-----------------------------------------------------------------------
C
C     ATTENTION : ICI F LINEAIRE MAIS DISCONTINUE ENTRE LES ELEMENTS
C
      ELSEIF(IELMF.EQ.10.AND.DISCF.EQ.11) THEN
C
C  COORDONNEE X
C
      IF(ICOORD.EQ.1) THEN
C
      DO 3 IELEM = 1 , NELEM
C
        F1  = F(IELEM)
        F2  = F(IELEM+NELMAX)
        F3  = F(IELEM+2*NELMAX)
C
        Y2  =  YEL(IELEM,2)
        Y3  =  YEL(IELEM,3)
C
        W1(IELEM) = 0.D0
        W2(IELEM) = 0.D0
        W3(IELEM) = 0.D0
        W4(IELEM) = (Y2*(F1-F3)+Y3*(F2-F1))*XSUR6
        W5(IELEM) = W4(IELEM) 
        W6(IELEM) = W4(IELEM)
C
3     CONTINUE
C
      ELSEIF(ICOORD.EQ.2) THEN
C
C  COORDONNEE Y
C
      DO 4 IELEM = 1 , NELEM
C
        F1  = F(IELEM)
        F2  = F(IELEM+NELMAX)
        F3  = F(IELEM+2*NELMAX)
C
        X2  =  XEL(IELEM,2)
        X3  =  XEL(IELEM,3)
C
        W1(IELEM) = 0.D0
        W2(IELEM) = 0.D0
        W3(IELEM) = 0.D0
        W4(IELEM) = (X2*(F3-F1)+X3*(F1-F2))*XSUR6
        W5(IELEM) = W4(IELEM) 
        W6(IELEM) = W4(IELEM)
C
4     CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
          CALL PLANTE(1)
          STOP
        ENDIF
C
C-----------------------------------------------------------------------
C     F P2
C-----------------------------------------------------------------------
C
      ELSEIF(IELMF.EQ.13) THEN
C
C================================
C  CAS DE LA DERIVEE SUIVANT X  =
C================================
C
      IF(ICOORD.EQ.1) THEN
C
        DO 5 IELEM = 1 , NELEM
C
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C        
        F1 = F(IKLE1(IELEM))
        F2 = F(IKLE2(IELEM))
        F3 = F(IKLE3(IELEM))
        F4 = F(IKLE4(IELEM))
        F5 = F(IKLE5(IELEM))
        F6 = F(IKLE6(IELEM))
C
        W1(IELEM) = ((2.D0*F1+F5-F4+F3-3.D0*F6)*Y2
     *            +  (F6-F5+3.D0*F4-F2-2.D0*F1)*Y3) * XSUR30
        W2(IELEM) = ((2.D0*F2-3.D0*F4+F1+F6-F5)*Y3
     *            +  (F3-F1-2.D0*F5+2.D0*F4   )*Y2) * XSUR30
        W3(IELEM) = ((3.D0*F6+F5-F4-F1-2.D0*F3)*Y2
     *            +  (F1-F2+2.D0*F5-2.D0*F6   )*Y3) * XSUR30
        W4(IELEM) = ((-4.D0*F6-8.D0*(F5-F4)+3.D0*F1+F3)*Y2
     *            +  (3.D0*(F2-F1)+4.D0*(F5-F6)       )*Y3) * XSUR30
        W5(IELEM) = ((4.D0*F6-8.D0*(F5-F4)-F1-3.D0*F3 )*Y2
     *            +  (3.D0*F2-4.D0*F4+F1+8.D0*(F5-F6) )*Y3) * XSUR30
        W6(IELEM) = ((4.D0*(F4-F5)+3.D0*(F1-F3)       )*Y2
     *            +  (4.D0*F4-3.D0*F1+8.D0*(F5-F6)-F2 )*Y3) * XSUR30
C
5     CONTINUE
C
C================================
C  CAS DE LA DERIVEE SUIVANT Y  =
C================================
C
      ELSEIF(ICOORD.EQ.2) THEN
C
        DO 6 IELEM = 1 , NELEM
C
        X2 = XEL(IELEM,2)
        X3 = XEL(IELEM,3)
C        
        F1  = F(IKLE1(IELEM))
        F2  = F(IKLE2(IELEM))
        F3  = F(IKLE3(IELEM))
        F4  = F(IKLE4(IELEM))
        F5  = F(IKLE5(IELEM))
        F6  = F(IKLE6(IELEM))
C
        W1(IELEM) = ((F4-F5-2.D0*F1+3.D0*F6-F3 )*X2
     *            +  (-3.D0*F4-F6+2.D0*F1+F2+F5)*X3) * XSUR30
        W2(IELEM) = ((F5-F1-2.D0*F2+3.D0*F4-F6)*X3
     *            +  (-F3+F1-2.D0*(F4-F5)     )*X2 ) * XSUR30
        W3(IELEM) = ((-3.D0*F6+2.D0*F3+F1+F4-F5)*X2
     *            +  (2.D0*(F6-F5)-F1+F2       )*X3) * XSUR30  
        W4(IELEM) = ((4.D0*F6-F3-3.D0*F1-8.D0*(F4-F5))*X2
     *            +  (4.D0*(F6-F5)+3.D0*(F1-F2)      )*X3) * XSUR30
        W5(IELEM) = ((3.D0*F3+F1-8.D0*(F4-F5)-4.D0*F6)*X2
     *            +  (8.D0*(F6-F5)+4.D0*F4-F1-3.D0*F2)*X3) * XSUR30
        W6(IELEM) = ((3.D0*(F3-F1)-4.D0*(F4-F5)      )*X2
     *            +  (8.D0*(F6-F5)+3.D0*F1+F2-4.D0*F4)*X3) * XSUR30
C
6     CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
          CALL PLANTE(1)
          STOP
        ENDIF
C
C-----------------------------------------------------------------------
C     F P2 MAIS DISCONTINUE ENTRE LES ELEMENTS
C-----------------------------------------------------------------------
C
      ELSEIF(IELMF.EQ.10.AND.DISCF.EQ.13) THEN
C
C================================
C  CAS DE LA DERIVEE SUIVANT X  =
C================================
C
      IF(ICOORD.EQ.1) THEN
C
        DO 7 IELEM = 1 , NELEM
C        
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        F1  = F(IELEM         )
        F2  = F(IELEM+  NELMAX)
        F3  = F(IELEM+2*NELMAX)
        F4  = F(IELEM+3*NELMAX)
        F5  = F(IELEM+4*NELMAX)
        F6  = F(IELEM+5*NELMAX)
C
        W1(IELEM) = ((2.D0*F1+F5-F4+F3-3.D0*F6)*Y2
     *            +  (F6-F5+3.D0*F4-F2-2.D0*F1)*Y3) * XSUR30
        W2(IELEM) = ((2.D0*F2-3.D0*F4+F1+F6-F5)*Y3
     *            +  (F3-F1-2.D0*F5+2.D0*F4   )*Y2) * XSUR30
        W3(IELEM) = ((3.D0*F6+F5-F4-F1-2.D0*F3)*Y2
     *            +  (F1-F2+2.D0*F5-2.D0*F6   )*Y3) * XSUR30
        W4(IELEM) = ((-4.D0*F6-8.D0*(F5-F4)+3.D0*F1+F3)*Y2
     *            +  (3.D0*(F2-F1)+4.D0*(F5-F6)       )*Y3) * XSUR30
        W5(IELEM) = ((4.D0*F6-8.D0*(F5-F4)-F1-3.D0*F3 )*Y2
     *            +  (3.D0*F2-4.D0*F4+F1+8.D0*(F5-F6) )*Y3) * XSUR30
        W6(IELEM) = ((4.D0*(F4-F5)+3.D0*(F1-F3)       )*Y2
     *            +  (4.D0*F4-3.D0*F1+8.D0*(F5-F6)-F2 )*Y3) * XSUR30
C
7     CONTINUE
C
C================================
C  CAS DE LA DERIVEE SUIVANT Y  =
C================================
C
      ELSEIF(ICOORD.EQ.2) THEN
C
      DO 8 IELEM = 1 , NELEM
C
        X2 = XEL(IELEM,2)
        X3 = XEL(IELEM,3)
C        
        F1  = F(IELEM         )
        F2  = F(IELEM+  NELMAX)
        F3  = F(IELEM+2*NELMAX)
        F4  = F(IELEM+3*NELMAX)
        F5  = F(IELEM+4*NELMAX)
        F6  = F(IELEM+5*NELMAX)
C
        W1(IELEM) = ((F4-F5-2.D0*F1+3.D0*F6-F3 )*X2
     *            +  (-3.D0*F4-F6+2.D0*F1+F2+F5)*X3) * XSUR30
        W2(IELEM) = ((F5-F1-2.D0*F2+3.D0*F4-F6)*X3
     *            +  (-F3+F1-2.D0*(F4-F5)     )*X2 ) * XSUR30
        W3(IELEM) = ((-3.D0*F6+2.D0*F3+F1+F4-F5)*X2
     *            +  (2.D0*(F6-F5)-F1+F2       )*X3) * XSUR30  
        W4(IELEM) = ((4.D0*F6-F3-3.D0*F1-8.D0*(F4-F5))*X2
     *            +  (4.D0*(F6-F5)+3.D0*(F1-F2)      )*X3) * XSUR30
        W5(IELEM) = ((3.D0*F3+F1-8.D0*(F4-F5)-4.D0*F6)*X2
     *            +  (8.D0*(F6-F5)+4.D0*F4-F1-3.D0*F2)*X3) * XSUR30
        W6(IELEM) = ((3.D0*(F3-F1)-4.D0*(F4-F5)      )*X2
     *            +  (8.D0*(F6-F5)+3.D0*F1+F2-4.D0*F4)*X3) * XSUR30
C
8     CONTINUE
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
          CALL PLANTE(1)
          STOP
        ENDIF
C
C-----------------------------------------------------------------------
C      AUTRES
C      ELSEIF
C-----------------------------------------------------------------------
C
      ELSE
C
       IF(LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
       IF(LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
100    FORMAT(1X,'VC13CC (BIEF) :',/,
     *        1X,'DISCRETISATION DE F NON PREVUE : ',1I6,
     *        1X,'NOM REEL : ',A6)
101    FORMAT(1X,'VC13CC (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F NOT AVAILABLE:',1I6,
     *        1X,'REAL NAME: ',A6)
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
200       FORMAT(1X,'VC13CC (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'VC13CC (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
C
C-----------------------------------------------------------------------
C
      RETURN
      END
