C                       *****************
                        SUBROUTINE VC01AA
C                       *****************
C
     *( XMUL,SF,F,SURFAC,
     *  IKLE1,IKLE2,IKLE3,NELEM,NELMAX,
     *  W1,W2,W3 )
C
C***********************************************************************
C BIEF VERSION 5.1           29/10/99    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C                    /
C    VEC(I) = XMUL  /    PSI(I) * F  D(OMEGA)
C                  /OMEGA
C
C
C    PSI(I) EST UNE BASE DE TYPE TRIANGLE P1
C
C    F EST UN VECTEUR DE DISCRETISATION P0, P1 OU P1 DISCONTINU
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
      USE BIEF, EX_VC01AA => VC01AA
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
      DOUBLE PRECISION, INTENT(INOUT) :: W1(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W2(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W3(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: SURFAC(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL
C
C     STRUCTURE DE F ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN)   :: SF
      DOUBLE PRECISION, INTENT(IN) :: F(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER IELEM,IELMF
      DOUBLE PRECISION XSUR03,XSUR12,F1,F2,F3,F123,COEF
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
C
C-----------------------------------------------------------------------
C
C     F LINEAIRE
C
      IF(IELMF.EQ.11) THEN
C
      XSUR12 = XMUL / 12.D0
C
      DO 3 IELEM = 1 , NELEM
C
        F1  = F(IKLE1(IELEM))
        F2  = F(IKLE2(IELEM))
        F3  = F(IKLE3(IELEM))
        F123  = F1 + F2 + F3
C
        COEF = XSUR12 * SURFAC(IELEM)
C
        W1(IELEM) = COEF * ( F123 + F1 )
        W2(IELEM) = COEF * ( F123 + F2 )
        W3(IELEM) = COEF * ( F123 + F3 )
C
3     CONTINUE
C
C-----------------------------------------------------------------------
C
C     F CONSTANTE PAR ELEMENT
C
      ELSEIF(IELMF.EQ.10.AND.SF%DIM2.EQ.1) THEN
C
      XSUR03 = XMUL / 3.D0
C
      DO 4 IELEM = 1 , NELEM
C
        W1(IELEM) = XSUR03 * SURFAC(IELEM) * F(IELEM)
        W2(IELEM) = W1(IELEM)
        W3(IELEM) = W1(IELEM)
C
4     CONTINUE
C
C-----------------------------------------------------------------------
C
C     F P1 DISCONTINUE
C
      ELSEIF(IELMF.EQ.10.AND.SF%DIM2.EQ.3.AND.SF%DIMDISC.EQ.11) THEN
C
      XSUR12 = XMUL / 12.D0
C
      DO 5 IELEM = 1 , NELEM
C
        F1  = F(IELEM)
        F2  = F(IELEM+NELEM)
        F3  = F(IELEM+2*NELEM)
        F123  = F1 + F2 + F3
C
        COEF = XSUR12 * SURFAC(IELEM)
C
        W1(IELEM) = COEF * ( F123 + F1 )
        W2(IELEM) = COEF * ( F123 + F2 )
        W3(IELEM) = COEF * ( F123 + F3 )
C
5     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
       IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
100    FORMAT(1X,'VC01AA (BIEF) :',/,
     *        1X,'DISCRETISATION DE F NON PREVUE : ',1I6,
     *        1X,'NOM REEL : ',A6)
101    FORMAT(1X,'VC01AA (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F NOT AVAILABLE:',1I6,
     *        1X,'REAL NAME: ',A6)
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
