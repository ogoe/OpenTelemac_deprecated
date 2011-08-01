C                       *****************
                        SUBROUTINE VC01FF
C                       *****************
C
     *( XMUL,SF,F,X,Y,Z,
     *  IKLE1,IKLE2,IKLE3,IKLE4,NBOR,NELEM,NELMAX,W1,W2,W3,W4 )
C
C***********************************************************************
C BIEF VERSION 5.1           05/02/91    J-M HERVOUET (LNH) 30 87 80 18
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
C    PSI(I) EST UNE BASE DE TYPE SEGMENT P1
C
C    F EST UN VECTEUR DE TYPE IELMF
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
C***********************************************************************
C
      USE BIEF, EX_VC01FF => VC01FF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: NBOR(*)
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX)
      INTEGER, INTENT(IN) :: IKLE3(NELMAX),IKLE4(NELMAX)
C
      DOUBLE PRECISION, INTENT(IN)    :: X(*),Y(*),Z(*)
      DOUBLE PRECISION, INTENT(INOUT) :: W1(NELMAX),W2(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W3(NELMAX),W4(NELMAX)      
C
C     STRUCTURE DE F ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN)   :: SF
      DOUBLE PRECISION, INTENT(IN) :: F(*)
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER IELEM,IELMF,I1,I2,I3,I4
      DOUBLE PRECISION XSUR72,H1,H2,HT,AL,F1,F2,F3,F4
C
      INTRINSIC SQRT
C
C***********************************************************************
C
      IELMF=SF%ELM
C
C-----------------------------------------------------------------------
C
C     F LINEAIRE PAR FACE DE BORD
C
      IF(IELMF.EQ.71) THEN
C
         XSUR72 = XMUL/72.D0
C
C   BOUCLE SUR LES FACES DE BORD
C
         DO 1 IELEM = 1,NELEM
C
C  NUMEROTATION GLOBALE DES SOMMETS DE LA FACE
C
            I1 = IKLE1(IELEM)
            I2 = IKLE2(IELEM)
            I3 = IKLE3(IELEM)
            I4 = IKLE4(IELEM)
C
            AL = SQRT((X(NBOR(I2))-X(NBOR(I1)))**2
     *               +(Y(NBOR(I2))-Y(NBOR(I1)))**2) * XSUR72
C
            H1 = Z(NBOR(I4)) - Z(NBOR(I1))
            H2 = Z(NBOR(I3)) - Z(NBOR(I2))
            HT = H1 + H2
            H1 = H1 + H1 + HT
            H2 = H2 + H2 + HT
C
            F1 = F(I1) + F(I1) + F(I4)
            F2 = F(I2) + F(I2) + F(I3)
            F3 = F(I2) + F(I3) + F(I3)
            F4 = F(I1) + F(I4) + F(I4)
C
            W1(IELEM) = (F1*H1+F2*HT)*AL
            W2(IELEM) = (F1*HT+F2*H2)*AL
            W3(IELEM) = (F4*HT+F3*H2)*AL
            W4(IELEM) = (F4*H1+F3*HT)*AL
C
1        CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
       IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
100    FORMAT(1X,'VC01FF (BIEF) :',/,
     *        1X,'DISCRETISATION DE F NON PREVUE : ',1I6,
     *        1X,'NOM REEL : ',A6)
101    FORMAT(1X,'VC01FF (BIEF) :',/,
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
