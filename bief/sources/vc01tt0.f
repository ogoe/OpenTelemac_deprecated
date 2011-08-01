C                       *****************
                        SUBROUTINE VC01TT0
C                       *****************
C
     *( XMUL,SF,F,X,Y,Z,
     *  IKLE1,IKLE2,IKLE3,IKLE4,NELEM,NELMAX,W)
C
C***********************************************************************
C BIEF VERSION 5.3           22/03/02    J-M HERVOUET (LNH) 30 87 80 18
C                                        
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C                    /
C    VEC(I) = XMUL  /    PSI(I) * F  D(OMEGA)
C                  /OMEGA
C
C
C    PSI(I) EST UNE BASE DE TYPE TETRAEDRE P1
C
C    F EST UN VECTEUR DE TYPE IELMF
C
C
C    ATTENTION : LE RESULTAT EST DANS W  SOUS FORME NON ASSEMBLEE.
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
C |      IKLE1,... | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE.
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE.
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE.
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      W1,2,3    |<-- |  VECTEUR RESULTAT SOUS FORME NON ASSEMBLEE.
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
      USE BIEF, EX_VC01TT0 => VC01TT0
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX)
      INTEGER, INTENT(IN) :: IKLE3(NELMAX),IKLE4(NELMAX)
C
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*),Z(*)
      DOUBLE PRECISION, INTENT(INOUT) :: W(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
C     STRUCTURE DE F ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF
      DOUBLE PRECISION, INTENT(IN) :: F(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER IELEM,IELMF,DISCF
      DOUBLE PRECISION XSUR24,COEF,X2,X3,X4,Y2,Y3,Y4,Z2,Z3,Z4
      INTEGER I1,I2,I3,I4
C
C***********************************************************************
C
      IELMF=SF%ELM
      DISCF = SF%DIMDISC
C
C-----------------------------------------------------------------------
C
C   F LINEAIRE
C
      IF(IELMF.EQ.30) THEN
C
         XSUR24 = XMUL / 24.D0
C
         DO 3 IELEM = 1 , NELEM
C
           I1 = IKLE1(IELEM)
           I2 = IKLE2(IELEM)
           I3 = IKLE3(IELEM)
           I4 = IKLE4(IELEM)
C
           X2 = X(I2)-X(I1)
           X3 = X(I3)-X(I1)
           X4 = X(I4)-X(I1)
C
           Y2 = Y(I2)-Y(I1)
           Y3 = Y(I3)-Y(I1)
           Y4 = Y(I4)-Y(I1)
C
           Z2 = Z(I2)-Z(I1)
           Z3 = Z(I3)-Z(I1)
           Z4 = Z(I4)-Z(I1)
C
           COEF = XSUR24*
     #           (X2*Y3*Z4-X2*Y4*Z3-Y2*X3*Z4+Y2*X4*Z3+Z2*X3*Y4-Z2*X4*Y3)
C
           
           W(IELEM) = COEF * F(IELEM) 
           
C
3        CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSE
C
       IF (LNG.EQ.1) WRITE(LU,101) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,102) IELMF,SF%NAME
101    FORMAT(1X,'VC01TT0 (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,' CAS NON PREVU',/,
     *        1X,'NOM REEL DE F : ',A6)
102    FORMAT(1X,'VC01TT0 (BIEF):',/,
     *        1X,'DISCRETISATION OF F : ',1I6,' NOT IMPLEMENTED',/,
     *        1X,'REAL NAME OF F: ',A6)
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
