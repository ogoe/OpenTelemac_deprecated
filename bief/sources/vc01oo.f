C                       *****************
                        SUBROUTINE VC01OO
C                       *****************
C
     *( XMUL,SF,F,LGSEG,IKLE1,IKLE2,NBOR,NELEM,NELMAX,W1,W2 )
C
C***********************************************************************
C BIEF VERSION 5.9        20/03/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C                                     F  LEPEINTRE (LNH) 30 87 78 54
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
C**********************************************************************
C
      USE BIEF, EX_VC01OO => VC01OO
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX)
      INTEGER, INTENT(IN) :: NBOR(*)
C
      DOUBLE PRECISION, INTENT(INOUT) :: W1(NELMAX),W2(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: LGSEG(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL
C
C     STRUCTURE DE F ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF
      DOUBLE PRECISION, INTENT(IN) :: F(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER IELEM,IELMF
      DOUBLE PRECISION XSUR3,XSUR6,F1,F2,V1,V2
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
C
C-----------------------------------------------------------------------
C
C     F CONSTANTE PAR SEGMENTS
C
      IF(IELMF.EQ.0) THEN
C
      DO IELEM = 1,NELEM
        W1(IELEM) = 0.5D0*XMUL*F(IELEM)*LGSEG(IELEM)
        W2(IELEM) = W1(IELEM)
      ENDDO
C
C-----------------------------------------------------------------------
C
C     F LINEAIRE PAR SEGMENTS
C
      ELSEIF(IELMF.EQ.1) THEN
C
      XSUR3 = XMUL/3.D0
      XSUR6 = XMUL/6.D0
C
      DO IELEM = 1,NELEM
        F1 = F(IKLE1(IELEM))
        F2 = F(IKLE2(IELEM))
        V1 = ( F1*XSUR3 + F2*XSUR6 )
        V2 = ( F2*XSUR3 + F1*XSUR6 )
        W1(IELEM) = V1 * LGSEG(IELEM)
        W2(IELEM) = V2 * LGSEG(IELEM)
      ENDDO
C
C-----------------------------------------------------------------------
C
C     F LINEAIRE PAR TRIANGLES OU QUADRILATERES OU QUASI-BULLE
C
      ELSEIF(IELMF.EQ.11.OR.IELMF.EQ.12.OR.IELMF.EQ.21) THEN
C
      XSUR3 = XMUL/3.D0
      XSUR6 = XMUL/6.D0
C
      DO IELEM = 1,NELEM
        F1 = F(NBOR(IKLE1(IELEM)))
        F2 = F(NBOR(IKLE2(IELEM)))
        V1 = ( F1*XSUR3 + F2*XSUR6 )
        V2 = ( F2*XSUR3 + F1*XSUR6 )
        W1(IELEM) = V1 * LGSEG(IELEM)
        W2(IELEM) = V2 * LGSEG(IELEM)
      ENDDO
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
       IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
100    FORMAT(1X,'VC01OO (BIEF) :',/,
     *        1X,'DISCRETISATION DE F NON PREVUE : ',1I6,
     *        1X,'NOM REEL : ',A6)
101    FORMAT(1X,'VC01OO (BIEF) :',/,
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
