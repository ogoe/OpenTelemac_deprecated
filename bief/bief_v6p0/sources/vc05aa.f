C                       *****************
                        SUBROUTINE VC05AA
C                       *****************
C
     *( XMUL,SW,W,SURFAC,IKLE1,IKLE2,IKLE3,NELEM,NELMAX,W1,W2,W3 )
C
C***********************************************************************
C BIEF VERSION 5.1           05/02/91    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C                    /         ->
C    VEC(I) = XMUL  /    (U,V).N  PSI(I) D(GAMMA)
C                  /GAMMA
C
C
C    PSI(I) EST UNE BASE DE TYPE TRIANGLE P1
C
C    ATTENTION : LE RESULTAT EST DANS W1 SOUS FORME NON ASSEMBLEE.
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
C |      X,Y,Z     | -->|  COORDONNEES DES POINTS
C |      SURFAC    | -->|  SURFACE DES ELEMENTS.
C |      IKLE1,... | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE.
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE.
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE.
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      W1,2,3    |<-- |  VECTEUR RESULTAT SOUS FORME NON ASSEMBLEE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C  PROGRAMMES APPELES :
C
C-----------------------------------------------------------------------
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_VC05AA => VC05AA
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
      DOUBLE PRECISION, INTENT(INOUT):: W1(NELMAX),W2(NELMAX),W3(NELMAX)
      DOUBLE PRECISION, INTENT(IN)   :: SURFAC(NELMAX)
      DOUBLE PRECISION, INTENT(IN)   :: XMUL
C
C-----------------------------------------------------------------------
C
C     STRUCTURES DE W ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN)   :: SW
      DOUBLE PRECISION, INTENT(IN) :: W(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELMW,IELEM,N1,N2,N3
      DOUBLE PRECISION XSUR12,A,WT
C
C-----------------------------------------------------------------------
C
      IELMW=SW%ELM
C
C-----------------------------------------------------------------------
C
C     W LINEAIRE PAR TRIANGLES
C
      IF(IELMW.EQ.11) THEN
C
         XSUR12 = XMUL/12.D0
C
C   BOUCLE SUR LES FACES DE BORD
C
         DO 1 IELEM = 1,NELEM
C
C  NUMEROTATION LOCALE DES SOMMETS DE LA FACE
C
            N1 = IKLE1(IELEM)
            N2 = IKLE2(IELEM)
            N3 = IKLE3(IELEM)
C
            WT = W(N1) + W(N2) + W(N3)
            A = SURFAC(IELEM) * XSUR12
C
            W1(IELEM) = A * (WT + W(N1))
            W2(IELEM) = A * (WT + W(N2))
            W3(IELEM) = A * (WT + W(N3))
C
1        CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
         IF (LNG.EQ.1) WRITE(LU,100) IELMW,SW%NAME
         IF (LNG.EQ.2) WRITE(LU,101) IELMW,SW%NAME
100      FORMAT(1X,'VC05AA (BIEF) :',/,
     *          1X,'DISCRETISATION DE W NON PREVUE : ',1I6,
     *          1X,'NOM REEL : ',A6)
101      FORMAT(1X,'VC05AA (BIEF) :',/,
     *          1X,'DISCRETIZATION OF W NOT AVAILABLE:',1I6,
     *          1X,'REAL NAME: ',A6)
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
