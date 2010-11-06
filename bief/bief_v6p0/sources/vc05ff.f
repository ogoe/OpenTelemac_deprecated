C                       *****************
                        SUBROUTINE VC05FF
C                       *****************
C
     *( XMUL,SU,SV,U,V,X,Y,Z,
     *  IKLE1,IKLE2,IKLE3,IKLE4,NBOR,NELEM,NELMAX,W1,W2,W3,W4)
C
C***********************************************************************
C BIEF VERSION 6.0      24/07/2009    J-M HERVOUET (LNHE) 01 30 87 80 18
C                                        
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C                    /         ->
C    VEC(I) = XMUL  /    (U,V).N  PSI(I)  D(GAMMA)
C                  /GAMMA
C
C
C    PSI(I) EST UNE BASE DE TYPE QUADRILATERE P1
C
C    ATTENTION : LE RESULTAT EST DANS W SOUS FORME NON ASSEMBLEE.
C
C                MAILLAGE REEL
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
C |      X,Y,Z     | -->|  COORDONNEES DES POINTS
C |      SURFAC    | -->|  SURFACE DES ELEMENTS.
C |      IKLE1,... | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE.
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE.
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE.
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      W1,2,3,4  |<-- |  VECTEUR RESULTAT SOUS FORME NON ASSEMBLEE.
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
      USE BIEF, EX_VC05FF => VC05FF
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
      DOUBLE PRECISION, INTENT(IN)    :: XMUL
C
C     STRUCTURES DE U,V ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN)   :: SU,SV
      DOUBLE PRECISION, INTENT(IN) :: U(*),V(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELMU,IELMV,IELEM,N1,N2,N3,N4,I1,I2,I3,I4
C
      DOUBLE PRECISION XSUR72,H1,H2,HT,AX,AY
      DOUBLE PRECISION U1,U2,U3,U4,V1,V2,V3,V4
C
C-----------------------------------------------------------------------
C
      IELMU=SU%ELM
      IELMV=SV%ELM
C
C-----------------------------------------------------------------------
C
      XSUR72 = XMUL/72.D0
C
C     U LINEAIRE PAR PRISMES
C
C-----------------------------------------------------------------------
C
      IF(IELMU.EQ.71.AND.IELMV.EQ.71) THEN
C
C-----------------------------------------------------------------------
C
C        BOUCLE SUR LES FACES DE BORD
C
         DO IELEM = 1,NELEM
C
C           NUMEROTATION LOCALE DES SOMMETS DE LA FACE
C
            I1 = IKLE1(IELEM)
            I2 = IKLE2(IELEM)
            I3 = IKLE3(IELEM)
            I4 = IKLE4(IELEM)
C
C           NUMEROTATION GLOBALE DES SOMMETS DE LA FACE
C
            N1 = NBOR(I1)
            N2 = NBOR(I2)
            N3 = NBOR(I3)
            N4 = NBOR(I4)
C
            H1 = Z(N4) - Z(N1)
            H2 = Z(N3) - Z(N2)
            HT = H1 + H2
            H1 = H1 + H1 + HT
            H2 = H2 + H2 + HT
C
            U1 = U(I1) + U(I1) + U(I4)
            U2 = U(I2) + U(I2) + U(I3)
            U3 = U(I2) + U(I3) + U(I3)
            U4 = U(I1) + U(I4) + U(I4)
C
            AX = (Y(N2)-Y(N1)) * XSUR72
C
            V1 = V(I1) + V(I1) + V(I4)
            V2 = V(I2) + V(I2) + V(I3)
            V3 = V(I2) + V(I3) + V(I3)
            V4 = V(I1) + V(I4) + V(I4)
C
            AY = (X(N1)-X(N2)) * XSUR72
C
            W1(IELEM) = (U1*H1+U2*HT)*AX + (V1*H1+V2*HT)*AY
            W2(IELEM) = (U1*HT+U2*H2)*AX + (V1*HT+V2*H2)*AY
            W3(IELEM) = (U4*HT+U3*H2)*AX + (V4*HT+V3*H2)*AY
            W4(IELEM) = (U4*H1+U3*HT)*AX + (V4*H1+V3*HT)*AY
C
         ENDDO
C
C-----------------------------------------------------------------------
C

C
C-----------------------------------------------------------------------
C
      ELSEIF(IELMU.EQ.41.AND.IELMV.EQ.41) THEN
C
C-----------------------------------------------------------------------
C
C        BOUCLE SUR LES FACES DE BORD
C
         DO IELEM = 1,NELEM
C
C  NUMEROTATION GLOBALE DES SOMMETS DE LA FACE
C
            N1 = NBOR(IKLE1(IELEM))
            N2 = NBOR(IKLE2(IELEM))
            N3 = NBOR(IKLE3(IELEM))
            N4 = NBOR(IKLE4(IELEM))
C
            H1 = Z(N4) - Z(N1)
            H2 = Z(N3) - Z(N2)
            HT = H1 + H2
            H1 = H1 + H1 + HT
            H2 = H2 + H2 + HT
C
            U1 = U(N1) + U(N1) + U(N4)
            U2 = U(N2) + U(N2) + U(N3)
            U3 = U(N2) + U(N3) + U(N3)
            U4 = U(N1) + U(N4) + U(N4)
            AX = (Y(N2)-Y(N1)) * XSUR72
C
            V1 = V(N1) + V(N1) + V(N4)
            V2 = V(N2) + V(N2) + V(N3)
            V3 = V(N2) + V(N3) + V(N3)
            V4 = V(N1) + V(N4) + V(N4)
            AY = (X(N1)-X(N2)) * XSUR72
C
            W1(IELEM) = (U1*H1+U2*HT)*AX + (V1*H1+V2*HT)*AY
            W2(IELEM) = (U1*HT+U2*H2)*AX + (V1*HT+V2*H2)*AY
            W3(IELEM) = (U4*HT+U3*H2)*AX + (V4*HT+V3*H2)*AY
            W4(IELEM) = (U4*H1+U3*HT)*AX + (V4*H1+V3*HT)*AY
C
         ENDDO
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
         IF (LNG.EQ.1) WRITE(LU,100) IELMU,SU%NAME
         IF (LNG.EQ.2) WRITE(LU,101) IELMU,SU%NAME
100      FORMAT(1X,'VC05FF (BIEF) :',/,
     *          1X,'DISCRETISATION DE U NON PREVUE : ',1I6,
     *          1X,'NOM REEL : ',A6)
101      FORMAT(1X,'VC05FF (BIEF) :',/,
     *          1X,'DISCRETIZATION OF U NOT AVAILABLE:',1I6,
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
