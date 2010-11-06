C                       *****************
                        SUBROUTINE VC05OO
C                       *****************
C
     *(XMUL,SU,SV,U,V,XNOR,YNOR,LGSEG,IKLE,NBOR,NELEM,NELMAX,W1,W2 )
C
C***********************************************************************
C BIEF VERSION 5.9        29/05/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C                    /             ->   ->
C    VEC(I) = XMUL  /    PSI(I) *  U  . N  D(OMEGA)
C                  /OMEGA
C
C
C    PSI(I) EST UNE BASE DE TYPE SEGMENT P1
C
C    ->
C    U EST UN VECTEUR DE COMPOSANTES U ET V
C
C    ->
C    N EST LE VECTEUR NORMAL EXTERIEUR A L'ELEMENT.
C
C    ATTENTION : LE RESULTAT EST DANS W SOUS FORME NON ASSEMBLEE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      XMUL      | -->|  COEFFICIENT MULTIPLICATEUR.
C |      SU,SV,SW  | -->|  STRUCTURES DES FONCTIONS U,V ET W
C |      U,V,W     | -->|  COMPOSANTES D'UN VECTEUR
C |                |    |  INTERVENANT DANS LA FORMULE.
C |      XEL,YEL,..| -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |      XNOR,..   | -->|  NORMALES AUX SEGMENTS.
C |      LGSEG..   | -->|  LONGUEUR DES SEGMENTS.
C |      IKLE1,... | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE.
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE.
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE.
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      W1,2      | -->|  VECTEUR RESULTAT SOUS FORME NON ASSEMBLEE.
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
      USE BIEF, EX_VC05OO => VC05OO
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE(NELMAX,*)
      INTEGER, INTENT(IN) :: NBOR(*)
C
      DOUBLE PRECISION, INTENT(IN)    :: XNOR(NELMAX),YNOR(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W1(NELMAX),W2(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: LGSEG(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL
C
C     STRUCTURES DE U,V ET DONNEES REELLES
C
      TYPE(BIEF_OBJ)  , INTENT(IN) :: SU,SV
      DOUBLE PRECISION, INTENT(IN) :: U(*),V(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER N1,N2,NG1,NG2,IELEM,IELMU,IELMV
      DOUBLE PRECISION XSUR06,U1,U2,V1,V2,VX1,VY1,VX2,VY2
C
C-----------------------------------------------------------------------
C
      XSUR06 = XMUL/6.D0
C
C-----------------------------------------------------------------------
C
      IELMU=SU%ELM
      IELMV=SV%ELM
C
C-----------------------------------------------------------------------
C        ->
C   F ET U  FONCTIONS LINEAIRES SUR DES TRIANGLES OU DES QUADRILATERES
C
      IF( (IELMU.EQ.11.OR.IELMU.EQ.12.OR.IELMU.EQ.21) .AND.
     *    (IELMV.EQ.11.OR.IELMV.EQ.12.OR.IELMV.EQ.21)       ) THEN
C
      DO IELEM =1,NELEM
C
C     NUMEROTATION DES POINTS DE BORD
C
C     NUMEROTATION GLOBALE
C
      NG1= NBOR(IKLE(IELEM,1))
      NG2= NBOR(IKLE(IELEM,2))
C
      U1 = U(NG1)
      U2 = U(NG2)
      V1 = V(NG1)
      V2 = V(NG2)
C
C   DETERMINATION DES FONCTIONS DE BASES SUR LA FRONTIERE :
C
      VX1 = XSUR06 * ( U2 + U1 + U1 )
      VY1 = XSUR06 * ( V2 + V1 + V1 )
      VX2 = XSUR06 * ( U1 + U2 + U2 )
      VY2 = XSUR06 * ( V1 + V2 + V2 )
C
      W1(IELEM) = LGSEG(IELEM) * ( VX1*XNOR(IELEM) + VY1*YNOR(IELEM) )
      W2(IELEM) = LGSEG(IELEM) * ( VX2*XNOR(IELEM) + VY2*YNOR(IELEM) )
C
      ENDDO
C
C-----------------------------------------------------------------------
C   ->
C   U  FONCTIONS LINEAIRES SUR DES SEGMENTS
C
      ELSEIF(IELMU.EQ.1.AND.IELMV.EQ.1) THEN
C
      DO IELEM =1,NELEM
C
C     NUMEROTATION DES POINTS DE BORD
C
      N1 = IKLE(IELEM,1)
      N2 = IKLE(IELEM,2)
C
C     NUMEROTATION GLOBALE
C
      NG1= NBOR(N1)
      NG2= NBOR(N2)
C
      U1 = U(N1)
      U2 = U(N2)
      V1 = V(N1)
      V2 = V(N2)
C
C   DETERMINATION DES FONCTIONS DE BASES SUR LA FRONTIERE :
C
      VX1 = XSUR06 * ( U2 + U1 + U1 )
      VY1 = XSUR06 * ( V2 + V1 + V1 )
      VX2 = XSUR06 * ( U1 + U2 + U2 )
      VY2 = XSUR06 * ( V1 + V2 + V2 )
C
      W1(IELEM) = LGSEG(IELEM) * ( VX1*XNOR(IELEM) + VY1*YNOR(IELEM) )
      W2(IELEM) = LGSEG(IELEM) * ( VX2*XNOR(IELEM) + VY2*YNOR(IELEM) )
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
        IF (LNG.EQ.1) WRITE(LU,100)
        IF (LNG.EQ.1) WRITE(LU,102) IELMU,SU%NAME
        IF (LNG.EQ.1) WRITE(LU,103) IELMV,SV%NAME
        IF (LNG.EQ.1) WRITE(LU,104)
        IF (LNG.EQ.2) WRITE(LU,110)
        IF (LNG.EQ.2) WRITE(LU,112) IELMU,SU%NAME
        IF (LNG.EQ.2) WRITE(LU,113) IELMV,SV%NAME
        IF (LNG.EQ.2) WRITE(LU,114)
100     FORMAT(1X,'VC05OO (BIEF) :')
102     FORMAT(1X,'DISCRETISATION DE U : ',1I6,
     *         1X,'NOM REEL : ',A6)
103     FORMAT(1X,'DISCRETISATION DE V : ',1I6,
     *         1X,'NOM REEL : ',A6)
104     FORMAT(1X,'CAS NON PREVU')
110     FORMAT(1X,'VC05OO (BIEF):')
112     FORMAT(1X,'DISCRETIZATION OF U:',1I6,
     *         1X,'REAL NAME: ',A6)
113     FORMAT(1X,'DISCRETIZATION OF V:',1I6,
     *         1X,'REAL NAME: ',A6)
114     FORMAT(1X,'CASE NOT IMPLEMENTED')
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
