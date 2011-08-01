C                       *****************
                        SUBROUTINE VC14AA
C                       *****************
C
     *( XMUL,SU,SV,U,V,XEL,YEL,SURFAC,IKLE1,IKLE2,IKLE3,NELEM,NELMAX,
     *  W1,W2,W3)
C
C***********************************************************************
C BIEF VERSION 5.9           09/07/2008    C. MOULIN    (LNH)
C                                        L VAN HAREN  (LNH)
C                                        A FROEHLY    (MATMECA)  
C***********************************************************************
C
C FONCTION : CALCUL DU TERME DE PRODUCTION TURBULENTE ( A NUT PRES)
C
C                   /
C  VEC(I) = XMUL   /     PSI(I) *
C                 /OMEGA
C
C
C              DU           DV          DU   DV
C     *  (  2*(--)^2  +  2*(--)^2  +  ( -- + -- )^2  )  D(OMEGA)
C              DX           DY          DY   DX
C
C
C    U ET V ETANT DES VECTEURS
C
C    PSI(I) EST UNE BASE DE TYPE IELM
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
      USE BIEF !, EX_VC14AA => VC14AA
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
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,*),YEL(NELMAX,*)
      DOUBLE PRECISION, INTENT(INOUT) ::W1(NELMAX),W2(NELMAX),W3(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: SURFAC(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
C     STRUCTURES DE U,V ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SU,SV
      DOUBLE PRECISION, INTENT(IN) :: U(*),V(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMU,IELMV
      DOUBLE PRECISION FACT,XSUR12,U21,U31,V21,V31,X2,X3,Y2,Y3
C
C-----------------------------------------------------------------------
C
      XSUR12 = XMUL / 12.D0
C
C-----------------------------------------------------------------------
C
C     ATTENTION U QUASI-BULLE ET U P2 SONT TRAITES ICI COMME LINEAIRES
C
C-----------------------------------------------------------------------
C
      IELMU=SU%ELM
      IELMV=SV%ELM
C
      IF(      (IELMU.EQ.11.AND.IELMV.EQ.11)
     *     .OR.(IELMU.EQ.12.AND.IELMV.EQ.12)  
     *     .OR.(IELMU.EQ.13.AND.IELMV.EQ.13)) THEN
C
      DO 1 IELEM = 1 , NELEM
C
        X2  = XEL(IELEM,2)
        X3  = XEL(IELEM,3)
        Y2  = YEL(IELEM,2)
        Y3  = YEL(IELEM,3)
C
        U21 = U(IKLE2(IELEM)) - U(IKLE1(IELEM))
        U31 = U(IKLE3(IELEM)) - U(IKLE1(IELEM))
        V21 = V(IKLE2(IELEM)) - V(IKLE1(IELEM))
        V31 = V(IKLE3(IELEM)) - V(IKLE1(IELEM))
C
        FACT = (       ( X2*U31-X3*U21-Y2*V31+Y3*V21)**2
     *           + 2 *(( X2*V31-X3*V21 )**2+( Y3*U21-Y2*U31 )**2 )
     *         )  * XSUR12 / SURFAC(IELEM)
C
        W1(IELEM) = FACT
        W2(IELEM) = FACT
        W3(IELEM) = FACT
C
1     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
       IF (LNG.EQ.1) WRITE(LU,100) IELMU,SU%NAME
       IF (LNG.EQ.1) WRITE(LU,200) IELMV,SV%NAME
       IF (LNG.EQ.1) WRITE(LU,300)
       IF (LNG.EQ.2) WRITE(LU,101) IELMU,SU%NAME
       IF (LNG.EQ.2) WRITE(LU,201) IELMV,SV%NAME
       IF (LNG.EQ.2) WRITE(LU,301)
100    FORMAT(1X,'VC14AA (BIEF) :',/,
     *        1X,'DISCRETISATION DE U : ',1I6,
     *        1X,'NOM REEL : ',A6)
200    FORMAT(1X,'DISCRETISATION DE V : ',1I6,
     *        1X,'NOM REEL : ',A6)
300    FORMAT(1X,'CAS NON PREVU')
101    FORMAT(1X,'VC14AA (BIEF) :',/,
     *        1X,'DISCRETIZATION OF U:',1I6,
     *        1X,'REAL NAME: ',A6)
201    FORMAT(1X,'DISCRETIZATION OF V:',1I6,
     *        1X,'REAL NAME: ',A6)
301    FORMAT(1X,'CASE NOT IMPLEMENTED')
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
