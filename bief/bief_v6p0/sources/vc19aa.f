C                       *****************
                        SUBROUTINE VC19AA
C                       *****************
C
     *( XMUL,SF,SG,SH,SU,SV,F,G,H,U,V,
     *  XEL,YEL,SURFAC,IKLE1,IKLE2,IKLE3,NELEM,NELMAX,W1,W2,W3,FORMUL)
C
C***********************************************************************
C BIEF VERSION 5.8           09/01/08    C-T PHAM (LNHE) 01 30 87 85 93
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C
C                    /          DPSII      DPSII
C      V  =  XMUL   /   F * ( U -----  + V ----- )   D(OMEGA)
C       I          /OMEGA        DX         DY
C
C
C    PSI(I) EST UNE BASE DE TYPE TRIANGLE P1
C
C    F, U ET V SONT DES VECTEURS
C
C
C    ATTENTION : LE RESULTAT EST DANS W SOUS FORME NON ASSEMBLEE.
C
C                                                     ->      
C    ATTENTION : SI FORMUL='HUGRADP ' LA VITESSE EST: U   
C                                                     ->      --->
C                SI FORMUL='HUGRADP2' LA VITESSE EST: U + G * GRAD(H) 
C
C
C
C    ATTENTION : U ET V QUASI-BULLE TRAITES COMME LINEAIRES !!!!!!!!!!!!
C
C
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      XMUL      | -->|  COEFFICIENT MULTIPLICATEUR.
C |      SF        | -->|  STRUCTURE DE LA FONCTION H
C |      SU,SV     | -->|  STRUCTURES DES FONCTIONS U ET V
C |      H         | -->|  FONCTIONS INTERVENANT DANS LA FORMULE.
C |      U,V       | -->|  COMPOSANTES D'UN VECTEUR
C |                |    |  INTERVENANT DANS LA FORMULE.
C |      XEL,YEL,..| -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |      SURFAC    | -->|  SURFACE DES ELEMENTS.
C |      IKLE1,... | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE.
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE.
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE.
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      W1,W2,W3  |<-- |  VECTEUR RESULTAT SOUS FORME NON ASSEMBLEE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C  PAS DE PROGRAMMES APPELES
C
C-----------------------------------------------------------------------
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF !, EX_VC19AA => VC19AA
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX),IKLE3(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL
      DOUBLE PRECISION, INTENT(IN)    :: XEL(NELMAX,*),YEL(NELMAX,*)
      DOUBLE PRECISION, INTENT(IN)    :: SURFAC(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W1(*),W2(*),W3(*)
C
C     STRUCTURES DE H,U,V ET DONNEES REELLES
C
      TYPE(BIEF_OBJ), INTENT(IN)    :: SF,SH,SG,SU,SV
      DOUBLE PRECISION, INTENT(IN)  :: F(*),G(*),H(*),U(*),V(*)
      CHARACTER(LEN=16), INTENT(IN) :: FORMUL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMU,IELMV,IELMF,IELMH,IELMG,I1,I2,I3
C
      DOUBLE PRECISION X2,Y2,X3,Y3,Z1,Z2,Z3,ZX,ZY,DET
      DOUBLE PRECISION H1,H2,H3,U1,U2,U3,V1,V2,V3
      DOUBLE PRECISION H123,U123,V123,HU123,HV123
      DOUBLE PRECISION XSUR24
C
C-----------------------------------------------------------------------
C
      XSUR24 = XMUL / 24.D0
C
C-----------------------------------------------------------------------
C
      IELMU=SU%ELM
      IELMV=SV%ELM
      IELMF=SF%ELM
      IELMG=SG%ELM
      IELMH=SH%ELM
C
C-----------------------------------------------------------------------
C
      IF(FORMUL(1:8).EQ.'HUGRADP ') THEN
C
C     F, U, V LINEAIRES
C
      IF(        IELMF.EQ.11
     *     .AND.(IELMU.EQ.11.OR.IELMU.EQ.12)
     *     .AND.(IELMV.EQ.11.OR.IELMV.EQ.12)  ) THEN
C
        DO IELEM = 1 , NELEM
C
          X2 = XEL(IELEM,2)
          X3 = XEL(IELEM,3)
          Y2 = YEL(IELEM,2)
          Y3 = YEL(IELEM,3)
C
          I1 = IKLE1(IELEM)
          I2 = IKLE2(IELEM)
          I3 = IKLE3(IELEM)
          H1 = F(I1)
          H2 = F(I2)
          H3 = F(I3)
          U1 = U(I1)
          U2 = U(I2)
          U3 = U(I3)
          V1 = V(I1)
          V2 = V(I2)
          V3 = V(I3)
C
          H123 = H1+H2+H3
          U123 = U1+U2+U3
          V123 = V1+V2+V3
C
          HU123 = H1*U1+H2*U2+H3*U3
          HV123 = H1*V1+H2*V2+H3*V3
C
          W1(IELEM) = ( (Y2-Y3)*(H123*U123+HU123)
     #                 +(X3-X2)*(H123*V123+HV123) )*XSUR24
          W2(IELEM) = (     Y3 *(H123*U123+HU123)
     #                     -X3 *(H123*V123+HV123) )*XSUR24
          W3(IELEM) = (    -Y2 *(H123*U123+HU123)
     #                     +X2 *(H123*V123+HV123) )*XSUR24
C
       ENDDO
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
        IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
        IF (LNG.EQ.1) WRITE(LU,400) IELMU,SU%NAME
        IF (LNG.EQ.1) WRITE(LU,500) IELMV,SV%NAME
        IF (LNG.EQ.1) WRITE(LU,600)
        IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
        IF (LNG.EQ.2) WRITE(LU,401) IELMU,SU%NAME
        IF (LNG.EQ.2) WRITE(LU,501) IELMV,SV%NAME
        IF (LNG.EQ.2) WRITE(LU,601)
100     FORMAT(1X,'VC19AA (BIEF) :',/,
     *         1X,'DISCRETISATION DE F : ',1I6,
     *         1X,'NOM REEL : ',A6)
200     FORMAT(1X,'DISCRETISATION DE G : ',1I6,
     *         1X,'NOM REEL : ',A6)
300     FORMAT(1X,'DISCRETISATION DE H : ',1I6,
     *         1X,'NOM REEL : ',A6)
400     FORMAT(1X,'DISCRETISATION DE U : ',1I6,
     *         1X,'NOM REEL : ',A6)
500     FORMAT(1X,'DISCRETISATION DE V : ',1I6,
     *         1X,'NOM REEL : ',A6)
600     FORMAT(1X,'CAS NON PREVU')
101     FORMAT(1X,'VC19AA (BIEF) :',/,
     *         1X,'DISCRETIZATION OF F:',1I6,
     *         1X,'REAL NAME: ',A6)
201     FORMAT(1X,'DISCRETIZATION OF G:',1I6,
     *         1X,'REAL NAME: ',A6)
301     FORMAT(1X,'DISCRETIZATION OF H:',1I6,
     *         1X,'REAL NAME: ',A6)
401     FORMAT(1X,'DISCRETIZATION OF U:',1I6,
     *         1X,'REAL NAME: ',A6)
501     FORMAT(1X,'DISCRETIZATION OF V:',1I6,
     *         1X,'REAL NAME: ',A6)
601     FORMAT(1X,'CASE NOT IMPLEMENTED')
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(FORMUL(1:8).EQ.'HUGRADP2') THEN
C
C-----------------------------------------------------------------------
C
C
C     F, G, U, V LINEAIRES   H LINEAIRE PAR MORCEAU
C
      IF(       IELMF.EQ.11
     *     .AND.IELMG.EQ.11
     *     .AND.IELMH.EQ.10
     *     .AND.(IELMU.EQ.11.OR.IELMU.EQ.12)
     *     .AND.(IELMV.EQ.11.OR.IELMV.EQ.12)  ) THEN
C
        DO IELEM = 1 , NELEM
C           
          X2 = XEL(IELEM,2)
          X3 = XEL(IELEM,3)
          Y2 = YEL(IELEM,2)
          Y3 = YEL(IELEM,3)
C
          I1 = IKLE1(IELEM)
          I2 = IKLE2(IELEM)
          I3 = IKLE3(IELEM)
          H1 = F(I1)
          H2 = F(I2)
          H3 = F(I3)
          DET=X2*Y3-X3*Y2
          Z1=H(IELEM)
          Z2=H(IELEM+  NELMAX)-Z1
          Z3=H(IELEM+2*NELMAX)-Z1
          ZX=(Z2*Y3-Z3*Y2)/DET
          ZY=(X2*Z3-X3*Z2)/DET
          U1 = U(I1) + G(I1)*ZX
          U2 = U(I2) + G(I2)*ZX
          U3 = U(I3) + G(I3)*ZX
          V1 = V(I1) + G(I1)*ZY
          V2 = V(I2) + G(I2)*ZY
          V3 = V(I3) + G(I3)*ZY
C
          H123 = H1+H2+H3
          U123 = U1+U2+U3
          V123 = V1+V2+V3
C
          HU123 = H1*U1+H2*U2+H3*U3
          HV123 = H1*V1+H2*V2+H3*V3
C
          W1(IELEM) = ( (Y2-Y3)*(H123*U123+HU123)
     #                 +(X3-X2)*(H123*V123+HV123) )*XSUR24
          W2(IELEM) = (     Y3 *(H123*U123+HU123)
     #                     -X3 *(H123*V123+HV123) )*XSUR24
          W3(IELEM) = (    -Y2 *(H123*U123+HU123)
     #                     +X2 *(H123*V123+HV123) )*XSUR24
C
       ENDDO
C
C     F, G, H, U, V LINEAIRES
C
      ELSEIF(       IELMF.EQ.11
     *         .AND.IELMG.EQ.11
     *         .AND.IELMH.EQ.11
     *         .AND.(IELMU.EQ.11.OR.IELMU.EQ.12)
     *         .AND.(IELMV.EQ.11.OR.IELMV.EQ.12)  ) THEN
C
        DO IELEM = 1 , NELEM
C           
          X2 = XEL(IELEM,2)
          X3 = XEL(IELEM,3)
          Y2 = YEL(IELEM,2)
          Y3 = YEL(IELEM,3)
C
          I1 = IKLE1(IELEM)
          I2 = IKLE2(IELEM)
          I3 = IKLE3(IELEM)
          H1 = F(I1)
          H2 = F(I2)
          H3 = F(I3)
          DET=X2*Y3-X3*Y2
          Z1=H(I1)
          Z2=H(I2)-Z1
          Z3=H(I3)-Z1
          ZX=(Z2*Y3-Z3*Y2)/DET
          ZY=(X2*Z3-X3*Z2)/DET
          U1 = U(I1) + G(I1)*ZX
          U2 = U(I2) + G(I2)*ZX
          U3 = U(I3) + G(I3)*ZX
          V1 = V(I1) + G(I1)*ZY
          V2 = V(I2) + G(I2)*ZY
          V3 = V(I3) + G(I3)*ZY
C
          H123 = H1+H2+H3
          U123 = U1+U2+U3
          V123 = V1+V2+V3
C
          HU123 = H1*U1+H2*U2+H3*U3
          HV123 = H1*V1+H2*V2+H3*V3
C
          W1(IELEM) = ( (Y2-Y3)*(H123*U123+HU123)
     #                 +(X3-X2)*(H123*V123+HV123) )*XSUR24
          W2(IELEM) = (     Y3 *(H123*U123+HU123)
     #                     -X3 *(H123*V123+HV123) )*XSUR24
          W3(IELEM) = (    -Y2 *(H123*U123+HU123)
     #                     +X2 *(H123*V123+HV123) )*XSUR24
C
       ENDDO
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
        IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
        IF (LNG.EQ.1) WRITE(LU,200) IELMG,SG%NAME
        IF (LNG.EQ.1) WRITE(LU,300) IELMH,SH%NAME
        IF (LNG.EQ.1) WRITE(LU,400) IELMU,SU%NAME
        IF (LNG.EQ.1) WRITE(LU,500) IELMV,SV%NAME
        IF (LNG.EQ.1) WRITE(LU,600)
        IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
        IF (LNG.EQ.2) WRITE(LU,201) IELMG,SG%NAME
        IF (LNG.EQ.2) WRITE(LU,301) IELMH,SH%NAME
        IF (LNG.EQ.2) WRITE(LU,401) IELMU,SU%NAME
        IF (LNG.EQ.2) WRITE(LU,401) IELMV,SV%NAME
        IF (LNG.EQ.2) WRITE(LU,601)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(FORMUL(1:8).EQ.'HUGRADP3') THEN
C
C-----------------------------------------------------------------------
C
C
C     F, G,   H LINEAIRE PAR MORCEAU    U ET V PAS PRIS EN COMPTE
C            OU LINEAIRE
C
      IF(       IELMF.EQ.11
     *     .AND.IELMG.EQ.11
     *     .AND.IELMH.EQ.10  ) THEN
C
        DO IELEM = 1 , NELEM
C           
          X2 = XEL(IELEM,2)
          X3 = XEL(IELEM,3)
          Y2 = YEL(IELEM,2)
          Y3 = YEL(IELEM,3)
C
          I1 = IKLE1(IELEM)
          I2 = IKLE2(IELEM)
          I3 = IKLE3(IELEM)
          H1 = F(I1)
          H2 = F(I2)
          H3 = F(I3)
          DET=X2*Y3-X3*Y2
          Z1=H(IELEM)
          Z2=H(IELEM+  NELMAX)-Z1
          Z3=H(IELEM+2*NELMAX)-Z1
          ZX=(Z2*Y3-Z3*Y2)/DET
          ZY=(X2*Z3-X3*Z2)/DET
          U1 = G(I1)*ZX
          U2 = G(I2)*ZX
          U3 = G(I3)*ZX
          V1 = G(I1)*ZY
          V2 = G(I2)*ZY
          V3 = G(I3)*ZY
C
          H123 = H1+H2+H3
          U123 = U1+U2+U3
          V123 = V1+V2+V3
C
          HU123 = H1*U1+H2*U2+H3*U3
          HV123 = H1*V1+H2*V2+H3*V3
C
          W1(IELEM) = ( (Y2-Y3)*(H123*U123+HU123)
     #                 +(X3-X2)*(H123*V123+HV123) )*XSUR24
          W2(IELEM) = (     Y3 *(H123*U123+HU123)
     #                     -X3 *(H123*V123+HV123) )*XSUR24
          W3(IELEM) = (    -Y2 *(H123*U123+HU123)
     #                     +X2 *(H123*V123+HV123) )*XSUR24
C
       ENDDO
C
       ELSEIF(       IELMF.EQ.11
     *          .AND.IELMG.EQ.11
     *          .AND.IELMH.EQ.11  ) THEN
C
        DO IELEM = 1 , NELEM
C           
          X2 = XEL(IELEM,2)
          X3 = XEL(IELEM,3)
          Y2 = YEL(IELEM,2)
          Y3 = YEL(IELEM,3)
C
          I1 = IKLE1(IELEM)
          I2 = IKLE2(IELEM)
          I3 = IKLE3(IELEM)
          H1 = F(I1)
          H2 = F(I2)
          H3 = F(I3)
          DET=X2*Y3-X3*Y2
          Z1=H(I1)
          Z2=H(I2)-Z1
          Z3=H(I3)-Z1
          ZX=(Z2*Y3-Z3*Y2)/DET
          ZY=(X2*Z3-X3*Z2)/DET
          U1 = G(I1)*ZX
          U2 = G(I2)*ZX
          U3 = G(I3)*ZX
          V1 = G(I1)*ZY
          V2 = G(I2)*ZY
          V3 = G(I3)*ZY
C
          H123 = H1+H2+H3
          U123 = U1+U2+U3
          V123 = V1+V2+V3
C
          HU123 = H1*U1+H2*U2+H3*U3
          HV123 = H1*V1+H2*V2+H3*V3
C
          W1(IELEM) = ( (Y2-Y3)*(H123*U123+HU123)
     #                 +(X3-X2)*(H123*V123+HV123) )*XSUR24
          W2(IELEM) = (     Y3 *(H123*U123+HU123)
     #                     -X3 *(H123*V123+HV123) )*XSUR24
          W3(IELEM) = (    -Y2 *(H123*U123+HU123)
     #                     +X2 *(H123*V123+HV123) )*XSUR24
C
       ENDDO
C
       ELSE
        IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
        IF (LNG.EQ.1) WRITE(LU,200) IELMG,SG%NAME
        IF (LNG.EQ.1) WRITE(LU,300) IELMH,SH%NAME
        IF (LNG.EQ.1) WRITE(LU,600)
        IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
        IF (LNG.EQ.2) WRITE(LU,201) IELMG,SG%NAME
        IF (LNG.EQ.2) WRITE(LU,301) IELMH,SH%NAME
        IF (LNG.EQ.2) WRITE(LU,601)
        CALL PLANTE(1)
        STOP       
       ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
        IF (LNG.EQ.2) WRITE(LU,2000) FORMUL
1000    FORMAT(1X,'VC19AA (BIEF) :',/,
     *         1X,'FORMULE : ',A16,' NON PREVUE')
2000    FORMAT(1X,'VC19AA (BIEF):',/,
     *         1X,'FORMULA: ',A16,' UNEXPECTED')
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
