C                       *****************
                        SUBROUTINE CHAR13
C                       *****************
C
     * ( U , V , DT , NRK , X , Y , IKLE , IFABOR ,
     *   XPLOT , YPLOT , DX , DY , SHP , ELT , NSP ,
     *   NPLOT , NPOIN , NELEM , NELMAX , SURDET , SENS , TEST)
C
C***********************************************************************
C BIEF VERSION 5.9        08/08/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C  FONCTION :
C
C     REMONTEE OU DESCENTE
C     DES COURBES CARACTERISTIQUES
C     SUR DES QUADRILATERES P1
C     DANS L'INTERVALLE DE TEMPS DT
C     AVEC UNE DISCRETISATION ELEMENTS FINIS
C
C
C  DISCRETISATION :
C
C     LE DOMAINE EST APPROCHE PAR UNE DISCRETISATION ELEMENTS FINIS
C     UNE APPROXIMATION LOCALE EST DEFINIE POUR LE VECTEUR VITESSE :
C     LA VALEUR EN UN POINT D'UN ELEMENT NE DEPEND QUE DES VALEURS
C     AUX NOEUDS DE CET ELEMENT
C
C
C  RESTRICTIONS ET HYPOTHESES :
C
C     LE CHAMP CONVECTEUR U EST SUPPOSE INDEPENDANT DU TEMPS
C     LE DERIVANT EST SUPPOSE PONCTUEL DONC NON DISPERSIF
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    U,V         | -->| COMPOSANTE DE LA VITESSE DU CONVECTEUR       |
C |    DT          | -->| PAS DE TEMPS.                                |
C |    NRK         | -->| NOMBRE DE SOUS-PAS DE RUNGE-KUTTA.           |
C |    X,Y         | -->| COORDONNEES DES POINTS DU MAILLAGE.          |
C |    IKLE        | -->| TRANSITION ENTRE LES NUMEROTATIONS LOCALE    |
C |                |    | ET GLOBALE.                                  |
C |    IFABOR      | -->| NUMEROS DES ELEMENTS AYANT UNE FACE COMMUNE  |
C |                |    | AVEC L'ELEMENT .  SI IFABOR<0 OU NUL         |
C |                |    | ON A UNE FACE LIQUIDE,SOLIDE,OU PERIODIQUE   |
C |  XPLOT,YPLOT   |<-->| POSITIONS SUCCESSIVES DES DERIVANTS.         |
C |    DX,DY       | -- | STOCKAGE DES SOUS-PAS . |
C |    SHP         |<-->| COORDONNEES BARYCENTRIQUES 2D AU PIED DES    |
C |                |    | COURBES CARACTERISTIQUES.                    |
C |    ELT         |<-->| NUMEROS DES ELEMENTS 2D AU PIED DES COURBES  |
C |                |    | CARACTERISTIQUES.                            |
C |    NSP         | -- | NOMBRE DE SOUS-PAS DE RUNGE KUTTA.           |
C |    NPLOT       | -->| NOMBRE DE DERIVANTS.                         |
C |    NPOIN       | -->| NOMBRE DE POINTS DU MAILLAGE.                |
C |    NELEM       | -->| NOMBRE D'ELEMENTS.                           |
C |    NELMAX      | -->| NOMBRE MAXIMAL D'ELEMENTS DANS LE MAILLAGE 2D|
C |    SURDET      | -->| VARIABLE UTILISEE PAR LA TRANSFORMEE ISOPARAM.
C |    SENS        | -->| DESCENTE OU REMONTEE DES CARACTERISTIQUES.   |
C |________________|____|______________________________________________|
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C     - APPELE PAR : CARACT , DERIVE , DERLAG
C     - PROGRAMMES APPELES : NEANT
C
C***********************************************************************
C
      USE BIEF, EX_CHAR13 => CHAR13
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)    :: SENS
      INTEGER         , INTENT(IN)    :: NPOIN,NELEM,NELMAX,NPLOT,NRK
      INTEGER         , INTENT(IN)    :: IKLE(NELMAX,6),IFABOR(NELMAX,3)
      INTEGER         , INTENT(INOUT) :: ELT(NPLOT),NSP(NPLOT)
      DOUBLE PRECISION, INTENT(IN)    :: U(NPOIN),V(NPOIN),SURDET(NELEM)
      DOUBLE PRECISION, INTENT(INOUT) :: XPLOT(NPLOT),YPLOT(NPLOT)
      DOUBLE PRECISION, INTENT(INOUT) :: SHP(3,NPLOT)
      DOUBLE PRECISION, INTENT(IN)    :: DT
      DOUBLE PRECISION, INTENT(IN)    :: X(*),Y(*)
      DOUBLE PRECISION, INTENT(INOUT) :: DX(NPLOT),DY(NPLOT),TEST(NPLOT)
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IPLOT,NSPMAX,ISP,I1,I2,I3,I4,I5,I6,IEL,ISO,IFA,ISUI(3)
C
      DOUBLE PRECISION PAS,EPSILO,A1,DX1,DY1,DXP,DYP,XP,YP,DENOM
C
      DATA ISUI   / 2 , 3 , 1 /
      DATA EPSILO / -1.D-6 /
C
      INTRINSIC INT,MAX,MIN,SQRT
C
C***********************************************************************
C  CALCUL DU NOMBRE DE SOUS PAS, LE MEME POUR TOUS LES POINTS
C-----------------------------------------------------------------------
C
      NSPMAX = 1
C
      DO 10 IPLOT = 1 ,NPLOT
C
         NSP(IPLOT) = 0
         IEL = ELT(IPLOT)
C
         IF(IEL.GT.0) THEN
C
            I1 = IKLE(IEL,1)
            I2 = IKLE(IEL,2)
            I3 = IKLE(IEL,3)
            I4 = IKLE(IEL,4)
            I5 = IKLE(IEL,5)
            I6 = IKLE(IEL,6)
C
            DXP = U(I1)*(2.D0*SHP(1,IPLOT)-1.D0)*SHP(1,IPLOT)
     *          + U(I2)*(2.D0*SHP(2,IPLOT)-1.D0)*SHP(2,IPLOT)
     *          + U(I3)*(2.D0*SHP(3,IPLOT)-1.D0)*SHP(3,IPLOT)
     *          + U(I4)*4.D0*SHP(1,IPLOT)*SHP(2,IPLOT)
     *          + U(I5)*4.D0*SHP(2,IPLOT)*SHP(3,IPLOT)
     *          + U(I6)*4.D0*SHP(3,IPLOT)*SHP(1,IPLOT)

            DYP = V(I1)*(2.D0*SHP(1,IPLOT)-1.D0)*SHP(1,IPLOT)
     *          + V(I2)*(2.D0*SHP(2,IPLOT)-1.D0)*SHP(2,IPLOT)
     *          + V(I3)*(2.D0*SHP(3,IPLOT)-1.D0)*SHP(3,IPLOT)
     *          + V(I4)*4.D0*SHP(1,IPLOT)*SHP(2,IPLOT)
     *          + V(I5)*4.D0*SHP(2,IPLOT)*SHP(3,IPLOT)
     *          + V(I6)*4.D0*SHP(3,IPLOT)*SHP(1,IPLOT)
C
           NSP(IPLOT) = MAX ( 1 ,
     *         INT(NRK*DT*SQRT((DXP*DXP+DYP*DYP)*SURDET(IEL))) )
            NSPMAX = MAX ( NSPMAX , NSP(IPLOT) )
C
         ENDIF         
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C  POUR TOUT PAS DE R-K REPETER
C-----------------------------------------------------------------------
C
       DO 20 ISP = 1 , NSPMAX
C
C-----------------------------------------------------------------------
C  LOCALISATION DU POINT D'ARRIVEE DE TOUTES LES CARACTERISTIQUES
C-----------------------------------------------------------------------
C
        DO 30 IPLOT = 1 ,NPLOT
C
         IF (ISP.LE.NSP(IPLOT)) THEN
C
               IEL = ELT(IPLOT)
               I1 = IKLE(IEL,1)
               I2 = IKLE(IEL,2)
               I3 = IKLE(IEL,3)
               I4 = IKLE(IEL,4)
               I5 = IKLE(IEL,5)
               I6 = IKLE(IEL,6)
               PAS = SENS * DT / NSP(IPLOT)
C
          DX(IPLOT) = ( U(I1)*(2.D0*SHP(1,IPLOT)-1.D0)*SHP(1,IPLOT)
     *                + U(I2)*(2.D0*SHP(2,IPLOT)-1.D0)*SHP(2,IPLOT)
     *                + U(I3)*(2.D0*SHP(3,IPLOT)-1.D0)*SHP(3,IPLOT)
     *                + U(I4)*4.D0*SHP(1,IPLOT)*SHP(2,IPLOT)
     *                + U(I5)*4.D0*SHP(2,IPLOT)*SHP(3,IPLOT)
     *                + U(I6)*4.D0*SHP(3,IPLOT)*SHP(1,IPLOT)) * PAS
          DY(IPLOT) = ( V(I1)*(2.D0*SHP(1,IPLOT)-1.D0)*SHP(1,IPLOT)
     *                + V(I2)*(2.D0*SHP(2,IPLOT)-1.D0)*SHP(2,IPLOT)
     *                + V(I3)*(2.D0*SHP(3,IPLOT)-1.D0)*SHP(3,IPLOT)
     *                + V(I4)*4.D0*SHP(1,IPLOT)*SHP(2,IPLOT)
     *                + V(I5)*4.D0*SHP(2,IPLOT)*SHP(3,IPLOT)
     *                + V(I6)*4.D0*SHP(3,IPLOT)*SHP(1,IPLOT)) * PAS
C
               XP = XPLOT(IPLOT) + DX(IPLOT)
               YP = YPLOT(IPLOT) + DY(IPLOT)
C
               SHP(1,IPLOT) = ((X(I3)-X(I2))*(YP-Y(I2))
     *                        -(Y(I3)-Y(I2))*(XP-X(I2)))* SURDET(IEL) 
               SHP(2,IPLOT) = ((X(I1)-X(I3))*(YP-Y(I3))
     *                        -(Y(I1)-Y(I3))*(XP-X(I3)))* SURDET(IEL)
               SHP(3,IPLOT) = ((X(I2)-X(I1))*(YP-Y(I1))
     *                        -(Y(I2)-Y(I1))*(XP-X(I1)))* SURDET(IEL)
C
               XPLOT(IPLOT) = XP
               YPLOT(IPLOT) = YP
C
            ENDIF
C
30        CONTINUE
C           
C-----------------------------------------------------------------------
C  TRAITEMENT PARTICULIER POUR LES CARACTERISTIQUES SORTIES
C  DE L'ELEMENT DE DEPART
C-----------------------------------------------------------------------
C
         DO 40 IPLOT = 1 , NPLOT
C
         IF(ISP.LE.NSP(IPLOT)) THEN
C
50       CONTINUE
C
            ISO = 0
            IF (SHP(1,IPLOT).LT.EPSILO) ISO = 1
            IF (SHP(2,IPLOT).LT.EPSILO) ISO = ISO + 2
            IF (SHP(3,IPLOT).LT.EPSILO) ISO = ISO + 4
C
            IF (ISO.NE.0) THEN
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QU'ON EST SORTI DE L'ELEMENT
C-----------------------------------------------------------------------
C
               IEL =ELT(IPLOT)               
               XP = XPLOT(IPLOT)
               YP = YPLOT(IPLOT)
C
               IF     (ISO.EQ.1) THEN
                  IFA = 2
               ELSEIF (ISO.EQ.2) THEN
                  IFA = 3
               ELSEIF (ISO.EQ.4) THEN
                  IFA = 1
               ELSEIF (ISO.EQ.3) THEN
                  IFA = 2
                  IF (DX(IPLOT)*(Y(IKLE(IEL,3))-YP).LT.
     *                DY(IPLOT)*(X(IKLE(IEL,3))-XP)) IFA = 3
               ELSEIF (ISO.EQ.6) THEN
                  IFA = 3
                  IF (DX(IPLOT)*(Y(IKLE(IEL,1))-YP).LT.
     *                DY(IPLOT)*(X(IKLE(IEL,1))-XP)) IFA = 1
               ELSE
                  IFA = 1
                  IF (DX(IPLOT)*(Y(IKLE(IEL,2))-YP).LT.
     *                DY(IPLOT)*(X(IKLE(IEL,2))-XP)) IFA = 2
               ENDIF
C
               IEL = IFABOR(IEL,IFA)  
C
               IF (IEL.GT.0) THEN
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE EST INTERNE AU DOMAINE
C  ON SE RELOCALISE DANS L'ELEMENT ADJACENT
C-----------------------------------------------------------------------
C
               I1 = IKLE(IEL,1)
               I2 = IKLE(IEL,2)
               I3 = IKLE(IEL,3)
               I4 = IKLE(IEL,4)
               I5 = IKLE(IEL,5)
               I6 = IKLE(IEL,6)
C
               ELT(IPLOT) = IEL
               SHP(1,IPLOT) = ((X(I3)-X(I2))*(YP-Y(I2))
     *                     -(Y(I3)-Y(I2))*(XP-X(I2)))* SURDET(IEL)
               SHP(2,IPLOT) = ((X(I1)-X(I3))*(YP-Y(I3))
     *                     -(Y(I1)-Y(I3))*(XP-X(I3)))* SURDET(IEL)
               SHP(3,IPLOT) = ((X(I2)-X(I1))*(YP-Y(I1))
     *                     -(Y(I2)-Y(I1))*(XP-X(I1)))* SURDET(IEL)
C
               GOTO 50
C
               ENDIF
C
               DXP = DX(IPLOT)
               DYP = DY(IPLOT)
               I1  = IKLE(ELT(IPLOT),IFA)
               I2  = IKLE(ELT(IPLOT),ISUI(IFA))
               DX1 = X(I2) - X(I1)
               DY1 = Y(I2) - Y(I1)
C
               IF(IEL.EQ.-1) THEN
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE EST UNE FRONTIERE SOLIDE
C  ON PROJETTE LE RELIQUAT SUR LA FRONTIERE PUIS ON SE RELOCALISE
C-----------------------------------------------------------------------
C
                  A1 = (DXP*DX1 + DYP*DY1) / (DX1*DX1 + DY1*DY1)
                  DX(IPLOT) = A1 * DX1
                  DY(IPLOT) = A1 * DY1
C

                A1 = ((XP-X(I1))*DX1+(YP-Y(I1))*DY1)/(DX1*DX1+DY1*DY1)
                  SHP(IFA,IPLOT) = 1.D0-A1
                  SHP(     ISUI(IFA)   ,IPLOT) = A1
                  SHP(ISUI(ISUI(IFA))  ,IPLOT) = 0.D0
    
                  XPLOT(IPLOT) = X(I1) + A1 * DX1
                  YPLOT(IPLOT) = Y(I1) + A1 * DY1
C
                 GOTO 50
C
              ENDIF
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE EST UNE FRONTIERE LIQUIDE
C  ON ARRETE LA REMONTEE DES CARACTERISTIQUE (SIGNE DE ELT)
C
C     OU QUE
C
C  LA, ON SAIT QUE LA FACE DE SORTIE EST UNE INTERFACE DE SOUS-DOMAINES
C  POINT D'INTERFACE QUI SERA TRAITE PAR LE SOUS-DOMAINE VOISIN
C  ON SE CONTENTE DE METTRE ICI TEST A ZERO
C-----------------------------------------------------------------------
C
               DENOM = DXP*DY1-DYP*DX1
               IF(DENOM.NE.0.D0) THEN
                 A1  = (DXP*(YP-Y(I1))-DYP*(XP-X(I1))) / DENOM
               ELSE
                 A1  = 0.D0
               ENDIF
               A1 = MAX(MIN(A1,1.D0),0.D0)
                SHP(IFA,IPLOT) = 1.D0-A1
                SHP(     ISUI(IFA)   ,IPLOT) = A1
                SHP(ISUI(ISUI(IFA))  ,IPLOT) = 0.D0
               XPLOT(IPLOT) = X(I1) + A1 * DX1
               YPLOT(IPLOT) = Y(I1) + A1 * DY1
               ELT(IPLOT) = - SENS * ELT(IPLOT)
               NSP(IPLOT) = ISP
C
               IF(IEL.EQ.-2.AND.NCSIZE.GT.1) THEN
                 IF(LNG.EQ.1) THEN
                 WRITE(LU,*) 'CHAR13 : PARALLELISME NON TRAITE'
                 ENDIF
                 IF(LNG.EQ.2) THEN
                 WRITE(LU,*) 'CHAR13: PARALLELISM NOT TREATED'
                 ENDIF
                 CALL PLANTE(1)
                 STOP
               ENDIF
C
            ENDIF 
C
         ENDIF
C
40       CONTINUE
C
C
20       CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
