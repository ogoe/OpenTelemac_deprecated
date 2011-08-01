C                       *****************
                        SUBROUTINE CHAR41
C                       *****************
C
     * ( U , V , W , DT , NRK , X , Y , ZSTAR , Z , IKLE2 , IBOR ,
     *   XPLOT , YPLOT , ZPLOT , DX , DY , DZ , SHP , SHZ , ELT , ETA ,
     *   NSP , NPLOT , NPOIN2 , NELEM2 , NPLAN , SURDET ,
     *   SENS , ISO_USELESS , TEST )
C
C***********************************************************************
C BIEF VERSION 6.0     16/02/2010     J-M HERVOUET (LNHE) 01 30 87 80 18
C                        28/04/93     J-M JANIN (LNH) 30 87 72 84
C
C 08/11/2004 : ADAPTATION A LA TRANSFORMEE SIGMA GENERALISEE
C 12/10/2005 : BUG CORRIGE, VOIR VARIABLE IELE QUI ETAIT AVANT IEL
C              ET EFFACAIT UN AUTRE IEL.
C 28/08/2008 : INVERSION DES BOUCLES SUR NSP ET IPLOT
C 23/01/2009 : CORRECTION D'UN BUG QUAND ON REFAIT LE POINT AU
C              FRANCHISSEMENT D'UN PLAN (VOIR PAS2)
C
C***********************************************************************
C
C  FONCTION :
C
C     REMONTEE OU DESCENTE
C     DES COURBES CARACTERISTIQUES
C     SUR DES PRISMES DE TELEMAC-3D
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
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    U,V,W       | -->| COMPOSANTE DE LA VITESSE DU CONVECTEUR       |
C |    DT          | -->| PAS DE TEMPS.                                |
C |    NRK         | -->| NOMBRE DE SOUS-PAS DE RUNGE-KUTTA.           |
C |    X,Y,ZSTAR   | -->| COORDONNEES DES POINTS DU MAILLAGE.          |
C |    Z           | -->| COTE DANS LE MAILLAGE REEL                   |
C |    IKLE2       | -->| TRANSITION ENTRE LES NUMEROTATIONS LOCALE    |
C |                |    | ET GLOBALE DU MAILLAGE 2D.                   |
C |    IBOR        | -->| NUMEROS 2D DES ELEMENTS AYANT UNE FACE COMMUNE
C |                |    | AVEC L'ELEMENT .  SI IFABOR<0 OU NUL         |
C |                |    | ON A UNE FACE LIQUIDE,SOLIDE,OU PERIODIQUE   |
C |  X..,Y..,ZPLOT |<-->| POSITIONS SUCCESSIVES DES DERIVANTS.         |
C |    DX,DY,DZ    | -- | STOCKAGE DES SOUS-PAS .                      |
C |    SHP         |<-->| COORDONNEES BARYCENTRIQUES 2D AU PIED DES    |
C |                |    | COURBES CARACTERISTIQUES.                    |
C |    SHZ         |<-->| COORDONNEES BARYCENTRIQUES SUIVANT Z DES     |
C |                |    | NOEUDS DANS LEURS ETAGES "ETA" ASSOCIES.     |
C |    ELT         |<-->| NUMEROS DES ELEMENTS 2D CHOISIS POUR CHAQUE  |
C |                |    | NOEUD.                                       |
C |    ETA         |<-->| NUMEROS DES ETAGES CHOISIS POUR CHAQUE NOEUD.|
C |    NSP         | -- | NOMBRE DE SOUS-PAS DE RUNGE KUTTA.           |
C |    NPLOT       | -->| NOMBRE DE DERIVANTS.                         |
C |    NPOIN2      | -->| NOMBRE DE POINTS DU MAILLAGE 2D.             |
C |    NELEM2      | -->| NOMBRE D'ELEMENTS DU MAILLAGE 2D.            |
C |    NPLAN       | -->| NOMBRE DE PLANS.                             |
C |    SURDET      | -->| VARIABLE UTILISEE PAR LA TRANSFORMEE ISOPARAM.
C |    SENS        | -->| DESCENTE OU REMONTEE DES CARACTERISTIQUES.   |
C |    ISO         | -- | STOCKAGE BINAIRE DE LA FACE DE SORTIE.       |
C |________________|____|______________________________________________|
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C     - APPELE PAR : CARACT
C     - PROGRAMMES APPELES : NEANT
C
C***********************************************************************
C
      USE BIEF, EX_CHAR41 => CHAR41
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)    :: SENS,NPLAN
      INTEGER         , INTENT(IN)    :: NPOIN2,NELEM2,NPLOT,NRK
      INTEGER         , INTENT(IN)    :: IKLE2(NELEM2,3)
      INTEGER         , INTENT(INOUT) :: ELT(NPLOT),NSP(NPLOT)
      DOUBLE PRECISION, INTENT(IN)    :: U(NPOIN2,NPLAN),V(NPOIN2,NPLAN)
      DOUBLE PRECISION, INTENT(IN)    :: W(NPOIN2,NPLAN),SURDET(NELEM2)
      DOUBLE PRECISION, INTENT(INOUT) :: XPLOT(NPLOT),YPLOT(NPLOT)
      DOUBLE PRECISION, INTENT(INOUT) :: ZPLOT(NPLOT)
      DOUBLE PRECISION, INTENT(INOUT) :: SHP(3,NPLOT),SHZ(NPLOT)
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN2),Y(NPOIN2),DT
      DOUBLE PRECISION, INTENT(IN)    :: Z(NPOIN2,NPLAN),ZSTAR(NPLAN)
      DOUBLE PRECISION, INTENT(INOUT) :: DX(NPLOT),DY(NPLOT),TEST(NPLOT)
      DOUBLE PRECISION, INTENT(INOUT) :: DZ(NPLOT)
      INTEGER         , INTENT(IN)    :: IBOR(NELEM2,5,NPLAN-1)
      INTEGER         , INTENT(INOUT) :: ETA(NPLOT),ISO_USELESS(NPLOT)
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELE,ISO
      INTEGER IPLOT,ISP,I1,I2,I3,IEL,IET,IET2,ISOH,ISOV,IFA,ISUI(3)
C
      DOUBLE PRECISION PAS,EPSILO,A1,DX1,DY1,DXP,DYP,XP,YP,ZP,DENOM
      DOUBLE PRECISION DELTAZ,EPSDZ,PAS2
C
      INTRINSIC ABS , INT , MAX , SQRT
C
      DATA ISUI   / 2 , 3 , 1 /
      DATA EPSILO / -1.D-6 /
      DATA EPSDZ  / 1.D-4 /
C
C-----------------------------------------------------------------------
C  CALCUL DU NOMBRE DE SOUS PAS, LE MEME POUR TOUS LES POINTS
C-----------------------------------------------------------------------
C
      DO 10 IPLOT = 1 , NPLOT
C
        IEL = ELT(IPLOT)
        IET = ETA(IPLOT)
C
        I1 = IKLE2(IEL,1)
        I2 = IKLE2(IEL,2)
        I3 = IKLE2(IEL,3)
C
        DXP = U(I1,IET  )*SHP(1,IPLOT)*(1.D0-SHZ(IPLOT))
     *      + U(I2,IET  )*SHP(2,IPLOT)*(1.D0-SHZ(IPLOT))
     *      + U(I3,IET  )*SHP(3,IPLOT)*(1.D0-SHZ(IPLOT))
     *      + U(I1,IET+1)*SHP(1,IPLOT)*SHZ(IPLOT)
     *      + U(I2,IET+1)*SHP(2,IPLOT)*SHZ(IPLOT)
     *      + U(I3,IET+1)*SHP(3,IPLOT)*SHZ(IPLOT)
C
        DYP = V(I1,IET  )*SHP(1,IPLOT)*(1.D0-SHZ(IPLOT))
     *      + V(I2,IET  )*SHP(2,IPLOT)*(1.D0-SHZ(IPLOT))
     *      + V(I3,IET  )*SHP(3,IPLOT)*(1.D0-SHZ(IPLOT))
     *      + V(I1,IET+1)*SHP(1,IPLOT)*SHZ(IPLOT)
     *      + V(I2,IET+1)*SHP(2,IPLOT)*SHZ(IPLOT)
     *      + V(I3,IET+1)*SHP(3,IPLOT)*SHZ(IPLOT)
C
        NSP(IPLOT)=MAX(1,INT(NRK*DT*SQRT((DXP**2+DYP**2)*SURDET(IEL))))
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C  POUR TOUT PAS DE R-K REPETER
C-----------------------------------------------------------------------
C
      DO IPLOT=1,NPLOT
C
      PAS = SENS * DT / NSP(IPLOT)
C
      DO ISP = 1 , NSP(IPLOT)
C
C-----------------------------------------------------------------------
C  LOCALISATION DU POINT D'ARRIVEE DE TOUTES LES CARACTERISTIQUES
C-----------------------------------------------------------------------
C
        ISO = 0
C       
        PAS2=PAS
C
        IEL = ELT(IPLOT)
        IET = ETA(IPLOT)
C
        I1 = IKLE2(IEL,1)
        I2 = IKLE2(IEL,2)
        I3 = IKLE2(IEL,3)
C
        DX(IPLOT) = ( U(I1,IET  )*SHP(1,IPLOT)*(1.D0-SHZ(IPLOT))
     *              + U(I2,IET  )*SHP(2,IPLOT)*(1.D0-SHZ(IPLOT))
     *              + U(I3,IET  )*SHP(3,IPLOT)*(1.D0-SHZ(IPLOT))
     *              + U(I1,IET+1)*SHP(1,IPLOT)*SHZ(IPLOT)
     *              + U(I2,IET+1)*SHP(2,IPLOT)*SHZ(IPLOT)
     *              + U(I3,IET+1)*SHP(3,IPLOT)*SHZ(IPLOT) ) * PAS
C
               DY(IPLOT) = ( V(I1,IET  )*SHP(1,IPLOT)*(1.D0-SHZ(IPLOT))
     *                     + V(I2,IET  )*SHP(2,IPLOT)*(1.D0-SHZ(IPLOT))
     *                     + V(I3,IET  )*SHP(3,IPLOT)*(1.D0-SHZ(IPLOT))
     *                     + V(I1,IET+1)*SHP(1,IPLOT)*SHZ(IPLOT)
     *                     + V(I2,IET+1)*SHP(2,IPLOT)*SHZ(IPLOT)
     *                     + V(I3,IET+1)*SHP(3,IPLOT)*SHZ(IPLOT) ) * PAS
C
               DELTAZ =  (Z(I1,IET+1)-Z(I1,IET))*SHP(1,IPLOT)
     *                 + (Z(I2,IET+1)-Z(I2,IET))*SHP(2,IPLOT)
     *                 + (Z(I3,IET+1)-Z(I3,IET))*SHP(3,IPLOT)
!
               IF(DELTAZ.GT.EPSDZ) THEN
               DZ(IPLOT) = ( W(I1,IET  )*SHP(1,IPLOT)*(1.D0-SHZ(IPLOT))
     *                     + W(I2,IET  )*SHP(2,IPLOT)*(1.D0-SHZ(IPLOT))
     *                     + W(I3,IET  )*SHP(3,IPLOT)*(1.D0-SHZ(IPLOT))
     *                     + W(I1,IET+1)*SHP(1,IPLOT)*SHZ(IPLOT)
     *                     + W(I2,IET+1)*SHP(2,IPLOT)*SHZ(IPLOT)
     *                     + W(I3,IET+1)*SHP(3,IPLOT)*SHZ(IPLOT) )
     *                     * PAS * (ZSTAR(IET+1)-ZSTAR(IET)) / DELTAZ
               ELSE
                 DZ(IPLOT) = 0.D0
               ENDIF
C
               XP = XPLOT(IPLOT) + DX(IPLOT)
               YP = YPLOT(IPLOT) + DY(IPLOT)
               ZP = ZPLOT(IPLOT) + DZ(IPLOT)
C
               SHP(1,IPLOT) = ((X(I3)-X(I2))*(YP-Y(I2))
     *                        -(Y(I3)-Y(I2))*(XP-X(I2))) * SURDET(IEL)
               SHP(2,IPLOT) = ((X(I1)-X(I3))*(YP-Y(I3))
     *                        -(Y(I1)-Y(I3))*(XP-X(I3))) * SURDET(IEL)
               SHP(3,IPLOT) = ((X(I2)-X(I1))*(YP-Y(I1))
     *                        -(Y(I2)-Y(I1))*(XP-X(I1))) * SURDET(IEL)
               SHZ(IPLOT) = (ZP-ZSTAR(IET)) / (ZSTAR(IET+1)-ZSTAR(IET))
C
               XPLOT(IPLOT) = XP
               YPLOT(IPLOT) = YP
               ZPLOT(IPLOT) = ZP
C
               IF(SHP(1,IPLOT).LT.EPSILO) ISO=IBSET(ISO,2)
               IF(SHP(2,IPLOT).LT.EPSILO) ISO=IBSET(ISO,3)
               IF(SHP(3,IPLOT).LT.EPSILO) ISO=IBSET(ISO,4)
C
               IF(SHZ(IPLOT).LT.EPSILO)      ISO=IBSET(ISO,0)
               IF(SHZ(IPLOT).GT.1.D0-EPSILO) ISO=IBSET(ISO,1)
C
C-----------------------------------------------------------------------
C  TRAITEMENT PARTICULIER POUR LES CARACTERISTIQUES SORTIES
C  DE L'ELEMENT DE DEPART
C-----------------------------------------------------------------------
C
50          CONTINUE
C  
C        
           IF(ISO.NE.0) THEN
C            
C-----------------------------------------------------------------------
C  LA, ON SAIT QU'ON EST SORTI DE L'ELEMENT
C-----------------------------------------------------------------------
C
               ISOH = IAND(ISO,28)
               ISOV = IAND(ISO, 3)
               IEL = ELT(IPLOT)
               IET = ETA(IPLOT)
               XP = XPLOT(IPLOT)
               YP = YPLOT(IPLOT)
               ZP = ZPLOT(IPLOT)
C
               IF(ISOH.NE.0) THEN
                  IF (ISOH.EQ.4) THEN
                     IFA = 2
                  ELSEIF (ISOH.EQ.8) THEN
                     IFA = 3
                  ELSEIF (ISOH.EQ.16) THEN
                     IFA = 1
                  ELSEIF (ISOH.EQ.12) THEN
                     IFA = 2
                     IF (DX(IPLOT)*(Y(IKLE2(IEL,3))-YP).LT.
     *                   DY(IPLOT)*(X(IKLE2(IEL,3))-XP)) IFA = 3
                  ELSEIF (ISOH.EQ.24) THEN
                     IFA = 3
                     IF (DX(IPLOT)*(Y(IKLE2(IEL,1))-YP).LT.
     *                   DY(IPLOT)*(X(IKLE2(IEL,1))-XP)) IFA = 1
                  ELSE
                     IFA = 1
                     IF (DX(IPLOT)*(Y(IKLE2(IEL,2))-YP).LT.
     *                   DY(IPLOT)*(X(IKLE2(IEL,2))-XP)) IFA = 2
                  ENDIF
                  IF(ISOV.GT.0) THEN
                     IF(ABS(DZ(IPLOT)).GT.EPSDZ) THEN
                       A1 = (ZP-ZSTAR(IET+ISOV-1)) / DZ(IPLOT)
                     ELSE
                       A1 = 0.D0
                     ENDIF
                     I1 = IKLE2(IEL,IFA)
                     I2 = IKLE2(IEL,ISUI(IFA))
                     IF ((X(I2)-X(I1))*(YP-A1*DY(IPLOT)-Y(I1)).GT.
     *                 (Y(I2)-Y(I1))*(XP-A1*DX(IPLOT)-X(I1))) IFA=ISOV+3
                  ENDIF
               ELSE
                  IFA = ISOV + 3
               ENDIF
C
               IEL = IBOR(IEL,IFA,IET)
C
               IF (IFA.LE.3) THEN
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE DU PRISME EST UNE FACE QUADRAN.
C  =================================================================
C-----------------------------------------------------------------------
C
                  IF(IEL.GT.0) THEN
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE EST INTERNE AU DOMAINE
C  ON SE RELOCALISE DANS L'ELEMENT ADJACENT
C-----------------------------------------------------------------------
C
                    I1 = IKLE2(IEL,1)
                    I2 = IKLE2(IEL,2)
                    I3 = IKLE2(IEL,3)
C
                    ELT(IPLOT) = IEL
                    SHP(1,IPLOT) = ((X(I3)-X(I2))*(YP-Y(I2))
     *                            -(Y(I3)-Y(I2))*(XP-X(I2)))*SURDET(IEL)
                    SHP(2,IPLOT) = ((X(I1)-X(I3))*(YP-Y(I3))
     *                            -(Y(I1)-Y(I3))*(XP-X(I3)))*SURDET(IEL)
                    SHP(3,IPLOT) = ((X(I2)-X(I1))*(YP-Y(I1))
     *                            -(Y(I2)-Y(I1))*(XP-X(I1)))*SURDET(IEL)
C
                    ISO = ISOV
C
                    IF(SHP(1,IPLOT).LT.EPSILO) ISO=IBSET(ISO,2)
                    IF(SHP(2,IPLOT).LT.EPSILO) ISO=IBSET(ISO,3)
                    IF(SHP(3,IPLOT).LT.EPSILO) ISO=IBSET(ISO,4)
C
                    GOTO 50
C
                  ENDIF
C
                  DXP = DX(IPLOT)
                  DYP = DY(IPLOT)
                  I1  = IKLE2(ELT(IPLOT),IFA)
                  I2  = IKLE2(ELT(IPLOT),ISUI(IFA))
                  DX1 = X(I2) - X(I1)
                  DY1 = Y(I2) - Y(I1)
C
                  IF(IEL.EQ.-1) THEN
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE EST UNE FRONTIERE SOLIDE
C  ON PROJETTE LE RELIQUAT SUR LA FRONTIERE PUIS SE RELOCALISE
C-----------------------------------------------------------------------
C
                    A1 = (DXP*DX1+DYP*DY1) / (DX1**2+DY1**2)
C                   DEPLACEMENT LE LONG DE LA FRONTIERE A PARTIR DU
C                   POINT DE SORTIE
                    DX(IPLOT) = A1 * DX1
                    DY(IPLOT) = A1 * DY1
C
                    A1=((XP-X(I1))*DX1+(YP-Y(I1))*DY1)/(DX1**2+DY1**2)
                    SHP(          IFA  ,IPLOT) = 1.D0 - A1
                    SHP(     ISUI(IFA) ,IPLOT) = A1
                    SHP(ISUI(ISUI(IFA)),IPLOT) = 0.D0
                    XPLOT(IPLOT) = X(I1) + A1 * DX1
                    YPLOT(IPLOT) = Y(I1) + A1 * DY1
C
                    ISO = ISOV
C
                    IF(SHP(1,IPLOT).LT.EPSILO) ISO=IBSET(ISO,2)
                    IF(SHP(2,IPLOT).LT.EPSILO) ISO=IBSET(ISO,3)
                    IF(SHP(3,IPLOT).LT.EPSILO) ISO=IBSET(ISO,4)
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
C>>>>
C                 A1 = (DXP*(YP-Y(I1))-DYP*(XP-X(I1)))/(DXP*DY1-DYP*DX1)
                  DENOM = DXP*DY1-DYP*DX1
                  IF(ABS(DENOM).GT.1.D-8) THEN
                     A1 = (DXP*(YP-Y(I1))-DYP*(XP-X(I1))) / DENOM
                  ELSE
                     A1 = 0.D0
                  ENDIF
C<<<<
                  IF (A1.GT.1.D0) A1 = 1.D0
                  IF (A1.LT.0.D0) A1 = 0.D0
                  SHP(          IFA  ,IPLOT) = 1.D0 - A1
                  SHP(     ISUI(IFA) ,IPLOT) = A1
                  SHP(ISUI(ISUI(IFA)),IPLOT) = 0.D0
                  XPLOT(IPLOT) = X(I1) + A1 * DX1
                  YPLOT(IPLOT) = Y(I1) + A1 * DY1
                  IF(ABS(DXP).GT.ABS(DYP)) THEN
                     A1 = (XP-XPLOT(IPLOT))/DXP
                  ELSE
                     A1 = (YP-YPLOT(IPLOT))/DYP
                  ENDIF
                  ZPLOT(IPLOT) = ZP - A1*DZ(IPLOT)
                  SHZ(IPLOT) = (ZPLOT(IPLOT)-ZSTAR(IET))
     *                       / (ZSTAR(IET+1)-ZSTAR(IET))
                  ELT(IPLOT) = - SENS * ELT(IPLOT)
                  NSP(IPLOT) = ISP
C
                  IF(IEL.EQ.-2.AND.NCSIZE.GT.1) TEST(IPLOT) = 0.D0
C
               ELSE
C
C-----------------------------------------------------------------------
C  CAS OU IFA = 4 OU 5 
C  LA, ON SAIT QUE LA FACE DE SORTIE DU PRISME EST UNE FACE TRIANGULAIRE
C  =====================================================================
C-----------------------------------------------------------------------
C
                  IFA = IFA - 4
C                 HENCE IFA NOW EQUALS 0 OR 1
C
                  IF (IEL.EQ.1) THEN
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE EST INTERNE AU DOMAINE
C  ET ON N'A PAS BESOIN DE RECALCULER LES VITESSES
C  ON SE RELOCALISE DANS L'ELEMENT ADJACENT
C-----------------------------------------------------------------------
C
                     ETA(IPLOT) = IET + IFA + IFA - 1
                     SHZ(IPLOT) = (ZP-ZSTAR(ETA(IPLOT)))
     *                   / (ZSTAR(ETA(IPLOT)+1)-ZSTAR(ETA(IPLOT)))
                     ISO = ISOH
                     IF(SHZ(IPLOT).LT.EPSILO)      ISO=IBSET(ISO,0)
                     IF(SHZ(IPLOT).GT.1.D0-EPSILO) ISO=IBSET(ISO,1)
                     GOTO 50
C
                  ENDIF
C
                  IF(IEL.EQ.-1) THEN
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE EST UNE FRONTIERE SOLIDE
C  ON PROJETTE LE RELIQUAT SUR LA FRONTIERE PUIS ON SE RELOCALISE
C-----------------------------------------------------------------------
C
                     ZPLOT(IPLOT) = ZSTAR(IET+IFA)
                     DZ   (IPLOT) = 0.D0
                     SHZ  (IPLOT) = IFA
C
                     ISO = ISOH
                     IF(ISOH.NE.0) GOTO 50
C
                  ELSE
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE EST UNE FRONTIERE LIQUIDE (CAS 0)
C      ON ARRETE ALORS LA REMONTEE DES CARACTERISTIQUES (SIGNE DE ELT)
C  OU, QUE L'ON VIENT DE TRAVERSER UN PLAN AVEC RECALCUL DES VITESSES
C  DEMANDE (CAS 2)
C-----------------------------------------------------------------------
C
                     IF(ABS(DZ(IPLOT)).GT.EPSDZ) THEN
                       A1 = (ZP-ZSTAR(IET+IFA)) / DZ(IPLOT)
                     ELSE
                       A1 = 0.D0
                     ENDIF
C                    RECUL JUSQU'AU POINT DE FRANCHISSEMENT
                     XP = XP - A1*DX(IPLOT)
                     YP = YP - A1*DY(IPLOT)
                     ZP = ZSTAR(IET+IFA)
                     IELE = ELT(IPLOT)
                     I1 = IKLE2(IELE,1)
                     I2 = IKLE2(IELE,2)
                     I3 = IKLE2(IELE,3)
C
                     SHP(1,IPLOT) = ((X(I3)-X(I2))*(YP-Y(I2))
     *                           -(Y(I3)-Y(I2))*(XP-X(I2)))*SURDET(IELE)
                     SHP(2,IPLOT) = ((X(I1)-X(I3))*(YP-Y(I3))
     *                           -(Y(I1)-Y(I3))*(XP-X(I3)))*SURDET(IELE)
                     SHP(3,IPLOT) = ((X(I2)-X(I1))*(YP-Y(I1))
     *                           -(Y(I2)-Y(I1))*(XP-X(I1)))*SURDET(IELE)
C
                     IF(IEL.EQ.2) THEN
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE SE SITUE SUR UN PLAN OU ON DEMANDE
C  UN RECALCUL DES VITESSES
C-----------------------------------------------------------------------
C
C                       IF IFA = 1 EXIT THROUGH THE TOP
C                       IF IFA = 0 EXIT THROUGH THE BOTTOM
C                       THEN NEW IET IS  IET+1 IF IFA = 1
C                                    AND IET-1 IF IFA = 0
C                       THIS IS SUMMARISED BY IET=IET+2*IFA-1
C
C                       RECOMPUTED VELOCITIES MUST BE TAKEN AT IET2=IET+IFA
C                       I.E. BOTTOM IF EXIT THROUGH THE BOTTOM
C                           AND TOP IF EXIT THROUGH THE TOP
C
                        IET2   = IET + IFA
                        IET    = IET + IFA + IFA - 1
C                       REMAINING TIME IS REDUCED
                        PAS2 = PAS2 * A1 
C
                        DX(IPLOT) = ( U(I1,IET2)*SHP(1,IPLOT)
     *                              + U(I2,IET2)*SHP(2,IPLOT)
     *                              + U(I3,IET2)*SHP(3,IPLOT) ) * PAS2
C
                        DY(IPLOT) = ( V(I1,IET2)*SHP(1,IPLOT)
     *                              + V(I2,IET2)*SHP(2,IPLOT)
     *                              + V(I3,IET2)*SHP(3,IPLOT) ) * PAS2
C
                        DELTAZ =  (Z(I1,IET+1)-Z(I1,IET))*SHP(1,IPLOT)
     *                          + (Z(I2,IET+1)-Z(I2,IET))*SHP(2,IPLOT)
     *                          + (Z(I3,IET+1)-Z(I3,IET))*SHP(3,IPLOT)
C 
                        IF(DELTAZ.GT.EPSDZ) THEN
                          DZ(IPLOT)=(W(I1,IET2)*SHP(1,IPLOT)
     *                             + W(I2,IET2)*SHP(2,IPLOT)
     *                             + W(I3,IET2)*SHP(3,IPLOT))*PAS2
     *                             *(ZSTAR(IET+1)-ZSTAR(IET))/DELTAZ
                        ELSE
                          DZ(IPLOT) = 0.D0
                        ENDIF
C
                        XP = XP + DX(IPLOT)
                        YP = YP + DY(IPLOT)
                        ZP = ZP + DZ(IPLOT)
C
                        SHP(1,IPLOT) = ((X(I3)-X(I2))*(YP-Y(I2))
     *                        -(Y(I3)-Y(I2))*(XP-X(I2))) * SURDET(IELE)
                        SHP(2,IPLOT) = ((X(I1)-X(I3))*(YP-Y(I3))
     *                        -(Y(I1)-Y(I3))*(XP-X(I3))) * SURDET(IELE)
                        SHP(3,IPLOT) = ((X(I2)-X(I1))*(YP-Y(I1))
     *                        -(Y(I2)-Y(I1))*(XP-X(I1))) * SURDET(IELE)
                        SHZ(IPLOT)=(ZP-ZSTAR(IET))/
     *                                         (ZSTAR(IET+1)-ZSTAR(IET))
C
                        XPLOT(IPLOT) = XP
                        YPLOT(IPLOT) = YP
                        ZPLOT(IPLOT) = ZP
                        ETA(IPLOT) = IET
C
                        ISO = 0
C
                        IF(SHP(1,IPLOT).LT.EPSILO) ISO=IBSET(ISO,2)
                        IF(SHP(2,IPLOT).LT.EPSILO) ISO=IBSET(ISO,3)
                        IF(SHP(3,IPLOT).LT.EPSILO) ISO=IBSET(ISO,4)
C
                        IF(SHZ(IPLOT).LT.EPSILO)      ISO=IBSET(ISO,0)
                        IF(SHZ(IPLOT).GT.1.D0-EPSILO) ISO=IBSET(ISO,1)
C
                        GOTO 50
C
                     ENDIF
C
                     XPLOT(IPLOT) = XP
                     YPLOT(IPLOT) = YP
                     ZPLOT(IPLOT) = ZP
                     SHZ  (IPLOT) = IFA
                     ELT  (IPLOT) = - SENS * ELT(IPLOT)
                     EXIT
C
                  ENDIF
C
               ENDIF
C
            ENDIF
C
      ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
