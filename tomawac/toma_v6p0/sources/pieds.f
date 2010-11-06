C                       *****************
                        SUBROUTINE PIEDS
C                       *****************
C
     *  (U , V , W , DT , NRK , X , Y , TETA , IKLE2 , IFABOR , ETAS ,
     *   XPLOT , YPLOT , ZPLOT , DX , DY , DZ , SHP1 , SHP2 , SHP3 ,
     *   SHZ , ELT , ETA , NSP , NPLOT , NPOIN2 , NELEM2 , NPLAN ,
     *   IFF , SURDET , SENS , ISO )
C
C***********************************************************************
C  TOMAWAC  VERSION 1.0       01/02/95        F MARCOS (LNH) 30 87 72 66
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
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    U,V,W       ! -->! COMPOSANTE DE LA VITESSE DU CONVECTEUR       !
C !    DT          ! -->! PAS DE TEMPS.                                !
C !    NRK         ! -->! NOMBRE DE SOUS-PAS DE RUNGE-KUTTA.           !
C !    X,Y,TETA    ! -->! COORDONNEES DES POINTS DU MAILLAGE.          !
C !    IKLE2       ! -->! TRANSITION ENTRE LES NUMEROTATIONS LOCALE    !
C !                !    ! ET GLOBALE DU MAILLAGE 2D.                   !
C !    IFABOR      ! -->! NUMEROS 2D DES ELEMENTS AYANT UNE FACE COMMUNE
C !                !    ! AVEC L'ELEMENT .  SI IFABOR<0 OU NUL         !
C !                !    ! ON A UNE FACE LIQUIDE,SOLIDE,OU PERIODIQUE   !
C !    ETAS        !<-->! TABLEAU DE TRAVAIL DONNANT LE NUMERO DE      !
C !                !    ! L'ETAGE SUPERIEUR                            !
C !  X..,Y..,ZPLOT !<-->! POSITIONS SUCCESSIVES DES DERIVANTS.         !
C !    DX,DY,DZ    ! -- ! STOCKAGE DES SOUS-PAS .                      !
C !    SHP1-2-3    !<-->! COORDONNEES BARYCENTRIQUES 2D AU PIED DES    !
C !                !    ! COURBES CARACTERISTIQUES.                    !
C !    SHZ         !<-->! COORDONNEES BARYCENTRIQUES SUIVANT Z DES     !
C !                !    ! NOEUDS DANS LEURS ETAGES "ETA" ASSOCIES.     !
C !    ELT         !<-->! NUMEROS DES ELEMENTS 2D CHOISIS POUR CHAQUE  !
C !                !    ! NOEUD.                                       !
C !    ETA         !<-->! NUMEROS DES ETAGES CHOISIS POUR CHAQUE NOEUD.!
C !    NSP         ! -- ! NOMBRE DE SOUS-PAS DE RUNGE KUTTA.           !
C !    NPLOT       ! -->! NOMBRE DE DERIVANTS.                         !
C !    NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE 2D.             !
C !    NELEM2      ! -->! NOMBRE D'ELEMENTS DU MAILLAGE 2D.            !
C !    NPLAN       ! -->! NOMBRE DE DIRECTIONS                         !
C !    SURDET      ! -->! VARIABLE UTILISEE PAR LA TRANSFORMEE ISOPARAM.
C !    SENS        ! -->! DESCENTE OU REMONTEE DES CARACTERISTIQUES.   !
C !    ISO         !<-->! INDIQUE PAR BIT LA FACE DE SORTIE DE l'ELEMEN!
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C     - APPELE PAR : WAC
C     - PROGRAMMES APPELES : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
      INTEGER NPOIN2,NELEM2,NPLAN,NPLOT,NSPMAX,NRK,SENS,IFF
C
      DOUBLE PRECISION U(NPOIN2,NPLAN),V(NPOIN2,NPLAN)
      DOUBLE PRECISION W(NPOIN2,NPLAN)
      DOUBLE PRECISION XPLOT(NPLOT),YPLOT(NPLOT),ZPLOT(NPLOT)
      DOUBLE PRECISION SURDET(NELEM2),SHZ(NPLOT)
      DOUBLE PRECISION SHP1(NPLOT),SHP2(NPLOT),SHP3(NPLOT)
      DOUBLE PRECISION X(NPOIN2),Y(NPOIN2),TETA(NPLAN+1)
      DOUBLE PRECISION DX(NPLOT),DY(NPLOT),DZ(NPLOT)
      DOUBLE PRECISION PAS,DT,A1,DX1,DY1,DXP,DYP,DZP,XP,YP,ZP
      DOUBLE PRECISION EPSILO, EPSI, EPM1
C
      INTEGER IKLE2(NELEM2,3),IFABOR(NELEM2,5),ETAS(NPLAN)
      INTEGER ELT(NPLOT),ETA(NPLOT),NSP(NPLOT),ISO(NPLOT)
      INTEGER IPLOT,ISP,I1,I2,I3,IEL,IET,ISOH,ISOV,IFA,ISUI(3)
C
      INTRINSIC ABS , INT , MAX , SQRT
C
      DATA ISUI   / 2 , 3 , 1 /
      DATA EPSILO / -1.D-6 /
      DATA EPSI   / 1.D-12 /
C
C-----------------------------------------------------------------------
C  CALCUL DU NOMBRE DE SOUS PAS, LE MEME POUR TOUS LES POINTS POUR UNE
C    FREQUENCE DONNEE
C-----------------------------------------------------------------------
C
      NSPMAX = 1
      EPM1=1.D0-EPSI
C
      DO 10 IPLOT = 1 , NPLOT
C
         NSP(IPLOT) = 0
         IEL = ELT(IPLOT)
C
         IF (IEL.GT.0) THEN
C
            IET = ETA(IPLOT)
C
            I1 = IKLE2(IEL,1)
            I2 = IKLE2(IEL,2)
            I3 = IKLE2(IEL,3)
C
         DXP = U(I1,IET  )*SHP1(IPLOT)*(1.D0-SHZ(IPLOT))
     *       + U(I2,IET  )*SHP2(IPLOT)*(1.D0-SHZ(IPLOT))
     *       + U(I3,IET  )*SHP3(IPLOT)*(1.D0-SHZ(IPLOT))
     *       + U(I1,ETAS(IET))*SHP1(IPLOT)*SHZ(IPLOT)
     *       + U(I2,ETAS(IET))*SHP2(IPLOT)*SHZ(IPLOT)
     *       + U(I3,ETAS(IET))*SHP3(IPLOT)*SHZ(IPLOT)
C
         DYP = V(I1,IET  )*SHP1(IPLOT)*(1.D0-SHZ(IPLOT))
     *       + V(I2,IET  )*SHP2(IPLOT)*(1.D0-SHZ(IPLOT))
     *       + V(I3,IET  )*SHP3(IPLOT)*(1.D0-SHZ(IPLOT))
     *       + V(I1,ETAS(IET))*SHP1(IPLOT)*SHZ(IPLOT)
     *       + V(I2,ETAS(IET))*SHP2(IPLOT)*SHZ(IPLOT)
     *       + V(I3,ETAS(IET))*SHP3(IPLOT)*SHZ(IPLOT)
C
         DZP = W(I1,IET  )*SHP1(IPLOT)*(1.D0-SHZ(IPLOT))
     *       + W(I2,IET  )*SHP2(IPLOT)*(1.D0-SHZ(IPLOT))
     *       + W(I3,IET  )*SHP3(IPLOT)*(1.D0-SHZ(IPLOT))
     *       + W(I1,ETAS(IET))*SHP1(IPLOT)*SHZ(IPLOT)
     *       + W(I2,ETAS(IET))*SHP2(IPLOT)*SHZ(IPLOT)
     *       + W(I3,ETAS(IET))*SHP3(IPLOT)*SHZ(IPLOT)
C
         NSP(IPLOT)= MAX(INT(NRK*DT*ABS(DZP/(TETA(IET)-TETA(IET+1)))),
     *         INT(NRK*DT*SQRT((DXP*DXP+DYP*DYP)*SURDET(IEL))) )
C
            NSP(IPLOT) = MAX (1,NSP(IPLOT))
C
            NSPMAX = MAX ( NSPMAX , NSP(IPLOT) )
C
         ENDIF 
C
10    CONTINUE
      IF (LNG.EQ.1) THEN
        WRITE(LU,*)
     *     '   FREQUENCE',IFF,', NOMBRE DE SOUS PAS RUNGE KUTTA :'
     *        ,NSPMAX
      ELSE
        WRITE(LU,*)
     *     '   FREQUENCY',IFF,', NUMBER OF RUNGE KUTTA SUB TIME-STEP
     *S :',NSPMAX
      ENDIF
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
         DO 30 IPLOT = 1 , NPLOT
C
            ISO(IPLOT) = 0
            IF (ISP.LE.NSP(IPLOT)) THEN
C
               IEL = ELT(IPLOT)
               IET = ETA(IPLOT)
C
               I1 = IKLE2(IEL,1)
               I2 = IKLE2(IEL,2)
               I3 = IKLE2(IEL,3)
               PAS = SENS * DT / NSP(IPLOT)
C
               DX(IPLOT) =
     * ( U(I1,IET  )*SHP1(IPLOT)*(1.D0-SHZ(IPLOT))
     * + U(I2,IET  )*SHP2(IPLOT)*(1.D0-SHZ(IPLOT))
     * + U(I3,IET  )*SHP3(IPLOT)*(1.D0-SHZ(IPLOT))
     * + U(I1,ETAS(IET))*SHP1(IPLOT)*SHZ(IPLOT)
     * + U(I2,ETAS(IET))*SHP2(IPLOT)*SHZ(IPLOT)
     * + U(I3,ETAS(IET))*SHP3(IPLOT)*SHZ(IPLOT) ) * PAS
C
               DY(IPLOT) =
     * ( V(I1,IET  )*SHP1(IPLOT)*(1.D0-SHZ(IPLOT))
     * + V(I2,IET  )*SHP2(IPLOT)*(1.D0-SHZ(IPLOT))
     * + V(I3,IET  )*SHP3(IPLOT)*(1.D0-SHZ(IPLOT))
     * + V(I1,ETAS(IET))*SHP1(IPLOT)*SHZ(IPLOT)
     * + V(I2,ETAS(IET))*SHP2(IPLOT)*SHZ(IPLOT)
     * + V(I3,ETAS(IET))*SHP3(IPLOT)*SHZ(IPLOT) ) * PAS
C
               DZ(IPLOT) =
     * ( W(I1,IET  )*SHP1(IPLOT)*(1.D0-SHZ(IPLOT))
     * + W(I2,IET  )*SHP2(IPLOT)*(1.D0-SHZ(IPLOT))
     * + W(I3,IET  )*SHP3(IPLOT)*(1.D0-SHZ(IPLOT))
     * + W(I1,ETAS(IET))*SHP1(IPLOT)*SHZ(IPLOT)
     * + W(I2,ETAS(IET))*SHP2(IPLOT)*SHZ(IPLOT)
     * + W(I3,ETAS(IET))*SHP3(IPLOT)*SHZ(IPLOT) ) * PAS
C
               XP = XPLOT(IPLOT) + DX(IPLOT)
               YP = YPLOT(IPLOT) + DY(IPLOT)
               ZP = ZPLOT(IPLOT) + DZ(IPLOT)
C
               SHP1(IPLOT) = ((X(I3)-X(I2))*(YP-Y(I2))
     *                        -(Y(I3)-Y(I2))*(XP-X(I2))) * SURDET(IEL)
               SHP2(IPLOT) = ((X(I1)-X(I3))*(YP-Y(I3))
     *                        -(Y(I1)-Y(I3))*(XP-X(I3))) * SURDET(IEL)
               SHP3(IPLOT) = ((X(I2)-X(I1))*(YP-Y(I1))
     *                        -(Y(I2)-Y(I1))*(XP-X(I1))) * SURDET(IEL)
               SHZ(IPLOT) = (ZP-TETA(IET)) / (TETA(IET+1)-TETA(IET))
C               IF (ABS(SHZ(IPLOT)).GT.2.5D0 ) THEN
C                  WRITE(LU,*)'SHZ***',IPLOT,IET,SHZ(IPLOT)
C                  WRITE(LU,*)TETA(IET),TETA(IET+1),ZP
C                  WRITE(LU,*)DZ(IPLOT),ZPLOT(IPLOT)
C                  STOP
C              ENDIF
C
               XPLOT(IPLOT) = XP
               YPLOT(IPLOT) = YP
               ZPLOT(IPLOT) = ZP
C
               IF (SHP1(IPLOT).LT.EPSILO)
     *              ISO(IPLOT)=IBSET(ISO(IPLOT),2)
               IF (SHP2(IPLOT).LT.EPSILO)
     *              ISO(IPLOT)=IBSET(ISO(IPLOT),3)
               IF (SHP3(IPLOT).LT.EPSILO)
     *              ISO(IPLOT)=IBSET(ISO(IPLOT),4)
C
               IF  (SHZ(IPLOT).LT.EPSILO)
     *              ISO(IPLOT)=IBSET(ISO(IPLOT),0)
               IF  (SHZ(IPLOT).GT.1.D0-EPSILO)
     *              ISO(IPLOT)=IBSET(ISO(IPLOT),1)
C
            ENDIF 
C
30       CONTINUE
C
C-----------------------------------------------------------------------
C  TRAITEMENT PARTICULIER POUR LES CARACTERISTIQUES SORTIES
C  DE L'ELEMENT DE DEPART
C-----------------------------------------------------------------------
C
         DO 40 IPLOT = 1 , NPLOT
C
50          CONTINUE
C
            IF (ISO(IPLOT).NE.0) THEN
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QU'ON EST SORTI DE L'ELEMENT
C-----------------------------------------------------------------------
C
               ISOH = IAND(ISO(IPLOT),28)
               ISOV = IAND(ISO(IPLOT), 3)
               IEL = ELT(IPLOT)
               IET = ETA(IPLOT)
               XP = XPLOT(IPLOT)
               YP = YPLOT(IPLOT)
               ZP = ZPLOT(IPLOT)
C
               IF (ISOH.NE.0) THEN
C
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
C
                  IF (ISOV.GT.0) THEN
                     A1 = (ZP-TETA(IET+ISOV-1)) / DZ(IPLOT)
                     I1 = IKLE2(IEL,IFA)
                     I2 = IKLE2(IEL,ISUI(IFA))
                     IF ((X(I2)-X(I1))*(YP-A1*DY(IPLOT)-Y(I1)).GT.
     *                 (Y(I2)-Y(I1))*(XP-A1*DX(IPLOT)-X(I1))) IFA=ISOV+3
                  ENDIF
C
               ELSE
C
                  IFA = ISOV + 3
C
               ENDIF
C
               IEL = IFABOR(IEL,IFA)
C
               IF (IFA.LE.3) THEN

C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE DU PRISME EST UNE FACE QUADRAN.
C  =================================================================
C-----------------------------------------------------------------------
C
                  IF (IEL.GT.0) THEN
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
                     SHP1(IPLOT) = ((X(I3)-X(I2))*(YP-Y(I2))
     *                           -(Y(I3)-Y(I2))*(XP-X(I2)))*SURDET(IEL)
                     SHP2(IPLOT) = ((X(I1)-X(I3))*(YP-Y(I3))
     *                           -(Y(I1)-Y(I3))*(XP-X(I3)))*SURDET(IEL)
                     SHP3(IPLOT) = ((X(I2)-X(I1))*(YP-Y(I1))
     *                           -(Y(I2)-Y(I1))*(XP-X(I1)))*SURDET(IEL)
C
                     ISO(IPLOT) = ISOV
C
         IF (SHP1(IPLOT).LT.EPSILO) ISO(IPLOT)=IBSET(ISO(IPLOT),2)
         IF (SHP2(IPLOT).LT.EPSILO) ISO(IPLOT)=IBSET(ISO(IPLOT),3)
         IF (SHP3(IPLOT).LT.EPSILO) ISO(IPLOT)=IBSET(ISO(IPLOT),4)
C
                     GOTO 50
C
                  ENDIF
                  DXP = DX(IPLOT)
                  DYP = DY(IPLOT)
                  I1  = IKLE2(ELT(IPLOT),IFA)
                  I2  = IKLE2(ELT(IPLOT),ISUI(IFA))
                  DX1 = X(I2) - X(I1)
                  DY1 = Y(I2) - Y(I1)
C
                  IF (IEL.EQ.-1) THEN
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE EST UNE FRONTIERE SOLIDE
C  ON MET LES SHP A 0, FIN DE LA REMONTEE
C-----------------------------------------------------------------------
C
                     SHP1(IPLOT) = 0.D0
                     SHP2(IPLOT) = 0.D0
                     SHP3(IPLOT) = 0.D0
                     ELT(IPLOT) = - SENS * ELT(IPLOT)
                     NSP(IPLOT) = ISP
                     GOTO 40
C
                  ENDIF 
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE EST UNE FRONTIERE LIQUIDE
C  ON ARRETE LA REMONTEE DES CARACTERISTIQUE (SIGNE DE ELT)
C-----------------------------------------------------------------------
C
                  A1 = (DXP*(YP-Y(I1))-DYP*(XP-X(I1)))/(DXP*DY1-DYP*DX1)
CFBG                  IF (A1.GT.1.D0) A1 = 1.D0
CFBG                  IF (A1.LT.0.D0) A1 = 0.D0
                  IF (A1.GT.EPM1) A1 = 1.D0
                  IF (A1.LT.EPSI) A1 = 0.D0
CFGB
                  IF (IFA.EQ.1) THEN
                    SHP1(IPLOT) = 1.D0 - A1
                    SHP2(IPLOT) = A1
                    SHP3(IPLOT) = 0.D0
                  ELSEIF (IFA.EQ.2) THEN
                    SHP2(IPLOT) = 1.D0 - A1
                    SHP3(IPLOT) = A1
                    SHP1(IPLOT) = 0.D0
                  ELSE
                    SHP3(IPLOT) = 1.D0 - A1
                    SHP1(IPLOT) = A1
                    SHP2(IPLOT) = 0.D0
                  ENDIF
                  XPLOT(IPLOT) = X(I1) + A1 * DX1
                  YPLOT(IPLOT) = Y(I1) + A1 * DY1
                  IF (ABS(DXP).GT.ABS(DYP)) THEN
                     A1 = (XP-XPLOT(IPLOT))/DXP
                  ELSE
                     A1 = (YP-YPLOT(IPLOT))/DYP
                  ENDIF
                  IF (A1.GT.EPM1) A1 = 1.D0
                  IF (A1.LT.EPSI) A1 = 0.D0
                  ZPLOT(IPLOT) = ZP - A1*DZ(IPLOT)
                  SHZ(IPLOT) = (ZPLOT(IPLOT)-TETA(IET))
     *                       / (TETA(IET+1)-TETA(IET))
                  ELT(IPLOT) = - SENS * ELT(IPLOT)
                  NSP(IPLOT) = ISP
C
               ELSE
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE DU PRISME EST UNE FACE TRIAN.
C  ===============================================================
C-----------------------------------------------------------------------
C
                  IFA = IFA - 4
C
                  IF (IEL.EQ.1) THEN
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE EST INTERNE AU DOMAINE
C  ON SE RELOCALISE DANS L'ELEMENT ADJACENT
C-----------------------------------------------------------------------
C
                     ETA(IPLOT) = IET + IFA + IFA - 1
                     IF (ETA(IPLOT).EQ.NPLAN+1) THEN
                         ETA(IPLOT)=1
                         ZP=ZP-2*3.14159265D0
                         ZPLOT(IPLOT)=ZP
                     ENDIF
                     IF (ETA(IPLOT).EQ.0) THEN
                         ETA(IPLOT) = NPLAN
                         ZP=ZP+2*3.14159265D0
                         ZPLOT(IPLOT)=ZP
                     ENDIF
                     SHZ(IPLOT) = (ZP-TETA(ETA(IPLOT)))
     *                   / (TETA(ETA(IPLOT)+1)-TETA(ETA(IPLOT)))
C
                     ISO(IPLOT) = ISOH
C
               IF (SHZ(IPLOT).LT.EPSILO)
     *             ISO(IPLOT)=IBSET(ISO(IPLOT),0)
               IF (SHZ(IPLOT).GT.1.D0-EPSILO)
     *             ISO(IPLOT)=IBSET(ISO(IPLOT),1)
C
                     GOTO 50
C
                  ELSE
C
C         WRITE(LU,*)'YA UN PROBLEME',IEL,IPLOT
C         WRITE(LU,*)'SHP',SHP1(IPLOT),SHP2(IPLOT),SHP3(IPLOT)
C         WRITE(LU,*)'SHZ',SHZ(IPLOT)
C         WRITE(LU,*)'DXYZ',DX(IPLOT),DY(IPLOT),DZ(IPLOT)
C         WRITE(LU,*)'XYZ',XPLOT(IPLOT),YPLOT(IPLOT),ZPLOT(IPLOT)

                  STOP
                  ENDIF
C
               ENDIF
C
            ENDIF
C
40       CONTINUE
C
20    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
