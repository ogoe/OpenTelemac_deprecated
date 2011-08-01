C                       *****************
                        SUBROUTINE PIED4D
C                       *****************
C
     *  (U , V , T , W , DT , NRK , X , Y , TETA , FREQ , IKLE2 , 
     *   IFABOR , ETAS , XPLOT , YPLOT , TPLOT , FPLOT , DX , DY , DW ,
     *   DF , SHP1 , SHP2 , SHP3 , SHT , SHF , ELT , ETA , FRE , NSP ,
     *   NPLOT , NPOIN2 , NELEM2 , NPLAN , NF , SURDET , SENS , 
     *   ISO )
C
C***********************************************************************
C  TOMAWAC  VERSION 1.0       01/02/95        F MARCOS (LNH) 30 87 72 66
C***********************************************************************
C
C  FONCTION :
C
C     REMONTEE OU DESCENTE
C     DES COURBES CARACTERISTIQUES
C     SUR LES HYPER PRISMES DE TOMAWAC
C     DANS L'INTERVALLE DE TEMPS DT
C     AVEC UNE DISCRETISATION HYBRIDE ELEMENTS FINIS+DIFF FINIS 2D
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
C !    U,V,T,W     ! -->! COMPOSANTE DE LA VITESSE DU CONVECTEUR       !
C !    DT          ! -->! PAS DE TEMPS.                                !
C !    NRK         ! -->! NOMBRE DE SOUS-PAS DE RUNGE-KUTTA.           !
C !  X,Y,TETA,FREQ ! -->! COORDONNEES DES POINTS DU MAILLAGE.          !
C !    IKLE2       ! -->! TRANSITION ENTRE LES NUMEROTATIONS LOCALE    !
C !                !    ! ET GLOBALE DU MAILLAGE 2D.                   !
C !    IFABOR      ! -->! NUMEROS 2D DES ELEMENTS AYANT UNE FACE COMMUNE
C !                !    ! AVEC L'ELEMENT .  SI IFABOR<0 OU NUL         !
C !                !    ! ON A UNE FACE LIQUIDE,SOLIDE,OU PERIODIQUE   !
C !    ETAS        !<-->! TABLEAU DE TRAVAIL DONNANT LE NUMERO DE      !
C !                !    ! L'ETAGE SUPERIEUR                            !
C ! X.,Y.,T.,FPLOT !<-->! POSITIONS SUCCESSIVES DES DERIVANTS.         !
C !  DX,DY,DW,DF   ! -- ! STOCKAGE DES SOUS-PAS .                      !
C !    SHP1-2-3    !<-->! COORDONNEES BARYCENTRIQUES 2D AU PIED DES    !
C !                !    ! COURBES CARACTERISTIQUES.                    !
C !    SHT         !<-->! COORDONNEES BARYCENTRIQUES SUIVANT TETA DES  !
C !                !    ! NOEUDS DANS LEURS ETAGES "ETA" ASSOCIES.     !
C !    SHF         !<-->! COORDONNEES BARYCENTRIQUES SUIVANT F DES     !
C !                !    ! NOEUDS DANS LEURS FREQUENCES "FRE" ASSOCIEES.!
C !    ELT         !<-->! NUMEROS DES ELEMENTS 2D CHOISIS POUR CHAQUE  !
C !                !    ! NOEUD.                                       !
C !    ETA         !<-->! NUMEROS DES ETAGES CHOISIS POUR CHAQUE NOEUD.!
C !    FRE         !<-->! NUMEROS DES FREQ. CHOISIES POUR CHAQUE NOEUD.!
C !    NSP         ! -- ! NOMBRE DE SOUS-PAS DE RUNGE KUTTA.           !
C !    NPLOT       ! -->! NOMBRE DE DERIVANTS.                         !
C !    NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE 2D.             !
C !    NELEM2      ! -->! NOMBRE D'ELEMENTS DU MAILLAGE 2D.            !
C !    NPLAN       ! -->! NOMBRE DE DIRECTIONS                         !
C !    NF          ! -->! NOMBRE DE FREQUENCES                         !
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
      INTEGER NPOIN2,NELEM2,NPLAN,NPLOT,NSPMAX,NRK,SENS,NF
C
      DOUBLE PRECISION U(NPOIN2,NPLAN,NF),V(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION T(NPOIN2,NPLAN,NF),W(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION XPLOT(NPLOT),YPLOT(NPLOT)
      DOUBLE PRECISION TPLOT(NPLOT),FPLOT(NPLOT)
      DOUBLE PRECISION SURDET(NELEM2),SHT(NPLOT),SHF(NPLOT)
      DOUBLE PRECISION SHP1(NPLOT),SHP2(NPLOT),SHP3(NPLOT)
      DOUBLE PRECISION X(NPOIN2),Y(NPOIN2),TETA(NPLAN+1),FREQ(NF)
      DOUBLE PRECISION DX(NPLOT),DY(NPLOT),DW(NPLOT),DF(NPLOT)
      DOUBLE PRECISION PAS,DT,EPSILO,A1,A2
      DOUBLE PRECISION DX1,DY1,DXP,DYP,DTP,DFP,XP,YP,TP,FP
C
      INTEGER IKLE2(NELEM2,3),IFABOR(NELEM2,7),ETAS(NPLAN)
      INTEGER ELT(NPLOT),ETA(NPLOT),FRE(NPLOT),NSP(NPLOT),ISO(NPLOT)
      INTEGER IPLOT,ISP,I1,I2,I3,IEL,IET,IFR,IFA,ISUI(3)
      INTEGER ISOH,ISOT,ISOF,ISOV
C
      INTRINSIC ABS , INT , MAX , SQRT
C
      DATA ISUI   / 2 , 3 , 1 /
      DATA EPSILO / -1.D-6 /
C
C-----------------------------------------------------------------------
C  CALCUL DU NOMBRE DE SOUS PAS MAXIMUM
C-----------------------------------------------------------------------
C
      NSPMAX = 1
C
      DO 10 IPLOT = 1 , NPLOT
C
         NSP(IPLOT) = 0
         IEL = ELT(IPLOT)
C
         IF (IEL.GT.0) THEN
C
            IET = ETA(IPLOT)
          IFR = FRE(IPLOT)
C
            I1 = IKLE2(IEL,1)
            I2 = IKLE2(IEL,2)
            I3 = IKLE2(IEL,3)
C
         DXP =(1.D0-SHF(IPLOT))*
     *              ( U(I1,IET  ,IFR)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *          + U(I2,IET  ,IFR)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *          + U(I3,IET  ,IFR)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *          + U(I1,ETAS(IET),IFR)*SHP1(IPLOT)*SHT(IPLOT)
     *          + U(I2,ETAS(IET),IFR)*SHP2(IPLOT)*SHT(IPLOT)
     *          + U(I3,ETAS(IET),IFR)*SHP3(IPLOT)*SHT(IPLOT))
     *        + SHF(IPLOT)*
     *              ( U(I1,IET  ,IFR+1)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *          + U(I2,IET  ,IFR+1)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *          + U(I3,IET  ,IFR+1)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *          + U(I1,ETAS(IET),IFR+1)*SHP1(IPLOT)*SHT(IPLOT)
     *          + U(I2,ETAS(IET),IFR+1)*SHP2(IPLOT)*SHT(IPLOT)
     *          + U(I3,ETAS(IET),IFR+1)*SHP3(IPLOT)*SHT(IPLOT))
C
         DYP =(1.D0-SHF(IPLOT))*
     *              ( V(I1,IET  ,IFR)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *          + V(I2,IET  ,IFR)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *          + V(I3,IET  ,IFR)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *          + V(I1,ETAS(IET),IFR)*SHP1(IPLOT)*SHT(IPLOT)
     *          + V(I2,ETAS(IET),IFR)*SHP2(IPLOT)*SHT(IPLOT)
     *          + V(I3,ETAS(IET),IFR)*SHP3(IPLOT)*SHT(IPLOT))
     *        + SHF(IPLOT)*
     *              ( V(I1,IET  ,IFR+1)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *          + V(I2,IET  ,IFR+1)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *          + V(I3,IET  ,IFR+1)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *          + V(I1,ETAS(IET),IFR+1)*SHP1(IPLOT)*SHT(IPLOT)
     *          + V(I2,ETAS(IET),IFR+1)*SHP2(IPLOT)*SHT(IPLOT)
     *          + V(I3,ETAS(IET),IFR+1)*SHP3(IPLOT)*SHT(IPLOT))
C
         DTP =(1.D0-SHF(IPLOT))*
     *              ( T(I1,IET  ,IFR)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *          + T(I2,IET  ,IFR)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *          + T(I3,IET  ,IFR)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *          + T(I1,ETAS(IET),IFR)*SHP1(IPLOT)*SHT(IPLOT)
     *          + T(I2,ETAS(IET),IFR)*SHP2(IPLOT)*SHT(IPLOT)
     *          + T(I3,ETAS(IET),IFR)*SHP3(IPLOT)*SHT(IPLOT))
     *        + SHF(IPLOT)*
     *              ( T(I1,IET  ,IFR+1)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *          + T(I2,IET  ,IFR+1)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *          + T(I3,IET  ,IFR+1)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *          + T(I1,ETAS(IET),IFR+1)*SHP1(IPLOT)*SHT(IPLOT)
     *          + T(I2,ETAS(IET),IFR+1)*SHP2(IPLOT)*SHT(IPLOT)
     *          + T(I3,ETAS(IET),IFR+1)*SHP3(IPLOT)*SHT(IPLOT))
C
         DFP =(1.D0-SHF(IPLOT))*
     *              ( W(I1,IET  ,IFR)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *          + W(I2,IET  ,IFR)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *          + W(I3,IET  ,IFR)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *          + W(I1,ETAS(IET),IFR)*SHP1(IPLOT)*SHT(IPLOT)
     *          + W(I2,ETAS(IET),IFR)*SHP2(IPLOT)*SHT(IPLOT)
     *          + W(I3,ETAS(IET),IFR)*SHP3(IPLOT)*SHT(IPLOT))
     *        + SHF(IPLOT)*
     *              ( W(I1,IET  ,IFR+1)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *          + W(I2,IET  ,IFR+1)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *          + W(I3,IET  ,IFR+1)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *          + W(I1,ETAS(IET),IFR+1)*SHP1(IPLOT)*SHT(IPLOT)
     *          + W(I2,ETAS(IET),IFR+1)*SHP2(IPLOT)*SHT(IPLOT)
     *          + W(I3,ETAS(IET),IFR+1)*SHP3(IPLOT)*SHT(IPLOT))
C
         NSP(IPLOT)= MAX( INT(NRK*DT*ABS(DTP/(TETA(IET)-TETA(IET+1)))),
     *                   INT(NRK*DT*ABS(DFP/(FREQ(IFR)-FREQ(IFR+1)))) )
         NSP(IPLOT)= MAX( NSP(IPLOT),
     *               INT(NRK*DT*SQRT((DXP*DXP+DYP*DYP)*SURDET(IEL))) )
C
          NSP(IPLOT) = MAX (1,NSP(IPLOT))
C
            NSPMAX = MAX ( NSPMAX , NSP(IPLOT) )
C
         ENDIF 
C
10    CONTINUE
      IF(LNG.EQ.1) THEN
         WRITE(LU,*) 'NOMBRE MAX DE SOUS PAS :',NSPMAX
      ELSE
         WRITE(LU,*) 'NUMBER OF SUB-ITERATIONS :',NSPMAX
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
               IFR = FRE(IPLOT)
C
               I1 = IKLE2(IEL,1)
               I2 = IKLE2(IEL,2)
               I3 = IKLE2(IEL,3)
               PAS = SENS * DT / NSP(IPLOT)
C
         DX(IPLOT) = ( (1.D0-SHF(IPLOT))*
     *          ( U(I1,IET  ,IFR)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *      + U(I2,IET  ,IFR)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *      + U(I3,IET  ,IFR)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *      + U(I1,ETAS(IET),IFR)*SHP1(IPLOT)*SHT(IPLOT)
     *      + U(I2,ETAS(IET),IFR)*SHP2(IPLOT)*SHT(IPLOT)
     *      + U(I3,ETAS(IET),IFR)*SHP3(IPLOT)*SHT(IPLOT))
     *        + SHF(IPLOT)*
     *          ( U(I1,IET  ,IFR+1)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *      + U(I2,IET  ,IFR+1)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *      + U(I3,IET  ,IFR+1)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *      + U(I1,ETAS(IET),IFR+1)*SHP1(IPLOT)*SHT(IPLOT)
     *      + U(I2,ETAS(IET),IFR+1)*SHP2(IPLOT)*SHT(IPLOT)
     *      + U(I3,ETAS(IET),IFR+1)*SHP3(IPLOT)*SHT(IPLOT)) )*PAS
C
         DY(IPLOT) = ( (1.D0-SHF(IPLOT))*
     *          ( V(I1,IET  ,IFR)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *      + V(I2,IET  ,IFR)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *      + V(I3,IET  ,IFR)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *      + V(I1,ETAS(IET),IFR)*SHP1(IPLOT)*SHT(IPLOT)
     *      + V(I2,ETAS(IET),IFR)*SHP2(IPLOT)*SHT(IPLOT)
     *      + V(I3,ETAS(IET),IFR)*SHP3(IPLOT)*SHT(IPLOT))
     *        + SHF(IPLOT)*
     *          ( V(I1,IET  ,IFR+1)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *      + V(I2,IET  ,IFR+1)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *      + V(I3,IET  ,IFR+1)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *      + V(I1,ETAS(IET),IFR+1)*SHP1(IPLOT)*SHT(IPLOT)
     *      + V(I2,ETAS(IET),IFR+1)*SHP2(IPLOT)*SHT(IPLOT)
     *      + V(I3,ETAS(IET),IFR+1)*SHP3(IPLOT)*SHT(IPLOT)) )*PAS
C
         DW(IPLOT) = ( (1.D0-SHF(IPLOT))*
     *          ( T(I1,IET  ,IFR)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *      + T(I2,IET  ,IFR)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *      + T(I3,IET  ,IFR)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *      + T(I1,ETAS(IET),IFR)*SHP1(IPLOT)*SHT(IPLOT)
     *      + T(I2,ETAS(IET),IFR)*SHP2(IPLOT)*SHT(IPLOT)
     *      + T(I3,ETAS(IET),IFR)*SHP3(IPLOT)*SHT(IPLOT))
     *        + SHF(IPLOT)*
     *          ( T(I1,IET  ,IFR+1)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *      + T(I2,IET  ,IFR+1)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *      + T(I3,IET  ,IFR+1)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *      + T(I1,ETAS(IET),IFR+1)*SHP1(IPLOT)*SHT(IPLOT)
     *      + T(I2,ETAS(IET),IFR+1)*SHP2(IPLOT)*SHT(IPLOT)
     *      + T(I3,ETAS(IET),IFR+1)*SHP3(IPLOT)*SHT(IPLOT)) )*PAS
C
         DF(IPLOT) = ( (1.D0-SHF(IPLOT))*
     *          ( W(I1,IET  ,IFR)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *      + W(I2,IET  ,IFR)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *      + W(I3,IET  ,IFR)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *      + W(I1,ETAS(IET),IFR)*SHP1(IPLOT)*SHT(IPLOT)
     *      + W(I2,ETAS(IET),IFR)*SHP2(IPLOT)*SHT(IPLOT)
     *      + W(I3,ETAS(IET),IFR)*SHP3(IPLOT)*SHT(IPLOT))
     *        + SHF(IPLOT)*
     *          ( W(I1,IET  ,IFR+1)*SHP1(IPLOT)*(1.D0-SHT(IPLOT))
     *      + W(I2,IET  ,IFR+1)*SHP2(IPLOT)*(1.D0-SHT(IPLOT))
     *      + W(I3,IET  ,IFR+1)*SHP3(IPLOT)*(1.D0-SHT(IPLOT))
     *      + W(I1,ETAS(IET),IFR+1)*SHP1(IPLOT)*SHT(IPLOT)
     *      + W(I2,ETAS(IET),IFR+1)*SHP2(IPLOT)*SHT(IPLOT)
     *      + W(I3,ETAS(IET),IFR+1)*SHP3(IPLOT)*SHT(IPLOT)) )*PAS
C
               XP = XPLOT(IPLOT) + DX(IPLOT)
               YP = YPLOT(IPLOT) + DY(IPLOT)
               TP = TPLOT(IPLOT) + DW(IPLOT)
               FP = FPLOT(IPLOT) + DF(IPLOT)
C
               SHP1(IPLOT) = ((X(I3)-X(I2))*(YP-Y(I2))
     *                        -(Y(I3)-Y(I2))*(XP-X(I2))) * SURDET(IEL)
               SHP2(IPLOT) = ((X(I1)-X(I3))*(YP-Y(I3))
     *                        -(Y(I1)-Y(I3))*(XP-X(I3))) * SURDET(IEL)
               SHP3(IPLOT) = ((X(I2)-X(I1))*(YP-Y(I1))
     *                        -(Y(I2)-Y(I1))*(XP-X(I1))) * SURDET(IEL)
               SHT(IPLOT) = (TP-TETA(IET)) / (TETA(IET+1)-TETA(IET))
               SHF(IPLOT) = (FP-FREQ(IFR)) / (FREQ(IFR+1)-FREQ(IFR))
C             IF (ABS(SHT(IPLOT)).GT.2.5D0 ) THEN
C      	  WRITE(LU,*) 'SHT***',IPLOT,IET,SHT(IPLOT)
C      	  WRITE(LU,*) TETA(IET),TETA(IET+1),ZP
C      	  WRITE(LU,*) DZ(IPLOT),ZPLOT(IPLOT)
C      	  STOP
C              ENDIF
C
               XPLOT(IPLOT) = XP
               YPLOT(IPLOT) = YP
               TPLOT(IPLOT) = TP
               FPLOT(IPLOT) = FP
C
               IF (SHP1(IPLOT).LT.EPSILO)
     *              ISO(IPLOT)=IBSET(ISO(IPLOT),4)
               IF (SHP2(IPLOT).LT.EPSILO)
     *              ISO(IPLOT)=IBSET(ISO(IPLOT),5)
               IF (SHP3(IPLOT).LT.EPSILO)
     *              ISO(IPLOT)=IBSET(ISO(IPLOT),6)
C
               IF  (SHT(IPLOT).LT.EPSILO)
     *              ISO(IPLOT)=IBSET(ISO(IPLOT),0)
               IF  (SHT(IPLOT).GT.1.D0-EPSILO)
     *              ISO(IPLOT)=IBSET(ISO(IPLOT),1)
C
               IF  (SHF(IPLOT).LT.EPSILO)
     *              ISO(IPLOT)=IBSET(ISO(IPLOT),2)
               IF  (SHF(IPLOT).GT.1.D0-EPSILO)
     *              ISO(IPLOT)=IBSET(ISO(IPLOT),3)
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
              ISOT = IAND(ISO(IPLOT), 3)
              ISOF = IAND(ISO(IPLOT),12)/4
            ISOV = IAND(ISO(IPLOT),15)
              ISOH = IAND(ISO(IPLOT),112)
              IEL = ELT(IPLOT)
              IET = ETA(IPLOT)
              IFR = FRE(IPLOT)
              XP = XPLOT(IPLOT)
              YP = YPLOT(IPLOT)
              TP = TPLOT(IPLOT)
              FP = FPLOT(IPLOT)
C
              IF (ISOH.NE.0) THEN
C
                IF (ISOH.EQ.16) THEN
                   IFA = 2
                ELSEIF (ISOH.EQ.32) THEN
                   IFA = 3
                ELSEIF (ISOH.EQ.64) THEN
                   IFA = 1
                ELSEIF (ISOH.EQ.48) THEN
                   IFA = 2
                   IF (DX(IPLOT)*(Y(IKLE2(IEL,3))-YP).LT.
     *                 DY(IPLOT)*(X(IKLE2(IEL,3))-XP)) IFA = 3
                ELSEIF (ISOH.EQ.96) THEN
                   IFA = 3
                   IF (DX(IPLOT)*(Y(IKLE2(IEL,1))-YP).LT.
     *                 DY(IPLOT)*(X(IKLE2(IEL,1))-XP)) IFA = 1
                ELSE
                   IFA = 1
                   IF (DX(IPLOT)*(Y(IKLE2(IEL,2))-YP).LT.
     *                 DY(IPLOT)*(X(IKLE2(IEL,2))-XP)) IFA = 2
                ENDIF
C
                IF (ISOV.GT.0) THEN
                  I1 = IKLE2(IEL,IFA)
                  I2 = IKLE2(IEL,ISUI(IFA))
                  IF (ISOF.GT.0) THEN
      	      IF (ISOT.GT.0) THEN
      	        A1=(FP-FREQ(IFR+ISOF-1))/DF(IPLOT)
      	        A2=(TP-TETA(IET+ISOT-1))/DW(IPLOT)
      	        IF (A1.LT.A2) THEN
                          IF ((X(I2)-X(I1))*(YP-A1*DY(IPLOT)-Y(I1)).GT.
     *             (Y(I2)-Y(I1))*(XP-A1*DX(IPLOT)-X(I1))) IFA=ISOF+5
      	        ELSE
                          IF ((X(I2)-X(I1))*(YP-A2*DY(IPLOT)-Y(I1)).GT.
     *             (Y(I2)-Y(I1))*(XP-A2*DX(IPLOT)-X(I1))) IFA=ISOT+3
      		ENDIF
      	      ELSE
                        A1 = (FP-FREQ(IFR+ISOF-1)) / DF(IPLOT)
                        IF ((X(I2)-X(I1))*(YP-A1*DY(IPLOT)-Y(I1)).GT.
     *             (Y(I2)-Y(I1))*(XP-A1*DX(IPLOT)-X(I1))) IFA=ISOF+5
      	       ENDIF
      	    ELSE
                      A1 = (TP-TETA(IET+ISOT-1)) / DW(IPLOT)
                      IF ((X(I2)-X(I1))*(YP-A1*DY(IPLOT)-Y(I1)).GT.
     *             (Y(I2)-Y(I1))*(XP-A1*DX(IPLOT)-X(I1))) IFA=ISOT+3
      	    ENDIF
                ENDIF
C
               ELSEIF (ISOT.GT.0) THEN
C
      	  IFA = ISOT + 3
C
                  IF (ISOF.GT.0) THEN
      	    A1=(FP-FREQ(IFR+ISOF-1))/DF(IPLOT)
      	    A2=(TP-TETA(IET+ISOT-1))/DW(IPLOT)
      	    IF (A1.LT.A2) IFA = ISOF + 5
      	  ENDIF
               ELSE
      	  IFA = ISOF + 5
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
         IF (SHP1(IPLOT).LT.EPSILO) ISO(IPLOT)=IBSET(ISO(IPLOT),4)
         IF (SHP2(IPLOT).LT.EPSILO) ISO(IPLOT)=IBSET(ISO(IPLOT),5)
         IF (SHP3(IPLOT).LT.EPSILO) ISO(IPLOT)=IBSET(ISO(IPLOT),6)
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
                  IF (A1.GT.1.D0) A1 = 1.D0
                  IF (A1.LT.0.D0) A1 = 0.D0
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
                  TPLOT(IPLOT) = TP - A1*DW(IPLOT)
                  SHT(IPLOT) = (TPLOT(IPLOT)-TETA(IET))
     *                       / (TETA(IET+1)-TETA(IET))
                  FPLOT(IPLOT) = FP - A1*DF(IPLOT)
                  SHF(IPLOT) = (FPLOT(IPLOT)-FREQ(IFR))
     *                       / (FREQ(IFR+1)-FREQ(IFR))
                  ELT(IPLOT) = - SENS * ELT(IPLOT)
                  NSP(IPLOT) = ISP
C
               ELSEIF (IFA.LE.5) THEN
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE DU PRISME EST UNE FACE TRIAN. TETA
C  =====================================================================
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
      		 TP=TP-2*3.14159265D0
      		 TPLOT(IPLOT)=TP
                     ENDIF
      	     IF (ETA(IPLOT).EQ.0) THEN
      		 ETA(IPLOT) = NPLAN
      		 TP=TP+2*3.14159265D0
      		 TPLOT(IPLOT)=TP
                     ENDIF
                     SHT(IPLOT) = (TP-TETA(ETA(IPLOT)))
     *                   / (TETA(ETA(IPLOT)+1)-TETA(ETA(IPLOT)))
C
                     ISO(IPLOT) = ISOH+ISOF*4
C
               IF (SHT(IPLOT).LT.EPSILO)
     *             ISO(IPLOT)=IBSET(ISO(IPLOT),0)
               IF (SHT(IPLOT).GT.1.D0-EPSILO)
     *             ISO(IPLOT)=IBSET(ISO(IPLOT),1)
C
                     GOTO 50
C
                  ELSE
C
        IF(LNG.EQ.1) THEN
         WRITE(LU,*) 'PROBLEME DANS PIED4D',IEL,IPLOT
        ELSE
         WRITE(LU,*) 'PROBLEM IN PIED4D',IEL,IPLOT
        ENDIF
        WRITE(LU,*) 'SHP',SHP1(IPLOT),SHP2(IPLOT),SHP3(IPLOT)
        WRITE(LU,*) 'SHT',SHT(IPLOT)
        WRITE(LU,*) 'DXYZ',DX(IPLOT),DY(IPLOT),DW(IPLOT)
        WRITE(LU,*) 'XYZ',XPLOT(IPLOT),YPLOT(IPLOT),TPLOT(IPLOT)

        STOP
                  ENDIF
C
               ELSE
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE DU PRISME EST UNE FACE TRIAN. FREQ
C  =====================================================================
C-----------------------------------------------------------------------
C
                  IFA = IFA - 6
C
                  IF ((IFA.EQ.1).AND.(IFR.EQ.NF-1)) IEL=-1
                  IF ((IFA.EQ.0).AND.(IFR.EQ.1)) IEL=-1
                  IF (IEL.EQ.1) THEN
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE EST INTERNE AU DOMAINE
C  ON SE RELOCALISE DANS L'ELEMENT ADJACENT
C-----------------------------------------------------------------------
C
                     FRE(IPLOT) = IFR + IFA + IFA - 1
                     SHF(IPLOT) = (FP-FREQ(FRE(IPLOT)))
     *                   / (FREQ(FRE(IPLOT)+1)-FREQ(FRE(IPLOT)))
C
                     ISO(IPLOT) = ISOH+ISOT
C
               IF (SHF(IPLOT).LT.EPSILO)
     *             ISO(IPLOT)=IBSET(ISO(IPLOT),2)
               IF (SHF(IPLOT).GT.1.D0-EPSILO)
     *             ISO(IPLOT)=IBSET(ISO(IPLOT),3)
C
                     GOTO 50
C
                  ELSE
C
C-----------------------------------------------------------------------
C  LA, ON SAIT QUE LA FACE DE SORTIE EST LA FREQUENCE MIN OU MAX
C  ON PROJETE LE RELICAT SUR LA FRONTIERE ET ON CONTINUE
C-----------------------------------------------------------------------
C
                    FPLOT(IPLOT)=FREQ(IFR+IFA)
      	    DF(IPLOT)=0.D0
      	    SHF(IPLOT)=IFA
      	    ISO(IPLOT) = ISOH +ISOT
      	    IF(ISO(IPLOT).NE.0) GOTO 50
C
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
