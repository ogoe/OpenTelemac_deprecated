C                       *****************
                        SUBROUTINE SPEINI
C                       *****************
C
     *( F     , SPEC  , FRA   , UV    , VV    , FREQ  , TETA  , GRAVIT,
     *  FREMAX, FETCH , SIGMAA, SIGMAB, GAMMA , FPIC  , HM0   , ALPHIL,
     *  TETA1 , SPRED1, TETA2 , SPRED2, XLAMDA, NPOIN2, NPLAN , NF    ,
     *  INISPE, E2FMIN, DEPTH , FRABI )
C
C**********************************************************************
C  TOMAWAC - V1.0    M. BENOIT               (EDF/DER/LNH)  -  13/07/95
C**********************************************************************
C
C  FONCTION : AFFECTATION DU SPECTRE DE VARIANCE SUR LE DOMAINE.
C  ********** PLUSIEURS OPTIONS SONT PROPOSEES SELON LA VALEUR
C             DE INISPE :
C               0. SPECTRE NUL EN TOUS LES POINTS.
C               1. SPECTRE DE TYPE JONSWAP FONCTION DU VENT ET NUL SI
C                  LA VITESSE DU VENT EST NULLE.
C               2. SPECTRE DE TYPE JONSWAP FONCTION DU VENT ET
C                  PARAMETRIQUE SI LA VITESSE DU VENT EST NULLE.
C               3. SPECTRE DE TYPE JONSWAP PARAMETRIQUE
C
C
C  ARGUMENTS :
C  ***********
C  +-------------+----+--------------------------------------------+
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C  +-------------+----+--------------------------------------------+
C  ! F(-,-,-)    !<-- ! SPECTRE DIRECTIONNEL DE VARIANCE           !
C  ! SPEC(-)     !<-- ! VECTEUR DU SPECTRE EN FREQUENCE            !
C  ! FRA(-)      !<-- ! VECTEUR DE LA FONCTION DE REPARTITION ANG. !
C  ! UV(-)       ! -->! TABLEAU DES COMPOSANTES OUEST-EST DU VENT  !
C  ! VV(-)       ! -->! TABLEAU DES COMPOSANTES SUD-NORD  DU VENT  !
C  ! FREQ(-)     ! -->! VECTEUR DES FREQUENCES DE DISCRETISATION   !
C  ! TETA(-)     ! -->! VECTEUR DES DIRECTIONS DE DISCRETISATION   !
C  ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !
C  ! FREMAX      ! -->! VALEUR MAXIMUM DE LA FREQUENCE DE PIC      !
C  ! FETCH       ! -->! FETCH MOYEN                                !
C  ! SIGMAA      ! -->! VALEUR DE SIGMA JONSWAP POUR F < FP        !
C  ! SIGMAB      ! -->! VALEUR DE SIGMA JONSWAP POUR F > FP        !
C  ! GAMMA       ! -->! FACTEUR DE FORME DE PIC JONSWAP            !
C  ! FPIC        ! -->! FREQUENCE DE PIC JONSWAP                   !
C  ! HM0         ! -->! HAUTEUR SIGNIFICATIVE JONSWAP              !
C  ! ALPHIL      ! -->! CONSTANTE DE PHILLIPS (ALPHA)              !
C  ! TETA1       ! -->! DIRECTION PRINCIPALE 1 POUR FRA            !
C  ! SPRED1      ! -->! ETALEMENT DIRECTIONNEL 1 POUR FRA          !
C  ! TETA2       ! -->! DIRECTION PRINCIPALE 2 POUR FRA            !
C  ! SPRED2      ! -->! ETALEMENT DIRECTIONNEL 2 POUR FRA          !
C  ! XLAMDA      ! -->! FACTEUR DE PONDERATION POUR LA FRA         !
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIOSN DE DISCRETISATION     !
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !
C  ! INISPE      ! -->! INDICATEUR D'INITIALISATION DU SPECTRE     !
C  ! E2FMIN      ! -->! SEUIL MINIMUM DE VARIANCE CONSIDERE        !
C  +-------------+----+--------------------------------------------+
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C  +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  CONDIW
C  ********    - PROGRAMME(S) APPELE(S) :  -      
C
C  REMARQUES :
C  ***********
C  - 
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NPOIN2, NPLAN , NF    , INISPE, FRABI
      DOUBLE PRECISION GRAVIT, FREMAX, FETCH , SIGMAA, SIGMAB, GAMMA
      DOUBLE PRECISION FPIC  , HM0   , ALPHIL, TETA1 , SPRED1, TETA2
      DOUBLE PRECISION SPRED2, XLAMDA, E2FMIN
      DOUBLE PRECISION FREQ(NF) , TETA(NPLAN), SPEC(NF), FRA(NPLAN)
      DOUBLE PRECISION UV(NPOIN2) , VV(NPOIN2), DEPTH(NPOIN2)
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF)
C
C.....VARIABLES LOCALES 
C     """""""""""""""""
      INTEGER  NPOIN4, IP    , JF    , JP
      DOUBLE PRECISION GX    , GXU   , UG     , AL    , FP    , DEUPI
      DOUBLE PRECISION UVMIN  , COEFA , COEFB , COEFD
      DOUBLE PRECISION COEFE , UVENT , FPMIN  , SPR1  , SPR2  , XLAM
      DOUBLE PRECISION TET1  , TET2  , COEF
C
C
      DEUPI = 6.283185307D0
      NPOIN4= NPOIN2*NPLAN*NF
      UVMIN = 1.D-6
      COEFA = 2.84D0
      COEFB = 0.033D0
      COEFD =-3.D0/10.D0
      COEFE = 2.D0/3.D0
      GX    = GRAVIT*FETCH
      FPMIN = 1.D-4
C
C
C     ===========================================================
C     SPECTRE INITIAL NUL EN TOUS POINTS (INISPE=0)
C     (FAIT EGALEMENT POUR INITIALISATION POUR LES AUTRES VALEURS)
C     ===========================================================
      CALL OV ( 'X=C     ' , F      , F , F , 0.D0 , NPOIN4 )
      IF (INISPE.EQ.0) RETURN
C
C
C     ==/ INISPE = 1 /===========================================
C     SI VENT NON NUL -E(F): JONSWAP FONCTION DU VENT (AL,FP)
C                     -FRA : UNIMODALE AUTOUR DE TETA(VENT)
C     SI VENT NUL     -E(F): NUL
C                     -FRA : NUL
C     ===========================================================
      IF (INISPE.EQ.1) THEN
C
        DO 100 IP=1,NPOIN2
          UVENT=SQRT(UV(IP)**2+VV(IP)**2)
          IF (UVENT.GT.UVMIN) THEN
C
C...........CALCUL DU SPECTRE EN FREQUENCE (JONSWAP).
C           """""""""""""""""""""""""""""""""""""""""
            GXU=GX/(UVENT*UVENT)
            UG = UVENT/GRAVIT
            FP = MAX(0.13D0,COEFA*GXU**COEFD)
            FP = MIN(FP,FREMAX*UG)
            AL = MAX(0.0081D0, COEFB*FP**COEFE)
            FP = FP/UG
            CALL SPEJON
     &( SPEC  , FREQ  , NF    , AL    , FP     , GAMMA , SIGMAA, SIGMAB,
     &  DEUPI , GRAVIT, E2FMIN, FPMIN )
C
C...........CALCUL DE LA FONCTION DE REPARTITION ANGULAIRE.
C           """""""""""""""""""""""""""""""""""""""""""""""
            SPR1=SPRED1
            TET1=ATAN2(UV(IP),VV(IP))
            SPR2=1.D0
            TET2=0.D0
            XLAM=1.D0
            IF(FRABI.EQ.2) THEN
              CALL FSPRD2
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
            ELSEIF(FRABI.EQ.3) THEN
              CALL FSPRD3
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
            ELSE
              CALL FSPRD1
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
            ENDIF
C
C...........CALCUL DU SPECTRE DIRECTIONNEL.
C           """""""""""""""""""""""""""""""
            DO 140 JF=1,NF
              DO 150 JP=1,NPLAN
                F(IP,JP,JF)=SPEC(JF)*FRA(JP)
  150         CONTINUE
  140       CONTINUE
          ENDIF
C
  100   CONTINUE
C
C     ==/ INISPE = 2 /===========================================
C     SI VENT NON NUL -E(F): JONSWAP FONCTION DU VENT (AL,FP)
C                     -FRA : UNIMODALE AUTOUR DE TETA(VENT)
C     SI VENT NUL     -E(F): JONSWAP PARAMETRE (AL,FP)
C                     -FRA : UNIMODALE PARAMETREE
C     ===========================================================
      ELSEIF (INISPE.EQ.2) THEN
C
        DO 200 IP=1,NPOIN2
          UVENT=SQRT(UV(IP)**2+VV(IP)**2)
C
C.........CALCUL DU SPECTRE EN FREQUENCE (JONSWAP).
C         """""""""""""""""""""""""""""""""""""""""
          IF (UVENT.GT.UVMIN) THEN
            GXU=GX/(UVENT*UVENT)
            UG = UVENT/GRAVIT
            FP = MAX(0.13D0,COEFA*GXU**COEFD)
            FP = MIN(FP,FREMAX*UG)
            AL = MAX(0.0081D0, COEFB*FP**COEFE)
            FP = FP/UG
          ELSE
            AL=ALPHIL
            FP=FPIC
          ENDIF
          CALL SPEJON
     &( SPEC  , FREQ  , NF    , AL    , FP     , GAMMA , SIGMAA, SIGMAB,
     &  DEUPI , GRAVIT, E2FMIN, FPMIN )
C
C.........CALCUL DE LA FONCTION DE REPARTITION ANGULAIRE.
C         """""""""""""""""""""""""""""""""""""""""""""""
          IF (UVENT.GT.UVMIN) THEN
            TET1=ATAN2(UV(IP),VV(IP))
          ELSE
            TET1=TETA1
          ENDIF
          SPR1=SPRED1
          SPR2=1.D0
          TET2=0.D0
          XLAM=1.D0
          IF(FRABI.EQ.2) THEN
            CALL FSPRD2
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
          ELSEIF(FRABI.EQ.3) THEN
            CALL FSPRD3
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
          ELSE
            CALL FSPRD1
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
          ENDIF
C
C.........CALCUL DU SPECTRE DIRECTIONNEL.
C         """""""""""""""""""""""""""""""
          DO 240 JF=1,NF
            DO 250 JP=1,NPLAN
              F(IP,JP,JF)=SPEC(JF)*FRA(JP)
  250       CONTINUE
  240     CONTINUE
C
  200   CONTINUE
C
C     ==/ INISPE = 3 /===========================================
C     SI VENT NON NUL -E(F): JONSWAP PARAMETRE (AL,FP)
C                     -FRA : UNIMODALE AUTOUR DE TETA(VENT)
C     SI VENT NUL     -E(F): NUL
C                     -FRA : NUL
C     ===========================================================
      ELSEIF (INISPE.EQ.3) THEN
C
        DO 300 IP=1,NPOIN2
          UVENT=SQRT(UV(IP)**2+VV(IP)**2)
          IF (UVENT.GT.UVMIN) THEN
C
C...........CALCUL DU SPECTRE EN FREQUENCE (JONSWAP).
C           """""""""""""""""""""""""""""""""""""""""
            AL = ALPHIL
            FP = FPIC
            CALL SPEJON
     &( SPEC  , FREQ  , NF    , AL    , FP     , GAMMA , SIGMAA, SIGMAB,
     &  DEUPI , GRAVIT, E2FMIN, FPMIN )
C
C...........CALCUL DE LA FONCTION DE REPARTITION ANGULAIRE.
C           """""""""""""""""""""""""""""""""""""""""""""""
            SPR1=SPRED1
            TET1=ATAN2(UV(IP),VV(IP))
            SPR2=1.D0
            TET2=0.D0
            XLAM=1.D0
            IF(FRABI.EQ.2) THEN
              CALL FSPRD2
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
            ELSEIF(FRABI.EQ.3) THEN
              CALL FSPRD3
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
            ELSE
              CALL FSPRD1
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
            ENDIF
C
C...........CALCUL DU SPECTRE DIRECTIONNEL.
C           """""""""""""""""""""""""""""""
            DO 340 JF=1,NF
              DO 350 JP=1,NPLAN
                F(IP,JP,JF)=SPEC(JF)*FRA(JP)
  350         CONTINUE
  340       CONTINUE
          ENDIF
C
  300   CONTINUE
C
C     ==/ INISPE = 4 /===========================================
C     SI VENT NON NUL -E(F): JONSWAP PARAMETRE (AL,FP)
C                     -FRA : UNIMODALE PARAMETREE
C     SI VENT NUL     -E(F): JONSWAP PARAMETRE (AL,FP)
C                     -FRA : UNIMODALE PARAMETREE
C     ===========================================================
      ELSEIF (INISPE.EQ.4) THEN
C
        DO 400 IP=1,NPOIN2
C
C.........CALCUL DU SPECTRE EN FREQUENCE (JONSWAP).
C         """""""""""""""""""""""""""""""""""""""""
          AL = ALPHIL
          FP = FPIC
          CALL SPEJON
     &( SPEC  , FREQ  , NF    , AL    , FP     , GAMMA , SIGMAA, SIGMAB,
     &  DEUPI , GRAVIT, E2FMIN, FPMIN )
C
C.........CALCUL DE LA FONCTION DE REPARTITION ANGULAIRE.
C         """""""""""""""""""""""""""""""""""""""""""""""
          SPR1=SPRED1
          TET1=TETA1
          SPR2=SPRED2
          TET2=TETA2
          XLAM=XLAMDA
          IF(FRABI.EQ.2) THEN
            CALL FSPRD2
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
          ELSEIF(FRABI.EQ.3) THEN
            CALL FSPRD3
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
          ELSE
            CALL FSPRD1
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
          ENDIF
C
C.........CALCUL DU SPECTRE DIRECTIONNEL.
C         """""""""""""""""""""""""""""""
          DO 440 JF=1,NF
            DO 450 JP=1,NPLAN
              F(IP,JP,JF)=SPEC(JF)*FRA(JP)
  450       CONTINUE
  440     CONTINUE
C
  400   CONTINUE
C
C     ==/ INISPE = 5 /===========================================
C     SI VENT NON NUL -E(F): JONSWAP PARAMETRE (HM0,FP)
C                     -FRA : UNIMODALE AUTOUR DE TETA(VENT)
C     SI VENT NUL     -E(F): NUL
C                     -FRA : NUL
C     ===========================================================
      ELSEIF (INISPE.EQ.5) THEN
C
        COEF=0.0624D0/(0.230D0+0.0336D0*GAMMA-0.185D0/(1.9D0+GAMMA))
     &      *(DEUPI*FPIC)**4*HM0*HM0/GRAVIT**2
C
        DO 500 IP=1,NPOIN2
          UVENT=SQRT(UV(IP)**2+VV(IP)**2)
          IF (UVENT.GT.UVMIN) THEN
C
C...........CALCUL DU SPECTRE EN FREQUENCE (JONSWAP).
C           """""""""""""""""""""""""""""""""""""""""
            AL=COEF
            FP = FPIC
            CALL SPEJON
     &( SPEC  , FREQ  , NF    , AL    , FP     , GAMMA , SIGMAA, SIGMAB,
     &  DEUPI , GRAVIT, E2FMIN, FPMIN )
C
C...........CALCUL DE LA FONCTION DE REPARTITION ANGULAIRE.
C           """""""""""""""""""""""""""""""""""""""""""""""
            SPR1=SPRED1
            TET1=ATAN2(UV(IP),VV(IP))
            SPR2=1.D0
            TET2=0.D0
            XLAM=1.D0
            IF(FRABI.EQ.2) THEN
              CALL FSPRD2
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
            ELSEIF(FRABI.EQ.3) THEN
              CALL FSPRD3
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
            ELSE
              CALL FSPRD1
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
            ENDIF
C
C...........CALCUL DU SPECTRE DIRECTIONNEL.
C           """""""""""""""""""""""""""""""
            DO 540 JF=1,NF
              DO 550 JP=1,NPLAN
                F(IP,JP,JF)=SPEC(JF)*FRA(JP)
  550         CONTINUE
  540       CONTINUE
          ENDIF
C
  500   CONTINUE
C
C     ==/ INISPE = 6 /===========================================
C     SI VENT NON NUL -E(F): JONSWAP PARAMETRE (HM0,FP)
C                     -FRA : UNIMODALE PARAMETREE
C     SI VENT NUL     -E(F): JONSWAP PARAMETRE (HM0,FP)
C                     -FRA : UNIMODALE PARAMETREE
C     ===========================================================
      ELSEIF (INISPE.EQ.6) THEN
C
        COEF=0.0624D0/(0.230D0+0.0336D0*GAMMA-0.185D0/(1.9D0+GAMMA))
     &      *(DEUPI*FPIC)**4*HM0*HM0/GRAVIT**2
C
        DO 600 IP=1,NPOIN2
C
C.........CALCUL DU SPECTRE EN FREQUENCE (JONSWAP).
C         """""""""""""""""""""""""""""""""""""""""
          AL = COEF
          FP = FPIC
          CALL SPEJON
     &( SPEC  , FREQ  , NF    , AL    , FP     , GAMMA , SIGMAA, SIGMAB,
     &  DEUPI , GRAVIT, E2FMIN, FPMIN )
C
C.........CALCUL DE LA FONCTION DE REPARTITION ANGULAIRE.
C         """""""""""""""""""""""""""""""""""""""""""""""
          SPR1=SPRED1
          TET1=TETA1
          SPR2=SPRED2
          TET2=TETA2
          XLAM=XLAMDA
          IF(FRABI.EQ.2) THEN
            CALL FSPRD2
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
          ELSEIF(FRABI.EQ.3) THEN
            CALL FSPRD3
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
          ELSE
            CALL FSPRD1
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
          ENDIF
C
C.........CALCUL DU SPECTRE DIRECTIONNEL.
C         """""""""""""""""""""""""""""""
          DO 640 JF=1,NF
            DO 650 JP=1,NPLAN
              F(IP,JP,JF)=SPEC(JF)*FRA(JP)
  650       CONTINUE
  640     CONTINUE
C
  600   CONTINUE
C
C     ==/ INISPE = 7 /===========================================
C     SI VENT NON NUL -E(F): TMA PARAMETRE (HM0,FP)
C                     -FRA : UNIMODALE PARAMETREE
C     SI VENT NUL     -E(F): TMA PARAMETRE (HM0,FP)
C                     -FRA : UNIMODALE PARAMETREE
C     ===========================================================
      ELSEIF (INISPE.EQ.7) THEN
C
        COEF=0.0624D0/(0.230D0+0.0336D0*GAMMA-0.185D0/(1.9D0+GAMMA))
     &      *(DEUPI*FPIC)**4*HM0*HM0/GRAVIT**2
C
        DO 700 IP=1,NPOIN2
C
C.........CALCUL DU SPECTRE EN FREQUENCE (JONSWAP).
C         """""""""""""""""""""""""""""""""""""""""
          AL = COEF
          FP = FPIC
C
          CALL SPETMA
     &( SPEC  , FREQ  , NF    , AL    , FP     , GAMMA , SIGMAA, SIGMAB,
     &  DEUPI , GRAVIT, E2FMIN, FPMIN , DEPTH(IP) )
C
C.........CALCUL DE LA FONCTION DE REPARTITION ANGULAIRE.
C         """""""""""""""""""""""""""""""""""""""""""""""
          SPR1=SPRED1
          TET1=TETA1
          SPR2=SPRED2
          TET2=TETA2
          XLAM=XLAMDA
          IF(FRABI.EQ.2) THEN
            CALL FSPRD2
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
          ELSEIF(FRABI.EQ.3) THEN
            CALL FSPRD3
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
          ELSE
            CALL FSPRD1
     &( FRA   , TETA  , NPLAN , SPR1  , TET1  , SPR2  , TET2  , XLAM  ,
     &  DEUPI )
          ENDIF
C
C.........CALCUL DU SPECTRE DIRECTIONNEL.
C         """""""""""""""""""""""""""""""
          DO 740 JF=1,NF
            DO 750 JP=1,NPLAN
              F(IP,JP,JF)=SPEC(JF)*FRA(JP)
  750       CONTINUE
  740     CONTINUE
C
  700   CONTINUE
      ENDIF
C
      RETURN
      END
