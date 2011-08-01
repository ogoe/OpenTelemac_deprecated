C                       *****************
                        SUBROUTINE TAUTOT
C                       *****************
C
     *( TAUT  , UVENT , TAUW  , CDRAG , ALPHA , XKAPPA, ZVENT , SEUIL ,
     *  GRAVIT, ITR   , ITRMIN, ITRMAX)
C
C**********************************************************************
C  TOMAWAC - V1.0    M. BENOIT               (EDF/DER/LNH)  -  25/04/95
C**********************************************************************
C
C  FONCTION : CALCUL DE LA CONTRAINTE TOTALE A PARTIR DE LA VITESSE DU
C  ********** VENT UVENT A LA COTE ZVENT (EN PRINCIPE ZVENT=10 M) ET
C             DE LA CONTRAINTE DUE A LA HOULE TAUW.
C             THEORIE DEVELOPPEE PAR JANSSEN (1989 ET 1991) ET UTILISEE
C             DANS WAM-CYCLE 4 (PROCEDURE STRESS DE PREPROC).
C
C  ARGUMENTS :
C  ***********
C  +-------------+----+--------------------------------------------+
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C  +-------------+----+--------------------------------------------+
C  ! TAUT        !<-- ! CONTRAINTE TOTALE                          !
C  ! UVENT       ! -->! VITESSE DU VENT A LA COTE ZVENT (M/S)      !
C  ! TAUW        ! -->! CONTRAINTE DUE A LA HOULE                  !
C  ! CDRAG       ! -->! COEFFICIENT DE TRAINEE                     !
C  ! ALPHA       ! -->! CONSTANTE DE LA LOI DE CHARNOCK            !
C  ! XKAPPA      ! -->! CONSTANTE DE VON KARMAN                    !
C  ! ZVENT       ! -->! COTE A LAQUELLE EST MESURE LE VENT (M)     !
C  ! SEUIL       ! -->! SEUIL DE CONVERGENCE METHODE DE NEWTON     !
C  ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !
C  ! ITR         !<-- ! NOMBRE D'ITERATIONS EFFECTUES              !
C  ! ITRMIN      ! -->! NOMBRE MINIMAL D'ITERATIONS SOUHAITE       !
C  ! ITRMAX      ! -->! NOMBRE MAXIMAL D'ITERATIONS SOUHAITE       !
C  +-------------+----+--------------------------------------------+
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C  +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  AIRSEA
C  ********    - PROGRAMME(S) APPELE(S) :    -
C
C  REMARQUES :
C  ***********
C  - LA DETERMINATION DE TAUT A PARTIR DE UVENT ET TAUW SE FAIT EN
C    RESOLVANT UNE EQUATION IMPLICITE PAR UNE METHODE ITERATIVE DE
C    NEWTON.
C
C  REFRENCES : - JANSSEN P.A.E.M (1989) : WIND-INDUCED STRESS AND THE
C              DRAG OF AIR FLOW OVER SEA WAVES. JPO, VOl 19, PP 745-754
C              - JANSSEN P.A.E.M (1991) : QUASI-LINEAR THEORY OF WIND-
C              WAVE GENERATION APPLIED TO WAVE FORECASTING. JPO, VOL 21
C              PP 1631-1642.
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  ITRMIN, ITRMAX, ITR
      DOUBLE PRECISION TAUT  , UVENT , TAUW  , ALPHA , XKAPPA , ZVENT
      DOUBLE PRECISION SEUIL , CDRAG , GRAVIT
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      DOUBLE PRECISION TAUMIN, XNUAIR, AUX   , USTO  , TAUO   , TAUN
      DOUBLE PRECISION USTN  , X     , ZE    , DIFF  , FOLD   , DFOLD
C
C
      ITR   =0
      TAUMIN=1.D-5
      XNUAIR=1.D-5
C
C.....VALEURS INITIALES.
C     """"""""""""""""""
      USTO  =UVENT*SQRT(CDRAG)
      TAUO  =MAX(USTO**2,TAUW+TAUMIN)
C
  190 CONTINUE
      ITR   = ITR+1
C
C.....ITERATION PAR LA METHODE DE NEWTON.
C     """""""""""""""""""""""""""""""""""
      USTO  = SQRT(TAUO)
      X     = TAUW/TAUO
      ZE    = MAX(0.1D0*XNUAIR/USTO,ALPHA*TAUO/(GRAVIT*SQRT(1.D0-X)))
      AUX   = DLOG(ZVENT/ZE)
      FOLD  = USTO-XKAPPA*UVENT/AUX
      DFOLD = 1.D0-2.D0*XKAPPA*UVENT*(1.D0-1.5D0*X)/AUX**2/USTO/(1.D0-X)
      USTN  = USTO-FOLD/DFOLD
      TAUN  = MAX(USTN**2,TAUW+TAUMIN)
C
C.....CRITERES DE CONVERGENCE.
C     """"""""""""""""""""""""
      DIFF=ABS(TAUN-TAUO)/TAUO
      TAUO=TAUN
      IF (ITR.LT.ITRMIN) GOTO 190
      IF ((DIFF.GT.SEUIL).AND.(ITR.LT.ITRMAX)) GOTO 190
C
C.....AFFECTATION DE LA SOLUTION TROUVEE.
C     """""""""""""""""""""""""""""""""""
      TAUT=TAUN
C
      RETURN
      END
