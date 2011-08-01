C                       *****************
                        SUBROUTINE STRESS
C                       *****************
C
     *( TAUWAV, TSTOT , F     , USNEW , TWNEW , Z0NEW , FREQ  , DFREQ ,
     *  TETA  , SINTET, COSTET, ROAIR , ROEAU , XKAPPA, BETAM , DECAL ,
     *  GRAVIT, NPOIN2, NPLAN , NF    , XTAUW , YTAUW , TAUHF )
C
C**********************************************************************
C  TOMAWAC - V1.0    M. BENOIT               (EDF/DER/LNH)  -  03/05/95
C**********************************************************************
C
C  FONCTION : CALCUL DE LA CONTRAINTE DE HOULE EN TOUS LES POINTS DU
C  ********** MAILLAGE SPATIAL.
C
C  ARGUMENTS :
C  ***********
C  +-------------+----+--------------------------------------------+
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C  +-------------+----+--------------------------------------------+
C  ! TAUWAV(-)   !<-- ! TABLEAU DES CONTRAINTES DUES A LA HOULE    !
C  ! TSTOT(-,-,-)! -->! CONTRIBUTION TERME SOURCE - PARTIE TOTALE  !
C  ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL                       !
C  ! USNEW(-)    ! -->! TABLEAU DES VITESSES DE FROTTEMENT         !
C  ! TWNEW(-)    ! -->! TABLEAU DES DIRECTIONS DU VENT             !
C  ! Z0NEW(-)    ! -->! TABLEAU DES LONGUEURS DE RUGOSITE          !
C  ! FREQ(-)     ! -->! VECTEUR DES FREQUENCES DE DISCRETISATION   !
C  ! DFREQ(-)    ! -->! VECTEUR DES INTERVALLES FREQUENTIELS       !
C  ! TETA(-)     ! -->! VECTEUR DES DIRECTIONS DE DISCRETISATION   !
C  ! SINTET(-)   ! -->! VECTEUR DES   SINUS DES DIRECTIONS         !
C  ! COSTET(-)   ! -->! VECTEUR DES COSINUS DES DIRECTIONS         !
C  ! ROAIR       ! -->! MASSE VOLUMIQUE DE L AIR                   !
C  ! ROEAU       ! -->! MASSE VOLUMIQUE DE L EAU                   !
C  ! XKAPPA      ! -->! CONSTANTE DE VON KARMAN                    !
C  ! BETAM       ! -->! CONSTANTE BETAMAX DE LA FORMULE DU VENT    !
C  ! DECAL       ! -->! CONSTANTE DE DECALAGE DE CROISSANCE VENT   !
C  ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !
C  ! XTAUW(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C  ! YTAUW(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C  ! TAUHF(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C  +-------------+----+--------------------------------------------+
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C  +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  SEMIMP
C  ********    - PROGRAMME(S) APPELE(S) :    -   
C
C  REMARQUES :
C  ***********
C  - LA DETERMINATION DE TAUHF A PARTIR DE USTAR ET ALFA SE FAIT
C    DIRECTEMENT DANS CETTE PROCEDURE (AUPARAVANT, ELLE SE FAISAIT PAR
C    APPEL A LA PROCEDURE TAUWHF).
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NPOIN2, NPLAN , NF
      DOUBLE PRECISION ROAIR , ROEAU , XKAPPA , BETAM , DECAL , GRAVIT
      DOUBLE PRECISION TAUWAV(NPOIN2), USNEW(NPOIN2), TWNEW(NPOIN2)
      DOUBLE PRECISION  Z0NEW(NPOIN2), FREQ(NF) , DFREQ(NF)
      DOUBLE PRECISION TSTOT(NPOIN2,NPLAN,NF), F(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION TETA(NPLAN), SINTET(NPLAN), COSTET(NPLAN)
      DOUBLE PRECISION XTAUW(NPOIN2) , YTAUW(NPOIN2), TAUHF(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  IP    , JP    , JF    , JTOT  , J
      DOUBLE PRECISION DEUPI , DTETAR, FRMAX , COEF1 , COEF2 , TTAUHF
      DOUBLE PRECISION USTAR , ALFA  , COSTMP, C1    , C2    , DIREC
      DOUBLE PRECISION Y     , OMEGA , CM    , ZX    , ZARG  , ZMU
      DOUBLE PRECISION ZLOG  , ZBETA , UST   , Z0    , UST2  , ALF2
      DOUBLE PRECISION CONST1, OMEGAM, X0    , YC    , DELY  , AUX
C
C
      DEUPI = 6.283185307D0
      DTETAR= DEUPI/DBLE(NPLAN)
      FRMAX = FREQ(NF)
      COEF1 = DTETAR*DEUPI**4*FRMAX**5/GRAVIT**2
      COEF2 = DEUPI*ROEAU/ROAIR*DTETAR
C
      DO 12 IP=1,NPOIN2
        XTAUW(IP)=0.D0
        YTAUW(IP)=0.D0
   12 CONTINUE
C
C.....INTEGRATION DU TERME SOURCE SUR LES FREQUENCES ET DIRECTIONS.
C     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      DO 110 JF=1,NF
        AUX=COEF2*FREQ(JF)*DFREQ(JF)
        DO 120 JP=1,NPLAN
          C1=AUX*SINTET(JP)
          C2=AUX*COSTET(JP)
C
          DO 100 IP=1,NPOIN2
            XTAUW(IP)=XTAUW(IP)+TSTOT(IP,JP,JF)*C1
            YTAUW(IP)=YTAUW(IP)+TSTOT(IP,JP,JF)*C2
  100     CONTINUE
  120   CONTINUE
  110 CONTINUE
C
C.....CALCUL DE LA PARTIE HAUTES-FREQUENCES PARAMETRISEE.
C     """""""""""""""""""""""""""""""""""""""""""""""""""
      DO 170 IP=1,NPOIN2
        TAUHF(IP)=0.D0
  170 CONTINUE
C
      DO 200 JP=1,NPLAN
        DIREC=TETA(JP)
        DO 171 IP=1,NPOIN2
          COSTMP=MAX(COS(DIREC-TWNEW(IP)),0.D0)
          TAUHF(IP)=TAUHF(IP)+F(IP,JP,NF)*COSTMP**3
  171   CONTINUE
  200 CONTINUE
C
      JTOT  = 50
      CONST1= BETAM/XKAPPA**2
      OMEGAM= 6.283185307D0*FRMAX
      X0    = 0.05D0
C
      DO 173 IP=1,NPOIN2
        USTAR=USNEW(IP)
        ALFA =Z0NEW(IP)*GRAVIT/USTAR**2
C
C----------------------------------------------ANCIENNE SUB TAUWHF
CC      CALL TAUWHF
CC   *( TAUHF , USTAR , ALFA  , BETAM , XKAPPA , DECAL , FRMAX , GRAVIT)
C----------------------------------------------
C
CMB.....LIMITATIONS POUR REPRODUIRE WAM4 (PROCEDURE TAUHF DE PREPROC)
C.......(CE "BRIDAGE" NE PARAIT PAS JUSTIFIE A PRIORI. IL EST PRODUIT
C.......PAR LE FAIT QUE TAUHF EST DISCRETISE SUR UNE GRILLE.)
        UST2  = MIN(USTAR,5.D0)
        ALF2  = MIN(ALFA,0.11D0)
CMB.....FIN DU BRIDAGE.
C
        UST   = MAX(UST2,0.000001D0)
        Z0    = ALF2*UST**2/GRAVIT
C
        YC    = MAX(OMEGAM,X0*GRAVIT/UST)*SQRT(Z0/GRAVIT)
        DELY  = MAX((1.D0-YC)/FLOAT(JTOT),0.D0)
        TTAUHF = 0.D0
        DO 102 J=1,JTOT
          Y     = YC+DBLE(J-1)*DELY
          OMEGA = Y*SQRT(GRAVIT/Z0)
          CM    = GRAVIT/OMEGA
          ZX    = UST/CM +DECAL
          ZARG  = MIN(XKAPPA/ZX,20.D0)
          ZMU   = MIN(GRAVIT*Z0/CM**2*DEXP(ZARG),1.D0)
          ZLOG  = MIN(DLOG(ZMU),0.D0)
          ZBETA = CONST1*ZMU*ZLOG**4
          TTAUHF= TTAUHF+ZBETA/Y*DELY
  102   CONTINUE
C----------------------------------------------ANCIENNE SUB TAUWHF
C
        TAUHF(IP) = TTAUHF*COEF1*USTAR**2*TAUHF(IP)
  173 CONTINUE
C
C.....PRISE EN COMPTE DE LA PARTIE HAUTES-FREQUENCES PARAMETRISEE.
C     """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
C
      DO 190 IP=1,NPOIN2
        XTAUW(IP) = XTAUW(IP) + TAUHF(IP)*SIN(TWNEW(IP))
        YTAUW(IP) = YTAUW(IP) + TAUHF(IP)*COS(TWNEW(IP))
        TAUWAV(IP)= SQRT(XTAUW(IP)**2+YTAUW(IP)**2)
  190 CONTINUE
C
      RETURN
      END
