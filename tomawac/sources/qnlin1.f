C                       *****************
                        SUBROUTINE QNLIN1
C                       *****************
C
     *( TSTOT , TSDER , IANGNL, COEFNL, NF    , NPLAN , F1    , RAISF ,
     *  TAILF , PROINF, NPOIN2, F     , DEPTH , XKMOY , TAUX1 , TAUX2 ,
     *  TAUX3 , TAUX4 , TAUX5 , DFINI )
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  QNLIN1 :  CALCUL DE LA CONTRIBUTION DU TERME SOURCE D'INTERACTIONS C
C  ********  NON-LINEAIRES ENTRE QUADRUPLETS DE FREQUENCES A L'AIDE   C
C            DE LA METHODE "DISCRETE INTERACTION APPROXIMATION (DIA)" C
C            PROPOSEE PAR HASSELMANN ET HASSELMANN (1985).            C
C            PROCEDURE SPECIFIQUE AU CAS OU LES FREQUENCES SONT EN    C
C            PROGRESSION GEOMETRIQUE ET LES DIRECTIONS REGULIEREMENT  C
C            ESPACEES SUR [0;2.PI].                                   C
C                                                                     C
C   - CREE POUR VERSION 1.0  LE 07/06/95 PAR M. BENOIT                C
C   - MOD. POUR VERSION 1.2  LE 26/06/96 PAR M. BENOIT                C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! TSTOT(-,-,-)!<-->! CONTRIBUTION TERME SOURCE - PARTIE TOTALE  !  C
C  ! TSDER(-,-,-)!<-->! CONTRIBUTION TERME SOURCE - PARTIE DERIVEE !  C
C  ! IANGNL(-,-) ! -->! TABLEAU DES INDICES ANGULAIRES POUR DIA    !  C
C  ! COEFNL(-)   ! -->! VECTEUR DES COEFFICIENTS DE CALCUL POUR DIA!  C
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! F1          ! -->! PREMIERE FREQUENCE DE DISCRETISATION       !  C
C  ! RAISF       ! -->! RAISON FREQUENTIELLE DE DISCRETISATION     !  C
C  ! TAILF       ! -->! FACTEUR DE QUEUE                           !  C
C  ! PROINF      ! -->! INDICATEUR DE PROFONDEUR INFINIE           !  C
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !  C
C  ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL DE VARIANCE           !  C
C  ! DEPTH(-)    ! -->! TABLEAU DES PROFONDEURS (METRES)           !  C
C  ! XKMOY(-)    ! -->! TABLEAU DES NOMBRES D'ONDE MOYENS          !  C
C  ! TAUX1(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUX2(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUX3(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUX4(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUX5(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! DFINI(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  SEMIMP                     C
C  ********    - PROGRAMME(S) APPELE(S) :  CQUEUE                     C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  REMARQUES :                                                        C
C  ***********                                                        C
C                                                                     C
C  - CETTE PROCEDURE UTILISE LES RESULTATS DU SOUS-PROGRAMME PRENL1   C
C    POUR OPTIMISER LES CALCULS DE LA METHODE DIA.                    C
C                                                                     C
C  - REFERENCES PRINCIPALES SUR LA METHODE DIA :                      C
C         * HASSELMANN S., HASSELMANN K. (1985) : COMPUTATIONS        C
C               AND PARAMETERIZATIONS OF THE NONLINEAR ENERGY         C
C               TRANSFER IN GRAVITY-WAVE SPECTRUM. PART1 : A NEW      C
C               METHOD FOR EFFICIENT COMPUTATION OF THE EXACT NON-    C
C               LINEAR TRANSFER INTEGRAL. JPO, VOL 15, PP 1369-1377.  C
C         * HASSELMANN S., HASSELMANN K. ET AL.(1985) : COMPUTATIONS  C
C               AND PARAMETERIZATIONS OF THE NONLINEAR ENERGY         C
C               TRANSFER IN GRAVITY-WAVE SPECTRUM. PART2 : PARAME-    C
C               TERIZATIONS OF THE NONLINEAR ENERGY TRANSFER FOR      C
C               APPLICATION IN WAVE MODELS. JPO, VOL 15, PP 1378-1391 C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NPOIN2, NPLAN , NF
      INTEGER  IANGNL(NPLAN,8)
      DOUBLE PRECISION F1  , RAISF , TAILF
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF), COEFNL(16)
      DOUBLE PRECISION TSTOT(NPOIN2,NPLAN,NF), TSDER(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION TAUX1(NPOIN2), TAUX2(NPOIN2), TAUX3(NPOIN2)
      DOUBLE PRECISION TAUX4(NPOIN2), TAUX5(NPOIN2), XKMOY(NPOIN2)
      DOUBLE PRECISION DFINI(NPOIN2), DEPTH(NPOIN2)
      LOGICAL  PROINF
C
C.....VARIABLES LOCALES 
C     """""""""""""""""
      INTEGER  JBP0  , JFP0  , JFP1  , JFM0  , JFM1  , JFP   , JFM
      INTEGER  JBP1  , JB    , JBM0  , JBM1  , IMAGE , JP 
      INTEGER  JPP0  , JPP1  , JPM0  , JPM1  , IP    , KAUX  , JF    
      INTEGER  JFMIN , JFMAX
      DOUBLE PRECISION COEFP0, COEFP1, COEFM0, COEFM1, COEFJF, XXFAC
      DOUBLE PRECISION FMOIN , FPLUS , TERM1 , TERM2 , US1PL4, US1ML4
      DOUBLE PRECISION C1    , C2    , C3    , C4    , C5    , C6
      DOUBLE PRECISION D1    , D2    , D3    , D4    , D5    , D6
      DOUBLE PRECISION C1SQ  , C2SQ  , C3SQ  , C4SQ  , C5SQ  , C6SQ
      DOUBLE PRECISION C7    , C8    , D7    , D8    , C7SQ  , C8SQ
      DOUBLE PRECISION TERM3 , FDEJF , FREQ
C
C
C.....RECUPERATION DES COEFFICIENTS CALCULES PAR PRENL1.
C     """"""""""""""""""""""""""""""""""""""""""""""""""
      C1    = COEFNL( 1)
      C2    = COEFNL( 2)
      C3    = COEFNL( 3)
      C4    = COEFNL( 4)
      C5    = COEFNL( 5)
      C6    = COEFNL( 6)
      C7    = COEFNL( 7)
      C8    = COEFNL( 8)
      JFP   = IDINT(COEFNL( 9)+1.D-7)
      JFM   = IDINT(COEFNL(10)-1.D-7)
      US1PL4= COEFNL(11)
      US1ML4= COEFNL(12)
      JFMIN = NINT(COEFNL(13))
      JFMAX = NINT(COEFNL(14))
      C1SQ  = C1*C1
      C2SQ  = C2*C2
      C3SQ  = C3*C3
      C4SQ  = C4*C4
      C5SQ  = C5*C5
      C6SQ  = C6*C6
      C7SQ  = C7*C7
      C8SQ  = C8*C8
C
C.....FACTEUR CORRECTIF POUR UNE PROFONDEUR FINIE.
C     """"""""""""""""""""""""""""""""""""""""""""
      IF (.NOT.PROINF) THEN
C
        DO 100 IP=1,NPOIN2
          TERM1 = MAX(0.75D0*DEPTH(IP)*XKMOY(IP),0.5D0)
          DFINI(IP) = 1.D0+(5.5D0/TERM1)*(1.D0-0.833D0*TERM1)
     *               *DEXP(-1.25D0*TERM1)
  100   CONTINUE
      ENDIF
C
C.....PREMIERE BOUCLE SUR LES FREQUENCES.
C     """""""""""""""""""""""""""""""""""
      DO 200 JF=JFMIN,JFMAX
C
C.......CALCUL DE LA FREQUENCE CONSIDEREE.
C       """"""""""""""""""""""""""""""""""
        FREQ = F1*RAISF**(JF-1)
C
C.......INDICES FREQUENTIELS ENTROURANT LA FREQUENCE 'SUPERIEURE'
C       FREQ(JFP0) < (1+LAMBDA).FREQ(JF) < FREQ(JFP1)
C       """""""""""""""""""""""""""""""""""""""""""""""""""""""""
        JFP0=JF+JFP
        JFP1=JFP0+1
C
C.......INDICES FREQUENTIELS ENTROURANT LA FREQUENCE 'INFERIEURE'
C       FREQ(JFM0) < (1-LAMBDA).FREQ(JF) < FREQ(JFM1)
C       """""""""""""""""""""""""""""""""""""""""""""""""""""""""
        JFM0=JF+JFM-1
        JFM1=JFM0+1
C
C.......LIMITATION DES INDICES A NF ET PRISE EN COMPTE DE LA QUEUE
C       DU SPECTRE DE FACON ANALYTIQUE (DECROISSANCE EN -TAILF).
C       """"""""""""""""""""""""""""""""""""""""""""""""""""""""""
        CALL CQUEUE( NF , RAISF , TAILF , JFP1 , JBP1 , COEFP1 )
        CALL CQUEUE( NF , RAISF , TAILF , JFP0 , JBP0 , COEFP0 )
        CALL CQUEUE( NF , RAISF , TAILF , JF   , JB   , COEFJF )
        CALL CQUEUE( NF , RAISF , TAILF , JFM1 , JBM1 , COEFM1 )
        CALL CQUEUE( NF , RAISF , TAILF , JFM0 , JBM0 , COEFM0 )
C
C.......COEFFICIENTS D'INTERPOLATION DU SPECTRE MODIFIES.
C       """""""""""""""""""""""""""""""""""""""""""""""""
        D1=C1*COEFP0*US1PL4
        D2=C2*COEFP0*US1PL4
        D3=C3*COEFP1*US1PL4
        D4=C4*COEFP1*US1PL4
        D5=C5*COEFM0*US1ML4
        D6=C6*COEFM0*US1ML4
        D7=C7*COEFM1*US1ML4
        D8=C8*COEFM1*US1ML4
C
C.......CALCUL DU COEFFICIENT MULTIPLICATIF (EN F**11) ET PRISE
C       EN COMPTE DU TERME CORRECTIF EN PROFONDEUR FINIE.
C       """""""""""""""""""""""""""""""""""""""""""""""""""""""
        XXFAC= 3000.D0*FREQ**11
        IF (PROINF) THEN
          DO 210 IP=1,NPOIN2
            TAUX1(IP) = XXFAC
  210     CONTINUE
        ELSE
          DO 220 IP=1,NPOIN2
            TAUX1(IP) = DFINI(IP)*XXFAC
  220     CONTINUE
        ENDIF
C
C.......DEUXIEME BOUCLE SUR SYMETRIE ANGULAIRE.
C       """""""""""""""""""""""""""""""""""""""
        DO 300 IMAGE=1,2
       KAUX=(IMAGE-1)*4
C
C........TROISIEME BOUCLE SUR LES DIRECTIONS.
C        """"""""""""""""""""""""""""""""""""
         DO 400 JP=1,NPLAN
          JPP0 = IANGNL(JP,KAUX+1)
          JPP1 = IANGNL(JP,KAUX+2)
          JPM0 = IANGNL(JP,KAUX+3)
          JPM1 = IANGNL(JP,KAUX+4)
C
C
          IF (JFM0.LT.1) THEN
C........./-------------------------------------------------------/
C........./ AU MOINS UNE DES FREQUENCES EST INFERIEURE A FREQ(1)  /
C........./ (LE SPECTRE F- A LA FREQUENCE (1-XLAMD).FREQ EST NUL) /
C........./-------------------------------------------------------/
C
          DO 500 IP=1,NPOIN2
           FDEJF = F(IP,JP,JB )*COEFJF
           FPLUS = F(IP,JPP0,JBP0)*D1 + F(IP,JPP1,JBP0)*D2
     *           + F(IP,JPP0,JBP1)*D3 + F(IP,JPP1,JBP1)*D4
C
           TERM1 = FDEJF*FPLUS
           TERM3 = TAUX1(IP)*FDEJF
C
           TAUX2(IP) = TERM1*TERM3
           TAUX3(IP) = 2.D0*TERM1*TAUX1(IP)
           TAUX5(IP) = FDEJF*US1PL4*TERM3
  500     CONTINUE
C
          IF (JB.EQ.JF) THEN
C
          DO 510 IP=1,NPOIN2
           TSTOT(IP,JP  ,JF  )=TSTOT(IP,JP  ,JF  )-TAUX2(IP)*2.D0
           TSDER(IP,JP  ,JF  )=TSDER(IP,JP  ,JF  )-TAUX3(IP)*2.D0
  510     CONTINUE
C
          IF (JBP0.EQ.JFP0) THEN
C
          DO 520 IP=1,NPOIN2
           TSTOT(IP,JPP0,JFP0)=TSTOT(IP,JPP0,JFP0)+TAUX2(IP)*C1
           TSTOT(IP,JPP1,JFP0)=TSTOT(IP,JPP1,JFP0)+TAUX2(IP)*C2
           TSDER(IP,JPP0,JFP0)=TSDER(IP,JPP0,JFP0)+TAUX5(IP)*C1SQ
           TSDER(IP,JPP1,JFP0)=TSDER(IP,JPP1,JFP0)+TAUX5(IP)*C2SQ
  520     CONTINUE
C
          IF (JBP1.EQ.JFP1) THEN
C
          DO 530 IP=1,NPOIN2
           TSTOT(IP,JPP0,JFP1)=TSTOT(IP,JPP0,JFP1)+TAUX2(IP)*C3
           TSTOT(IP,JPP1,JFP1)=TSTOT(IP,JPP1,JFP1)+TAUX2(IP)*C4
           TSDER(IP,JPP0,JFP1)=TSDER(IP,JPP0,JFP1)+TAUX5(IP)*C3SQ
           TSDER(IP,JPP1,JFP1)=TSDER(IP,JPP1,JFP1)+TAUX5(IP)*C4SQ
  530     CONTINUE
          ENDIF
          ENDIF
          ENDIF
C
          ELSE
C........./--------------------------------------------------------/
C........./ LES FREQUENCES F-, F, F+ PEUVENT CONTENIR DE L'ENERGIE /
C........./--------------------------------------------------------/
C
          DO 600 IP=1,NPOIN2
           FDEJF = F(IP,JP,JB )*COEFJF
           FPLUS = F(IP,JPP0,JBP0)*D1 + F(IP,JPP1,JBP0)*D2
     *           + F(IP,JPP0,JBP1)*D3 + F(IP,JPP1,JBP1)*D4
           FMOIN = F(IP,JPM0,JBM0)*D5 + F(IP,JPM1,JBM0)*D6
     *           + F(IP,JPM0,JBM1)*D7 + F(IP,JPM1,JBM1)*D8
C
           TERM1 = FDEJF*(FPLUS+FMOIN)
           TERM2 = 2.D0*FPLUS*FMOIN
           TERM3 = TAUX1(IP)*FDEJF
C
           TAUX2(IP) = (TERM1-TERM2)*TERM3
           TAUX3(IP) = (2.D0*TERM1-TERM2)*TAUX1(IP)
           TAUX5(IP) = (FDEJF-2.D0*FMOIN)*US1PL4*TERM3
           TAUX4(IP) = (FDEJF-2.D0*FPLUS)*US1ML4*TERM3
  600     CONTINUE
C
          IF (JBM0.EQ.JFM0) THEN
C
          DO 710 IP=1,NPOIN2
           TSTOT(IP,JPM0,JFM0)=TSTOT(IP,JPM0,JFM0)+TAUX2(IP)*C5
           TSTOT(IP,JPM1,JFM0)=TSTOT(IP,JPM1,JFM0)+TAUX2(IP)*C6
           TSDER(IP,JPM0,JFM0)=TSDER(IP,JPM0,JFM0)+TAUX4(IP)*C5SQ
           TSDER(IP,JPM1,JFM0)=TSDER(IP,JPM1,JFM0)+TAUX4(IP)*C6SQ
  710     CONTINUE
C
          IF (JBM1.EQ.JFM1) THEN
C
          DO 720 IP=1,NPOIN2
           TSTOT(IP,JPM0,JFM1)=TSTOT(IP,JPM0,JFM1)+TAUX2(IP)*C7
           TSTOT(IP,JPM1,JFM1)=TSTOT(IP,JPM1,JFM1)+TAUX2(IP)*C8
           TSDER(IP,JPM0,JFM1)=TSDER(IP,JPM0,JFM1)+TAUX4(IP)*C7SQ
           TSDER(IP,JPM1,JFM1)=TSDER(IP,JPM1,JFM1)+TAUX4(IP)*C8SQ
  720     CONTINUE
C
          IF (JB.EQ.JF) THEN
C
          DO 730 IP=1,NPOIN2
           TSTOT(IP,JP  ,JF  )=TSTOT(IP,JP  ,JF  )-TAUX2(IP)*2.D0
           TSDER(IP,JP  ,JF  )=TSDER(IP,JP  ,JF  )-TAUX3(IP)*2.D0
  730     CONTINUE
C
          IF (JBP0.EQ.JFP0) THEN
C
          DO 740 IP=1,NPOIN2
           TSTOT(IP,JPP0,JFP0)=TSTOT(IP,JPP0,JFP0)+TAUX2(IP)*C1
           TSTOT(IP,JPP1,JFP0)=TSTOT(IP,JPP1,JFP0)+TAUX2(IP)*C2
           TSDER(IP,JPP0,JFP0)=TSDER(IP,JPP0,JFP0)+TAUX5(IP)*C1SQ
           TSDER(IP,JPP1,JFP0)=TSDER(IP,JPP1,JFP0)+TAUX5(IP)*C2SQ
  740     CONTINUE
C
          IF (JBP1.EQ.JFP1) THEN
C
          DO 750 IP=1,NPOIN2
           TSTOT(IP,JPP0,JFP1)=TSTOT(IP,JPP0,JFP1)+TAUX2(IP)*C3
           TSTOT(IP,JPP1,JFP1)=TSTOT(IP,JPP1,JFP1)+TAUX2(IP)*C4
           TSDER(IP,JPP0,JFP1)=TSDER(IP,JPP0,JFP1)+TAUX5(IP)*C3SQ
           TSDER(IP,JPP1,JFP1)=TSDER(IP,JPP1,JFP1)+TAUX5(IP)*C4SQ
  750     CONTINUE
          ENDIF
          ENDIF
          ENDIF
          ENDIF
          ENDIF
C
          ENDIF
C
  400    CONTINUE
C
  300   CONTINUE
C
  200 CONTINUE
C
      RETURN
      END
