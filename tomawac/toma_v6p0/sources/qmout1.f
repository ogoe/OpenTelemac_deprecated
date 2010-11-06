C                       *****************
                        SUBROUTINE QMOUT1
C                       *****************
C
     *( TSTOT , TSDER , F     , XK    , ENRJ  , FREQ  , FMOY  , XKMOY ,
     *  PROINF, CMOUT1, CMOUT2, GRAVIT, NF    , NPLAN , NPOIN2, TAUX1 ,
     *  BETA  )
C
C**********************************************************************
C  TOMAWAC - V1.0    P. THELLIER & M. BENOIT (EDF/DER/LNH)  -  06/04/95
C**********************************************************************
C
C  FONCTION : CALCUL DE LA CONTRIBUTION DU TERME SOURCE DE DISSIPATION
C  ********** PAR MOUTONNEMENT. UTILISATION DE LA PARAMETRISATION DE
C             KOMEN ET AL. (1984).
C
C  ARGUMENTS :
C  ***********
C +-------------+----+--------------------------------------------+
C ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C +-------------+----+--------------------------------------------+
C ! TSTOT(-,-,-)!<-->! CONTRIBUTION TERME SOURCE - PARTIE TOTALE  !
C ! TSDER(-,-,-)!<-->! CONTRIBUTION TERME SOURCE - PARTIE DERIVEE !
C ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL                       !
C ! XK(-,-)     ! -->! NOMBRE D'ONDE                              !
C ! ENRJ(-)     ! -->! ENERGIE TOTALE DU SPECTRE                  !
C ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES                     !
C ! FMOY(-)     ! -->! TABLEAU DES FREQUENCES MOYENNES            !
C ! XKMOY(-)    ! -->! TABLEAU DES NOMBRES D'ONDE MOYENS          !
C ! PROINF      ! -->! INDICATEUR DE PROFONDEUR INFINIE           !
C ! CMOUT1      ! -->! CONSTANTE DE L'EXPRESSION DE MOUTONEMENT   !
C ! CMOUT2      ! -->! CONSTANTE DE L'EXPRESSION DE MOUTONEMENT   !
C ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !
C ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !
C ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !
C ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C ! TAUX1(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C ! BETA(-)     !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C +-------------+----+--------------------------------------------+
C ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  : SEMIMP
C  ********    - PROGRAMME(S) APPELE(S) :    -
C
C  REMARQUES :
C  ***********
C  - LA CONSTANTE CMOUT1 UTILISEE DANS WAM-CYCLE 4 VAUT 4.5.
C  - LA CONSTANTE CMOUT2 UTILISEE DANS WAM-CYCLE 4 VAUT 0.5.
C
C  REFERENCE : KOMEN G.J., HASSELMANN S., HASSELMANN K. (1984) : ON
C              THE EXISTENCE OF A FULLY DEVELOPED WINDSEA SPECTRUM.
C              JPO, VOL 14, PP 1271-1285.
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NF    , NPLAN , NPOIN2
      DOUBLE PRECISION CMOUT1, CMOUT2, GRAVIT
      DOUBLE PRECISION XKMOY(NPOIN2), ENRJ(NPOIN2) ,  BETA(NPOIN2)
      DOUBLE PRECISION      FREQ(NF), FMOY(NPOIN2) , TAUX1(NPOIN2)
      DOUBLE PRECISION TSTOT(NPOIN2,NPLAN,NF), TSDER(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION     F(NPOIN2,NPLAN,NF),    XK(NPOIN2,NF)
      LOGICAL PROINF
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  JP    , JF    , IP
      DOUBLE PRECISION AUX   , DEUPI , C1    , C2
C
C
      DEUPI = 6.283185307D0
      C1 = - CMOUT1*DEUPI**9.D0/GRAVIT**4.D0
      C2 = - CMOUT1*DEUPI
C
      IF (PROINF) THEN
C     ---------------- CAS DE PROFONDEUR INFINIE (FORMULATION AVEC F).
C
C.......TABLEAU DE TRAVAIL (TERME QUI NE DEPEND QUE DU POINT D'ESPACE)
C       """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        DO IP=1,NPOIN2
          TAUX1(IP) = C1 * ENRJ(IP)**2.D0 * FMOY(IP)**9.D0
        ENDDO
C
C.......BOUCLE SUR LES FREQUENCES DE DISCRETISATION
C       """""""""""""""""""""""""""""""""""""""""""
        DO JF=1,NF
C
C.........CALCUL DU COEFFICIENT BETA : QMOUT1 = BETA * F
C         """"""""""""""""""""""""""""""""""""""""""""""
          DO IP=1,NPOIN2
            AUX = (FREQ(JF)/FMOY(IP))**2.D0
            BETA(IP)=TAUX1(IP)*((1.D0-CMOUT2)*AUX+CMOUT2*AUX**2.D0)
          ENDDO
C
C.........PRISE EN COMPTE DU TERME SOURCE.
C         """"""""""""""""""""""""""""""""
          DO JP=1,NPLAN
            DO IP=1,NPOIN2
              TSTOT(IP,JP,JF) = TSTOT(IP,JP,JF)+BETA(IP)*F(IP,JP,JF)
              TSDER(IP,JP,JF) = TSDER(IP,JP,JF)+BETA(IP)
            ENDDO
          ENDDO
        ENDDO
C
      ELSE
C     ---------------- CAS DE PROFONDEUR FINIE (FORMULATION AVEC K).
C
C.......TABLEAU DE TRAVAIL (TERME QUI NE DEPEND QUE DU POINT D'ESPACE)
C       """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        DO IP=1,NPOIN2
          TAUX1(IP) = C2 * ENRJ(IP)**2.D0 * FMOY(IP) * XKMOY(IP)**4.D0
        ENDDO 
C
C.......BOUCLE SUR LES FREQUENCES DE DISCRETISATION
C       """""""""""""""""""""""""""""""""""""""""""
        DO JF=1,NF
C
C.........CALCUL DU COEFFICIENT BETA : QMOUT1 = BETA * F
C         """"""""""""""""""""""""""""""""""""""""""""""
          DO IP=1,NPOIN2
            AUX = XK(IP,JF) / XKMOY(IP)
            BETA(IP)=TAUX1(IP)*((1.D0-CMOUT2)*AUX+CMOUT2*AUX**2.D0)
          ENDDO 
C
C.........PRISE EN COMPTE DU TERME SOURCE.
C         """"""""""""""""""""""""""""""""
          DO JP=1,NPLAN
            DO IP=1,NPOIN2
              TSTOT(IP,JP,JF) = TSTOT(IP,JP,JF)+BETA(IP)*F(IP,JP,JF)
              TSDER(IP,JP,JF) = TSDER(IP,JP,JF)+BETA(IP)
            ENDDO
          ENDDO
        ENDDO
C
      ENDIF
C
      RETURN
      END
