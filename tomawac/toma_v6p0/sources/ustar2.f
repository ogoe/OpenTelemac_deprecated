C                       *****************
                        SUBROUTINE USTAR2
C                       *****************
C
     *( USTAR , UV    , VV    , NPOIN2)
C
C**********************************************************************
C  TOMAWAC - V1.1    M. BENOIT               (EDF/DER/LNH)  -  26/03/96
C**********************************************************************
C
C  FONCTION : CALCUL DE LA VITESSE DE FROTTEMENT U* EN TOUS LES POINTS
C  ********** DU MAILLAGE SPATIAL.
C             UTILISATION D'UN COEFFICIENT DE TRAINEE VARIANT DE FACON
C             LINEAIRE AVEC LA VITESSE DU VENT. EXPRESSION IDENTIQUE A
C             CELLE UTILISEE DANS WAM-CYCLE 3 (WAMDI GROUP, 1988).
C
C  ARGUMENTS :
C  ***********
C +-------------+----+--------------------------------------------+
C ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C +-------------+----+--------------------------------------------+
C ! USTAR(-)    !<-- ! TABLEAU DES VITESSES DE FROTTEMENT         !
C ! UV(-)       ! -->! TABLEAU DES COMPOSANTES OUEST-EST DU VENT  !
C ! VV(-)       ! -->! TABLEAU DES COMPOSANTES SUD-NORD  DU VENT  !
C ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C +-------------+----+--------------------------------------------+
C ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  SEMIMP
C  ********    - PROGRAMME(S) APPELE(S) :  -
C
C  REMARQUES :
C  ***********
C
C  REFRENCES : - WAMDI GROUP (1988) : THE WAM MODEL - A THIRD GENERA-
C              TION OCEAN WAVE PREDICTION MODEL. JPO, VOL 18, PP 1775-
C              1810.
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NPOIN2
      DOUBLE PRECISION USTAR(NPOIN2) , UV(NPOIN2) , VV(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  IP
      DOUBLE PRECISION UVENT , CDRAG
C
C
C.....BOUCLE PRINCIPALE SUR LES POINTS DU MAILLAGE SPATIAL.
C     """""""""""""""""""""""""""""""""""""""""""""""""""""
      DO 100 IP=1,NPOIN2
C
C.......CALCUL DE LA VITESSE DU VENT A 10 METRES.
C       """""""""""""""""""""""""""""""""""""""""
        UVENT=DSQRT(UV(IP)**2+VV(IP)**2)
C
C.......CALCUL DU COEFFICIENT DE TRAINEE.
C       """""""""""""""""""""""""""""""""
        CDRAG = 6.5D-5*UVENT + 8.D-4
        IF (UVENT.LT.7.5D0) CDRAG = 1.2875D-3
C
C.......CALCUL DE LA VITESSE DE FROTTEMENT.
C       """""""""""""""""""""""""""""""""""
        USTAR(IP)=DSQRT(CDRAG)*UVENT
C
  100 CONTINUE
C
      RETURN
      END
