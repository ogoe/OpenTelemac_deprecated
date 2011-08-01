C                       *****************
                        SUBROUTINE WNSCOU
C                       *****************
C
     *( CK2   , FREQ  , DEPTH )
C
C**********************************************************************
C  TOMAWAC - V1.0    M. BENOIT               (EDF/DER/LNH)  -  07/02/95
C**********************************************************************
C
C  FONCTION : CALCUL DU NOMBRE D'ONDE DE LA HOULE PAR RESOLUTION DE
C  ********** L'EQUATION DE DISPERSION DANS LE CAS SANS COURANT.
C
C  ARGUMENTS :
C  ***********
C  +-------------+----+--------------------------------------------+
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C  +-------------+----+--------------------------------------------+
C  ! CK2         !<-- ! NOMBRE D'ONDE                              !
C  ! FREQ        ! -->! FREQUENCE                                  !
C  ! DEPTH       ! -->! PROFONDEUR                                 !
C  +-------------+----+--------------------------------------------+
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C  +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :
C  ********    - PROGRAMME(S) APPELE(S) : NEANT
C
C  REMARQUES :
C  ***********
C  -> 3 METHODES SONT UTILISEES SUIVANT LES VALEURS DE K0*D
C     (K0 ETANT LE NOMBRE D'ONDE EN PROFONDEUR INFINIE)
C
C                      3.2                  5.6
C      -----------------!--------------------!-----------------> K0*D
C                       !                    !
C   METHODE EXPLICITE   ! METHODE ITERATIVE  !     K=K0
C   (HUNT A L'ORDRE 9)  !                    ! (PROF. INFINIE)
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      DOUBLE PRECISION CK2   , FREQ  , DEPTH
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  I
      DOUBLE PRECISION P(9)  , EPS   , XG    , DEUPI , XK0   , XK0D
      DOUBLE PRECISION AUX   , A     , Y     , YI    , OM
C
      DATA DEUPI/6.28318531D0/
      DATA EPS/0.0001D0/
      DATA XG/9.806D0/
      DATA P/0.66667D0,0.35550D0,0.16084D0,0.06320D0,0.02174D0,
     &       0.00654D0,0.00171D0,0.00039D0,0.00011D0/
C
C
C
C.....CALCUL DE LA PULSATION (OM), K0 ET K0D
      OM=FREQ*DEUPI
      XK0=OM*OM/XG
      XK0D=XK0*DEPTH
C
C.....TEST SUR LA VALEUR DE XK0D POUR DEFINIR LA METHODE DE RESOLUTION
C     ================================================================
C
      IF (XK0D.LE.3.2) THEN
C.......METHODE EXPLICITE DE HUNT A L'ORDRE 9.
        Y=XK0*DEPTH
        AUX=1.D0
        YI=1.D0
        DO 12 I=1,9
          YI=YI*Y
          AUX=AUX+P(I)*YI
   12   CONTINUE
        AUX=Y+1.D0/AUX
        CK2=OM/SQRT(XG*DEPTH/AUX)
C
      ELSEIF (XK0D.LE.5.6) THEN
C.......METHODE ITERATIVE A PARTIR DE LA METHODE DE HUNT A L'ORDRE 9.
        Y=XK0*DEPTH
        AUX=1.D0
        YI=1.D0
        DO 11 I=1,9
          YI=YI*Y
          AUX=AUX+P(I)*YI
   11   CONTINUE
        AUX=Y+1.D0/AUX
        CK2=OM/SQRT(XG*DEPTH/AUX)
    2   CONTINUE
        A=CK2
        CK2=XK0/TANH(A*DEPTH)
        IF (ABS(CK2-A)/CK2.LT.EPS) GOTO 3
        GOTO 2
    3   CONTINUE
C
      ELSE
C.......APPROXIMATION DE PROFONDEUR INFINIE
        CK2=XK0
      ENDIF
C
      RETURN
      END
