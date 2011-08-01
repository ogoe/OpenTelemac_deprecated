C                       *****************
                        SUBROUTINE SPETMA
C                       *****************
C
     &( SPEC  , FREQ  , NF    , AL    , FP     , GAMMA , SIGMAA, SIGMAB,
     &  DEUPI , GRAVIT, E2FMIN, FPMIN , DEPTH  )
C
C**********************************************************************
C  TOMAWAC - V1.1    M. BENOIT               (EDF/DER/LNH)  -  10/01/96
C**********************************************************************
C
C  FONCTION : CALCUL DU SPECTRE EN FREQUENCE DE TYPE TMA POUR UNE
C  ********** SERIE DE FREQUENCES DONNEES.
C
C  ARGUMENTS :
C  ***********
C  +-------------+----+--------------------------------------------+
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C  +-------------+----+--------------------------------------------+
C  ! SPEC(-)     !<-- ! VALEURS DU SPECTRE TMA                     !
C  ! FREQ(-)     ! -->! FREQUENCES DE DISCRETISATION               !
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !
C  ! AL          ! -->! CONSTANTE DE PHILLIPS (ALPHA)              !
C  ! FP          ! -->! FREQUENCE DE PIC DU SPECTRE JONSWAP        !
C  ! GAMMA       ! -->! FACTEUR DE FORME DE PIC JONSWAP            !
C  ! SIGMAA      ! -->! VALEUR DE SIGMA JONSWAP POUR F < FP        !
C  ! SIGMAB      ! -->! VALEUR DE SIGMA JONSWAP POUR F > FP        !
C  ! DEUPI       ! -->! 2.PI                                       !
C  ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !
C  ! E2FMIN      ! -->! SEUIL MINIMUM DE SPECTRE CONSIDERE         !
C  ! FPMIN       ! -->! VALEUR MINIMUM DE LA FREQUENCE DE PIC      !
C  ! DEPTH       ! -->! PROFONDEUR D'EAU AU POINT CONSIDERE (M)    !
C  +-------------+----+--------------------------------------------+
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C  +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  SPEINI
C  ********    - PROGRAMME(S) APPELE(S) :  -      
C
C  REMARQUES :
C  ***********
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NF
      DOUBLE PRECISION GRAVIT, SIGMAA, SIGMAB, GAMMA , DEUPI , FPMIN
      DOUBLE PRECISION FP    , E2FMIN, AL    , DEPTH
      DOUBLE PRECISION SPEC(NF)      , FREQ(NF)
C
C.....VARIABLES LOCALES 
C     """""""""""""""""
      INTEGER  JF
      DOUBLE PRECISION COEF  , ARG1   , ARG2  , ARG3  , SIG   , FF
      DOUBLE PRECISION ARG4  , OMEGH
C
C
      IF (FP.GT.FPMIN) THEN
        COEF=AL*GRAVIT**2/DEUPI**4
        DO 100 JF=1,NF
          FF=FREQ(JF)
          IF (FF.LT.FP) THEN 
            SIG=SIGMAA
          ELSE
            SIG=SIGMAB
          ENDIF
          ARG1=0.5D0*((FF-FP)/(SIG*FP))**2
          IF (ARG1.LT.99.D0) THEN
            ARG1=GAMMA**EXP(-ARG1)
          ELSE
            ARG1=1.D0
          ENDIF
          ARG2=1.25D0*(FP/FF)**4
          IF (ARG2.LT.99.D0) THEN
            ARG2=EXP(-ARG2)
          ELSE
            ARG2=0.D0
          ENDIF
          ARG3=COEF/FF**5
          OMEGH=DEUPI*FF*SQRT(DEPTH/GRAVIT)
          IF (OMEGH.LT.1.D0) THEN
            ARG4=0.5D0*OMEGH*OMEGH
          ELSEIF (OMEGH.LT.2.D0) THEN
            ARG4=1.D0-0.5D0*(2.D0-OMEGH)**2
          ELSE
            ARG4=1.D0
          ENDIF
          SPEC(JF)=ARG1*ARG2*ARG3*ARG4
          IF (SPEC(JF).LT.E2FMIN) SPEC(JF)=0.D0
  100   CONTINUE
      ELSE
        DO 150 JF=1,NF
          SPEC(JF)=0.D0
  150   CONTINUE
      ENDIF
C
      RETURN
      END
