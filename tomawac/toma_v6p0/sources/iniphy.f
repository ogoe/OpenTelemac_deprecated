C                       *****************
                        SUBROUTINE INIPHY
C                       *****************
C
     *( XK    , CG    , B     , DEPTH , FREQ  , COSPHI, NPOIN2, NF    ,
     *  PROINF, SPHE  )
C
C**********************************************************************
C  TOMAWAC - V1.0    M. BENOIT               (EDF/DER/LNH)  -  07/02/95
C**********************************************************************
C
C  FONCTION : CALCUL DES GRANDEURS DE HOULE INDEPENDANTES DU TEMPS
C  ********** (NOMBRE D'ONDE, VITESSE DE GROUPE,...)
C
C  ARGUMENTS :
C  ***********
C  +-------------+----+--------------------------------------------+
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C  +-------------+----+--------------------------------------------+
C  ! XK(-,-)     !<-- ! MATRICE DES NOMBRES D'ONDE EN FR           !
C  ! CG(-,-)     !<-- ! MATRICE DE VITESSES DE GROUPE EN FR        !
C  ! B(-,-)      !<-- ! MATRICE DE PASSAGE DE (KX,KY) A (FR,TETA)  !
C  ! DEPTH(-)    ! -->! VECTEUR DES PROFONDEURS                    !
C  ! FREQ(-)     ! -->! VECTEUR DES FREQUENCES DE DISCRETISATION   !
C  ! COSPHI(-)   ! -->! VECTEUR DES COSINUS DES LATITUDES          !
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL 2D    !
C  ! NF          ! -->! NOMBRE DE FREQUENCES                       !
C  ! PROINF      ! -->! INDICATEUR CALCUL EN PROFONDEUR INFINIE    !
C  ! SPHE        ! -->! INDICATEUR CALCUL EN SPHERIQUE             !
C  +-------------+----+--------------------------------------------+
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C  +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  : WAC
C  ********    - PROGRAMME(S) APPELE(S) : WNSCOU
C
C  REMARQUES :
C  ***********
C   - TOUTES LES DIRECTIONS SONT EN RADIAN ET COMPRISES ENTRE 0 - 2PI
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER          NF    , NPOIN2
      DOUBLE PRECISION DEPTH(NPOIN2)    , COSPHI(NPOIN2), FREQ(NF)
      DOUBLE PRECISION B(NPOIN2,NF)  , XK(NPOIN2,NF) , CG(NPOIN2,NF)
      LOGICAL          PROINF, SPHE
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER          IP    , JF
      DOUBLE PRECISION DEUPI2, DEUPI , XG    , DPDSUG, AUX2
      DOUBLE PRECISION AUX1  , AUX3  , DEUKD , R2
C
C
      XG=9.81D0
      DEUPI=2.D0*3.14159265D0
      DEUPI2=DEUPI*DEUPI
      DPDSUG=DEUPI2/XG
      R2=(6400.D3)**2
C
      IF (PROINF) THEN
C                               +---------------+
C.....................CALCUL EN ! PROF. INFINIE !
C                               +---------------+
        DO 310 JF=1,NF
          AUX1=DPDSUG*(FREQ(JF))**2
          AUX3=0.5D0*XG/(DEUPI*FREQ(JF))
          DO 320 IP=1,NPOIN2
            XK(IP,JF)=AUX1
            CG(IP,JF)=AUX3
  320     CONTINUE
  310   CONTINUE
      ELSE
C                               +---------------+
C.....................CALCUL EN ! PROF.   FINIE !
C                               +---------------+
        DO 410 JF=1,NF
          AUX2=DEUPI*FREQ(JF)
          DO 430 IP=1,NPOIN2
            CALL WNSCOU(AUX1,FREQ(JF),DEPTH(IP))
            DEUKD=2.D0*AUX1*DEPTH(IP)
            IF (DEUKD.GT.7.D2) THEN
              AUX3=0.5D0*AUX2/AUX1
            ELSE
              AUX3=0.5D0*(1.D0+DEUKD/SINH(DEUKD))*AUX2/AUX1
            ENDIF
            XK(IP,JF)=AUX1
            CG(IP,JF)=AUX3
  430     CONTINUE
  410   CONTINUE
      ENDIF
C
C
C.....CALCUL DE  B POUR LE PASSAGE DE (KX,KY) A (FR,TETA)
C     ===================================================
      IF (.NOT.SPHE) THEN
C                               +-----------+
C.....................CALCUL EN ! CARTESIEN !
C                               +-----------+
        DO 710 JF=1,NF
          AUX1=DEUPI2*FREQ(JF)
          DO 720 IP=1,NPOIN2
            B(IP,JF)= CG(IP,JF)/(AUX1*XK(IP,JF))
  720     CONTINUE
  710   CONTINUE
C
      ELSE
C                               +-----------+
C.....................CALCUL EN ! SPHERIQUE !
C                               +-----------+
        DO 810 JF=1,NF
          AUX1=DEUPI2*FREQ(JF)*R2
          DO 820 IP=1,NPOIN2
            B(IP,JF)= CG(IP,JF)/(AUX1*XK(IP,JF)*COSPHI(IP))
  820     CONTINUE
  810   CONTINUE
      ENDIF
C
      RETURN
      END
