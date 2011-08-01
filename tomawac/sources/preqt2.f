C                       *****************
                        SUBROUTINE PREQT2
C                       *****************
C
     *( TETA  , NPLAN , BDISPB, BDSSPB, NBD , INDI ) 
C
C**********************************************************************
C  TOMAWAC - V5P0                           (EDF/DER/LNH)  -  11/06/98
C**********************************************************************
C
C  FONCTION : TERME SOURCE LIE AUX INTERACTIONS NON-LINEAIRES ENTRE
C  ********** TRIPLETS DE FREQUENCES.
C             MODELE ISSU DES EQUATIONS DE BOUSSINESQ.
C
C  ARGUMENTS :
C  ***********
C  +-------------+----+--------------------------------------------+
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C  +-------------+----+--------------------------------------------+
C  ! F(-,-,-)    !<-->! SPECTRE DIRECTIONNEL DE VARIANCE           !
C  ! XK(-,-)     ! -->! TABLEAU DES NOMBRES D'ONDE                 !
C  ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !
C  ! DFREQ(-)    ! -->! TABLEAU DES PAS DE FREQUENCE               !
C  ! DEPTH(-)    ! -->! TABLEAU DES PROFONDEURS (METRES)           !
C  ! TETA(-)     ! -->! VECTEUR DES DIRECTIONS DE DISCRETISATION   !
C  ! SINTET(-)   ! -->! VECTEUR DES   SINUS DES DIRECTIONS         !
C  ! COSTET(-)   ! -->! VECTEUR DES COSINUS DES DIRECTIONS         !
C  ! RAISF       ! -->! RAISON FREQUENTIELLE POUR DISCRETISATION   !
C  ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C  ! TSTOT(-,-,-)!<-- ! CONTRIBUTION TERME SOURCE - PARTIE TOTALE  !
C  ! TSDER(-,-,-)!<-- ! CONTRIBUTION TERME SOURCE - PARTIE DERIVEE !
C  +-------------+----+--------------------------------------------+
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C  +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  SEMIMP
C  ********    - PROGRAMME(S) APPELE(S) :  
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
      INTEGER           NPLAN
      DOUBLE PRECISION  BDISPB , BDSSPB
      DOUBLE PRECISION  TETA(NPLAN)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER           IPL
      DOUBLE PRECISION  AP2 , EPS , DTETA
C
      INTEGER           NBPL , NBPU, NBD, NB1
      INTEGER           INDI(NPLAN)
C
C     """"""""""""""""""""""""""""""""""""""""""""""""""""""""""
C
      DTETA = TETA(2)-TETA(1)
      EPS  = 1.D-5
      IF(BDSSPB.GE.BDISPB) THEN
         AP2  = (BDISPB-TETA(1))/DTETA
         NBPL = IDINT(AP2)
         AP2  = AP2 - DBLE(NBPL)
         IF(AP2.GT.EPS) THEN
            NBPL = NBPL + 2
         ELSE
            NBPL = NBPL + 1
         ENDIF
         AP2  = (BDSSPB-TETA(1))/DTETA
         NBPU = IDINT(AP2) + 1
         NBD=NBPU-NBPL+1
C         ALLOCATE(INDI(1:NBD))
         DO IPL=1,NBD
            INDI(IPL)=NBPL+IPL-1
         END DO
      ELSE
         AP2  = (BDSSPB-TETA(1))/DTETA
         NBPU = IDINT(AP2) + 1
         AP2  = (BDISPB-TETA(1))/DTETA
         NBPL = IDINT(AP2)
         AP2  = AP2 - DBLE(NBPL)
         IF(AP2.GT.EPS) THEN
            NBPL = NBPL + 2
         ELSE
            NBPL = NBPL + 1
         ENDIF
         IF(NBPL.GT.NPLAN) THEN
            NBPL = 1
            INDI(1) = 1
            NBD  = NBPU - NBPL + 1
C            ALLOCATE(INDI(1:NBD))
            DO IPL = 2,NBD
               INDI(IPL)=IPL
            END DO
         ELSE
            NB1 = NPLAN - NBPL + 1
            NBD = NB1 + NBPU
C            ALLOCATE(INDI(1:NBD))
            DO IPL = 1,NB1
               INDI(IPL)=NBPL+IPL-1
            END DO
            DO IPL = 1,NBPU
               INDI(IPL+NB1)=IPL
            END DO
         ENDIF
      ENDIF
C
      RETURN
      END
