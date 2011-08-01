C                       *****************
                        SUBROUTINE SPREAD
C                       *****************
C
     *( DIRSPR, F     , COSTET, SINTET, NPLAN , FREQ  , DFREQ , NF    ,
     *  NPOIN2, TAILF , COSMOY, SINMOY, VARIAN, TAUXC , TAUXS , TAUXE )
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  SPREAD :  CALCUL DE L'ETALEMENT ANGULAIRE MOYEN SUR LE SPECTRE     C
C  ********  (DIRECTIONAL WIDTH) EN DEGRES.                           C
C                                                                     C
C   - CREE POUR VERSION 1.1  LE 28/12/95 PAR M. BENOIT                C
C   - MOD. POUR VERSION 1.2  LE 05/07/96 PAR M. BENOIT                C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! DIRSPR(-)   !<-- ! VECTEUR DES ETALEMENTS DIRECTIONNELS       !  C
C  ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL DE VARIANCE           !  C
C  ! COSTET(-)   ! -->! VECTEUR DES COSINUS DES DIRECTIONS         !  C
C  ! SINTET(-)   ! -->! VECTEUR DES SINUS   DES DIRECTIONS         !  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !  C
C  ! DFREQ(-)    ! -->! VECTEUR DES PAS DE FREQUENCE               !  C
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !  C
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL 2D    !  C
C  ! TAILF       ! -->! FACTEUR DE QUEUE DU SPECTRE                !  C
C  ! COSMOY(-)   !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! SINMOY(-)   !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! VARIAN(-)   !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUXC(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUXS(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUXE(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  DUMP2D                     C
C  ********    - PROGRAMME(S) APPELE(S) :                             C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NF    , NPLAN , NPOIN2
      DOUBLE PRECISION TAILF
      DOUBLE PRECISION DIRSPR(NPOIN2), SINMOY(NPOIN2), COSMOY(NPOIN2)
      DOUBLE PRECISION TAUXS (NPOIN2), TAUXC (NPOIN2), TAUXE (NPOIN2)
      DOUBLE PRECISION COSTET(NPLAN) , SINTET(NPLAN)
      DOUBLE PRECISION FREQ(NF), DFREQ(NF)
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF), VARIAN(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  IP    , JP    , JF
      DOUBLE PRECISION AUXC  , AUXS  , DEUPI , DFDTET, DTETAR, AUXI
      DOUBLE PRECISION CNVD  , SEUIL , COEFT
C
C
      SEUIL=1.D-20
      DEUPI=2.D0*3.14159265D0
      CNVD =360.D0/DEUPI
      DTETAR=DEUPI/DBLE(NPLAN)
C
      DO 10 IP=1,NPOIN2
        COSMOY(IP)=0.D0
        SINMOY(IP)=0.D0
        VARIAN(IP)=0.D0
   10 CONTINUE
C
C-----C-------------------------------------------------------C
C-----C  SOMMATION SUR LA PARTIE DISCRETISEE DU SPECTRE       C
C-----C-------------------------------------------------------C
      DO 30 JF=1,NF
C
        DFDTET=DFREQ(JF)*DTETAR
C
        DO 35 IP=1,NPOIN2
          TAUXC(IP)=0.D0
          TAUXS(IP)=0.D0
          TAUXE(IP)=0.D0
   35   CONTINUE
C
        DO 20 JP=1,NPLAN
          AUXC=COSTET(JP)*DFDTET
          AUXS=SINTET(JP)*DFDTET
          DO 40 IP=1,NPOIN2
            TAUXC(IP)=TAUXC(IP)+F(IP,JP,JF)*AUXC
            TAUXS(IP)=TAUXS(IP)+F(IP,JP,JF)*AUXS
            TAUXE(IP)=TAUXE(IP)+F(IP,JP,JF)*DFDTET
   40     CONTINUE
   20   CONTINUE
C
        DO 45 IP=1,NPOIN2
          COSMOY(IP)=COSMOY(IP)+TAUXC(IP)
          SINMOY(IP)=SINMOY(IP)+TAUXS(IP)
          VARIAN(IP)=VARIAN(IP)+TAUXE(IP)
   45   CONTINUE
C
   30 CONTINUE
C
C-----C-------------------------------------------------------------C
C-----C  PRISE EN COMPTE EVENTUELLE DE LA PARTIE HAUTES-FREQUENCES  C
C-----C-------------------------------------------------------------C
      IF (TAILF.GT.1.D0) THEN
        COEFT=FREQ(NF)/((TAILF-1.D0)*DFREQ(NF))
        DO 55 IP=1,NPOIN2
          COSMOY(IP)=COSMOY(IP)+TAUXC(IP)*COEFT
          SINMOY(IP)=SINMOY(IP)+TAUXS(IP)*COEFT
          VARIAN(IP)=VARIAN(IP)+TAUXE(IP)*COEFT
   55   CONTINUE
      ENDIF
C
C-----C-------------------------------------------------------------C
C-----C  CALCUL DE L'ETALEMENT ANGULAIRE PROPREMENT DIT             C
C-----C-------------------------------------------------------------C
      DO 60 IP=1,NPOIN2
        IF (VARIAN(IP).GT.SEUIL) THEN
          AUXS=SINMOY(IP)/VARIAN(IP)
          AUXC=COSMOY(IP)/VARIAN(IP)
          AUXI=MIN(DSQRT(AUXS*AUXS+AUXC*AUXC),1.D0)
          DIRSPR(IP)=DSQRT(2.D0*(1.D0-AUXI))*CNVD
        ELSE
          DIRSPR(IP)=SEUIL
        ENDIF
   60 CONTINUE
C
      RETURN
      END
