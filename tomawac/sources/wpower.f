C                       *****************
                        SUBROUTINE WPOWER
C                       *****************
C
     *( POWER , F     , FREQ  , DFREQ , CG    , TAILF , NF    , NPLAN ,
     *  NPOIN2, ROEAU , GRAVIT)
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 5.3 ----                         C
C                                                                     C
C  TOTNRJ : CALCUL DE LA VARIANCE DU SPECTRE DIRECTIONNEL EN TOUS     C
C  ******** LES POINTS DU MAILLAGE SPATIAL 2D. CETTE VARIANCE EST     C
C           OBTENUE PAR INTEGRATION SUR LES FREQUENCES ET DIRECTIONS  C
C           ET PREND EVENTUELLEMENT EN COMPTE LA PARTIE HAUTES        C
C           FREQUENCES DU SPECTRE.                                    C
C                                                                     C
C   - CREE POUR VERSION 5.3  LE 20/05/2003 PAR M. BENOIT              C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! POWER(-)    !<-- ! TABLEAU DE PUISSANCE LINEIQUE              !  C
C  ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL DE VARIANCE           !  C
C  ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !  C
C  ! DFREQ(-)    ! -->! TABLEAU DES PAS DE FREQUENCE               !  C
C  ! CG(-,-)     ! -->! TABLEAU DES VITESSES DE GROUPE             !  C
C  ! TAILF       ! -->! FACTEUR DE QUEUE DU SPECTRE                !  C
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL 2D    !  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  DUMP2D, TERSOU             C
C  ********    - PROGRAMME(S) APPELE(S) :                             C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  REMARQUES :                                                        C
C  ***********                                                        C
C  - LA PARTIE HAUTES-FREQUENCES DU SPECTRE N'EST PRISE EN COMPTE QUE C
C    SI LE FACTEUR DE QUEUE (TAILF) EST STRICTEMENT SUPERIEUR A 1.    C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER          NF    , NPLAN , NPOIN2
      DOUBLE PRECISION TAILF , POWER(NPOIN2), FREQ(NF), DFREQ(NF)
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF), CG(NPOIN2,NF)
      DOUBLE PRECISION GRAVIT, ROEAU
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER          IP    , JP    , JF
      DOUBLE PRECISION AUX1  , DTETAR, ROGER
C
C
      DTETAR=2.D0*3.14159265D0/DBLE(NPLAN)
      ROGER=ROEAU*GRAVIT/1000.D0
      DO IP=1,NPOIN2
        POWER(IP)=0.D0
      ENDDO
C
C-----C-------------------------------------------------------C
C-----C  SOMMATION SUR LA PARTIE DISCRETISEE DU SPECTRE       C
C-----C-------------------------------------------------------C
      DO JF=1,NF
        AUX1=DFREQ(JF)*DTETAR
        DO JP=1,NPLAN
          DO IP=1,NPOIN2
            POWER(IP) = POWER(IP) + F(IP,JP,JF)*CG(IP,JF)*AUX1
          ENDDO
        ENDDO
      ENDDO
C
C-----C-------------------------------------------------------------C
C-----C  PRISE EN COMPTE EVENTUELLE DE LA PARTIE HAUTES-FREQUENCES  C
C-----C-------------------------------------------------------------C
      IF (TAILF.GT.1.D0) THEN
        AUX1=DTETAR*GRAVIT/(4.D0*3.14159265D0*TAILF)
        DO JP=1,NPLAN
          DO IP=1,NPOIN2
            POWER(IP)=POWER(IP) + F(IP,JP,NF)*AUX1
          ENDDO
        ENDDO
      ENDIF
C
C-----C-------------------------------------------------------------C
C-----C  PASSAGE EN KW/m  (en mulitpliant par ro.g/1000)            C
C-----C-------------------------------------------------------------C
      DO IP=1,NPOIN2
        POWER(IP)=POWER(IP)*ROGER
      ENDDO
C
      RETURN
      END
