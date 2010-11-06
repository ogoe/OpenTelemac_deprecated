      ! ******************** !
        SUBROUTINE SIS_ARRET ! 
      ! ******************** !

     &(ESM,EMAX,HN,VARSOR,NPOIN,MN,NRES,FMTRES,MAXVAR,AT0,RC,HIST,
     & BINRESSIS,TEXTE,SORLEO,SORIMP,T1,T2)

C**********************************************************************C
C SISYPHE VERSION 5.8  05/01/04  F.    HUVELIN                         C
C SISYPHE VERSION 5.1  11/09/95  E.    PELTIER                         C
C SISYPHE VERSION 5.1  11/09/95  C.    LENORMANT                       C
C SISYPHE VERSION 5.1  11/09/95  J.-M. HERVOUET                        C
C**********************************************************************C


             ! ======================================= !
             ! Stop test in case of too high evolution !
             ! ======================================= !


C**********************************************************************C
C                                                                      C
C                 SSSS I   SSSS Y   Y PPPP  H   H EEEEE                C
C                S     I  S      Y Y  P   P H   H E                    C
C                 SSS  I   SSS    Y   PPPP  HHHHH EEEE                 C
C                    S I      S   Y   P     H   H E                    C
C                SSSS  I  SSSS    Y   P     H   H EEEEE                C
C                                                                      C
C----------------------------------------------------------------------C
C                             ARGUMENTS                                C
C .________________.____.______________________________________________C
C |      NOM       |MODE|                   ROLE                       C
C |________________|____|______________________________________________C
C |________________|____|______________________________________________C
C                    <=  Can't be change by the user                   C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY SISYPHE                                                    !
!                                                                      !
! CALL DESIMP, PREDES, PLANTE
! CALL OS    
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE INTERFACE_SISYPHE,EX_SIS_ARRET => SIS_ARRET
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ),    INTENT(IN)    :: ESM, EMAX, HN, VARSOR
      INTEGER,           INTENT(IN)    :: NPOIN, MN, NRES, MAXVAR
      DOUBLE PRECISION,  INTENT(IN)    :: AT0, RC, HIST(1)
      CHARACTER(LEN=3),  INTENT(IN)    :: BINRESSIS
      CHARACTER(LEN=32), INTENT(IN)    :: TEXTE(MAXVAR)
      CHARACTER(LEN=8),  INTENT(IN)    :: FMTRES
      LOGICAL,           INTENT(IN)    :: SORLEO(*), SORIMP(*)
      TYPE(BIEF_OBJ),    INTENT(INOUT) :: T1, T2
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER          :: IMIN
      DOUBLE PRECISION :: XMIN
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!

      ! ************************************ !
      ! I - CALCUL DE LA VALEUR ABSOLUE DE F ! 
      ! ************************************ !
      CALL OS('X=ABS(Y)', X=T1, Y=ESM)

      ! ************************************************** !
      ! II - CALCUL DE LA DIFFERENCE ENTRE LA VALEUR       ! (_IMP_)
      !      ADMISSIBLE DU DEPOT ET LA VALEUR ABSOLUE DE F ! 
      ! ************************************************** !
      CALL OS('X=Y-Z   ', X=T2, Y=EMAX, Z=T1)

      ! ************************************************************** !
      ! III - CALCUL DE LA VALEUR MINIMALE DE LA DIFFERENCE PRECEDENTE ! 
      ! ************************************************************** !
      CALL MINI(XMIN, IMIN, T2%R, NPOIN)

      ! ************************************************************* !
      ! IV - SI LA VALEUR MINIMALE EST NEGATIVE, LE CALCUL EST ARRETE ! 
      ! ************************************************************* !
      IF (XMIN < 0.D0) THEN

         ! IV.1 - Printing the values 
         ! --------------------------
         IF(LNG.EQ.1) THEN
            WRITE(LU,400) MN
            WRITE(LU,*) ' '
            WRITE(LU,402) IMIN
            WRITE(LU,404) HN%R(IMIN)
            WRITE(LU,406) RC
            WRITE(LU,408) EMAX%R(IMIN)
            WRITE(LU,410) ESM%R(IMIN)
            WRITE(LU,412) AT0
            WRITE(LU,*) ' '
            WRITE(LU,*) 'DERNIER RESULTAT SAUVEGARDE'
         ELSE IF(LNG.EQ.2) THEN
            WRITE(LU,401) MN
            WRITE(LU,*) ' '
            WRITE(LU,403) IMIN
            WRITE(LU,405) HN%R(IMIN)
            WRITE(LU,407) RC
            WRITE(LU,409) EMAX%R(IMIN)
            WRITE(LU,411) ESM%R(IMIN)
            WRITE(LU,413) AT0
            WRITE(LU,*) ' '
            WRITE(LU,*) 'LAST RESULT SAVED'
         ENDIF


         ! IV.2 - Saving the last result
         ! -----------------------------
         CALL PREDES(1,AT0)
         CALL BIEF_DESIMP(FMTRES,VARSOR,HIST,0,NPOIN,NRES,BINRESSIS,AT0,
     &                    1,1,1,SORLEO,SORIMP,MAXVAR,TEXTE,1,1)
         CALL PLANTE(1)
         STOP
      ENDIF


      !----------------------------------------------------------------!
400   FORMAT(1X,/,' EVOLUTION TROP FORTE AU CALCUL  : ',1I6)
402   FORMAT(' NOEUD NUMERO                    : ',1I6)
404   FORMAT(' HAUTEUR D''EAU                  : ',G16.7)
406   FORMAT(' RAPPORT D''EVOLUTION CRITIQUE   : ',G16.7)
408   FORMAT(' EVOLUTION MAXIMALE ADMISSIBLE   : ',G16.7)
410   FORMAT(' EVOLUTION CUMULEE CALCULEE      : ',G16.7)
412   FORMAT(' TEMPS                           : ',G16.7)
      !----------------------------------------------------------------!
401   FORMAT(1X,/,' TOO MUCH EVOLUTION AT COMPUTATION: ',1I6)
403   FORMAT(' NODE NUMBER                     : ',1I6)
405   FORMAT(' WATER DEPTH                     : ',G16.7)
407   FORMAT(' CRITICAL EVOLUTION RATIO        : ',G16.7)
409   FORMAT(' MAXIMAL ALLOWED EVOLUTION       : ',G16.7)
411   FORMAT(' COMPUTED CUMULATED EVOLUTION    : ',G16.7)
413   FORMAT(' TIME                            : ',G16.7)
      !----------------------------------------------------------------!

!======================================================================!
!======================================================================!

      RETURN
      END

