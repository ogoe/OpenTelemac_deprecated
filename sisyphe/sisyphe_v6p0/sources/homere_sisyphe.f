C                       **********************
                        PROGRAM HOMERE_SISYPHE
C                       **********************
C
C***********************************************************************
C  SISYPHE VERSION 6.0    09/04/2009    C LE NORMANT (LNH) 30 87 78 54
C
C***********************************************************************
C
C     FONCTIONS:
C     ==========
C
C 1)  ACQUISITION DE TOUTES LES DONNEES NECESSAIRES
C     AU CALCUL DES POINTEURS: FICHIER CAS + PARTIELLEMENT LA GEOMETRIE
C
C 2)  APPEL DU SOUS-PROGRAMME SISYPHE.
C
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PRINCI
C
C SOUS-PROGRAMMES APPELES : LECDON , POINT , SISYPHE
C
C**********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_SISYPHE
      USE INTERFACE_SISYPHE
C
      IMPLICIT NONE
      INTEGER     LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER TDEB,TFIN,IFLOT,NCAR,DUMINT
      LOGICAL DUMLOG
C
      CHARACTER(LEN=24), PARAMETER :: CODE='SISYPHE                 '
C
      CHARACTER(LEN=250) PATH
      CHARACTER(LEN=144) MOTCAR(300),FILE_DESC(4,300)
C
      INTEGER  TIME_IN_SECONDS
      EXTERNAL TIME_IN_SECONDS
C
C======================================================================
C
      CALL BIEF_INIT(CODE,PATH,NCAR,.TRUE.)
c
      TDEB = TIME_IN_SECONDS()
C
C  ENTETE SUR LISTING
C
      IF(LNG.EQ.1) WRITE(LU,100)
      IF(LNG.EQ.2) WRITE(LU,101)
      WRITE(LU,102)
100   FORMAT(/////,1X,'LISTING DE SISYPHE ',78('-'))
101   FORMAT(/////,1X,'LISTING OF SISYPHE ',78('-'))
102   FORMAT(/////,
     *14X,'    SSSS I   SSSS Y   Y PPPP  H   H EEEEE',/,
     *14X,'   S     I  S      Y Y  P   P H   H E    ',/,
     *14X,'    SSS  I   SSS    Y   PPPP  HHHHH EEEE  ',/,
     *14X,'       S I      S   Y   P     H   H E     ',/,
     *14X,'   SSSS  I  SSSS    Y   P     H   H EEEEE',/,
     *14X,'                                          ',/,
     *14X,'                 VERSION 6.0              ',/,
     *14X,/////)
C
C-----------------------------------------------------------------------
C
C LECTURE DU FICHIER CAS
C
      CALL LECDON_SISYPHE(MOTCAR,FILE_DESC,PATH,NCAR,CODE)
C
C-----------------------------------------------------------------------
C
      IFLOT = 0
      CALL BIEF_OPEN_FILES(CODE,SIS_FILES,MAXLU_SIS,
     *                     PATH,NCAR,.FALSE.,IFLOT,2) 
C
C-----------------------------------------------------------------------
C
C ALLOCATION DES VECTEURS, MATRICES ET BLOCS
C
      CALL POINT_SISYPHE
C
C-----------------------------------------------------------------------
C
C  APPEL DU SOUS-PROGRAMME TELEMAC
C  -1 VEUT DIRE QU'ON PASSE DANS TOUTE LA SUBROUTINE  SISYPHE CAR ON NE
C  COUPLE PAS. LES AUTRES VARIABLES NE SERVENT QUE POUR LE COUPLAGE
C
C     INOUT VARIABLES IN SISYPHE CANNOT BE HARDCODED
      DUMINT=1
      DUMLOG=.FALSE.
C
      CALL SISYPHE(-1,0,0,0,0,T1,T1,T1,T1,T1,T1,T1,
     *             DUMLOG,DUMINT,DUMLOG,CODE,1,
     *             T1,T1,0.D0,T1,0.D0,DUMLOG,DUMLOG,
     *             T1,1,T1,T1,T1,T1)
C
C-----------------------------------------------------------------------
C
      CALL BIEF_CLOSE_FILES(CODE,SIS_FILES,MAXLU_SIS,.TRUE.)
C
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1) WRITE(LU,10)
      IF(LNG.EQ.2) WRITE(LU,11)
10    FORMAT(1X,///,1X,'FIN NORMALE DU PROGRAMME',///)
11    FORMAT(1X,///,1X,'CORRECT END OF RUN',///)
C
C-----------------------------------------------------------------------
C
      TFIN = TIME_IN_SECONDS()
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'DUREE DU CALCUL : ',TFIN-TDEB,' SECONDES'
      ENDIF
      IF(LNG.EQ.2) THEN
        WRITE(LU,*) 'COMPUTER TIME: ',TFIN-TDEB,' SECONDS'
      ENDIF
C
C-----------------------------------------------------------------------
C
      STOP
      END PROGRAM HOMERE_SISYPHE
