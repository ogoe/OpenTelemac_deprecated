C                       **********************
                        PROGRAM HOMERE_ARTEMIS
C                       **********************
C
C***********************************************************************
C
C  ARTEMIS VERSION 6.0       19/04/99  D. AELBRECHT (LNH) 01 30 87 74 12
C
C  LINKED TO BIEF VERS. 5.5  24/04/97  J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C     FUNCTIONS:
C     ==========
C
C 1)  ACQUISITION OF DATA REQUIRED FOR ALL ALLOCATING MEMORY
C     (PARAMETER FILE + GEOMETRY)
C
C 2)  CALL TO THE REAL MAIN PROGRAM ARTEMIS 
C
C
C**********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_ARTEMIS
C
      IMPLICIT NONE
      INTEGER     LNG,LU,I,J
      COMMON/INFO/LNG,LU
C
      INTEGER TDEB,TFIN,NCAR,IFLOT
C
      CHARACTER(LEN=24), PARAMETER :: CODE='ARTEMIS                 '
C
      INTEGER  TIME_IN_SECONDS
      EXTERNAL TIME_IN_SECONDS
C
C======================================================================
C
      CHARACTER(LEN=250) PATH
      CHARACTER(LEN=144) FILE_DESC(4,300)
C
C-----------------------------------------------------------------------
C      
      CALL BIEF_INIT(CODE,PATH,NCAR,.TRUE.)
 
C
      TDEB = TIME_IN_SECONDS()
C
C     HEADING
C
      IF(LNG.EQ.1) WRITE(LU,100)
      IF(LNG.EQ.2) WRITE(LU,101)
      WRITE(LU,102)
100   FORMAT(/////,1X,'LISTING DE ARTEMIS ',78('-'))
101   FORMAT(/////,1X,'LISTING OF ARTEMIS ',78('-'))
102   FORMAT(/////,
     *14X,'    AAA  RRRR  TTTTT EEEEE M   M IIIII  SSSS',/,
     *14X,'   A   A R   R   T   E     MM MM   I   S    ',/,
     *14X,'   AAAAA RRRR    T   EEEEE M M M   I    SSS ',/,
     *14X,'   A   A R   R   T   E     M   M   I       S',/,
     *14X,'   A   A R   R   T   EEEEE M   M IIIII SSSS ',/,
     *14X,'                                            ',/,
     *14X,'          VERSION 6.0      FORTRAN 90 ',/,
     *14X,/////)
C
C-----------------------------------------------------------------------
C
C     READING THE PARAMETER FILE
C
      CALL LECDON_ARTEMIS(FILE_DESC,PATH,NCAR,CODE)
  
C-----------------------------------------------------------------------
C
C     OPENING FILES
C
      IFLOT = 0
      CALL BIEF_OPEN_FILES(CODE,ART_FILES,44,PATH,NCAR,.FALSE.,IFLOT,1)
      
C-----------------------------------------------------------------------
C
C     ALLOCATING MEMORY IN BIEF_OBJ STRUCTURES (VECTORS, MATRICES)
     

C
      CALL POINT_ARTEMIS
C
C-----------------------------------------------------------------------
C
C     CALL TO REAL MAIN PROGRAM
C
      CALL ARTEMIS
C
C-----------------------------------------------------------------------
C
      CALL BIEF_CLOSE_FILES(CODE,ART_FILES,44,.TRUE.)
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
      END
