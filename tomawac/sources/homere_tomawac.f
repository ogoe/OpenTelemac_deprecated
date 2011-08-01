C                       **********************
                        PROGRAM HOMERE_TOMAWAC
C                       **********************
C
C**********************************************************************
C  TOMAWAC VERSION 6.0       12/01/01           OPTIMER  02 98 44 24 51
C**********************************************************************
C
C                       PROGRAMME PRINCIPAL DE
C
C           TTTTT  OOOOO  M   M  AAAAA  W   W  AAAAA  CCCC 
C             T    O   O  MM MM  A   A  W   W  A   A  C    
C             T    O   O  M M M  AAAAA  W W W  AAAAA  C    
C             T    O   O  M   M  A   A  WW WW  A   A  C    
C             T    OOOOO  M   M  A   A  W   W  A   A  CCCC
C
C                    RESOLUTION DE l'EQUATION DE
C                 DU SPECTRE DE HOULE DIRECTIONNELLE
C
C-----------------------------------------------------------------------
C
C    1) LIT LES DONNEES NECESSAIRES POUR L ALLOCATION DE MEMOIRE
C       READING ALL THE NECESSARY DATA FOR ALLOCATING MEMORY
C    2) ALLOUE LA MEMOIRE  - ALLOCATING MEMORY
C    3) APPELLE LE VERITABLE PROGRAMME PRINCIPAL WAC
C       CALLING THE REAL MAIN PROGRAM WAC
C
C**********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TOMAWAC
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
      INTEGER NCAR,IFLOT
C
      CHARACTER(LEN=24), PARAMETER :: CODE='TOMAWAC                '
C
      CHARACTER(LEN=250) PATH
      CHARACTER(LEN=144) FILE_DESC(4,300)
C
      CALL BIEF_INIT(CODE,PATH,NCAR,.TRUE.)
C
C-----------------------------------------------------------------------
C     EN TETE   -  HEADING
C-----------------------------------------------------------------------
C
      WRITE(LU,100)
      WRITE(LU,110)
100   FORMAT(100(1H-),////////,
     *16X,
     *'TTTTT  OOOOO  M   M  AAAAA  W   W  AAAAA  CCCCC '
     *,/,16X,
     *'  T    O   O  MM MM  A   A  W   W  A   A  C     '
     *,/,16X,
     *'  T    O   O  M W M  AAAAA  W W W  AAAAA  C     '
     *,/,16X,
     *'  T    O   O  M   M  A   A  WW WW  A   A  C     '
     *,/,16X,
     *'  T    OOOOO  M   M  A   A  W   W  A   A  CCCCC '
     *,//)
110   FORMAT(15X,
     *'               |    |    |                 '
     *,/,15X,
     *'              )_)  )_)  )_) _              '
     *,/,15X,
     *'             )___))___))___)\              '
     *,/,15X,
     *'             )____)____)_____)\\           '
     *,/,15X,
     *'           _____|____|____|____\\\__       '
     *,/,15X,
     *'  ---------\               6.0  /---------  '
     *,/,15X,
     *'    ^^^^^^^^^^^^^^^^^^^^^^^^^^^             '
     *,/,15X,
     *'         ^^^^      ^^^^     ^^^    ^^      '
     *,/,15X,
     *'             ^^^^      ^^^                 '
     *,///)
C
C-----------------------------------------------------------------------
C     LECTURE DU FICHIER DES PARAMETRES
C     READING THE PARAMETER FILE
C
      CALL LECDON_TOMAWAC(FILE_DESC,PATH,NCAR,CODE)
C
C-----------------------------------------------------------------------
C     OUVERTURE DES FICHIERS
C     OPENING FILES
C
      IFLOT = 0
      CALL BIEF_OPEN_FILES(CODE,WAC_FILES,44,PATH,NCAR,.FALSE.,IFLOT,1)
C
C-----------------------------------------------------------------------
C     ALLOCATION DE LA MEMOIRE
C     ALLOCATING MEMORY
C
      CALL POINT_TOMAWAC
C
C-----------------------------------------------------------------------
C     APPEL DU VRAI PROGRAMME PRINCIPAL
C     CALLING THE REAL MAIN PROGRAM
C
      CALL WAC
C
C-----------------------------------------------------------------------
C     FERMETURE DES FICHIERS
C     CLOSING THE FILES
C
      CALL BIEF_CLOSE_FILES(CODE,WAC_FILES,44,.TRUE.)
C
C-----------------------------------------------------------------------
C      
      IF (LNG.EQ.1) WRITE(LU,10)
      IF (LNG.EQ.2) WRITE(LU,20)
10    FORMAT(1X,////,1X,'FIN NORMALE DU PROGRAMME',/////)
20    FORMAT(1X,////,1X,'CORRECT END OF RUN',/////)
C
C-----------------------------------------------------------------------
C
      STOP
      END
