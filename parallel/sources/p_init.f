C                       *****************
                        SUBROUTINE P_INIT
C                       *****************
C
     *(CHAINE,NCAR,IPID,NCSIZE)
C
C***********************************************************************
C  PARALLEL  VERSION 5.9       /06/96           HANS HERRMANN (HANOVRE)
C             MODIFIED        08/06/96     REINHARD HINKELMANN (HANOVRE)
C             MODIFIED        17/12/96            J-M HERVOUET (LNH)
C  ADAPTED FOR MPI              /10/99       RAINER JOHANNI (SGI MUNICH)
C  VERSION 5.0 MODIFIED       28/12/99    J.A. JANKOWSKI (BAW KARLSRUHE)
C        MODIFIED (SIZE OF EXTENSION)       16/05/08    P. VEZOLLE (IBM)
C***********************************************************************
C
C      FONCTIONS: INITIALISATIONS.
C      ==========
C
C      MELDET PROGRAMM BEI PARASTATION AN.
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   CHAINE       | <--| CHAINE DE CARACTERES CONTENANT LE NOM DU
C |                |    | REPERTOIRE DE TRAVAIL
C |   NCAR         | <--| NOMBRE DE LETTRES DE LA CHAINE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PRINCI
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON /INFO/ LNG,LU
C
      INCLUDE 'mpif.h'
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER, INTENT(OUT)            :: NCAR,IPID,NCSIZE
      CHARACTER(LEN=250), INTENT(OUT) :: CHAINE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER MYTID,IER,I,NPREAD
C
      CHARACTER*255 EXE
C
      LOGICAL YAPARA
      CHARACTER*11 PNUMBER
      CHARACTER*13 MYNAM
C
C-----------------------------------------------------------------------
C MPI IS SILENT WHEN EVERYTHING IS GOING ON PROPERLY
C I SET DEFAULT LANGUAGE 2 AND STANDARD OUTPUT 6
C IN ORDER TO SEE ERROR MESSAGES (TO THE MASTER!)
C THE SUBROUTINE CALLED NEXT IS READ_CONFIG !
C THIS IS NOT PRETTY...
C
      LNG=2
      LU=6
C
C     ALL WRITE STATEMENTS BEFORE OPENING A FILE ON LU=6 SEEM TO RAISE
C     A PROBLEM ON WINDOWS COMPAQ COMPILER
C
C     PROPOSITION DE SOGREAH (MAIS DENYNONE PAS DANS LA NORME)
C     OPEN(UNIT=LU,FILE="PARALLEL.LOG",FORM='FORMATTED',ACTION='WRITE',
C    *     SHARE='DENYNONE')
C
c$$$      IF(LNG.EQ.1) WRITE(LU,*) 'ENTREE DANS P_INIT'
c$$$      IF(LNG.EQ.2) WRITE(LU,*) 'ENTERING P_INIT'
c$$$      IF(LNG.EQ.1) WRITE(LU,*)
c$$$     *    'LOGGING OUTPUT DIRECTED INTO FILES FOR EACH PROCESSOR'
c$$$      IF(LNG.EQ.2) WRITE(LU,*)
c$$$     *    'LOGGING OUTPUT DIRECTED INTO FILES FOR EACH PROCESSOR'
C
C INITIALISE MPI
C
      CALL MPI_INIT(IER)
      IF(IER.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'P_INIT: ERREUR DANS MPI_INIT'
        IF(LNG.EQ.2) WRITE(LU,*) 'P_INIT: ERROR IN MPI_INIT'
        WRITE(LU,*) 'MPI ERROR ',IER
        STOP
      ENDIF
C
C OBTAIN MYTID, IT IS VIRTUALLY THE PROCESSOR NUMBER (RANK)
C
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IPID,IER)
      IF(IER.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'P_INIT: ERREUR DANS MPI_COMM_RANK'
        IF(LNG.EQ.2) WRITE(LU,*) 'P_INIT: ERROR IN MPI_COMM_RANK'
        WRITE(LU,*) 'MPI ERROR ',IER
        STOP
      ENDIF
C
C OBTAIN NCSIZE
C
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NCSIZE,IER)
      IF(IER.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'P_INIT: ERREUR DANS MPI_COMM_SIZE'
        IF(LNG.EQ.2) WRITE(LU,*) 'P_INIT: ERROR IN MPI_COMM_SIZE'
        WRITE(LU,*) 'MPI ERROR ',IER
        STOP
      ENDIF
C
c$$$      IF(LNG.EQ.1) WRITE(LU,*)
c$$$     &    'MPI: CALCUL AVEC ',NCSIZE,' PROCESSEURS'
c$$$      IF(LNG.EQ.2) WRITE(LU,*)
c$$$     &    'MPI: COMPUTATION WITH ',NCSIZE,' PROCESSORS'
C
C MANIPULATE MASTER'S AND SLAVES' STANDART OUTPUT TO FILES PE#LOG.TXT
C SLAVES WRITE TO CHANNEL 95 (?)
C WORKS FOR PE# 0-999 (WE HAVE A DECENT NUMBER OF PROCESSORS!)
C
      IF(NCSIZE.GT.1) THEN
        PNUMBER = '00000-00000'
C
        IF((NCSIZE-1).LT.10) THEN
          WRITE(PNUMBER(05:05),'(I1)') NCSIZE-1
        ELSEIF((NCSIZE-1).LT.100) THEN
          WRITE(PNUMBER(04:05),'(I2)') NCSIZE-1
        ELSEIF((NCSIZE-1).LT.1000) THEN
          WRITE(PNUMBER(03:05),'(I3)') NCSIZE-1
        ELSEIF((NCSIZE-1).LT.10000) THEN
          WRITE(PNUMBER(02:05),'(I4)') NCSIZE-1
        ELSE
          WRITE(PNUMBER(01:05),'(I5)') NCSIZE-1
        ENDIF
C
        IF(IPID.LT.10) THEN
          WRITE(PNUMBER(11:11),'(I1)') IPID
        ELSEIF(IPID.LT.100) THEN
          WRITE(PNUMBER(10:11),'(I2)') IPID
        ELSEIF(IPID.LT.1000) THEN
          WRITE(PNUMBER(9:11),'(I3)') IPID
        ELSEIF(IPID.LT.10000) THEN
          WRITE(PNUMBER(8:11),'(I4)') IPID
        ELSE
          WRITE(PNUMBER(7:11),'(I5)') IPID
        ENDIF
        WRITE(MYNAM,'("PE", A11)') PNUMBER
C
      ENDIF
C
      IF(IPID.EQ.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'MAITRE PROCESSEUR NUMERO ',
     &                  IPID,' OF THE GROUP OF ',NCSIZE
        IF(LNG.EQ.2) WRITE(LU,*) 'MASTER PROCESSOR NUMBER ',
     &                  IPID,' OF THE GROUP OF ',NCSIZE
      ELSE
        OPEN(UNIT=LU,FILE=MYNAM//'.LOG', FORM='FORMATTED',
     &         STATUS='UNKNOWN')
        IF(LNG.EQ.1) WRITE(LU,*) 'ESCLAVE PROCESSEUR NUMERO ',
     &                           IPID,' IN THE GROUP OF ',NCSIZE
        IF(LNG.EQ.2) WRITE(LU,*) 'SLAVE  PROCESSOR NUMBER ',
     &                           IPID,' IN THE GROUP OF ',NCSIZE
      ENDIF
C
C LECTURE DU NOMBRE DE PROCESSEURS ET NOM DE L'EXECUTABLE
C
      NCAR=0
      NPREAD=1
      YAPARA=.FALSE.
      INQUIRE(FILE='./PARAL',EXIST=YAPARA)

      IF(YAPARA) THEN
        OPEN(40,FILE='PARAL',FORM='FORMATTED',ACTION='READ')
        READ(40,*) NPREAD
        IF(NPREAD.NE.NCSIZE) THEN
          WRITE (LU,*)
     &      'P_INIT: FILE PARAL IS INCONSISTENT WITH MPI PARAMETERS'
          WRITE (LU,*) 'MPI NCSIZE   = ',NCSIZE
          WRITE (LU,*) 'PARAL NCSIZE = ',NPREAD
        ENDIF
CC        IF(LNG.EQ.1) WRITE(LU,*)'CALCUL AVEC ',NPREAD,' PROCESSEURS'
CC        IF(LNG.EQ.2) WRITE(LU,*)'COMPUTATION WITH ',NPREAD,' PROCESSORS'
        READ(40,*) NCAR
        READ(40,100) CHAINE
100     FORMAT(A250)
        EXE(1:NCAR+5)=CHAINE(1:NCAR) // 'A.EXE'
        IF(LNG.EQ.1) WRITE(LU,*) 'FICHIER EXECUTABLE : ',EXE(1:NCAR+5)
        IF(LNG.EQ.2) WRITE(LU,*) 'EXECUTABLE FILE: ',EXE(1:NCAR+5)
        CLOSE(40)
      ELSEIF(IPID.EQ.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'P_INIT: FICHIER PARAL NON TROUVE'
        IF(LNG.EQ.2) WRITE(LU,*) 'P_INIT: FILE PARAL NOT FOUND'
        STOP
      ENDIF
C
C THE BARRIER COMES USUALLY UNEXPECTED.
C
      CALL MPI_BARRIER(MPI_COMM_WORLD,IER)
      IF (IER.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'P_INIT: ERREUR DANS MPI_BARRIER'
        IF(LNG.EQ.2) WRITE(LU,*) 'P_INIT: ERROR MPI_BARRIER'
        WRITE(LU,*) 'MPI ERROR ',IER
        STOP
      ELSE
        IF(LNG.EQ.1) WRITE(LU,*) 'BARRIERE PASSEE'
        IF(LNG.EQ.2) WRITE(LU,*) 'BARRIER PASSED'
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END









