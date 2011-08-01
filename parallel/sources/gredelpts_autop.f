C
C 04/06/2008
C
C Warning by JMH: there is a CALL EXIT(ICODE) which is a Fortran extension
C                 it will not work with some compilers, like Nag
C
C
C
C
C
cjaj 2001/2
c     slightly changed to deal with: 
c     (1) arbitrary number of subdomains
c     (2) arbitrary names of the geometry and result files
c     (3) automatic parallel runs
c
CHW   IMPROVED READING OF DATASETS, 20.02.2003, BAW-HAMBURG
C
cjaj  added exit codes Fri Mar 14 15:47:51 MET 2003
C
      PROGRAM GREDELPTS
C
C     MERGES THE RESULTS OF A PARALLEL COMPUTATION WITH COUPLING WITH DELWAQ
C     TO WRITE A SINGLE FILE IN DELWAQ FORMAT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      INTEGER LI
      COMMON/INFO/LNG,LU
C
      CHARACTER(LEN=30) GEO   
C
      INTEGER IPID,ERR,FU
      INTEGER NELEM,ECKEN,NDUM,I,J,NBV1,NBV2,PARAM(10)
      INTEGER NPLAN,NPOIN2,NPOIN2LOC,NPLANLOC
      INTEGER NPROC,NRESU,NPOINMAX
      integer i_s, i_sp, i_len
      INTEGER IT
C
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: NPOIN,VERIF
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KNOLG
C      
C
      REAL   , DIMENSION(:)  , ALLOCATABLE :: GLOBAL_VALUE
      REAL   , DIMENSION(:)  , ALLOCATABLE :: LOCAL_VALUE      
C
      LOGICAL IS,ENDE
C
      CHARACTER*30 RES
      CHARACTER*50 RESPAR
      CHARACTER*11 EXTENS
      EXTERNAL    EXTENS
      INTRINSIC MAXVAL           

      
C-------------------------------------------------------------------------
C
      LI=5
      LU=6
      LNG=2
chw
cjaj introduce yourself with the version date
c
      write(lu,*) 'I am Gredelpts, cousin of Gretel from BAW Hamburg' 
      write(lu,*)
c
c reading file names and the number of processors / partitions
c
      write (lu, advance='no', 
     &    fmt='(/,'' Global geometry file: '')')
!      REWIND(LI)      
      read(li,*) geo
      write(lu,*) geo
c
      write (lu, advance='no', fmt='(/,'' Result file: '')')
      read(li,*) res   
      write(lu,*) res
c
      write (lu,advance='no',fmt='(/,'' Number of processors: '')')
      read (li,*) nproc
      write(lu,*) nproc
c      
      inquire (file=geo,exist=is)
      if (.not.is) then 
        write (lu,*) 'file does not exist: ', geo
        call plante (-1)
        stop
      end if     
c
      i_s  = len (res)
      i_sp = i_s + 1
      do i=1,i_s
         if(RES(i_sp-i:i_sp-i) .ne. ' ') exit
      enddo
      i_len=i_sp - i

C
C     FICHIER DE GEOMETRIE DU CALCUL, LU JUSQU'AUX 10 PARAMETRES:
C
      OPEN(2,FILE=GEO,FORM='UNFORMATTED',STATUS='OLD',ERR=990)  
      READ(2,ERR=990)
      READ(2,ERR=990) NBV1,NBV2
      DO 10 I=1,NBV1+NBV2
        READ(2,ERR=990)
10    CONTINUE
      GO TO 992
990   WRITE(LU,*) 'ERROR WHEN OPENING OR READING FILE: ',GEO
      call plante(-1)
      STOP
992   CONTINUE
C     LECTURE DES 10 PARAMETRES ET DE LA DATE
      READ(2) (PARAM(I),I=1,10)
      IF(PARAM(10).EQ.1) READ(2) (PARAM(I),I=1,6)
C
C     FICHIER  DE RESULTATS :
C
      OPEN(3,FILE=RES,FORM='UNFORMATTED',ERR=991)    
      GO TO 993
991   WRITE(LU,*) 'ERROR WHEN OPENING FILE: ',RES
      call plante(-1)
      STOP
993   CONTINUE
C
C     1) LECTURE DU DEBUT DU PREMIER FICHIER DE RESULTATS.
C
ccc      RESPAR=RES // EXTENS(2**IDIMS-1,0)
c
      respar=res(1:i_len) // extens(nproc-1,0)
c
      inquire (file=respar,exist=is)
      if (.not.is) then 
        write (lu,*) 'file does not exist: ', respar
        write (lu,*) 'check the number of processors'
        write (lu,*) 'and the result file core name'
        call plante(-1)
        stop
      end if  
c
      OPEN(4,FILE=RESPAR,FORM='UNFORMATTED',ERR=994)
      GO TO 995
994   WRITE(LU,*) 'ERROR WHEN OPENING FILE: ',RESPAR
      call plante(-1)
      STOP
995   CONTINUE
C
      READ(4) NPOIN2
      READ(4) NPLAN
      IF(NPLAN.EQ.1) NPLAN = 0
C
      CLOSE(4)
C
C  5 : 4 parametres
C      
      READ(2) NELEM,NPOIN2,ECKEN,NDUM
      WRITE(LU,*) '4 PARAMETERS IN GEOMETRY FILE'
      WRITE(LU,*) 'NELEM=',  NELEM
      WRITE(LU,*) 'NPOIN2=', NPOIN2
      WRITE(LU,*) 'ECKEN=',  ECKEN
      WRITE(LU,*) 'NDUM=',   NDUM
C
C  ALLOCATIONS DYNAMIQUES DES TABLEAUX
C
      ALLOCATE(NPOIN(NPROC),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'npoin')
      IF(NPLAN.EQ.0) THEN
        ALLOCATE(VERIF(NPOIN2)    ,STAT=ERR)
      ELSE
        ALLOCATE(VERIF(NPOIN2*NPLAN)    ,STAT=ERR)
      ENDIF
      IF(ERR.NE.0) call ALLOER (lu, 'verif')
C  GLOBAL_VALUES, STORING THE WHOLE DATASET (NBV1-VALUES)
      IF(NPLAN.EQ.0) THEN
        ALLOCATE(GLOBAL_VALUE(NPOIN2)       ,STAT=ERR)
      ELSE
        ALLOCATE(GLOBAL_VALUE(NPOIN2*NPLAN) ,STAT=ERR) 
      ENDIF
      IF(ERR.NE.0) call ALLOER (lu, 'global_value')
C
C  FIN ALLOCATION ...
C
C------------------------------------------------------------------------------
C
C OPENING FILES AND READING/SKIPPING HEADERS -> NPOIN(NPROC), NPOINMAX
C
      DO IPID = 0,NPROC-1
         FU = IPID +10
         RESPAR=RES(1:I_LEN) // EXTENS(NPROC-1,IPID)
         OPEN (FU,FILE=RESPAR,FORM='UNFORMATTED',ERR=998)
         GO TO 999
998      WRITE(LU,*) 'ERROR WHEN OPENING FILE: ',RESPAR,
     &                      ' USING FILE UNIT: ', FU
         call plante(-1)
         STOP
999      REWIND(FU)
         READ(FU) NPOIN(IPID+1)
         READ(FU) NPLANLOC
      END DO
C
      NPOINMAX = MAXVAL(NPOIN)
C TABLE FOR LOCAL-GLOBAL NUMBERS, 2D-FIELD
      IF(NPLAN.EQ.0) THEN
         ALLOCATE (KNOLG(NPOINMAX,NPROC),STAT=ERR)
      ELSE
         ALLOCATE (KNOLG(NPOINMAX/NPLAN,NPROC),STAT=ERR)
      ENDIF
      IF(ERR.NE.0) call ALLOER (lu, 'knolg')
C  LOCAL_VALUES, STORING THE WHOLE DATASET (NBV1-VALUES)
      ALLOCATE(LOCAL_VALUE(NPOINMAX),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'local_value')
C
C READING KNOLG(NPOIN,NPROC)
C
      DO IPID = 0,NPROC-1
         FU = IPID +10
         IF(NPLAN.EQ.0) THEN
            READ(FU) (KNOLG(I,IPID+1),I=1,NPOIN(IPID+1))
         ELSE
            READ(FU) (KNOLG(I,IPID+1),I=1,NPOIN(IPID+1)/NPLAN)
         ENDIF
      END DO
C
C READING DATASETS
C
      NRESU = 0
C
2000  NRESU = NRESU + 1
C
      IF(NPLAN.EQ.0) THEN
        DO I=1,NPOIN2
          VERIF(I)=0
        ENDDO
      ELSE
        DO I=1,NPOIN2*NPLAN
          VERIF(I)=0
        ENDDO
      ENDIF
C
      WRITE(LU,*)'TRY TO READ DATASET NO.',NRESU
C
      DO IPID = 0,NPROC-1
         FU = IPID +10
         CALL READ_DATASET
     +   (LOCAL_VALUE,NPOINMAX,NPOIN(IPID+1),IT,FU,ENDE)
         IF (ENDE) GOTO 3000
C STORE EACH DATASET
         IF(NPLAN.EQ.0) THEN
            DO I=1,NPOIN(IPID+1)
              GLOBAL_VALUE(KNOLG(I,IPID+1)) = LOCAL_VALUE(I)
              VERIF(KNOLG(I,IPID+1))   = 1
            END DO
         ELSE
            NPOIN2LOC = NPOIN(IPID+1)/NPLAN
            DO I=1,NPOIN2LOC
            DO J=1,NPLAN
            GLOBAL_VALUE(KNOLG(I,IPID+1) + NPOIN2   *(J-1)) = 
     +       LOCAL_VALUE(      I         + NPOIN2LOC*(J-1))
                   VERIF(KNOLG(I,IPID+1) + NPOIN2   *(J-1)) = 1
            END DO
            END DO
         ENDIF       
      END DO
C WRITING GLOBAL DATASET
      WRITE(LU,*)'WRITING DATASET NO.',NRESU,' TIME =',IT
C
      IF(NPLAN.EQ.0) THEN
         WRITE(3) IT, (GLOBAL_VALUE(I),I=1,NPOIN2)
      ELSE
         WRITE(3) IT, (GLOBAL_VALUE(I),I=1,NPOIN2*NPLAN)
      ENDIF
C VERIFICATIONS ...
      IF(NPLAN.EQ.0) THEN
        DO I=1,NPOIN2
          IF(VERIF(I).EQ.0) THEN
            WRITE(LU,*) 'ERROR, POINT I=',I,' FALSE FOR NRESU=',NRESU
          ENDIF
        END DO
      ELSE
        DO I=1,NPOIN2*NPLAN
          IF(VERIF(I).EQ.0) THEN
            WRITE(LU,*) 'ERROR, POINT I=',I,' FALSE FOR NRESU=',NRESU
          ENDIF
        END DO
      ENDIF     
C
      GO TO 2000
C
3000  WRITE(LU,*) 'END OF PROGRAM, ',NRESU-1,' DATASETS FOUND'
C
      CLOSE(2)       
      CLOSE(3)
C
      DO IPID = 0,NPROC-1
         FU = IPID +10
         CLOSE (FU)
      END DO      
C   
      STOP

      END PROGRAM GREDELPTS


C                       ***************************
                        CHARACTER*11 FUNCTION EXTENS
C                       ***************************
     *(N,IPID)
C
C***********************************************************************
C  PARA       VERSION 4.0         08/01/97        J-M HERVOUET (LNH)
C
C***********************************************************************
C
C      FONCTIONS: EXTENSION DES FICHIERS SUR CHAQUE PROCESSEUR.
C      ==========
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________|
C |     N          | -->| NOMBRE DE PROCESSEURS MOINS UN = NCSIZE-1
C |     IPID       | -->| NUMERO DU PROCESSEUR
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR :
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER IPID,N
C
C-----------------------------------------------------------------------
C
      IF(N.GT.0) THEN
C
        EXTENS='00000-00000'
C
        IF(N.LT.10) THEN
          WRITE(EXTENS(05:05),'(I1)') N
        ELSEIF(N.LT.100) THEN
          WRITE(EXTENS(04:05),'(I2)') N
        ELSEIF(N.LT.1000) THEN
          WRITE(EXTENS(03:05),'(I3)') N
        ELSEIF(N.LT.10000) THEN
          WRITE(EXTENS(02:05),'(I4)') N
        ELSE
          WRITE(EXTENS(01:05),'(I5)') N
        ENDIF
C
        IF(IPID.LT.10) THEN
          WRITE(EXTENS(11:11),'(I1)') IPID
        ELSEIF(IPID.LT.100) THEN
          WRITE(EXTENS(10:11),'(I2)') IPID
        ELSEIF(IPID.LT.1000) THEN
          WRITE(EXTENS(09:11),'(I3)') IPID
        ELSEIF(IPID.LT.10000) THEN
          WRITE(EXTENS(08:11),'(I4)') IPID
        ELSE
          WRITE(EXTENS(07:11),'(I5)') IPID
        ENDIF
C
      ELSE
C
        EXTENS='       '
C
      ENDIF  
C
C-----------------------------------------------------------------------
C
      RETURN
      END
C                         ***********************
                          SUBROUTINE READ_DATASET
C                         ***********************
     *(LOCAL_VALUE,NPOINMAX,NPOIN,IT,FU,ENDE)
C
      IMPLICIT NONE
C
      INTEGER NPOINMAX,NPOIN,FU
      INTEGER IPOIN
      INTEGER IT
C
      REAL LOCAL_VALUE(NPOINMAX)
C
      LOGICAL ENDE
C
      ENDE = .TRUE.
C
      READ(FU,END=999) IT, (LOCAL_VALUE(IPOIN),IPOIN=1,NPOIN)
C
      ENDE = .FALSE.
C
 999  RETURN
      END
 

      subroutine ALLOER (n, chfile)
      implicit none
      integer, intent(in) :: n
      character*(*), intent(in) :: chfile
      write(n,*) 'error by allocation of ',chfile
      call plante(-1)
      stop
      end subroutine ALLOER


      subroutine PLANTE(ival)
      implicit none
      integer, intent(in) :: ival
      integer icode      
      if (ival < 0) then      ! this indicates a controlled error
        icode = 1 
      else if (ival==0) then  ! this indicates a program failure
        icode = -1
      else                    ! this indicates a normal stop
        icode = 0
      endif 
      CALL EXIT(ICODE)
      stop    ! which is usually equivalent to call EXIT(0)
      end subroutine PLANTE
