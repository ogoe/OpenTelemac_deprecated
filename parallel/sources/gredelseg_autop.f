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
      PROGRAM GREDELSEG
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
      INTEGER NELEM,ECKEN,NDUM,I,J,K,NBV1,NBV2,PARAM(10)
      INTEGER NPLAN,NPOIN2,NPOIN2LOC,NOQ2,NPLANLOC,NSEG2LOC,NOQ2LOC
      INTEGER MBNDLOC,NPTFRLOC
      INTEGER NPROC,NRESU,NPOINMAX,NSEGMAX,NOQMAX,NPTFRMAX
      integer i_s, i_sp, i_len
      INTEGER IT
      INTEGER IDUM, NPTFR
      INTEGER IELM,NELEM2,NELMAX2,NPTFR2,NSEG2,KLOG,MBND2
      INTEGER MAXNVOIS,ISEG,IG1,IG2,IGTEMP,IVOIS,IL1,IL2
C
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: NPOIN,VERIF,NOQ,NSEG
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: MBND,NODENRS,NPTFRL
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KNOLG,KSEGLG
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NODENRSLOC,NBORLOC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: LIHBORLOC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IKLESA
C      
C
      REAL   , DIMENSION(:)  , ALLOCATABLE :: GLOBAL_VALUE
      REAL   , DIMENSION(:)  , ALLOCATABLE :: LOCAL_VALUE      
C
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IKLE       ! IKLE(SIZIKL,*) OU IKLE(NELMAX,*)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IFABOR     ! IFABOR(NELMAX,*) OU IFABOR(NELMAX2,*)
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: NVOIS,IADR ! NVOIS(NPOIN),IADR(NPOIN)
C
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: NELBOR,LIHBOR      ! NELBOR(NPTFR),LIHBOR(NPTFR)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NULONE             ! NULONE(NPTFR,2) OU NULONE(NPTFR)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KP1BOR             ! KP1BOR(NPTFR,2) OU KP1BOR(NPTFR)
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: NBOR               ! NBOR(*)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IKLBOR             ! IKLBOR(NPTFR,2)
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: T3                 ! T3(NPOIN) 
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: NBOR0,LIHBOR0      ! NBOR0(NPTFR),LIHBOR0(NPTFR)
C
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: GLOSEG         ! GLOSEG(MAXSEG,2)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ELTSEG,ORISEG  ! ELTSEG(NELMAX,*),ORISEG(NELMAX,3)     
C
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: GLOSEGLOC
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: SEGMENT
C
      REAL RDUM
C
      LOGICAL IS,ENDE
C
      CHARACTER*30 RES
      CHARACTER*50 RESPAR
      CHARACTER*11 EXTENS
      CHARACTER*30 CONLIM
      CHARACTER*7  FILETYPE
      EXTERNAL    EXTENS
      INTRINSIC MAXVAL           

      
      LI=5
      LU=6
      LNG=2
chw
cjaj introduce yourself with the version date
c
      write(lu,*) 'I am Gredelseg, cousin of Gretel from BAW Hamburg' 
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
      READ(4) FILETYPE
      READ(4) NPOIN2
      READ(4) NSEG2LOC
      READ(4) MBNDLOC
      READ(4) NOQ2LOC
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
      ALLOCATE(NOQ(NPROC),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'noq')
      ALLOCATE(NSEG(NPROC),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'nseg')
      ALLOCATE(MBND(NPROC),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'mbnd')
      ALLOCATE(IKLESA(3,NELEM),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'iklesa')
      ALLOCATE(NODENRS(NPOIN2),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'nodenrs')
      ALLOCATE(NPTFRL(NPROC),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'nptfr2loc')
C
      ALLOCATE(IFABOR(NELEM,3),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'ifabor')
      ALLOCATE(IKLE(NELEM,3),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'ikle')
      ALLOCATE(IADR(NPOIN2),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'iadr')
      ALLOCATE(NVOIS(NPOIN2),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'nvois')
      ALLOCATE(T3(NPOIN2),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 't3')
C
C  FIN ALLOCATION ...
C
C  6 : IKLE 
C 
      READ(2)  ((IKLESA(I,J),I=1,ECKEN),J=1,NELEM)
C
C----------------------------------------------------------------------
C
c
      IF(NPLAN.LE.1) THEN
        conlim = "T2DCLI"
      ELSE
        conlim = "T3DCLI"
      ENDIF
C
      OPEN(4,FILE=CONLIM,FORM='FORMATTED',ERR=996)
      GO TO 997
 996  WRITE(LU,*) 'ERROR WHEN OPENING FILE: ',CONLIM
      call plante(-1)
      STOP
 997  CONTINUE
C
      ALLOCATE(LIHBOR0(NPOIN2),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'lihbor')
      ALLOCATE(NBOR0(NPOIN2),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'nbor')
      DO I=1,NPOIN2
        READ(4,*,END=989) LIHBOR0(I),IDUM,IDUM,RDUM,RDUM,RDUM,RDUM,
     1                    IDUM,RDUM,RDUM,RDUM,NBOR0(I),IDUM
      ENDDO
C
      CLOSE(4) 
 989  NPTFR=I-1
C
      ALLOCATE(LIHBOR(NPTFR),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'lihbor')
      ALLOCATE(NBOR(NPTFR),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'nbor')
      ALLOCATE(NELBOR(NPTFR),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'nelbor')
      ALLOCATE(NULONE(NPTFR,2),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'nulone')
      ALLOCATE(KP1BOR(NPTFR,2),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'kp1bor')
      ALLOCATE(IKLBOR(NPTFR,2),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'iklbor')
      ALLOCATE(ELTSEG(NELEM,3),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'eltseg')
      ALLOCATE(ORISEG(NELEM,3),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'oriseg')
C
      MBND2=0
C
      DO I=1,NPOIN2
        NODENRS(I) = I
      ENDDO
C
      DO I=1,NPTFR
        NBOR(I)   = NBOR0(I)
        LIHBOR(I) = LIHBOR0(I)
        IF (LIHBOR(I).NE.2) THEN
          MBND2 = MBND2 + 1
          NODENRS(NBOR(I)) = -MBND2
        ENDIF
      ENDDO
C
C------------------------------------------------------------------------------
C
C LOCAL CONSTRUCTION OF GLOSEG
C 
C------------------------------------------------------------------------------
C
C     WITH PRISMS, DIFFERENT FROM 2D VALUES, OTHERWISE
C
      IELM = 11 ! WARNING EN DUR !!!
        NELEM2  =NELEM
        NELMAX2 =NELEM
        NPTFR2  =NPTFR
C
C     CALCUL DES VOISINS DES FACES DE BORD POUR LE MAILLAGE DE TRIANGLES
C
        DO J=1,NELEM
          DO I=1,3
            IKLE(J,I)=IKLESA(I,J)
          ENDDO
        ENDDO

      IF(IELM.EQ.11.OR.IELM.EQ.41.OR.IELM.EQ.51) THEN
        CALL VOISIN(IFABOR,NELEM2,NELEM,IELM,IKLE,
     *              NELEM,
     *              NPOIN2,IADR,NVOIS)
        MAXNVOIS = MAXVAL(NVOIS)/2
      ELSE
        WRITE(LU,*) 'UNEXPECTED ELEMENT IN INBIEF:',IELM
        CALL PLANTE(1)
        STOP
      ENDIF

      KLOG = 2 ! CONVENTION CONDITION A LA LIMITE DE PAROI, EN DUR !!!

      IF(IELM.EQ.11.OR.IELM.EQ.41.OR.IELM.EQ.51) THEN
        CALL ELEBD(NELBOR,NULONE,KP1BOR,
     *             IFABOR,NBOR,IKLE,NELEM,
     *             IKLBOR,NELEM2,NELMAX2,
     *             NPOIN2,NPTFR2,IELM,
     *             LIHBOR,KLOG,
     *             IADR,NVOIS,T3)
      ELSE 
        WRITE(LU,*) 'UNEXPECTED ELEMENT IN INBIEF:',IELM
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C  DATA STRUCTURE FOR EDGE-BASED STORAGE (FROM 5.9 ON ALWAYS DONE IN 2D)
C  SEE CALL TO COMP_SEG BELOW FOR COMPLETING THE STRUCTURE
C
      IF(IELM.EQ.11) THEN
C
         NSEG2 = (3*NELEM+NPTFR)/2
         NOQ2=NPLAN*(NSEG2+MBND2)+(NPLAN-1)*NPOIN2
         IF(NPLAN.EQ.0) THEN
           ALLOCATE(VERIF(NSEG2+MBND2),STAT=ERR)
         ELSE
           ALLOCATE(VERIF(NOQ2) ,STAT=ERR)
         ENDIF
         IF(ERR.NE.0) call ALLOER (lu, 'verifseg')
C
C  GLOBAL_VALUES, STORING THE WHOLE DATASET (NBV1-VALUES)
         IF(NPLAN.EQ.0) THEN
           ALLOCATE(GLOBAL_VALUE(NSEG2+MBND2),STAT=ERR)
         ELSE
           ALLOCATE(GLOBAL_VALUE(NOQ2),STAT=ERR) 
         ENDIF
         IF(ERR.NE.0) call ALLOER (lu, 'global_value')
C
         ALLOCATE(GLOSEG(NSEG2,2),STAT=ERR)
         IF(ERR.NE.0) call ALLOER (lu, 'gloseg')
C
      CALL STOSEG(IFABOR,NELEM,NELMAX2,NELMAX2,IELM,IKLE,
     *            NBOR,NPTFR,
     *            GLOSEG,NSEG2,    ! GLOSEG%MAXDIM1,
     *            ELTSEG,ORISEG,NSEG2,
     *            KP1BOR,NELBOR,NULONE)
      ENDIF
C
      ALLOCATE(SEGMENT(NPOIN2,MAXNVOIS,2),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'segment')
C
C INITIALISATION OF SEGMENT
      DO K=1,2
        DO J=1,MAXNVOIS
          DO I=1,NPOIN2
            SEGMENT(I,J,K) = 0
          ENDDO
        ENDDO
      ENDDO
C
      DO ISEG=1,NSEG2
        IG1 = GLOSEG(ISEG,1)
        IG2 = GLOSEG(ISEG,2)
C GLOBAL NUMBER IN INCREASING ORDER
        IF(IG1.GT.IG2) THEN
          IGTEMP = IG1
          IG1 = IG2
          IG2 = IGTEMP
        ENDIF
        IVOIS=1
        DO WHILE ((SEGMENT(IG1,IVOIS,1).NE.0).AND.(IVOIS.LE.MAXNVOIS))
          IVOIS = IVOIS + 1
        ENDDO
        SEGMENT(IG1,IVOIS,1) = IG2
        SEGMENT(IG1,IVOIS,2) = ISEG
      ENDDO
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
         READ(FU) FILETYPE
         READ(FU) NPOIN(IPID+1)
         READ(FU) NSEG(IPID+1)
         READ(FU) MBND(IPID+1)
         READ(FU) NOQ(IPID+1)
         READ(FU) NPLANLOC
         READ(FU) NPTFRL(IPID+1)
      END DO
C
      NPOINMAX = MAXVAL(NPOIN)
      NSEGMAX = MAXVAL(NSEG)
      NOQMAX = MAXVAL(NOQ)
      NPTFRMAX = MAXVAL(NPTFRL)
C TABLE FOR LOCAL-GLOBAL NUMBERS, 2D-FIELD
      ALLOCATE (GLOSEGLOC(NSEGMAX,2,NPROC),STAT=ERR)
      IF(NPLAN.EQ.0) THEN
         ALLOCATE (KNOLG(NPOINMAX,NPROC),STAT=ERR)
         ALLOCATE (KSEGLG(NSEGMAX,NPROC),STAT=ERR)
         ALLOCATE (NODENRSLOC(NPOINMAX,NPROC),STAT=ERR)
         ALLOCATE (NBORLOC(NPTFRMAX,NPROC),STAT=ERR)
         ALLOCATE (LIHBORLOC(NPTFRMAX,NPROC),STAT=ERR)
      ELSE
         ALLOCATE (KNOLG(NPOINMAX/NPLAN,NPROC),STAT=ERR)
         ALLOCATE (KSEGLG(NOQMAX,NPROC),STAT=ERR)
         ALLOCATE (NODENRSLOC(NPOINMAX/NPLAN,NPROC),STAT=ERR)
         ALLOCATE (NBORLOC(NPTFRMAX,NPROC),STAT=ERR)
         ALLOCATE (LIHBORLOC(NPTFRMAX,NPROC),STAT=ERR)
      ENDIF
      IF(ERR.NE.0) call ALLOER (lu, 'knolg')
      IF(ERR.NE.0) call ALLOER (lu, 'kseglg')
      IF(ERR.NE.0) call ALLOER (lu, 'nodenrsloc')
      IF(ERR.NE.0) call ALLOER (lu, 'nborloc')
C  LOCAL_VALUES, STORING THE WHOLE DATASET (NBV1-VALUES)
      ALLOCATE(LOCAL_VALUE(NOQMAX),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'local_value')
C
C READING KNOLG(NPOIN,NPROC)
C
C
      IF(NPLAN.EQ.0) THEN
        DO I=1,NSEG2+MBND2
          VERIF(I)=0
        ENDDO
      ELSE
        DO I=1,NOQ2
          VERIF(I)=0
        ENDDO
      ENDIF
C
      DO IPID = 0,NPROC-1
         FU = IPID +10
C CHECK
         IF(NPLAN.EQ.0) THEN
            READ(FU) (KNOLG(I,IPID+1),I=1,NPOIN(IPID+1))
            READ(FU) ((GLOSEGLOC(I,J,IPID+1),J=1,2),I=1,NSEG(IPID+1))
            READ(FU) (NODENRSLOC(I,IPID+1),I=1,NPOIN(IPID+1))
            READ(FU) (NBORLOC(I,IPID+1),I=1,NPTFRL(IPID+1))
            READ(FU) (LIHBORLOC(I,IPID+1),I=1,NPTFRL(IPID+1))
         ELSE
            READ(FU) (KNOLG(I,IPID+1),I=1,NPOIN(IPID+1)/NPLAN)
            READ(FU) ((GLOSEGLOC(I,J,IPID+1),J=1,2),I=1,NSEG(IPID+1))
            READ(FU) (NODENRSLOC(I,IPID+1),I=1,NPOIN(IPID+1)/NPLAN)
            READ(FU) (NBORLOC(I,IPID+1),I=1,NPTFRL(IPID+1))
            READ(FU) (LIHBORLOC(I,IPID+1),I=1,NPTFRL(IPID+1))
         ENDIF
C
C INITIALISATION OF SEGMENT
C
         DO ISEG=1,NSEG(IPID+1)
           IL1 = GLOSEGLOC(ISEG,1,IPID+1)
           IL2 = GLOSEGLOC(ISEG,2,IPID+1)
           IG1 = KNOLG(IL1,IPID+1)
           IG2 = KNOLG(IL2,IPID+1)
C GLOBAL NUMBER IN INCREASING ORDER
           IF(IG1.GT.IG2) THEN
             IGTEMP = IG1
             IG1 = IG2
             IG2 = IGTEMP
           ENDIF
           IVOIS=1
           DO WHILE ((SEGMENT(IG1,IVOIS,1).NE.IG2)
     1               .AND.(IVOIS.LE.MAXNVOIS))
             IVOIS = IVOIS + 1
           ENDDO
           IF(IVOIS.LE.MAXNVOIS) THEN
             KSEGLG(ISEG,IPID+1) = SEGMENT(IG1,IVOIS,2)
           ENDIF
         ENDDO
C
      END DO
C
C FURTHER VERIFICATIONS
C
C READING DATASETS
C
      NRESU = 0
C
2000  NRESU = NRESU + 1
C
      IF(NPLAN.EQ.0) THEN
        DO I=1,NSEG2+MBND2
          VERIF(I)=0
        ENDDO
      ELSE
        DO I=1,NOQ2
          VERIF(I)=0
        ENDDO
      ENDIF
C
      WRITE(LU,*)'TRY TO READ DATASET NO.',NRESU
C
      IF(NPLAN.EQ.0) THEN
        DO I=1,NSEG2+MBND2
          GLOBAL_VALUE(I) = 0.D0
        ENDDO
      ELSE
        DO I=1,NOQ2
          GLOBAL_VALUE(I) = 0.D0
        ENDDO
      ENDIF
C
      DO IPID = 0,NPROC-1
         FU = IPID +10
C FOR READING LOCAL X INSTEAD OF READ_DATASET
         CALL READ_DATASET
     +   (LOCAL_VALUE,NOQMAX,NOQ(IPID+1),IT,FU,ENDE)
         IF (ENDE) GOTO 3000
C STORE EACH DATASET
         IF(NPLAN.EQ.0) THEN
            NSEG2LOC  = NSEG(IPID+1)
            NPTFRLOC  = NPTFRL(IPID+1)
            DO I=1,NSEG2LOC
              GLOBAL_VALUE(KSEGLG(I,IPID+1)) = 
     *        GLOBAL_VALUE(KSEGLG(I,IPID+1)) + LOCAL_VALUE(I)
              VERIF(KSEGLG(I,IPID+1)) =   VERIF(KSEGLG(I,IPID+1))
     *                                     + 1
            END DO
C
           DO I=1,NPTFRLOC
             IF(LIHBORLOC(I,IPID+1).NE.2) THEN
               IF(FILETYPE(1:7).EQ.'SUMAREA') THEN
                 GLOBAL_VALUE(-NODENRS(KNOLG(NBORLOC(I,IPID+1),IPID+1))
     *                        + NSEG2) =
     *           LOCAL_VALUE(-NODENRSLOC(NBORLOC(I,IPID+1),IPID+1)
     *                       + NSEG2LOC)
                 VERIF( -NODENRS(KNOLG(NBORLOC(I,IPID+1),IPID+1))
     *                 + NSEG2) = 1
               ELSEIF(FILETYPE(1:7).EQ.'SUMFLOW') THEN
                 GLOBAL_VALUE(-NODENRS(KNOLG(NBORLOC(I,IPID+1),IPID+1))
     *                        + NSEG2) =
     *           GLOBAL_VALUE(-NODENRS(KNOLG(NBORLOC(I,IPID+1),IPID+1))
     *                        + NSEG2) +
     *           LOCAL_VALUE(-NODENRSLOC(NBORLOC(I,IPID+1),IPID+1)
     *                       + NSEG2LOC)
                 VERIF( -NODENRS(KNOLG(NBORLOC(I,IPID+1),IPID+1))
     *                 + NSEG2) =
     *           VERIF( -NODENRS(KNOLG(NBORLOC(I,IPID+1),IPID+1))
     *                 + NSEG2) + 1
               ELSE
                 WRITE(LU,*) 'CAS NON PREVU'
                 STOP
               ENDIF
             ENDIF
           END DO
C
         ELSE
           NPOIN2LOC = NPOIN(IPID+1)/NPLAN
           NSEG2LOC  = NSEG(IPID+1)
           MBNDLOC   = MBND(IPID+1)
           NPTFRLOC  = NPTFRL(IPID+1)
           DO I=1,NSEG2LOC
             DO J=1,NPLAN
               GLOBAL_VALUE(KSEGLG(I,IPID+1) + (NSEG2+MBND2)*(J-1)) =
     *         GLOBAL_VALUE(KSEGLG(I,IPID+1) + (NSEG2+MBND2)*(J-1)) +
     *         LOCAL_VALUE(       I      + (NSEG2LOC+MBNDLOC)*(J-1))
               VERIF(KSEGLG(I,IPID+1) + (NSEG2+MBND2)*(J-1)) =
     *       + VERIF(KSEGLG(I,IPID+1) + (NSEG2+MBND2)*(J-1)) + 1
             END DO
           END DO
C
           DO I=1,NPTFRLOC
             IF(LIHBORLOC(I,IPID+1).NE.2) THEN
               DO J=1,NPLAN
                 IF(FILETYPE(1:7).EQ.'SUMAREA') THEN
                  GLOBAL_VALUE(-NODENRS(KNOLG(NBORLOC(I,IPID+1),IPID+1))
     *                         + NSEG2 + (NSEG2+MBND2)*(J-1)) =
     *             LOCAL_VALUE(-NODENRSLOC(NBORLOC(I,IPID+1),IPID+1)
     *                         + NSEG2LOC + (NSEG2LOC+MBNDLOC)*(J-1))
                   VERIF( -NODENRS(KNOLG(NBORLOC(I,IPID+1),IPID+1))
     *                   + NSEG2 + (NSEG2+MBND2)*(J-1)) = 1
                 ELSEIF(FILETYPE(1:7).EQ.'SUMFLOW') THEN
                  GLOBAL_VALUE(-NODENRS(KNOLG(NBORLOC(I,IPID+1),IPID+1))
     *                         + NSEG2 + (NSEG2+MBND2)*(J-1)) =
     *            GLOBAL_VALUE(-NODENRS(KNOLG(NBORLOC(I,IPID+1),IPID+1))
     *                         + NSEG2 + (NSEG2+MBND2)*(J-1)) +
     *             LOCAL_VALUE(-NODENRSLOC(NBORLOC(I,IPID+1),IPID+1)
     *                         + NSEG2LOC + (NSEG2LOC+MBNDLOC)*(J-1))
                   VERIF( -NODENRS(KNOLG(NBORLOC(I,IPID+1),IPID+1))
     *                   + NSEG2 + (NSEG2+MBND2)*(J-1)) =
     *             VERIF( -NODENRS(KNOLG(NBORLOC(I,IPID+1),IPID+1))
     *                   + NSEG2 + (NSEG2+MBND2)*(J-1)) + 1
                 ELSE
                   WRITE(LU,*) 'CAS NON PREVU'
                   STOP
                 ENDIF
               END DO
             ENDIF
           END DO
C
           DO I=1,NPOIN2LOC
             DO J=1,NPLAN-1
               IF(FILETYPE(1:7).EQ.'SUMAREA') THEN
                 GLOBAL_VALUE(  KNOLG(I,IPID+1) + NPOIN2*(J-1)
     *                        + (NSEG2+MBND2)*NPLAN) =
     *         LOCAL_VALUE(I+NPOIN2LOC*(J-1)+(NSEG2LOC+MBNDLOC)*NPLAN)
                 VERIF( KNOLG(I,IPID+1) + NPOIN2*(J-1)
     *               + (NSEG2+MBND2)*NPLAN) = 1
               ELSEIF(FILETYPE(1:7).EQ.'SUMFLOW') THEN
                 GLOBAL_VALUE( KNOLG(I,IPID+1) + NPOIN2*(J-1)
     *                        + (NSEG2+MBND2)*NPLAN) =
     *           GLOBAL_VALUE( KNOLG(I,IPID+1) + NPOIN2*(J-1)
     *                        + (NSEG2+MBND2)*NPLAN) +
     *         LOCAL_VALUE(I+NPOIN2LOC*(J-1)+(NSEG2LOC+MBNDLOC)*NPLAN)
                 VERIF( KNOLG(I,IPID+1) + NPOIN2*(J-1)
     *               + (NSEG2+MBND2)*NPLAN) =
     *           VERIF( KNOLG(I,IPID+1) + NPOIN2*(J-1)
     *               + (NSEG2+MBND2)*NPLAN) + 1
               ELSE
                 WRITE(LU,*) 'CAS NON PREVU'
                 STOP
               ENDIF
             END DO
           END DO
         ENDIF       
      END DO
C WRITING GLOBAL DATASET
      WRITE(LU,*)'WRITING DATASET NO.',NRESU,' TIME =',IT
C
      IF(NPLAN.EQ.0) THEN
         WRITE(3) IT, (GLOBAL_VALUE(I),I=1,NSEG2+MBND2)
      ELSE
         WRITE(3) IT, (GLOBAL_VALUE(I),I=1,NOQ2)
      ENDIF
C VERIFICATIONS ...
      IF(NPLAN.EQ.0) THEN
        DO I=1,NSEG2+MBND2
          IF(VERIF(I).EQ.0) THEN
            WRITE(LU,*) 'ERROR, SEGMENT I=',I,' FALSE FOR NRESU=',NRESU
          ENDIF
        END DO
      ELSE
        DO I=1,NOQ2
          IF(VERIF(I).EQ.0) THEN
            WRITE(LU,*) 'ERROR, SEGMENT I=',I,' FALSE FOR NRESU=',NRESU
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
        !!!Fabs
      
      STOP

      END PROGRAM GREDELSEG


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

C                       *****************
                        SUBROUTINE VOISIN
C                       *****************
C
     *(IFABOR,NELEM,NELMAX,IELM,IKLE,SIZIKL,
     * NPOIN,IADR,NVOIS)
C
C***********************************************************************
C BIEF VERSION 5.9        19/02/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C    FONCTION : CONSTRUCTION DU TABLEAU IFABOR, OU IFABOR(IELEM,IFACE)
C               EST LE NUMERO GLOBAL DU VOISIN DE LA FACE IFACE DE
C               L'ELEMENT IELEM SI CE VOISIN EXISTE ET 0 SI LA FACE EST
C               SUR LA FRONTIERE DU DOMAINE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    IFABOR      |<-- | TABLEAU DES VOISINS DES FACES.
C |    NELEM       | -->| NOMBRE D'ELEMENTS DANS LE MAILLAGE.
C |    NELMAX      | -->| NOMBRE MAXIMUM D'ELEMENTS DANS LE MAILLAGE.
C |                |    | (CAS DES MAILLAGES ADAPTATIFS)
C |    IELM        | -->| 11: TRIANGLES
C |                |    | 21: QUADRILATERES
C |    IKLE        | -->| NUMEROS GLOBAUX DES POINTS DE CHAQUE ELEMENT
C |    NPOIN       | -->| NOMBRE TOTAL DE POINTS DU DOMAINE
C |________________|____|_______________________________________________
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT DANS TELEMAC 2D : PREDAT
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER SIZIKL,NELEM,NELMAX,IELM,NPOIN
      INTEGER IKLE(SIZIKL,*)
      INTEGER IFABOR(NELMAX,*)
      INTEGER NVOIS(NPOIN),IADR(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NFACE,NDP,KEL,IMAX,IFACE,IELEM,M1,M2,IV,IELEM2,IFACE2
      INTEGER I,ERR,I1,I2,IDIMAT
C
      INTEGER SOMFAC(2,4,2)
      DATA SOMFAC / 1,2 , 2,3 , 3,1 , 0,0   ,  1,2 , 2,3 , 3,4 , 4,1 /
C
C     TABLEAUX DE TRAVAIL ALLOUES DYNAMIQUEMENT
C
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT1,MAT2,MAT3
C
C-----------------------------------------------------------------------
C
      IF(IELM.EQ.21) THEN
C       QUADRILATERES
        NFACE = 4
C       NOMBRE DE POINTS PAR ELEMENT
        NDP = 4
C       ADRESSE DANS SOMFAC
        KEL = 2
      ELSEIF(IELM.EQ.11.OR.IELM.EQ.41.OR.IELM.EQ.51) THEN
C       TRIANGLES
        NFACE = 3
C       NOMBRE DE POINTS PAR ELEMENT
        NDP = 3
C       ADRESSE DANS SOMFAC
        KEL = 1
      ELSE
        IF(LNG.EQ.1) WRITE(LU,98) IELM
        IF(LNG.EQ.2) WRITE(LU,99) IELM
98      FORMAT(1X,'VOISIN: IELM=',1I6,' TYPE D''ELEMENT NON PREVU')
99      FORMAT(1X,'VOISIN: IELM=',1I6,' TYPE OF ELEMENT NOT AVAILABLE')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C     IDIMAT EST UNE MAJORATION DE LA SOMME DES NOMBRES DE VOISINS DE
C     TOUS LES POINTS (VOISIN = RELIE PAR UN SEGMENT)
C
      IDIMAT = NDP*2*NELEM
C
      ALLOCATE(MAT1(IDIMAT),STAT=ERR)
      ALLOCATE(MAT2(IDIMAT),STAT=ERR)
      ALLOCATE(MAT3(IDIMAT),STAT=ERR)
C
      IF(ERR.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,1000) ERR
        IF(LNG.EQ.2) WRITE(LU,2000) ERR
1000    FORMAT(1X,'VOISIN : ERREUR A L''ALLOCATION DE MEMOIRE : ',/,1X,
     *            'CODE D''ERREUR : ',1I6)
2000    FORMAT(1X,'VOISIN: ERROR DURING ALLOCATION OF MEMORY: ',/,1X,
     *            'ERROR CODE: ',1I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C  CALCUL DU TABLEAU NVOIS POUR CHAQUE POINT
C  ATTENTION : NVOIS N'EST QU'UNE MAJORATION DU NOMBRE DE VOISINS
C              DONT LA SOMME VA FAIRE IDIMAT
C
      DO I=1,NPOIN
        NVOIS(I) = 0
      ENDDO
C
      DO IFACE = 1,NFACE
        DO IELEM=1,NELEM
          I1 = IKLE( IELEM , SOMFAC(1,IFACE,KEL) )
          I2 = IKLE( IELEM , SOMFAC(2,IFACE,KEL) )
          NVOIS(I1) = NVOIS(I1) + 1
          NVOIS(I2) = NVOIS(I2) + 1
        ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C
C  CALCUL DES ADRESSES DE CHAQUE POINT DANS UNE STRUCTURE DE TYPE
C  MATRICE COMPACTE
C
      IADR(1) = 1
      DO 50 I= 2,NPOIN
        IADR(I) = IADR(I-1) + NVOIS(I-1)
50    CONTINUE
C
      IMAX = IADR(NPOIN) + NVOIS(NPOIN) - 1
      IF(IMAX.GT.IDIMAT) THEN
        IF(LNG.EQ.1) WRITE(LU,51) IDIMAT,IMAX
        IF(LNG.EQ.2) WRITE(LU,52) IDIMAT,IMAX
51      FORMAT(1X,'VOISIN: TAILLE DE MAT1,2,3 (',1I9,') INSUFFISANTE',/,
     *         1X,'IL FAUT AU MOINS : ',1I9)
52      FORMAT(1X,'VOISIN: SIZE OF MAT1,2,3 (',1I9,') TOO SHORT',/,
     *         1X,'MINIMUM SIZE: ',1I9)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C  INITIALISATION A ZERO DE LA MATRICE COMPACTE
C
      DO I=1,IMAX
        MAT1(I) = 0
      ENDDO
C
C-----------------------------------------------------------------------
C
C  BOUCLE SUR LES FACES DE CHAQUE ELEMENT :
C
      DO 60 IFACE = 1 , NFACE
      DO 70 IELEM = 1 , NELEM
C
      IFABOR(IELEM,IFACE) = -1
C
C        NUMEROS GLOBAUX DES POINTS DE LA FACE :
C
         I1 = IKLE( IELEM , SOMFAC(1,IFACE,KEL) )
         I2 = IKLE( IELEM , SOMFAC(2,IFACE,KEL) )
C
C        NUMEROS GLOBAUX ORDONNES :
C
         M1 = MIN0(I1,I2)
         M2 = MAX0(I1,I2)
C
         DO 80 IV = 1,NVOIS(M1)
C
           IF(MAT1(IADR(M1)+IV-1).EQ.0) THEN
              MAT1(IADR(M1)+IV-1)=M2
              MAT2(IADR(M1)+IV-1)=IELEM
              MAT3(IADR(M1)+IV-1)=IFACE
              GO TO 81
           ELSEIF(MAT1(IADR(M1)+IV-1).EQ.M2) THEN
              IELEM2 = MAT2(IADR(M1)+IV-1)
              IFACE2 = MAT3(IADR(M1)+IV-1)
              IFABOR(IELEM,IFACE) = IELEM2
              IFABOR(IELEM2,IFACE2) = IELEM
              GO TO 81
           ENDIF
C
80       CONTINUE
C
         IF(LNG.EQ.1) WRITE(LU,82)
         IF(LNG.EQ.2) WRITE(LU,83)
82       FORMAT(1X,'VOISIN : ERREUR DANS LE MAILLAGE       ',/,1X,
     *             '         PEUT-ETRE DES POINTS CONFONDUS')
83       FORMAT(1X,'VOISIN : ERROR IN THE MESH             ',/,1X,
     *             '         MAYBE SUPERIMPOSED POINTS     ')
         CALL PLANTE(1)
         STOP
C
81       CONTINUE
C
70    CONTINUE
60    CONTINUE
C
C-----------------------------------------------------------------------
C
      DEALLOCATE(MAT1)
      DEALLOCATE(MAT2)
      DEALLOCATE(MAT3)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
C                       ****************
                        SUBROUTINE ELEBD
C                       ****************
C
     *(NELBOR,NULONE,KP1BOR,IFABOR,NBOR,IKLE,SIZIKL,IKLBOR,NELEM,NELMAX,
     * NPOIN,NPTFR,IELM,LIHBOR,KLOG,T1,T2,T3)
C
C***********************************************************************
C BIEF VERSION 5.9        20/03/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C COPYRIGHT 1999
C***********************************************************************
C
C    PRISMES DECOUPES EN TETRAEDRES
C
C    FONCTION : 1) CONSTRUCTION DES TABLEAUX NELBOR ET NULONE
C               2) CONSTRUCTION DU TABLEAU KP1BOR
C               3) DISTINCTION DANS LE TABLEAU IFABOR ENTRE
C                  LES FACES DE BORD SOLIDES ET LES FACES LIQUIDES
C               4) CALCUL DE IKLBOR. 
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    NELBOR      |<-- | NUMERO DE L'ELEMENT ADJACENT AU KIEME SEGMENT|
C |    NULONE      |<-- | NUMERO LOCAL D'UN POINT DE BORD DANS         |
C |                |    | L'ELEMENT ADJACENT DONNE PAR NELBOR          |
C |    KP1BOR      |<-- | NUMERO DU POINT SUIVANT LE POINT DE BORD K.  |
C |    IFABOR      | -->| TABLEAU DES VOISINS DES FACES.
C |    NBOR        | -->| NUMERO GLOBAL DU POINT DE BORD K.
C |    IKLE        | -->| NUMEROS GLOBAUX DES POINTS DE CHAQUE ELEMENT.
C |    NELEM       | -->| NOMBRE TOTAL D'ELEMENTS DANS LE MAILLAGE.
C |    T1,2,3      | -->| TABLEAUX DE TRAVAIL ENTIERS.
C |    NPOIN       | -->| NOMBRE TOTAL DE POINTS DU DOMAINE.
C |    NPTFR       | -->| NOMBRE DE POINTS FRONTIERES.
C |    IELM        | -->| TYPE D'ELEMENT.
C |                |    | 11 : TRIANGLES.
C |                |    | 21 : QUADRILATERES.
C |    LIHBOR      | -->| TYPES DE CONDITIONS AUX LIMITES SUR H
C |    KLOG        | -->| CONVENTION POUR LA CONDITION LIMITE DE PAROI
C |    MXPTVS      | -->| NOMBRE MAXIMUM DE VOISINS D'UN POINT
C |    MXELVS      | -->| NOMBRE MAXIMUM D'ELEMENTS AUTOUR D'UN POINT
C |________________|____|______________________________________________|
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : 
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER KLOG,NELMAX,NELEM,SIZIKL
      INTEGER NPOIN,NPTFR,IELM
      INTEGER NELBOR(NPTFR),LIHBOR(NPTFR)
      INTEGER NULONE(NPTFR,2),KP1BOR(NPTFR,2)
      INTEGER NBOR(NPTFR)
      INTEGER IFABOR(NELMAX,3)
      INTEGER IKLE(SIZIKL,3)
      INTEGER IKLBOR(NPTFR,2)
      INTEGER T1(NPOIN),T2(NPOIN),T3(NPOIN) 
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,NFACE,NPT,KEL,IPOIN
      INTEGER K,IFACE,I1,I2,N1,N2,IPT,IEL,I,K1,K2
C
      INTEGER SOMFAC(2,4,2)
C
      DATA SOMFAC / 1,2 , 2,3 , 3,1 , 0,0   ,  1,2 , 2,3 , 3,4 , 4,1 /
C
C-----------------------------------------------------------------------
C
      IF(IELM.EQ.11.OR.IELM.EQ.41.OR.IELM.EQ.51) THEN
C       TRIANGLES
        NFACE = 3
        NPT = 3
        KEL = 1
      ELSE
        IF(LNG.EQ.1) WRITE(LU,900) IELM
        IF(LNG.EQ.2) WRITE(LU,901) IELM
900     FORMAT(1X,'ELEBD : IELM=',1I6,' TYPE D''ELEMENT INCONNU')
901     FORMAT(1X,'ELEBD: IELM=',1I6,' UNKNOWN TYPE OF ELEMENT')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C  INITIALISATION DE T1,2,3 A ZERO
C
      DO IPOIN=1,NPOIN
        T1(IPOIN) = 0
        T2(IPOIN) = 0
        T3(IPOIN) = 0
      ENDDO
C
C  ON STOCKE K DANS TRAV(*,3) A L'ADRESSE NBOR(K)
C  CE QUI PERMET DE PASSER DE NUMERO GLOBAL A NUMERO DE BORD
C
      DO K = 1, NPTFR
         T3(NBOR(K)) = K
      ENDDO
C
C  BOUCLE SUR TOUTES LES FACES DE TOUS LES ELEMENTS :
C
      DO 20 IFACE = 1 , NFACE
      DO 10 IELEM = 1 , NELEM
C
      IF(IFABOR(IELEM,IFACE).EQ.-1) THEN
C
C      C'EST UNE VRAIE FACE DE BORD (LES FACES INTERNES EN PARALLELISME
C                                    SONT MARQUEES AVEC DES -2).
C      NUMEROS GLOBAUX DES POINTS DE LA FACE :
C
       I1 = IKLE( IELEM , SOMFAC(1,IFACE,KEL) )
       I2 = IKLE( IELEM , SOMFAC(2,IFACE,KEL) )
C
C      ON STOCKE DANS T1 ET T2 A L'ADRESSE I1 : I2 ET IELEM
C
       T1(I1) = I2
       T2(I1) = IELEM
C
C      UNE FACE LIQUIDE EST RECONNUE AVEC LA CONDITION LIMITE SUR H
C
C      07/02/03 IF(NPTFR...  COURTESY OLIVER GOETHEL, HANNOVER UNIVERSITY
       IF(NPTFR.GT.0) THEN
       IF(LIHBOR(T3(I1)).NE.KLOG.AND.LIHBOR(T3(I2)).NE.KLOG) THEN
C        FACE LIQUIDE : IFABOR=0  FACE SOLIDE : IFABOR=-1
         IFABOR(IELEM,IFACE)=0
       ENDIF
       ENDIF
C
      ENDIF
C
10    CONTINUE
20    CONTINUE
C
CPARA
C
C VALEURS NULLES MISES POUR NELBOR QUAND L'ELEMENT EST DANS UN AUTRE
C SOUS-DOMAINE.
C
C
C      LIGNES SUIVANTES DEJA FAITES PLUS HAUT
C      IF(NCSIZE.GT.1) THEN
C        DO 39 K1=1,NPTFR
C          NELBOR(K)=0
C39      CONTINUE
C      ENDIF
C
CPARAFIN
C
C  BOUCLE SUR TOUS LES POINTS:
C
C     07/02/03 IF(NPTFR...  CORRECTION BY OLIVER GOETHELS, HANNOVER
      IF(NPTFR.GT.0) THEN
      DO I = 1 , NPOIN
         IF(T1(I).NE.0) THEN
C          POINT SUIVANT
           KP1BOR(T3(I),1)=T3(T1(I))
C          POINT PRECEDENT
           KP1BOR(T3(T1(I)),2)=T3(I)
           NELBOR(T3(I))=T2(I)
         ENDIF
      ENDDO
      ENDIF
C
CPARAFIN
C
C CALCUL DU TABLEAU NULONE
C
      DO 50 K1=1,NPTFR
C
CPARAFIN
      K2=KP1BOR(K1,1)
      IEL = NELBOR(K1)
      N1  = NBOR(K1)
      N2  = NBOR(K2)
C
      I1 = 0
      I2 = 0
C
      DO 60 IPT=1,NPT
C
        IF(IKLE(IEL,IPT).EQ.N1) THEN
          NULONE(K1,1) = IPT
          I1 = 1
        ENDIF
        IF(IKLE(IEL,IPT).EQ.N2) THEN
          NULONE(K1,2) = IPT
          I2 = 1
        ENDIF
C
60    CONTINUE
C
      IF(I1.EQ.0.OR.I2.EQ.0) THEN
        IF(LNG.EQ.1) WRITE(LU,810) IEL
        IF(LNG.EQ.2) WRITE(LU,811) IEL
810     FORMAT(1X,'ELEBD: ERREUR DE NUMEROTATION DANS L''ELEMENT:',I6,/,
     *         1X,'       CAUSE POSSIBLE :                       '   ,/,
     *         1X,'       LE FICHIER DES CONDITIONS AUX LIMITES NE'  ,/,
     *         1X,'       CORRESPOND PAS AU FICHIER DE GEOMETRIE  ')
811     FORMAT(1X,'ELEBD: ERROR OF NUMBERING IN THE ELEMENT:',I6,
     *         1X,'       POSSIBLE REASON:                       '   ,/,
     *         1X,'       THE BOUNDARY CONDITION FILE IS NOT      '  ,/,
     *         1X,'       RELEVANT TO THE GEOMETRY FILE           ')
        CALL PLANTE(1)
        STOP
      ENDIF
C
50    CONTINUE
C
C  CALCUL DE IKLBOR : LIKE IKLE FOR BOUNDARY POINTS, WITH BOUNDARY
C                     POINTS NUMBERING
C
      DO K=1,NPTFR
        IKLBOR(K,1) = K
        IKLBOR(K,2) = KP1BOR(K,1)
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
C                       *****************
                        SUBROUTINE STOSEG
C                       *****************
C
     *(IFABOR,NELEM,NELMAX,NELMAX2,IELM,IKLE,NBOR,NPTFR,
     * GLOSEG,MAXSEG,ELTSEG,ORISEG,NSEG,KP1BOR,NELBOR,NULONE)
C
C***********************************************************************
C BIEF VERSION 5.9         02/10/08    J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C    FONCTION : BUILDING THE DATA STRUCTURE FOR EDGE-BASED STORAGE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    IFABOR      |<-- | TABLEAU DES VOISINS DES FACES.
C |    NELEM       | -->| NOMBRE D'ELEMENTS DANS LE MAILLAGE.
C |    NELMAX      | -->| NOMBRE MAXIMUM D'ELEMENTS DANS LE MAILLAGE.
C |                |    | (CAS DES MAILLAGES ADAPTATIFS)
C |    NELMAX2     | -->| PREMIERE DIMENSION DE IFABOR
C |                |    | (EN 3D LE NOMBRE D'ELEMENTS 2D)
C |    IELM        | -->| 11: TRIANGLES.
C |                |    | 21: QUADRILATERES.
C |    IKLE        | -->| NUMEROS GLOBAUX DES POINTS DE CHAQUE ELEMENT.
C |    NBOR        | -->| GLOBAL NUMBERS OF BOUNDARY POINTS.
C |    NPTFR       | -->| NUMBER OF BOUNDARY POINTS.
C |    GLOSEG      |<-- | GLOBAL NUMBERS OF POINTS OF SEGMENTS.
C |    MAXSEG      |<-- | 1st DIMENSION OF MAXSEG.
C |    ELTSEG      |<-- | SEGMENTS OF EVERY TRIANGLE.
C |    ORISEG      |<-- | ORIENTATION OF SEGMENTS OF EVERY TRIANGLE.
C |    NSEG        |<-- | NUMBER OF SEGMENTS OF THE MESH.
C |    KP1BOR      | -->| NUMBER OF POINT FOLLOWING BOUNDARY POINT K
C |                |    | (I.E. K+1 MOST OF THE TIME BUT NOT ALWAYS).
C |    NELBOR      | -->| NUMBER OF ELEMENT CONTAINING SEGMENT K OF
C |                |    | THE BOUNDARY.
C |    NULONE      | -->| LOCAL NUMBER OF BOUNDARY POINTS IN A BOUNDARY
C |                |    | ELEMENT.
C |________________|____|_______________________________________________
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C CALLING PROGRAMME : INBIEF
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER NELMAX,NELMAX2,NPTFR,NSEG,MAXSEG,IELM,NELEM
      INTEGER NBOR(NPTFR),KP1BOR(NPTFR)
      INTEGER IFABOR(NELMAX2,3),IKLE(NELMAX,3)
      INTEGER NELBOR(NPTFR),NULONE(NPTFR)
      INTEGER GLOSEG(MAXSEG,2)
      INTEGER ELTSEG(NELMAX,3),ORISEG(NELMAX,3)     
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IPTFR,NSE
C
      INTEGER NEL,IFA,I1,I2,J1,J2,IFACE,JFACE,IG1,IG2
      INTEGER IELEM,IELEM1,IELEM2
C
      INTEGER NEXT(3)
      DATA NEXT / 2,3,1 /
C
C-----------------------------------------------------------------------
C
      IF(IELM.NE.11.AND.IELM.NE.12.AND.IELM.NE.13.AND.IELM.NE.14) THEN
        IF (LNG.EQ.1) WRITE(LU,500) IELM
        IF (LNG.EQ.2) WRITE(LU,501) IELM
500     FORMAT(1X,'STOSEG (BIEF) : ELEMENT NON PREVU : ',1I6)
501     FORMAT(1X,'STOSEG (BIEF) : UNEXPECTED ELEMENT: ',1I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C     INITIALISING ELTSEG
C
      DO IELEM = 1 , NELEM
        ELTSEG(IELEM,1) = 0
        ELTSEG(IELEM,2) = 0
        ELTSEG(IELEM,3) = 0
      ENDDO
C
C-----------------------------------------------------------------------
C
C     LOOP ON BOUNDARY POINTS : 
C
      NSE = 0
      DO IPTFR = 1 , NPTFR
C
C       IN PARALLELISM, IF THE BOUNDARY POINT FOLLOWING IPTFR IS IN
C       ANOTHER SUB-DOMAIN, KP1BOR(IPTFR)=IPTFR 
C       IN THIS CASE THE SEGMENT
C       BASED ON IPTFR AND THIS POINT IS NOT IN THE LOCAL DOMAIN
C       A CONSEQUENCE IS THAT NSE IS NOT EQUAL TO IPTFR
C
        IF(KP1BOR(IPTFR).NE.IPTFR) THEN
C
          NSE = NSE + 1
          GLOSEG(NSE,1) = NBOR(IPTFR)
          GLOSEG(NSE,2) = NBOR(KP1BOR(IPTFR))
          NEL = NELBOR(IPTFR)
          IFA = NULONE(IPTFR)
          ELTSEG(NEL,IFA) = NSE
          ORISEG(NEL,IFA) = 1
C
        ENDIF
C
      ENDDO
C
C-----------------------------------------------------------------------
C
C     LOOP ON ELEMENTS FOR NUMBERING INTERNAL SEGMENTS AND FILLING:
C     GLOSEG, ELTSEG, ORISEG
C
      DO IELEM1 = 1 , NELEM
        DO IFACE = 1 , 3
          IF(ELTSEG(IELEM1,IFACE).EQ.0) THEN
C           NEW SEGMENT (HENCE INTERNAL SO IFABOR<>0)
            NSE = NSE + 1
C           BOTH NEIGHBOURING ELEMENTS ARE TREATED FOR THIS SEGMENT
            I1 = IKLE(IELEM1,     IFACE)
            I2 = IKLE(IELEM1,NEXT(IFACE))
            IF(I1.EQ.I2) THEN
              IF(LNG.EQ.1) THEN
               WRITE(LU,*) 'STOSEG : SEGMENT AVEC UN SEUL POINT'
               WRITE(LU,*) '         ELEMENT ',IELEM1,' FACE ',IFACE               
              ENDIF
              IF(LNG.EQ.2) THEN
               WRITE(LU,*) 'STOSEG: EDGE MADE OF ONLY ONE POINT'
               WRITE(LU,*) '        ELEMENT ',IELEM1,' FACE ',IFACE               
              ENDIF
              CALL PLANTE(1)
              STOP
            ENDIF
            ELTSEG(IELEM1,IFACE) = NSE
            IG1=I1
            IG2=I2
C           SEGMENT ORIENTED LOWER RANK TO HIGHER RANK
            IF(IG1.LT.IG2) THEN
              GLOSEG(NSE,1) = I1
              GLOSEG(NSE,2) = I2
              ORISEG(IELEM1,IFACE) = 1
            ELSE
              GLOSEG(NSE,1) = I2
              GLOSEG(NSE,2) = I1
              ORISEG(IELEM1,IFACE) = 2
            ENDIF
C           OTHER ELEMENT NEIGHBOURING THIS SEGMENT
            IELEM2 = IFABOR(IELEM1,IFACE)
C           IELEM2 = 0 OR -1 MAY OCCUR IN PARALLELISM
            IF(IELEM2.GT.0) THEN
C             LOOKING FOR THE RIGHT FACE OF ELEMENT IELEM2
              DO JFACE = 1,3
                J1 = IKLE(IELEM2,     JFACE)
                J2 = IKLE(IELEM2,NEXT(JFACE))
C               ALL ELEMENTS HAVE A COUNTER-CLOCKWISE NUMBERING
                IF(I1.EQ.J2.AND.I2.EQ.J1) THEN
                  ELTSEG(IELEM2,JFACE) = NSE
                  ORISEG(IELEM2,JFACE) = 3-ORISEG(IELEM1,IFACE)
C                 FACE FOUND, NO NEED TO GO ON
                  GO TO 1000
                ELSEIF(I1.EQ.J1.AND.I2.EQ.J2) THEN
C                 FACE MAL ORIENTEE
                  IF(LNG.EQ.1) THEN
                    WRITE(LU,*) 'STOSEG : MAILLAGE DEFECTUEUX'
                    WRITE(LU,*) '         LA FACE ',JFACE
                    WRITE(LU,*) '         DE L''ELEMENT ',IELEM2
                    WRITE(LU,*) '         EST MAL ORIENTEE' 
                    WRITE(LU,*) '         (POINTS ',I1,' ET ',I2,')'           
                  ENDIF
                  IF(LNG.EQ.2) THEN
                    WRITE(LU,*) 'STOSEG: WRONG MESH'
                    WRITE(LU,*) '        FACE ',JFACE
                    WRITE(LU,*) '        OF ELEMENT ',IELEM2
                    WRITE(LU,*) '        IS NOT WELL ORIENTED' 
                    WRITE(LU,*) '         (POINTS ',I1,' AND ',I2,')'            
                  ENDIF
                  CALL PLANTE(1)
                  STOP                  
                ENDIF
              ENDDO
C             FACE NOT FOUND, THIS IS AN ERROR
              IF(LNG.EQ.1) THEN
                WRITE(LU,*) 'STOSEG : MAILLAGE DEFECTUEUX'
                WRITE(LU,*) '         ELEMENTS ',IELEM1,' ET ',IELEM2 
                WRITE(LU,*) '         LIES PAR LES POINTS ',I1,' ET ',I2
                WRITE(LU,*) '         MAIS CES POINTS NE FONT PAS UNE'
                WRITE(LU,*) '         FACE DE L''ELEMENT ',IELEM2               
              ENDIF
              IF(LNG.EQ.2) THEN
                WRITE(LU,*) 'STOSEG: WRONG MESH'
                WRITE(LU,*) '        ELEMENTS ',IELEM1,' AND ',IELEM2 
                WRITE(LU,*) '        LINKED BY POINTS ',I1,' AND ',I2
                WRITE(LU,*) '        BUT THESE POINTS ARE NOT AN EDGE'
                WRITE(LU,*) '        OF ELEMENT ',IELEM2               
              ENDIF
              CALL PLANTE(1)
              STOP
            ENDIF
1000        CONTINUE
          ENDIF
        ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C
C     VERIFICATION
C
      IF(NSEG.NE.NSE) THEN
        IF (LNG.EQ.1) WRITE(LU,502) NSE,NSEG
        IF (LNG.EQ.2) WRITE(LU,503) NSE,NSEG
502     FORMAT(1X,'STOSEG (BIEF) : MAUVAIS NOMBRE DE SEGMENTS : ',1I6,
     *            '                AU LIEU DE ',1I6,' ATTENDUS')
503     FORMAT(1X,'STOSEG (BIEF): WRONG NUMBER OF SEGMENTS : ',1I6,
     *            '               INSTEAD OF ',1I6,' EXPECTED')
        CALL PLANTE(1)
        STOP     
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
