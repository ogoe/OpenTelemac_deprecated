C
C 24/06/2009
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
      PROGRAM GREDELHYD
C
C     MERGES THE RESULTS OF A PARALLEL COMPUTATION WITH COUPLING WITH DELWAQ
C     TO WRITE A SINGLE FILE IN DELWAQ FORMAT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
      INTEGER LI
C
      CHARACTER(LEN=30) GEO   
C
      INTEGER ERR
      INTEGER NELEM,ECKEN,NDUM,I,J,NBV1,NBV2,PARAM(10)
      INTEGER NPLAN,NPOIN2
      INTEGER NPROC
      integer i_s, i_sp, i_len
      INTEGER IDUM, NPTFR,NSEG2,MBND
C
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: LIHBOR             ! LIHBOR(NPTFR)
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: NBOR               ! NBOR(*)
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: NBOR0,LIHBOR0      ! NBOR0(NPTFR),LIHBOR0(NPTFR)
C
      REAL RDUM
      REAL,    DIMENSION(:)  , ALLOCATABLE :: F
C
      LOGICAL IS
C
      CHARACTER*30 RES
      CHARACTER*50 RESPAR
      CHARACTER*11 EXTENS
      CHARACTER*30 CONLIM
      EXTERNAL    EXTENS
      INTRINSIC MAXVAL
C
      INTEGER ITSTRT,ITSTOP,ITSTEP,NSTEPA
      INTEGER MARDAT(3),MARTIM(3)
      CHARACTER*72  TITRE
      CHARACTER*144 NOMGEO,NOMLIM
      CHARACTER*144 NOMSOU,NOMMAB,NOMCOU,NOMSAL,NOMTEM
      CHARACTER*144 NOMINI,NOMVEB,NOMMAF,NOMVEL,NOMVIS
      LOGICAL SALI_DEL,TEMP_DEL
      LOGICAL VELO_DEL,DIFF_DEL
C
      LI=5
      LU=6
      LNG=2
chw
cjaj introduce yourself with the version date
c
      write(lu,*) 'I am Gredelhyd, cousin of Gretel from BAW Hamburg' 
      write(lu,*)
C
      write (lu, advance='no', 
     &    fmt='(/,'' Global geometry file: '')')
!      REWIND(LI)      
      read(li,*) geo
      write(lu,*) geo
c
c reading file names and the number of processors / partitions
c
      write (lu, advance='no', fmt='(/,'' Result file: '')')
      read(li,*) res   
      write(lu,*) res
c
      write (lu,advance='no',fmt='(/,'' Number of processors: '')')
      read (li,*) nproc
      write(lu,*) nproc

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
      OPEN(3,FILE=RES,FORM='FORMATTED',ERR=991)    
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
      OPEN(4,FILE=RESPAR,FORM='FORMATTED',ERR=994)
      GO TO 995
994   WRITE(LU,*) 'ERROR WHEN OPENING FILE: ',RESPAR
      call plante(-1)
      STOP
995   CONTINUE
C
      READ(4,'(I6)')NPLAN
      CLOSE(4)
C
      ALLOCATE(F(NPLAN),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'f')
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
C----------------------------------------------------------------------
C
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
C
      MBND=0
C
      DO I=1,NPTFR
        NBOR(I)   = NBOR0(I)
        LIHBOR(I) = LIHBOR0(I)
        IF (LIHBOR(I).NE.2) THEN
          MBND = MBND + 1
        ENDIF
      ENDDO
C
C     WITH PRISMS, DIFFERENT FROM 2D VALUES, OTHERWISE
C
      NSEG2 = (3*NELEM+NPTFR)/2
C
c
      OPEN(4,FILE=RESPAR,FORM='FORMATTED',ERR=984)
      GO TO 985
984   WRITE(LU,*) 'ERROR WHEN OPENING FILE: ',RESPAR
      call plante(-1)
      STOP
985   CONTINUE
C
      READ(4,'(I6)')NPLAN
      READ(4,'(I3)')J
      READ(4,'(A)')TITRE(1:J)
      READ(4,'(I4)')MARDAT(1)
      READ(4,'(I2)')MARDAT(2)
      READ(4,'(I2)')MARDAT(3)
      READ(4,'(I2)')MARTIM(1)
      READ(4,'(I2)')MARTIM(2)
      READ(4,'(I2)')MARTIM(3)
      READ(4,'(I14)')ITSTRT
      READ(4,'(I14)')ITSTOP
      READ(4,'(I14)')NSTEPA
      READ(4,'(I6)')NPLAN
C
      WRITE(3, '(A)' )
     *    "task      full-coupling                              "
      WRITE(3, '(A)' )
     *    "                                                     "
      WRITE(3, '(A)' )
     *    "#                                                    "
      WRITE(3, '(A)' )
     *    "# Telemac data                                       "
      WRITE(3, '(A)' )
     *    "#                                                    "
      WRITE(3, '(A)' )
     *    "                                                     "
      WRITE(3, '(A)' )
     *    "geometry  finite-elements                            "
      WRITE(3, '(A)' )
     *    "                                                     "
      WRITE(3, '(A)' )
     *    "horizontal-aggregation       no                      "
      WRITE(3, '(A)' )
     *    "minimum-vert-diffusion-used  no                      "
      WRITE(3, '(A)' )
     *    "vertical-diffusion           calculated              "
      WRITE(3, '(A)' )
     *    "description                                          "
      WRITE(3, '(A,A,A)' )
     *    "   '",TITRE(1:J),"'"
      WRITE(3, '(A)' )
     *    "   '                                    '            "
      WRITE(3, '(A)' )
     *    "   '                                    '            "
      WRITE(3, '(A)' )
     *    "end-description                                      "
      WRITE(3, '(A,I4,I2,I2,I2,I2,I2,A)' )
     *"reference-time           '",MARDAT(1),MARDAT(2),MARDAT(3),
     *                             MARTIM(1),MARTIM(2),MARTIM(3),"'"
      WRITE(3, '(A,I14,A)' )
     *    "hydrodynamic-start-time  '",ITSTRT,"'"
      WRITE(3, '(A,I14,A)' )
     *    "hydrodynamic-stop-time   '",ITSTOP,"'"
      WRITE(3, '(A,I14,A)' )
     *    "hydrodynamic-timestep    '",NSTEPA,"'"
      WRITE(3, '(A,I14,A)' )
     *    "conversion-ref-time      '",ITSTRT,"'"
      WRITE(3, '(A,I14,A)' )
     *    "conversion-start-time    '",ITSTRT,"'"
      WRITE(3, '(A,I14,A)' )
     *    "conversion-stop-time     '",ITSTOP,"'"
      WRITE(3, '(A,I14,A)' )
     *    "conversion-timestep      '",NSTEPA,"'"
      WRITE(3, '(A,I6)'  )
     *    "grid-cells-first-direction ",NPOIN2
      WRITE(3, '(A,I6,A)')
     *    "grid-cells-second-direction",NSEG2+MBND," # nr of exchanges!"
      WRITE(3, '(A,I6)' )
     *    "number-hydrodynamic-layers ",NPLAN
      WRITE(3, '(A,I6)' )
     *    "number-water-quality-layers",NPLAN
      READ(4,'(I3)')J
      READ(4,'(A)')NOMGEO(1:J)
      WRITE(3, '(A,A,A)' )
     *    "hydrodynamic-file        '",NOMGEO(1:J),"'"
      WRITE(3, '(A)' )
     *    "aggregation-file         none                        "
      WRITE(3, '(A,A,A)' )
     *    "grid-indices-file        '",NOMGEO(1:J),"'"
      READ(4,'(I3)')J
      READ(4,'(A)')NOMLIM(1:J)
      WRITE(3, '(A,A,A)' )
     *    "grid-coordinates-file    '",NOMLIM(1:J),"'"
      READ(4,'(I3)')J
      READ(4,'(A)')NOMSOU(1:J)
      WRITE(3, '(A,A,A)' )
     *    "volumes-file             '",NOMSOU(1:J),"'"
      READ(4,'(I3)')J
      READ(4,'(A)')NOMMAB(1:J)
      WRITE(3, '(A,A,A)' )
     *    "areas-file               '",NOMMAB(1:J),"'"
      READ(4,'(I3)')J
      READ(4,'(A)')NOMCOU(1:J)
      WRITE(3, '(A,A,A)' )
     *    "flows-file               '",NOMCOU(1:J),"'"
      READ(4,'(I3)')J
      READ(4,'(A)')NOMVEB(1:J)
      WRITE(3, '(A,A,A)' )
     *    "pointers-file            '",NOMVEB(1:J),"'"
      READ(4,'(I3)')J
      READ(4,'(A)')NOMMAF(1:J)
      WRITE(3, '(A,A,A)' )
     *    "lengths-file             '",NOMMAF(1:J),"'"
      READ(4,'(L)')SALI_DEL
      IF(SALI_DEL) THEN
        READ(4,'(I3)')J
        READ(4,'(A)')NOMSAL(1:J)
        WRITE(3, '(A,A,A)' )
     *    "salinity-file            '",NOMSAL(1:J),"'"
      ELSE
        WRITE(3, '(A)' )
     *    "salinity-file            none                        "
      ENDIF
      READ(4,'(L)')TEMP_DEL
      IF(TEMP_DEL) THEN
        READ(4,'(I3)')J
        READ(4,'(A)')NOMTEM(1:J)
        WRITE(3, '(A,A,A)' )
     *    "temperature-file         '",NOMTEM(1:J),"'"
      ELSE
        WRITE(3, '(A)' )
     *    "temperature-file         none                        "
      ENDIF
      READ(4,'(L)')DIFF_DEL
      IF(DIFF_DEL) THEN
        READ(4,'(I3)')J
        READ(4,'(A)')NOMVIS(1:J)
        WRITE(3, '(A,A,A)' )
     *    "vert-diffusion-file      '",NOMVIS(1:J),"'"
      ELSE
        WRITE(3, '(A)' )
     *    "vert-diffusion-file      none                        "
      ENDIF
      READ(4,'(L1)') VELO_DEL
      IF(VELO_DEL) THEN
        READ(4,'(I3)')J
        READ(4,'(A)')NOMVEL(1:J)
        WRITE(3, '(A,A,A)' )
     *    "velocity-file            '",NOMVEL(1:J),"'"
      ELSE
        WRITE(3, '(A)' )
     *    "velocity-file            none                        "
      ENDIF
      READ(4,'(I3)')J
      READ(4,'(A)')NOMINI(1:J)
      WRITE(3, '(A,A,A)' )
     *    "surfaces-file            '",NOMINI(1:J),"'"
C
      WRITE(3, '(A)' )
     *    "total-grid-file          none                        "
      WRITE(3, '(A)' )
     *    "discharges-file          none                        "
      WRITE(3, '(A)' )
     *    "chezy-coefficients-file  none                        "
      WRITE(3, '(A)' )
     *    "shear-stresses-file      none                        "
      WRITE(3, '(A)' )
     *    "walking-discharges-file  none                        "
      if ( NPLAN .gt. 1 ) then
         WRITE(3, '(A)' )
     *       "minimum-vert-diffusion                            "
         WRITE(3, '(A)' )
     *       "   upper-layer       0.0000E+00                   "
         WRITE(3, '(A)' )
     *       "   lower-layer       0.0000E+00                   "
         WRITE(3, '(A)' )
     *       "   interface-depth   0.0000E+00                   "
         WRITE(3, '(A)' )
     *       "end-minimum-vert-diffusion                        "
      endif
      WRITE(3, '(A)' )
     *    "constant-dispersion                                  "
      WRITE(3, '(A)' )
     *    "   first-direction    0.0000                         "
      WRITE(3, '(A)' )
     *    "   second-direction   0.0000                         "
      WRITE(3, '(A)' )
     *    "   third-direction    0.0000                         "
      WRITE(3, '(A)' )
     *    "end-constant-dispersion                              "
      WRITE(3, '(A)' )
     *    "hydrodynamic-layers                               "
      DO I=1,NPLAN
        READ(4,'(F10.4)')F(I)
      ENDDO
      do I=1,NPLAN
         WRITE(3, '(F10.4)' ) F(I)
      enddo
      WRITE(3, '(A)' )
     *    "end-hydrodynamic-layers                           "
      WRITE(3, '(A)' )
     *    "water-quality-layers                              "
      do I=1,NPLAN
         WRITE(3, '(F10.4)' ) 1.0
      enddo
      WRITE(3, '(A)' )
     *    "end-water-quality-layers                          "
      WRITE(3, '(A)' )
     *    "discharges                                           "
      WRITE(3, '(A)' )
     *    "end-discharges                                       "
C
      WRITE(LU,*) 'END OF PROGRAM '
C
      CLOSE(2)       
      CLOSE(3)
      CLOSE(4)
C
      STOP

      END PROGRAM GREDELHYD


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
