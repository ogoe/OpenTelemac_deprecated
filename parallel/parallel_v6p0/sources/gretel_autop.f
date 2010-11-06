C
C 08/08/2007
C
C Warning by JMH: there is a CALL EXIT(ICODE) which is a Fortran extension
C                 it will not work with some compilers, like Nag
C
C
C Note JMH 28/04/2009: there is a strange treatment here of NPLAN=0
C                      which seems to correspond to 2D whereas 
C                      NPLAN.NE.0 would be 3D. It corresponds to
C                      a special programming of ECRGEO
C                      and WRITE_MESH_SERAFIN. This is to be changed.    
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
      PROGRAM GRETEL
C
C     REGROUPE LES RESULTATS D'UN CALCUL PARALLELE POUR EN FAIRE
C     UN FICHIER AU FORMAT SELAFIN NORMAL.
C
      IMPLICIT NONE
      INTEGER LNG,LU
      INTEGER LI
      COMMON/INFO/LNG,LU
C     
C
      CHARACTER(LEN=30) GEO   
C      	
C-------------------------------------------------------------------------
C
      LI=5
      LU=6
      LNG=2
chw
cjaj introduce yourself with the version date
c
      write(lu,*) 'I am Gretel from BAW Hamburg' 
      write(lu,*) 'reincarnated by Holger Weilbeer' 
      write(lu,*) 'on 20th February 2003'
      write(lu,*)
      
        !!!!Fabs: Test pour savoir dans quelle configuration on se trouve:
        !!!!      ** Si le premier nom est T2DGEO, alors ===> decomposition de domaine TELEMAC-2D
        !!!!      OU
        !!!!      ** Si le premier nom est different de T2DGEO, et egal a E2DSERA ou E2DVOL ou E2DSCAL,
        !!!!      alors ===> decomposition particulaire ESTEL-2D.

        !!!!Fabs: Il faudrait effectuer des tests de rapidite sur les machines, quelle difference
        !!!!    : entre SCAL et PARAL... Le mode SCAL semble etre plus rapide et prendre moins de
        !!!!    : place. 

        !!!!Fabs: La decomposition particulaire ne modifie pas le SERAPHIN, donc tou interet a le
        !!!!Fabs: declarer juste en SCAL plutot qu en SELAFIN => gani de memoire et gani de temps au
        !!!!Fabs: regrouepement des resultats. 
C	
C
      write (lu, advance='no', 
     &    fmt='(/,'' Global geometry file: '')')
!      REWIND(LI)      
      read(li,*) geo
 !       write(lu,*) geo     
C      	
      IF ((GEO.EQ.'E2DSERA').OR.(GEO.EQ.'E2DVOL')
     &       .OR.(GEO.EQ.'E2DSCAL')) THEN 
C              
C|=======================================================================/     
C|     	                                                                 /
C| DEBUT DU REGROUPEMENT DES FICHIERS ISSUS DU PARALLELE PARTICULAIRE    /
C|                                                                       /
C| SERAFIN  = INCHANGE => SIMPLE RECOPIE D UN DES NPROC FICHIERS         /
C| VOLFIN   = CHANGE => SOMME AU NIVEAU DES RESULTATS PART_INS,PART_CUM  /
C| SCALAIRE = CHANGE => SOMME AU NIVEAU DES RESULTATS NBPART_LOST...     /
C|                                                                       /     
C|=======================================================================/
C
C
      CALL RECOMPOSITION_PARTICULAIRE(GEO)
C      
C      
C|==================================================================|     
C|     	                                                            |
C| FIN DU REGROUPEMENT DES FICHIERS ISSUS DU PARALLELE PARTICULAIRE |
C|                                                                  |
C|==================================================================|
C
C
        ELSE 
C   
C   
C|==================================================================|     
C|     	                                                            |
C| DEBUT DU REGROUPEMENT DES FICHIERS ISSUS DU PARALLELE DE DOMAINE |
C|                                                                  |
C|==================================================================|
C
        CALL RECOMPOSITION_DECOMP_DOMAINE (GEO)
C   
C|==================================================================|     
C|     	                                                            |
C| DEBUT DU REGROUPEMENT DES FICHIERS ISSUS DU PARALLELE DE DOMAINE |
C|                                                                  |
C|==================================================================|
C	


      ENDIF       
      
      STOP

      END PROGRAM GRETEL
!
!=============================================================================      
!

C
C			  ****************************************
                          SUBROUTINE RECOMPOSITION_PARTICULAIRE (GEO)
C                         ****************************************
C
C***********************************************************************
C BIEF VERSION ?           18/07/05    Fabien DECUNG (Stagiaire MATMECA) 
C***********************************************************************
C
C FONCTION : Regroupement des resultats issus du calcul en decomposition 
C	     particulaire dans ESTEL-2D. Les seuls fichiers regroupes sont:
C
C	     1/ Serafin (E2DSERA): aucune modification. Simple recopie
C	        d'un seul fichier.        
C	     2/ Volfin (E2DVOL): Variables particulaires modifiees.
C	        Il s'agit de sommer les valeurs stockees sur les fichiers de
C	        chaque processeur. 		
C	     3/ Resultats scalaires (E2DSCAL): Variables particulaires modifiees.
C	        Il s'agit de sommer les valeurs 
C	        
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      IKLE      |<-->|  TABLEAU DES CONNECTIVITES                   |
C |      NELEM     | -->|  NOMBRE D'ELEMENTS
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS
C |      NPOIN     | -->|  NOMBRE DE SOMMETS DU MAILLAGE
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : gretel_autop.f
C
C  SOUS-PROGRAMME APPELE : extens, alloer, cpikle2, SKIP_HEADER, READ_DATASET, READ_DATASET_ELEM
C
C**********************************************************************
C                         ****************************************

      IMPLICIT NONE
      INTEGER LNG,LU
      INTEGER LI
      COMMON/INFO/LNG,LU
!
!=>Fabs
      CHARACTER(LEN=30), INTENT(IN) :: GEO
!<=Fabs
      INTEGER IPID,ERR,FU
      INTEGER NELEM,ECKEN,NDUM,I,J,K,NBV1,NBV2,PARAM(10)
      INTEGER NPLAN,NPOIN2,NPOIN2LOC
      INTEGER NPROC,NRESU,NPOINMAX
      INTEGER I_S, I_SP, I_LEN
!
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: NPOIN,IPOBO,VERIF,IPOBO3D
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KNOLG
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IKLESA,IKLE3D
      
      !Fabs------------------------------------------------------!
      INTEGER, DIMENSION(:)    , ALLOCATABLE   :: PART
!=>Fabs : bug Nag
      INTEGER, DIMENSION(:)    , ALLOCATABLE   :: PART_REP
!<=Fabs
      REAL   , DIMENSION(:,:)  , ALLOCATABLE   :: GLOBAL_VALUE
      REAL   , DIMENSION(:,:)  , ALLOCATABLE   :: LOCAL_VALUE
      REAL   , DIMENSION(:)    , ALLOCATABLE   :: XORIG,YORIG      
      REAL   , DIMENSION(:)    , ALLOCATABLE   :: SOMMEPART      
      REAL   , DIMENSION(:,:,:), ALLOCATABLE   :: LOCAL_VALUELEM      
      DOUBLE PRECISION, DIMENSION (:)    , ALLOCATABLE :: SOMMERESU
      DOUBLE PRECISION, DIMENSION (:,:,:), ALLOCATABLE :: VALUESCP
      !Fabs------------------------------------------------------!      
!
      !Fabs---------------------!	
      DOUBLE PRECISION SOMME
      !!!!Fabs: Attention a ne pas mettre en double precision, sinon buuuuug!!!   
      REAL AT
      INTEGER NUM_PROC      
      INTEGER NBVAR,LINE,NBLINE,TEMPS
      LOGICAL IS,ENDE,ENDEOFFILE
      !Fabs---------------------!      
!
!      CHARACTER*32 GEO,GEOM
!=>Fabs
      CHARACTER*32 GEOM
!<=Fabs

      CHARACTER*32 RUB,RUBENS      
      CHARACTER*32 TEXTLU(200)      
      CHARACTER*72 TITRE
      CHARACTER*72 PTEXCL, TITLE      
      CHARACTER*80 TITSEL
!	
      CHARACTER*11 EXTENS
      EXTERNAL    EXTENS
      INTEGER, INTRINSIC ::  MAXVAL
!            
!
!   
!
        !!!!Fabs: Li = Numero du fichier de lecture des noms de fichiers et 
        !!!!    : du nombre de processeurs => gretel.par (repertoire temporaire)           
        LI = 5
!	
        !!!!Fabs: Information obligatoire car en decomposition particulaire, le fichier de geometrie	
        !!!!    : n est pas declare en "selafin-geom" mais en "paral" et donc, n est pas ecrit dans gretel.par
!	
        write(lu, advance='no', fmt=
     &         '(/,''  PARTIAL RESULT SELAFIN FILE: <input_name>: '')')
        read(LI,*) RUBENS
!		
        if (RUBENS.eq.' ') then
          write (lu,'('' no filename'')') 
        else
          write(lu,*) 'input: ',RUBENS 
        end if        
!	
        WRITE (lu,advance='no',fmt='(/,'' Number of processors: '')')
        READ (LI,*) Nproc
        WRITE (LU,*) Nproc

        !!!!Fabs: Il faut tester si les fichiers sont des seraphin ou volfin ou scp
        !!!!    : car on a divers methodes de regroupement selon le type de fichiers.
      
        IF ((RUBENS.EQ.'E2DSERA').OR.(RUBENS.EQ.'E2DVOL')) THEN
C
C     FICHIER DE GEOMETRIE DU CALCUL, LU JUSQU'AUX 10 PARAMETRES:
C

!!!!Fabs: Uniquement si le fichier GEO est declare en PARAL dans le dictionnaire
!Fabs-----------------------------------------!
!	i_s  = len (GEO)
!        i_sp = i_s + 1
!        do i=1,i_s
!         if(GEO(i_sp-i:i_sp-i) .ne. ' ') exit
!        enddo
!        i_len=i_sp - i
!	
!	GEOM=GEO(1:i_len) // extens(nproc-1,0)
!!!!Fabs: Sinon on prend la racine E2DGEO
        GEOM = GEO
!Fabs-----------------------------------------!

      OPEN(2,FILE=GEOM,FORM='UNFORMATTED',STATUS='OLD',ERR=9990)
      READ(2,ERR=9990)
      READ(2,ERR=9990) NBV1,NBV2
      DO 110 I=1,NBV1+NBV2
        READ(2,ERR=9990)
110    CONTINUE
      GO TO 9992
9990   WRITE(LU,*) 'ERROR WHEN OPENING OR READING FILE: ',GEOM
      call plante(-1)
      STOP
9992   CONTINUE
C     LECTURE DES 10 PARAMETRES ET DE LA DATE
      READ(2) (PARAM(I),I=1,10)
      IF(PARAM(10).EQ.1) READ(2) (PARAM(I),I=1,6)
C
C     FICHIER  DE RESULTATS :
C
      OPEN(3,FILE=RUBENS,FORM='UNFORMATTED',ERR=9991)
      GO TO 9993
9991   WRITE(LU,*) 'ERROR WHEN OPENING FILE: ',RUBENS
      call plante(-1)
      STOP
9993   CONTINUE
C
C     1) LECTURE DU DEBUT DU PREMIER FICHIER DE RESULTATS.
C
C
        i_s  = len (RUBENS)
        i_sp = i_s + 1
        do i=1,i_s
         if(RUBENS(i_sp-i:i_sp-i) .ne. ' ') exit
        enddo
        i_len=i_sp - i
!	
        RUB=RUBENS(1:i_len) // extens(nproc-1,0)
!
      !!!!Fabs: Pour tester au moins si le premier fichier temporaire parallele existe
      inquire (file=RUB,exist=is)
      if (.not.is) then 
        write (lu,*) 'file does not exist: ', RUB
        write (lu,*) 'check the number of processors'
        write (lu,*) 'and the result file core name'
        call plante(-1)
        stop
      end if  
!
      OPEN(4,FILE=RUB,FORM='UNFORMATTED',ERR=9994)
      GO TO 9995
9994  WRITE(LU,*) 'ERROR WHEN OPENING FILE: ',RUB
      call plante(-1)
      STOP
9995   CONTINUE
!
!       
!
C
C  1 : TITRE
C
      READ(4) TITRE
      WRITE(LU,*) 'TITLE=',TITRE
      TITSEL=TITRE // 'SERAPHIN'
      WRITE(3) TITSEL
C
C  2 : NBV1,NBV2
C	
      READ(4) NBV1,NBV2
      WRITE(LU,*) 'NBV1=',NBV1,'   NBV2=',NBV2
      WRITE(3) NBV1,NBV2
C
C  3 : NOMS ET UNITES DES VARIABLES
C
!
      !!!!Fabs: Allocation du pointeur particulaire
!
      ALLOCATE (PART(1:2))
      PART = 0
!	
      DO 5500 I=1,NBV1     
        READ(4) TEXTLU(I)
        WRITE(LU,*) 'VARIABLE ',I,' : ',TEXTLU(I)

        !!!!FABS: REPERAGE DU NUMERO DES VARIABLES A SOMMER
        !!!!      POUR PART_INS => PART(1)
        !!!!      POUR PART_CUM => PART(2)

      IF (RUBENS.EQ.'E2DVOL') THEN     
        IF ( (TEXTLU(I).EQ.'PARTICULES INST. -')
     &  .OR.(TEXTLU(I).EQ.'PARTICLES INST. -') ) THEN
           PART(1) = I
        ELSE
           IF ( (TEXTLU(I).EQ.'PARTICULES CUM.  -' ) 
     &  .OR.(TEXTLU(I).EQ.'PARTICLES CUM.  -') ) THEN
              PART(2) = I
           ENDIF
       ENDIF
      ENDIF
      WRITE(3) TEXTLU(I)
            
5500  CONTINUE
!
!      IF ((RUBENS.EQ.'E2DVOL').AND.
!     &   ((PART(1).EQ.0).OR.(PART(2).EQ.0))) THEN
!      	WRITE(LU,*) 'PAS RESULTATS PART INS OU CUM'
!	WRITE(LU,*) 'VERIFIER VOS SORTIES PARTICULAIRES VOLFIN'
!	CALL PLANTE(-1)
!      ENDIF
C
C  4 : 10 PARAMETRES
C
      READ(4) (PARAM(I),I=1,10)
      WRITE(LU,*) '10 PARAMETERS : ',PARAM
      PARAM(9)=0
      PARAM(8)=0
      NPLAN=PARAM(7)
      WRITE(3) (PARAM(I),I=1,10)
C LECTURE EVENTUELLE DE LA DATE ET ECRITURE
      IF(PARAM(10).EQ.1) THEN
        READ(4)  (PARAM(I),I=1,6)
        WRITE(3) (PARAM(I),I=1,6)
      ENDIF
      
      CLOSE(4) 
      
C
C  5 : LECTURE DES VARIABLES NELEM :4 PARAMETRES
C       
      READ(2) NELEM,NPOIN2,ECKEN,NDUM
      WRITE(LU,*) '4 PARAMETERS IN GEOMETRY FILE'
      WRITE(LU,*) 'NELEM=',  NELEM
      WRITE(LU,*) 'NPOIN2=', NPOIN2
      WRITE(LU,*) 'ECKEN=',  ECKEN
      WRITE(LU,*) 'NDUM=',   NDUM
C
      IF(NPLAN.EQ.0) THEN
        WRITE(3) NELEM,NPOIN2,ECKEN,NDUM
      ELSE
        WRITE(3) NELEM*(NPLAN-1),NPOIN2*NPLAN,6,NDUM
      ENDIF
C
C  ALLOCATIONS DYNAMIQUES DES TABLEAUX
C
      ALLOCATE(NPOIN(NPROC),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'npoin')
      ALLOCATE(IKLESA(3,NELEM),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'iklesa')
      ALLOCATE(IPOBO(NPOIN2)      ,STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'ipobo')
      IF(NPLAN.EQ.0) THEN
        ALLOCATE(VERIF(NPOIN2)    ,STAT=ERR)
      ELSE
        ALLOCATE(VERIF(NPOIN2*NPLAN)    ,STAT=ERR)
      ENDIF
      IF(ERR.NE.0) call ALLOER (lu, 'verif')
C  GLOBAL_VALUES, STORING THE WHOLE DATASET (NBV1-VALUES)
      IF(NPLAN.EQ.0) THEN
       ! ALLOCATE(GLOBAL_VALUE(NPOIN2,NBV1)       ,STAT=ERR)
       ALLOCATE(GLOBAL_VALUE(NELEM,NBV1))
      ELSE
        ALLOCATE(GLOBAL_VALUE(NPOIN2*NPLAN,NBV1) ,STAT=ERR) 
      ENDIF
      IF(ERR.NE.0) call ALLOER (lu, 'global_value')
C  X AND Y
      ALLOCATE(XORIG(NPOIN2)    ,STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'xorig')
      ALLOCATE(YORIG(NPOIN2)    ,STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'yorig')
C  CONCERNE UN CAS 3D
      IF(NPLAN.NE.0) THEN
      ALLOCATE(IKLE3D(NELEM*(NPLAN-1),6),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'ikle3d')
      ALLOCATE(IPOBO3D(NPOIN2*NPLAN),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'ipobo3d')
      ENDIF
C
C  FIN ALLOCATION ...
C
C  6 : IKLE 
C 
      READ(2)  ((IKLESA(I,J),I=1,ECKEN),J=1,NELEM)
      WRITE(LU,*) 'WRITING IKLE'
      IF(NPLAN.EQ.0) THEN
        WRITE(3) ((IKLESA(I,J),I=1,ECKEN),J=1,NELEM)
      ELSE
C       ICI ECRITURE DE IKLE3D (AVEC INVERSION DES DIMENSIONS)
        CALL CPIKLE2(IKLE3D,IKLESA,NELEM,NELEM,NPOIN2,NPLAN)
        WRITE(3) ((IKLE3D(I,J),J=1,6),I=1,NELEM*(NPLAN-1))
      ENDIF
C
C  7 : IPOBO
C      
      READ(2)  (IPOBO(I),I=1,NPOIN2)
      WRITE(LU,*) 'WRITING IPOBO'
      IF(NPLAN.EQ.0) THEN
        WRITE(3) (IPOBO(I),I=1,NPOIN2)
      ELSE
C       DUMMY VALUES
        DO I=1,NPOIN2*NPLAN
          IPOBO3D(I) = 0
        ENDDO
        WRITE(3) (IPOBO3D(I),I=1,NPOIN2*NPLAN)
      ENDIF
C
C  8 : X et Y, WILL BE CHECKED LATER ...
C
      READ(2)  (XORIG(I),I=1,NPOIN2)
      READ(2)  (YORIG(I),I=1,NPOIN2)   
           
C 
C------------------------------------------------------------------------------
C
C OPENING FILES AND READING/SKIPPING HEADERS -> NPOIN(NPROC), NPOINMAX
C
      DO IPID = 0,NPROC-1
         !!!!Fabs: FU est le numero du fichier de resultats temporaires.
         !!!!    : Cela peut occasionner des erreurs dans des calculs sur 
         !!!!    : beaucoup de processeurs si deux processeurs ont la meme valeur FU.
         FU = IPID + 10
         RUB=RUBENS(1:I_LEN) // EXTENS(NPROC-1,IPID)
         OPEN (FU,FILE=RUB,FORM='UNFORMATTED',ERR=9998)
         GO TO 9999
9998     WRITE(LU,*) 'ERROR WHEN OPENING FILE: ',RUB,
     &                'USING FILE UNIT: ', FU
         CALL PLANTE(-1)
         STOP
9999      REWIND(FU)
         CALL SKIP_HEADER(FU,NPOIN(IPID+1),NBV1,ERR,lu)
         IF(ERR.NE.0) then 
           WRITE(LU,*) 'ERROR READING FILE'
           CALL PLANTE(-1)
           STOP
         ENDIF
      END DO
!      
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
        ALLOCATE(LOCAL_VALUE(NPOINMAX,NBV1)   ,STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'local_value')
      !!!!Fabs: Allocation d un nouveau vecteur temporaire de stockage des resultats
        ALLOCATE(LOCAL_VALUELEM(0:NPROC-1,1:NELEM,1:NBV1)   ,STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'LOCAL_VALUELEM')
      !!!!Fabs: Allocation d un vecteur temporaire de stockage des sommes
        ALLOCATE(SOMMEPART(1:NELEM)   ,STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'SOMMEPART')
!      
      !!!! JUSQU ICI OK POUR LE SERAFIN ET LE VOLFIN
C
C READING KNOLG(NPOIN,NPROC)
C
      !!!!Fabs: Utilite de stocker dans notre decomposition... Non!
C
      DO IPID = 0,NPROC-1
         FU = IPID + 10
         IF(NPLAN.EQ.0) THEN
            READ(FU) (KNOLG(I,IPID+1),I=1,NPOIN(IPID+1))
         ELSE
            READ(FU) (KNOLG(I,IPID+1),I=1,NPOIN(IPID+1)/NPLAN)
         ENDIF
      END DO
C
C READING LOCAL X
C
!
!
!
        !!!!Fabs: On lit et on avance sur tous les fichiers temporaires...
        !!!!    : Lecture des abscisses X
        DO IPID = 0,NPROC-1
        FU = IPID +10
        READ(FU) (LOCAL_VALUE(I,1),I=1,NPOIN(IPID+1)) 
        ENDDO
!	
!
        !!!!Fabs: On se fixe par rapport au processeur 0 (+1)	 
        !!!!    : et on ecrit les coordonnees des abscisses X.
        WRITE(LU,*) 'WRITING X-COORDINATES'
        IF(NPLAN.EQ.0) THEN
          WRITE(3) (LOCAL_VALUE(I,1),I=1,NPOIN(1))        
        ENDIF             
!
        !!!!Fabs: On lit et on avance sur tous les fichiers temporaires...
        !!!!    : Lecture des abscisses Y
        DO IPID = 0,NPROC-1
        FU = IPID +10
        READ(FU) (LOCAL_VALUE(I,1),I=1,NPOIN(IPID+1)) 
        ENDDO
!	
        !!!!Fabs: On se fixe par rapport au processeur 0 (+1)	 
        !!!!    : et on ecrit les coordonnees des abscisses X.
        WRITE(LU,*) 'WRITING Y-COORDINATES'
        IF(NPLAN.EQ.0) THEN
          WRITE(3) (LOCAL_VALUE(I,1),I=1,NPOIN(1))       
        ENDIF 
C
C READING DATASETS
C
      NRESU = 0
C
20000 NRESU = NRESU + 1
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
!
      !!!!Fabs: Test du type de fichier  :=> SERAFIN = E2DSERA
      !!!!    :                          :=> VOLFIN  = E2DVOL
      !!!!    :                          :=> SCP     = E2DSCAL
!                             
!	
      IF (RUBENS.EQ.'E2DSERA') THEN

                !!! Fabs: Le fichier en question est un Serafin.
                !!!     : Les valeurs sont alors È considÅrer par noeuds.
                !!!     : On ne considÉre alors qu'un seul processeur car 
                !!!     : elles sont les memes sur chaque processeur

                !!! Lecture des donnÅes sur le dernier fichier ouvert. 
                CALL READ_DATASET
     +   (LOCAL_VALUE,NPOINMAX,NPOIN(1),NBV1,AT,FU,ENDE)     
                IF (ENDE) GOTO 30000
                !!!Fabs: Ecriture des rÅsultats dans le fichier final
                WRITE(LU,*)'WRITING DATASET NO.',NRESU,' TIME =',AT
                WRITE(3) AT
                DO K = 1,NBV1
                        IF(NPLAN.EQ.0) THEN
                        WRITE(3) (LOCAL_VALUE(I,K),I=1,NPOIN(1))
                        ELSE
                        !Fabs: Utile dans notre cas...non?
                        !WRITE(3) (LOCAL_VALUE(I,K),I=1,NPOIN(1)*NPLAN)
                        ENDIF
                END DO
                GO TO 20000

        ELSE
          IF (RUBENS.EQ.'E2DVOL') THEN
            !!! Fabs: 
            !!! Le fichier en question est un Volfin, Part_Cum et Part_Ins
            !!! sont modifiÅes par la parallÅlisation.
            !!! Les valeurs sont ici donnÅes par mailles.
            !!! On doit selon les cas considÅrer un ou tous les processeurs.

            DO IPID = 0,NPROC-1 
              !!!! Fabs:
              !!!! On lit ici les donnÅes de chaque pas de temps
              !!!! pour toutes les variables et tous les processeurs. 	           
              FU = IPID +10
              CALL READ_DATASET_ELEM
     +   (LOCAL_VALUELEM,NPROC,NELEM,NBV1,AT,FU,IPID,ENDE)  
              IF (ENDE) GOTO 30000
            ENDDO
            !!!Fabs: Ecriture dans le fichier rÅsultat final.
            WRITE(LU,*)'WRITING DATASET NO.',NRESU,' TIME =',AT
            WRITE(3) AT
!	    
            !!!!Fabs: Ecriture dans le fichier rÅsultat final des rÅsultats. 
            DO K = 1,NBV1
              IF (NPLAN.EQ.0) THEN
                !!! Fabs: 
                !!! Dans le cas de donnÅes particulaires issues de
                !!! la parallÅlisation, il nous faut sommer celles-ci 
                !!! sur tous les processeurs. 
                IF ( (K.NE.PART(1).AND.K.NE.PART(2)).OR.AT.EQ.0) THEN
                !!! Fabs: Le paramÉtre K n'est pas dÅpendant de la 
                !!! parallÅlisation, il n'est alors utile que de prendre
                !!! les valeurs du processeur 0.  		            					
                  WRITE(3) (LOCAL_VALUELEM(0,I,K),I=1,NELEM)
                ELSE !( (K.NE.PART(1)...
                  !!! Fabs: DonnÅes issues de la parallÅlisation, on les somme 
                  !!! alors sur le nombre de processeurs.	
                  SOMMEPART = 0.
                  DO I= 1, NELEM
                    DO NUM_PROC = 0, NPROC-1
                  SOMMEPART(I)=SOMMEPART(I)+LOCAL_VALUELEM(NUM_PROC,I,K)
                    ENDDO
                  ENDDO
!		  IF (K.EQ.PART(2)) THEN 
                  !!! Fabs:
                  !!! il s agit du fichier PARTICULES CUM, il faut alors retirer
                  !!! les particules initialement injectÅes qui ne font pas partie de la
                  !!! parallelisation. 
!		  WRITE(3) ( SOMMEPART(I), I=1,NELEM )
!		  ELSE
                  WRITE(3) ( SOMMEPART(I), I=1,NELEM )
!		  ENDIF		  		  
                ENDIF !( (K.NE.PART(1)....
              ENDIF ! (NPLAN.EQ.0)
            ENDDO !(K = 1,NBV1)		    					  
          ENDIF
        GO TO 20000
      ENDIF ! (RUBENS.EQ.'E2DSERA')	    
!		
!		
30000  WRITE(LU,*) 'END OF PROGRAM, ',NRESU-1,' DATASETS FOUND'
!
      !!!!Fabs: Fermeture du fichier de lecture (2)
      !!!!    : et de celui de resultat final regroupe (3)
      !!!!    : ainsi que ceux temporaires sur chaque Proc.
      
      !Fabs--------------!
      CLOSE(2)       
      CLOSE(3)
!
      DO IPID = 0,NPROC-1
         FU = IPID +10
         CLOSE (FU)
      END DO      
      !Fabs--------------!      
!                              
      DEALLOCATE (PART)
      DEALLOCATE (LOCAL_VALUELEM)
!      
      
      ELSE !!!! IF (RUBENS.EQ.'E2DVOL').OR.(RUBENS.EQ.'E2DSERA')
      
      !!!!Fabs: Le fichier a lire est alors un fichier de resultats scalaires.
!
      OPEN(3,FILE=RUBENS,FORM='FORMATTED',ERR=99991)
      GO TO 99993
99991   WRITE(LU,*) 'ERROR WHEN OPENING FILE: ',RUBENS
      CALL PLANTE(-1)
      STOP
99993 CONTINUE
C
C     1) LECTURE DU DEBUT DU PREMIER FICHIER DE RESULTATS.
C
C
        i_s  = len (RUBENS)
        i_sp = i_s + 1
        do i=1,i_s
         if(RUBENS(i_sp-i:i_sp-i) .ne. ' ') exit
        enddo
        i_len=i_sp - i
!			
        RUB=RUBENS(1:i_len) // extens(nproc-1,0)
!	
      inquire (file=RUB,exist=is)
      if (.not.is) then 
        write (lu,*) 'file does not exist: ', RUB
        write (lu,*) 'check the number of processors'
        write (lu,*) 'and the result file core name'
        call plante(-1)
        stop
      end if  
C
C OPENING FILES AND READING/SKIPPING HEADERS -> NPOIN(NPROC), NPOINMAX
C
       ALLOCATE (PART_REP(1:8))       
       DO IPID = 0,NPROC-1
         FU = IPID + 10
         RUB=RUBENS(1:I_LEN) // EXTENS(NPROC-1,IPID)
         OPEN (FU,FILE=RUB,FORM='FORMATTED',ERR=99998)
         GO TO 99999
99998      WRITE(LU,*) 'ERROR WHEN OPENING FILE: ',RUB,
     &                 'USING FILE UNIT: ', FU
         call plante(-1)
         STOP
99999    REWIND(FU)
         !!!Fabs: On alloue ce tableau de 1 a 8 en raison des huits valeurs particulaires 
         !!!      a reconstituer.         	
         !!!!!! CHAUD !!!!!!
         ENDEOFFILE = .FALSE.
         NBVAR  = 0
         NBLINE = 0
         NBVAR  = 0
         DO WHILE (.NOT.ENDEOFFILE)
           READ(FU,1981,END = 9799) PTEXCL,TITLE
           IF (PTEXCL.EQ."'") THEN
             IF(IPID.EQ.0) WRITE(3,1981) "'",TITLE
             NBLINE = NBLINE + 1
             NBVAR  = NBVAR  + 1
             !!!! Fabs; il faut reperer les places des variables particulaires.
             IF (TITLE.EQ."NBPART          -'" ) PART_REP(1)=NBVAR-3
             IF (TITLE.EQ."NBPART_OUT      -'" ) PART_REP(2)=NBVAR-3
             IF (TITLE.EQ."NBPART_NEW      -'" ) PART_REP(3)=NBVAR-3
             IF (TITLE.EQ."NBPART_LOST     -'" ) PART_REP(4)=NBVAR-3
             IF (TITLE.EQ."NBPART_AT       -'" ) PART_REP(5)=NBVAR-3
             IF (TITLE.EQ."NBPART_OUT_AT   -'" ) PART_REP(6)=NBVAR-3
             IF (TITLE.EQ."NBPART_NEW_AT   -'" ) PART_REP(7)=NBVAR-3
             IF (TITLE.EQ."NBPART_LOST_AT  -'" ) PART_REP(8)=NBVAR-3

           ELSE
             !!!!Fabs:Le nombre de variables a recopier est trouve. 
             !!!!     On se replace au debut du fichier pour les ecrire.
             !!!!     NBVAR - 3: Il y a trois lignes de titre a chaque fois. 
             !!!!     On pourrait juste utiliser NBLINE - 3...	
             WRITE(LU,*) 'NUMBER OF VARIABLES', NBVAR-3
             IF ((PART_REP(1).NE.0.).AND.(PART_REP(8).NE.0.)) THEN    
               WRITE(LU,*) 'PARTICULAR VARIABLES FOUND'
             ELSE
               WRITE(LU,*) 'PARTICULAR VARIABLES NOT FOUND'
               GO TO 9799
             ENDIF
             REWIND(FU)
             GO TO 10190
           ENDIF
         END DO
         !!!!! Fabs: A modifier pitet bien que oui, pitet bien que non.
         !!!!! il y a actuellement trois lignes de titres!!!!	 
         !!!!! peut etre a modifier pour que ca prenne a partir de time?
         !!!!! pitet pas besoin vu que ces trois lignes sont ds le code: h2d_resscp.f

         !!!! On est oblige de stocker selon le pas de temps.

10190    TEMPS = 0
         !!!! Compte du nombre de pas de temps.
         DO WHILE (.NOT.ENDEOFFILE)
           READ(FU,*,END=6996) TITLE
           TEMPS = TEMPS + 1
         ENDDO 

6996     TEMPS = TEMPS - NBLINE
         WRITE(LU,*) 'NUMBER OF TIME STEPS', TEMPS

         IF (IPID.EQ.0) THEN
         ALLOCATE (VALUESCP(1:TEMPS+1,0:NPROC-1,1:NBVAR-3))
         ALLOCATE (SOMMERESU(1:8))
         VALUESCP=0.
         ENDIF

         !!!!Pour passer les lignes de titre et en venir au valeurs scalaires.
         REWIND(FU)
         DO LINE=1,NBLINE
           READ(FU,*) TITLE
         ENDDO
         !!!!Lecture de toutes les valeurs scalaires par processeurs et par pas de temps.
         TEMPS = 0
         DO WHILE (.NOT.ENDEOFFILE)
           READ(FU,*,END=6969) (VALUESCP(TEMPS+1,IPID,I),I=1,NBVAR-3)
           TEMPS = TEMPS + 1
         ENDDO
6969     CONTINUE
       ENDDO
       !!!! Fabs: Il nous faut maintenant reconstituer les bonnes sommes!!!
       !!!! Fabs: C chaud!!!!
       !!!! Recopiage des donnees non particulaires a partir du processeur 0	     	 
       DO LINE=1,TEMPS         
         SOMMERESU = 0.
         DO IPID=0,NPROC-1                   
           !!!!Les variables scalaires 1,3,5,7 ne changent pas
           !!!!Tandis qu il faut sommer les 2,4,6,8  
           DO I=1,7,2
             SOMMERESU(I) =  VALUESCP(LINE,0,PART_REP(I))
           ENDDO
           DO I=2,8,2
             SOMMERESU(I) = SOMMERESU(I)+VALUESCP(LINE,IPID,PART_REP(I))
           ENDDO
         ENDDO                       
         !!!! Condition necessaire s il reste des variables apres celles issues du particulaire sinon, else...	   	   	  	 
         IF ((PART_REP(8)+1).LE.NBVAR-3) THEN
         !!!! A tester!!!!
         WRITE(3,1010) (VALUESCP(LINE,0,I),I=1,PART_REP(1)-1),
     &   (SOMMERESU(I),I=1,8),
     &   (VALUESCP(LINE,0,I),I=PART_REP(8)+1,SIZE(VALUESCP,3))
         ELSE
         WRITE(3,1010) (VALUESCP(LINE,0,I),I=1,PART_REP(1)-1),
     &   (SOMMERESU(I),I=1,8)
         ENDIF
       ENDDO  
       DO LINE=1,TEMPS
       WRITE(LU,*) 'Time',line
       DO IPID=0,NPROC-1
       ENDDO
       ENDDO
       GO TO 97909

9799    PRINT*, 'ERROR'        
        CALL PLANTE(-1)
        STOP

97909   PRINT*, 'DATA SETS FOUND'
        DO IPID = 0,NPROC-1
        FU = IPID +10
        CLOSE (FU)
        END DO 
        CLOSE (3)
        DEALLOCATE (VALUESCP)
!
1981    FORMAT (A1,A70)
1010    FORMAT(E14.6,1X,30(E14.6,1X))
      
      END IF !!!! IF (RUBENS.EQ.'E2DVOL').OR.(RUBENS.EQ.'E2DSERA')
!        
      STOP  
!
      END SUBROUTINE RECOMPOSITION_PARTICULAIRE  
!
!=============================================================================      
!
!=============================================================================      
!
!=============================================================================      
!      
      
C			  *************************************************
                          SUBROUTINE RECOMPOSITION_DECOMP_DOMAINE (GEO)
C                         *************************************************
      IMPLICIT NONE
      INTEGER LNG,LU
      INTEGER LI
      COMMON/INFO/LNG,LU    
C     
!=>Fabs
      CHARACTER(LEN=30), INTENT(IN) :: GEO
!<=Fabs

      INTEGER IPID,ERR,FU
      INTEGER NELEM,ECKEN,NDUM,I,J,K,NBV1,NBV2,PARAM(10)
      INTEGER NPLAN,NPOIN2,NPOIN2LOC
      INTEGER NPROC,NRESU,NPOINMAX
      integer i_s, i_sp, i_len
C
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: NPOIN,IPOBO,VERIF,IPOBO3D
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KNOLG
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IKLESA,IKLE3D
C      
C
      REAL   , DIMENSION(:,:), ALLOCATABLE :: GLOBAL_VALUE
      REAL   , DIMENSION(:,:), ALLOCATABLE :: LOCAL_VALUE      
      REAL   , DIMENSION(:)  , ALLOCATABLE :: XORIG,YORIG      
C
      REAL AT
C
      LOGICAL IS,ENDE
C
      CHARACTER*30 RES
      CHARACTER*50 RESPAR
      CHARACTER*72 TITRE
      CHARACTER*80 TITSEL
      CHARACTER*32 TEXTLU(200)
      CHARACTER*11 EXTENS
      EXTERNAL    EXTENS
      INTEGER, INTRINSIC ::  MAXVAL           
!
!      
        LI = 5      
c
c reading file names and the number of processors / partitions
c
!=>Fabs
!      read(li,*) geo
!        write(lu,*) geo
!      write (lu, advance='no', 
!     &    fmt='(/,'' Global geometry file: '')')
!<=Fabs
c
      write (lu, advance='no', fmt='(/,'' Result file: '')')
      read(li,*) res   
c
      write (lu,advance='no',fmt='(/,'' Number of processors: '')')
      read (li,*) nproc
      write (lu,*) ' '
      
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
C  1 : TITRE
C
      READ(4) TITRE
      WRITE(LU,*) 'TITLE=',TITRE
      TITSEL=TITRE // 'SERAFIN'
      WRITE(3) TITSEL
C
C  2 : NBV1,NBV2
C
      READ(4) NBV1,NBV2
      WRITE(LU,*) 'NBV1=',NBV1,'   NBV2=',NBV2
      WRITE(3) NBV1,NBV2
C
C  3 : NOMS ET UNITES DES VARIABLES
C
      DO 500 I=1,NBV1     
        READ(4) TEXTLU(I)
        WRITE(LU,*) 'VARIABLE ',I,' : ',TEXTLU(I)
        WRITE(3) TEXTLU(I)
500   CONTINUE
C
C  4 : 10 PARAMETRES
C
      READ(4) (PARAM(I),I=1,10)
      WRITE(LU,*) '10 PARAMETERS : ',PARAM
      PARAM(9)=0
      PARAM(8)=0
      NPLAN=PARAM(7)
      WRITE(3) (PARAM(I),I=1,10)
C LECTURE EVENTUELLE DE LA DATE ET ECRITURE
      IF(PARAM(10).EQ.1) THEN
        READ(4)  (PARAM(I),I=1,6)
        WRITE(3) (PARAM(I),I=1,6)
      ENDIF
      CLOSE(4)
C
C  2) LECTURE DU FICHIER DE GEOMETRIE
C
C  5 : 4 parametres
C      
      READ(2) NELEM,NPOIN2,ECKEN,NDUM
      WRITE(LU,*) '4 PARAMETERS IN GEOMETRY FILE'
      WRITE(LU,*) 'NELEM=',NELEM
      WRITE(LU,*) 'NPOIN2=',NPOIN2
      WRITE(LU,*) 'ECKEN=',ECKEN
      WRITE(LU,*) 'NDUM=',NDUM
C
      IF(NPLAN.EQ.0) THEN
        WRITE(3) NELEM,NPOIN2,ECKEN,NDUM
      ELSE
        WRITE(3) NELEM*(NPLAN-1),NPOIN2*NPLAN,6,NDUM
      ENDIF
C
C  ALLOCATIONS DYNAMIQUES DES TABLEAUX
C
      ALLOCATE(NPOIN(NPROC),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'npoin')
      ALLOCATE(IKLESA(3,NELEM),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'iklesa')
      ALLOCATE(IPOBO(NPOIN2)      ,STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'ipobo')
      IF(NPLAN.EQ.0) THEN
        ALLOCATE(VERIF(NPOIN2)    ,STAT=ERR)
      ELSE
        ALLOCATE(VERIF(NPOIN2*NPLAN)    ,STAT=ERR)
      ENDIF
      IF(ERR.NE.0) call ALLOER (lu, 'verif')
C  GLOBAL_VALUES, STORING THE WHOLE DATASET (NBV1-VALUES)
      IF(NPLAN.EQ.0) THEN
        ALLOCATE(GLOBAL_VALUE(NPOIN2,NBV1)       ,STAT=ERR)
      ELSE
        ALLOCATE(GLOBAL_VALUE(NPOIN2*NPLAN,NBV1) ,STAT=ERR) 
      ENDIF
      IF(ERR.NE.0) call ALLOER (lu, 'global_value')
C  X AND Y
      ALLOCATE(XORIG(NPOIN2)    ,STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'xorig')
      ALLOCATE(YORIG(NPOIN2)    ,STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'yorig')
C  3D
      IF(NPLAN.NE.0) THEN
      ALLOCATE(IKLE3D(NELEM*(NPLAN-1),6),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'ikle3d')
      ALLOCATE(IPOBO3D(NPOIN2*NPLAN),STAT=ERR)
      IF(ERR.NE.0) call ALLOER (lu, 'ipobo3d')
      ENDIF
C
C  FIN ALLOCATION ...
C
C  6 : IKLE 
C 
      READ(2)  ((IKLESA(I,J),I=1,ECKEN),J=1,NELEM)
      WRITE(LU,*) 'WRITING IKLE'
      IF(NPLAN.EQ.0) THEN
        WRITE(3) ((IKLESA(I,J),I=1,ECKEN),J=1,NELEM)
      ELSE
C       ICI ECRITURE DE IKLE3D (AVEC INVERSION DES DIMENSIONS)
        CALL CPIKLE2(IKLE3D,IKLESA,NELEM,NELEM,NPOIN2,NPLAN)
        WRITE(3) ((IKLE3D(I,J),J=1,6),I=1,NELEM*(NPLAN-1))
      ENDIF
C
C  7 : IPOBO
C      
      READ(2)  (IPOBO(I),I=1,NPOIN2)
      WRITE(LU,*) 'WRITING IPOBO'
      IF(NPLAN.EQ.0) THEN
        WRITE(3) (IPOBO(I),I=1,NPOIN2)
      ELSE
C       DUMMY VALUES
        DO I=1,NPOIN2*NPLAN
          IPOBO3D(I) = 0
        ENDDO
        WRITE(3) (IPOBO3D(I),I=1,NPOIN2*NPLAN)
      ENDIF
C
C  8 : X et Y, WILL BE CHECKED LATER ...
C
      READ(2)  (XORIG(I),I=1,NPOIN2)
      READ(2)  (YORIG(I),I=1,NPOIN2)
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
         CALL SKIP_HEADER(FU,NPOIN(IPID+1),NBV1,ERR,lu)
         IF(ERR.NE.0) then 
           write(lu,*) 'ERROR READING FILE '
           call plante(-1)
           STOP
         ENDIF
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
        ALLOCATE(LOCAL_VALUE(NPOINMAX,NBV1)    ,STAT=ERR)
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
C READING LOCAL X
C
      DO IPID = 0,NPROC-1
         FU = IPID +10
         READ(FU) (LOCAL_VALUE(I,1),I=1,NPOIN(IPID+1))
         IF(NPLAN.EQ.0) THEN
          DO I=1,NPOIN(IPID+1)
            GLOBAL_VALUE(KNOLG(I,IPID+1),1) =
     +       LOCAL_VALUE(      I        ,1)
                   VERIF(KNOLG(I,IPID+1))   = 1
          ENDDO
         ELSE
          NPOIN2LOC = NPOIN(IPID+1)/NPLAN
          DO I=1,NPOIN2LOC
          DO J=1,NPLAN
            GLOBAL_VALUE(KNOLG(I,IPID+1) + NPOIN2   *(J-1) , 1)= 
     +       LOCAL_VALUE(      I         + NPOIN2LOC*(J-1) , 1)
                   VERIF(KNOLG(I,IPID+1) + NPOIN2   *(J-1))  = 1
          ENDDO
          ENDDO
         ENDIF    
      END DO
C
C COMPARISON WITH GLOBAL VALUES
C
C     IN 3D, CHECKING THE FIRST PLANE ONLY
      DO I=1,NPOIN2
         IF(ABS(XORIG(I)-GLOBAL_VALUE(I,1)).GT.0.1) THEN
            WRITE(LU,*) 'POINT ',I,' XORIG=',XORIG(I),
     +                ' GLOBAL_VALUE=',GLOBAL_VALUE(I,1)
            WRITE(LU,*) 'GEO IS PROBABLY NOT THE RIGHT ORIGINAL FILE'
         ENDIF
      ENDDO
C FURTHER VERIFICATIONS
      IF(NPLAN.EQ.0) THEN
        DO I=1,NPOIN2
          IF(VERIF(I).EQ.0) THEN
            WRITE(LU,*) 'ERROR, POINT I=',I,' FALSE FOR X-COORDINATES'
          ENDIF
        ENDDO
      ELSE
        DO I=1,NPOIN2*NPLAN
          IF(VERIF(I).EQ.0) THEN
            WRITE(LU,*) 'ERROR, POINT I=',I,' FALSE FOR X-COORDINATES'
          ENDIF
        ENDDO
      ENDIF     
C WRITING X
      WRITE(LU,*) 'WRITING X-COORDINATES'
      IF(NPLAN.EQ.0) THEN
         WRITE(3) (GLOBAL_VALUE(I,1),I=1,NPOIN2)
      ELSE
         WRITE(3) (GLOBAL_VALUE(I,1),I=1,NPOIN2*NPLAN)
      ENDIF
C
C READING LOCAL Y
C
      DO IPID = 0,NPROC-1
         FU = IPID +10
         READ(FU) (LOCAL_VALUE(I,1),I=1,NPOIN(IPID+1))
         IF(NPLAN.EQ.0) THEN
          DO I=1,NPOIN(IPID+1)
            GLOBAL_VALUE(KNOLG(I,IPID+1),1) =
     +       LOCAL_VALUE(      I        ,1)
                   VERIF(KNOLG(I,IPID+1))   = 1
          ENDDO
         ELSE
          NPOIN2LOC = NPOIN(IPID+1)/NPLAN
          DO I=1,NPOIN2LOC
          DO J=1,NPLAN
            GLOBAL_VALUE(KNOLG(I,IPID+1) + NPOIN2   *(J-1) , 1)= 
     +       LOCAL_VALUE(      I         + NPOIN2LOC*(J-1) , 1)
                   VERIF(KNOLG(I,IPID+1) + NPOIN2   *(J-1))  = 1
          ENDDO
          ENDDO
         ENDIF    
      END DO
C
C COMPARISON WITH GLOBAL VALUES
C
C IN 3D, CHECKING THE FIRST PLANE ONLY
      DO I=1,NPOIN2
         IF(ABS(YORIG(I)-GLOBAL_VALUE(I,1)).GT.0.1) THEN
            WRITE(LU,*) 'POINT ',I,' YORIG=',YORIG(I),
     +                      ' GLOBAL_VALUE=',GLOBAL_VALUE(I,1)
            WRITE(LU,*) 'GEO IS PROBABLY NOT THE RIGHT ORIGINAL FILE'
         ENDIF
      END DO
C FURTHER VERIFICATIONS
      IF(NPLAN.EQ.0) THEN
        DO I=1,NPOIN2
          IF(VERIF(I).EQ.0) THEN
            WRITE(LU,*) 'ERROR, POINT I=',I,' FALSE FOR Y-COORDINATES'
          ENDIF
        ENDDO
      ELSE
        DO I=1,NPOIN2*NPLAN
          IF(VERIF(I).EQ.0) THEN
            WRITE(LU,*) 'ERROR, POINT I=',I,' FALSE FOR Y-COORDINATES'
          ENDIF
        ENDDO
      ENDIF   
C WRITING Y
      WRITE(LU,*) 'WRITING Y-COORDINATES'
      IF(NPLAN.EQ.0) THEN
         WRITE(3) (GLOBAL_VALUE(I,1),I=1,NPOIN2)
      ELSE
         WRITE(3) (GLOBAL_VALUE(I,1),I=1,NPOIN2*NPLAN)
      ENDIF
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
     +   (LOCAL_VALUE,NPOINMAX,NPOIN(IPID+1),NBV1,AT,FU,ENDE)
         IF (ENDE) GOTO 3000
C STORE EACH DATASET
         IF(NPLAN.EQ.0) THEN
            DO I=1,NPOIN(IPID+1)
            DO K=1,NBV1
              GLOBAL_VALUE(KNOLG(I,IPID+1),K) = LOCAL_VALUE(I,K)
            END DO
              VERIF(KNOLG(I,IPID+1))   = 1
            END DO
         ELSE
            NPOIN2LOC = NPOIN(IPID+1)/NPLAN
            DO I=1,NPOIN2LOC
            DO J=1,NPLAN
            DO K=1,NBV1
            GLOBAL_VALUE(KNOLG(I,IPID+1) + NPOIN2   *(J-1) , K)= 
     +       LOCAL_VALUE(      I         + NPOIN2LOC*(J-1) , K)
            END DO
                   VERIF(KNOLG(I,IPID+1) + NPOIN2   *(J-1)) = 1
            END DO
            END DO
         ENDIF       
      END DO
C WRITING GLOBAL DATASET
      WRITE(LU,*)'WRITING DATASET NO.',NRESU,' TIME =',AT
C
      WRITE(3) AT
      DO K = 1,NBV1
         IF(NPLAN.EQ.0) THEN         
            WRITE(3) (GLOBAL_VALUE(I,K),I=1,NPOIN2)
         ELSE
            WRITE(3) (GLOBAL_VALUE(I,K),I=1,NPOIN2*NPLAN)
         ENDIF
      END DO
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
        !!!Fabs
      
      END SUBROUTINE RECOMPOSITION_DECOMP_DOMAINE


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
C                         **********************
                          SUBROUTINE SKIP_HEADER
C                         **********************
     +(FU,NPOIN,NVALUE,ERR,lu)
     
      !!!!!Fabs
      !!!!! Utilisee pour sauter l entete des fichiers issus du parallele
      !!!!!FinFabs
C
      IMPLICIT NONE
C
      INTEGER NPOIN,NELEM,ECKEN,NDUM,NBV1,NVALUE,I,NPLAN
      INTEGER FU,ERR,lu
      INTEGER PARAM(10)
C
C  1 : SKIPPING TITRE
C
      READ(FU,ERR=999)
C
C  2 : READINGING NBV1
C
      READ(FU,ERR=999) NBV1
      IF (NBV1.NE.NVALUE) then 
        write(lu,*)  'NBV1.NE.NVALUE! CHECK OUTPUT FILES ...'
        call plante(-1)
        STOP
      endif 
C
C  3 : SKIPPING NOMS ET UNITES DES VARIABLES
C
      DO I=1,NBV1     
        READ(FU,ERR=999)
      END DO
C
C  4 : 10 PARAMETRES
C
      READ(FU,ERR=999) (PARAM(I),I=1,10)
      NPLAN=PARAM(7)
C  LECTURE EVENTUELLE DE LA DATE ET ECRITURE
      IF(PARAM(10).EQ.1) THEN
        READ(FU,ERR=999)  (PARAM(I),I=1,6)
      ENDIF
C
C  5 : 4 parametres
C      
      READ(FU,ERR=999) NELEM,NPOIN,ECKEN,NDUM
C
C  6 : IKLE 
C 
      READ(FU,ERR=999)
C
 999  RETURN
      END
C                         ***********************
                          SUBROUTINE READ_DATASET
C                         ***********************
     +(LOCAL_VALUE,NPOINMAX,NPOIN,NVALUE,AT,FU,ENDE)
C
      IMPLICIT NONE
C
      INTEGER NPOINMAX,NPOIN,NVALUE,FU
      INTEGER IPOIN,IVALUE
C
      REAL AT
      REAL LOCAL_VALUE(NPOINMAX,NVALUE)
C
      LOGICAL ENDE
C
      ENDE = .TRUE.
C
      READ(FU,END=999) AT      
      DO IVALUE = 1,NVALUE
         READ(FU,END=999) (LOCAL_VALUE(IPOIN,IVALUE),IPOIN=1,NPOIN)
      END DO
C
      ENDE = .FALSE.
C
 999  RETURN
      END
      
!!!!FABS : G RAJOUTE CETTE SUBROUTINE....      
C      			  ****************************
                          SUBROUTINE READ_DATASET_ELEM
C                         ****************************
     +(LOCAL_VALUELEM,NPROC,NELEM,NBV1,AT,FU,IPID,ENDE)          
C
      IMPLICIT NONE
C
      INTEGER NPROC,NELEM,NBV1,FU,IPID
      INTEGER IELEM,IVALUE
C
      REAL AT     
      REAL LOCAL_VALUELEM(0:NPROC-1,1:NELEM,1:NBV1)
C
      LOGICAL ENDE
C
      ENDE = .TRUE.
C
      READ(FU,END=9099) AT      
      DO IVALUE = 1,NBV1
         READ(FU,END=9099) (LOCAL_VALUELEM(IPID,IELEM,IVALUE)  
     &   ,IELEM=1,NELEM)
      END DO
C
      ENDE = .FALSE.
C
 9099  RETURN
      END                            
      
C                       ******************
                        SUBROUTINE CPIKLE2
C                       ******************
C
     *(IKLE3,IKLES,NELEM2,NELMAX2,NPOIN2,NPLAN)
C
C***********************************************************************
C BIEF VERSION 5.1           23/08/99    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : EXTENSION DU TABLEAU DES CONNECTIVITES
C            CAS DE L'EXTENSION A UN ELEMENT QUASI-BULLE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      IKLE      |<-->|  TABLEAU DES CONNECTIVITES                   |
C |      NELEM     | -->|  NOMBRE D'ELEMENTS
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS
C |      NPOIN     | -->|  NOMBRE DE SOMMETS DU MAILLAGE
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : INBIEF
C
C  SOUS-PROGRAMME APPELE : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NELEM2,NELMAX2,NPOIN2,NPLAN
      INTEGER, INTENT(INOUT) :: IKLES(3,NELEM2)
      INTEGER, INTENT(INOUT) :: IKLE3(NELMAX2,NPLAN-1,6)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,I
C
C-----------------------------------------------------------------------
C
C     BOTTOM AND TOP OF ALL LAYERS
C
      IF(NPLAN.GE.2) THEN
        DO I = 1,NPLAN-1
          DO IELEM = 1,NELEM2          
            IKLE3(IELEM,I,1) = IKLES(1,IELEM) + (I-1)*NPOIN2
            IKLE3(IELEM,I,2) = IKLES(2,IELEM) + (I-1)*NPOIN2
            IKLE3(IELEM,I,3) = IKLES(3,IELEM) + (I-1)*NPOIN2
            IKLE3(IELEM,I,4) = IKLES(1,IELEM) +  I   *NPOIN2
            IKLE3(IELEM,I,5) = IKLES(2,IELEM) +  I   *NPOIN2
            IKLE3(IELEM,I,6) = IKLES(3,IELEM) +  I   *NPOIN2  
          ENDDO
        ENDDO
      ELSE
        IF(LNG.EQ.1) WRITE(LU,*) 'CPIKLE2 : IL FAUT AU MOINS 2 PLANS'
        IF(LNG.EQ.2) WRITE(LU,*) 'CPIKLE2 : MINIMUM OF 2 PLANES NEEDED'
        call plante(-1)
        STOP
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
      !!! write(*,*) 'Returning exit code: ', icode
      CALL EXIT(ICODE)
      stop    ! which is usually equivalent to call EXIT(0)
      end subroutine PLANTE

