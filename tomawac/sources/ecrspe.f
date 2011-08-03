C                       *****************
                        SUBROUTINE ECRSPE
C                       *****************
C
     *( F     , B     , TETA  , NPLAN , FREQ  , NF    , NK    ,
     *  NPOIN2, AT    , LT    , AUXIL , INUTIL, NOLEO , NLEO  , NSCO  ,
     *  BINSCO, DEBRES, TITCAS, DATE  , TIME ,ISLEO ,KNOLG)
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 6.0 ----                         C
C                                                                     C
C  ECRSPE :  PREPARATION DES IMPRESSIONS DES SPECTRES DIRECTIONNELS   C
C  ********  DE VARIANCE AUX POINTS D'ESPACE SELECTIONNES             C
C            (FORMAT SERAPHIN BINAIRE)                                C
C                                                                     C
C   - CREE    POUR VERSION 5.0  LE 28/08/00   OPTIMER  02 98 44 24 51 C
C   - MODIFIE POUR VERSION 5.2  LE 07/06/01                           C
C   - MODIFIE POUR VERSION 5.5  LE 13/07/2004 M. Benoit (Correction   C
C     d'un bug sur la declaration de IPOBO tels qu'il etait passe en  C
C     argument. Le tableau passe n'est plus utilise et renomme en     C
C     INUTIL. IPOBO est declare en ALLOCATABLE dans la subroutine)    C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL DE VARIANCE           !  C
C  ! TETA(-)     ! -->! VECTEUR DES DIRECTIONS DE DISCRETISATION   !  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !  C
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !  C
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL 2D    !  C
C  ! AT          ! -->! DATE COURANTE DU CALCUL  
!BD_ICNKA
!  ! LT          ! -->! NUMERO D IMRPESSION                        !  C   
!BD_INCKA
C  ! AUXIL(-,-)  !<-->! TABLEAU DE TRAVAIL POUR SPECTRE DIRECT.    !  C
C  ! INUTIL(-)   !<-->! TABLEAU NON-UTILISE (ANCIEN IPOBO)         !  C
C  ! NOLEO(-)    ! -->! NUMEROS DE NOEUD DES POINTS DE SORTIE      !  C
C  ! NLEO        ! -->! NOMBRE DE POINTS DE SORTIE                 !  C
C  ! NSCO        ! -->! NUM. DU FICHIER DE SORTIE DES SPECTRES     !  C
C  ! BINSCO      ! -->! BINAIRE DU FICHIER DE SORTIE DES SPECTRES  !  C
C  ! DEBRES      ! -->! INDICATEUR DE PREMIERE DATE A SAUVER       !  C
C  ! TITCAS      ! -->! TITRE DU CAS DE CALCUL                     !  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  WAC-COW                    C
C  ********    - PROGRAMME(S) APPELE(S) :  ECRIT                      C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
      USE BIEF
      USE TOMAWAC_MPI
      USE TOMAWAC_MPI_TOOLS
C
      IMPLICIT NONE

      COMMON/ECRSPE_MPI/SPE_SEND
      INTEGER ::SPE_SEND
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NPOIN2, NLEO  , NSCO  , NF , NK , NPLAN
      INTEGER  NOLEO(NLEO)   , INUTIL(1)
      INTEGER  DATE(3),TIME(3)
      DOUBLE PRECISION AT    , AUXIL(NPLAN,NK), AAT(1)
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF) , TETA(NPLAN), FREQ(NF)
      DOUBLE PRECISION B(NPOIN2,NPLAN)
      INTEGER  LT
C
      LOGICAL DEBRES
      CHARACTER*72 TITCAS
      CHARACTER(LEN=*)   BINSCO
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  ISTAT , II    , JF    , K     , IB(10)
      INTEGER  KAMP1 , KAMP2 , KAMP3 , KAMP4 , KAMP5 , KAMP6 , ILEO
      INTEGER  IBID(1), NELEM, NPSPE
      CHARACTER*72 C
CMB...Modif MB/GM nb de points de sortie
      CHARACTER*32 TEXTE(99)
CMB...Fin de modif
      CHARACTER*6  NUM
      CHARACTER*2  CC
      CHARACTER*1  C0    , C1    , C2    , C3    , C4    , C5    , C6
!BD_INCKA ajout d un maillage sur les frequences et les plans
      TYPE(BIEF_MESH) :: MESHF
CMB...Modif MB/GM nb de points de sortie      
      LOGICAL         :: SORLEO(99)
CMB...Fin de modif
      INTEGER :: I,IER
      INTEGER, ALLOCATABLE :: IKLE(:) ! GLOBAL CONNECTIVITY
      TYPE(BIEF_OBJ)  :: BVARSOR
!BD_INCKA modif dans le cas parallel
      INTEGER, DIMENSION(NLEO) :: NRECV_LEO
      LOGICAL         :: ISLEO(NLEO) !
      INTEGER, DIMENSION(NCSIZE) :: NSPE_RECV
      INTEGER :: KNOLG(NPOIN2)
!BD_INCKA
!BD_INCKA fin modif
C
      NPSPE=NF*NPLAN
      NELEM=(NF-1)*NPLAN
      SORLEO = .FALSE.
      DO ILEO=1,NLEO
          KAMP1=NOLEO(ILEO)
          IF (NCSIZE.GT.1) KAMP1=KNOLG(NOLEO(ILEO))
          KAMP2=MOD(KAMP1,100000)
          KAMP3=MOD(KAMP2,10000)
          KAMP4=MOD(KAMP3,1000)
          KAMP5=MOD(KAMP4,100)
          KAMP6=MOD(KAMP5,10)
CMB...Modif MB/GM nb de points de sortie      
CMB       IF(ILEO.GT.9) THEN
CMB         C0='1'
CMB       ELSE
CMB         C0='0'
CMB       ENDIF
CMB       CC=C0//CHAR(48+MOD(ILEO,10))
          CC=CHAR(48+INT(ILEO/10))//CHAR(48+MOD(ILEO,10))
CMB...Fin de modif	  
          C1=CHAR(48+INT(KAMP1/100000))
          C2=CHAR(48+INT(KAMP2/10000))
          C3=CHAR(48+INT(KAMP3/1000))
          C4=CHAR(48+INT(KAMP4/100))
          C5=CHAR(48+INT(KAMP5/10))
          C6=CHAR(48+KAMP6)
          NUM=C1//C2//C3//C4//C5//C6
          TEXTE(ILEO)='F'//CC//' PT2D'//NUM//'  UNITE SI       '
          IF (.not.ISLEO(ILEO)) TEXTE(ILEO) = 
     *     'POINT HORS MAILLAGE             '
          SORLEO(ILEO) = .TRUE.
      ENDDO
C
C=====C
C  1  C AU PREMIER PAS DE TEMPS IMPRIME, ON ECRIT L'EN-TETE DU FICHIER.
C=====C================================================================
      IF (DEBRES) THEN
C
C.......2.1 CONSTRUCTION DU NOM DES VARIABLES
C       """""""""""""""""""""""""""""""""""""
!BD_INCKA Creation de meshf, maillage associe a la discretisation
!frequence direction 
        ALLOCATE(MESHF%TYPELM)
        ALLOCATE(MESHF%NELEM)
        ALLOCATE(MESHF%NPOIN)
        ALLOCATE(MESHF%IKLE)
        ALLOCATE(MESHF%IKLE%I(4*NELEM))
        ALLOCATE(MESHF%X)
        ALLOCATE(MESHF%Y)
        ALLOCATE(MESHF%NPTFR)
        ALLOCATE(MESHF%NBOR)
        ALLOCATE(MESHF%NBOR%I(NPSPE))
        ALLOCATE(MESHF%DIM)
        ALLOCATE(MESHF%KNOLG)
        ALLOCATE(MESHF%KNOLG%I(NPSPE))
        MESHF%NAME = 'MESH'
        MESHF%TYPELM = 20 !maillage quadrangle 2D
        MESHF%NELEM  = NELEM
        MESHF%NPOIN  = NPSPE
        MESHF%DIM    = 2
        ALLOCATE(IKLE(4*NELEM))
        II=0
        DO JF=1,NF-1
          DO K=1,NPLAN
           II=II+1
           IKLE(II)=MOD(II,NPLAN)+1+(JF-1)*NPLAN
          ENDDO
        ENDDO
        DO II=1,NELEM
          IKLE(II+NELEM)=II
          IKLE(II+2*NELEM)=II+NPLAN
          IKLE(II+3*NELEM)=IKLE(II)+NPLAN
        ENDDO
        MESHF%IKLE%I=IKLE
C
!        DEALLOCATE(IKLE)
!        DEALLOCATE(IPOBO)
C
C.......2.9 ECRITURE DES TABLEAUX X ET Y
C       """"""""""""""""""""""""""""""""
        ALLOCATE(MESHF%X%R(NPLAN*NF))
        ALLOCATE(MESHF%Y%R(NPLAN*NF))
        MESHF%NPTFR = 2*NPLAN!+2*(NF-2)
        DO JF=1,NF
          DO II=1,NPLAN
            MESHF%X%R(II+NPLAN*(JF-1))=FREQ(JF)*SIN(TETA(II))
          ENDDO
        ENDDO
        DO JF=1,NF
          DO II=1,NPLAN
            MESHF%Y%R(II+NPLAN*(JF-1))=FREQ(JF)*COS(TETA(II))
          ENDDO
        ENDDO
        MESHF%NBOR%I=0
        DO II = 1,NPLAN
           MESHF%NBOR%I(II) = II
        ENDDO
        DO II = NPLAN+1,2*NPLAN
           MESHF%NBOR%I(II)=NPLAN+1+NPSPE-II
        ENDDO
        MESHF%KNOLG%I = 0
!BD_INCKA dans le cas // on initialise les valeurs des parametres MPI
      IF (NCSIZE.GT.1)         CALL GET_MPI_PARAMETERS(MPI_INTEGER,
     *                          MPI_REAL8,MPI_UB,
     *                          MPI_COMM_WORLD,MPI_SUCCESS)
!BD_INCKA fin modif
!BD_INCKA dnas le cas // on attribue les numero des noeud des sous
! domaines aux points d'ecriture de spectre
        IF (NCSIZE.GT.1) THEN
        CALL SPECTRE_SEND(SPE_SEND,NSPE_RECV,NLEO,ISLEO,
     *                           NRECV_LEO)
        CALL TEXTE_SENDRECV(TEXTE(1:NLEO),NLEO,NPSPE,ISLEO,NRECV_LEO)
        ENDIF
!BD_INCKA fin modif
        IF (((NCSIZE.GT.1).AND.(IPID==0)).OR.NCSIZE.LE.1) THEN
        ! CREATION DU JEU DE DONNEES POUR UN FORMAT DE FICHIER
        ! FORMAT_RES.
        ! LE JEU DE DONNEES EST CREE DANS LE FICHIER NRES, ET EST
        ! DEFINIT PAR UN TITRE ET LES VARIABLES A ECRIRE. 
        CALL CREATE_DATASET(BINSCO, ! FORMAT FICHIER RESULTAT
     *                      NSCO,  ! LU FICHIER RESULTAT
     *                      TITCAS,     ! TITRE DE L'ETUDE
     *                      NLEO,        ! MAX VARIABLES SORTIE
     *                      TEXTE,      ! NOMS VARIABLES SORTIE
     *                      SORLEO)     ! SORTIE OU PAS DES VARIABLES
        ! ECRITURE DU MAILLAGE DANS LE FICHIER SORTIE :
        ! SI ON EST ON PARALLEL, FAUT L'INDIQUER VIA NCSIZE ET NPTIR.
        ! LES AUTRES INFORMATIONS SONT DANS MESH.
        ! EN PLUS : DATE/TEMPS DE DEPART ET LES COORDONNEES DE
        ! L'ORIGINE.
        CALL WRITE_MESH(BINSCO, ! FORMAT FICHIER RESULTAT     
     *                  NSCO,  ! LU FICHIER RESULTAT
     *                  MESHF,          ! DESCRIPTEUR MAILLAGE
     *                  1,             ! NOMBRE DE PLAN /NA/
     *                  DATE,          ! DATE DEBUT
     *                  TIME,          ! HEURE DEBUT
     *                  0,0)  ! COORDONNEES ENTIERES DE L'ORIGINE.
        ENDIF
C
      ENDIF
!BD_INCKA dnas le cas // on attribue les numero des noeud des sous
! domaines aux points d'ecriture de spectre
        IF (NCSIZE.GT.1) THEN
        CALL SPECTRE_SEND(SPE_SEND,NSPE_RECV,NLEO,ISLEO,
     *                           NRECV_LEO)
        CALL TEXTE_SENDRECV(TEXTE,NLEO,NPSPE,ISLEO,NRECV_LEO)
        ENDIF
!BD_INCKA fin modif
C=====C
C  3  C ENREGISTREMENT DU PAS DE TEMPS COURANT.
C=====C========================================
C
C.....3.1 ECRITURE DU TEMPS AT
C     """"""""""""""""""""""""
      AAT(1) = AT
      ALLOCATE(BVARSOR%ADR(NLEO))
      DO II=1,NLEO
        ALLOCATE(BVARSOR%ADR(II)%P)
        ALLOCATE(BVARSOR%ADR(II)%P%R(NPSPE))
        BVARSOR%ADR(II)%P%DIM1 = NPSPE
        BVARSOR%ADR(II)%P%ELM  = 21
      ENDDO
      DO ILEO=1,NLEO
        II=NOLEO(ILEO)
        DO JF=1,NF
          DO K=1,NPLAN
            BVARSOR%ADR(ILEO)%P%R(K+(JF-1)*NPLAN)=F(II,K,JF)
          ENDDO
        ENDDO
      ENDDO
!
!BD_INCKA modif cas parallel
      IF (NCSIZE.GT.1) THEN
        CALL BVARSOR_SENDRECV(BVARSOR,NLEO,NPSPE,ISLEO,NRECV_LEO)
CMB...Modif MB/GM nb de points de sortie
        IF (IPID==0) CALL WRITE_DATA(BINSCO,NSCO,99,AT,LT,SORLEO,TEXTE,
     *                BVARSOR,NPSPE)
CMB...Fin de modif
      ELSE
CMB...Modif MB/GM nb de points de sortie
        CALL WRITE_DATA(BINSCO,NSCO,99,AT,LT,SORLEO,TEXTE,
     *                BVARSOR,NPSPE)
CMB...Fin de modif
      ENDIF 
!
!
      RETURN
      END
