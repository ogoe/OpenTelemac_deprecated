C                       *****************
                        SUBROUTINE ECRGEO
C                       *****************
C
     *(X,Y,NPOIN,NBOR,NFIC,NVAR,TEXTE,VARCLA,NVARCL,
     * TITRE,SORLEO,NSOR,IKLE,NELEM,NPTFR,NDP,DATE,TIME,
     * NCSIZE,NPTIR,KNOLG,NPLAN,I3,I4)
C
C
C   NOTE JMH 01/12/2003 : VARCLA,NVARCL ARE NO LONGER USED
C
C
C***********************************************************************
C BIEF VERSION 5.6         27/12/05    J-M HERVOUET (LNH) 01 30 71 80 18
C
C***********************************************************************
C
C     FONCTIONS :  ECRITURE DU FICHIER GEOMETRIQUE AU STANDARD SELAFIN
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   X,Y          |<-->| COORDONNEES DU MAILLAGE.
C |   NPOIN        |<-->| NOMBRE DE POINTS DU MAILLAGE.
C |   NBOR         | -->| NUMEROTAION GLOBALE DES POINTS DE BORD.
C |   NFIC         | -->| NUMERO DE CANAL DU FICHIER A LIRE OU ECRIRE.
C |   STAND        | -->| NON UTILISE
C |   STD          | -->| BINAIRE DU FICHIER (STD, IBM, I3E)
C |   NVAR         |<-->| NOMBRE DE VARIABLES DANS LE FICHIER
C |   TEXTE        |<-->| NOMS ET UNITES DES VARIABLES.
C |   VARCLA       | -->| TABLEAU AVEC LES NOMS DES VARIABLES CLANDESTINES.
C |   NVARCL       | -->| NOMBRE DE VARIABLES CLANDESTINES.
C |   TITRE        |<-->| TITRE DU FICHIER.
C |   SORLEO       | -->| VARIABLES QUE L'ON SOUHAITE ECRIRE DANS LE
C |                |    | FICHIER (TABLEAU DE 26 LOGIQUES)
C |   NSOR         | -->| DIMENSION DE SORLEO ET SORIMP
C |   W            | -->| TABLEAU DE TRAVAIL CONSIDERE ICI COMME REEL
C |                |    | DE TAILLE NPOIN.
C |   IKLE         |<-->| TABLE DE CONNECTIVITE (I.E. PASSAGE DE LA
C |                |    | NUMEROTATION LOCALE DES POINTS D'UN ELEMENT
C |                |    | A LA NUMEROTATION GLOBALE
C |   NELEM        |<-->| NOMBRE D'ELEMENTS DU MAILLAGE.
C |   NPTFR        |<-->| NOMBRE DE POINTS FRONTIERE DU DOMAINE.
C |   NDP          |<-->| NOMBRE DE SOMMETS PAR ELEMENT.
C |   DATE,TIME    | -->| DATE (3 INTEGERS) AND TIME (3 INTEGERS)
C |   NCSIZE       | -->| NUMBER OF PROCESSORS
C |   NPTIR        | -->| NUMBER OF INTERFACE POINTS IN PARALLEL
C |   KNOLG        | -->| GLOBAL NUMBERS OF LOCAL POONTS IN PARALLEL
C |   NPLAN        | -->| NUMBER OF PLANES (3D MESHES IN PRISMS)
C |   I3,I4        | -->| INTEGERS, WILL BE PUT IN FILE IN POSITION 3
C |                |    | AND 4 OF THE ARRAY OF 10 INTEGERS
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : LIT , ECRIT
C
C***********************************************************************
C
C    LISTE DES ENREGISTREMENTS DU FICHIER GEOMETRIQUE:
C
C      1    : TITRE DE L'ETUDE
C      2    : NOMBRE DE FONCTIONS LUES SUR LA GRILLE 1 ET LA GRILLE 2.
C      3    : NOM ET UNITE DES VARIABLES
C      4    : 1,0,0,0,0,0,0,0,0,0
C      5    : NELEM,NPOIN,NDP,1
C      6    : IKLE
C      7    : IPOBO TABLEAU DE DIMENSION NPOIN, 0 POUR LES POINTS
C             INTERIEURS, UN NUMERO SINON.
C      8    : X
C      9    : Y
C
C    CE QUI SUIT N'EST PAS FAIT DANS FM3SEL.
C
C     10    : TEMPS
C     11    : VARIABLES DECLAREES EN 3 (DANS L'ORDRE DES DECLARATIONS)
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPOIN,NFIC,NVARCL,NSOR,NELEM,NPTFR,NDP
      INTEGER, INTENT(OUT) :: NVAR
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*)
C                                    IKLE(NELEM,NDP)
      INTEGER, INTENT(IN) :: NBOR(*),IKLE(*)
      CHARACTER(LEN=32), INTENT(IN) :: TEXTE(*),VARCLA(NVARCL)
C                                            NSOR      NSOR+NVARCL
      CHARACTER(LEN=72), INTENT(IN) :: TITRE
      LOGICAL, INTENT(IN) :: SORLEO(*)
      INTEGER, INTENT(IN) :: NCSIZE,NPTIR
      INTEGER, INTENT(IN) :: TIME(3),DATE(3)
      INTEGER, INTENT(IN) :: KNOLG(NPOIN)
      INTEGER, INTENT(IN), OPTIONAL :: NPLAN,I3,I4 
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION XBID(2)
C
      INTEGER IB(10),ISTAT,I,IBID(1),IELEM,ERR
C
      INTEGER, ALLOCATABLE :: IPOBO(:),IKLES(:)
C
      LOGICAL YA_IPOBO,YA_IKLES
C
      CHARACTER*2 CBID
      CHARACTER*80 TITSEL
C
C-----------------------------------------------------------------------
C
      YA_IPOBO = .FALSE.
      YA_IKLES = .FALSE.
C
C   ON SE PLACE AU DEBUT DU FICHIER
C
      REWIND NFIC
C
C   LEC/ECR 1   : NOM DU FICHIER GEOMETRIQUE.
C
      TITSEL = TITRE // 'SERAPHIN'
      CALL ECRI2(XBID,IBID,TITSEL,80,'CH',NFIC,'STD',ISTAT)
C
C   LEC/ECR 2   : NOMBRE DE FONCTIONS DE DISCRETISATION 1 ET 2
C
      IB(1)=0
      IB(2)=0
      DO 91 I=1,NSOR
        IF(SORLEO(I)) IB(1) = IB(1) + 1
91    CONTINUE
      CALL ECRI2(XBID,IB,CBID,2,'I ',NFIC,'STD',ISTAT)
      NVAR =  IB(1)  +  IB(2)
C
C   LEC/ECR 3 : NOMS ET UNITES DES VARIABLES
C
      IF(NVAR.GE.1) THEN
        DO I=1,NSOR
          IF(SORLEO(I)) THEN
           CALL ECRI2(XBID,IBID,TEXTE(I)(1:32),32,'CH',NFIC,'STD',ISTAT)
          ENDIF
        ENDDO
C       IF(NVARCL.NE.0) THEN
C         DO I=1,NVARCL
C         CALL ECRI2(XBID,IBID,VARCLA(I)(1:32),32,'CH',NFIC,'STD',ISTAT)
C         ENDDO
C       ENDIF
      ENDIF
C
C   LEC/ECR 4   : LISTE DE 10 PARAMETRES ENTIERS
C
        IB(1) = 1
        DO 29 I = 2,10
         IB(I) = 0
29      CONTINUE
C
C       ORIGIN COORDINATES IN METRES
C
        IF(PRESENT(I3)) IB(3)=I3
        IF(PRESENT(I4)) IB(4)=I4
C
C       NUMBER OF PLANES IN 3D
C
        IF(PRESENT(NPLAN)) IB(7)=NPLAN
C
CPARA   MARQUAGE POUR ANNONCER LA LECTURE DE KNOLG
        IF(NCSIZE.GT.1) THEN
          IB(8)=NPTFR
          IB(9)=NPTIR
        ENDIF
CPARAFIN
C   Y-A-T-IL PASSAGE DE LA DATE ?
        IF(DATE(1)+DATE(2)+DATE(3)+TIME(1)+TIME(2)+TIME(3).NE.0) THEN
         IB(10) = 1
        ENDIF
C   ECRITURE DU TABLEAU DE 10 PARAMETRES
        CALL ECRI2(XBID,IB,CBID,10,'I ',NFIC,'STD',ISTAT)
C   PASSAGE DE LA DATE
        IF(IB(10).EQ.1) THEN
          IB(1)=DATE(1)
          IB(2)=DATE(2)
          IB(3)=DATE(3)
          IB(4)=TIME(1)
          IB(5)=TIME(2)
          IB(6)=TIME(3)
          CALL ECRI2(XBID,IB,CBID,6,'I ',NFIC,'STD',ISTAT)
        ENDIF
C
C   LEC/ECR 5 : 4 ENTIERS
C
      IF(NDP.NE.4) THEN
        IB(1) = NELEM
      ELSE
C       TETRAHEDRONS REGROUPED INTO PRISMS
        IB(1)=NELEM/3
      ENDIF
      IB(2) = NPOIN
      IF(NDP.NE.4) THEN
        IB(3) = NDP
      ELSE
C       TETRAHEDRONS REGROUPED INTO PRISMS
        IB(3) = 6
      ENDIF
      IB(4) = 1
      CALL ECRI2(XBID,IB,CBID,4,'I ',NFIC,'STD',ISTAT)
C
C   LEC/ECR 6 : IKLE
C
      IF(NDP.NE.4) THEN
        ALLOCATE(IKLES(NELEM*NDP),STAT=ERR)
      ELSE
C       TETRAHEDRONS REGROUPED INTO PRISMS
        ALLOCATE(IKLES(NELEM*2)  ,STAT=ERR)
      ENDIF
      IF(ERR.NE.0) STOP 'ECRGEO : ALLOCATION DE IKLES'
      YA_IKLES = .TRUE.
C     INVERSION DE IKLE  EN IKLES POUR SELAFIN
      IF(NDP.NE.4) THEN
        DO I      = 1,NDP
          DO IELEM  = 1,NELEM
            IKLES((IELEM-1)*NDP+I) = IKLE((I-1)*NELEM+IELEM)
          ENDDO
        ENDDO
      ELSE
C     TETRAHEDRONS REGROUPED INTO PRISMS
        DO IELEM  = 1,NELEM/3
          IKLES((IELEM-1)*6+1) = IKLE(      IELEM)
          IKLES((IELEM-1)*6+2) = IKLE(NELEM+IELEM)
          IKLES((IELEM-1)*6+3) = IKLE(NELEM+IELEM)
          IKLES((IELEM-1)*6+4) = IKLE(      IELEM)+NPOIN/NPLAN
          IKLES((IELEM-1)*6+5) = IKLE(NELEM+IELEM)+NPOIN/NPLAN
          IKLES((IELEM-1)*6+6) = IKLE(NELEM+IELEM)+NPOIN/NPLAN
        ENDDO
      ENDIF
C  
      IF(NDP.NE.4) THEN
      CALL ECRI2(XBID,IKLES,CBID,NELEM*NDP,'I ',NFIC,'STD',ISTAT)
      ELSE
C     TETRAHEDRONS REGROUPED INTO PRISMS
      CALL ECRI2(XBID,IKLES,CBID,NELEM*2,'I ',NFIC,'STD',ISTAT)
      ENDIF
C
C   LEC/ECR 7 : IPOBO (CAS DES FICHIERS SANS PARALLELISME)
C
      IF(IB(8).EQ.0.AND.IB(9).EQ.0) THEN
C
        ALLOCATE(IPOBO(NPOIN),STAT=ERR)
        IF(ERR.NE.0) STOP 'ECRGEO : ALLOCATION DE IPOBO'
        YA_IPOBO = .TRUE.
        DO 40 I=1,NPOIN
         IPOBO(I) = 0
40      CONTINUE
C       ONLY LATERAL BOUNDARY POINTS WITH PRISMS
        DO 41 I =1,NPTFR
         IPOBO(NBOR(I)) = I
41      CONTINUE
        CALL ECRI2(XBID,IPOBO,CBID,NPOIN,'I ',NFIC,'STD',ISTAT)
C
      ENDIF
C
      IF(IB(8).NE.0.OR.IB(9).NE.0) THEN
C
C   LEC/ECR  7.1 KNOLG (SEULEMENT EN CAS DE PARALLELISME)
C
      CALL ECRI2(XBID,KNOLG,CBID,NPOIN,'I ',NFIC,'STD',ISTAT)
C
      ENDIF
C
C   LEC/ECR  8 ET 9 : X ET Y  COORDONNEES DES POINTS DU MAILLAGE
C
      CALL ECRI2(X   ,IBID,CBID,NPOIN,'R4',NFIC,'STD',ISTAT)
      CALL ECRI2(Y   ,IBID,CBID,NPOIN,'R4',NFIC,'STD',ISTAT)
C
C-----------------------------------------------------------------------
C
      IF(YA_IPOBO) DEALLOCATE(IPOBO)
      IF(YA_IKLES) DEALLOCATE(IKLES)
C
C-----------------------------------------------------------------------
C
      RETURN
      END

