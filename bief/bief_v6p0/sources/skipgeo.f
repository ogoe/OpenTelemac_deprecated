C                       ******************
                        SUBROUTINE SKIPGEO
C                       ******************
C
     *(NFIC,TITFIC,NPOIN,NVAR,TEXTLU,NPLAN)
C
C***********************************************************************
C BIEF VERSION 5.5         18/11/04    J-M HERVOUET (LNH) 01 30 71 80 18
C
C***********************************************************************
C
C     FONCTION : SAUTE LA GEOMETRIE DANS UN FICHIER AU FORMAT SELAFIN
C
C                SKIPS THE GEOMETRY IN A SELAFIN FILE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   NFIC         | -->| LOGICAL UNIT OF FILE TO READ
C |   TITFIC       |<-- | TITLE OF FILE (FIRST RECORD)
C |   NPOIN        |<-- | NUMBER OF POINTS IN THE MESH
C |   NVAR         |<-- | NUMBER OF VARIABLES IN THE FILE
C |   TEXTLU       |<-- | NAMES OF VARIABLES (32 CHARACTERS FOR EACH)
C |                |<-- | 16 FIRST : NAME  16 LAST : UNIT
C |   NPLAN        |<-- | NUMBER OF PLANES (3D CAS)  OPTIONAL !!!!!!!
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : LIT 
C
C***********************************************************************
C
C    LISTE DES ENREGISTREMENTS DU FICHIER GEOMETRIQUE:
C
C      1    : TITRE DE L'ETUDE
C      2    : NOMBRE DE FONCTIONS LUES SUR LA GRILLE 1 ET LA GRILLE 2.
C      3    : NOM ET UNITE DES VARIABLES
C      4    : 1,0,0,0,0,0,0,0,0,N
C      4.1  : DATE(3 INTEGERS) AND TIME(3 INTEGERS) IF N=1
C      5    : NELEM,NPOIN,NDP,1
C      6    : IKLE
C      7    : IPOBO TABLEAU DE DIMENSION NPOIN, 0 POUR LES POINTS
C             INTERIEURS, UN NUMERO SINON.
C      8    : X
C      9    : Y
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)            :: NFIC
      INTEGER, INTENT(OUT), OPTIONAL :: NPLAN
      INTEGER, INTENT(OUT)           :: NPOIN,NVAR
      CHARACTER(LEN=72), INTENT(OUT) :: TITFIC
      CHARACTER(LEN=32), INTENT(OUT) :: TEXTLU(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION XBID(1)
      REAL W(1)
      INTEGER IB(10),ISTAT,I,IBID(1)
      CHARACTER*1 CBID
C
C-----------------------------------------------------------------------
C
C   ON SE PLACE AU DEBUT DU FICHIER
C
      REWIND NFIC
C
C   LEC/ECR 1   : NOM DU FICHIER GEOMETRIQUE.
C
      CALL LIT(XBID,W,IBID,TITFIC,72,'CH',NFIC,'STD',ISTAT)
C
C   LEC/ECR 2   : NOMBRE DE FONCTIONS DE DISCRETISATION 1 ET 2
C
      CALL LIT(XBID,W,IB,CBID,2,'I ',NFIC,'STD',ISTAT)
      NVAR = IB(1)+IB(2)
C
C   LEC/ECR 3 : NOMS ET UNITES DES VARIABLES
C
      IF(NVAR.GE.1) THEN
        DO 10 I=1,NVAR
          CALL LIT(XBID,W,IBID,TEXTLU(I),32,'CH',NFIC,'STD',ISTAT)
10      CONTINUE
      ENDIF
C
C   LEC/ECR 4   : LISTE DE 10 PARAMETRES ENTIERS
C
      CALL LIT(XBID,W,IB,CBID,10,'I ',NFIC,'STD',ISTAT)
      IF(PRESENT(NPLAN)) NPLAN=IB(7)
      IF(IB(10).EQ.1) THEN
        CALL LIT(XBID,W,IB,CBID,1,'I ',NFIC,'STD',ISTAT)
      ENDIF
C
C   LEC/ECR 5 : 4 ENTIERS
C
      CALL LIT(XBID,W,IB,CBID,4,'I ',NFIC,'STD',ISTAT)
      NPOIN = IB(2)
C
C   LEC/ECR 6 : IKLE
C
      CALL LIT(XBID,W,IB,CBID,1,'I ',NFIC,'STD',ISTAT)
C
C   LEC/ECR 7 : IPOBO (CAS DES FICHIERS SANS PARALLELISME)
C
      CALL LIT(XBID,W,IB,CBID,1,'I ',NFIC,'STD',ISTAT)
C
C   LEC/ECR  8 ET 9 : X ET Y  COORDONNEES DES POINTS DU MAILLAGE
C
      CALL LIT(XBID,W,IBID,CBID,1,'R4',NFIC,'STD',ISTAT)
      CALL LIT(XBID,W,IBID,CBID,1,'R4',NFIC,'STD',ISTAT)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
