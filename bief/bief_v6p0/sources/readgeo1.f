C                       *******************
                        SUBROUTINE READGEO1
C                       *******************
C
     *(NPOIN,NELEM,NPTFR,NDP,IB,NFIC,NELEBD)
C
C***********************************************************************
C BIEF VERSION 5.5           29/04/04  J-M HERVOUET (LNH) 01 30 71 80 18
C                                      REGINA NEBAUER
C                                      LAM MINH PHUONG
C***********************************************************************
C
C   USER SUBROUTINE READGEO1
C
C   FUNCTION: 
C     
C   READS OR COMPUTES THE VALUES OF NPOIN, NELEM, NPTFR, MXPTVS, MXELVS
C   IN THE GEOMETRY FILE IN THE CHANNEL NGEO.
C
C   MAY BE REWRITTEN FOR ANOTHER FILE FORMAT
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   NPOIN        |<-- | NOMBRE DE POINTS DU MAILLAGE.
C |   NELEM        |<-- | NOMBRE D'ELEMENTS DU MAILLAGE.
C |   NPTFR        |<-- | NOMBRE DE POINTS FRONTIERE DU DOMAINE.
C |   NDP          |<-- | NOMBRE DE NOEUD PAR ELEMENT
C |   IB           |<-- | 10 INTEGERS, SEE SELAFIN FILE STANDARD
C |   NFIC         | -->| UNITE LOGIQUE FICHIER GEO
C |   NELEBD       |<-- | NOMBRE D'ELEMENTS DE BORD 
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : LIT 
C
C***********************************************************************
C
      USE BIEF, EX_READGEO1 => READGEO1
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C     
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(OUT)         :: NPOIN  ! Number of mesh nodes
      INTEGER, INTENT(OUT)         :: NELEM  ! Number of elements
      INTEGER, INTENT(OUT)         :: NDP    ! Number of element faces 
      INTEGER, INTENT(OUT)         :: IB(10) ! Integer array
      INTEGER, INTENT(OUT)         :: NPTFR  ! Number of border nodes
      INTEGER, INTENT(IN)          :: NFIC   ! File to read
      INTEGER,OPTIONAL,INTENT(OUT) :: NELEBD ! Number of border elements
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION  :: XB(2)
      REAL              :: RB(2)
      INTEGER           :: ISTAT
      INTEGER           :: NVAR
      INTEGER           :: I,IB6(6)
      CHARACTER(LEN=2)  :: CB
      CHARACTER(LEN=72) :: TITRE
C
C-----------------------------------------------------------------------
C
C   ON SE PLACE AU DEBUT DU FICHIER
C
      REWIND NFIC
C
C     1: TITLE.
      CALL LIT(XB,RB,IB,TITRE,72,'CH',NFIC,'STD',ISTAT)      
C
C     2: NUMBER OF ARRAYS IN THE RESULT FILE
      CALL LIT(XB,RB,IB,CB,2,'I ',NFIC,'STD',ISTAT)
      NVAR =  IB(1)  +  IB(2)      
C     3: NAMES AND UNITS OF VARIABLES
      IF(NVAR.GE.1) THEN
        DO I=1,NVAR
          CALL LIT(XB,RB,IB,CB,2,'CH',NFIC,'STD',ISTAT)                
        ENDDO
      ENDIF
C
C     4: LIST OF 10 INTEGER PARAMETERS
      CALL LIT(XB,RB,IB,CB,10,'I ',NFIC,'STD',ISTAT)
C      
C     CASE WHERE DATE AND TIME IN THE FILE
      IF(IB(10).EQ.1) CALL LIT(XB,RB,IB6,CB,6,'I ',NFIC,'STD',ISTAT)     
C
C     READ THE NUMBER OF BORDER ELEMENTS FOR 3D MESH
      IF(IB(7).NE.0.AND.PRESENT(NELEBD)) THEN
        NELEBD = IB(7)
      END IF     
C     CASE WHERE KNOLG IS GIVEN INSTEAD IPOBO (PARALLELISM)
      IF(IB(8).NE.0) THEN 
        NPTFR=IB(8)
C       NOTE JMH : NEXT LINE MOVED AFTER ENDIF ON 22/07/02
C                  SUBDOMAINS MAY HAVE NPTFR=0
C       NPTIR=IB(9)
      ENDIF
      NPTIR=IB(9)
C
C     5: 4 INTEGERS
      CALL LIT(XB,RB,IB6,CB,4,'I ',NFIC,'STD',ISTAT)  
C
      NELEM = IB6(1)
      NPOIN = IB6(2)
      NDP   = IB6(3)
C
C-----------------------------------------------------------------------
C
C  FORMATS D'IMPRESSION :
C
      IF(LNG.EQ.1) WRITE(LU,300) TITRE
      IF(LNG.EQ.1) WRITE(LU,500) NELEM,NPOIN
      IF(LNG.EQ.2) WRITE(LU,301) TITRE
      IF(LNG.EQ.2) WRITE(LU,501) NELEM,NPOIN
C
      IF(NPOIN.LT.3) THEN
        IF(LNG.EQ.1) WRITE(LU,23) NPOIN
        IF(LNG.EQ.2) WRITE(LU,24) NPOIN
        CALL PLANTE(1)
        STOP
      ENDIF
C
23    FORMAT(1X,'READGEO1 : NOMBRE DE POINTS DU MAILLAGE : ',1I6,/,1X,
     *          '           NOMBRE DE POINTS DE FRONTIERE: ',1I6,/,1X,
     *          '           DONNEES ERRONEES, ARRET DU PROGRAMME')
24    FORMAT(1X,'READGEO1 : NUMBER OF POINTS IN THE MESH: ',1I6,/,1X,
     *          '           NUMBER OF BOUNDARY POINTS: ',1I6,/,1X,
     *          '           WRONG DATA, PROGRAMME STOPPED')
300   FORMAT(1X,//,1X,'READGEO1 : TITRE= ',A72,/)
301   FORMAT(1X,//,1X,'READGEO1: TITLE= ',A72,/)
500   FORMAT(1X,'NOMBRE D''ELEMENTS:',1I6,/,
     *       1X,'NOMBRE REEL DE POINTS:',1I6)
501   FORMAT(1X,'NUMBER OF ELEMENTS:',1I6,/,
     *       1X,'NUMBER OF POINTS:',1I6)
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
