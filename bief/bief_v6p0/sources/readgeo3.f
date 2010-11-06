C                       *******************
                        SUBROUTINE READGEO3
C                       *******************
C
     *(KNOLG,X,Y,NPOIN,NFIC,IB,Z)
C
C***********************************************************************
C BIEF VERSION 5.6          19/10/03  J-M HERVOUET (LNHE) 01 30 71 80 18
C                                     REGINA NEBAUER
C                                     LAM MINH PHUONG
C                           19/10/05  EMILE RAZAFINDRAKOTO
C***********************************************************************
C
C   USER SUBROUTINE READGEO3
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
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : LIT 
C
C***********************************************************************
C
      USE BIEF, EX_READGEO3 => READGEO3
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C     
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)           :: NPOIN,NFIC
      INTEGER, INTENT(INOUT)        :: IB(10)
      INTEGER, INTENT(OUT)          :: KNOLG(NPOIN)
      DOUBLE PRECISION, INTENT(OUT) :: X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: Z(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION XB(2)
      REAL, ALLOCATABLE :: RB(:)
      REAL RBID(1)
      INTEGER ISTAT,ERR
      CHARACTER(LEN=1)  :: CB
C
C-----------------------------------------------------------------------
C
C     LE DEBUT DU FICHIER EST DEJA LU PAR READGEO1
C
C     REWIND NFIC
C
C     7 : IPOBO REMPLACE PAR KNOLG (CAS DES FICHIERS AVEC PARALLELISME)
C
      IF(IB(8).NE.0.OR.IB(9).NE.0) THEN
C       PARALLELISME
C       CAS OU IPOBO EST REMPLACE PAR KNOLG 
        CALL LIT(XB,RBID,KNOLG,CB,NPOIN,'I ',NFIC,'STD',ISTAT)
      ENDIF
C
C     8 ET 9 : COORDONNEES X ET Y
C
      ALLOCATE(RB(NPOIN),STAT=ERR)
      IF(ERR.NE.0) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'READGEO3 : ALLOCATION DE RB DEFECTUEUSE'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'READGEO3 : WRONG ALLOCATION OF RB'
        ENDIF
        STOP
      ENDIF
C     
      CALL LIT(X   ,RB,IB,CB,NPOIN,'R4',NFIC,'STD',ISTAT)
      CALL LIT(Y   ,RB,IB,CB,NPOIN,'R4',NFIC,'STD',ISTAT)
C
C     SPECIAL FORMAT FOR TETRAHEDRONS : Z AFTER X AND Y
C     A RECORD FOR TIME IS PRESENT WITH THE SELAFIN FORMAT
C     WHEN Z IS GIVEN AS VARIABLE IN TIME, BUT THIS IS NEVER USED.
C
      IF(PRESENT(Z)) THEN
C       RECORD FOR TIME
C       CALL LIT(Z,RB,IB,CB,1,'R4',NFIC,'STD',ISTAT)
C       RECORD FOR Z (FIRST VARIABLE IN SELAFIN FORMAT)
        CALL LIT(Z,RB,IB,CB,NPOIN,'R4',NFIC,'STD',ISTAT)
      ENDIF
C
C-----------------------------------------------------------------------
C
CMODIF ER 19/10/2005
      CALL CORRXY(X,Y,NPOIN)
CFIN MODIF ER 19/10/2005
C
      DEALLOCATE(RB)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
