C                       **********************
                        SUBROUTINE FIND_IN_SEL
C                       **********************
C
     *(RES,NAME,NFIC,W,OK,RECORD,NP,TIME)
C
C***********************************************************************
C  BIEF VERSION 5.2           08/08/98    J-M HERVOUET (LNH) 30 71 80 18
C                                          
C***********************************************************************
C
C  FONCTION  :  RECHERCHE D'UN TABLEAU DE RESULTATS DANS UN FICHIER
C               AU FORMAT SELAFIN.
C
C
C
C  ATTENTION : IT SEEMS THAT THERE IS A BUG IN COMPILER NAG WITH
C              OPTIMIZATION O4 : TIME IS THEN MANDATORY.
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   RES          |<-- | WHERE TO PUT THE RESULT
C |   NAME         | -->| NAME OF VARIABLE (16 CHARACTERS)
C |   NFIC         | -->| NUMERO DU CANAL DU FICHIER
C |   W            |    | TABLEAU DE TRAVAIL REEL DE DIMENSION NPOIN.
C |   OK           |<-- | TRUE IF ARRAY IS FOUND
C |   
C |   PARAMETRES OPTIONNELS:  
C |  
C |   RECORD       | -->| NUMBER OF THE REQUESTED RECORD
C |   NP           |<-- | NUMBER OF POINTS
C |   TIME         |<-- | TIME OF THE RECORD
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : PREDON
C PROGRAMMES APPELES : LIT , OV
C
C***********************************************************************
C
      USE BIEF, EX_FIND_IN_SEL => FIND_IN_SEL
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: RES
      CHARACTER(LEN=16), INTENT(IN) :: NAME
      LOGICAL, INTENT(OUT)          :: OK
      REAL, INTENT(INOUT)           :: W(*)
      INTEGER, INTENT(IN) :: NFIC
      INTEGER, INTENT(IN),  OPTIONAL          :: RECORD
      INTEGER, INTENT(OUT), OPTIONAL          :: NP
      DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: TIME
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NPOIN,ISTAT,I,NVAR,IB(10),REC,IREC
C
      DOUBLE PRECISION XB(2)
      REAL RB(2)
C
      CHARACTER*1 CB
      CHARACTER*32 TEXTLU(36)
C
C-----------------------------------------------------------------------
C
      IF(PRESENT(RECORD)) THEN
        REC = RECORD
      ELSE
        REC = 1
      ENDIF
C
      OK = .FALSE.
C
C-----------------------------------------------------------------------
C
C     LECTURE 'RAPIDE' JUSQU'A L'ENREGISTREMENT DU TEMPS
C
C
C     ON SE PLACE AU DEBUT DU FICHIER
C
      REWIND NFIC
C
C     1: TITLE.
      CALL LIT(XB,RB,IB,CB,1,'CH',NFIC,'STD',ISTAT)
C
C     2: NUMBER OF ARRAYS IN THE RESULT FILE
      CALL LIT(XB,RB,IB,CB,2,'I ',NFIC,'STD',ISTAT)
      NVAR =  IB(1)  +  IB(2)
C
C     3: NAMES AND UNITS OF VARIABLES
      IF(NVAR.GE.1) THEN
        DO I=1,NVAR
           CALL LIT(XB,RB,IB,TEXTLU(I),32,'CH',NFIC,'STD',ISTAT)
        ENDDO
      ENDIF
C
C     4: LIST OF 10 INTEGER PARAMETERS
      CALL LIT(XB,RB,IB,CB,10,'I ',NFIC,'STD',ISTAT)
C     CASE WHERE DATE AND TIME IN THE FILE
      IF(IB(10).EQ.1) CALL LIT(XB,RB,IB,CB,6,'I ',NFIC,'STD',ISTAT)
C
C     5: 4 INTEGERS
      CALL LIT(XB,RB,IB,CB,4,'I ',NFIC,'STD',ISTAT)
      NPOIN = IB(2)
C
      IF(PRESENT(NP)) NP = NPOIN
C
C     6: IKLES (LIKE IKLE BUT INDICES EXCHANGED)
      CALL LIT(XB,RB,IB,CB,1,'I ',NFIC,'STD',ISTAT)
C
C     7: IPOBO OU KNOLG
      CALL LIT(XB,RB,IB,CB,1,'I ',NFIC,'STD',ISTAT)
C
C     8 ET 9: X ET Y
      CALL LIT(XB,W,IB,CB,1,'R4',NFIC,'STD',ISTAT)
      CALL LIT(XB,W,IB,CB,1,'R4',NFIC,'STD',ISTAT)
C
C-----------------------------------------------------------------------
C
      IREC = 0
500   IREC = IREC + 1
      IF (NVAR.GE.1) THEN
C
C       ENREGISTREMENT DU TEMPS
C
        CALL LIT(XB,W,IB,CB,1,'R4',NFIC,'STD',ISTAT)
C       NOTE JMH : THE FOLLOWING INSTRUCTION RAISES PROBLEMS
C       WHEN TIME IS NOT PRESENT, WITH NAG COMPILER WITH OPTION -O4
        IF(PRESENT(TIME)) TIME=XB(1)
C
        DO I=1,NVAR
C
C         LECTURE DE LA VARIABLE, OU SAUT DE L'ENREGISTREMENT
          IF(TEXTLU(I)(1:16).EQ.NAME.AND.REC.EQ.IREC) THEN
            CALL LIT(RES%R,W,IB,CB,NPOIN,'R4',NFIC,'STD',ISTAT)
            OK=.TRUE.
          ELSE
            CALL LIT(XB,W,IB,CB,1,'R4',NFIC,'STD',ISTAT)
          ENDIF
C
        ENDDO
C
      ENDIF
      IF(IREC.NE.REC) GO TO 500
C
C-----------------------------------------------------------------------
C
      RETURN
      END
