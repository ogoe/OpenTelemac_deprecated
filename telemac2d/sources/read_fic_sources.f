C                       ***************************
                        SUBROUTINE READ_FIC_SOURCES
C                       ***************************
C
     *( Q , WHAT , AT , NFIC , LISTIN , STAT )
C
C***********************************************************************
C  TELEMAC 2D VERSION 6.0  28/06/2010  J-M HERVOUET (LNH) 01 30 87 80 18
C                            
C***********************************************************************
C
C IMPORTANT NOTE : THIS SUBROUTINE IS A COPY OF 
C                  SUBROUTINE READ_FIC_FRLIQ BECAUSE IT USES THE SAME
C                  FILE FORMAT (ONLY MESSAGES ON LISTING ARE CHANGED). 
C                  THE ONLY DIFFERENCE IS THAT 
C                  THE ALLOCATABLE ARRAYS TIME AND INFIC WILL STORE
C                  DIFFERENT DATA.
C
C                  THE PROBLEM IS : WHERE TO STORE THESE DATA BECAUSE
C                  THESE ROUTINES MAY BE CALLED BY TELEMAC-2D OR 3D
C                  A SPECIFIC MODULE COULD BE DONE
C
C
C 28/06/2010 : SIZE OF LINE PARAMETERIZED (SEE SIZELIGN)
C
C-----------------------------------------------------------------------                  
C
C FUNCTION : READ AND INTERPOLATE VALUES IN THE SOURCES FILE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   WHAT         | -->| VARIABLE TO LOOK FOR IN 8 CHARACTERS
C |     AT         | -->| TIME
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : BORD
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER*8     , INTENT(IN)       :: WHAT
      DOUBLE PRECISION, INTENT(IN)       :: AT
      DOUBLE PRECISION, INTENT(INOUT)    :: Q
      INTEGER         , INTENT(IN)       :: NFIC
      LOGICAL         , INTENT(IN)       :: LISTIN
      LOGICAL         , INTENT(OUT)      :: STAT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      LOGICAL DEJA
      DATA DEJA /.FALSE./
C
C     MAXIMUM NUMBER OF CHARACTERS PER LIGN (MAY BE CHANGED)
C
      INTEGER, PARAMETER :: SIZELIGN = 300
C
      INTEGER IVALUE,NVALUE,ILIG,NLIG,OK,J,IWHAT,IDEB,IFIN,IL1,IL2
      INTEGER, PARAMETER :: MAXVAL=50
      DOUBLE PRECISION TL1,TL2,TETA,TOL,LASTAT
C
      CHARACTER(LEN=SIZELIGN) LIGNE
      CHARACTER*8 CHOIX(MAXVAL),LASTWHAT
C      
      DATA TOL /1.D-3/
C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: INFIC
      DOUBLE PRECISION, DIMENSION(:)  , ALLOCATABLE :: TIME
C
      SAVE DEJA,INFIC,TIME,CHOIX,IL1,IL2,TL1,TL2,NVALUE,LASTWHAT,LASTAT
      SAVE NLIG
C
      INTRINSIC ABS
C
C-----------------------------------------------------------------------
C
C     1) FIRST PART (AT FIRST CALL)
C        READING THE FILE OF LIQUID BOUNDARIES
C        INITIALISING CURRENT LINES AND INTERVAL OF TIME
C
      IF(.NOT.DEJA) THEN
        REWIND(NFIC)
C       SKIPPING COMMENTS
1       READ(NFIC,FMT='(A)') LIGNE
        IF(LIGNE(1:1).EQ.'#') GO TO 1
C
C       FINDING WHAT AND HOW MANY VALUES ARE GIVEN IN THE FILE
C
        NVALUE = -1
        IFIN = 1
40      IDEB = IFIN                
C
C       LOOKING FOR FIRST CHARACTER OF NAME
50      IF(LIGNE(IDEB:IDEB).EQ.' '.AND.IDEB.LT.SIZELIGN) THEN
          IDEB=IDEB+1
          GO TO 50
        ENDIF
C       LOOKING FOR LAST CHARACTER OF NAME
        IFIN = IDEB
60      IF(LIGNE(IFIN:IFIN).NE.' '.AND.IFIN.LT.SIZELIGN) THEN
          IFIN=IFIN+1
          GO TO 60
        ENDIF
C
        IF(IDEB.EQ.IFIN) GO TO 4
C
        NVALUE = NVALUE + 1
        IF(NVALUE.EQ.0) THEN
          IF(LIGNE(IDEB:IFIN-1).NE.'T') THEN
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'LA PREMIERE VARIABLE DOIT ETRE LE TEMPS T'
            WRITE(LU,*) 'DANS LE FICHIER DES SOURCES'
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'FIRST VALUE MUST BE TIME, DENOTED T'
            WRITE(LU,*) 'IN SOURCES FILE'
          ENDIF
          CALL PLANTE(1)
          STOP
          ENDIF
        ELSEIF(NVALUE.LE.MAXVAL) THEN
          CHOIX(NVALUE)='        '
          CHOIX(NVALUE)(1:IFIN-IDEB+1)=LIGNE(IDEB:IFIN-1)
        ELSE
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'AUGMENTER MAXVAL DANS READ_FIC_SOURCES'
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'INCREASE MAXVAL IN READ_FIC_SOURCES'
          ENDIF
          CALL PLANTE(1)
          STOP          
        ENDIF
        IF(IFIN.LT.SIZELIGN) GO TO 40
C
C       SKIPPING THE LINE WITH UNITS OR NAMES
4       READ(NFIC,FMT='(A)') LIGNE
        IF(LIGNE(1:1).EQ.'#') GO TO 4
C
C       COUNTING LINES OF DATA
        NLIG = 0
998     READ(NFIC,*,END=1000,ERR=999) LIGNE
        IF(LIGNE(1:1).NE.'#') NLIG=NLIG+1
        GO TO 998
999     CONTINUE
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'ERREUR DE LECTURE SUR LE FICHIER DES SOURCES'
          WRITE(LU,*) 'LIQUIDES A LA LIGNE DE DONNEES : ',NLIG
          WRITE(LU,*) '(COMMENTAIRES NON COMPTES)'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'READING ERROR ON THE SOURCES FILE'
          WRITE(LU,*) 'AT LINE OF DATA : ',NLIG
          WRITE(LU,*) '(COMMENTS EXCLUDED)'
        ENDIF
        CALL PLANTE(1)
        STOP
1000    CONTINUE    
C
C       DYNAMIC ALLOCATION OF TIME AND INFIC
C
        ALLOCATE(TIME(NLIG),STAT=OK)
        IF(OK.NE.0) WRITE(LU,*) 'MEMORY ALLOCATION ERROR FOR TIME'
        ALLOCATE(INFIC(NVALUE,NLIG),STAT=OK)
        IF(OK.NE.0) WRITE(LU,*) 'MEMORY ALLOCATION ERROR FOR INFIC'
C
C       FINAL READING OF TIME AND INFIC
C
        REWIND(NFIC)
C       SKIPPING COMMENTS AND FIRST TWO MANDATORY LINES
2       READ(NFIC,FMT='(A)') LIGNE
        IF(LIGNE(1:1).EQ.'#') GO TO 2
        READ(NFIC,FMT='(A)') LIGNE
C        
        DO ILIG=1,NLIG       
3         READ(NFIC,FMT='(A)') LIGNE
          IF(LIGNE(1:1).EQ.'#') THEN
            GO TO 3
          ELSE
            BACKSPACE(NFIC)
            READ(NFIC,*) TIME(ILIG),(INFIC(IVALUE,ILIG),IVALUE=1,NVALUE)
          ENDIF
        ENDDO 
C
        CLOSE(NFIC)
        DEJA = .TRUE.
C
        IL1 = 1
        IL2 = 2
        TL1 = TIME(1)
        TL2 = TIME(2) 
C             
      ENDIF
C
C-----------------------------------------------------------------------
C
C     2) INTERPOLATING THE DATA TO GET THE CORRECT TIME
C 
C     2.a) FINDING THE ADDRESS IN THE ARRAY OF STORED DATA     
C
C     2.b) INTERPOLATING DATA OF THE ARRAY INFIC           
C
C-----------------------------------------------------------------------
C
C
C     WHAT VARIABLE IS ASKED ?
      IWHAT = 0
      DO J=1,NVALUE
        IF(WHAT.EQ.CHOIX(J)) IWHAT=J
      ENDDO
      IF(IWHAT.EQ.0) THEN
        STAT=.FALSE.
        RETURN
      ENDIF
C
70    IF(AT.GE.TL1-TOL.AND.AT.LE.TL2+TOL) THEN
        TETA = (AT-TL1)/(TL2-TL1)
      ELSE
         DO J=1,NLIG-1
            IF(AT.GE.TIME(J)-TOL.AND.AT.LE.TIME(J+1)+TOL) THEN
               TL1=TIME(J)
               TL2=TIME(J+1)
               IL1=J
               IL2=J+1
               GO TO 70
            ENDIF
         ENDDO
         IL1=IL2
         IL2=IL2+1
        IF(IL2.GT.NLIG) THEN
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'T=',AT,' HORS LIMITES'
            WRITE(LU,*) 'DU FICHIER DES SOURCES'
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'T=',AT,' OUT OF RANGE'
            WRITE(LU,*) 'OF THE SOURCES FILE'
          ENDIF
          CALL PLANTE(1)
          STOP        
        ENDIF
        TL1=TIME(IL1)
        TL2=TIME(IL2)
        GO TO 70
      ENDIF
C
      Q = (1.D0-TETA)*INFIC(IWHAT,IL1)
     *  +       TETA *INFIC(IWHAT,IL2)
C
      STAT=.TRUE.
C
C     PRINT ONLY IF NEW TIME OR NEW VALUE IS ASKED
C
      IF(LISTIN) THEN
        IF(ABS(AT-LASTAT).GT.TOL.OR.LASTWHAT.NE.WHAT) THEN
          WRITE(LU,*) 'SOURCES : ',WHAT,'=',Q
        ENDIF
      ENDIF
      LASTAT=AT
      LASTWHAT=WHAT
C
C-----------------------------------------------------------------------
C
      RETURN
      END
