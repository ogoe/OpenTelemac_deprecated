C                       ***********************
                        INTEGER FUNCTION DIMENS
C                       ***********************
C
     *( IELM )
C
C***********************************************************************
C BIEF VERSION 5.9      13/02/2008    J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C  FONCTION  : DONNE LA DIMENSION D'UN ELEMENT
C
C  NOTE : ON POURRAIT FAIRE UN TABLEAU EN DATA POUR ACCELERER
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   IELM         | -->| TYPE D'ELEMENT.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: IELM
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      IF(IELM.EQ.0 .OR.
     *   IELM.EQ.1 .OR. 
     *   IELM.EQ.2) THEN
C
        DIMENS = 1
C
      ELSEIF(IELM.EQ.10.OR.
     *       IELM.EQ.11.OR.
     *       IELM.EQ.12.OR.
     *       IELM.EQ.13.OR.
     *       IELM.EQ.14.OR.
     *       IELM.EQ.70.OR.
     *       IELM.EQ.71.OR.
     *       IELM.EQ.80.OR.
     *       IELM.EQ.81.OR.     
     *       IELM.EQ.61.OR.
     *       IELM.EQ.60.OR.
     *       IELM.EQ.20.OR.
     *       IELM.EQ.21) THEN
C
        DIMENS = 2
C
      ELSEIF(IELM.EQ.30.OR.
     *       IELM.EQ.31.OR.
     *       IELM.EQ.40.OR.
     *       IELM.EQ.41.OR.
     *       IELM.EQ.50.OR.
     *       IELM.EQ.51    ) THEN
C
        DIMENS = 3
C
      ELSE
        IF(LNG.EQ.1) WRITE(LU,100) IELM
        IF(LNG.EQ.2) WRITE(LU,101) IELM
100     FORMAT(1X,'DIMENS (BIEF) : ',1I6,' ELEMENT NON PREVU')
101     FORMAT(1X,'DIMENS (BIEF) : ',1I6,' ELEMENT NOT IMPLEMENTED')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
