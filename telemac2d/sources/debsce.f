C                       ********************************
                        DOUBLE PRECISION FUNCTION DEBSCE
C                       ********************************
C
     *( TIME , I , DISCE )
C
C***********************************************************************
C  TELEMAC 2D VERSION 6.0   03/04/08  J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C FONCTION  : DONNE LA VALEUR DU DEBIT POUR TOUTES LES SOURCES
C
C             PERMET DE PROGRAMMER DES VARIATIONS EN FONCTION DU TEMPS
C             ET DE LA PROFONDEUR.
C
C-----------------------------------------------------------------------
C
C FUNCTION: GIVES THE PRESCRIBED DISCHARGE OF EVERY POINT SOURCE
C           VARIATIONS FUNCTION OF TIME AND SPACE MAY BE IMPLEMENTED
C
C-----------------------------------------------------------------------
C
C NOTE: T2DVEF IS THE SOURCES FILE IN TELEMAC-2D
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   TIME         | -->| TIME
C |   I            | -->| NUMBER OF THE SOURCE
C |   DISCE        | -->| ARRAY OF DISCHARGES OF SOURCES.
C |                |    | READ IN THE PARAMETER FILE.
C |                |    | NAME OF DISCE IS DSCE IN TELEMAC-2D.
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : DIFSOU ET PROSOU
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC2D, ONLY: MAXSCE,AT,ENTET,NREJET,
     *                                  T2D_FILES,T2DVEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(IN) :: TIME,DISCE(*)
      INTEGER         , INTENT(IN) :: I
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER*8 FCT
      INTEGER N
      LOGICAL DEJA,OK(MAXSCE)
      DATA    DEJA /.FALSE./
      SAVE    OK,DEJA
C
C     FIRST CALL, OK INITIALISED TO .TRUE.
C
      IF(.NOT.DEJA) THEN
        DO N=1,NREJET
          OK(N)=.TRUE.
        ENDDO
        DEJA=.TRUE.
      ENDIF
C
C     IF SOURCES FILE EXISTING, ATTEMPT TO FIND
C     THE VALUE IN IT. IF YES, OK REMAINS TO .TRUE. FOR NEXT CALLS
C                      IF  NO, OK SET     TO .FALSE.
C
      IF(OK(I).AND.T2D_FILES(T2DVEF)%NAME(1:1).NE.' ') THEN
C
C       FCT WILL BE Q(1), Q(2), ETC, Q(99), DEPENDING ON I
        FCT(1:2)='Q('
        IF(I.LT.10) THEN 
          WRITE(FCT(3:3),FMT='(I1)') I
          FCT(4:8)=')    '
        ELSEIF(I.LT.100) THEN
          WRITE(FCT(3:4),FMT='(I2)') I
          FCT(5:8)=')   '
        ELSE
          WRITE(LU,*) 'DEBSCE NOT PROGRAMMED FOR MORE THAN 99 SOURCES'
          CALL PLANTE(1)
          STOP 
        ENDIF
        CALL READ_FIC_SOURCES(DEBSCE,FCT,AT,T2D_FILES(T2DVEF)%LU,
     *                        ENTET,OK(I)) 
C
      ENDIF
C
C     BEWARE, AN ERROR IN THE SOURCES FILE MAY REMAIN UNNOTICED
C     BECAUSE WE RESORT HERE TO THE PARAMETER FILE
C
      IF(.NOT.OK(I).OR.T2D_FILES(T2DVEF)%NAME(1:1).EQ.' ') THEN
C 
C       PROGRAMMABLE PART                              
C       DISCE IS TAKEN IN THE PARAMETER FILE 
C
C       GLOBAL NUMBER OF SOURCE I IS ISCE(I) IN TELEMAC-2D
        DEBSCE = DISCE(I)                                              
C 
      ENDIF
C          
C-----------------------------------------------------------------------
C
      RETURN
      END
