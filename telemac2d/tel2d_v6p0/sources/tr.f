C                       ****************************
                        DOUBLE PRECISION FUNCTION TR
C                       ****************************
C
     *( I , ITRAC , N , IERR )
C
C***********************************************************************
C TELEMAC 2D VERSION 6.0  02/04/2009  J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C FONCTION  : DONNE LA VALEUR DU TRACEUR POUR TOUTES LES ENTREES OU IL
C             EST IMPOSE.
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C FUNCTION  : GIVES THE PRESCRIBED VALUES OF TRACERS AT
C             THE LIQUID BOUNDARIES.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   I            | -->| RANG DE LA FRONTIERE A DEBIT IMPOSE          |
C |                |    | (1 S'IL N'Y EN A QU'UNE)                     |
C |   N            | -->| NUMERO GLOBAL DU POINT   
C |                |    | (LU DANS LE FICHIER DES PARAMETRES)          |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : BORD
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D
      USE INTERFACE_TELEMAC2D, EX_TR => TR
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: I,N,ITRAC
      INTEGER, INTENT(INOUT) :: IERR
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER*8 FCT
      INTEGER J,IRANK
      LOGICAL DEJA,OK(MAXFRO*MAXTRA)
      DATA    DEJA /.FALSE./
      SAVE    OK,DEJA
C
C     WE A PRIORI ASSUME THAT TR WILL BE FOUND
C
      IERR=0
C
C     FIRST CALL, OK INITIALISED TO .TRUE.
C
      IF(.NOT.DEJA) THEN
        DO J=1,MAXFRO
          OK(J)=.TRUE.
        ENDDO
        DEJA=.TRUE.
      ENDIF
C
C     IF FILE OF LIQUID BOUNDARIES EXISTING, ATTEMPT TO FIND
C     THE VALUE IN IT. IF YES, OK REMAINS TO .TRUE. FOR NEXT CALLS
C                      IF  NO, OK SET     TO .FALSE.
C     RANK OF VALUE IN ARRAY TRACER OR IN LIQUID BOUNDARIES FILE
C
      IRANK=ITRAC+(I-1)*NTRAC
      IF(OK(IRANK).AND.T2D_FILES(T2DIMP)%NAME(1:1).NE.' ') THEN
C
C       FCT WILL BE TR(1), TR(2), ETC, TR(99), DEPENDING ON I
        FCT(1:3)='TR('
        IF(IRANK.LT.10) THEN 
          WRITE(FCT(4:4),FMT='(I1)') IRANK
          FCT(5:8)=')   '
        ELSEIF(IRANK.LT.100) THEN
          WRITE(FCT(4:5),FMT='(I2)') IRANK
          FCT(6:8)=')  '
        ELSE
          WRITE(LU,*) 'TR NOT PROGRAMMED FOR MORE THAN 99 BOUNDARIES'
          CALL PLANTE(1)
          STOP
        ENDIF
        CALL READ_FIC_FRLIQ(TR,FCT,AT,T2D_FILES(T2DIMP)%LU,
     *                      ENTET,OK(IRANK))
C
      ENDIF
C
C     IF VALUE NOT FOUND IN LIQUID BOUNDARIES FILE
C     OR NO LIQUID BOUNDARIES FILE
C     ATTEMPT TO FIND IT IN THE PARAMETER FILE
C
      IF(.NOT.OK(IRANK).OR.T2D_FILES(T2DIMP)%NAME(1:1).EQ.' ') THEN
C                                                                           
        IF(NTRACE.GE.IRANK) THEN
          TR = TRACER(IRANK)
          OK(IRANK)=.TRUE.
        ELSEIF(NTRACE.NE.0) THEN
          IF(LNG.EQ.1) WRITE(LU,300) IRANK
300       FORMAT(1X,/,1X,'TR : VALEURS IMPOSEES DU TRACEUR'
     *             ,/,1X,'     EN NOMBRE INSUFFISANT'
     *             ,/,1X,'     DANS LE FICHIER DES PARAMETRES'
     *             ,/,1X,'     IL EN FAUT AU MOINS : ',1I6)
          IF(LNG.EQ.2) WRITE(LU,301) IRANK
301       FORMAT(1X,/,1X,'TR : MORE PRESCRIBED TRACER VALUES'
     *             ,/,1X,'     ARE REQUIRED IN THE PARAMETER FILE'
     *             ,/,1X,'     AT LEAST ',1I6,' MUST BE GIVEN')
          CALL PLANTE(1)
          STOP
C
        ENDIF
C 
      ENDIF            
C
C     NOTHING FOUND: VALUES OF BOUNDARY CONDITIONS FILE WILL BE TAKEN
C
      IF(.NOT.OK(IRANK)) THEN
        TR=0.D0
        IERR=1
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
