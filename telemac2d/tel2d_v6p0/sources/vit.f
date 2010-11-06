C                       *****************************
                        DOUBLE PRECISION FUNCTION VIT
C                       *****************************
C
     *( I , N )
C
C***********************************************************************
C  TELEMAC 2D VERSION 6.0    17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C FONCTION  : DONNE LA VALEUR DE LA VITESSE POUR TOUTES LES ENTREES A
C             VITESSE IMPOSEE.
C
C-----------------------------------------------------------------------
C
C FUNCTION  : GIVES THE PRESCRIBED VELOCITY AT ALL LIQUID BOUNDARIES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   I            | -->| NUMBER OF LIQUID BOUNDARY
C |   N            | -->| GLOBAL NUMBER OF POINT  
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
      USE INTERFACE_TELEMAC2D, EX_VIT => VIT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: I,N
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER*8 FCT
      INTEGER J
      LOGICAL DEJA,OK(MAXFRO)
      DATA    DEJA /.FALSE./
      SAVE    OK,DEJA
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
C
      IF(OK(I).AND.T2D_FILES(T2DIMP)%NAME(1:1).NE.' ') THEN
C
C       FCT WILL BE VIT(1), VIT(2), ETC, VIT(99), DEPENDING ON I
        FCT(1:4)='VIT('
        IF(I.LT.10) THEN 
          WRITE(FCT(5:5),FMT='(I1)') I
          FCT(6:8)=')  '
        ELSEIF(I.LT.100) THEN
          WRITE(FCT(5:6),FMT='(I2)') I
          FCT(7:8)=') '
        ELSE
          STOP 'VIT NOT PROGRAMMED FOR MORE THAN 99 BOUNDARIES'
        ENDIF
        CALL READ_FIC_FRLIQ(VIT,FCT,AT,T2D_FILES(T2DIMP)%LU,ENTET,OK(I))
C
      ENDIF
C
      IF(.NOT.OK(I).OR.T2D_FILES(T2DIMP)%NAME(1:1).EQ.' ') THEN
C 
C     PROGRAMMABLE PART                              
C     SL IS TAKEN IN THE PARAMETER FILE, BUT MAY BE CHANGED 
C                                                                                     
        IF(NVITES.GE.I) THEN
          VIT = VITES(I)
        ELSE
          IF(LNG.EQ.1) WRITE(LU,200) I
200       FORMAT(1X,/,1X,'VIT : VITESSES IMPOSEES EN NOMBRE INSUFFISANT'
     *             ,/,1X,'      DANS LE FICHIER DES PARAMETRES'
     *             ,/,1X,'      IL EN FAUT AU MOINS : ',1I6)
          IF(LNG.EQ.2) WRITE(LU,201) I
201       FORMAT(1X,/,1X,'VIT : MORE PRESCRIBED VELOCITIES ARE REQUIRED'
     *             ,/,1X,'      IN THE PARAMETER FILE'
     *             ,/,1X,'      AT LEAST ',1I6,' MUST BE GIVEN')
          CALL PLANTE(1)
          STOP
        ENDIF
C 
      ENDIF             
C
C-----------------------------------------------------------------------
C
      RETURN
      END
