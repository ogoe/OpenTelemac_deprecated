C                       *******************************
                        DOUBLE PRECISION FUNCTION TRSCE
C                       *******************************
C
     *( TIME , I , ITRAC )
C
C***********************************************************************
C  TELEMAC 2D VERSION 6.0   08/04/08  J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C FONCTION  : DONNE LA VALEUR DU TRACEUR POUR TOUTES LES SOURCES
C
C             PERMET DE FAIRE VARIER LA VALEUR AUX SOURCES EN FONCTION
C             DU TEMPS.
C
C-----------------------------------------------------------------------
C
C FUNCTION  : GIVES THE PRESCRIBED VALUE OF TRACERS AT THE SOURCES
C
C             THIS VALUE MAY VARY IN TIME.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   TIME         | -->| TIME
C |   I            | -->| SOURCE RANK
C |   TSCE         | -->| ARRAY OF PRESCRIBED VALUES OF THE TRACER
C |                |    | (READ IN THE PARAMETER FILE) 
C |   ITRAC        | -->| TRACER RANK
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
      USE DECLARATIONS_TELEMAC2D, ONLY: AT,ENTET,NTRAC,TSCE,NREJET,
     *                                  T2D_FILES,T2DVEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(IN) :: TIME
      INTEGER         , INTENT(IN) :: I,ITRAC
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER*8 FCT
      INTEGER N,IRANK
      LOGICAL DEJA,OK(99)  ! 99.GE.NREJET*NTRAC
      DATA    DEJA /.FALSE./
      SAVE    OK,DEJA
C
C     FIRST CALL, OK INITIALISED TO .TRUE.
C
      IF(.NOT.DEJA) THEN
        IF(NREJET*NTRAC.GT.99) THEN
          WRITE(LU,*) 'CHANGE DIMENSION OF OK IN TRSCE, ',NREJET*NTRAC,
     *                ' REQUIRED'
          CALL PLANTE(1)
          STOP 
        ENDIF
        DO N=1,NREJET*NTRAC
          OK(N)=.TRUE.
        ENDDO
        DEJA=.TRUE.
      ENDIF
C
C     IF SOURCES FILE EXISTING, ATTEMPT TO FIND
C     THE VALUE IN IT. IF YES, OK REMAINS TO .TRUE. FOR NEXT CALLS
C                      IF  NO, OK SET     TO .FALSE.
C
C     IRANK CORRESPONDS TO TELEMAC2D DOCUMENTATION
C     TRACER 1 OF SOURCE 1, TRACER 2 OF SOURCE 1, ETC.
      IRANK=ITRAC+NTRAC*(I-1)
      IF(OK(IRANK).AND.T2D_FILES(T2DVEF)%NAME(1:1).NE.' ') THEN
C
C       FCT WILL BE T(1), T(2), ETC, T(99), DEPENDING ON I AND ITRAC
        FCT='TR(     '
        IF(IRANK.LT.10) THEN 
          WRITE(FCT(4:4),FMT='(I1)') IRANK
          FCT(5:5)=')'
        ELSEIF(IRANK.LT.100) THEN
          WRITE(FCT(4:5),FMT='(I2)') IRANK
          FCT(6:6)=')'
        ELSE
          WRITE(LU,*) 'TRSCE NOT PROGRAMMED FOR MORE THAN 99 DATA'
          CALL PLANTE(1)
          STOP 
        ENDIF
        CALL READ_FIC_SOURCES(TRSCE,FCT,AT,T2D_FILES(T2DVEF)%LU,
     *                        ENTET,OK(IRANK)) 
C
      ENDIF
C
C     BEWARE, AN ERROR IN THE SOURCES FILE MAY REMAIN UNNOTICED
C     BECAUSE WE RESORT HERE TO THE PARAMETER FILE
C
      IF(.NOT.OK(I).OR.T2D_FILES(T2DVEF)%NAME(1:1).EQ.' ') THEN
C 
C       PROGRAMMABLE PART                              
C       TSCE IS TAKEN IN THE PARAMETER FILE 
C
C       GLOBAL NUMBER OF SOURCE I IS ISCE(I) IN TELEMAC-2D
        TRSCE = TSCE(I,ITRAC)                                              
C 
      ENDIF       
C
C-----------------------------------------------------------------------
C
      RETURN
      END
