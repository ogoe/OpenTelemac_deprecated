C                       **********************
                        SUBROUTINE READ_CONFIG
C                       **********************
C
     *(LNG,LU,CHAINE,NCAR)
C
C***********************************************************************
C BIEF VERSION 5.1             J-M HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C      FONCTIONS:
C      ==========
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |                | -->|
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : HOMERE
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER      , INTENT(INOUT) :: LNG,LU
      CHARACTER*250, INTENT(IN)    :: CHAINE
      INTEGER      , INTENT(IN)    :: NCAR
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      LOGICAL YACONFIG
      INTEGER NC
      CHARACTER*257 CONFIG
C
C-----------------------------------------------------------------------
C
      IF(NCAR.GT.0) THEN
        CONFIG(1:NCAR+6)=CHAINE(1:NCAR) // 'CONFIG'
        NC=NCAR+6
      ELSE
        CONFIG(1:6)='CONFIG'
        NC=6
      ENDIF
C   
      YACONFIG=.FALSE.
      INQUIRE(FILE=CONFIG(1:NC),EXIST=YACONFIG)
      IF(YACONFIG) THEN
C
        OPEN(40,FILE=CONFIG(1:NC), FORM='FORMATTED')
        READ(40,*) LNG
C$DC$ : ne pas surcharger LU en mode PARALLEL sous WinNT
C       (garder la redirection output faite par P_INIT sur
C       le canal 95)
        IF(LU.NE.95) READ(40,*) LU
        CLOSE(40)
C
      ELSE
C
        WRITE(LU,*) 'READ_CONFIG : FICHIER CONFIG NON TROUVE : ',CONFIG
        WRITE(LU,*) 'VALEURS PAR DEFAUT DE LU ET LNG : ',LU,' ET ',LNG
        WRITE(LU,*) ' '
        WRITE(LU,*) 'READ_CONFIG: FILE CONFIG NOT FOUND: ',CONFIG
        WRITE(LU,*) 'DEFAULTS VALUES OF LU AND LNG: ',LU,' AND ',LNG
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
