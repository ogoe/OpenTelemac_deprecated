C                       ******************
                        SUBROUTINE WAITFOR
C                       ******************
C
     * (DOSSIER,FICHIER)                                            
C
C***********************************************************************
C  BIEF VERSION 5.2    08/02/2001    NATHALY BARBRY (UNIVERSITE DE CAEN)
C                                        J-M HERVOUET (LNHE) 30 87 80 18
C***********************************************************************
C
C  FONCTION : 'FICHIER' EST UN NOM DE FICHIER ATTENDU DANS LE DOSSIER.
C             SA PRESENCE EST INDIQUEE PAR L'EXISTENCE D'UN FICHIER VIDE
C             DONT LE NOM EST LE MEME QUE 'FICHIER' PRECEDE DE YA.
C
C             QUAND LE FICHIER YAFICHIER EXISTE, ON LE DETRUIT ET ON SORT
C
C
C  FUNCTION : FICHIER IS THE NAME OF A FILE WHICH IS WAITED IN DIRECTOTY
C             DOSSIER. THIS FILE EXISTS IF THERE EXISTS ANOTHER EMPTY
C             FILE CALLED YAFICHIER (NAME OF THE FILE WITH 'YA' BEFORE).
C
C             WHEN FILE YAFICHIER EXISTS, IT IS DELETED AND WE RETURN TO
C             THE CALLING PROGRAM.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   DOSSIER      | -->| DOSSIER OU SE TROUVENT LES FICHIERS A LIRE   |
C |   FICHIER      | -->| FICHIER A LIRE                               |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C----------------------------------------------------------------------
C
C  APPELE PAR : COUPLAGE
C
C**********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      CHARACTER(LEN=250), INTENT(IN) :: DOSSIER
      CHARACTER(LEN=*)  , INTENT(IN) :: FICHIER
C
      CHARACTER(LEN=270) NOMFIC
      LOGICAL OUI
      INTEGER ERR1
C
C     TEMPS D'ATTENTE PARAMETRABLE
C
      INTEGER, PARAMETER :: LAPS = 3
C
C-----------------------------------------------------------------------
C
      INTRINSIC TRIM
C
C-----------------------------------------------------------------------
C
10    CONTINUE
C
      NOMFIC=TRIM(DOSSIER)//'YA'//FICHIER
      INQUIRE(FILE=NOMFIC,EXIST=OUI,ERR=84,IOSTAT=ERR1)
C
C-----------------------------------------------------------------------
C  
      IF(OUI) THEN
C
        OPEN(94,FILE=NOMFIC,
     *          STATUS='OLD',FORM='UNFORMATTED',ERR=85,IOSTAT=ERR1)
        CLOSE(94,STATUS='DELETE',ERR=86,IOSTAT=ERR1)
C
        GO TO 1000
C
      ELSE
C
        INQUIRE(FILE=TRIM(DOSSIER)//'STOP',EXIST=OUI,ERR=84,IOSTAT=ERR1)
        IF(OUI) THEN
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'WAITFOR : ARRET DU PROGRAMME, ATTENTE INUTILE'
            WRITE(LU,*) '          CAR UN FICHIER STOP A ETE CREE'
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'WAITFOR : PROGRAM STOPPED, WAITING IS USELESS'
            WRITE(LU,*) '          BECAUSE A FILE STOP HAS BEEN CREATED'
          ENDIF
          CALL PLANTE(1)
          STOP
        ELSE
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'ATTENTE DE ',LAPS,' SECONDES'
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'WAITING ',LAPS,' SECONDS'
          ENDIF
          CALL ATTEND(LAPS) 
        ENDIF
C
        GO TO 10
C
      ENDIF 
C
C-----------------------------------------------------------------------
C     MESSAGES D'ERREUR
C-----------------------------------------------------------------------
C
84    CONTINUE
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'WAITFOR : ERREUR DE INQUIRE SUR LE FICHIER :'
        WRITE(LU,*) '         ',TRIM(DOSSIER)//'YA'//FICHIER
        WRITE(LU,*) '          ERREUR NUMERO : ',ERR1
      ENDIF
      IF(LNG.EQ.2) THEN
        WRITE(LU,*) 'WAITFOR : ERROR OF COMMAND INQUIRE ON FILE :'
        WRITE(LU,*) '         ',TRIM(DOSSIER)//'YA'//FICHIER
        WRITE(LU,*) '          ERROR NUMBER : ',ERR1
      ENDIF
      CALL PLANTE(1)
      STOP
85    CONTINUE
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'WAITFOR : ERREUR ',ERR1
        WRITE(LU,*) 'A L''OUVERTURE DU FICHIER ',NOMFIC
      ENDIF
      IF(LNG.EQ.2) THEN
        WRITE(LU,*) 'WAITFOR : ERROR ',ERR1
        WRITE(LU,*) 'WHEN OPENING THE FILE ',NOMFIC
      ENDIF
      CALL PLANTE(1)
      STOP
86    CONTINUE
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'WAITFOR : ERREUR ',ERR1
        WRITE(LU,*) 'A LA FERMETURE DU FICHIER ',NOMFIC
      ENDIF
      IF(LNG.EQ.2) THEN
        WRITE(LU,*) 'WAITFOR : ERROR ',ERR1
        WRITE(LU,*) 'WHEN CLOSING THE FILE ',NOMFIC
      ENDIF
      CALL PLANTE(1)
      STOP
C
C-----------------------------------------------------------------------
C
1000  CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
