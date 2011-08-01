C                       **************************
                        SUBROUTINE BIEF_OPEN_FILES
C                       **************************
C
     &(CODE,FILES,NFILES,PATH,NCAR,FLOT,IFLOT,ICODE)
C
C***********************************************************************
C BIEF VERSION 6.0     12/10/2009     J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C      FONCTIONS: OPENING FILES DECLARED IN THE PARAMETER FILE
C      ==========
C
C      NOTE : PARAMETER FILE AND DICTIONARY ARE OPENED AND CLOSED
C             IN LECDON
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    CODE        | -->| NAME OF CALLING PROGRAMME
C |    FILES       | -->| STRUCTURES OF CODE FILES
C |    NFILES      | -->| NUMBER OF FILES
C |    PATH        | -->| FULL NAME OF THE PATH WHERE THE CASE IS
C |    NCAR        | -->| NUMBER OF CHARACTERS IN THE PATH
C |    FLOT        | -->| LOGICAL, IF YES LOGICAL UNITS DECIDED BY
C |                |    | THIS SUBROUTINE, IF NO, TAKEN IN SUBMIT
C |    IFLOT       | -->| IF FLOT=YES, START NEW LOGICAL UNIT NUMBERS
C |                |    | AT IFLOT+1
C |    ICODE       |    | NUMERO DU CODE EN CAS DE COUPLAGE
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
      USE BIEF, EX_BIEF_OPEN_FILES => BIEF_OPEN_FILES
      USE DECLARATIONS_TELEMAC
      USE M_MED
C
      IMPLICIT NONE
      INTEGER     LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER           , INTENT(IN)    :: NFILES
      CHARACTER(LEN=24) , INTENT(IN)    :: CODE
      TYPE(BIEF_FILE)   , INTENT(INOUT) :: FILES(NFILES)
      CHARACTER(LEN=250), INTENT(IN)    :: PATH
      INTEGER           , INTENT(IN)    :: NCAR,ICODE
      INTEGER           , INTENT(INOUT) :: IFLOT
      LOGICAL           , INTENT(IN)    :: FLOT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
      CHARACTER(LEN=11) :: FORME
C
C-----------------------------------------------------------------------
C     
C     MESSAGE
C
      IF(LNG.EQ.1) WRITE(LU,*) 'OUVERTURE DES FICHIERS POUR ',CODE
      IF(LNG.EQ.2) WRITE(LU,*) 'OPENING FILES FOR ',CODE
C
C     DECRYPTAGE DE LA CHAINE SUBMIT POUR LES FICHIERS
C     TROUVES DANS LE FICHIER CAS.
C
      DO I=1,NFILES
C
        IF(FILES(I)%NAME(1:1).NE.' ') THEN
C
C         LOGICAL UNIT MODIFIED WHEN COUPLING
C
          IF(FLOT) THEN
            IFLOT=IFLOT+1
C           2 AND 3 SKIPPED (DICTIONARY AND PARAMETER FILE)
            IF(IFLOT.EQ.2) IFLOT=4
C           5 AND 6 SKIPPED (STANDARD INPUT AND OUTPUT)
            IF(IFLOT.EQ.5) IFLOT=7
            FILES(I)%LU=IFLOT
          ENDIF
C
          IF(FILES(I)%BINASC.EQ.'ASC') THEN
            FORME='FORMATTED  '
          ELSE
            FORME='UNFORMATTED'
          ENDIF
C 
C         OUVERTURE DU FICHIER
C
          IF(FILES(I)%FMT.EQ.'MED     ') THEN
C
            IF(NCSIZE.LE.1) THEN
              CALL OPEN_FILE_MED(FILES(I)%TELNAME,FILES(I)%LU,
     *                           FILES(I)%ACTION)
            ELSE
C             PARALLELE, FICHIER DE TYPE SCAL
              IF(FILES(I)%TYPE(1:4).EQ.'SCAL') THEN
              CALL OPEN_FILE_MED(PATH(1:NCAR)//TRIM(FILES(I)%TELNAME)
     *               ,FILES(I)%LU,FILES(I)%ACTION)
C             PARALLELE, FICHIER D'AUTRE TYPE
              ELSE
              CALL OPEN_FILE_MED(PATH(1:NCAR)//TRIM(FILES(I)%TELNAME)
     *               //EXTENS(NCSIZE-1,IPID)//'.med',FILES(I)%LU,
     *               FILES(I)%ACTION)
              ENDIF
            ENDIF              
C
          ELSE
C
            IF(NCSIZE.LE.1) THEN
C             SCALAIRE
              OPEN(FILES(I)%LU,FILE=FILES(I)%TELNAME,
     *             FORM=FORME,ACTION=FILES(I)%ACTION)
            ELSE
C             PARALLELE, FICHIER DE TYPE SCAL
              IF(FILES(I)%TYPE(1:4).EQ.'SCAL') THEN
                OPEN(FILES(I)%LU,
     *               FILE=PATH(1:NCAR)//TRIM(FILES(I)%TELNAME),
     *               FORM=FORME,ACTION=FILES(I)%ACTION)
C             PARALLELE, FICHIER D'AUTRE TYPE
              ELSE
                OPEN(FILES(I)%LU,
     *               FILE=PATH(1:NCAR)//TRIM(FILES(I)%TELNAME)
     *               //EXTENS(NCSIZE-1,IPID),
     *               FORM=FORME,ACTION=FILES(I)%ACTION)
              ENDIF
            ENDIF      
C
          ENDIF
C
        ENDIF
C
      ENDDO
C
C     SETTING AND STORING NAME OF CODE
C
      NAMECODE = CODE
      NNAMECODE(ICODE) = CODE              
C
C-----------------------------------------------------------------------
C
      RETURN
      END
