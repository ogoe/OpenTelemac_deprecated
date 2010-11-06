C                       **********************
                        SUBROUTINE READ_SUBMIT
C                       **********************
C
     &(FILES,NFILES,CODE,SUBMIT,NMOT)
C
C***********************************************************************
C BIEF VERSION 6.0      27/03/2009    J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C      FONCTIONS: OPENING FILES DECLARED IN THE PARAMETER FILE
C      ==========
C
C      NOTE : PARAMETER FILE AND DICTIONARY ARE OPENED
C             AND CLOSED IN LECDON
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    FILES       | -->| FILES STRUCTURES
C |    NFILES      | -->| NUMBER OF FILES IN ARRAY FILES
C |    CODE        | -->| NAME OF CALLING PROGRAMME
C |    SUBMIT      | -->| CHARACTER STRINGS STEMMING FROM DICTIONARY
C |    NMOT        | -->| SECOND DIMENSION OF SUBMIT AND MOTCAR 
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
      USE BIEF, EX_READ_SUBMIT => READ_SUBMIT
      USE DECLARATIONS_TELEMAC
C
      IMPLICIT NONE
      INTEGER     LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER           , INTENT(IN) :: NFILES,NMOT
      TYPE(BIEF_FILE), INTENT(INOUT) :: FILES(NFILES)
      CHARACTER(LEN=24) , INTENT(IN) :: CODE
      CHARACTER(LEN=144), INTENT(IN) :: SUBMIT(4,NMOT)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,ICOL,I1,CANAL
C
      CHARACTER(LEN=7) :: NOMCANAL
      CHARACTER(LEN=9) :: LITECR
C
      INTEGER  PREVAL,INTLU
      EXTERNAL PREVAL,INTLU
C
C-----------------------------------------------------------------------
C
      DO I=1,NFILES
        FILES(I)%LU=0
        FILES(I)%TELNAME='      '
        FILES(I)%NAME(1:1)=' '
      ENDDO
C
C-----------------------------------------------------------------------
C
C     DECRYPTAGE DE LA CHAINE SUBMIT POUR LES FICHIERS
C     TROUVES DANS LE FICHIER CAS.
C
      DO I=1,NMOT
C
C EXEMPLE DE SUBMIT : 'NGEO-READ-01;T2DGEO;OBLIG;BIN;LIT;SELAFIN-GEOM'
C
        IF(     SUBMIT(4,I).NE.' '
C       IF(     SUBMIT(4,I).NE.' '.AND.MOTCAR(I)(1:1).NE.' '
     *     .AND.SUBMIT(4,I)(1:7).NE.'INUTILE'  ) THEN
C         RECHERCHE DU NOM FORTRAN DU CANAL (PAR EXEMPLE NGEO)
          ICOL=PREVAL(1,SUBMIT(4,I),'-','-','-')
          NOMCANAL=SUBMIT(4,I)(1:ICOL-1)
C         RECHERCHE DE LA CHAINE READ OU WRITE OU READWRITE
C         QUI SE TROUVE AVANT LE SIGNE - SUIVANT
          I1=ICOL+1
          ICOL=PREVAL(I1,SUBMIT(4,I),'-','-','-')
          LITECR=SUBMIT(4,I)(I1:ICOL-1)
C         LECTURE DU CANAL APRES LE SIGNE -
          CANAL=INTLU(ICOL,SUBMIT(4,I))
          IF(CANAL.GT.NFILES) THEN
            IF(LNG.EQ.1) THEN
              WRITE(LU,*) 'READ_SUBMIT : NFILES TROP PETIT : ',NFILES
              WRITE(LU,*) '              IL FAUT AU MOINS ',CANAL
            ENDIF
            IF(LNG.EQ.2) THEN
              WRITE(LU,*) 'READ_SUBMIT: NFILES TOO SMALL : ',NFILES
              WRITE(LU,*) '             IT SHOULD BE AT LEAST ',CANAL
            ENDIF
            CALL PLANTE(1)
            STOP
          ENDIF
          FILES(CANAL)%LU=CANAL
          FILES(CANAL)%ACTION=LITECR
C
          ICOL=PREVAL(ICOL,SUBMIT(4,I),';',';',';')
C         LECTURE DU NOM DU FICHIER A METTRE DANS LE DOSSIER TMP
          I1=PREVAL(ICOL+1,SUBMIT(4,I),';',';',';')
          FILES(CANAL)%TELNAME=SUBMIT(4,I)(ICOL+1:I1-1)
C         ON SAUTE ;FACUL; OU ;OBLIG;
          ICOL=PREVAL(I1+1,SUBMIT(4,I),';',';',';')
C         BINAIRE OU ASCII
          FILES(CANAL)%BINASC=SUBMIT(4,I)(ICOL+1:ICOL+3)
C         REMARQUE : SUBMIT(4,I)(ICOL+5:ICOL+7) CONTAINS LIT OU ECR
C                    NOT USED HERE
C         MODE SELAFIN-GEOM, PARAL, SCAL, ETC.
          FILES(CANAL)%TYPE=TRIM(SUBMIT(4,I)(ICOL+9:MIN(144,ICOL+20)))          
C        
        ENDIF
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
