C                       ***************************
                        SUBROUTINE FLUXPR_TELEMAC2D
C                       ***************************
C
     *(NSEC,CTRLSC,FLX,VOLNEG,VOLPOS,INFO,TPS,NSEG,NCSIZE,CUMFLO)
C
C***********************************************************************
C  BIEF VERSION 5.5    25/03/99    J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C  FONCTION : CALCUL DES FLUX A TRAVERS DES SECTIONS DE CONTROLE
C             ET CUMUL DE CES FLUX POUR OBTENIR LES VOLUMES OSCILLANTS.
C
C  PRINTOUTS OF DISCHARGES THROUGH CONTROL SECTIONS ARE IN THIS ROUTINE
C  YOU CAN REWRITE IT TO DIVERT THESE PRINTOUTS TO A FILE OR TO CHANGE
C  THE FORMAT
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   NSEC         | -->|  NUMBER OF CONTROL SECTIONS
C |   CTRLSC       | -->|  NUMBERS OF POINTS IN THE CONTROL SECTIONS
C |   FLX          | -->|  FLUXES THROUGH CONTROL SECTIONS
C |   VOLNEG       | -->|  CUMULATED NEGATIVE VOLUME THROUGH SECTIONS
C |   VOLPOS       | -->|  CUMULATED POSITIVE VOLUME THROUGH SECTIONS
C |   INFO         | -->|  IF YES : INFORMATION IS PRINTED
C |   TPS          | -->|  TEMPS
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR :
C
C SOUS-PROGRAMMES APPELES :
C
C***********************************************************************
C
      USE BIEF_DEF, ONLY: IPID
      USE DECLARATIONS_TELEMAC2D, ONLY: 
     &          T2D_FILES,T2DSEC,T2DSEO,CHAIN,TITCAS
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)          :: NSEC,NCSIZE
      INTEGER, INTENT(IN)          :: CTRLSC(*)
      INTEGER, INTENT(IN)          :: NSEG(NSEC)
      LOGICAL, INTENT(IN)          :: INFO,CUMFLO
      DOUBLE PRECISION, INTENT(IN) :: FLX(NSEC)
      DOUBLE PRECISION, INTENT(IN) :: VOLNEG(NSEC),VOLPOS(NSEC)
      DOUBLE PRECISION, INTENT(IN) :: TPS
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: WORK(:)
      DOUBLE PRECISION P_DMAX,P_DMIN, P_DSUM
      INTEGER                        P_IMIN
      EXTERNAL         P_DMAX,P_DMIN,P_DSUM,P_IMIN
C
      INTEGER ISEC,II,ERR
      CHARACTER(LEN=16) :: FMTZON='(4(1X,1PG21.14))'
      LOGICAL :: CLASSIC=.FALSE. 
      LOGICAL, SAVE :: INIT=.TRUE.
      INTEGER, SAVE :: NSEO
C
C-----------------------------------------------------------------------
C 
      IF (.NOT.ALLOCATED(CHAIN)) CLASSIC=.TRUE.
C
      IF(INFO) THEN
C
      IF (CLASSIC) THEN ! follow fluxpr.f of BIEF blindly 
C
      IF(NCSIZE.LE.1) THEN
C
      DO ISEC = 1,NSEC
C
      IF(CUMFLO) THEN
      IF(LNG.EQ.1) WRITE(LU,130) ISEC,CTRLSC(1+2*(ISEC-1)),
     *                                CTRLSC(2+2*(ISEC-1)),
     *                                FLX(ISEC),
     *                                VOLNEG(ISEC),
     *                                VOLPOS(ISEC)
      IF(LNG.EQ.2) WRITE(LU,131) ISEC,CTRLSC(1+2*(ISEC-1)),
     *                                CTRLSC(2+2*(ISEC-1)),
     *                                FLX(ISEC),
     *                                VOLNEG(ISEC),
     *                                VOLPOS(ISEC)
      ELSE
      IF(LNG.EQ.1) WRITE(LU,136) ISEC,CTRLSC(1+2*(ISEC-1)),
     *                                CTRLSC(2+2*(ISEC-1)),
     *                                FLX(ISEC)
      IF(LNG.EQ.2) WRITE(LU,137) ISEC,CTRLSC(1+2*(ISEC-1)),
     *                                CTRLSC(2+2*(ISEC-1)),
     *                                FLX(ISEC)
      ENDIF
130   FORMAT(1X,/,1X,'SECTION DE CONTROLE ',1I2,
     *               ' (ENTRE LES POINTS ',1I5,' ET ',1I5,')',//,5X,
     *               'DEBIT : '                    ,G16.7,/,5X,
     *               'CUMUL DES DEBITS NEGATIFS : ',G16.7,/,5X,
     *               'CUMUL DES DEBITS POSITIFS : ',G16.7)
131   FORMAT(1X,/,1X,'CONTROL SECTION NUMBER ',1I2,
     *               ' (BETWEEN POINTS ',1I5,' AND ',1I5,')',//,5X,
     *               'DISCHARGE: '                 ,G16.7,/,5X,
     *               'NEGATIVE VOLUME THROUGH THE SECTION: ',G16.7,/,5X,
     *               'POSITIVE VOLUME THROUGH THE SECTION: ',G16.7)
136   FORMAT(1X,/,1X,'SECTION DE CONTROLE ',1I2,
     *               ' (ENTRE LES POINTS ',1I5,' ET ',1I5,')',//,5X,
     *               'DEBIT : '                    ,G16.7)
137   FORMAT(1X,/,1X,'CONTROL SECTION NUMBER ',1I2,
     *               ' (BETWEEN POINTS ',1I5,' AND ',1I5,')',//,5X,
     *               'DISCHARGE: '                 ,G16.7)
C
      ENDDO
C
      ELSE
C
      DO ISEC = 1,NSEC
C     SECTIONS ACROSS 2 SUB-DOMAINS WILL HAVE NSEG=0 OR -1
C     AND -1 WANTED HERE FOR RELEVANT MESSAGE.
      II=P_IMIN(NSEG(ISEC))
C
      IF(II.GE.0) THEN
C
      IF(LNG.EQ.1) WRITE(LU,132) ISEC,CTRLSC(1+2*(ISEC-1)),
     *                                CTRLSC(2+2*(ISEC-1)),
     *              P_DMIN(FLX(ISEC))+P_DMAX(FLX(ISEC)),
     *                                P_DMIN(VOLNEG(ISEC)),
     *                                P_DMAX(VOLPOS(ISEC))
      IF(LNG.EQ.2) WRITE(LU,133) ISEC,CTRLSC(1+2*(ISEC-1)),
     *                                CTRLSC(2+2*(ISEC-1)),
     *              P_DMIN(FLX(ISEC))+P_DMAX(FLX(ISEC)),
     *                                P_DMIN(VOLNEG(ISEC)),
     *                                P_DMAX(VOLPOS(ISEC))
132   FORMAT(1X,/,1X,'SECTION DE CONTROLE ',1I2,
     *               ' (ENTRE LES POINTS ',1I5,' ET ',1I5,')',//,5X,
     *               'DEBIT : '                    ,G16.7,/,5X,
     *               'CUMUL DES DEBITS NEGATIFS : ',G16.7,/,5X,
     *               'CUMUL DES DEBITS POSITIFS : ',G16.7)
133   FORMAT(1X,/,1X,'CONTROL SECTION NUMBER ',1I2,
     *               ' (BETWEEN POINTS ',1I5,' AND ',1I5,')',//,5X,
     *               'DISCHARGE: '                 ,G16.7,/,5X,
     *               'NEGATIVE VOLUME THROUGH THE SECTION: ',G16.7,/,5X,
     *               'POSITIVE VOLUME THROUGH THE SECTION: ',G16.7)
C
      ELSE
C
      IF(LNG.EQ.1) WRITE(LU,134) ISEC,CTRLSC(1+2*(ISEC-1)),
     *                                CTRLSC(2+2*(ISEC-1))
      IF(LNG.EQ.2) WRITE(LU,135) ISEC,CTRLSC(1+2*(ISEC-1)),
     *                                CTRLSC(2+2*(ISEC-1))
134   FORMAT(1X,/,1X,'SECTION DE CONTROLE ',1I2,
     *               ' (ENTRE LES POINTS ',1I5,' ET ',1I5,')',//,5X,
     *               'A CHEVAL SUR DEUX SOUS-DOMAINES, PAS DE CALCUL')
135   FORMAT(1X,/,1X,'CONTROL SECTION NUMBER ',1I2,
     *               ' (BETWEEN POINTS ',1I5,' AND ',1I5,')',//,5X,
     *               'ACROSS TWO SUB-DOMAINS, NO COMPUTATION')
      ENDIF
C
      ENDDO
C
      ENDIF
C
C-----------------------------------------------------------------------
C CHAIN allocated, i.e. serial or parallel case from SECTIONS INPUT FILE
C       we can apply co-ordinates instead and/or names of sections
C
      ELSE 
        IF(NCSIZE.LE.1) THEN ! serial
          DO ISEC = 1,NSEC
            IF(CUMFLO) THEN
              IF(LNG.EQ.1) WRITE(LU,230) ISEC,TRIM(CHAIN(ISEC)%DESCR),
     *                           FLX(ISEC),VOLNEG(ISEC),VOLPOS(ISEC)
              IF(LNG.EQ.2) WRITE(LU,231) ISEC,TRIM(CHAIN(ISEC)%DESCR),
     *                           FLX(ISEC),VOLNEG(ISEC),VOLPOS(ISEC)
            ELSE
              IF(LNG.EQ.1) WRITE(LU,236) ISEC,TRIM(CHAIN(ISEC)%DESCR),
     *                           FLX(ISEC)
              IF(LNG.EQ.2) WRITE(LU,237) ISEC,TRIM(CHAIN(ISEC)%DESCR),
     *                           FLX(ISEC)
            ENDIF

230   FORMAT(1X,/,1X,'SECTION DE CONTROLE ',1I2,
     *               ' (NOM ',A,')',//,5X,
     *               'DEBIT : '                    ,G16.7,/,5X,
     *               'CUMUL DES DEBITS NEGATIFS : ',G16.7,/,5X,
     *               'CUMUL DES DEBITS POSITIFS : ',G16.7)
231   FORMAT(1X,/,1X,'CONTROL SECTION NUMBER ',1I2,
     *               ' (NAME ',A,')',//,5X,
     *               'DISCHARGE: '                 ,G16.7,/,5X,
     *               'NEGATIVE VOLUME THROUGH THE SECTION: ',G16.7,/,5X,
     *               'POSITIVE VOLUME THROUGH THE SECTION: ',G16.7)
236   FORMAT(1X,/,1X,'SECTION DE CONTROLE ',1I2,
     *               ' (NOM ',A,')',//,5X,
     *               'DEBIT : '                    ,G16.7)
237   FORMAT(1X,/,1X,'CONTROL SECTION NUMBER ',1I2,
     *               ' (NAME ',A,')',//,5X,
     *               'DISCHARGE: '                 ,G16.7)
          ENDDO
C
        ELSE
C
          DO ISEC = 1,NSEC
C
            IF(LNG.EQ.1) WRITE(LU,232) ISEC,TRIM(CHAIN(ISEC)%DESCR),
     *                                P_DSUM(FLX(ISEC)),
     *                                P_DSUM(VOLNEG(ISEC)),
     *                                P_DSUM(VOLPOS(ISEC))
            IF(LNG.EQ.2) WRITE(LU,233) ISEC,TRIM(CHAIN(ISEC)%DESCR),
     *                                P_DSUM(FLX(ISEC)),
     *                                P_DSUM(VOLNEG(ISEC)),
     *                                P_DSUM(VOLPOS(ISEC))

232         FORMAT(1X,/,1X,'SECTION DE CONTROLE ',1I2,
     *               ' (NOM ',A,')',//,5X,
     *               'DEBIT : '                    ,G16.7,/,5X,
     *               'CUMUL DES DEBITS NEGATIFS : ',G16.7,/,5X,
     *               'CUMUL DES DEBITS POSITIFS : ',G16.7)
233         FORMAT(1X,/,1X,'CONTROL SECTION NUMBER ',1I2,
     *               ' (NAME ',A,')',//,5X,
     *               'DISCHARGE: '                 ,G16.7,/,5X,
     *               'NEGATIVE VOLUME THROUGH THE SECTION: ',G16.7,/,5X,
     *               'POSITIVE VOLUME THROUGH THE SECTION: ',G16.7)
C
          ENDDO
        ENDIF
C
      ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C master writes a nice sections output file, the header only once
C
      IF ( (.NOT.CLASSIC) .AND. 
     &      (TRIM(T2D_FILES(T2DSEO)%NAME).NE.'') ) THEN 
        IF (INIT) THEN 
          INIT=.FALSE.
          IF ((NCSIZE.GT.1 .AND. IPID.EQ.0).OR.(NCSIZE.LE.1)) THEN 
            NSEO=T2D_FILES(T2DSEO)%LU
            WRITE(NSEO,*) 'TITLE = "Fluxes for ',TRIM(TITCAS),'"' 
            WRITE(NSEO,*) 'VARIABLES = t', 
     &         (' '//TRIM(CHAIN(ISEC)%DESCR), ISEC=1,NSEC)        
          ENDIF
          IF (NCSIZE.GT.1) THEN 
            ALLOCATE (WORK(NSEC), STAT=ERR)
            IF (ERR.NE.0) THEN
              WRITE(LU,*) 
     &          'FLUXPR_TELEMAC2D: error allocating WORK:',ERR
              CALL PLANTE(1)
              STOP
            ENDIF
          ENDIF
        ENDIF
        ! deadlock with write and p_dsum in an implied WRITE loop 
        ! because it is only master to write the message...
        IF (NCSIZE.GT.1) THEN
          DO ISEC=1,NSEC
            WORK(ISEC) = P_DSUM(FLX(ISEC))
          END DO 
          IF (IPID.EQ.0) 
     &      WRITE (NSEO, FMT=FMTZON) TPS, (WORK(ISEC), ISEC=1,NSEC)
        ELSE
          WRITE (NSEO, FMT=FMTZON) TPS, (FLX(ISEC), ISEC=1,NSEC)
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C 
      RETURN
      END 
