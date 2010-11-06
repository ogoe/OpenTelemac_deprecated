C                       *************************
                        SUBROUTINE FLUXPR_SISYPHE
C                       *************************
C
     *(NSEC,CTRLSC,FLX,VOLNEG,VOLPOS,INFO,TPS,NSEG,NCSIZE,
     * FLXS,VOLNEGS,VOLPOSS,SUSP,FLXC,VOLNEGC,VOLPOSC,CHARR)
C
C***********************************************************************
C  SISYPHE VERSION 5.7    27/12/06    J-M HERVOUET (LNHE) 01 30 87 80 18
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
      USE DECLARATIONS_SISYPHE, ONLY:
     &          SIS_FILES,SISSEO,CHAIN,TITCA
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)          :: NSEC,NCSIZE
      INTEGER, INTENT(IN)          :: CTRLSC(*)
      INTEGER, INTENT(IN)          :: NSEG(NSEC)
      LOGICAL, INTENT(IN)          :: INFO,SUSP,CHARR
      DOUBLE PRECISION, INTENT(IN) :: FLX(NSEC),TPS
      DOUBLE PRECISION, INTENT(IN) :: VOLNEG(NSEC),VOLPOS(NSEC)
      DOUBLE PRECISION, INTENT(IN) :: FLXS(NSEC),FLXC(NSEC)
      DOUBLE PRECISION, INTENT(IN) :: VOLNEGS(NSEC),VOLPOSS(NSEC)
      DOUBLE PRECISION, INTENT(IN) :: VOLNEGC(NSEC),VOLPOSC(NSEC)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: WORK(:)
      DOUBLE PRECISION P_DMAX,P_DMIN,P_DSUM
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
      IF (CLASSIC) THEN !jaj #### follow fluxpr.f of BIEF blindly
C
      IF(NCSIZE.LE.1) THEN
C
      DO ISEC = 1,NSEC
C
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
130   FORMAT(1X,/,1X,'SECTION DE CONTROLE ',1I2,
     *               ' (ENTRE LES POINTS ',1I5,' ET ',1I5,')',//,5X,
     *               'DEBIT :                     ',G16.7,/,5X,
     *               'CUMUL DES DEBITS NEGATIFS : ',G16.7,/,5X,
     *               'CUMUL DES DEBITS POSITIFS : ',G16.7)
131   FORMAT(1X,/,1X,'CONTROL SECTION NUMBER ',1I2,
     *               ' (BETWEEN POINTS ',1I5,' AND ',1I5,')',//,5X,
     *               'DISCHARGE:                 ',G16.7,/,5X,
     *               'CUMULATED NEGATIVE VOLUME: ',G16.7,/,5X,
     *               'CUMULATED POSITIVE VOLUME: ',G16.7)
      IF(SUSP) THEN
        IF(LNG.EQ.1) WRITE(LU,1301) FLXS(ISEC),
     *                              VOLNEGS(ISEC),
     *                              VOLPOSS(ISEC)
        IF(LNG.EQ.2) WRITE(LU,1302) FLXS(ISEC),
     *                              VOLNEGS(ISEC),
     *                              VOLPOSS(ISEC)
1301    FORMAT(5X,'DEBIT EN SUSPENSION :       ',G16.7,/,5X,
     *            'CUMUL DES DEBITS NEGATIFS : ',G16.7,/,5X,
     *            'CUMUL DES DEBITS POSITIFS : ',G16.7)
1302    FORMAT(5X,'DISCHARGE IN SUSPENSION:   ',G16.7,/,5X,
     *            'CUMULATED NEGATIVE VOLUME: ',G16.7,/,5X,
     *            'CUMULATED POSITIVE VOLUME: ',G16.7)
      ENDIF
      IF(CHARR) THEN
        IF(LNG.EQ.1) WRITE(LU,1303) FLXC(ISEC),
     *                              VOLNEGC(ISEC),
     *                              VOLPOSC(ISEC)
        IF(LNG.EQ.2) WRITE(LU,1304) FLXC(ISEC),
     *                              VOLNEGC(ISEC),
     *                              VOLPOSC(ISEC)
1303    FORMAT(5X,'DEBIT EN CHARRIAGE :        ',G16.7,/,5X,
     *            'CUMUL DES DEBITS NEGATIFS : ',G16.7,/,5X,
     *            'CUMUL DES DEBITS POSITIFS : ',G16.7)
1304    FORMAT(5X,'BEDLOAD DISCHARGE:         ',G16.7,/,5X,
     *            'CUMULATED NEGATIVE VOLUME: ',G16.7,/,5X,
     *            'CUMULATED POSITIVE VOLUME: ',G16.7)
      ENDIF
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
      IF(LNG.EQ.1) WRITE(LU,130) ISEC,CTRLSC(1+2*(ISEC-1)),
     *                                CTRLSC(2+2*(ISEC-1)),
     *                 P_DMIN(FLX(ISEC))+P_DMAX(FLX(ISEC)),
     *                                P_DMIN(VOLNEG(ISEC)),
     *                                P_DMAX(VOLPOS(ISEC))
      IF(LNG.EQ.2) WRITE(LU,131) ISEC,CTRLSC(1+2*(ISEC-1)),
     *                                CTRLSC(2+2*(ISEC-1)),
     *                 P_DMIN(FLX(ISEC))+P_DMAX(FLX(ISEC)),
     *                                P_DMIN(VOLNEG(ISEC)),
     *                                P_DMAX(VOLPOS(ISEC))
      IF(SUSP) THEN
        IF(LNG.EQ.1) WRITE(LU,1301) 
     *              P_DMIN(FLXS(ISEC))+P_DMAX(FLXS(ISEC)),
     *                              P_DMIN(VOLNEGS(ISEC)),
     *                              P_DMAX(VOLPOSS(ISEC))
        IF(LNG.EQ.2) WRITE(LU,1302) 
     *              P_DMIN(FLXS(ISEC))+P_DMAX(FLXS(ISEC)),
     *                              P_DMIN(VOLNEGS(ISEC)),
     *                              P_DMAX(VOLPOSS(ISEC))
      ENDIF
      IF(CHARR) THEN
        IF(LNG.EQ.1) WRITE(LU,1303) 
     *              P_DMIN(FLXC(ISEC))+P_DMAX(FLXC(ISEC)),
     *                              P_DMIN(VOLNEGC(ISEC)),
     *                              P_DMAX(VOLPOSC(ISEC))     
        IF(LNG.EQ.2) WRITE(LU,1304) 
     *              P_DMIN(FLXC(ISEC))+P_DMAX(FLXC(ISEC)),
     *                              P_DMIN(VOLNEGC(ISEC)),
     *                              P_DMAX(VOLPOSC(ISEC))
      ENDIF
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
      ENDIF ! ncsize
C
C-----------------------------------------------------------------------
C CHAIN allocated, i.e. serial or parallel case from SECTIONS INPUT FILE
C       we can apply co-ordinates instead and/or names of sections

      ELSE ! .not.classic 
C
        DO ISEC = 1,NSEC
C
          IF(LNG.EQ.1) WRITE(LU,230) ISEC,TRIM(CHAIN(ISEC)%DESCR),
     *                 P_DSUM(FLX(ISEC)),P_DSUM(VOLNEG(ISEC)),
     *                                   P_DSUM(VOLPOS(ISEC))
          IF(LNG.EQ.2) WRITE(LU,231) ISEC,TRIM(CHAIN(ISEC)%DESCR),
     *                 P_DSUM(FLX(ISEC)),P_DSUM(VOLNEG(ISEC)),
     *                                   P_DSUM(VOLPOS(ISEC))
230       FORMAT(1X,/,1X,'SECTION DE CONTROLE ',1I2,
     *               ' (NOM ',A,')',//,5X,
     *               'DEBIT :                     ',G16.7,/,5X,
     *               'CUMUL DES DEBITS NEGATIFS : ',G16.7,/,5X,
     *               'CUMUL DES DEBITS POSITIFS : ',G16.7)
231       FORMAT(1X,/,1X,'CONTROL SECTION NUMBER ',1I2,
     *               ' (NAME ',A,')',//,5X,
     *               'DISCHARGE:                 ',G16.7,/,5X,
     *               'CUMULATED NEGATIVE VOLUME: ',G16.7,/,5X,
     *               'CUMULATED POSITIVE VOLUME: ',G16.7)
          IF(SUSP) THEN
            IF(LNG.EQ.1) WRITE(LU,2301) 
     *              P_DSUM(FLXS(ISEC)),P_DSUM(VOLNEGS(ISEC)),
     *                                 P_DSUM(VOLPOSS(ISEC))
            IF(LNG.EQ.2) WRITE(LU,2302) 
     *              P_DSUM(FLXS(ISEC)),P_DSUM(VOLNEGS(ISEC)),
     *                                 P_DSUM(VOLPOSS(ISEC))
2301        FORMAT(5X,'DEBIT EN SUSPENSION :       ',G16.7,/,5X,
     *            'CUMUL DES DEBITS NEGATIFS : ',G16.7,/,5X,
     *            'CUMUL DES DEBITS POSITIFS : ',G16.7)
2302        FORMAT(5X,'DISCHARGE IN SUSPENSION:   ',G16.7,/,5X,
     *            'CUMULATED NEGATIVE VOLUME: ',G16.7,/,5X,
     *            'CUMULATED POSITIVE VOLUME: ',G16.7)
          ENDIF
          IF(CHARR) THEN
            IF(LNG.EQ.1) WRITE(LU,2303) 
     *              P_DSUM(FLXC(ISEC)),P_DSUM(VOLNEGC(ISEC)),
     *                                 P_DSUM(VOLPOSC(ISEC))
            IF(LNG.EQ.2) WRITE(LU,2304) 
     *              P_DSUM(FLXC(ISEC)),P_DSUM(VOLNEGC(ISEC)),
     *                                 P_DSUM(VOLPOSC(ISEC))
2303        FORMAT(5X,'DEBIT EN CHARRIAGE :        ',G16.7,/,5X,
     *            'CUMUL DES DEBITS NEGATIFS : ',G16.7,/,5X,
     *            'CUMUL DES DEBITS POSITIFS : ',G16.7)
2304        FORMAT(5X,'BEDLOAD DISCHARGE:         ',G16.7,/,5X,
     *            'CUMULATED NEGATIVE VOLUME: ',G16.7,/,5X,
     *            'CUMULATED POSITIVE VOLUME: ',G16.7)
          ENDIF
C
        ENDDO
C
      ENDIF ! if classic 
C
      ENDIF ! if info 
C !jaj ####
C-----------------------------------------------------------------------
C master writes a nice sections output file, the header only once
C NOTE: programmed for the bedl load discharge only (if charr implied)
C
      IF ( (.NOT.CLASSIC) .AND. CHARR .AND. 
     &     (TRIM(SIS_FILES(SISSEO)%NAME).NE.'') ) THEN 
        IF (INIT) THEN 
          INIT=.FALSE.
          IF ((NCSIZE.GT.1 .AND. IPID.EQ.0).OR.(NCSIZE.LE.1)) THEN
            NSEO=SIS_FILES(SISSEO)%LU
            WRITE(NSEO,*) 'TITLE = "Bedload discharges for ',
     &                     TRIM(TITCA),'"' 
            WRITE(NSEO,*) 'VARIABLES = t', 
     &           (' '//TRIM(CHAIN(ISEC)%DESCR), ISEC=1,NSEC)
          ENDIF 
          IF (NCSIZE.GT.1) THEN 
            ALLOCATE (WORK(NSEC), STAT=ERR)
            IF (ERR.NE.0) THEN
              WRITE(LU,*) 'FLUXPR_SISYPHE: error allocating WORK:',ERR
              CALL PLANTE(1)
              STOP
            ENDIF
          ENDIF
        ENDIF 
        ! deadlock with write and p_dsum in an implied WRITE loop 
        ! because it is only master to write the message...
        IF (NCSIZE.GT.1) THEN 
          DO ISEC=1,NSEC
            WORK(ISEC)=P_DSUM(FLXC(ISEC))
          END DO
          IF (IPID.EQ.0) 
     &      WRITE (NSEO, FMT=FMTZON) TPS, (WORK(ISEC), ISEC=1,NSEC)
        ELSE 
          WRITE (NSEO, FMT=FMTZON) TPS, (FLXC(ISEC), ISEC=1,NSEC)
        ENDIF 
      ENDIF
C
C-----------------------------------------------------------------------
C 
      RETURN
      END SUBROUTINE FLUXPR_SISYPHE
