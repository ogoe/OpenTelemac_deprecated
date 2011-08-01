C                       *****************
                        SUBROUTINE FLUXPR
C                       *****************
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
      DOUBLE PRECISION, INTENT(IN) :: FLX(NSEC),TPS
      DOUBLE PRECISION, INTENT(IN) :: VOLNEG(NSEC),VOLPOS(NSEC)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION P_DMAX,P_DMIN
      INTEGER                        P_IMIN
      EXTERNAL         P_DMAX,P_DMIN,P_IMIN
C
      INTEGER ISEC,II
C
C-----------------------------------------------------------------------
C 
      IF(INFO) THEN
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
      ENDIF
C
C-----------------------------------------------------------------------
C 
      RETURN
      END 
