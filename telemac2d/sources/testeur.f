C                              ******************
                               SUBROUTINE TESTEUR
C                              ******************
C
     *(NS,NSEG,NPTFR,NUBO,DT,NBOR,NORDRE,AIRS,AIRST,HSTOK,
     * HCSTOK,FLUXT,FLUXTEMP,FLUHBOR,FLUHBTEMP,LOGFR,TEST,NTRAC)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.8    05/09/07                      INRIA
C
C***********************************************************************
C        
C   TEST POUR LA POSITIVITE DU TRACEUR
C
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
      USE BIEF
C	
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NS,NSEG,NPTFR,NORDRE,NTRAC
      INTEGER, INTENT(IN)             :: NUBO(2,NSEG),NBOR(NPTFR)
      INTEGER, INTENT(IN)             :: LOGFR(*)
      DOUBLE PRECISION, INTENT(IN)    :: AIRS(NS),AIRST(2,*)
      DOUBLE PRECISION, INTENT(IN)    :: DT
      DOUBLE PRECISION, INTENT(INOUT) :: TEST
      DOUBLE PRECISION, INTENT(IN)    :: HSTOK(*)
      DOUBLE PRECISION, INTENT(IN)    :: HCSTOK(2,*)
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: FLUHBOR,FLUHBTEMP,FLUXT
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: FLUXTEMP
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IS,K,NSG,NUBO1,NUBO2,ERR,ITRAC
C   
      DOUBLE PRECISION AUX 
C
C------------------------------------------------------------------------
C
C     EX TABLEAU AUTOMATIQUE !!!
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: FLUXTEST(:)
C
      LOGICAL DEJA
      DATA DEJA/.FALSE./
C
      IF(.NOT.DEJA) THEN
        ALLOCATE(FLUXTEST(NS),STAT=ERR)
        IF(ERR.NE.0) GO TO 1001
        GO TO 1002
1001    CONTINUE
        IF(LNG.EQ.1) WRITE(LU,1000) ERR
        IF(LNG.EQ.2) WRITE(LU,2000) ERR
1000    FORMAT(1X,'TESTEUR : ERREUR A L''ALLOCATION DE MEMOIRE : ',/,1X,
     *         'CODE D''ERREUR : ',1I6)
2000    FORMAT(1X,'TESTEUR: ERROR DURING ALLOCATION OF MEMORY: ',/,1X,
     *        'ERROR CODE: ',1I6)
        CALL PLANTE(1)
        STOP
1002    CONTINUE
        DEJA=.TRUE.
      ENDIF
C
C------------------------------------------------------------------------
C
C     TEST SUR TOUS LES TRACEURS, IL SUFFIT D'UN TEST<0 POUR SORTIR
C
      DO ITRAC=1,NTRAC
C
      DO IS=1,NS
        FLUXTEST(IS) = 0.D0
      ENDDO
C
C   TEST POUR  ORDRE 1
C
      IF(NORDRE.EQ.1) THEN
C
      DO NSG=1,NSEG
         NUBO1=NUBO(1,NSG)
         NUBO2=NUBO(2,NSG)
         AUX = FLUXT%ADR(ITRAC)%P%R(NSG)+DT*FLUXTEMP%ADR(ITRAC)%P%R(NSG)
         IF(AUX.GE.0.D0) THEN
           FLUXTEST(NUBO1) = FLUXTEST(NUBO1) + AUX
         ELSE
           FLUXTEST(NUBO2) = FLUXTEST(NUBO2) - AUX
         ENDIF
      ENDDO
      DO K=1,NPTFR
         IS =NBOR(K)
         AUX = FLUHBOR%ADR(ITRAC)%P%R(K)+DT*FLUHBTEMP%ADR(ITRAC)%P%R(K)
         IF(AUX.GE.0.D0) FLUXTEST(IS)=FLUXTEST(IS)+AUX
      ENDDO
      DO IS=1,NS
         TEST=AIRS(IS)*HSTOK(IS)-FLUXTEST(IS)
         IF(TEST.LT.0.D0) RETURN 
      ENDDO
C
      ELSE
C
C   TEST POUR ORDRE 2
C
         DO NSG=1,NSEG 
         NUBO1=NUBO(1,NSG)
         NUBO2=NUBO(2,NSG)
C
         AUX = FLUXT%ADR(ITRAC)%P%R(NSG)+DT*FLUXTEMP%ADR(ITRAC)%P%R(NSG)
C
         IF(AUX.GE.0.D0) THEN
C
         TEST=AIRST(1,NSG)*HCSTOK(1,NSG)-AUX 
         IF(LOGFR(NUBO1).NE.0.AND.
     &     AIRST(1,NSG)*HCSTOK(1,NSG).GT.0.D0) THEN
         FLUXTEST(NUBO1)= MAX(FLUXTEST(NUBO1),
     &     AUX/(AIRST(1,NSG)*HCSTOK(1,NSG)))
         ENDIF
C
         ELSE
C
         TEST=AIRST(2,NSG)*HCSTOK(2,NSG)+AUX          
         IF(LOGFR(NUBO2).NE.0.AND.
     &     AIRST(2,NSG)*HCSTOK(2,NSG).GT.0.D0) THEN
         FLUXTEST(NUBO2)= MAX(FLUXTEST(NUBO2),
     &     -AUX/(AIRST(2,NSG)*HCSTOK(2,NSG)))
         ENDIF
C
         ENDIF
         IF(TEST.LT.0.D0) RETURN 
      ENDDO
C
C   TEST POUR LES NOEUDS FRONTIERE
C
        DO K=1,NPTFR
         IS =NBOR(K)
         AUX = FLUHBOR%ADR(ITRAC)%P%R(K)+DT*FLUHBTEMP%ADR(ITRAC)%P%R(K)
         IF(AUX.GE.0.D0) THEN 
           TEST=AIRS(IS)*HSTOK(IS)-AUX
           IF(AIRS(IS)*HSTOK(IS).GT.0.D0.AND.
     &     (1.D0-FLUXTEST(IS)).LT. AUX/(AIRS(IS)*HSTOK(IS)))
     &     TEST =-1.D0
         ENDIF
         IF(TEST.LT.0.D0) RETURN 
C
       ENDDO
C
      ENDIF
C
C     DO ITRAC=1,NTRAC
      ENDDO
C
C------------------------------------------------------------------------
C
      RETURN
      END
