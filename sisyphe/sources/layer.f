C                       ****************
                        SUBROUTINE LAYER 
C                       ****************
C
     *(ZFCL_W,NLAYER,ZR,ZF,ESTRAT,ELAY,MASBAS,ACLADM,NSICLA,NPOIN,
     * ELAY0,VOLTOT,ES,AVAIL,CONST_ALAYER,DTS,ESTRATNEW,NLAYNEW)
C
C***********************************************************************
C SISYPHE VERSION 6.0                       MATTHIEU GONZALES DE LINARES
C                                                2002       
C
C 16/09/2009 : JMH  ->  AVAIL(NPOIN,10,NSICLA) 
C 10/05/2010 : JMH  ->  CASE WITH DEPOSIT REWRITTEN, TESTS CHANGED
C                       OTHER PARTS REWRITTEN AND/OR OPTIMISED
C                       CLEAN STOP IN PARALLEL IF PROBLEM
C 
C COPYRIGHT EDF   
C***********************************************************************
C                                   
C FONCTION : CALCULATION OF AVAIL FOR EACH CLASS AND EACH LAYER                                     
C            NEW STRATUM THICKNESS ESTRAT 
C
C            ACTIVE LAYER IS LAYER 1, IT IS KEPT AT A PRESCRIBED HEIGHT
C            OF ELAY0.
C
C            STRATUM IS LAYER 2 OF HEIGHT ESTRAT
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |    ZFCL_W      | -->| EVOLUTION FOR EACH SEDIMENT CLASS
C |    NLAYER      |<-- | NUMBER OF LAYER FOR EACH POINT
C |    ELAY        |<-- | ACTIVE LAYER THICKNESS FOR EACH POINT
C |    ESTRAT      |<-- | ACTIVE STRATUM THICKNESS FOR EACH POINT
C |    AVAIL       |<-- | SEDIMENT FRACTION FOR EACH LAYER, CLASS, POINT
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C PROGRAMME APPELANT : SISYPHE
C PROGRAMMES APPELES : 
C***********************************************************************
C
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE (BIEF_OBJ),  INTENT(IN)    :: ZFCL_W,ZR,ZF
      TYPE (BIEF_OBJ),  INTENT(IN)    :: MASBAS,ACLADM
      INTEGER,          INTENT(IN)    :: NSICLA,NPOIN
      DOUBLE PRECISION, INTENT(IN)    :: DTS
      LOGICAL,          INTENT(IN)    :: CONST_ALAYER 
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: NLAYER,ESTRAT,ELAY
      DOUBLE PRECISION, INTENT(INOUT) :: ELAY0
      DOUBLE PRECISION, INTENT(INOUT) :: ES(NPOIN,10)
      DOUBLE PRECISION, INTENT(INOUT) :: AVAIL(NPOIN,10,NSICLA)
      DOUBLE PRECISION, INTENT(INOUT) :: VOLTOT(10),ESTRATNEW(NPOIN)
      INTEGER         , INTENT(INOUT) :: NLAYNEW(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION P_DSUM
      EXTERNAL         P_DSUM
      INTEGER  P_ISUM
      EXTERNAL P_ISUM
C
C-----------------------------------------------------------------------
C
      INTEGER I,J,K,ARRET,ARRET2
      DOUBLE PRECISION EVOL,HEIGH,TEST1,TEST2,AEVOL,AUX
C
C-----------------------------------------------------------------------
C
C     FOR CHECKING FRACTIONS IN THE RANGE [-ZERO,1+ZERO]
C
      DOUBLE PRECISION ZERO
      DATA             ZERO/1.D-10/
C
C-----------------------------------------------------------------------
C
      ARRET=0
C
C-----------------------------------------------------------------------
C      
      DO J=1,NPOIN
C
        IF(.NOT.CONST_ALAYER) ELAY0 = 3.D0 * ACLADM%R(J)
C
        NLAYNEW(J) = NLAYER%I(J)
C
C QUESTION JMH, EVOLUTION HAS BEEN COMPUTED BEFORE IN ARRAY E, WHY NOT
C               EVOL=E(J) ?????
C               ELAY(J) = ES(J,1) WHY IS IT AN EXTRA ARRAY ??
C
C     
        EVOL  = 0.D0
        HEIGH = ZF%R(J)-ZR%R(J)
!       HERE ELAY.NE.HEIGH BECAUSE ELAY IS THE ACTIVE LAYER THICKNESS         
        DO I=1,NSICLA
          EVOL = EVOL + ZFCL_W%ADR(I)%P%R(J)
        ENDDO
C
        IF(NLAYER%I(J).GT.1) THEN
C
          IF(EVOL.GE.0.D0) THEN
C
C           DEPOSIT
C
C           NEW HEIGHT OF LAYER 2 (IT RECEIVES EVOL TO KEEP LAYER 1 CONSTANT)
            ESTRATNEW(J) = ESTRAT%R(J) + EVOL
C
            DO I=1,NSICLA
C             JMH 28/04/2010. THE OLD IMPLEMENTATION CONSISTED OF FIRST
C             GIVING EVOL TO LAYER 2, WITH OLD AVAIL, THEN OF RECEIVING
C             THE DEPOSIT, BUT IF A CLASS DISAPPEARS IN LAYER 1,
C             IT IS NOT POSSIBLE TO GIVE IT FIRST TO LAYER 2, SO NEW
C             FRACTIONS MUST BE COMPUTED BEFORE GIVING EVOL TO LAYER 2 
C             THEN I SEE NO DIFFERENCE BETWEEN EVOL<ELAY0 AND EVOL>ELAY0 
C             SO BOTH ARE TREATED BELOW, UNLIKE VERSION 5.9.           
C
C             1) LAYER 1 RECEIVES ZFCL_W OF CLASS I, WE COMPUTE THE
C                PROVISIONAL NEW AVAIL(J,1,I) IN AUX.
              AUX=(AVAIL(J,1,I)*ELAY0+ZFCL_W%ADR(I)%P%R(J))/
     *                        (ELAY0+EVOL) 
C
C             2) LAYER 2 RECEIVES AUX*EVOL OF CLASS I (AUX MAY BE 0 HERE)                 
              AVAIL(J,2,I)=(AUX*EVOL+AVAIL(J,2,I)*ESTRAT%R(J))/
     *                               ESTRATNEW(J)
C
C             3) SEEN FROM LAYER 1, AUX*EVOL OF CLASS I HAS BEEN GIVEN
C                AND THE NEW LAYER THICKNESS IS ELAY0, HENCE THE NEW AVAIL
              AVAIL(J,1,I)=( AVAIL(J,1,I)*ELAY0-AUX*EVOL
     &                          +ZFCL_W%ADR(I)%P%R(J) )/ELAY0
!
! note CV: can be replaced by 
!              AVAIL(J,1,I)= AUX
!
!             OLD (AND WRONG) FORMULA (IN IT -AVAIL*EVOL SHOULD BE -AUX*EVOL)
!             AVAIL(J,1,I)=( AVAIL(J,1,I)*(ELAY0-EVOL)
!    &                       +ZFCL_W%ADR(I)%P%R(J)     )/ELAY0
!
C             IF(AVAIL(J,1,I).GT.1.D0+ZERO.OR.
C    *           AVAIL(J,1,I).LT.-ZERO) THEN
C               WRITE(LU,*) 'ERROR IN LAYER CASE 1'
C               STOP
C             ENDIF 
            ENDDO
C           NEW HEIGHT OF LAYER 1   
            ELAY%R(J) = ELAY0
C          
          ELSEIF(EVOL.GT.-ELAY0) THEN
C CV: je ne suis pas d'accord avec le commentaire ci-desssous 
C   On est dans le cas : -ELAY0 <EVOL<0 donc ELAY0> -EVOL>0
C   l'épaisseur de la premiere couche est donc suffisante 
C 
C           EROSION GREATER THAN LAYER 1, WE HAVE TO DESTROY A STRATUM
C               
            IF(-EVOL.GT.ESTRAT%R(J)) THEN
!
!  CV : dans le cas habituel, la couche 2 est très large et deux couches sont en général suffisantes
!       ici il y a destruction de la couche 2
!  
!             USUAL CASE (NOTE JMH : WHY NOT .GE.2 ?
! 
              IF(NLAYER%I(J).GT.2) THEN
!
                DO I=1,NSICLA
                  AVAIL(J,1,I) = (  AVAIL(J,1,I)*ELAY0
     &                            + ZFCL_W%ADR(I)%P%R(J)
     &                            + AVAIL(J,2,I)*ESTRAT%R(J)
     &                            - AVAIL(J,3,I)*(EVOL+ESTRAT%R(J))
     &                           )/ ELAY0
C                 IF(AVAIL(J,1,I).GT.1.D0+ZERO.OR.
C    &               AVAIL(J,1,I).LT.-ZERO) THEN
C                   WRITE(LU,*) 'ERROR IN LAYER CASE 2'
C                   STOP
C                 ENDIF
                  AVAIL(J,2,I) = AVAIL(J,3,I)
                  DO K=3,MIN(9,NLAYER%I(J))
                    AVAIL(J,K,I) = AVAIL(J,K+1,I)
                  ENDDO
                ENDDO
                ELAY%R(J) = ELAY0
                NLAYNEW(J) = NLAYER%I(J) - 1
                ESTRATNEW(J) = ESTRAT%R(J) + EVOL + ES(J,3)
                DO K=3,MIN(9,NLAYER%I(J))
                  ES(J,K) = ES(J,K+1)
                ENDDO
!
!             ONLY ONE LAYER LEFT (NOTE JMH : 1 OR 2 ?)
!             
              ELSE
                DO I=1,NSICLA
                  IF(HEIGH.GT.0.D0) THEN            
                    AVAIL(J,1,I) = (  AVAIL(J,1,I)*ELAY0
     &                              + ZFCL_W%ADR(I)%P%R(J)
     &                              + ESTRAT%R(J)*AVAIL(J,2,I)
     &                              )/(ELAY%R(J)+EVOL+ESTRAT%R(J)) 
C                   IF(AVAIL(J,1,I).GT.1.D0+ZERO.OR.
C    &                 AVAIL(J,1,I).LT.-ZERO) THEN
C                     WRITE(LU,*) 'ERROR IN LAYER CASE 3'
C                     STOP
C                   ENDIF
                  ELSE
                    AVAIL(J,1,I) = 0.D0
                  ENDIF
                  AVAIL(J,2,I) = 0.D0 
                ENDDO
                NLAYNEW(J) = NLAYER%I(J) - 1
                ELAY%R(J) = HEIGH
                ESTRATNEW(J) = 0.D0                   
              ENDIF
!
!           ONLY LAYER 1 ERODED
!
            ELSE
              DO I=1,NSICLA
                AVAIL(J,1,I) = (  AVAIL(J,1,I) * ELAY0
     &                          + ZFCL_W%ADR(I)%P%R(J)
     &                          - EVOL*AVAIL(J,2,I)    )/ELAY0
C               IF(AVAIL(J,1,I).GT.1.D0+ZERO.OR.
C    &             AVAIL(J,1,I).LT.-ZERO) THEN
C                 WRITE(LU,*) 'ERROR IN LAYER CASE 4'
C                 STOP
C               ENDIF
              ENDDO
              ELAY%R(J) = ELAY0
              ESTRATNEW(J) = ESTRAT%R(J) + EVOL
            ENDIF
!
          ELSE
C
C           TEST D'ARRET SI EROSION SUPERIEURE A 
C           EPAISSEUR DE LA COUCHE ACTIVE !
C 
            IF(LNG.EQ.1) THEN
              WRITE(LU,*) 'EROSION TROP FORTE AU NOEUD J=',J
              WRITE(LU,*) 'DIMINUER DT OU AUGMENTER ELAY0'
            ENDIF
            IF(LNG.EQ.2) THEN
              WRITE(LU,*) 'TOO MUCH EROSION AT POINT J=',J
              WRITE(LU,*) 'DECREASE DT OR INCREASE ELAY0'              
            ENDIF
            WRITE(LU,*) 'EVOL=', EVOL, 'ELAY0=',ELAY0
            CALL PLANTE(1)
            STOP
C
!         END OF TESTS ON EVOL
!
          ENDIF 
!
! THERE WAS ONLY ONE LAYER
! ------------------------
!
        ELSE
!
!         IT IS NOW BIG ENOUGH TO MAKE TWO LAYERS
!
          IF(HEIGH.GT.ELAY0) THEN
            NLAYNEW(J) = 2
            ESTRATNEW(J) = HEIGH - ELAY0
            ELAY%R(J) = ELAY0
            DO I=1,NSICLA
              AVAIL(J,2,I) = AVAIL(J,1,I)
              AVAIL(J,1,I) = (AVAIL(J,1,I) * (ELAY0-EVOL)
     &                     + ZFCL_W%ADR(I)%P%R(J) )/ELAY0
C             IF(AVAIL(J,1,I).GT.1.D0+ZERO.OR.
C    &           AVAIL(J,1,I).LT.-ZERO) THEN
C                WRITE(LU,*) 'ERROR IN LAYER CASE 5'
C                STOP
C             ENDIF
            ENDDO        
!
! IF THERE REMAINS ONLY ONE LAYER
! -------------------------------
!
          ELSE
!           NOTE JMH: THE TRICKIEST PART...
!           THE PROBLEM OF 0/0 CREATED BY THE CHOICE OF AVAIL
!           AS MAIN VARIABLE...
            IF(ELAY%R(J)+EVOL.GT.1.D-15) THEN
              DO I=1,NSICLA 
C               AUX=AVAIL(J,1,I)   
                AVAIL(J,1,I) = (AVAIL(J,1,I)*ELAY%R(J)+
     &                          ZFCL_W%ADR(I)%P%R(J))
     &                          / (ELAY%R(J)+EVOL)
C               IF(AVAIL(J,1,I).GT.1.D0+ZERO.OR.
C    &            AVAIL(J,1,I).LT.-ZERO) THEN
C                 WRITE(LU,*) 'ERROR IN LAYER CASE 6'
C                 WRITE(LU,*) 'INITIAL AVAIL=',AUX
C                 WRITE(LU,*) 'J=',J,' CLASS ',I
C                 WRITE(LU,*) 'EVOL=',EVOL,' ELAY=',ELAY%R(J)
C                 WRITE(LU,*) 'ZFCL=',ZFCL_W%ADR(I)%P%R(J)
C                 WRITE(LU,*) 'DENOMINATOR=',ELAY%R(J)+EVOL
C                 WRITE(LU,*) 'NUMERATOR=',AUX*ELAY%R(J)+
C    &                                     ZFCL_W%ADR(I)%P%R(J)
C               ENDIF                
                AVAIL(J,2,I) = 0.D0
              ENDDO
              IF(ELAY%R(J)+EVOL.LT.1.D-7) THEN
!               PLAYING WITH ZEROES, RISK OF SUM NOT EQUAL TO 1
!               ONLY BECAUSE OF TRUNCATION ERRORS, WE NORMALIZE
                TEST1=0.D0
                DO I=1,NSICLA 
                  AVAIL(J,1,I)=MAX(0.D0,MIN(1.D0,AVAIL(J,1,I)))
                  TEST1=TEST1+AVAIL(J,1,I)
                ENDDO
                IF((TEST1-1.D0)**2.GT.ZERO) THEN
                  DO I=1,NSICLA 
                    AVAIL(J,1,I)=AVAIL(J,1,I)/MAX(TEST1,1.D-21)
                  ENDDO
                ENDIF
              ENDIF
            ELSE
              DO I=1,NSICLA 
                AVAIL(J,1,I) = 0.D0
                AVAIL(J,2,I) = 0.D0
              ENDDO
            ENDIF                        
            ELAY%R(J) = HEIGH
            ESTRATNEW(J) = 0.D0
            NLAYNEW(J) = 1
          ENDIF
        ENDIF
!
      NLAYER%I(J) = NLAYNEW(J)
      ESTRAT%R(J) = ESTRATNEW(J)
      ES(J,1) = ELAY%R(J)
      IF(NLAYER%I(J).GT.1) ES(J,2) = ESTRAT%R(J)
!
      TEST1 = 0.D0
      TEST2 = 0.D0          
!
      DO I=1,NSICLA
        DO K=1,NLAYER%I(J)
!         TESTING THAT AVAIL IS IN THE RANGE (-ZERO,1+ZERO)  
          IF(AVAIL(J,K,I).GT.1.D0+ZERO.OR.AVAIL(J,K,I).LT.-ZERO) THEN 
            WRITE(LU,*) 'ERROR ON FRACTIONS'
            WRITE(LU,*) 'LAYER ',K,' CLASS ',I,' POINT ',J
            IF(AVAIL(J,K,I).LT.0.D0) THEN
              WRITE(LU,*) 'AVAIL=' ,AVAIL(J,K,I)
            ELSE
              WRITE(LU,*) 'AVAIL-1=' ,AVAIL(J,K,I)-1.D0
            ENDIF
            WRITE(LU,*) 'ZFCL=',ZFCL_W%ADR(I)%P%R(J)
            WRITE(LU,*) 'EVOL=',EVOL,' ELAY=',ELAY%R(J)
            ARRET=1
          ELSE 
!           ONCE CHECKED THAT WE HAVE ONLY TRUNCATION ERRORS, CLIPPING
            AVAIL(J,1,I)=MAX(AVAIL(J,1,I),0.D0)
            AVAIL(J,1,I)=MIN(AVAIL(J,1,I),1.D0)    
          ENDIF  
        ENDDO
        TEST1 = TEST1 + AVAIL(J,1,I)
        TEST2 = TEST2 + AVAIL(J,2,I)
      ENDDO
! 
!     TESTING THAT SUM OF AVAIL IS 1 FOR FIRST 2 LAYERS
!
      IF(TEST1.GT.ZERO.AND.(TEST1-1.D0)**2>ZERO) THEN
        WRITE(LU,*) ' PROBLEM IN LAYER: J,TEST1',J,TEST1
        WRITE(LU,*) ' IN LAYER 1 SUM OF FRACTIONS NOT 1'
        ARRET=1 
      ENDIF
      IF(TEST2.GT.ZERO.AND.(TEST2-1.D0)**2>ZERO) THEN
        WRITE(LU,*) ' PROBLEM IN LAYER: J,TEST2',J,TEST2
        WRITE(LU,*) ' IN LAYER 2 SUM OF FRACTIONS IS NOT 1'
        ARRET=1 
      ENDIF
!
!     END OF LOOP ON ALL POINTS
!      
      ENDDO    
C
C     CALCUL DU VOLUME TOTAL DE SEDIMENTS DANS LE DOMAINE
C          
      DO I = 1, NSICLA  
        VOLTOT(I) = 0.D0
      ENDDO
      DO I=1,NSICLA
        DO J=1,NPOIN
          DO K=1,NLAYER%I(J)        
            VOLTOT(I) = VOLTOT(I) + ES(J,K)*AVAIL(J,K,I)*MASBAS%R(J) 
          ENDDO
        ENDDO
      ENDDO
!
      IF(NCSIZE.GT.1) THEN
        DO I=1,NSICLA
          VOLTOT(I) = P_DSUM(VOLTOT(I))
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
!     CLEAN STOP FOR ALL PROCESSORS IF PROBLEM
!
      ARRET2=ARRET
      IF(NCSIZE.GT.1) ARRET2=P_ISUM(ARRET)     
      IF(ARRET2.GT.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'ARRET APRES ERREUR DANS LAYER'
        IF(LNG.EQ.2) WRITE(LU,*) 'STOP AFTER AN ERROR IN LAYER'
        IF(ARRET.EQ.0) THEN
          IF(LNG.EQ.1) WRITE(LU,*) 'DANS ',ARRET2,' PROCESSEUR(S)'
          IF(LNG.EQ.2) WRITE(LU,*) 'IN ',ARRET2,' PROCESSOR(S)'
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
!
!-----------------------------------------------------------------------
! 
      RETURN
      END
