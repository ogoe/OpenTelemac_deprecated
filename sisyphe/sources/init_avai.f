C                         ********************
                          SUBROUTINE INIT_AVAI 
C                         ********************
C
C
C***********************************************************************
C SISYPHE VERSION 6.0              
C
C                            BUI MINH DUC  2002  
C                            initial fraction distribution
C                            for non uniform bed materials 
C
C                            Matthieu GONZALES DE LINARES 2002/2003
C
C                                                
C COPYRIGHT EDF-BAW-IFH   
C***********************************************************************
C
C     FONCTION  : INITIAL FRACTION DISTRIBUTION AND LAYER THICKNESS
C
C     PHILOSOPHIE : INIT_COMPO PERMET DE DEFINIR LA STRATIFICATION
C     CORRESPONDANT A LA VARIATION DE COMPOSITION INITIALE DU FOND,
C     AINSI NCOUCHES CORRESPOND AU NOMBRE DE COUCHES INITIALES REELLES
C                   INIT_AVAI CORRIGE ET COMPLETE CETTE STRATIFICATION
C     EN CAS DE PROBLEMES, ET RAJOUTE LA COUCHE ACTIVE, AINSI NLAYER
C     CORRESPOND AU NOMBRE DE COUCHES UTILISEES DANS LE CALCUL
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |    ZF          | -->|  BOTTOM
C |    ZR          | -->|  RIGID BOTTOM
C |    ELAY0       | -->|  WANTED ACTIVE LAYER THICKNESS
C |    NLAYER      |<-- |  NUMBER OF LAYER FOR EACH POINT
C |    ELAY        |<-- |  ACTIVE LAYER THICKNESS FOR EACH POINT
C |    ESTRAT      |<-- |  ACTIVE STRATUM THICKNESS FOR EACH POINT
C |    AVAIL       |<-- |  SEDIMENT FRACTION FOR EACH LAYER, CLASS, NODE
C |    ES          |<-- |  THICKNESS FOR EACH LAYER AND NODE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C PROGRAMME APPELANT : SISYPHE
C PROGRAMMES APPELES : INIT_COMPO
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_SISYPHE
      USE INTERFACE_SISYPHE,EX_INIT_AVAI=> INIT_AVAI
C
      IMPLICIT NONE
C      
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION P_DSUM
      EXTERNAL         P_DSUM
C
      INTEGER I,J,K,NMAXI
C     
      DOUBLE PRECISION HAUTSED,TEST1
C
C-----------------------------------------------------------------------
C
      NMAXI = 0
C
C     USER SET THE INITIAL NUMBER OF LAYERS, THEIR THICKNESS 
C     AND THEIR COMPOSITION
C
C     NOTE: WHEN COMPUTATION CONTINUED INIT_COMPO MUST NOT
C           CHANGE ES AND AVAIL
C
      IF(DEBU) THEN
!       ADDED BY JMH 30/06/2010
        DO J=1,NPOIN
          I=1
          DO K=2,NOMBLAY
            IF(ES(J,K).GT.0.D0) I = I + 1
          ENDDO
          NLAYER%I(J)=I
!         CHECKING ALL LAYERS AND CORRECTING AVAIL
!         DUE TO POSSIBLE SHIFT OF SINGLE PRECISION STORAGE
          DO K=1,NLAYER%I(J)
            TEST1=0.D0
            DO I=1,NSICLA
              TEST1=TEST1+AVAIL(J,K,I)
            ENDDO
            IF(TEST1.GT.1.D-10.AND.(TEST1-1.D0)**2.GT.1.D-10) THEN
              DO I=1,NSICLA
                AVAIL(J,K,I)=AVAIL(J,K,I)/TEST1
              ENDDO
            ENDIF 
          ENDDO          
        ENDDO
      ELSE
        CALL INIT_COMPO(IT1%I)     
C
        DO 45 J=1,NPOIN
C  
C       10 IS THE MAXIMUM NUMBER OF LAYERS ALLOWED
        NLAYER%I(J) = IT1%I(J)
        IF(NLAYER%I(J).GT.10) THEN
          IF(LNG.EQ.1) WRITE(LU,1800)         
          IF(LNG.EQ.2) WRITE(LU,1815)
          CALL PLANTE(1)
          STOP
        ENDIF
C         
C       THE HEIGHT OF SEDIMENT (SUM OF ES) MUST NOT BE MORE THAN ZF-ZR
C       IF SO, THE HEIGHT OF THE LAST LAYER IS REDUCED
C       IF THERE ARE LAYERS UNDER ZR, THEY ARE NOT TAKEN INTO ACCOUNT
        HAUTSED = 0.D0
        DO K=1,IT1%I(J)
        IF(HAUTSED + ES(J,K) .GE. ZF%R(J) - ZR%R(J)) THEN
          ES(J,K) = ZF%R(J) - ZR%R(J) -  HAUTSED 
          NLAYER%I(J) = K
          HAUTSED = HAUTSED + ES(J,K)
          GOTO 144
        ENDIF
        HAUTSED = HAUTSED + ES(J,K)
        ENDDO
C
144     CONTINUE
C
C       THE THICKNESS OF THE LAST LAYER IS ENLARGED SO THAT
C       THE HEIGHT OF SEDIMENT (SUM OF ES) IS EQUAL TO ZF-ZR
        IF(HAUTSED.LT.ZF%R(J)-ZR%R(J)) THEN
          ES(J,NLAYER%I(J))=ES(J,NLAYER%I(J))+ZF%R(J)-ZR%R(J)-HAUTSED
        ENDIF
C
        IF(NLAYER%I(J).GT.1) THEN
C         On suppose que ELAY0 est plus petit que la premiere strate.
C         A rajouter : le cas ou ELAY0 est plus grand
          IF(ELAY0.GT.ES(J,1)) THEN
            IF(LNG.EQ.1) THEN
              WRITE(LU,*) 'COUCHE ACTIVE TROP GROSSE/STRATIFICATION'
            ENDIF
            IF(LNG.EQ.2) THEN
              WRITE(LU,*) ' ACTIVE LAYER TOO BIG/STRATIFICATION'
            ENDIF
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
C       ON SEPARE LA PREMIERE STRATE EN ACTIVE LAYER + ACTIVE STRATUM   
        IF(ELAY0.LT.ES(J,1)) THEN 
          NLAYER%I(J) = NLAYER%I(J) + 1
          IF(NLAYER%I(J).GT.10) THEN
            IF(LNG.EQ.1) WRITE(LU,1800)         
            IF(LNG.EQ.2) WRITE(LU,1815)
            CALL PLANTE(1)
            STOP
          ENDIF
C         IL FAUT DECALER LES INDICES POUR ES ET AVAIL           
          IF(NLAYER%I(J).GT.2) THEN 
            DO K=NLAYER%I(J),3,-1
              ES(J,K) = ES(J,K-1)
            ENDDO
          ENDIF
          ES(J,2) = ES(J,1) - ELAY0
          ES(J,1) = ELAY0
          DO I=1,NSICLA
            DO K=NLAYER%I(J),2,-1
              AVAIL(J,K,I) = AVAIL(J,K-1,I)
            ENDDO
          ENDDO
        ENDIF
C
45      CONTINUE
C
      ENDIF
C
      NMAXI=0
      DO J=1,NPOIN
        ELAY%R(J) = ES(J,1)
        IF(NLAYER%I(J).GT.1) THEN
          ESTRAT%R(J) = ES(J,2)           
        ENDIF
C       A REVOIR :
C       ON REMPLIT DE ZERO LES AVAIL INUTILES, CA VAUT MIEUX !!!???
        IF(NLAYER%I(J).LT.10) THEN
          DO I = 1, NSICLA
            DO K = NLAYER%I(J)+1,10
              AVAIL(J,K,I) = 0.D0
            ENDDO
          ENDDO
        ENDIF
        IF(NLAYER%I(J).GT.NMAXI) NMAXI = NLAYER%I(J)
      ENDDO
C 
C     CALCUL DU VOLUME TOTAL DE SEDIMENTS DE CHAQUE CLASSE
C
      DO I=1,NSICLA
        VOLTOT(I)=0.D0
        DO J=1,NPOIN
          DO K=1,NLAYER%I(J)
            VOLTOT(I) = VOLTOT(I) + ES(J,K)*AVAIL(J,K,I)*VOLU2D%R(J)
          ENDDO
        ENDDO
      ENDDO
C
      IF(NCSIZE.GT.1) THEN
        DO I=1,NSICLA
          VOLTOT(I) = P_DSUM(VOLTOT(I))
        ENDDO
      ENDIF
C
      WRITE(LU,*) 'MAXIMUM INITIAL NUMBER OF LAYERS :',NMAXI
      DO I=1,NSICLA
        IF(LNG.EQ.1) THEN
          WRITE(LU,*)'VOLUME TOTAL DE LA CLASSE ',I ,' :',VOLTOT(I)
        ENDIF
        IF(LNG.EQ.1) THEN
          WRITE(LU,*)'TOTAL VOLUME OF CLASS ',I ,' :',VOLTOT(I)
        ENDIF
      ENDDO
C
C-----------------------------------------------------------------------
C      
1800  FORMAT(1X,'IL Y A PLUS DE 10 COUCHES DANS LA STRATIFICATION')
1815  FORMAT(1X,'THERE ARE MORE THAN 10 LAYERS IN THE STRATIFICATION')
C
C-----------------------------------------------------------------------
C
      RETURN
      END
