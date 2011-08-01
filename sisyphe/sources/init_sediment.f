C                       ************************
                        SUBROUTINE INIT_SEDIMENT
C                       ************************
C
     *(NSICLA,ELAY,ZF,ZR,NPOIN,AVAIL,FRACSED_GF,AVA0,  
     * LGRAFED,CALWC,XMVS,XMVE,GRAV,VCE,XWC,FDM, 
     * CALAC,AC, SEDCO, ES, NCOUCH_TASS,CONC_VASE,
     * MS_SABLE, MS_VASE, ACLADM, UNLADM)
C
C***********************************************************************
C SISYPHE VERSION 6.0   30/12/2008    C. VILLARET (LNHE) 01 30 87 83 28
C
C 16/09/2009 JMH : AVAIL(NPOIN,10,NSICLA)
C
C***********************************************************************
C                                                                      
C                                                                      
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____._______________________________________________
C |      NOM       |MODE|                   ROLE                        
C |________________|____|_______________________________________________
C |   NSICLA       | -->| NUMBER OF SEDIMENT CLASSES
C |   ELAY         | -->|       
C |   HN           | -->| WATER DEPTH    
C |   ZF           | -->| BOTTOM                             
C |   ZR           | -->| NON ERODABLE BED                                              
C |   NPOIN        | -->| NUMBER OF POINTS
C |   AVAI         | -->|  
C |   AVAIL        | -->|
C |   FRACSED_GF   | -->|   
C |   AVA0         | -->| 
C |   LGRAFED      | -->|
C |   CALWC        | -->|
C |   XMVS         | -->|
C |   XMVE         | -->| 
C |   GRAV         | -->| GRAVITY ACCELERATION
C |   VCE          | -->|  
C |   Z            | -->|
C |   MESH         | -->|
C |   CHOIX        | -->|
C |   XWC          | -->| 
C |   FDM          | -->|
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C----------------------------------------------------------------------
C PROGRAMME APPELANT : SISYPHE
C PROGRAMMES APPELES : NOEROD, INIT_AVAI
C----------------------------------------------------------------------
C
      USE BIEF
      USE INTERFACE_SISYPHE, EX_INIT_SEDIMENT => INIT_SEDIMENT
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,           INTENT(IN)     :: NSICLA,NPOIN,NCOUCH_TASS
      TYPE(BIEF_OBJ),    INTENT(INOUT)  :: ELAY,ZF,ZR 
CV
      TYPE(BIEF_OBJ), INTENT(INOUT)     :: MS_SABLE, MS_VASE
      TYPE(BIEF_OBJ),    INTENT(INOUT)  :: ACLADM, UNLADM
CV         
      LOGICAL,           INTENT(IN)     :: LGRAFED,CALWC
CV      
      LOGICAL,           INTENT(IN)     :: CALAC
      DOUBLE PRECISION,  INTENT(IN)     :: XMVS,XMVE,GRAV,VCE
      DOUBLE PRECISION,  INTENT(INOUT)  :: AVA0(NSICLA)
      DOUBLE PRECISION,  INTENT(INOUT)  :: AVAIL(NPOIN,10,NSICLA)
      DOUBLE PRECISION,  INTENT(INOUT)  :: FRACSED_GF(NSICLA)    
      DOUBLE PRECISION,  INTENT(INOUT)  :: FDM(NSICLA),XWC(NSICLA) 
CV
      DOUBLE PRECISION,  INTENT(INOUT)  :: AC(NSICLA)
C modif CV
      LOGICAL,           INTENT(IN)     :: SEDCO(NSICLA)     
C
C sedco (1) ou sedco(2) = Yes --> modele de tassement
C
C      
      DOUBLE PRECISION, INTENT(IN)    :: CONC_VASE(10)
      DOUBLE PRECISION, INTENT(INOUT) :: ES(NPOIN,10)  
C end modif  CV
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER            :: I,J
      DOUBLE PRECISION   :: DENS,DSTAR
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
!  ------ COMPOSITION DU LIT 
! UNe seule classse
!
        CALL OS('X=Y-Z   ',X=ELAY,Y=ZF,Z=ZR)
!
      IF(NSICLA.EQ.1) THEN
         DO I=1,NPOIN
          AVAIL(I,1,1) = 1.D0
          ACLADM%R(I) = FDM(1)
        ENDDO             
C vase pure 
        IF(SEDCO(1)) CALL INIT_MIXTE(XMVS,NPOIN,AVAIL,NSICLA,ES,
     *                    ELAY%R, NCOUCH_TASS,CONC_VASE,MS_SABLE%R,
     *                     MS_VASE%R,ZF%R,ZR%R,AVA0)
C
      ELSE 
C  
C  multi-classes  non cohesifs
C
        IF(.NOT.SEDCO(2)) THEN
C
          CALL INIT_AVAI
C         CALL MEAN_GRAIN_SIZE
C this part can be integrated into init_avai
          DO J=1,NPOIN   
            ACLADM%R(J) = 0.D0
            UNLADM%R(J) = 0.D0
            DO I=1,NSICLA
              IF(AVAIL(J,1,I).GT.0.D0) THEN
                ACLADM%R(J) = ACLADM%R(J) + FDM(I)*AVAIL(J,1,I)
                UNLADM%R(J) = UNLADM%R(J) + FDM(I)*AVAIL(J,2,I)
              ENDIF
            ENDDO
            ACLADM%R(J)=MAX(ACLADM%R(J),0.D0)
            UNLADM%R(J)=MAX(UNLADM%R(J),0.D0)
          ENDDO
        ELSE
C        mixte (non cohesif /cohesif)
          CALL INIT_MIXTE(XMVS,NPOIN,AVAIL,NSICLA,ES,ELAY%R,
     *                     NCOUCH_TASS,CONC_VASE,MS_SABLE%R,
     *                     MS_VASE%R,ZF%R,ZR%R,AVA0)
          DO I=1,NPOIN
            ACLADM%R(I) = FDM(1)
          ENDDO 
        ENDIF
CV fin multi classes
      ENDIF
C
      IF(LGRAFED) THEN
        DO I=1, NSICLA
          FRACSED_GF(I)=AVA0(I)
        ENDDO
      ENDIF
C         
C
C ------ Vitesse de chute          
C
      IF(.NOT.CALWC) THEN
        DENS = (XMVS - XMVE) / XMVE 
        DO I = 1, NSICLA
          CALL VITCHU_SISYPHE(XWC(I),DENS,FDM(I),GRAV,VCE)
        ENDDO         
      ENDIF
C
C------ Parametre de Shields 
C
      IF(.NOT.CALAC) THEN      
        DENS  = (XMVS - XMVE )/ XMVE      
        DO I = 1, NSICLA
          DSTAR = FDM(I)*(GRAV*DENS/VCE**2)**(1.D0/3.D0) 
          IF (DSTAR <= 4.D0) THEN
            AC(I) = 0.24*DSTAR**(-1.0D0)
          ELSEIF (DSTAR <= 10.D0) THEN
            AC(I) = 0.14D0*DSTAR**(-0.64D0)
          ELSEIF (DSTAR <= 20.D0) THEN
            AC(I) = 0.04D0*DSTAR**(-0.1D0)
          ELSEIF (DSTAR <= 150.D0) THEN
            AC(I) = 0.013D0*DSTAR**(0.29D0)
          ELSE
            AC(I) = 0.055D0
          ENDIF
        ENDDO
      ENDIF
C

C
C-----------------------------------------------------------------------
C
      RETURN
      END
      
