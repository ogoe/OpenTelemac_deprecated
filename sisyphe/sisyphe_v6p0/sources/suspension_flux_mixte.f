      ! ******************************* !
        SUBROUTINE SUSPENSION_FLUX_MIXTE!
      ! ******************************* !

     &  (TAUP,HN,ACLADM,CS,NPOIN,
     &   CHARR,XMVE,XMVS,GRAV,HMIN,XWC,
     &   ZERO,KARMAN,PARTHENIADES,FLUER_SABLE,FLUER_VASE,ZREF,
     &   AC,CSTAEQ,QSC,ICQ,DEBUG,AVAIL,NSICLA,ES,
     &   TOCE_VASE,NCOUCH_TASS,DT,TOCE_MIXTE,MS_SABLE,MS_VASE)

          ! ================================================= !
          ! COMPUTATION OF THE FLUX OF DEPOSITION AND EROSION !
          ! ================================================= !

!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE INTERFACE_SISYPHE, EX_FLUX_MIXTE=>SUSPENSION_FLUX_MIXTE
      USE BIEF
      USE DECLARATIONS_SISYPHE, ONLY : FDM
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU         


      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE (BIEF_OBJ),  INTENT(IN)    :: TAUP,HN,ACLADM,CS
      INTEGER,          INTENT(IN)    :: NPOIN,DEBUG,NSICLA
      INTEGER,          INTENT(IN)    :: NCOUCH_TASS
      LOGICAL,          INTENT(IN)    :: CHARR
      DOUBLE PRECISION, INTENT(IN)    :: XMVE, XMVS, GRAV, HMIN
      DOUBLE PRECISION, INTENT(IN)    :: XWC
      DOUBLE PRECISION, INTENT(IN)    :: ZERO, KARMAN, PARTHENIADES
      TYPE (BIEF_OBJ),  INTENT(IN)    :: ZREF
      DOUBLE PRECISION, INTENT(INOUT) :: AC,AVAIL(NPOIN,10,NSICLA)
      DOUBLE PRECISION, INTENT(INOUT) :: ES(NPOIN,10)
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: CSTAEQ
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: FLUER_SABLE,FLUER_VASE
      DOUBLE PRECISION,  INTENT(INOUT) :: MS_SABLE(NPOIN,10) 
      DOUBLE PRECISION,  INTENT(INOUT) :: MS_VASE(NPOIN,10)
      DOUBLE PRECISION,  INTENT(INOUT) ::TOCE_MIXTE(NPOIN,10)
C
      DOUBLE PRECISION, INTENT(IN)      :: DT
C
      TYPE(BIEF_OBJ),   INTENT(IN)       ::  QSC
      INTEGER,          INTENT (IN)      :: ICQ
C
      DOUBLE PRECISION, INTENT(IN)     :: TOCE_VASE(10)
C      
      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER I, J,N
      DOUBLE PRECISION FLUERSABLE,FLUERVASE,FLUER_LOC(10)  
C
      DOUBLE PRECISION QE_MOY,TOCE_SABLE,TEMPS,QER_VASE,QER_SABLE     
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!       

      ! ******************************************** !
      ! I - COMPUTATION OF THE CRITICAL SHEAR STRESS
      !    --> TOCE_SABLE
      ! ******************************************** !
C
C CALCUL DU TOCE_SABLE VIA LE PARAMETRE DE SHIELDS
         TOCE_SABLE= AC*(XMVS-XMVE)*GRAV*FDM(1)         
C           
      ! **************************************** !
	!        II-CACUL DE L EROSION      
      ! **************************************** !      
C---------ON NE PASSE DANS LE CALCUL D'EROSION QU'UNE SEULE FOIS (SABLE PAR EXEMPLE
C         CAR LE FLUX CALCULE EST UN FLUX GLOBAL COMMEUN AU 2 SEDIMENT)
C---------JE CALCULE LES FLUX D'EROSION THEORIQUES POUR CHAQUE(SEDIMENT DIPONIBLE DE FACON INFINIE DANS CHAQUE COUCHE)
C
C---------CALCUL DE LA CONTRAINTE CRITIQUE DE CHAQUE COUCHE FONCTION DU % DE VASE
      DO J=1,NCOUCH_TASS
        DO I=1,NPOIN
          IF(AVAIL(I,J,2).LE.0.3D0)THEN
            TOCE_MIXTE(I,J)=TOCE_SABLE
          ELSEIF(AVAIL(I,J,2).GE.0.5D0)THEN
             TOCE_MIXTE(I,J)=TOCE_VASE(J)
          ELSE                  
                 TOCE_MIXTE(I,J)= TOCE_SABLE +
     &   (AVAIL(I,J,2)-0.3D0)*(TOCE_VASE(J)-TOCE_SABLE)/(0.5D0-0.3D0) 
          ENDIF
        ENDDO
      ENDDO
C MODIFS CV INTRODUIRE TOCE EN ARGUMENT 
C         AC(I) = TOCE_MIXTE(I,J)/((XMVS-XMVE)*GRAV*ACLADM%R(I))
C
        IF(ICQ.EQ.1) THEN 
          IF (DEBUG > 0) WRITE(LU,*) 'SUSPENSION_FREDSOE'
C  
           CALL SUSPENSION_FREDSOE(ACLADM,TAUP,NPOIN,
     &         GRAV,XMVE,XMVS,ZERO,AC,CSTAEQ)
C
          IF (DEBUG > 0) WRITE(LU,*) 'END SUSPENSION_FREDSOE'
C
        ELSEIF(ICQ.EQ.2) THEN
C
          IF (DEBUG > 0) WRITE(LU,*) 'SUSPENSION_BIJKER'
C
               CALL SUSPENSION_BIJKER(TAUP,HN,NPOIN,CHARR,QSC,ZREF,
     &                                ZERO,HMIN,CSTAEQ,XMVE)
C
          IF (DEBUG > 0) WRITE(LU,*) 'END SUSPENSION_BIJKER'
C     
        ENDIF
C 
C      DO J=NCOUCH_TASS,1,-1         
C        DO I=1,NPOIN
C           CSTAEQ_COUCHE(I,J)=CSTAEQ%R(I)
C        ENDDO
C      ENDDO
C
      DO I=1,NPOIN
C
        DO J=1,NCOUCH_TASS
C
C-----------CALCUL FLUER_SABLE_VASE FONCTION DU POURCENTAGE DE VASE
C 
          IF(AVAIL(I,J,2).LE.0.3D0)THEN
C-------------FRACTION DE VASE < 30%, LES FUX SONT SEMBLABLES A DU SABLE PUR
            IF(TAUP%R(I).GT.TOCE_MIXTE(I,J))THEN
                 FLUER_LOC(J)=CSTAEQ%R(I)*XWC
            ELSE
               FLUER_LOC(J)=0.D0
            ENDIF
C-------------FRACTION DE VASE > 50%, LES FLUX SONT SEMBLABLES A DE LA VASE PURE
          ELSEIF(AVAIL(I,J,2).GE.0.5D0)THEN
            IF(TAUP%R(I).GT.TOCE_MIXTE(I,J))THEN
               FLUER_LOC(J)=PARTHENIADES*
     &              ((TAUP%R(I)/TOCE_MIXTE(I,J))-1.D0)
            ELSE
               FLUER_LOC(J)=0.D0
            ENDIF              
C-------------FRACTION DE VASE >30% ET <50%, INTERPOLATION DES FLUX
C             ET DE LA LA CONTRAINTE CRITIQUE
          ELSE
            IF(TAUP%R(I).GT.TOCE_MIXTE(I,J))THEN
               FLUERSABLE=CSTAEQ%R(I)*XWC
               FLUERVASE=PARTHENIADES*
     *             ((TAUP%R(I)/TOCE_MIXTE(I,J))-1.D0)
            ELSE
               FLUERSABLE=0.D0
               FLUERVASE=0.D0
            ENDIF
               FLUER_LOC(J)=(AVAIL(I,J,2)-0.3D0)/
     &           (0.5D0-0.3D0)*(FLUERVASE-FLUERSABLE)+FLUERSABLE
          ENDIF          
        ENDDO
    
C
C CALCUL DE LA PROFONDEUR D EROSION ZER_MOY  
C et des masses erodees           
          QER_VASE = 0.D0
          QER_SABLE = 0.D0
C 
          TEMPS= DT
C          
          DO J= 1, NCOUCH_TASS 
           IF(ES(I,J).GE.1.D-6) THEN
C 
C CALCUL DE LA MASSE POTENTIELLEMENT ERODABLE DANS LA COUCHE J (KG/M2)
C       
             QE_MOY= FLUER_LOC(J) *XMVS * TEMPS
C          
             IF(QE_MOY.LT.(MS_SABLE(I,J)
     *            +MS_VASE(I,J))) THEN
C
                  QER_VASE = QER_VASE 
     *                  + QE_MOY*MS_VASE(I,J)/
     *                      (MS_VASE(I,J)+MS_SABLE(I,J))
                  QER_SABLE = QER_SABLE  
     *                    + QE_MOY*MS_SABLE(I,J)
     *                      /(MS_VASE(I,J)+MS_SABLE(I,J))
CV                   
                 GO TO 10
C             
              ELSE
C
                  QER_VASE = QER_VASE + MS_VASE(I,J)
                  QER_SABLE = QER_SABLE + MS_SABLE(I,J)
                 TEMPS= TEMPS -
     *             (MS_SABLE(I,J)+MS_VASE(I,J))
     *                      /FLUER_LOC(J)/XMVS
              ENDIF
          ENDIF
C
         ENDDO
          WRITE(LU,*) 'ATTENTION TOUTES LES COUCHES SONT VIDES'
C          STOP           

  10    CONTINUE   
C
      ! ************************************************ !
      ! II-CACUL DU FLUX D EROSION POUR LE SABLE/LA VASE !   
      ! ************************************************ ! 
C
C Q_VASE REPRESENTE LA MASSE SURFACIQUE DE VASE A ERODER POUR ATTEINDRE ZER_MOY
C Q_SABLE REPRESENTE LA MASSE SURFACIQUE DE VASE A ERODER POUR ATTEINDRE ZER_MOY
C                         
        FLUER_VASE%R(I)  = QER_VASE /(DT*XMVS) 
        FLUER_SABLE%R(I) = QER_SABLE/(DT*XMVS) 
C
      ENDDO                
C
C-----------------------------------------------------------------------
C      
      RETURN      
      END
