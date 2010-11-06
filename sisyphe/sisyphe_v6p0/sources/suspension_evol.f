C recalculer les ES à la fin et verifier le critère 
C ELAY=ZF-ZR
      ! ************************** !
        SUBROUTINE SUSPENSION_EVOL
      ! ************************** !

     *  (ZFCL_S,FLUDP,FLUER,DT, NPOIN,CSF,XMVS, QFLUX,MS,
     *   SEDCO,CONC_VASE,NCOUCH_TASS)
C 
C ---------------------------------------------------------------------C
!                                                                       !
! CALLED BY SUSPENSION_COMPUTATION                                      !
!      CALCUL DE L'EVOLUTION POUR LA VASE EN FONCTION DES FLUDP ET FLUER! 
!      ET MISE A JOUR DE LA MASSE DES COUCHE ET EPAISSEURS DE CHAQUE COUCHE
!             ET EPAISSEUR TOTALE                                       !    
!                                                                       !
!                                                                       !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU

      ! 2/ GLOBAL VARIABLES
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: ZFCL_S,FLUDP,FLUER,QFLUX
      DOUBLE PRECISION, INTENT(IN)    :: DT, XMVS, CSF
      INTEGER, INTENT(IN) :: NPOIN,NCOUCH_TASS
      LOGICAL, INTENT(IN) :: SEDCO 
      DOUBLE PRECISION, INTENT(IN) :: CONC_VASE(NCOUCH_TASS)
      DOUBLE PRECISION,  INTENT(INOUT) :: MS(NPOIN,NCOUCH_TASS)
 
C
      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER :: I,J 
C      
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION CONC,MER
C
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
       ZERO = 1.D-08
C
C CALCUL DU flux de sédiments à chadue pas de temps
C
           CALL OS('X=Y-Z   ', X=QFLUX, Y=FLUDP, Z=FLUER)
           CALL OS('X=CX    ', X=QFLUX, C=DT)
           IF(NCOUCH_TASS.EQ.1)   CALL OS('X=CY    ', 
     *         X=ZFCL_S, Y= QFLUX, C=1.D0/CSF)

           IF(NCOUCH_TASS.GT.1) THEN
 
             DO I = 1, NPOIN       
C
C DEPOT DANS LA PREMIERE COUCHE 
C
             IF (QFLUX%R(I).GT.ZERO) THEN
                ZFCL_S%R(I) = QFLUX%R(I) / CSF  
                MS(I,1) = MS (I,1) +QFLUX%R(I)*XMVS
!
              ELSEIF(QFLUX%R(I).LT.ZERO) THEN 
C
C EROSION DES COUCHES SUCCESSIVES
C               
C
                ZFCL_S%R(I) = 0.D0
                MER = - QFLUX%R(I) *XMVS
C            
                DO J = 1, NCOUCH_TASS
C           
                 IF(.NOT.SEDCO) CONC= XMVS * CSF
                 IF(SEDCO) CONC=XMVS*CONC_VASE(J)
C
                 IF (MER.LE.MS(I,J)) THEN            
                   MS(I,J)= MS(I,J) - MER
                   ZFCL_S%R(I)= ZFCL_S%R(I) - MER/CONC
                   GO TO 40
C
                ELSE
C
C EROSION DE LA TOTALITE DE LA SOUS-COUCHE
C
                   MER= MER - MS(I,J)
                   ZFCL_S%R(I)= ZFCL_S%R(I) - 
     *                MS(I,J)/CONC             
                   MS(I,J)=0.D0
C                                       
               ENDIF
C FIN DE LA BOUCLE SUR LES COUCHES
             ENDDO
C fin erosio
          ENDIF
C
  40      CONTINUE  
C
C FIN DE LA BOUCLE SUR LES POINTS
C
        ENDDO
      ENDIF
!======================================================================!
!======================================================================!
!
      RETURN      
      END
