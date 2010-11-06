      ! ************************* !
        SUBROUTINE BEDLOAD_DIFFIN ! 
      ! ************************* !

     &  (U, V, NBOR, XNEBOR, YNEBOR, KP1BOR, MASKEL, NELBOR, NPTFR,
     &   KENT, KSORT, KLOG, KINC, KDIR, KDDL, KNEU, MSK, CLT, LITBOR,
     &   MASKTR, LIMTRA)


C**********************************************************************
C SISYPHE VERSION 6.0  17/08/2004  FRANCOIS MENARD (STAGIAIRE LNHE)    
C**********************************************************************


               ! ===================================== !
               ! Initialisation des conditions limites !
               !                                       !
               ! ===================================== !


C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT
C**********************************************************************C
C                                                                      C
C                 SSSS I   SSSS Y   Y PPPP  H   H EEEEE                C
C                S     I  S      Y Y  P   P H   H E                    C
C                 SSS  I   SSS    Y   PPPP  HHHHH EEEE                 C
C                    S I      S   Y   P     H   H E                    C
C                SSSS  I  SSSS    Y   P     H   H EEEEE                C
C                                                                      C
C----------------------------------------------------------------------C
C                             ARGUMENTS                                C
C .____________.____.__________________________________________________C
C |    NOM     |MODE|                    ROLE                          C
C |____________|____|__________________________________________________C
C |____________|____|__________________________________________________C
C                                                                      C
C                =>  Can't be change                                   C
C                <=> Can be change                                     C
C                <=  Must be set                                       C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY BEDLOAD_MAIN                                               !
!                                                                      !
! CALL      ------                                                     !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE INTERFACE_SISYPHE, EX_BEDLOAD_DIFFIN => BEDLOAD_DIFFIN
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU


      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE(BIEF_OBJ), INTENT(IN)    :: U, V, NBOR, XNEBOR, YNEBOR
      TYPE(BIEF_OBJ), INTENT(IN)    :: KP1BOR, MASKEL, NELBOR
      INTEGER,        INTENT(IN)    :: NPTFR, KENT, KSORT, KLOG
      INTEGER,        INTENT(IN)    :: KINC, KDIR, KDDL, KNEU
      LOGICAL,        INTENT(IN)    :: MSK
      TYPE(BIEF_OBJ), INTENT(INOUT) :: CLT
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: LITBOR, MASKTR, LIMTRA


      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER            :: K, K1, K2
      DOUBLE PRECISION   :: USCALN,C
      INTEGER, PARAMETER :: DIR = 1
      INTEGER, PARAMETER :: DDL = 2
      INTEGER, PARAMETER :: NEU = 3
      INTEGER, PARAMETER :: OND = 4


!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!

      ! ****************************************************** !
      ! I - TYPES DE CONDITIONS AUX LIMITES DU TRACEUR         ! (_IMP_)
      !     EVENTUELLEMENT MODIFIE EN FONCTION DU SIGNE DE U.N ! (_IMP_)
      !     SUR LES FRONTIERES LIQUIDES (N : NORMALE SORTANTE) ! (_IMP_)
      ! ****************************************************** !
      DO K = 1, NPTFR

         CLT%I(K) = LITBOR%I(K)

         ! I.1 - Reperage des frontieres liquides (_IMP_)
         ! --------------------------------------
         IF (CLT%I(K) == KENT) THEN
            USCALN = U%R(NBOR%I(K))*XNEBOR%R(K)
     &             + V%R(NBOR%I(K))*YNEBOR%R(K)

            ! Vitesse sortante, traceur libre
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            IF (USCALN >= 0.D0) CLT%I(K) = KSORT

         ELSEIF(CLT%I(K) == KSORT) THEN

            USCALN = U%R(NBOR%I(K))*XNEBOR%R(K)
     &             + V%R(NBOR%I(K))*YNEBOR%R(K)

            ! Vitesse entrante, traceur libre
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            IF (USCALN <= 0.D0) CLT%I(K) = KENT 
         ENDIF

      ENDDO


      ! **************************************************************** !
      ! II - CONSTRUCTION DU TABLEAU MASKTR EN FONCTION DE CLT           ! (_IMP_)
      !      MASKTR VAUT 1 POUR UN SEGMENT DE TYPE NEUMANN ET ZERO SINON ! (_IMP_)
      !      UN SEGMENT EST DE TYPE NEUMANN SI UN DE SES POINTS AU MOINS ! (_IMP_)
      !      EST DONNE COMME NEUMANN PAR L'UTILISATEUR.                  ! (_IMP_)
      ! **************************************************************** !
      CALL OS('X=0     ', X=MASKTR)

      DO K1 = 1 , NPTFR

         K2 = KP1BOR%I(K1)

         ! II.1 - SEGMENTS DE TYPE NEUMANN
         ! -------------------------------
         IF (CLT%I(K1).EQ.KLOG.OR.CLT%I(K2).EQ.KLOG) THEN
            MASKTR%ADR(NEU)%P%R(K1) = 1.D0

         ! II.2 - SEGMENTS DE TYPE SORTIE (_IMP_)
         ! ------------------------------
         ELSEIF ((CLT%I(K1) == KENT) .AND. (CLT%I(K2) == KSORT)) THEN
            MASKTR%ADR(DDL)%P%R(K1) = 1.D0

         ELSEIF ((CLT%I(K1) == KSORT) .OR. (CLT%I(K2) == KSORT)) THEN
            MASKTR%ADR(DDL)%P%R(K1) = 1.D0

         ! II.3 - SEGMENTS DE TYPE SORTIE (_IMP_)
         ! ------------------------------
         ELSEIF ((CLT%I(K1) == KSORT) .AND. (CLT%I(K2) == KENT)) THEN
            MASKTR%ADR(DDL)%P%R(K1) = 1.D0
         ELSEIF ((CLT%I(K1) == KENT) .OR. (CLT%I(K2) == KENT)) THEN
            MASKTR%ADR(DIR)%P%R(K1) = 1.D0
         ELSEIF ((CLT%I(K1) == KINC) .OR. (CLT%I(K2) == KINC)) THEN
            MASKTR%ADR(OND)%P%R(K1)=1.D0
         ELSE
            IF (LNG == 1) WRITE(LU,101)
            IF (LNG == 2) WRITE(LU,102)
            CALL PLANTE(1)
         ENDIF
      ENDDO


      ! *********************** !
      ! III - MASQUAGE EVENTUEL !
      ! *********************** !
      IF(MSK) THEN                                                      
        DO K1 = 1 , NPTFR 
          C=MASKEL%R(NELBOR%I(K1))                                            
          MASKTR%ADR(DIR)%P%R(K1) = MASKTR%ADR(DIR)%P%R(K1)*C
          MASKTR%ADR(DDL)%P%R(K1) = MASKTR%ADR(DDL)%P%R(K1)*C
          MASKTR%ADR(NEU)%P%R(K1) = MASKTR%ADR(NEU)%P%R(K1)*C
          MASKTR%ADR(OND)%P%R(K1) = MASKTR%ADR(OND)%P%R(K1)*C
        ENDDO
      ENDIF


      ! ************************************************************** !
      ! IV - PASSAGE DES CONDITION PHYSIQUES AUX CONDITIONS TECHNIQUES !
      ! ************************************************************** !
      DO K = 1, NPTFR                                                    

         ! IV.1 - ENTREE DE DOMAINE : TRACEUR IMPOSE (_IMP_)
         ! -----------------------------------------
         IF(CLT%I(K).EQ.KENT) THEN
            LIMTRA%I(K) = KDIR

         ELSEIF(CLT%I(K).EQ.KSORT) THEN
            LIMTRA%I(K) = KDDL

         ! IV.2 - PAROI : CONDITIONS DE NEUMANN (_IMP_)
         ! ------------------------------------
         ELSEIF(CLT%I(K).EQ.KLOG ) THEN
            LIMTRA%I(K) = KNEU

         ! IV.3 - ERREUR, VALEUR DE LITBOR INCONNUE (_IMP_)
         ! ----------------------------------------
         ELSE
            IF (LNG == 1) WRITE(LU,11) K, LITBOR%I(K)
            IF (LNG == 2) WRITE(LU,12) K, LITBOR%I(K)
            CALL PLANTE(1)
            STOP      
         ENDIF 
                                                                   
      ENDDO

      !----------------------------------------------------------------! 
101   FORMAT(' DIFFIN_SISYPHE : CAS NON PREVU')
11    FORMAT(' DIFFIN_SISYPHE : POINT ',1I8,' LITBOR= ',1I8,' ?')
      !----------------------------------------------------------------! 
102   FORMAT(' DIFFIN_SISYPHE: UNEXPECTED CASE')
12    FORMAT(' DIFFIN_SISYPHE : POINT ',1I8,' LITBOR= ',1I8,' ?')
      !----------------------------------------------------------------! 


!======================================================================!
!======================================================================!

      RETURN
      END
