      ! *************************** !
        SUBROUTINE BEDLOAD_SOLVS_VF ! 
      ! *************************** !

     &(MESH,QSX,QSY,LIEBOR,UNSV2D,EBOR,BREACH,NSEG,NPTFR,
     & NPOIN,KENT,KSORT,DT,T10,ZFCL,FLUX)

C
C**********************************************************************
C SISYPHE VERSION 5.8  30/10/2007  J-M HERVOUET
C SISYPHE VERSION 5.5  07/05/2002  M. GONZALES DE LINARES            
C SISYPHE VERSION 5.5  14/09/2004  F. HUVELIN 
C
C JMH 15/09/09 : KENT KSORT ADDED (WERE HARDCODED BEFORE !!!)
C                       
C**********************************************************************
C
C FONCTION : RESOLUTION OF EXNER EQUATION WITH THE FINITE VOLUME METHOD
C
C 30/10/2007  J-M HERVOUET (UNSV2D +DIRICL SUPPRIME)
C
C
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
C .________________.____.______________________________________________C
C |      NOM       |MODE|                   ROLE                       C
C |________________|____|______________________________________________C
C |________________|____|______________________________________________C
C                    <=  Can't be change by the user                   C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY BEDLOAD_EVOL                                               !                                                                      !                                             !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE INTERFACE_SISYPHE, EX_BEDLOAD_SOLVS_VF => BEDLOAD_SOLVS_VF
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU


      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE(BIEF_MESH),  INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ),   INTENT(IN)    :: QSX, QSY
      TYPE(BIEF_OBJ),   INTENT(IN)    :: LIEBOR,UNSV2D, EBOR
      TYPE(BIEF_OBJ),   INTENT(IN)    :: BREACH
      INTEGER,          INTENT(IN)    :: NSEG,NPTFR,NPOIN,KENT,KSORT
      DOUBLE PRECISION, INTENT(IN)    :: DT
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: T10
      TYPE(BIEF_OBJ),   INTENT(INOUT)   :: ZFCL, FLUX


      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER          :: ISEGIN, K
      INTEGER          :: IEL, IEL1, IEL2
      DOUBLE PRECISION :: QSMOY1, QSMOY2
      DOUBLE PRECISION :: QSP
      DOUBLE PRECISION :: VNOIN1, VNOIN2, RNORM
      DOUBLE PRECISION :: XN, YN

!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!

      ! ***************** !
      ! I - INTIALIZATION ! 
      ! ***************** !
      
      CALL OS('X=0     ', X=FLUX)

      ! ************************************************* !
      ! II - DETERMINER LE FLUX SORTANT DE CHAQUE CELLULE ! 
      ! ************************************************* !
      
      DO ISEGIN = 1, NSEG

         IEL1 = MESH%NUBO%I(2*ISEGIN - 1)
         IEL2 = MESH%NUBO%I(2*ISEGIN    )

         ! II.1 - RNORM : LONGUEUR DU SEGMENT 
         ! ----------------------------------
         VNOIN1 = MESH%VNOIN%R(3*ISEGIN - 2)
         VNOIN2 = MESH%VNOIN%R(3*ISEGIN - 1)
         RNORM  = MESH%VNOIN%R(3*ISEGIN    )

         ! II.2 - QS AU SEGMENT DECOMPOSE SUIVANT X ET Y 
         ! ---------------------------------------------
         QSMOY1 = 0.5D0*(QSX%R(IEL1) + QSX%R(IEL2))
         QSMOY2 = 0.5D0*(QSY%R(IEL1) + QSY%R(IEL2))

         ! II.3 - PROJECTION DE QS AU SEGMENT SUR LA NORMALE DU SEGMENT 
         ! ------------------------------------------------------------
         QSP = VNOIN1*QSMOY1 + VNOIN2*QSMOY2

         ! II.4 - UPWIND SCHEME ON NODES WITH A "PROBLEM" 
         ! ----------------------------------------------
         IF(BREACH%I(IEL1).EQ.1.AND.QSP.GT.0.D0) THEN
           QSMOY1 = QSX%R(IEL1)
           QSMOY2 = QSY%R(IEL1)
         ENDIF
         IF(BREACH%I(IEL2).EQ.1.AND.QSP.LT.0.D0) THEN
           QSMOY1 = QSX%R(IEL2)
           QSMOY2 = QSY%R(IEL2)
         ENDIF
   
         QSP = VNOIN1*QSMOY1 + VNOIN2*QSMOY2

         ! II.5 - INTEGRATION PAR LA LONGUEUR DU SEGMENT 
         ! ---------------------------------------------
         FLUX%R(IEL1) = FLUX%R(IEL1) + RNORM*QSP
         FLUX%R(IEL2) = FLUX%R(IEL2) - RNORM*QSP

      ENDDO

      ! ******************************* !
      ! III - TRAITEMENT DES FRONTIERES ! (_IMP_)
      ! ******************************* !
      DO K = 1 , NPTFR

         IEL = MESH%NBOR%I(K)

         ! III.1 - EVOLUTION LIBRE : ON LAISSE SORTIR LES SEDIMENTS 
         ! --------------------------------------------------------
         IF (LIEBOR%I(K).EQ.KSORT) THEN

            ! XNEBOR(*+NPTFR) ET YNEBOR(*+NPTFR)
            ! CONTIENNENT LE VECTEUR NORMAL A UN POINT FRONTIERE
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            XN = MESH%XNEBOR%R(K+NPTFR)
            YN = MESH%YNEBOR%R(K+NPTFR)

            ! AJOUT DE LA CONTRIBUTION DU FLUX SUR LE SEGMENT FRONTIERE
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            FLUX%R(IEL) = FLUX%R(IEL) + QSX%R(IEL)*XN + QSY%R(IEL)*YN

         ENDIF

         ! III.2 - POUR LA PAROI SOLIDE IL N'Y A RIEN A PROGRAMMER 
         !         CAR LE FLUX DE SEDIMENTS Y EST NUL              
         ! -------------------------------------------------------

      ENDDO

      IF(NCSIZE.GT.1) CALL PARCOM(FLUX, 2, MESH)

      ! ************************** !
      ! IV - RESOLUTION DU SYSTEME ! 
      ! ************************** !

      ! Signe moins car flux sortant positif
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL OS('X=CYZ   ', X=ZFCL, Y=FLUX, Z=UNSV2D, C=-DT)
!
      DO K=1,NPTFR
        IF(LIEBOR%I(K).EQ.KENT) THEN
          ZFCL%R(MESH%NBOR%I(K)) = EBOR%R(K)
        ENDIF
      ENDDO
!
!======================================================================!
!======================================================================!
!
      RETURN
      END
