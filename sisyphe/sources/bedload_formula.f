
      ! ************************** !
        SUBROUTINE BEDLOAD_FORMULA 
      ! ************************** !

     &(U2D,V2D,UCMOY,HN,CF,MU,TOB,TOBW,UW,TW,THETAW,FW, 
     & ACLADM, UNLADM,KSP,KSR,AVA,NPOIN,ICF,HIDFAC,XMVS,XMVE,
     & DM,GRAV,VCE,XKV,HMIN,XWC,D90,KARMAN,ZERO,
     & PI,SUSP, AC, HIDING, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10,
     & T11,TETAP, QSC, QSS,IELMT,SECCURRENT,SLOPEFF, 
     & COEFPN,BIJK,HOULE)


C**********************************************************************C
C SISYPHE VERSION 5.6  12/01/2005  F. HUVELIN                          C
C SISYPHE VERSION 5.4  --/10/2003  C. VILLARET                         C
C SISYPHE VERSION 5.2  --/01/2002  BUI MINH DUC                        C
C**********************************************************************C


               ! ===================================== !
               ! Computation of the bed-load transport !
               ! ===================================== !


C COPYRIGHT EDF-BAW-IFH
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
C                    <=  Can't be changed by the user                  C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY BEDLOAD_SOLIDISCHARGE                                      !
!                                                                      !
! CALL      BEDLOAD_MEYER                                              !
!           BEDLOAD_EINST                                              !
!           BEDLOAD_ENGEL                                              !
!           BEDLOAD_ENGEL_OLD                                          !
!           BEDLOAD_BIJKER                                             !
!           BEDLOAD_SOULSBY                                            !
!           BEDLOAD_HUNZ_MEYER                                         !
!           BEDLOAD_VANRIJN                                            !
!           BEDLOAD_BAILARD                                            !
!           BEDLOAD_DIBWAT   
!           BEDLOAD_CHENG                                         !
!           QSFORM                                                     !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE INTERFACE_SISYPHE,EX_BEDLOAD_FORMULA => BEDLOAD_FORMULA
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU


      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE(BIEF_OBJ),   INTENT(IN)    :: U2D, V2D, UCMOY,HN, CF, TOB
      TYPE(BIEF_OBJ),   INTENT(IN)    :: MU,TOBW, UW, TW, THETAW, FW
      TYPE(BIEF_OBJ),   INTENT(IN)    :: ACLADM,UNLADM,KSR,KSP
      INTEGER,          INTENT(IN)    :: NPOIN, ICF, HIDFAC,IELMT
      DOUBLE PRECISION, INTENT(IN)    :: XMVS, XMVE, DM, GRAV, VCE
      DOUBLE PRECISION, INTENT(IN)    :: XKV, HMIN, XWC, D90
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN, ZERO, PI
      LOGICAL,          INTENT(IN)    :: SUSP,SECCURRENT,HOULE
      DOUBLE PRECISION, INTENT(INOUT) :: AC
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: HIDING
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: T1, T2, T3, T4, T5, T6, T7
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: T8, T9, T10,T11
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: TETAP ! work array T12
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: QSC, QSS
      TYPE(BIEF_OBJ),   INTENT(INOUT) ::  COEFPN
      INTEGER,          INTENT(IN)    :: SLOPEFF 
C
      DOUBLE PRECISION, INTENT (IN) :: BIJK,AVA(NPOIN)
      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER                     :: I
      DOUBLE PRECISION            :: DENS,DSTAR,ALPHA
      DOUBLE PRECISION, PARAMETER :: ZERO_LOCAL = 1.D-6
      DOUBLE PRECISION            :: C1


!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!


      ! *************************** !
      ! I - ADIMENSIONAL PARAMETERS !
      ! *************************** !
      
      DENS  = (XMVS - XMVE )/ XMVE      
      DSTAR = DM*(GRAV*DENS/VCE**2)**(1.D0/3.D0) 

      ! ************************ !
      ! II -  FROTTEMENT DE PEAU !
      ! ************************ !
! 
      C1 = 1.D0/(DENS*XMVE*GRAV*DM)      
      CALL OS('X=CYZ   ', X=TETAP, Y=TOB,Z=MU,  C=C1) 
      CALL OS('X=+(Y,C)', X= TETAP,Y=TETAP, C=ZERO_LOCAL)
!      
      ! *********************************************** !
      ! III - NOMBRE DE SHIELDS CRITIQUE D'ENTRAINEMENT !
      !       FORMULE DE VAN RIJN                       !
      ! *********************************************** !
! deplace dans init_sediment.f
!      
!      IF(ICF.EQ.7.OR.ICF.EQ.6.OR.AC.LE.0.D0) THEN
!         IF (DSTAR <= 4.D0) THEN
!            AC = 0.24*DSTAR**(-1.0D0)
!         ELSEIF (DSTAR <= 10.D0) THEN
!            AC = 0.14D0*DSTAR**(-0.64D0)
!         ELSEIF (DSTAR <= 20.D0) THEN
!            AC = 0.04D0*DSTAR**(-0.1D0)
!         ELSEIF (DSTAR <= 150.D0) THEN
!            AC = 0.013D0*DSTAR**(0.29D0)
!         ELSE
!            AC = 0.055D0
!         ENDIF
!      ENDIF

      IF(SECCURRENT) CALL BEDLOAD_SECCURRENT(IELMT)

      ! ****************************************** !
      ! IV - CALCUL DES 2 COMPOSANTES DU TRANSPORT !
      !      QSS : SUSPENSION                      ! 
      !      QSC : CHARRIAGE                       ! 
      ! ****************************************** !

      ! ===================================== !
      ! IV(1) - FORMULE DE MEYER-PETER-MULLER ! 
      !         CHARRIAGE SEUL                ! 
      ! ===================================== !
      
      IF(ICF == 1) THEN   

          CALL BEDLOAD_MEYER(TETAP,HIDING,HIDFAC,DENS,GRAV,DM,AC,
     &                       T1,QSC,SLOPEFF,COEFPN)
          DO I=1,NPOIN
            QSC%R(I)=XKV*QSC%R(I)*AVA(I)
          ENDDO
          ALPHA = -3.D0

      ! =========================== !
      ! IV(2) - FORMULE DE EINSTEIN ! 
      !         CHARRIAGE SEUL      !
      ! =========================== !
      
      ELSEIF(ICF == 2) THEN

         CALL BEDLOAD_EINST(TETAP,NPOIN,DENS,GRAV,DM,DSTAR,QSC)
         DO I=1,NPOIN
           QSC%R(I)=XKV*QSC%R(I)*AVA(I)*HIDING%R(I)
         ENDDO
         ALPHA = -6.D0

      ! =================================== !
      ! IV(30) - FORMULE DE ENGELUND-HANSEN !
      !          TRANSPORT TOTAL            ! 
      ! =================================== !
      
      ELSEIF(ICF == 30) THEN
C v6p0 MU remplace CF
C attention differences
C         CALL BEDLOAD_ENGEL(TETAP,DENS,GRAV,DM,QSC)
C retour version antérieure
         CALL BEDLOAD_ENGEL(TOB,CF,DENS,GRAV,DM,XMVE,T1,QSC)
C        Repartition arbitraire
         DO I=1,NPOIN
           QSC%R(I)=XKV*QSC%R(I)*AVA(I)*HIDING%R(I)
         ENDDO
         ALPHA = -5.D0

      ! ======================================== !
      ! IV(3) - FORMULE DE ENGELUND-HANSEN       !                 
      !         MODIFICATION DE CHOLLET ET CUNGE ! 
      !         TRANSPORT TOTAL                  !  
      ! ======================================== !
      
      ELSEIF(ICF == 3) THEN
C        KSP remplace CFP
         CALL BEDLOAD_ENGEL_OLD
     &        (TETAP,CF,NPOIN,GRAV,DM,DENS,T1,QSC)
C        Repartition arbitraire
         DO I=1,NPOIN
           QSC%R(I)=XKV*QSC%R(I)*AVA(I)*HIDING%R(I)
         ENDDO
         ALPHA = -5.D0

      ! ============================== !
      ! IV(4) - FORMULE DE BIJKER      ! 
      !         CHARRIAGE + SUSPENSION ! 
      ! ============================== !
      
      ELSEIF (ICF == 4) THEN
      
         CALL BEDLOAD_BIJKER
     &    (TOBW,TOB,MU,KSP,KSR,HN,NPOIN,DM,DENS,XMVE,GRAV,
     &     XWC,KARMAN,ZERO,T4,T7,T8,T9,QSC,QSS,BIJK,HOULE)
         DO I=1,NPOIN
           QSC%R(I)=XKV*QSC%R(I)*AVA(I)*HIDING%R(I)
           QSS%R(I)=XKV*QSS%R(I)*AVA(I)*HIDING%R(I)
         ENDDO
         ALPHA = -1.D0

      ! ============================== !
      ! IV(5) - FORMULE DE SOULSBY     ! 
      !         CHARRIAGE + SUSPENSION ! 
      ! ============================== !
      
      ELSEIF (ICF == 5) THEN

         CALL BEDLOAD_SOULSBY
     &        (UCMOY,HN,UW,NPOIN,DENS,GRAV,DM,DSTAR,HMIN,
     &         D90,QSC,QSS)
         DO I=1,NPOIN
           QSC%R(I)=XKV*QSC%R(I)*AVA(I)*HIDING%R(I)
           QSS%R(I)=XKV*QSS%R(I)*AVA(I)*HIDING%R(I)
         ENDDO
         ALPHA = -4.6D0

      ! ================================================== !
      ! IV(6) - FORMULE DE HUNZIKER / MEYER-PETER & MULLER ! 
      !         CHARRIAGE SEUL                             ! 
      ! ================================================== !
      
      ELSEIF (ICF == 6) THEN

         CALL BEDLOAD_HUNZ_MEYER
     &        (TOB, MU, ACLADM, UNLADM, NPOIN, DENS, XMVE, GRAV,
     &         DM, AC, T1, T2, T3, HIDING, QSC)
         DO I=1,NPOIN
           QSC%R(I)=XKV*QSC%R(I)*AVA(I)
         ENDDO
         ALPHA = -3.D0

      ! =========================== !
      ! IV(7) - FORMULE DE VAN RIJN ! 
      !         CHARRIAGE SEUL      ! 
      ! =========================== !
      
      ELSEIF (ICF == 7) THEN
C
         CALL BEDLOAD_VANRIJN
     &        (TOB,MU,NPOIN,DM,DENS,GRAV,DSTAR,AC,QSC)
         DO I=1,NPOIN
           QSC%R(I)=XKV*QSC%R(I)*AVA(I)*HIDING%R(I)
         ENDDO
         ALPHA = -4.2D0

      ! ============================== !
      ! IV(8) - FORMULE DE BAILARD     ! 
      !         CHARRIAGE + SUSPENSION ! 
      ! ============================== !
      
      ELSEIF (ICF == 8) THEN
C
         CALL BEDLOAD_BAILARD
     &        (U2D,V2D,UCMOY,TOB,TOBW,THETAW,UW,FW,CF,NPOIN,
     &         PI,XMVE,GRAV,DENS,XWC,T1,T2,T3,T4,T5,T6,T7,
     &         T8,T9,T10,T11,QSC,QSS,HOULE)
         DO I=1,NPOIN
           QSC%R(I)=XKV*QSC%R(I)*AVA(I)*HIDING%R(I)
           QSS%R(I)=XKV*QSS%R(I)*AVA(I)*HIDING%R(I)
         ENDDO
         ALPHA = -3.D0

      ! ======================================= !
      ! IV(9) - FORMULE DE DIBAJNIA ET WATANABE ! 
      !         TRANSPORT TOTAL                 ! 
      ! ======================================= !
      
      ELSEIF(ICF == 9) THEN
C
         CALL BEDLOAD_DIBWAT
     &        (U2D,V2D,UCMOY, CF, TOB, TOBW, UW, TW, FW, THETAW,
     &         NPOIN, XMVE, DENS, GRAV, DM, XWC, PI, T1, T2, T3, T4,
     &         T5, T6, T7, T8, T9, T10, T11, QSC,HOULE)
C        Repartition arbitraire
         DO I=1,NPOIN
           QSC%R(I)=XKV*QSC%R(I)*AVA(I)*HIDING%R(I)
         ENDDO
         ALPHA = -3.D0
      ! ======================================= !
      ! IV(20) - FORMULE DE DIBAJNIA ET WATANABE ! 
      !         TRANSPORT TOTAL                 ! 
      ! ======================================= !
      
      ELSEIF(ICF == 20) THEN

         CALL BEDLOAD_CHENG(TETAP,NPOIN,DENS,GRAV,DM,DSTAR,QSC)
         DO I=1,NPOIN
           QSC%R(I)=XKV*QSC%R(I)*AVA(I)*HIDING%R(I)
         ENDDO
         ALPHA = -6.D0

      ! ============================================ !
      ! IV(0) - FORMULE PROGRAMMEE PAR L'UTILISATEUR !
      ! ============================================ !
      
      ELSEIF (ICF == 0) THEN

         ALPHA = -1.D0 ! Initialisation de alpha
         CALL QSFORM                 
         DO I=1,NPOIN
           QSC%R(I)=XKV*QSC%R(I)*AVA(I)*HIDING%R(I)
           QSS%R(I)=XKV*QSS%R(I)*AVA(I)*HIDING%R(I)
         ENDDO

      ! ================= !
      ! IV(else) - ERREUR ! 
      ! ================= !
      
      ELSE
        IF(LNG == 1) WRITE(LU,200) ICF
        IF(LNG == 2) WRITE(LU,201) ICF
200     FORMAT(1X,'TRANSP : FORMULE DE TRANSPORT INCONNUE :',1I6)
201     FORMAT(1X,'TRANSP : TRANSPORT FORMULA UNKNOWN:',1I6)
        CALL PLANTE(1)
        STOP
      ENDIF
!
!-----------------------------------------------------------------------
!
!     WHEN SUSPENSION IS NOT ASKED SPECIFICALLY, SOME BEDLOAD TRANSPORT
!     FORMULAS GIVE A VALUE
!
      IF(.NOT.SUSP) THEN
        IF(ICF.EQ.4.OR.ICF.EQ.5.OR.ICF.EQ.8.OR.ICF.EQ.0) THEN
          DO I = 1,NPOIN
            QSC%R(I) = QSC%R(I) + QSS%R(I)
          ENDDO
        ELSE
!         NOTE JMH: IS THIS REALLY USEFUL ???
          DO I = 1,NPOIN
            QSS%R(I) = 0.D0
          ENDDO
        ENDIF
      ENDIF     
!
!=======================================================================
!=======================================================================
!
      RETURN
      END
