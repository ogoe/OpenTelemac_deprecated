      ! ************************* ! 
        SUBROUTINE BEDLOAD_EFFPNT ! (_IMP_) 
      ! ************************* ! 
 
     & (MASKEL,LIQBOR,S,ZF,U2D,V2D,UCMOY,NPOIN,NPTFR,IELMT,KENT, 
     &  BETA,PI,MSK,MESH,DZFDX,DZFDY,CTETA,STETA, 
     &  COEF,CALFA,SALFA,SLOPEFF,PHISED,DEVIA,BETA2, 
     &  TOB,XMVS,XMVE,DM,GRAV,UNSV2D) 
 
C**********************************************************************C 
C SISYPHE VERSION 5.1  11/09/1995  E. PELTIER                          C 
C SISYPHE VERSION 5.1  11/09/1995  C. LENORMANT                        C 
C SISYPHE VERSION 5.1  11/09/1995  J.-M. HERVOUET                      C 
C**********************************************************************C 
 
             ! ========================================= ! 
             ! Calcul des parametres de l'effet de pente ! 
             ! ========================================= ! 
 
 
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
! CALLED BY BEDLOAD_INIT                                               ! 
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
      USE INTERFACE_SISYPHE,EX_BEDLOAD_EFFPNT => BEDLOAD_EFFPNT 
      USE BIEF 
      IMPLICIT NONE 
      INTEGER LNG,LU 
      COMMON/INFO/LNG,LU 
 
 
      ! 2/ GLOBAL VARIABLES 
      ! ------------------- 
      TYPE(BIEF_OBJ),   INTENT(IN)    :: MASKEL,LIQBOR,S,UNSV2D 
      TYPE(BIEF_OBJ),   INTENT(IN)    :: ZF, U2D,V2D, UCMOY, TOB 
      INTEGER,          INTENT(IN)    :: NPOIN, NPTFR, IELMT, KENT 
      INTEGER,          INTENT(IN)    :: SLOPEFF,DEVIA 
      DOUBLE PRECISION, INTENT(IN)    :: BETA, PI, PHISED, BETA2 
      DOUBLE PRECISION, INTENT(IN)    :: XMVS, XMVE, GRAV, DM 
      LOGICAL,          INTENT(IN)    :: MSK 
      TYPE(BIEF_MESH),  INTENT(INOUT) :: MESH 
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: DZFDX, DZFDY 
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: CTETA,STETA 
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: COEF, CALFA, SALFA 
 
 
      ! 3/ LOCAL VARIABLES 
      ! ------------------ 
      INTEGER          :: I, K 
      DOUBLE PRECISION :: C,ZETA,C1,CALPHA,SALPHA,AA,BB
      DOUBLE PRECISION :: CPSI,SPSI,DZF,TANPHI,CZETA,SZETA,SURBETA2 
      DOUBLE PRECISION :: NORM ,TT1 
! 
!======================================================================! 
!======================================================================! 
!                               PROGRAMME                              ! 
!======================================================================! 
!======================================================================! 
! 
!     DETERMINATION DE COS ET SIN TETA 
!     TETA = ANGLE DE L'ECOULEMENT PAR RAPPORT A AXE X 
!
      DO I=1,NPOIN
        IF(UCMOY%R(I).GE.1.D-12) THEN
          CTETA%R(I)=U2D%R(I)/UCMOY%R(I)
          STETA%R(I)=V2D%R(I)/UCMOY%R(I)
        ELSE
          CTETA%R(I)=1.D0
          STETA%R(I)=0.D0
        ENDIF
      ENDDO
! 
!----------------------------------------------------------------------  
! 
!     CALCUL DE LA PENTE  : D(ZF)/DX ET D(ZF)/DY (VALEURS NODALES) 
! 
      CALL VECTOR(DZFDX, '=', 'GRADF          X',IELMT,1.D0,ZF,S,S, 
     *            S,S,S,MESH,MSK,MASKEL) 
      CALL VECTOR(DZFDY, '=', 'GRADF          Y',IELMT,1.D0,ZF,S,S, 
     *            S,S,S,MESH,MSK,MASKEL) 
C 
      IF(NCSIZE.GT.1) THEN
        CALL PARCOM(DZFDX,2,MESH)
        CALL PARCOM(DZFDY,2,MESH)
      ENDIF 
C 
      CALL OS('X=XY    ',X=DZFDX,Y=UNSV2D)  
      CALL OS('X=XY    ',X=DZFDY,Y=UNSV2D)    
! 
!====================================================================== 
!
!     CALCUL DE L'ANGLE DU TRANSPORT SOLIDE ALFA = TETA + DEVIATION 
!
!     1 : KOCH ET FLOKSTRA
!
      IF(DEVIA==1) THEN 
! 
      C = 2.D0*(XMVS-XMVE)*GRAV*DM/3.D0
      DO I=1,NPOIN
        TT1=C/MAX(TOB%R(I),1.D-10)
        AA=STETA%R(I)-TT1*DZFDY%R(I)
        BB=CTETA%R(I)-TT1*DZFDX%R(I)
        NORM=MAX(SQRT(AA**2+BB**2),1.D-10)
        SALFA%R(I)=AA/NORM
        CALFA%R(I)=BB/NORM 
      ENDDO 
! 
!     2 : TALMON ET AL. JHR 1995 33(4) 
! 
      ELSEIF(DEVIA==2) THEN 
! 
      SURBETA2=1.D0/BETA2
      C = (XMVS-XMVE)*GRAV*DM*SURBETA2**2
      DO I=1,NPOIN
        TT1=SQRT(C/MAX(TOB%R(I),1.D-10))
        AA=STETA%R(I)-TT1*DZFDY%R(I)
        BB=CTETA%R(I)-TT1*DZFDX%R(I)
        NORM=MAX(SQRT(AA**2+BB**2),1.D-10)
        SALFA%R(I)=AA/NORM
        CALFA%R(I)=BB/NORM 
      ENDDO 
!  
      ENDIF 
! 
!====================================================================== 
! 
!     CALCUL DE COEF POUR LA PRISE EN COMPTE DE L'EFFET DE PENTE 
!     SUR L'AMPLITUDE DU TRANSPORT SOLIDE                         
! 
!     METHODE 1 (EMPIRICAL METHOD) 
! 
      IF(SLOPEFF==1) THEN 
! 
        DO I=1,NPOIN
          COEF%R(I)=MAX(0.D0,
     *    1.D0-BETA*(DZFDX%R(I)*CTETA%R(I)+DZFDY%R(I)*STETA%R(I)) )
        ENDDO
! 
!     METHODE 2 : SOULSBY 1997 DYNAMICS OF MARINE SANDS p107-108 
!
      ELSEIF(SLOPEFF.EQ.2) THEN 
C 
        TANPHI = TAN(PHISED*PI/180.D0) 
C 
        DO I=1,NPOIN 
C
C         COSINUS ET SINUS DE LA DIRECTION DE LA PENTE
          DZF=SQRT(DZFDX%R(I)**2+DZFDY%R(I)**2)
          IF(DZF.GT.1.D-12) THEN
            CALPHA=DZFDX%R(I)/DZF
            SALPHA=DZFDY%R(I)/DZF
          ELSE
            CALPHA=1.D0
            SALPHA=0.D0
          ENDIF
C 
C         ZETA ANGLE QUE FAIT LA PENTE AVEC HORIZONTALE (BETA DE SOULSBY) 
          ZETA=ATAN(DZF)
          CZETA=COS(ZETA)
          SZETA=SIN(ZETA)    
C   
C         PSI ANGLE DU COURANT PAR RAPPORT A LA DIRECTION DE LA PENTE
C         PSI=TETA%R(I)-ALPHA
          CPSI=CTETA%R(I)*CALPHA+STETA%R(I)*SALPHA
          SPSI=STETA%R(I)*CALPHA-CTETA%R(I)*SALPHA  
          C1=(CZETA*TANPHI)**2-(SPSI*SZETA)**2 
          COEF%R(I)=MAX((CPSI*SZETA+SQRT(MAX(C1,0.D0)))/TANPHI,0.D0) 
          COEF%R(I)=MAX(COEF%R(I),0.D0) 
C 
        ENDDO 
! 
      ENDIF 
!
! ********************************************************************* ! 
!     V - TRAITEMENT DES POINTS FRONTIERES A DEBIT IMPOSE               ! 
!      PAS DE MODIFICATION DU QS QUAND IL EST DETERMINE PAR UTILISATEUR !  
! ********************************************************************* !
! 
      DO K = 1 , NPTFR 
         IF (LIQBOR%I(K) == KENT) COEF%R(MESH%NBOR%I(K)) = 1.D0 
!                           R.K. mai 2007 
!                           KSORT = 4
         IF (LIQBOR%I(K) == 4) COEF%R(MESH%NBOR%I(K)) = 1.D0 
      ENDDO 
! 
!====================================================================== 
!====================================================================== 
! 
      RETURN 
      END SUBROUTINE BEDLOAD_EFFPNT
