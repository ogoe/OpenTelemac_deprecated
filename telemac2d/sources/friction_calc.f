C                       ************************
                        SUBROUTINE FRICTION_CALC
C                       ************************
C
     &(N_START, N_END, KFROT, NDEF, VK, GRAV,
     & KARMAN, CHESTR, DW_MESH, HC, VRES, CF)
C
C***********************************************************************
C  TELEMAC-2D VERSION 5.5              J-M HERVOUET (LNH) 01 30 87 80 18
C                           20/04/04    F. HUVELIN 
C***********************************************************************
C
C FONCTION : SETTING THE FRICTION COEFFICIENT
C
C----------------------------------------------------------------------
C                             ARGUMENTS                                
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE                       
C |________________|____|______________________________________________
C | N_START,N_END  | => | Starting and ending point                    
C | KFROT          | => | Law used for the calculation                 
C | NDEF           | => | Default's Manning                            
C | VK             | => | Kinematic viscosity                          
C | GRAV           | => | Gravity acceleration                         
C | KARMAN         | => | Von Karman's constant                        
C | CHESTR         | => | Friction parameter                           
C | DW_MESH        | => | Distance to the boundary                    
C | HC             | => | Water depth : max(H,HMIN)                    
C | VRES           | => | Resultant velocity                           
C | CF             | <= | Friction coefficient (bottom or wall)        
C |________________|____|______________________________________________
C                    <=  input value                                   
C                    =>  output value                                   
C----------------------------------------------------------------------
C
C     FRICTION LAW PROGRAMMED :
C
C     KFROT = 0 :  NO FRICTION
C     KFROT = 1 :  LAW OF HAALAND
C     KFROT = 2 :  LAW OF CHEZY
C     KFROT = 3 :  LAW OF STRICKLER
C     KFROT = 4 :  LAW OF MANNING
C     KFROT = 5 :  LAW OF NIKURADSE
C     KFROT = 6 :  LOG LAW OF WALL
C     KFROT = 7 :  LAW OF COLEBROOK-WHITE
C
C
!=======================================================================!
!=======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS                !
!=======================================================================!
!=======================================================================!
!
      USE BIEF
!
      IMPLICIT NONE      
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,          INTENT(IN)    :: N_START, N_END, KFROT
      DOUBLE PRECISION, INTENT(IN)    :: NDEF, VK, GRAV, KARMAN
      TYPE(BIEF_OBJ),   INTENT(IN)    :: CHESTR, DW_MESH, HC, VRES
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: CF
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER                       :: I, ITER
      DOUBLE PRECISION, PARAMETER   :: TIERS = 1.D0/3.D0
      DOUBLE PRECISION              :: UNORM, INLOG, AUX
      DOUBLE PRECISION              :: OLDUST, OLDCF
      DOUBLE PRECISION              :: RE, UST, DW, DWPLUS
      DOUBLE PRECISION              :: TERM1, TERM2
!
!=======================================================================!
!=======================================================================!
!                               PROGRAMME                               !
!=======================================================================!
!=======================================================================!
!
      SELECT CASE (KFROT)
!
! NO FRICTION
! -----------
!
      CASE(0)
!
         DO I = N_START, N_END
            CF%R(I) = 0.D0
         ENDDO
!
! LAW OF HAALAND
! --------------
!
      CASE(1)
!
         DO I = N_START, N_END
            UNORM = MAX(VRES%R(I),1.D-6)
C                         1.D-6 : VISCOSITE LAMINAIRE DE L'EAU
            INLOG = (6.9D0*1.D-6/4.D0  /HC%R(I)/UNORM)**3
     &            + (CHESTR%R(I)/14.8D0/HC%R(I))**3.33
            INLOG = MIN(1.D0-1.D-6,INLOG)
            AUX   = -0.6D0*LOG(INLOG)/LOG(10.D0)
            CF%R(I) = 0.25D0 / AUX**2
         ENDDO
!
! LAW OF CHEZY
! ------------
!
      CASE(2)
!
         DO I = N_START, N_END
            CF%R(I) = 2.D0*GRAV/(CHESTR%R(I)**2)
         ENDDO
!
! LAW OF STRICKLER
! ----------------
!
      CASE(3)
!
         DO I = N_START, N_END
            CF%R(I) = 2.D0*GRAV/CHESTR%R(I)**2/HC%R(I)**TIERS
         ENDDO
!
! LAW OF MANNING
! --------------
!
      CASE(4)
!
         DO I = N_START, N_END
            CF%R(I) = 2.D0*GRAV*(CHESTR%R(I)**2)/HC%R(I)**TIERS
         ENDDO
!
! LAW OF NIKURADSE
! ----------------
!
      CASE(5)
!
         DO I = N_START, N_END
            CF%R(I) = 2.D0/(LOG( 11.D0*HC%R(I)/CHESTR%R(I))/KARMAN )**2
         ENDDO
!
! LOG LAW OF WALL
! ---------------
!
      CASE(6)
!
         DO I = N_START, N_END
!
            IF(VRES%R(I) < 1.0D-9) THEN
               CF%R(I) = 20.D0 ! Rismo2d = 10.D0 and Telemac2d = 2*10.D0
            ELSE
!
               DW = 0.33D0*DW_MESH%R(I)
!
               IF (CHESTR%R(I) < 1.0D-9) THEN
!
! iterative computation of friction velocity Ust
! ----------------------------------------------
!
                  UST    = 100.0*VK/DW
                  OLDUST = 0.D0
!
                  DO ITER = 1, 50
!
                     IF (ABS((UST-OLDUST)/UST)<=1.0D-6) EXIT
!
                     DWPLUS = DW*UST/VK
!
                     IF (DWPLUS < 11.D0) DWPLUS = 11.D0
!  
                     OLDUST = UST
                     UST    = KARMAN*VRES%R(I) / LOG(9.D0*DWPLUS)
! 
                  ENDDO
!
               ELSE
                  UST = KARMAN*VRES%R(I) / (LOG(DW/CHESTR%R(I))+8.5D0)
                  RE  = CHESTR%R(I)*UST  / VK
!
                  IF (RE < 70.D0) THEN
!
! iterative computation of friction velocity Ust
! ----------------------------------------------
!
                     OLDUST = 0.D0
!
                     DO ITER = 1, 50
!
                        IF (ABS((UST-OLDUST)/UST)<=1.0D-6) EXIT
!
                        DWPLUS = DW*UST/VK
!
                        IF (DWPLUS < 11.D0) DWPLUS = 11.D0
!
                        RE     = CHESTR%R(I)*UST/VK
                        OLDUST = UST
!         
                        IF (RE < 3.32D0) THEN
                           UST = KARMAN*VRES%R(I) / LOG(9.D0*DWPLUS)
                        ELSE
                           UST = KARMAN*VRES%R(I)
     &                         / (  LOG(DW/CHESTR%R(I))
     &                            + 3.32D0*LOG(RE)/RE
     &                            + KARMAN*(8.5D0-9.96D0/RE))
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
!
               DWPLUS = DW*UST/VK
!   
               IF (DWPLUS < 11.D0 ) THEN
                  UST = 11.0*VK / DW
                  UST = SQRT(VRES%R(I)*VK/DW)
               ENDIF
!
               CF%R(I) = 2.D0*(UST**2) / (VRES%R(I)**2)
            ENDIF
         ENDDO
!
! LAW OF COLEBROOK-WHITE
! ----------------------
!
      CASE(7)
!
         DO I = N_START, N_END
!
            RE = 4.D0*VRES%R(I)*HC%R(I)/VK
!
! the original condition for laminar/turbulent flow
! could not be hold (problems during NR iteration):
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
            IF(RE.GT.500.D0) THEN
!
! test condition: roughness less than flow depth
! ----------------------------------------------
!
               IF(CHESTR%R(I) < HC%R(I)) THEN
!                 NDEF : default Manning's n
                  CF%R(I) = 2.D0*(NDEF**2)*GRAV/HC%R(I)**TIERS
               ELSE
                  TERM1   = 4.4D0 / RE
                  TERM2   = CHESTR%R(I) / 14.84D0 / HC%R(I)
                  CF%R(I) = 2.5D0 ! initialize cf=1/sqrt(cf) for iteration
                  OLDCF   = 0.D0
!
                  DO ITER = 1, 50
                     IF (ABS((CF%R(I)-OLDCF)/CF%R(I))<=1.0D-6) EXIT
                     OLDCF = CF%R(I)
                     CF%R(I) = -2.03D0*LOG10(OLDCF*TERM1 + TERM2)
                  ENDDO
!        
                  IF (ITER.GE.50) THEN
                     CF%R(I) = -2.03D0*LOG10(TERM2)
                  ENDIF
!
                  CF%R(I) = 2.D0 / (CF%R(I)**2) / 8.D0
               ENDIF
!
            ELSEIF (RE.GT.100.D0 ) THEN
               CF%R(I) = 16.D0 / RE
            ELSE
               CF%R(I) = 0.16D0
            ENDIF
         ENDDO
!
! OTHER CASES
! -----------
!
      CASE DEFAULT
!
         IF (LNG.EQ.1) WRITE(LU,1) KFROT
1        FORMAT(I5,' : LOI DE FROTTEMENT INCONNUE')
         IF (LNG.EQ.2) WRITE(LU,2) KFROT
2        FORMAT(I5,' : UNKNOWN FRICTION LAW')
         CALL PLANTE(1)
         STOP
!
      END SELECT
!
!=======================================================================!
!=======================================================================!
!
      RETURN
      END
