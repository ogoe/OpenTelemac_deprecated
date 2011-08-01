

      ! *************************** !
        SUBROUTINE SUSPENSION_ROUSE 
      ! *************************** !

     &(USTAR,HN,NPOIN,KARMAN,HMIN,ZERO,XWC,ZREF,T2)
     
C**********************************************************************C
C SISYPHE VERSION 5.6  04/01/05  F. HUVELIN                            C
C SISYPHE VERSION 5.5  14/04/04  C. VILLARET   01 30 87 83 28          C
C SISYPHE VERSION 5.5  14/04/04  J-M HERVOUET  01 30 87 80 18          C
C SISYPHE VERSION 5.7  13/07/07  J-M HERVOUET  01 30 87 80 18          C
C**********************************************************************C

            ! ============================================ !
            !    Computation of the deposition flux and    !
            ! concentration thanks to the profile of Rouse !
            ! ============================================ !

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
C |   TOB          | => |
C |   CF           | => |
C |   HN           | => |
C |   Q            | => |
C |   CS           | => |
C |   NPOIN        | => |
C |   ZERO         | => |
C |   XMVE         | => |
C |   KARMAN       | => |
C |   HMIN         | => |
C |   XWC          | => |
C |   T1..T2       | <=>|  Work array                                   
C |   USTAR        | <=>| (Work array)
C |   HCHAR        | <=>| (Work array)
C |   IEIN         | <=>| (Work array)
C |   FLUDP        | <= |
C !________________|____|______________________________________________C
C                    <=  Can't be change by the user                   C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY SUSPENSION_FLUX                                            !
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
      USE INTERFACE_SISYPHE,EX_SUSPENSION_ROUSE => SUSPENSION_ROUSE
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU

      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE(BIEF_OBJ),   INTENT(IN)    :: USTAR,HN,ZREF
      INTEGER,          INTENT(IN)    :: NPOIN
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN,XWC,HMIN,ZERO
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: T2
     
      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER          :: I
      DOUBLE PRECISION :: B,EXP,ROUSE
!      
      INTRINSIC MAX,MIN,LOG
!     
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
!     ROUSE NUMBER AND MINIMUM BOUND OF THE EINSTEIN INTEGRAL  
!      
      DO I=1,NPOIN
!
!        Rouse number 
!
         ROUSE=  XWC / (KARMAN*MAX(USTAR%R(I),ZERO))
!
!        Minimum bound of the Einstein integral -->  B = KS/H 
! 
!        MODIF JMH 16/12/2009    B ALWAYS < 1 ? 
!        B = ZREF%R(I)/MAX(HN%R(I),HMIN)
         B = ZREF%R(I)/MAX(HN%R(I),ZREF%R(I))
!
!        RATIO BETWEEN reference CONC ON BOTTOM and MEAN CONC     
!        Assuming EXPONENTIAL PROFILE WITH EXPONENT ROUSE NUMBER --> T2
!     
         EXP=ROUSE-1.D0
         IF(ABS(EXP).GT.1.D-4) THEN
!          AJOUTE PAR JMH 12/07/2007 
           EXP=MIN(EXP,3.D0)
           T2%R(I)=B*(1.D0-B**EXP)/EXP
         ELSE
           T2%R(I)=-B*LOG(B)
         ENDIF 
         T2%R(I)=MAX(1.D0/MAX(T2%R(I),ZERO),1.D0) 
      ENDDO
!    
!======================================================================!
!======================================================================!
!
      RETURN     
      END
