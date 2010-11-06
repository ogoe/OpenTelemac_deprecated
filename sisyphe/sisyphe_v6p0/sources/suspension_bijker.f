
      ! *****************************     !
        SUBROUTINE SUSPENSION_BIJKER      ! 
      ! *****************************     !

     &  (TAUP, HN, NPOIN, CHARR, QSC, ZREF,
     &    ZERO, HMIN, CSTAEQ,XMVE)


C**********************************************************************C
C SISYPHE VERSION 5.6  04/01/05  F. HUVELIN                            C
C SISYPHE VERSION 5.5  14/04/04  C. VILLARET  01 30 87 83 28           C
C**********************************************************************C


         ! ==================================================== !
         !   Reference concentration calculation at z= 2*d50    !
         ! thanks to the formula of Zyserman and Fredsoe (1994) !
         ! ==================================================== !


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
C |   ACLADM       | => |
C |   CF           | => |
C |   TOB          | => |
C |   HCLIP        | => |
C |   AVA          | => |
C |   NPOIN        | => |
C |   CHARR        | => |
C |   KSPRATIO     | => |
C |   HMIN         | => |
C |   GRAV         | => |
C |   XMVE         | => |
C |   XMVS         | => |
C |   AC           | <=>|
C |   FLUER        | <= |
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
      USE INTERFACE_SISYPHE,EX_SUSPENSION_BIJKER => SUSPENSION_BIJKER
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU


      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE(BIEF_OBJ),   INTENT(IN)    :: TAUP,QSC,HN
      TYPE (BIEF_OBJ),  INTENT(IN)    :: ZREF
      INTEGER,          INTENT(IN)    :: NPOIN
      LOGICAL,          INTENT(IN)    :: CHARR
      DOUBLE PRECISION, INTENT(IN)    :: HMIN,ZERO,XMVE
C
      TYPE(BIEF_OBJ),   INTENT(INOUT)   ::  CSTAEQ


      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER                     :: I
      DOUBLE PRECISION            :: USTARP,CMAX
C
C     MAXIMUM CONCENTRATION CORRESPONDING TO DENSE PACKING
C
      DATA CMAX/0.6D0/
C
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!

      IF(.NOT.CHARR) THEN
        WRITE(LU,*) 'SUSPENSION_BIJKER ERROR ON CHARR'
        CALL PLANTE(1)
        STOP
      ENDIF     
C
      DO I=1,NPOIN
C           
        IF(TAUP%R(I).LE.ZERO) THEN
          CSTAEQ%R(I) = 0.D0
        ELSE
          USTARP=SQRT(TAUP%R(I)/XMVE)        
          CSTAEQ%R(I) = QSC%R(I)/(6.34D0*USTARP*ZREF%R(I))
          CSTAEQ%R(I) = MIN(CSTAEQ%R(I),CMAX)
        ENDIF
C     
      ENDDO

!======================================================================!
!======================================================================!

      RETURN      
      END SUBROUTINE SUSPENSION_BIJKER
