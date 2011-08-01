C                       *****************                          
                        SUBROUTINE NEWSTR                                    
C                       *****************                                    
C     
     *(SETSTR,SETSTR2,DESC,RO,RSTART,NPARAM,ESTIME,KFROT)
C
C***********************************************************************
C TELEMAC-2D VERSION 5.2  22/10/01  J-M HERVOUET TEL: 30 87 80 18
C                         22/03/94  E. BARROS   (IN TELEMAC 2.2)
C                         02/10/00  A. LEOPARDI (UNINA) (UPGRADE TO 5.1)                                                                      
C***********************************************************************
C     
C     FUNCTION : COMPUTE THE NEW SET OF FRICTION COEFFICIENTS 
C     
C-----------------------------------------------------------------------
C     ARGUMENTS                                         
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      SETSTR    | <--|  NEW SET
C |      SETSTR2   | -->|  OLD SET
C |      DESC      | -->|  
C |      RO        | -->|  SETSTR=SETSTR2+RO*DESC
C |      RSTART    | -->|  LOGICAL, RESTART COMPUTATION BECAUSE OUT OF LIMITS
C |      NZONE     | -->|  NUMBER OF ZONES
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C     
C     APPELE PAR :            HOMERE_PIT
C     
C     SOUS-PROGRAMME APPELE : OS
C     
C**********************************************************************
C
      USE BIEF
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION , INTENT(IN)    :: RO
      TYPE (BIEF_OBJ)  , INTENT(IN)    :: DESC
      TYPE (BIEF_OBJ)  , INTENT(IN)    :: SETSTR2
      TYPE (BIEF_OBJ)  , INTENT(INOUT) :: SETSTR
      LOGICAL          , INTENT(INOUT) :: RSTART
      INTEGER          , INTENT(IN)    :: NPARAM,KFROT
      CHARACTER(LEN=72), INTENT(IN)    :: ESTIME
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C        
C---------------------------------------------------------------------
C
      CALL OV( 'X=Y+CZ  ',SETSTR%R,SETSTR2%R,DESC%R,RO,NPARAM)
C
C     TEST ON LIMITS
C     LIMITS (1,100)
C
      RSTART=.FALSE.
      DO I=1,NPARAM 
         IF (SETSTR%R(I).LT.1.D0) THEN
           SETSTR%R(I)=1.D0
           RSTART=.TRUE.
         ENDIF        
         IF (SETSTR%R(I).GT.100.D0) THEN
           SETSTR%R(I)=100.D0
           RSTART=.TRUE.
         ENDIF        
      ENDDO
C        
C---------------------------------------------------------------------
C            
      RETURN
      END
