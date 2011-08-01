C                       ********************                                                        
                        SUBROUTINE ASSIGNSTR                                    
C                       ********************                                    
C     
     *(CHESTR,SETSTR,PZONE,NZONE,NPOIN)
C
C***********************************************************************
C TELEMAC-2D VERSION 5.2         22/10/01  J-M HERVOUET TEL: 30 87 80 18
C                                08/11/00  A. LEOPARDI (UNINA)                                                                       
C***********************************************************************
C     
C     FUNCTION : ASSIGN VALUES OF STRICKLERS
C          
C-----------------------------------------------------------------------
C     ARGUMENTS                                         
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      CHESTR    | <--|  STRICKLERS PER POINTS
C |      SETSTR    | -->|  SET OF STRICKLERS' (ZONES)
C |      PZONE     | -->|  TABLE OF ZONES
C |      NZONE     | -->|  NUMBER OF ZONES
C |      NPOIN     | -->|  NUMBER OF POINTS
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C     
C     APPELE PAR :  HOMERE_TELEMAC2D
C     
C     SOUS-PROGRAMME APPELE : OS
C     
C**********************************************************************
C
      USE BIEF
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: CHESTR 
      TYPE(BIEF_OBJ), INTENT(IN)      :: SETSTR
      INTEGER, INTENT(IN)             :: PZONE(*)
      INTEGER, INTENT(IN)             :: NZONE
      INTEGER, INTENT(IN)             :: NPOIN
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
      INTEGER I,J
C        
C---------------------------------------------------------------------
C
      IF(NZONE.GT.0) THEN
        DO J=1,NZONE
           DO I=1,NPOIN
              IF (PZONE(I).EQ.J) CHESTR%R(I)=SETSTR%R(J)
           ENDDO
        ENDDO
      ELSE
        CALL OS('X=Y     ',X=CHESTR,Y=SETSTR)
      ENDIF
C        
C-----------------------------------------------------------------------
C            
      RETURN
      END
