C                       ******************                                                       
                        SUBROUTINE INITSTR                                    
C                       ******************                                    
C     
     *(CHESTR,SETSTR,PZONE,NZONE,NPOIN,T1)
C
C***********************************************************************
C TELEMAC-2D VERSION 5.2         22/10/01  J-M HERVOUET TEL: 30 87 80 18                                                                  
C***********************************************************************
C     
C     FUNCTION : ASSIGN INITIAL VALUES OF STRICKLERS PER ZONE
C     
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
C     APPELE PAR :            HOMERE_PIT
C     
C     SOUS-PROGRAMME APPELE : OS
C     
C**********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(IN)      :: CHESTR
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: SETSTR,T1
      INTEGER, INTENT(IN)             :: PZONE(*)
      INTEGER, INTENT(IN)             :: NZONE
      INTEGER, INTENT(IN)             :: NPOIN
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,J
C        
C----------------------------------------------------------------------
C
      IF(NZONE.GT.0) THEN
C
C       ZONATION : SETSTR=AVERAGE PER ZONE OF CHESTR
C
        CALL OS('X=C     ',X=SETSTR,C=0.D0)
        CALL OS('X=Y     ',X=T1    ,Y=SETSTR)
        DO J=1,NZONE
          DO I=1,NPOIN
            IF(PZONE(I).EQ.J) THEN
              SETSTR%R(J)=SETSTR%R(J)+CHESTR%R(I)
              T1%R(J)=T1%R(J)+1.D0
            ENDIF
          ENDDO
          SETSTR%R(J)=SETSTR%R(J)/T1%R(J)
        ENDDO
C
      ELSE
C
C       NO ZONATION : SETSTR=CHESTR
C
        CALL OS('X=Y     ',X=SETSTR,Y=CHESTR)
C
      ENDIF
C        
C----------------------------------------------------------------------
C            
      RETURN
      END
