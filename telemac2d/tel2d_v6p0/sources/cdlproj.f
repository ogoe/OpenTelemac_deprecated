C                       ******************
                        SUBROUTINE CDLPROJ
C                       ******************
C
     *(NS,NPTFR,NBOR,LIMPRO,XNEBOR,YNEBOR,KNEU,UA)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.3               INRIA
C
C***********************************************************************
C
C     FONCTION  : PROJECTION OF THE SOLUTION ON THE BOUNDARY CONDITIONS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |                | -- |  
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NS,NPTFR,KNEU
      INTEGER, INTENT(IN) :: NBOR(NPTFR),LIMPRO(NPTFR,6)
      DOUBLE PRECISION, INTENT(IN) :: XNEBOR(2*NPTFR),YNEBOR(2*NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: UA(3,NS)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C       
      INTEGER IS,K  
C      
      DOUBLE PRECISION VNX,VNY,CEN
C                                                          
C-----------------------------------------------------------------------
C
C  CONDITION DE GLISSEMENT
C  ***********************
C
C
      DO K=1,NPTFR
C
       IF(LIMPRO(K,1).EQ.KNEU) THEN
C
         IS=NBOR(K)
         VNX=XNEBOR(K)
         VNY=YNEBOR(K)
         CEN=UA(2,IS)*VNX+UA(3,IS)*VNY
         UA(2,IS) = UA(2,IS) -CEN*VNX  
         UA(3,IS) = UA(3,IS) -CEN*VNY
C
       ENDIF
C
      ENDDO  
C                                                          
C-----------------------------------------------------------------------
C
      RETURN
      END
