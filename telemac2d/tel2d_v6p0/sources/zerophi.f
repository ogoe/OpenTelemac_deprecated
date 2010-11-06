C                       ******************
                        SUBROUTINE ZEROPHI
C                       ******************
C
     *(X0,X,NIT,CA1)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.4                                           INRIA
C
C***********************************************************************
C
C FONCTION  : ZERO DE PHI-CA1 PAR LA METHODE DE NEWTON
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
      INTEGER, INTENT(INOUT)          :: NIT
      DOUBLE PRECISION, INTENT(IN)    :: X0,CA1
      DOUBLE PRECISION, INTENT(INOUT) :: X
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NITEPS 
C
      DOUBLE PRECISION EPS,EPSX,SQ32,SQ3,AMPLUS,AMMOINS,DPHI,FDF
C                                                          
C-----------------------------------------------------------------------
C
      SQ3 =SQRT(3.D0)
      SQ32=SQRT(1.5D0)
      EPS=1.E-12
      EPSX=1.E-12
      NIT=0
      NITEPS=0
      X=X0
C
1     NIT=NIT+1
C
      IF(X.GE.SQ32) THEN
      X= SQ32 - EPSX
      NITEPS= NITEPS + 1
      ENDIF
      IF(NITEPS.EQ.3) THEN
      X = SQ32
      GOTO 10
      ENDIF
C
      AMPLUS =MIN(-X,+SQ32)
      AMMOINS=MIN(-X,-SQ32)
C
      DPHI= AMPLUS-AMMOINS
      FDF = (DPHI*(AMPLUS+AMMOINS+2.D0*X)-2.D0*SQ3*CA1)/(2.D0*DPHI)
C
      IF(ABS(FDF).LT.EPS) GOTO 10
      IF(NIT.EQ.20) GOTO 10
C
      X=X-FDF
      GOTO 1
C
10    CONTINUE
C                                                          
C-----------------------------------------------------------------------
C
      RETURN
      END
