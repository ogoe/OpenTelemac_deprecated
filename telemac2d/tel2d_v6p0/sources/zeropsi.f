C                       ******************
                        SUBROUTINE ZEROPSI
C                       ******************
C
     *(X0,X,NIT,CA1,A2)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.4                                           INRIA
C
C***********************************************************************
C
C     FONCTION  : ZERO DE PSI-A2 PAR LA METHODE DE NEWTON
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
      DOUBLE PRECISION, INTENT(IN)    :: X0,A2,CA1
      DOUBLE PRECISION, INTENT(INOUT) :: X
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NITEPS 
C
      DOUBLE PRECISION EPS,EPSX,SQ32,SQ3I
      DOUBLE PRECISION AMPLUS,AMMOINS,PHI1,PHI,DPHIPH3,CAPHI,FDF
C                                                          
C-----------------------------------------------------------------------
C
      SQ32=SQRT(1.5D0)
      SQ3I=1.D0/SQRT(3.D0)
      EPS=1.E-12
      EPSX=1.E-12
      NIT=0
      NITEPS=0
      X=X0
C
1     NIT=NIT+1
C
      IF(X.LE.-SQ32) THEN
      X= -SQ32 + EPSX
      NITEPS= NITEPS +1
      ENDIF
      IF(NITEPS.EQ.3) THEN
      X= -SQ32
      GOTO 10
      ENDIF
C     
      AMPLUS=MAX(-X,+SQ32)
      AMMOINS=MAX(-X,-SQ32)
C
      PHI1=0.5D0*(AMPLUS+AMMOINS+2.*X)
      PHI = PHI1*SQ3I*(AMPLUS-AMMOINS)
      DPHIPH3 = 1.D0/(3.D0*PHI1)
C
      CAPHI = CA1* PHI**(1.D0/3.D0)
      FDF= (X - 2.D0 - A2 *CAPHI)/(1.D0- DPHIPH3*(X-2.D0))
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
