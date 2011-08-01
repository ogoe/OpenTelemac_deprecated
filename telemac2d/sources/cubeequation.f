C                       ***********************
                        SUBROUTINE cubeEquation
C                       ***********************
C
     & (aCof, bCof, cCof, dCof, reals, x)
C
C***********************************************************************
C  TELEMAC-2D VERSION 5.5                 20/04/2004    F. HUVELIN 
C***********************************************************************
C                                                                      
C----------------------------------------------------------------------
C                             ARGUMENTS                                
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE                       
C |________________|____|______________________________________________
C | aCof           | => | Constant for X**3                            
C | bCof           | => | Constant for X**2                            
C | cCof           | => | Constant for X                               
C | dCof           | => | Constant of the equation                     
C | reals          | <= | Number of real roots                         
C | x              | <= | Value(s) of the real root(s)                 
C |________________|____|______________________________________________
C                    <=  output value                                   
C                    =>  input value                                   
C----------------------------------------------------------------------
C                                                                     
C CALLED BY FRICTION_LINDNER                                           
C                                                                                                                      !
C                                                                     
C======================================================================
C======================================================================
C                             EXPLANATION                              
C======================================================================
C======================================================================
C
C   Resolution of the equation of aX**3+bX**2+cX+d=0 (E1)                
C                                                                      
C   1/ Change the variable x by x = X -b/(3a)                            
C                              P = (c/a) - b*b/(3a*a)                    
C                              Q = 2b**3/(27a**3)+d/a-bc/(3a*a)         
C                                                                    
C      =>  (E1) <=> X**3+PX+Q = 0 (E2)                                  
C                                                                     
C   2/ Compute delta = (Q*Q)/4 + (P**3)/27                            
C                                                                      
C   3/ delta  > 0 : solution of Cardan      (1 real root)              
C                                                                      
C   4/ delta <= 0 : solution of Bombelli(?) (3 real roots)              
C                                                                      
C======================================================================
C======================================================================
C                    DECLARATION DES TYPES ET DIMENSIONS               
C======================================================================
C======================================================================
C
      IMPLICIT NONE
C       
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(IN)  :: aCof, bCof, cCof, dCof
      INTEGER,          INTENT(OUT) :: reals
      DOUBLE PRECISION, INTENT(OUT) :: x(3)
C       
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979323846D0
      DOUBLE PRECISION            :: ba, ca, P, Q, Q2P3, U, V
      DOUBLE PRECISION            :: expo, sign, tmp, phi
C
C=======================================================================
C======================================================================= 
C
      ba = bCof / aCof / 3.D0
      ca = cCof / aCof
C  
      P  = ca/3.D0 - ba**2
      Q  = ba**3 - ba*ca/2.D0 + dCof/aCof/2.D0
C  
      Q2P3 = Q**2 + P**3
C   
      IF (Q2P3 > 0.D0) THEN  
         reals = 1
         expo  = 1.D0/3.D0
         tmp   = -Q + SQRT(Q2P3)
         sign  = tmp / ABS(tmp)
         U     = sign * ABS(tmp)**(expo)
         tmp   = -Q - SQRT(Q2P3)
         sign  = tmp / ABS(tmp)
         V     = sign * ABS(tmp)**expo
         x(1)  = (U + V) - ba
C
      ELSE
C
        reals = 3
        tmp = -Q / (-P)**(1.5D0)
C
        IF (tmp >= 1.D0) THEN
           phi = 0.D0
        ELSE IF (tmp <= -1.D0) THEN
           phi = PI
        ELSE
          phi = acos (tmp)
        ENDIF
C
        x(1) = 2.D0* SQRT(-P)* COS(phi/3.D0)           -  ba
        x(2) = 2.D0* SQRT(-P)* COS((phi+2.D0*PI)/3.D0) -  ba
        x(3) = 2.D0* SQRT(-P)* COS((phi+4.D0*PI)/3.D0) -  ba
C
      ENDIF
C  
C=======================================================================
C=======================================================================
C
      RETURN
      END
