C                        ******************
                         SUBROUTINE COEFMAT
C                        ******************
C
     *(PERIAF,DT,M,AM,NPERIAF)
C     
C***********************************************************************
C  TELEMAC 2D VERSION 5.7    28/07/2006                        Chun WANG 
C    
C***********************************************************************
C
C      FONCTION:    Establish the coefficient matrice used for spectrum 
C                   analysis. The theory employed here is the Least Mean 
C                   Error Square method (refer to J. M. Janin et. al. 
C                   Simulation Des Courants de Maree en Manche et Proche 
C                   Atlantique,Page 27-28).
C
C      This subroutine was called by SPECTRE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !   PERIAF       !--> !  period of waves.
C !   AM           !<-->!  (2NPERIAF*2NPERIAF) coefficient matrix,  
C !   NPERIAF      !--> !  number of waves
C !   DT           !--> !  time interval.
C !   M            !--> !  number of sampling points
C !   W            !<-- !  Circular frequence.
C !________________!____!_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------     
C
      IMPLICIT NONE
      INTEGER LNG,LU             
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,          INTENT(IN   ) :: NPERIAF,M
      DOUBLE PRECISION, INTENT(IN   ) :: DT,PERIAF(NPERIAF)
      DOUBLE PRECISION, INTENT(INOUT) :: AM(2*NPERIAF,2*NPERIAF)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C                        500>NPERIAF
      DOUBLE PRECISION W(500),PI
      DOUBLE PRECISION A,B,C,D,AA,BB
      INTEGER I,J
C
      INTRINSIC COS,SIN,ACOS
C
C-----------------------------------------------------------------------     
C
      PI = ACOS(-1.D0)
C
      DO I = 1, NPERIAF
        W(I)=2.D0*PI/PERIAF(I)
      ENDDO
C      
      DO I = 1, NPERIAF
         DO J = 1, NPERIAF
            AA = (W(I)+W(J))*DT
            BB = (W(I)-W(J))*DT
            A = (-1.D0+COS(AA)+COS(M*AA)-COS((M+1)*AA))/(2-2*COS(AA))
            C = (SIN(AA)+SIN(M*AA)-SIN((M+1)*AA))/(2-2*COS(AA))
            IF(I.EQ.J) THEN
              B = M*1.D0
              D = 0.D0
            ELSE   
              B = (-1+COS(BB)+COS(M*BB)-COS((M+1)*BB))/(2-2*COS(BB))            
              D = (SIN(BB)+SIN(M*BB)-SIN((M+1)*BB))/(2-2*COS(BB))
            ENDIF   
            AM(I,J) = (A+B)/(2*M)
            AM(I,NPERIAF+J) = (C-D)/(2*M)
            AM(I+NPERIAF,J) = (C+D)/(2*M)
            AM(I+NPERIAF,NPERIAF+J) = (B-A)/(2*M)
         ENDDO
      ENDDO         
C
C-----------------------------------------------------------------------     
C
      RETURN
      END
