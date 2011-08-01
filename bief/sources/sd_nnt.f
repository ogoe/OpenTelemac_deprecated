C                          *****************                                                                           
                           SUBROUTINE SD_NNT
C                          *****************
C    
     *(N,R,C,IL,JL,L,D,IU,JU,U,Z,B,TMP)                     
C
C***********************************************************************
C BIEF VERSION 5.9     18/02/08   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C
C  FUNCTION: NUMERIC SOLUTION OF THE TRANSPOSE OF A SPARSE NONSYMMETRIC
C            SYSTEM OF LINEAR EQUATIONS GIVEN LDU-FACTORIZATION
C            (UNCOMPRESSED POINTER STORAGE)
C                                                                       
C       INPUT VARIABLES:   N, R,C, IL,JL,L, D, IU,JU,U, B               
C       OUTPUT VARIABLES:  Z                                            
C                                                                       
C       PARAMETERS USED INTERNALLY:                                     
C FIA   \ TMP   - HOLDS NEW RIGHT-HAND SIDE B' FOR SOLUTION OF THE      
C       \           EQUATION LX = B'.                                   
C       \           SIZE = N.                                           
C 
C                          
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |                | -- |        
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************            
C 
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                                                      
      INTEGER R(*),C(*),IL(*),JL(*),IU(*),JU(*),N                           
      DOUBLE PRECISION L(*),D(*),U(*),Z(*),B(*),TMP(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,J,K,JMIN,JMAX
      DOUBLE PRECISION TMPK 
C
C-----------------------------------------------------------------------
C                                                                       
C  ******  SOLVE  UT Y = B  BY FORWARD SUBSTITUTION  *******************
C
      DO K=1,N                                                      
        TMP(K) = B(C(K))
      ENDDO
C                                              
      DO 3 K=1,N                                                      
        TMPK = - TMP(K)                                               
        JMIN = IU(K)                                                  
        JMAX = IU(K+1) - 1                                            
        IF (JMIN.GT.JMAX)  GO TO 3                                    
        DO 2 J=JMIN,JMAX                                              
          TMP(JU(J)) = TMP(JU(J)) + U(J) * TMPK 
2       CONTINUE                      
3     CONTINUE                                                      
C                                                                       
C  ******  SOLVE  D LT X = Y  BY BACK SUBSTITUTION  ********************
C
      K = N                                                           
      DO I=1,N                                                      
        TMPK = - TMP(K) * D(K)                                        
        JMIN = IL(K)                                                  
        JMAX = IL(K+1) - 1                                            
        IF(JMIN.GT.JMAX) GO TO 5                                    
        DO 4 J=JMIN,JMAX                                              
          TMP(JL(J)) = TMP(JL(J)) + L(J) * TMPK
4       CONTINUE                       
5       Z(R(K)) = - TMPK                                              
        K = K-1
      ENDDO 
C
C-----------------------------------------------------------------------
C                                                       
      RETURN                                                          
      END
