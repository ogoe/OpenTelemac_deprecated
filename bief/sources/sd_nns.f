C                          *****************                                                                      
                           SUBROUTINE SD_NNS
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
C  FUNCTION: NUMERIC SOLUTION OF A SPARSE NONSYMMETRIC SYSTEM OF LINEAR         
C            EQUATIONS GIVEN LDU-FACTORIZATION
C            (UNCOMPRESSED POINTER STORAGE)
C                                                                       
C       INPUT VARIABLES:   N, R,C, IL,JL,L, D, IU,JU,U, B               
C       OUTPUT VARIABLES:  Z                                            
C                                                                       
C       PARAMETERS USED INTERNALLY:                                     
C FIA   \ TMP   - HOLDS NEW RIGHT-HAND SIDE B' FOR SOLUTION OF THE      
C       \           EQUATION UX = B'.                                   
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
      INTEGER N,R(*),C(*),IL(*),JL(*),IU(*),JU(*)                             
      DOUBLE PRECISION L(*),D(*),U(*),Z(*),B(*),TMP(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,J,K,JMIN,JMAX
      DOUBLE PRECISION BSUM    
C
C-----------------------------------------------------------------------
C                                                                       
C  ******  SOLVE LDY = B  BY FORWARD SUBSTITUTION  *********************
C
        DO 2 K=1,N                                                      
          BSUM = B(R(K))                                                 
          JMIN = IL(K)                                                  
          JMAX = IL(K+1) - 1                                            
          IF (JMIN.GT.JMAX)  GO TO 2                                    
          DO 1 J=JMIN,JMAX                                              
            BSUM = BSUM - L(J) * TMP(JL(J))
1         CONTINUE                               
2         TMP(K) = BSUM * D(K)                                           
C                                                                       
C  ******  SOLVE  UX = Y  BY BACK SUBSTITUTION  ************************
C
      K = N                                                           
      DO 5 I=1,N                                                      
        BSUM = TMP(K)                                                  
        JMIN = IU(K)                                                  
        JMAX = IU(K+1) - 1                                            
        IF (JMIN.GT.JMAX)  GO TO 4                                    
        DO 3 J=JMIN,JMAX                                              
          BSUM = BSUM - U(J) * TMP(JU(J)) 
3       CONTINUE                              
   4    TMP(K) = BSUM                                                  
        Z(C(K)) = BSUM                                                 
        K = K-1 
5     CONTINUE 
C
C------------------------------------------------------------------------
C                                                     
      RETURN                                                          
      END
