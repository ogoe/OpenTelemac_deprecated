C                          *****************                                                                      
                           SUBROUTINE SD_NNF
C                          *****************
C
     *(N,R,C,IC,IA,JA,A,Z,B,IL,JL,L,LMAX,D,IU,JU,U,UMAX,ROW,TMP,FLAG)                                             
C 
C***********************************************************************
C BIEF VERSION 5.9     18/02/08   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C 
C  FUNCTION: NUMERIC LDU-FACTORIZATION OF SPARSE NONSYMMETRIC MATRIX AND        
C            SOLUTION OF SYSTEM OF LINEAR EQUATIONS (UNCOMPRESSED POINTER     
C            STORAGE)
C                                                                     
C       INPUT VARIABLES:   N, R,C,IC, IA,JA,A, B, IL,JL,LMAX, IU,JU,UMAX
C       OUTPUT VARIABLES:  Z, L,D,U, FLAG                               
C                                                                       
C       PARAMETERS USED INTERNALLY:                                     
C FIA   \ ROW   - HOLDS INTERMEDIATE VALUES IN CALCULATION OF L, D, U.  
C       \           SIZE = N.                                           
C FIA   \ TMP   - HOLDS NEW RIGHT-HAND SIDE B' FOR SOLUTION OF THE      
C       \           EQUATION  UX = B'.                                  
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
      INTEGER R(*),C(*),IC(*),IA(*),JA(*),N                      
      INTEGER IL(*),JL(*),LMAX,IU(*),JU(*),UMAX,FLAG                                             
      DOUBLE PRECISION A(*),Z(*),B(*),L(*),D(*),U(*),ROW(*),TMP(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+  
C
      INTEGER LI,I,J,K,IMIN,IMAX,JMIN,JMAX
      DOUBLE PRECISION BSUM,DK
C
C-----------------------------------------------------------------------
C                                                                                                        
C     CHECK STORAGE 
C
      IF(IL(N+1)-1.GT.LMAX) THEN
C       ERROR: INSUFFICIENT STORAGE FOR L
        FLAG = 4*N + 1                                                  
        RETURN
      ENDIF                             
      IF(IU(N+1)-1.GT.UMAX) THEN
C       ERROR: INSUFFICIENT STORAGE FOR U                                 
        FLAG = 7*N + 1                                                  
        RETURN  
      ENDIF                             
C                                                                       
C     FOR EACH ROW 
C
      DO 10 K=1,N
C                                                     
C       SET THE INITIAL STRUCTURE OF ROW
C
        JMIN = IL(K)                                                  
        JMAX = IL(K+1) - 1                                            
        IF(JMAX.GE.JMIN) THEN                                    
C         IF L(K,M) .NE. 0, ROW(M)=0
          DO J=JMIN,JMAX                                              
            ROW(JL(J)) = 0 
          ENDDO
        ENDIF                                             
        ROW(K) = 0                                                    
        JMIN = IU(K)                                                  
        JMAX = IU(K+1) - 1                                            
        IF(JMAX.GE.JMIN) THEN                                    
C         IF U(K,M) .NE. 0, ROW(M)=0
          DO J=JMIN,JMAX                                              
            ROW(JU(J)) = 0
          ENDDO 
        ENDIF                                             
        JMIN = IA(R(K))                                               
        JMAX = IA(R(K)+1) - 1                                         
C       SET ROW TO KTH ROW OF REORDERED A
        DO J=JMIN,JMAX                                              
          ROW(IC(JA(J))) = A(J)
        ENDDO                                       
C       INITIALIZE BSUM
        BSUM = B(R(K))                                                 
C                                                                       
C       ASSIGN THE KTH ROW OF L AND ADJUST ROW, BSUM  
        IMIN = IL(K)                                                  
        IMAX = IL(K+1) - 1                                            
        IF(IMAX.GT.IMIN) THEN                                    
          DO I=IMIN,IMAX                                              
            LI = - ROW(JL(I))                                           
C           IF L IS NOT REQUIRED, THEN COMMENT OUT THE FOLLOWING LINE
            L(I) = - LI                                                 
            BSUM = BSUM + LI * TMP(JL(I))                                 
            JMIN = IU(JL(I))                                            
            JMAX = IU(JL(I)+1) - 1                                      
            IF(JMAX.GT.JMIN) THEN                                  
              DO J=JMIN,JMAX                                            
                ROW(JU(J)) = ROW(JU(J)) + LI * U(J) 
              ENDDO 
            ENDIF                     
          ENDDO
        ENDIF                                                    
C                                                                       
C       ASSIGN DIAGONAL D AND KTH ROW OF U, SET TMP(K)
C
        IF(ROW(K).EQ.0) THEN
C         ERROR:  ZERO PIVOT                                                 
          FLAG = 8*N + K                                                  
          RETURN   
        ENDIF                                   
        DK = 1 / ROW(K)                                               
        D(K) = DK                                                     
        TMP(K) = BSUM * DK                                             
        JMIN = IU(K)                                                  
        JMAX = IU(K+1) - 1                                            
        IF(JMAX.GE.JMIN) THEN                                   
          DO J=JMIN,JMAX                                              
            U(J) = ROW(JU(J)) * DK 
          ENDDO 
        ENDIF                                    
10    CONTINUE                                                      
C                                                                       
C     SOLVE  UX = TMP  BY BACK SUBSTITUTION
C
      K = N                                                           
      DO I=1,N                                                     
        BSUM = TMP(K)                                                  
        JMIN = IU(K)                                                  
        JMAX = IU(K+1) - 1                                            
        IF(JMAX.GE.JMIN) THEN                                   
          DO J=JMIN,JMAX                                             
            BSUM = BSUM - U(J) * TMP(JU(J))
          ENDDO 
        ENDIF                              
        TMP(K)  = BSUM                                                  
        Z(C(K)) = BSUM                                                 
        K = K-1
      ENDDO                                                       
C                                                                       
      FLAG = 0
C
C-----------------------------------------------------------------------
C                                                        
      RETURN                                                                                                                 
      END
