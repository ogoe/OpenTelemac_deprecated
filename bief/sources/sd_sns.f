C                            *****************
                             SUBROUTINE SD_SNS
C                            *****************
C
     *(N,P,D,IJU,JU,IU,U,Z,B,TMP)
C
C***********************************************************************
C BIEF VERSION 5.9     13/11/08   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C
C  FONCTION :  SOLUTION OF SPARSE SYMMETRIC POSITIVE DEFINITE SYSTEM OF      
C              LINEAR EQUATIONS  MX = B  GIVEN UT-D-U FACTORIZATION OF M              
C 
C  COMPRESSED STORAGE OF SPARSE MATRICES                                
C                                                                       
C    THE STRICT UPPER TRIANGULAR PORTION OF THE MATRIX U IS STORED IN   
C    (IA,JA,A) FORMAT USING THE ARRAYS IU, JU, AND U, EXCEPT THAT AN    
C    ADDITIONAL ARRAY IJU IS USED TO REDUCE THE STORAGE REQUIRED FOR JU 
C    BY ALLOWING SOME SEQUENCES OF COLUMN INDICES TO CORRESPOND TO MORE 
C    THAN ONE ROW.  FOR I < N, IJU(I) IS THE INDEX IN JU OF THE FIRST   
C    ENTRY FOR THE I-TH ROW;  IJU(N) IS THE NUMBER OF ENTRIES IN JU.    
C    THUS, THE NUMBER OF ENTRIES IN THE I-TH ROW IS  IU(I+1) - IU(I),   
C    THE NONZERO ENTRIES OF THE I-TH ROW ARE STORED CONSECUTIVELY IN    
C                                                                       
C        U(IU(I)),   U(IU(I)+1),   ..., U(IU(I+1)-1),                   
C                                                                       
C    AND THE CORRESPONDING COLUMN INDICES ARE STORED CONSECUTIVELY IN   
C                                                                       
C        JU(IJU(I)), JU(IJU(I)+1), ..., JU(IJU(I)+IU(I+1)-IU(I)-1).     
C                                                                       
C    COMPRESSION IN JU OCCURS IN TWO WAYS.  FIRST, IF A ROW I WAS MERGED
C    INTO ROW K, AND THE NUMBER OF ELEMENTS MERGED IN FROM (THE TAIL    
C    PORTION OF) ROW I IS THE SAME AS THE FINAL LENGTH OF ROW K, THEN   
C    THE KTH ROW AND THE TAIL OF ROW I ARE IDENTICAL AND IJU(K) POINTS  
C    TO THE START OF THE TAIL.  SECOND, IF SOME TAIL PORTION OF THE     
C    (K-1)ST ROW IS IDENTICAL TO THE HEAD OF THE KTH ROW, THEN IJU(K)   
C    POINTS TO THE START OF THAT TAIL PORTION.  FOR EXAMPLE, THE NONZERO
C    STRUCTURE OF THE STRICT UPPER TRIANGULAR PART OF THE MATRIX        
C                                                                       
C             ( D 0 0 0 X X X )                                         
C             ( 0 D 0 X X 0 0 )                                         
C             ( 0 0 D 0 X X 0 )                                         
C         U = ( 0 0 0 D X X 0 )                                         
C             ( 0 0 0 0 D X X )                                         
C             ( 0 0 0 0 0 D X )                                         
C             ( 0 0 0 0 0 0 D )                                         
C                                                                       
C    WOULD BE STORED AS                                                 
C                                                                       
C             \ 1  2  3  4  5  6  7  8                                  
C         ----+------------------------                                 
C          IU \ 1  4  6  8 10 12 13 13                                  
C          JU \ 5  6  7  4  5  6                                        
C         IJU \ 1  4  5  5  2  3  6           .                         
C                                                                       
C    THE DIAGONAL ENTRIES OF U ARE EQUAL TO ONE AND ARE NOT STORED.     
C                                                                    
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C | N              | -->| ORDER OF THE MATRIX
C | P              | -->| INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN 
C |                |    | THE PERMUTATION OF THE ROWS AND COLUMNS OF M 
C |                |    | CORRESPONDING TO THE MINIMUM DEGREE ORDERING;  
C |                |    | DIMENSION = N 
C | D              |    | REAL ONE-DIMENSIONAL ARRAY CONTAINING THE 
C |                |    | RECIPROCALS OF THE DIAGONAL ENTRIES OF THE 
C |                |    | MATRIX D;  DIMENSION = N  
C | IJU            | -->| INTEGER ONE-DIMENSIONAL ARRAY CONTAINING 
C |                |    | POINTERS TO THE START OF EACH ROW IN JU;  DIMENSION = N 
C | JU             | -->| INTEGER ONE-DIMENSIONAL ARRAY CONTAINING THE 
C |                |    | COLUMN INDICES CORRESPONDING TO THE ELEMENTS 
C |                |    | OF U;  DIMENSION = JUMAX
C | IU             | -->| INTEGER ONE-DIMENSIONAL ARRAY CONTAINING 
C |                |    | POINTERS TO DELIMIT ROWS IN U;  DIMENSION = N+1
C | U              |    | REAL ONE-DIMENSIONAL ARRAY CONTAINING THE 
C |                |    | NONZERO ENTRIES IN THE STRICT UPPER TRIANGLE
C |                |    | OF U, STORED BY ROWS; DIMENSION = UMAX 
C | Z              |    | REAL ONE-DIMENSIONAL ARRAY CONTAINING THE 
C |                |    | SOLUTION X;  Z AND B CAN BE THE SAME ARRAY;  
C |                |    | DIMENSION = N 
C | B              |    | REAL ONE-DIMENSIONAL ARRAY CONTAINING THE 
C |                |    | RIGHT-HAND SIDE B; B AND Z CAN BE THE SAME ARRAY;  
C |                |    | DIMENSION = N 
C | TMP            | -- | REAL ONE-DIMENSIONAL WORK ARRAY; DIMENSION N
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_SD_SNS => SD_SNS
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: N
      INTEGER, INTENT(INOUT)          :: P(N),IJU(*),JU(*),IU(N+1)
      DOUBLE PRECISION, INTENT(IN)    :: B(N)
      DOUBLE PRECISION, INTENT(INOUT) :: TMP(N),Z(N),D(N),U(*)
C                                                           UMAX
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,J,K,JMIN,JMAX,MU                                                                    
      DOUBLE PRECISION TMPK,SU             
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C----SET TMP TO PERMUTED B
C                                              
      DO K=1,N                                                      
        TMP(K) = B(P(K))
      ENDDO                                              
C                                                                       
C----SOLVE  UT D Y = B  BY FORWARD SUBSTITUTION
C                         
      DO K=1,N                                                      
        TMPK = TMP(K)                                                 
        JMIN = IU(K)                                                  
        JMAX = IU(K+1) - 1                                            
        IF(JMIN.GT.JMAX)  GO TO 3                                    
        MU = IJU(K) - JMIN                                            
        DO J=JMIN,JMAX                                              
          TMP(JU(MU+J)) = TMP(JU(MU+J)) + U(J) * TMPK
        ENDDO                 
3       TMP(K) = TMPK * D(K) 
      ENDDO                                           
C                                                                       
C----SOLVE  U X = Y  BY BACK SUBSTITUTION
C                               
      K = N                                                           
      DO I=1,N                                                      
        SU   = TMP(K)                                                  
        JMIN = IU(K)                                                  
        JMAX = IU(K+1) - 1                                            
        IF(JMIN.GT.JMAX) GO TO 5                                    
        MU = IJU(K) - JMIN                                            
        DO J=JMIN,JMAX                                              
          SU = SU + U(J) * TMP(JU(MU+J))
        ENDDO                            
5       TMP(K) = SU                                                  
        Z(P(K)) = SU                                                 
        K = K-1
      ENDDO                                                       
C                                                                       
C-----------------------------------------------------------------------
C                                                                        
      RETURN                                                          
      END
