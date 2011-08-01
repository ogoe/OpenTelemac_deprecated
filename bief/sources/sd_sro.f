C                              *****************
                               SUBROUTINE SD_SRO
C                              *****************
C
     *(N,IP,IA,JA,A,Q,R,DFLAG)
C
C***********************************************************************
C BIEF VERSION 5.7     20/11/06   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C
C  FONCTION : SYMMETRIC REORDERING OF SPARSE SYMMETRIC MATRIX 
C
C    THE NONZERO ENTRIES OF THE MATRIX M ARE ASSUMED TO BE STORED       
C    SYMMETRICALLY IN (IA,JA,A) FORMAT (I.E., NOT BOTH M(I,J) AND M(J,I)
C    ARE STORED IF I NE J).                                             
C                                                                       
C    SRO DOES NOT REARRANGE THE ORDER OF THE ROWS, BUT DOES MOVE        
C    NONZEROES FROM ONE ROW TO ANOTHER TO ENSURE THAT IF M(I,J) WILL BE 
C    IN THE UPPER TRIANGLE OF M WITH RESPECT TO THE NEW ORDERING, THEN  
C    M(I,J) IS STORED IN ROW I (AND THUS M(J,I) IS NOT STORED);  WHEREAS
C    IF M(I,J) WILL BE IN THE STRICT LOWER TRIANGLE OF M, THEN M(J,I) IS
C    STORED IN ROW J (AND THUS M(I,J) IS NOT STORED).                   
C                                                                     
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   Q            | -->| INTEGER ONE-DIMENSIONAL WORK ARRAY
C |   R            | -->| INTEGER ONE-DIMENSIONAL WORK ARRAY  
C |                |    | DIMENSION = NUMBER OF 
C |                |    | NONZERO ENTRIES IN THE UPPER TRIANGLE OF M
C |   DFLAG        | -->| LOGICAL VARIABLE;  IF DFLAG = .TRUE., THEN 
C |                |    | STORE NONZERO DIAGONAL ELEMENTS AT THE 
C |                |    | BEGINNING OF THE ROW      
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_SD_SRO => SD_SRO
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
      INTEGER, INTENT(IN)             :: N
      INTEGER, INTENT(IN)             :: IP(*)
      INTEGER, INTENT(INOUT)          :: JA(*),R(*),Q(N),IA(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A(*)
      LOGICAL, INTENT(IN)             :: DFLAG                               
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,J,K,JMIN,JMAX,ILAST,JDUMMY,JAK 
      DOUBLE PRECISION AK                                                                                     
C                                                                       
C-----------------------------------------------------------------------                                                 
C                                                                      
C--PHASE 1 -- FIND ROW IN WHICH TO STORE EACH NONZERO                   
C----INITIALIZE COUNT OF NONZEROES TO BE STORED IN EACH ROW
C             
      DO I=1,N                                                      
        Q(I) = 0
      ENDDO                                                      
C                                                                       
C----FOR EACH NONZERO ELEMENT A(J)
C                                      
      DO 3 I=1,N                                                      
        JMIN = IA(I)                                                  
        JMAX = IA(I+1) - 1                                            
        IF(JMIN.GT.JMAX)  GO TO 3                                    
        DO 2 J=JMIN,JMAX                                              
C                                                                       
C--------FIND ROW (=R(J)) AND COLUMN (=JA(J)) IN WHICH TO STORE A(J) 
C
          K = JA(J)                                                   
          IF (IP(K).LT.IP(I))  JA(J) = I                              
          IF (IP(K).GE.IP(I))  K = I                                  
          R(J) = K                                                    
C                                                                       
C--------AND INCREMENT COUNT OF NONZEROES (=Q(R(J)) IN THAT ROW
C     
          Q(K) = Q(K) + 1
2       CONTINUE                                             
3     CONTINUE                                                      
C                                                                      
C--PHASE 2 -- FIND NEW IA AND PERMUTATION TO APPLY TO (JA,A)            
C----DETERMINE POINTERS TO DELIMIT ROWS IN PERMUTED (JA,A)
C              
      DO 4 I=1,N                                                      
        IA(I+1) = IA(I) + Q(I)                                        
        Q(I) = IA(I+1)
4     CONTINUE                                                
C                                                                       
C----DETERMINE WHERE EACH (JA(J),A(J)) IS STORED IN PERMUTED (JA,A)     
C----FOR EACH NONZERO ELEMENT (IN REVERSE ORDER)
C                        
      ILAST = 0                                                       
      JMIN = IA(1)                                                    
      JMAX = IA(N+1) - 1                                              
      J = JMAX                                                        
      DO 6 JDUMMY=JMIN,JMAX                                           
        I = R(J)                                                      
        IF(.NOT.DFLAG .OR. JA(J).NE.I .OR. I.EQ.ILAST)  GO TO 5      
C                                                                       
C------IF DFLAG, THEN PUT DIAGONAL NONZERO AT BEGINNING OF ROW
C          
        R(J) = IA(I)                                                
        ILAST = I                                                   
        GO TO 6                                                     
C                                                                       
C------PUT (OFF-DIAGONAL) NONZERO IN LAST UNUSED LOCATION IN ROW
C        
5       Q(I) = Q(I) - 1                                             
        R(J) = Q(I)                                                 
C                                                                       
        J = J-1 
6     CONTINUE                                                      
C                                                                                                                                              
C--PHASE 3 -- PERMUTE (JA,A) TO UPPER TRIANGULAR FORM (WRT NEW ORDERING)
C
      DO 8 J=JMIN,JMAX                                                
7       IF (R(J).EQ.J)  GO TO 8                                       
        K = R(J)                                                    
        R(J) = R(K)                                                 
        R(K) = K                                                    
        JAK = JA(K)                                                 
        JA(K) = JA(J)                                               
        JA(J) = JAK                                                 
        AK = A(K)                                                   
        A(K) = A(J)                                                 
        A(J) = AK                                                   
        GO TO 7                                                     
8     CONTINUE                                                      
C                                                                       
C-----------------------------------------------------------------------
C                                                                        
      RETURN                                                          
      END
