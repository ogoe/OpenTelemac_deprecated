C                          *****************
                           SUBROUTINE SD_SSF
C                          *****************
C
     *(N,P,IP,IA,JA,IJU,JU,IU,JUMAX,Q,MARK,JL,FLAG)         
C
C***********************************************************************
C BIEF VERSION 5.7     20/11/06   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C
C  FONCTION : SYMBOLIC UT-D-U FACTORIZATION OF SPARSE SYMMETRIC MATRIX  
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
C | P              |<-- | INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN 
C |                |    | THE PERMUTATION OF THE ROWS AND COLUMNS OF M 
C |                |    | CORRESPONDING TO THE MINIMUM DEGREE ORDERING;  
C |                |    | DIMENSION = N 
C | IP             |<-- | INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN
C |                |    | THE INVERSE OF THE PERMUTATION RETURNED IN P;  
C |                |    | DIMENSION = N 
C | IA             |<-- | INTEGER ONE-DIMENSIONAL ARRAY CONTAINING 
C |                |    | POINTERS TO DELIMIT ROWS IN JA AND A;  
C |                |    | DIMENSION = N+1                        
C | JA             |<-- | INTEGER ONE-DIMENSIONAL ARRAY CONTAINING THE 
C |                |    | COLUMN INDICES CORRESPONDING TO THE ELEMENTS 
C |                |    | OF A;  DIMENSION = NUMBER OF NONZERO ENTRIES 
C |                |    | IN (THE UPPER TRIANGLE OF) M
C | IJU            |    | INTEGER ONE-DIMENSIONAL ARRAY CONTAINING 
C |                |    | POINTERS TO THE START OF EACH ROW IN JU;  DIMENSION = N 
C | JU             |    | INTEGER ONE-DIMENSIONAL ARRAY CONTAINING THE 
C |                |    | COLUMN INDICES CORRESPONDING TO THE ELEMENTS 
C |                |    | OF U;  DIMENSION = JUMAX
C | IU             |    | INTEGER ONE-DIMENSIONAL ARRAY CONTAINING 
C |                |    | POINTERS TO DELIMIT ROWS IN U;  DIMENSION = N+1
C | JUMAX          |    | DECLARED DIMENSION OF THE ONE-DIMENSIONAL 
C |                |    | ARRAY JU;  JUMAX MUST BE AT LEAST THE SIZE 
C |                |    | OF U MINUS COMPRESSION (IJU(N) AFTER THE 
C |                |    | CALL TO SSF)       
C | Q              | -->| INTEGER ONE-DIMENSIONAL WORK ARRAY, DIMENSION N
C |                |    | Q CONTAINS AN ORDERED LINKED LIST 
C |                |    | REPRESENTATION OF THE NONZERO    
C |                |    | STRUCTURE OF THE K-TH ROW OF U --                                
C |                |    | Q(K) IS THE FIRST COLUMN WITH A NONZERO ENTRY                  
C |                |    | Q(I) IS THE NEXT COLUMN WITH A NONZERO ENTRY 
C |                |    | AFTER COLUMN I    
C |                |    | IN EITHER CASE, Q(I) = N+1 INDICATES THE 
C |                |    | END OF THE LIST 
C | MARK           | -- | INTEGER ONE-DIMENSIONAL WORK ARRAY; DIMENSION N
C |                |    | THE LAST ROW STORED IN JU FOR WHICH U(MARK(I),I)
C |                |    | NE 0
C | JL             | -->| INTEGER ONE-DIMENSIONAL WORK ARRAY; DIMENSION N   
C |                |    | DIMENSION = NUMBER OF NONZERO ENTRIES IN THE 
C |                |    | UPPER TRIANGLE OF M. JL CONTAINS LISTS OF ROWS 
C |                |    | TO BE MERGED INTO UNELIMINATED ROWS --   
C |                |    | I GE K => JL(I) IS THE FIRST ROW TO BE 
C |                |    | MERGED INTO ROW I       
C |                |    | I LT K => JL(I) IS THE ROW FOLLOWING ROW I IN 
C |                |    | SOME LIST OF ROWS
C |                |    | IN EITHER CASE, JL(I) = 0 INDICATES THE 
C |                |    | END OF A LIST 
C | FLAG           |<-- | LOGICAL VARIABLE;  IF DFLAG = .TRUE., THEN 
C |                |    | STORE NONZERO DIAGONAL ELEMENTS AT THE 
C |                |    | BEGINNING OF THE ROW      
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_SD_SSF => SD_SSF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
      INTEGER, INTENT(IN) :: N,JUMAX
      INTEGER, INTENT(INOUT) :: P(N),IP(N),IA(N+1),JA(*),IJU(N),FLAG
      INTEGER, INTENT(INOUT) :: JU(JUMAX),IU(N+1),Q(N),MARK(N),JL(N)                                                                                                                             
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
      INTEGER I,J,M,TAG,VJ,QM,JUMIN,JUPTR,K,LUK,LUI,JMIN,JMAX,LMAX                  
      LOGICAL CLIQUE            
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C----INITIALIZATION
C                                                                       
C    JUMIN AND JUPTR ARE THE INDICES IN JU OF THE FIRST AND LAST        
C    ELEMENTS IN THE LAST ROW SAVED IN JU                             
C                                                                       
C    LUK IS THE NUMBER OF NONZERO ENTRIES IN THE K-TH ROW  
C                                                     
      JUMIN = 1                                                       
      JUPTR = 0                                                       
      IU(1) = 1                                                       
      DO K=1,N                                                      
        MARK(K) = 0                                                   
        JL(K) = 0 
      ENDDO                                                    
C                                                                       
C----FOR EACH ROW K
C                                                     
      DO 18 K=1,N                                                     
        LUK = 0                                                       
        Q(K) = N+1                                                    
C                                                                       
        TAG = MARK(K)                                                 
        CLIQUE = .FALSE.                                              
        IF(JL(K).NE.0)  CLIQUE = JL(JL(K)).EQ.0                      
C                                                                       
C------INITIALIZE NONZERO STRUCTURE OF K-TH ROW TO ROW P(K) OF M
C        
        JMIN = IA(P(K))                                               
        JMAX = IA(P(K)+1) - 1                                         
        IF (JMIN.GT.JMAX)  GO TO 4                                    
        DO 3 J=JMIN,JMAX                                              
          VJ = IP(JA(J))                                              
          IF(VJ.LE.K)  GO TO 3                                       
C                                                                       
          QM = K                                                    
2         M = QM                                                    
          QM = Q(M)                                                 
          IF (QM.LT.VJ)  GO TO 2                                    
          IF (QM.EQ.VJ)  GO TO 102                                  
          LUK = LUK+1                                             
          Q(M) = VJ                                               
          Q(VJ) = QM                                              
          IF (MARK(VJ).NE.TAG)  CLIQUE = .FALSE.                  
C                                                                       
3       CONTINUE                                                    
C                                                                       
C------IF EXACTLY ONE ROW IS TO BE MERGED INTO THE K-TH ROW AND THERE IS
C------A NONZERO ENTRY IN EVERY COLUMN IN THAT ROW IN WHICH THERE IS A  
C------NONZERO ENTRY IN ROW P(K) OF M, THEN DO NOT COMPUTE FILL-IN, JUST
C------USE THE COLUMN INDICES FOR THE ROW WHICH WAS TO HAVE BEEN MERGED
C 
4       IF(.NOT.CLIQUE)  GO TO 5                                     
        IJU(K) = IJU(JL(K)) + 1                                     
        LUK = IU(JL(K)+1) - (IU(JL(K))+1)                           
        GO TO 17                                                    
C                                                                       
C------MODIFY NONZERO STRUCTURE OF K-TH ROW BY COMPUTING FILL-IN        
C------FOR EACH ROW I TO BE MERGED IN
C                                   
5       LMAX = 0                                                      
        IJU(K) = JUPTR                                                
C                                                                       
        I = K                                                         
6       I = JL(I)                                                     
        IF (I.EQ.0)  GO TO 10                                         
C                                                                       
C--------MERGE ROW I INTO K-TH ROW
C                                      
        LUI = IU(I+1) - (IU(I)+1)                                   
        JMIN = IJU(I) +  1                                          
        JMAX = IJU(I) + LUI                                         
        QM = K                                                      
C                                                                       
        DO 8 J=JMIN,JMAX                                            
          VJ = JU(J)                                                
7         M = QM                                                    
          QM = Q(M)                                                 
          IF (QM.LT.VJ)  GO TO 7                                    
          IF (QM.EQ.VJ)  GO TO 8                                    
          LUK = LUK+1                                             
          Q(M) = VJ                                               
          Q(VJ) = QM                                              
          QM = VJ                                                 
8       CONTINUE                                                  
C                                                                       
C--------REMEMBER LENGTH AND POSITION IN JU OF LONGEST ROW MERGED
C       
        IF(LUI.LE.LMAX)  GO TO 6                                   
        LMAX = LUI                                                
        IJU(K) = JMIN                                             
C                                                                       
        GO TO 6                                                     
C                                                                       
C------IF THE K-TH ROW IS THE SAME LENGTH AS THE LONGEST ROW MERGED,    
C------THEN USE THE COLUMN INDICES FOR THAT ROW
C                         
10      IF (LUK.EQ.LMAX)  GO TO 17                                    
C                                                                       
C------IF THE TAIL OF THE LAST ROW SAVED IN JU IS THE SAME AS THE HEAD  
C------OF THE K-TH ROW, THEN OVERLAP THE TWO SETS OF COLUMN INDICES --  
C--------SEARCH LAST ROW SAVED FOR FIRST NONZERO ENTRY IN K-TH ROW ...
C  
        I = Q(K)                                                    
        IF (JUMIN.GT.JUPTR)  GO TO 12                               
        DO 11 JMIN=JUMIN,JUPTR                                      
          IF (JU(JMIN)-I)  11, 13, 12                               
11      CONTINUE                                                  
12      GO TO 15                                                    
C                                                                       
C--------... AND THEN TEST WHETHER TAIL MATCHES HEAD OF K-TH ROW
C        
13      IJU(K) = JMIN                                               
        DO J=JMIN,JUPTR                                          
          IF (JU(J).NE.I)  GO TO 15                                 
          I = Q(I)                                                  
          IF (I.GT.N)  GO TO 17                                     
        ENDDO                                                  
        JUPTR = JMIN - 1                                            
C                                                                       
C------SAVE NONZERO STRUCTURE OF K-TH ROW IN JU
C                         
15      I = K                                                         
        JUMIN = JUPTR +  1                                            
        JUPTR = JUPTR + LUK                                           
        IF (JUPTR.GT.JUMAX)  GO TO 106                                
        DO J=JUMIN,JUPTR                                           
          I = Q(I)                                                    
          JU(J) = I                                                   
          MARK(I) = K
        ENDDO                                                 
        IJU(K) = JUMIN                                                
C                                                                       
C------ADD K TO ROW LIST FOR FIRST NONZERO ELEMENT IN K-TH ROW
C          
17      IF(LUK.LE.1)  GO TO 18                                       
        I = JU(IJU(K))                                              
        JL(K) = JL(I)                                               
        JL(I) = K                                                   
C                                                                       
18    IU(K+1) = IU(K) + LUK                                         
C                                                                       
      FLAG = 0                                                  
      RETURN 
C                                                                                                                             
C ** ERROR -- DUPLICATE ENTRY IN A 
C                                     
102   FLAG = 2*N + P(K)                                               
      RETURN
C                                                          
C ** ERROR -- INSUFFICIENT STORAGE FOR JU
C                               
106   FLAG = 6*N + K                                                  
      RETURN                                                          
      END
