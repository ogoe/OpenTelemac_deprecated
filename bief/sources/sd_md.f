C                        ****************
                         SUBROUTINE SD_MD
C                        ****************
C
     *(N,IA,JA,MAX,V,L,HEAD,LAST,NEXT,MARK,FLAG)
C
C***********************************************************************
C BIEF VERSION 5.7     20/11/06   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C
C  FONCTION : MINIMUM DEGREE ALGORITHM (BASED ON ELEMENT MODEL) 
C
C    MD FINDS A MINIMUM DEGREE ORDERING OF THE ROWS AND COLUMNS OF A    
C    SYMMETRIC MATRIX M STORED IN (IA,JA,A) FORMAT.      
C
C                                                                       
C  ADDITIONAL PARAMETERS                                                
C                                                                       
C    MAX  - DECLARED DIMENSION OF THE ONE-DIMENSIONAL ARRAYS V AND L;   
C           MAX MUST BE AT LEAST  N+2K,  WHERE K IS THE NUMBER OF       
C           NONZEROES IN THE STRICT UPPER TRIANGLE OF M                 
C                                                                       
C    V    - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = MAX        
C                                                                       
C    L    - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = MAX        
C                                                                       
C    HEAD - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N          
C                                                                       
C    LAST - INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN THE PERMUTATION
C           OF THE ROWS AND COLUMNS OF M CORRESPONDING TO THE MINIMUM   
C           DEGREE ORDERING;  DIMENSION = N                             
C                                                                       
C    NEXT - INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN THE INVERSE OF 
C           THE PERMUTATION RETURNED IN LAST;  DIMENSION = N            
C                                                                       
C    MARK - INTEGER ONE-DIMENSIONAL WORK ARRAY (MAY BE THE SAME AS V);  
C           DIMENSION = N                                               
C                                                                       
C    FLAG - INTEGER ERROR FLAG;  VALUES AND THEIR MEANINGS ARE -        
C             0      NO ERRORS DETECTED                                 
C             11N+1  INSUFFICIENT STORAGE IN MD                         
C                                                                       
C                                                                       
C  DEFINITIONS OF INTERNAL PARAMETERS                                   
C                                                                       
C    ---------+---------------------------------------------------------
C    V(S)     \ VALUE FIELD OF LIST ENTRY                               
C    ---------+---------------------------------------------------------
C    L(S)     \ LINK FIELD OF LIST ENTRY  (0 => END OF LIST)            
C    ---------+---------------------------------------------------------
C    L(VI)    \ POINTER TO ELEMENT LIST OF UNELIMINATED VERTEX VI       
C    ---------+---------------------------------------------------------
C    L(EJ)    \ POINTER TO BOUNDARY LIST OF ACTIVE ELEMENT EJ           
C    ---------+---------------------------------------------------------
C    HEAD(D)  \ VJ => VJ HEAD OF D-LIST D                               
C             \  0 => NO VERTEX IN D-LIST D                                                                                                       
C             \          VI UNELIMINATED VERTEX                 
C             \          VI IN EK           \       VI NOT IN EK        
C    ---------+-----------------------------+---------------------------
C    NEXT(VI) \ UNDEFINED BUT NONNEGATIVE   \ VJ => VJ NEXT IN D-LIST   
C             \                             \  0 => VI TAIL OF D-LIST   
C    ---------+-----------------------------+---------------------------
C    LAST(VI) \ (NOT SET UNTIL MDP)         \ -D => VI HEAD OF D-LIST D 
C             \-VK => COMPUTE DEGREE        \ VJ => VJ LAST IN D-LIST   
C             \ EJ => VI PROTOTYPE OF EJ    \  0 => VI NOT IN ANY D-LIST
C             \  0 => DO NOT COMPUTE DEGREE \                           
C    ---------+-----------------------------+---------------------------
C    MARK(VI) \ MARK(VK)                    \ NONNEGATIVE TAG < MARK(VK)
C                                                                       
C                                                                       
C             \                   VI ELIMINATED VERTEX                  
C             \      EI ACTIVE ELEMENT      \           OTHERWISE       
C    ---------+-----------------------------+---------------------------
C    NEXT(VI) \ -J => VI WAS J-TH VERTEX    \ -J => VI WAS J-TH VERTEX  
C             \       TO BE ELIMINATED      \       TO BE ELIMINATED    
C    ---------+-----------------------------+---------------------------
C    LAST(VI) \  M => SIZE OF EI = M        \ UNDEFINED                 
C    ---------+-----------------------------+---------------------------
C    MARK(VI) \ -M => OVERLAP COUNT OF EI   \ UNDEFINED                 
C             \       WITH EK = M           \                           
C             \ OTHERWISE NONNEGATIVE TAG   \                           
C             \       < MARK(VK)            \                       
C                          
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |                | -->|    
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_SD_MD => SD_MD
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: N,MAX
      INTEGER, INTENT(INOUT) :: IA(*),JA(*),V(MAX),L(MAX),HEAD(N)
      INTEGER, INTENT(INOUT) :: LAST(N),NEXT(N),MARK(N),FLAG                                                                   
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
      INTEGER TAG,DMIN,VK,EK,TAIL,K                      
      EQUIVALENCE  (VK,EK)                  
C                                                                       
C-----------------------------------------------------------------------                                         
C                                                                       
C----INITIALIZATION
C                                                     
      TAG = 0                                                         
      CALL SD_MDI(N,IA,JA,MAX,V,L,HEAD,LAST,NEXT,MARK,TAG,FLAG)          
      IF(FLAG.NE.0)  RETURN                                          
C                                                                       
      K = 0                                                           
      DMIN = 1
C                                                                       
C----WHILE  K < N  DO
C                                                   
1     IF(K.GE.N) GO TO 4                                            
C                                                                       
C------SEARCH FOR VERTEX OF MINIMUM DEGREE
C                              
2     IF(HEAD(DMIN).GT.0)  GO TO 3                                 
      DMIN = DMIN + 1                                             
      GO TO 2                                                     
C                                                                       
C------REMOVE VERTEX VK OF MINIMUM DEGREE FROM DEGREE LIST
C              
3     VK = HEAD(DMIN)                                               
      HEAD(DMIN) = NEXT(VK)                                         
      IF (HEAD(DMIN).GT.0) LAST(HEAD(DMIN)) = -DMIN                
C                                                                       
C------NUMBER VERTEX VK, ADJUST TAG, AND TAG VK 
C                        
      K = K+1                                                       
      NEXT(VK) = -K                                                 
      LAST(EK) = DMIN - 1                                           
      TAG = TAG + LAST(EK)                                          
      MARK(VK) = TAG                                                
C                                                                       
C------FORM ELEMENT EK FROM UNELIMINATED NEIGHBORS OF VK
C                
      CALL SD_MDM(VK,TAIL,V,L,LAST,NEXT,MARK)                            
C                                                                       
C------PURGE INACTIVE ELEMENTS AND DO MASS ELIMINATION
C                  
      CALL SD_MDP(K,EK,TAIL,V,L,HEAD,LAST,NEXT,MARK)                     
C                                                                       
C------UPDATE DEGREES OF UNELIMINATED VERTICES IN EK
C                    
      CALL SD_MDU(EK,DMIN,V,L,HEAD,LAST,NEXT,MARK)                       
C                                                                       
      GO TO 1                                                       
C                                                                       
C----GENERATE INVERSE PERMUTATION FROM PERMUTATION
C                      
4     DO 5 K=1,N                                                      
      NEXT(K) = -NEXT(K)                                            
5     LAST(NEXT(K)) = K                                             
C                                                                       
C-----------------------------------------------------------------------
C                                                                        
      RETURN                                                          
      END
