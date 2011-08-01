C                        ******************
                         SUBROUTINE SD_ODRV
C                        ******************
C
     *(N,IA,JA,A,P,IP,NSP,ISP,PATH,FLAG)
C
C***********************************************************************
C BIEF VERSION 5.7     20/11/06   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C
C  FONCTION : DRIVER FOR SPARSE MATRIX REORDERING ROUTINE
C
C    ODRV FINDS A MINIMUM DEGREE ORDERING OF THE ROWS AND COLUMNS OF A  
C    SYMMETRIC MATRIX M STORED IN (IA,JA,A) FORMAT (SEE BELOW).  FOR THE
C    REORDERED MATRIX, THE WORK AND STORAGE REQUIRED TO PERFORM GAUSSIAN
C    ELIMINATION IS (USUALLY) SIGNIFICANTLY LESS.                       
C                                                                       
C    IF ONLY THE NONZERO ENTRIES IN THE UPPER TRIANGLE OF M ARE BEING   
C    STORED, THEN ODRV SYMMETRICALLY REORDERS (IA,JA,A), (OPTIONALLY)   
C    WITH THE DIAGONAL ENTRIES PLACED FIRST IN EACH ROW.  THIS IS TO    
C    ENSURE THAT IF M(I,J) WILL BE IN THE UPPER TRIANGLE OF M WITH      
C    RESPECT TO THE NEW ORDERING, THEN M(I,J) IS STORED IN ROW I (AND   
C    THUS M(J,I) IS NOT STORED);  WHEREAS IF M(I,J) WILL BE IN THE      
C    STRICT LOWER TRIANGLE OF M, THEN M(J,I) IS STORED IN ROW J (AND    
C    THUS M(I,J) IS NOT STORED).                                        
C                                                                       
C                                                                       
C  STORAGE OF SPARSE MATRICES                                           
C                                                                       
C    THE NONZERO ENTRIES OF THE MATRIX M ARE STORED ROW-BY-ROW IN THE   
C    ARRAY A.  TO IDENTIFY THE INDIVIDUAL NONZERO ENTRIES IN EACH ROW,  
C    WE NEED TO KNOW IN WHICH COLUMN EACH ENTRY LIES.  THESE COLUMN     
C    INDICES ARE STORED IN THE ARRAY JA;  I.E., IF  A(K) = M(I,J),  THEN
C    JA(K) = J.  TO IDENTIFY THE INDIVIDUAL ROWS, WE NEED TO KNOW WHERE 
C    EACH ROW STARTS.  THESE ROW POINTERS ARE STORED IN THE ARRAY IA;   
C    I.E., IF M(I,J) IS THE FIRST NONZERO ENTRY (STORED) IN THE I-TH ROW
C    AND  A(K) = M(I,J),  THEN  IA(I) = K.  MOREOVER, IA(N+1) POINTS TO 
C    THE FIRST LOCATION FOLLOWING THE LAST ELEMENT IN THE LAST ROW.     
C    THUS, THE NUMBER OF ENTRIES IN THE I-TH ROW IS  IA(I+1) - IA(I),   
C    THE NONZERO ENTRIES IN THE I-TH ROW ARE STORED CONSECUTIVELY IN    
C                                                                       
C            A(IA(I)),  A(IA(I)+1),  ..., A(IA(I+1)-1),                 
C                                                                       
C    AND THE CORRESPONDING COLUMN INDICES ARE STORED CONSECUTIVELY IN   
C                                                                       
C            JA(IA(I)), JA(IA(I)+1), ..., JA(IA(I+1)-1).                
C                                                                       
C    SINCE THE COEFFICIENT MATRIX IS SYMMETRIC, ONLY THE NONZERO ENTRIES
C    IN THE UPPER TRIANGLE NEED BE STORED.  FOR EXAMPLE, THE MATRIX     
C                                                                       
C             ( 1  0  2  3  0 )                                         
C             ( 0  4  0  0  0 )                                         
C         M = ( 2  0  5  6  0 )                                         
C             ( 3  0  6  7  8 )                                         
C             ( 0  0  0  8  9 )                                         
C                                                                       
C    COULD BE STORED AS                                                 
C                                                                       
C            \ 1  2  3  4  5  6  7  8  9 10 11 12 13                    
C         ---+--------------------------------------                    
C         IA \ 1  4  5  8 12 14                                         
C         JA \ 1  3  4  2  1  3  4  1  3  4  5  4  5                    
C          A \ 1  2  3  4  2  5  6  3  6  7  8  8  9                    
C                                                                       
C    OR (SYMMETRICALLY) AS                                              
C                                                                       
C            \ 1  2  3  4  5  6  7  8  9                                
C         ---+--------------------------                                
C         IA \ 1  4  5  7  9 10                                         
C         JA \ 1  3  4  2  3  4  4  5  5                                
C          A \ 1  2  3  4  5  6  7  8  9                   
C                          
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C | N              | -->| ORDER OF THE MATRIX    
C | IA             |<-- | INTEGER ONE-DIMENSIONAL ARRAY CONTAINING 
C |                |    | POINTERS TO DELIMIT ROWS IN JA AND A;  
C |                |    | DIMENSION = N+1                        
C | JA             |<-- | INTEGER ONE-DIMENSIONAL ARRAY CONTAINING THE 
C |                |    | COLUMN INDICES CORRESPONDING TO THE ELEMENTS 
C |                |    | OF A;  DIMENSION = NUMBER OF NONZERO ENTRIES 
C |                |    | IN (THE UPPER TRIANGLE OF) M 
C | A              |<-- | REAL ONE-DIMENSIONAL ARRAY CONTAINING THE 
C |                |    | NONZERO ENTRIES IN (THE UPPER TRIANGLE OF) M,
C |                |    | STORED BY ROWS;  DIMENSION =NUMBER OF NONZERO
C |                |    | ENTRIES IN (THE UPPER TRIANGLE OF) M 
C | P              |<-- | INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN 
C |                |    | THE PERMUTATION OF THE ROWS AND COLUMNS OF M 
C |                |    | CORRESPONDING TO THE MINIMUM DEGREE ORDERING;  
C |                |    | DIMENSION = N 
C | IP             |<-- | INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN
C |                |    | THE INVERSE OF THE PERMUTATION RETURNED IN P;  
C |                |    | DIMENSION = N 
C | NSP            | -->| DECLARED DIMENSION OF THE ONE-DIMENSIONAL 
C |                |    | ARRAY ISP;  NSP MUST BE AT LEAST  3N+4K,  
C |                |    | WHERE K IS THE NUMBER OF NONZEROES
C |                |    | IN THE STRICT UPPER TRIANGLE OF M   
C | ISP            |<-- | INTEGER ONE-DIMENSIONAL ARRAY USED FOR 
C |                |    | WORKING STORAGE; DIMENSION = NSP 
C | PATH           | -->| INTEGER PATH SPECIFICATION;  
C |                |    | VALUES AND THEIR MEANINGS ARE -
C |                |    | 1  FIND MINIMUM DEGREE ORDERING ONLY                      
C |                |    | 2  FIND MINIMUM DEGREE ORDERING AND 
C |                |    |    REORDER SYMMETRICALLY 
C |                |    |    STORED MATRIX (USED WHEN ONLY THE NONZERO 
C |                |    |    ENTRIES IN THE UPPER TRIANGLE OF M ARE 
C |                |    |    BEING STORED)            
C |                |    | 3  REORDER SYMMETRICALLY STORED MATRIX AS 
C |                |    |    SPECIFIED BY INPUT PERMUTATION (USED WHEN
C |                |    |    AN ORDERING HAS ALREADY BEEN DETERMINED 
C |                |    |    AND ONLY THE NONZERO ENTRIES IN THE  
C |                |    |    UPPER TRIANGLE OF M ARE BEING STORED)                
C |                |    | 4  SAME AS 2 BUT PUT DIAGONAL ENTRIES AT 
C |                |    |    START OF EACH ROW
C |                |    | 5  SAME AS 3 BUT PUT DIAGONAL ENTRIES AT 
C |                |    |    START OF EACH ROW
C | FLAG           |<-- | INTEGER ERROR FLAG;  VALUES AND THEIR MEANINGS ARE -        
C |                |    | 0    NO ERRORS DETECTED                                 
C |                |    | 9N+K  INSUFFICIENT STORAGE IN MD                         
C |                |    | 10N+1  INSUFFICIENT STORAGE IN ODRV                       
C |                |    | 11N+1  ILLEGAL PATH SPECIFICATION
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_SD_ODRV => SD_ODRV
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
      INTEGER, INTENT(IN)             :: N,NSP,PATH
      INTEGER, INTENT(INOUT)          :: FLAG
      INTEGER, INTENT(INOUT)          :: IA(N),JA(*),P(N),IP(N),ISP(NSP)
      DOUBLE PRECISION, INTENT(INOUT) :: A(*)                      
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
      INTEGER V,L,HEAD,TMP,Q,NEXT,MAX                                        
      LOGICAL DFLAG                                                
C                                                                       
C-----------------------------------------------------------------------                                                
C                                                                       
C----INITIALIZE ERROR FLAG AND VALIDATE PATH SPECIFICATION
C              
      FLAG = 0                                                        
      IF(PATH.LT.1.OR.5.LT.PATH) GO TO 111                        
C                                                                       
C----ALLOCATE STORAGE AND FIND MINIMUM DEGREE ORDERING
C                  
      IF((PATH-1)*(PATH-2)*(PATH-4).NE.0) GO TO 1             
      MAX = (NSP-N)/2                                               
      V    = 1                                                      
      L    = V     +  MAX                                           
      HEAD = L     +  MAX                                           
      NEXT = HEAD  +  N                                             
      IF(MAX.LT.N) GO TO 110                                      
C                                                                       
      CALL SD_MD(N,IA,JA,MAX,ISP(V),ISP(L),ISP(HEAD),P,IP,ISP(V),FLAG)
C
      IF(FLAG.NE.0) GO TO 100                                     
C                                                                       
C----ALLOCATE STORAGE AND SYMMETRICALLY REORDER MATRIX
C                  
1     IF ((PATH-2) * (PATH-3) * (PATH-4) * (PATH-5) .NE. 0)  GO TO 2
C  
      TMP = (NSP+1) -      N                                        
      Q   = TMP     - (IA(N+1)-1)                                   
      IF (Q.LT.1)  GO TO 110                                        
C                                                                       
      DFLAG = PATH.EQ.4 .OR. PATH.EQ.5                              
      CALL SD_SRO(N,IP,IA,JA,A,ISP(TMP),ISP(Q),DFLAG)           
C                                                                       
2     RETURN                                                          
C                                                                       
C ** ERROR -- ERROR DETECTED IN MD
C                                      
100   RETURN
C                                                          
C ** ERROR -- INSUFFICIENT STORAGE
C                                      
110   FLAG = 10*N + 1                                                 
      RETURN
C                                                          
C ** ERROR -- ILLEGAL PATH SPECIFIED
C                                    
111   FLAG = 11*N + 1
C                                                                       
C-----------------------------------------------------------------------
C                                                 
      RETURN                                                          
      END
