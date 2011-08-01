C                         ******************
                          SUBROUTINE SD_SDRV
C                         ******************
C
     *(N,P,IP,IA,JA,A,B,Z,NSP,ISP,RSP,ESP,PATH,FLAG)
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
C    SDRV SOLVES SPARSE SYMMETRIC POSITIVE DEFINITE SYSTEMS OF LINEAR   
C    EQUATIONS.  THE SOLUTION PROCESS IS DIVIDED INTO THREE STAGES --   
C                                                                       
C      SSF - THE COEFFICIENT MATRIX M IS FACTORED SYMBOLICALLY TO       
C            DETERMINE WHERE FILLIN WILL OCCUR DURING THE NUMERIC       
C            FACTORIZATION.                                             
C                                                                       
C      SNF - M IS FACTORED NUMERICALLY INTO THE PRODUCT UT-D-U, WHERE   
C            D IS DIAGONAL AND U IS UNIT UPPER TRIANGULAR.              
C                                                                       
C      SNS - THE LINEAR SYSTEM  MX = B  IS SOLVED USING THE UT-D-U      
C            FACTORIZATION FROM SNF.                                    
C                                                                       
C    FOR SEVERAL SYSTEMS WITH THE SAME COEFFICIENT MATRIX, SSF AND SNF  
C    NEED BE DONE ONLY ONCE (FOR THE FIRST SYSTEM);  THEN SNS IS DONE   
C    ONCE FOR EACH ADDITIONAL RIGHT-HAND SIDE.  FOR SEVERAL SYSTEMS     
C    WHOSE COEFFICIENT MATRICES HAVE THE SAME NONZERO STRUCTURE, SSF    
C    NEED BE DONE ONLY ONCE (FOR THE FIRST SYSTEM);  THEN SNF AND SNS   
C    ARE DONE ONCE FOR EACH ADDITIONAL SYSTEM.                          
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
C    IN THE UPPER TRIANGLE NEED BE STORED, FOR EXAMPLE, THE MATRIX      
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
C          A \ 1  2  3  4  5  6  7  8  9          .                     
C                                                                       
C                                                                       
C  REORDERING THE ROWS AND COLUMNS OF M                                 
C                                                                       
C    A SYMMETRIC PERMUTATION OF THE ROWS AND COLUMNS OF THE COEFFICIENT 
C    MATRIX M (E.G., WHICH REDUCES FILLIN OR ENHANCES NUMERICAL         
C    STABILITY) MUST BE SPECIFIED.  THE SOLUTION Z IS RETURNED IN THE   
C    ORIGINAL ORDER.                                                    
C                                                                       
C    TO SPECIFY THE TRIVIAL ORDERING (I.E., THE IDENTITY PERMUTATION),  
C    SET  P(I) = IP(I) = I,  I=1,...,N.  IN THIS CASE, P AND IP CAN BE  
C    THE SAME ARRAY.                                                    
C                                                                       
C    IF A NONTRIVIAL ORDERING (I.E., NOT THE IDENTITY PERMUTATION) IS   
C    SPECIFIED AND M IS STORED SYMMETRICALLY (I.E., NOT BOTH M(I,J) AND 
C    M(J,I) ARE STORED FOR I NE J), THEN ODRV SHOULD BE CALLED (WITH    
C    PATH = 3 OR 5) TO SYMMETRICALLY REORDER (IA,JA,A) BEFORE CALLING   
C    SDRV.  THIS IS TO ENSURE THAT IF M(I,J) WILL BE IN THE UPPER       
C    TRIANGLE OF M WITH RESPECT TO THE NEW ORDERING, THEN M(I,J) IS     
C    STORED IN ROW I (AND THUS M(J,I) IS NOT STORED);  WHEREAS IF M(I,J)
C    WILL BE IN THE STRICT LOWER TRIANGLE OF M, THEN M(J,I) IS STORED IN
C    ROW J (AND THUS M(I,J) IS NOT STORED).                             
C                                                 
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
C | A              |<-- | REAL ONE-DIMENSIONAL ARRAY CONTAINING THE 
C |                |    | NONZERO ENTRIES IN (THE UPPER TRIANGLE OF) M,
C |                |    | STORED BY ROWS;  DIMENSION =NUMBER OF NONZERO
C |                |    | ENTRIES IN (THE UPPER TRIANGLE OF) M 
C | B              |    | REAL ONE-DIMENSIONAL ARRAY CONTAINING THE 
C |                |    | RIGHT-HAND SIDE B; B AND Z CAN BE THE SAME ARRAY;  
C |                |    | DIMENSION = N
C | Z              |    | REAL ONE-DIMENSIONAL ARRAY CONTAINING THE 
C |                |    | SOLUTION X;  Z AND B CAN BE THE SAME ARRAY;  
C |                |    | DIMENSION = N  
C | NSP            | -->| DECLARED DIMENSION OF THE ONE-DIMENSIONAL 
C |                |    | ARRAY ISP;  NSP MUST BE AT LEAST  3N+4K,  
C |                |    | WHERE K IS THE NUMBER OF NONZEROES
C |                |    | IN THE STRICT UPPER TRIANGLE OF M   
C | ISP            |<-- | INTEGER ONE-DIMENSIONAL ARRAY USED FOR 
C |                |    | WORKING STORAGE; DIMENSION = NSP
C | RSP            | -- | REAL ONE-DIMENSIONAL ARRAY USED FOR WORKING 
C |                |    | STORAGE;  RSP AND ISP SHOULD BE EQUIVALENCED; 
C |                |    | DIMENSION = NSP
C | ESP            | -- | INTEGER VARIABLE;  IF SUFFICIENT STORAGE WAS 
C |                |    | AVAILABLE TO PERFORM THE SYMBOLIC 
C |                |    | FACTORIZATION (SSF), THEN ESP IS SET TO
C |                |    | THE AMOUNT OF EXCESS STORAGE PROVIDED 
C |                |    | (NEGATIVE IF INSUFFICIENT STORAGE WAS 
C |                |    | AVAILABLE TO PERFORM THE NUMERIC 
C |                |    | FACTORIZATION (SNF))   
C | PATH           | -->| INTEGER PATH SPECIFICATION;  
C |                |    | VALUES AND THEIR MEANINGS ARE -
C |                |    | 1  PERFORM SSF, SNF, AND SNS                              
C |                |    | 2  PERFORM SNF AND SNS (ISP/RSP IS ASSUMED 
C |                |    |    TO HAVE BEEN SET UP IN AN EARLIER CALL 
C |                |    |    TO SDRV (FOR SSF))         
C |                |    | 3  PERFORM SNS ONLY (ISP/RSP IS ASSUMED 
C |                |    |    TO HAVE BEEN SET UP IN AN EARLIER CALL 
C |                |    |    TO SDRV (FOR SSF AND SNF))     
C |                |    | 4  PERFORM SSF                                            
C |                |    | 5  PERFORM SSF AND SNF                                    
C |                |    | 6  PERFORM SNF ONLY (ISP/RSP IS ASSUMED TO 
C |                |    |    HAVE BEEN SET UP IN AN EARLIER CALL TO 
C |                |    |    SDRV (FOR SSF)) 
C | FLAG           |<-- | INTEGER ERROR FLAG; VALUES AND THEIR MEANINGS:        
C |                |    | 0     NO ERRORS DETECTED                                
C |                |    | 2N+K   DUPLICATE ENTRY IN A  --  ROW = K                 
C |                |    | 6N+K   INSUFFICIENT STORAGE IN SSF -- ROW = K          
C |                |    | 7N+1   INSUFFICIENT STORAGE IN SNF                       
C |                |    | 8N+K   ZERO PIVOT  --  ROW = K                           
C |                |    | 10N+1  INSUFFICIENT STORAGE IN SDRV                      
C |                |    | 11N+1  ILLEGAL PATH SPECIFICATION   
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_SD_SDRV => SD_SDRV
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: N,NSP,PATH
      INTEGER, INTENT(INOUT) :: FLAG,P(*),IP(*),IA(*),JA(*),ISP(*),ESP
      DOUBLE PRECISION, INTENT(IN)    :: B(N)
      DOUBLE PRECISION, INTENT(INOUT) :: A(1),Z(N),RSP(NSP)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
      INTEGER RATIO,Q,MARK,D,U,TMP,UMAX,IJU,IU,IL,JL,JU,JUMAX                        
C                                                                       
C-----------------------------------------------------------------------                                                                            
C                         
      DATA RATIO/2/                                                  
C                                                                       
C----VALIDATE PATH SPECIFICATION 
C                                       
      IF(PATH.LT.1.OR.6.LT.PATH) GO TO 111                        
C                                                                       
C----ALLOCATE STORAGE AND FACTOR M SYMBOLICALLY TO DETERMINE FILL-IN
C    
      IJU   = 1                                                       
      IU    = IJU     +  N                                            
      JL    = IU      +  N+1                                           
      JU    = JL      +  N                                            
      Q     = (NSP+1) -  N                                            
      MARK  = Q       -  N                                            
      JUMAX = MARK    - JU                                            
C                                                                       
      IF((PATH-1) * (PATH-4) * (PATH-5) .NE. 0)  GO TO 1             
      IF(JUMAX.LE.0) GO TO 110                                    
      CALL SD_SSF(N,P,IP,IA,JA,ISP(IJU),ISP(JU),ISP(IU),JUMAX,   
     *            ISP(Q),ISP(MARK),ISP(JL),FLAG)                         
      IF (FLAG.NE.0) GO TO 100                                     
C                                                                       
C----ALLOCATE STORAGE AND FACTOR M NUMERICALLY
C                          
1     IL   = JU      + ISP(IJU+(N-1))                                 
      TMP  = ((IL-1)+(RATIO-1)) / RATIO  +  1                         
      D    = TMP     + N                                              
      U    = D       + N                                              
      UMAX = (NSP+1) - U                                              
      ESP  = UMAX    - (ISP(IU+N)-1)                                  
C                                                                       
      IF ((PATH-1) * (PATH-2) * (PATH-5) * (PATH-6) .NE. 0)  GO TO 2  
      IF (UMAX.LE.0)  GO TO 110                                     
      CALL SD_SNF(N,P,IP,IA,JA,A,                                    
     *            RSP(D),ISP(IJU),ISP(JU),ISP(IU),RSP(U),UMAX,        
     *            ISP(IL),ISP(JL),FLAG)                                 
      IF(FLAG.NE.0)  GO TO 100                                     
C                                                                       
C----SOLVE SYSTEM OF LINEAR EQUATIONS  MX = B
C                           
2     IF((PATH-1) * (PATH-2) * (PATH-3) .NE. 0)  GO TO 3             
      IF (UMAX.LE.0)  GO TO 110                                     
      CALL SD_SNS(N,P,RSP(D),ISP(IJU),ISP(JU),
     *            ISP(IU),RSP(U),Z,B,RSP(TMP))                                                 
C                                                                       
3     RETURN                                                          
C                                                                       
C ** ERROR -- ERROR DETECTED IN SSF, SNF, OR SNS
C                        
100   RETURN 
C                                                         
C ** ERROR -- INSUFFICIENT STORAGE
C                                      
110   FLAG = 10*N + 1                                                 
      RETURN
C                                                          
C ** ERROR -- ILLEGAL PATH SPECIFICATION
C                                
111   FLAG = 11*N + 1
C                                                                       
C-----------------------------------------------------------------------
C                                                  
      RETURN                                                          
      END
