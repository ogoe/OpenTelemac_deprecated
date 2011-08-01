C                          ******************
                           SUBROUTINE SD_NDRV
C                          ******************
C
     *(N,R,C,IC,IA,JA,A,B,Z,NSP,ISP,RSP,ESP,PATH,FLAG) 
C
C***********************************************************************
C BIEF VERSION 5.9     18/02/08   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C
C  FUNCTION:  DRIVER FOR SUBROUTINES FOR SOLVING SPARSE NONSYMMETRIC
C             SYSTEMS OF LINEAR EQUATIONS (UNCOMPRESSED POINTER STORAGE)    
C
C    PARAMETERS                                                         
C    CLASS ABBREVIATIONS ARE --                                         
C       N - INTEGER VARIABLE                                            
C       F - REAL VARIABLE                                               
C       V - SUPPLIES A VALUE TO THE DRIVER                              
C       R - RETURNS A RESULT FROM THE DRIVER                            
C       I - USED INTERNALLY BY THE DRIVER                               
C       A - ARRAY                                                       
C                                                                       
C CLASS \ PARAMETER                                                     
C ------+----------                                                     
C       \                                                               
C         THE NONZERO ENTRIES OF THE COEFFICIENT MATRIX M ARE STORED    
C    ROW-BY-ROW IN THE ARRAY A.  TO IDENTIFY THE INDIVIDUAL NONZERO     
C    ENTRIES IN EACH ROW, WE NEED TO KNOW IN WHICH COLUMN EACH ENTRY    
C    LIES.  THE COLUMN INDICES WHICH CORRESPOND TO THE NONZERO ENTRIES  
C    OF M ARE STORED IN THE ARRAY JA;  I.E., IF  A(K) = M(I,J),  THEN   
C    JA(K) = J.  IN ADDITION, WE NEED TO KNOW WHERE EACH ROW STARTS AND 
C    HOW LONG IT IS.  THE INDEX POSITIONS IN JA AND A WHERE THE ROWS OF 
C    M BEGIN ARE STORED IN THE ARRAY IA;  I.E., IF M(I,J) IS THE FIRST  
C    NONZERO ENTRY (STORED) IN THE I-TH ROW AND A(K) = M(I,J),  THEN    
C    IA(I) = K.  MOREOVER, THE INDEX IN JA AND A OF THE FIRST LOCATION  
C    FOLLOWING THE LAST ELEMENT IN THE LAST ROW IS STORED IN IA(N+1).   
C    THUS, THE NUMBER OF ENTRIES IN THE I-TH ROW IS GIVEN BY            
C    IA(I+1) - IA(I),  THE NONZERO ENTRIES OF THE I-TH ROW ARE STORED   
C    CONSECUTIVELY IN                                                   
C            A(IA(I)),  A(IA(I)+1),  ..., A(IA(I+1)-1),                 
C    AND THE CORRESPONDING COLUMN INDICES ARE STORED CONSECUTIVELY IN   
C            JA(IA(I)), JA(IA(I)+1), ..., JA(IA(I+1)-1).                
C    FOR EXAMPLE, THE 5 BY 5 MATRIX                                     
C                ( 1. 0. 2. 0. 0.)                                      
C                ( 0. 3. 0. 0. 0.)                                      
C            M = ( 0. 4. 5. 6. 0.)                                      
C                ( 0. 0. 0. 7. 0.)                                      
C                ( 0. 0. 0. 8. 9.)                                      
C    WOULD BE STORED AS                                                 
C               \ 1  2  3  4  5  6  7  8  9                             
C            ---+--------------------------                             
C            IA \ 1  3  4  7  8 10                                      
C            JA \ 1  3  2  2  3  4  4  4  5                             
C             A \ 1. 2. 3. 4. 5. 6. 7. 8. 9.         .                  
C                                                                       
C NV    \ N     - NUMBER OF VARIABLES/EQUATIONS.                        
C FVA   \ A     - NONZERO ENTRIES OF THE COEFFICIENT MATRIX M, STORED   
C       \           BY ROWS.                                            
C       \           SIZE = NUMBER OF NONZERO ENTRIES IN M.              
C NVA   \ IA    - POINTERS TO DELIMIT THE ROWS IN A.                    
C       \           SIZE = N+1.                                         
C NVA   \ JA    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF A.    
C       \           SIZE = SIZE OF A.                                   
C FVA   \ B     - RIGHT-HAND SIDE B;  B AND Z CAN THE SAME ARRAY.       
C       \           SIZE = N.                                           
C FRA   \ Z     - SOLUTION X;  B AND Z CAN BE THE SAME ARRAY.           
C       \           SIZE = N.                                           
C                                                                       
C         THE ROWS AND COLUMNS OF THE ORIGINAL MATRIX M CAN BE          
C    REORDERED (E.G., TO REDUCE FILLIN OR ENSURE NUMERICAL STABILITY)   
C    BEFORE CALLING THE DRIVER.  IF NO REORDERING IS DONE, THEN SET     
C    R(I) = C(I) = IC(I) = I  FOR I=1,...,N.  THE SOLUTION Z IS RETURNED
C    IN THE ORIGINAL ORDER.                                             
C                                                                       
C NVA   \ R     - ORDERING OF THE ROWS OF M.                            
C       \           SIZE = N.                                           
C NVA   \ C     - ORDERING OF THE COLUMNS OF M.                         
C       \           SIZE = N.                                           
C NVA   \ IC    - INVERSE OF THE ORDERING OF THE COLUMNS OF M;  I.E.,   
C       \           IC(C(I)) = I  FOR I=1,...,N.                        
C       \           SIZE = N.                                           
C                                                                       
C         THE SOLUTION OF THE SYSTEM OF LINEAR EQUATIONS IS DIVIDED INTO
C    THREE STAGES --                                                    
C      NSF -- THE MATRIX M IS PROCESSED SYMBOLICALLY TO DETERMINE WHERE 
C              FILLIN WILL OCCUR DURING THE NUMERIC FACTORIZATION.      
C      NNF -- THE MATRIX M IS FACTORED NUMERICALLY INTO THE PRODUCT LDU 
C              OF A UNIT LOWER TRIANGULAR MATRIX L, A DIAGONAL MATRIX D,
C              AND A UNIT UPPER TRIANGULAR MATRIX U, AND THE SYSTEM     
C              MX = B  IS SOLVED.                                       
C      NNS -- THE LINEAR SYSTEM  MX = B  IS SOLVED USING THE LDU        
C  OR          FACTORIZATION FROM NNF.                                  
C      NNT -- THE TRANSPOSED LINEAR SYSTEM  MT X = B  IS SOLVED USING   
C              THE LDU FACTORIZATION FROM NNF.                          
C    FOR SEVERAL SYSTEMS WHOSE COEFFICIENT MATRICES HAVE THE SAME       
C    NONZERO STRUCTURE, NSF NEED BE DONE ONLY ONCE (FOR THE FIRST       
C    SYSTEM);  THEN NNF IS DONE ONCE FOR EACH ADDITIONAL SYSTEM.  FOR   
C    SEVERAL SYSTEMS WITH THE SAME COEFFICIENT MATRIX, NSF AND NNF NEED 
C    BE DONE ONLY ONCE (FOR THE FIRST SYSTEM);  THEN NNS OR NNT IS DONE 
C    ONCE FOR EACH ADDITIONAL RIGHT-HAND SIDE.                          
C                                                                       
C NV    \ PATH  - PATH SPECIFICATION;  VALUES AND THEIR MEANINGS ARE -- 
C       \           1  PERFORM NSF AND NNF.                             
C       \           2  PERFORM NNF ONLY  (NSF IS ASSUMED TO HAVE BEEN   
C       \               DONE IN A MANNER COMPATIBLE WITH THE STORAGE    
C       \               ALLOCATION USED IN THE DRIVER).                 
C       \           3  PERFORM NNS ONLY  (NSF AND NNF ARE ASSUMED TO    
C       \               HAVE BEEN DONE IN A MANNER COMPATIBLE WITH THE  
C       \               STORAGE ALLOCATION USED IN THE DRIVER).         
C       \           4  PERFORM NNT ONLY  (NSF AND NNF ARE ASSUMED TO    
C       \               HAVE BEEN DONE IN A MANNER COMPATIBLE WITH THE  
C       \               STORAGE ALLOCATION USED IN THE DRIVER).         
C       \           5  PERFORM NSF ONLY.                                
C                                                                       
C         VARIOUS ERRORS ARE DETECTED BY THE DRIVER AND THE INDIVIDUAL  
C    SUBROUTINES.                                                       
C                                                                       
C NR    \ FLAG  - ERROR FLAG;  VALUES AND THEIR MEANINGS ARE --         
C       \             0     NO ERRORS DETECTED                          
C       \             N+K   NULL ROW IN A  --  ROW = K                  
C       \            2N+K   DUPLICATE ENTRY IN A  --  ROW = K           
C       \            3N+K   INSUFFICIENT STORAGE IN NSF  --  ROW = K    
C       \            4N+1   INSUFFICIENT STORAGE IN NNF                 
C       \            5N+K   NULL PIVOT  --  ROW = K                     
C       \            6N+K   INSUFFICIENT STORAGE IN NSF  --  ROW = K    
C       \            7N+1   INSUFFICIENT STORAGE IN NNF                 
C       \            8N+K   ZERO PIVOT  --  ROW = K                     
C       \           10N+1   INSUFFICIENT STORAGE IN NDRV                
C       \           11N+1   ILLEGAL PATH SPECIFICATION                  
C                                                                       
C         WORKING STORAGE IS NEEDED FOR THE FACTORED FORM OF THE MATRIX 
C    M PLUS VARIOUS TEMPORARY VECTORS.  THE ARRAYS ISP AND RSP SHOULD BE
C    EQUIVALENCED;  INTEGER STORAGE IS ALLOCATED FROM THE BEGINNING OF  
C    ISP AND REAL STORAGE FROM THE END OF RSP.                          
C                                                                       
C NV    \ NSP   - DECLARED DIMENSION OF RSP;  NSP GENERALLY MUST        
C       \           BE LARGER THAN  5N+3 + 2K  (WHERE  K = (NUMBER OF   
C       \           NONZERO ENTRIES IN M)).                             
C NVIRA \ ISP   - INTEGER WORKING STORAGE DIVIDED UP INTO VARIOUS ARRAYS
C       \           NEEDED BY THE SUBROUTINES;  ISP AND RSP SHOULD BE   
C       \           EQUIVALENCED.                                       
C       \           SIZE = LRATIO*NSP                                   
C FVIRA \ RSP   - REAL WORKING STORAGE DIVIDED UP INTO VARIOUS ARRAYS   
C       \           NEEDED BY THE SUBROUTINES;  ISP AND RSP SHOULD BE   
C       \           EQUIVALENCED.                                       
C       \           SIZE = NSP.                                         
C NR    \ ESP   - IF SUFFICIENT STORAGE WAS AVAILABLE TO PERFORM THE    
C       \           SYMBOLIC FACTORIZATION (NSF), THEN ESP IS SET TO THE
C       \           AMOUNT OF EXCESS STORAGE PROVIDED (NEGATIVE IF      
C       \           INSUFFICIENT STORAGE WAS AVAILABLE TO PERFORM THE   
C       \           NUMERIC FACTORIZATION (NNF)).                       
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
C    CONVERSION TO DOUBLE PRECISION                                       
C                                                                       
C    TO CONVERT THESE ROUTINES FOR DOUBLE PRECISION ARRAYS, SIMPLY USE  
C    THE DOUBLE PRECISION DECLARATIONS IN PLACE OF THE REAL DECLARATIONS
C    IN EACH SUBPROGRAM;  IN ADDITION, THE DATA VALUE OF THE INTEGER    
C    VARIABLE LRATIO MUST BE SET AS INDICATED IN SUBROUTINE NDRV        
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
      INTEGER N,NSP                                                                       
      INTEGER R(*),C(*),IC(*),IA(*),JA(*),ISP(*)        
      INTEGER ESP,PATH,FLAG,Q,IM,D,U,ROW,TMP,UMAX,VMAX                                                  
      DOUBLE PRECISION A(*),B(*),Z(*),RSP(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
      INTEGER IL,JL,IU,JU,JLMAX,JUTMP,JUMAX,J,L,LMAX                  
C                                                                       
C  SET LRATIO EQUAL TO THE RATIO BETWEEN THE LENGTH OF FLOATING POINT   
C  AND INTEGER ARRAY DATA;  E. G., LRATIO = 1 FOR (REAL, INTEGER),      
C  LRATIO = 2 FOR (DOUBLE PRECISION, INTEGER)                           
C
      INTEGER LRATIO                                                                       
      DATA LRATIO/2/                                                  
C
C-----------------------------------------------------------------------
C                                                                       
      IF(PATH.LT.1.OR.5.LT.PATH) GO TO 111                        
C  ******  INITIALIZE AND DIVIDE UP TEMPORARY STORAGE  *****************
      IL = 1                                                          
      IU = IL + N+1                                                   
      JL = IU + N+1                                                   
C                                                                       
C  ******  CALL NSF IF FLAG IS SET  ************************************
      IF ((PATH-1) * (PATH-5) .NE. 0)  GO TO 2                        
      VMAX = (LRATIO*NSP + 1 - JL) - (N+1) - N                       
      JLMAX = VMAX/2                                                 
      Q = JL + JLMAX                                           
      IM = Q + (N+1)                                           
      JUTMP = IM  +   N                                             
      JUMAX = LRATIO*NSP + 1 - JUTMP                                
      ESP = VMAX/LRATIO                                              
      IF (JLMAX.LE.0 .OR. JUMAX.LE.0)  GO TO 110                    
      CALL SD_NSF(N,R,IC,IA,JA,                                       
     *        ISP(IL),ISP(JL),JLMAX,ISP(IU),ISP(JUTMP),JUMAX,     
     *        ISP(Q),ISP(IM),FLAG)                                  
      IF (FLAG.NE.0)  GO TO 100                                     
C  ******  MOVE JU NEXT TO JL  *****************************************
      JLMAX = ISP(IL+N)-1                                           
      JU    = JL + JLMAX                                            
      JUMAX = ISP(IU+N)-1                                           
      IF(JUMAX.GE.1) THEN                                     
        DO J=1,JUMAX                                                
          ISP(JU+J-1) = ISP(JUTMP+J-1)
        ENDDO 
      ENDIF                               
C                                                                       
C  ******  CALL REMAINING SUBROUTINES  *********************************
2     CONTINUE
      JLMAX = ISP(IL+N)-1                                             
      JU    = JL  + JLMAX                                             
      JUMAX = ISP(IU+N)-1                                             
      L     = (JU + JUMAX - 2 + LRATIO)  /  LRATIO    +    1          
      LMAX  = JLMAX                                                   
      D     = L   + LMAX                                              
      U     = D   + N                                                 
      ROW   = NSP + 1 - N                                             
      TMP   = ROW - N                                                 
      UMAX  = TMP - U                                                 
      ESP = UMAX - JUMAX                                              
C                                                                       
      IF((PATH-1)*(PATH-2).NE.0) GO TO 3                        
        IF(UMAX.LE.0) GO TO 110                                     
          CALL SD_NNF(N,R,C,IC,IA,JA,A,Z,B,                         
     *                ISP(IL),ISP(JL),RSP(L),LMAX,RSP(D),                 
     *                ISP(IU),ISP(JU),RSP(U),UMAX,                        
     *                RSP(ROW),RSP(TMP),FLAG)                               
          IF(FLAG.NE.0) GO TO 100                                     
          RETURN                                                        
C                                                                       
   3    IF ((PATH-3) .NE. 0)  GO TO 4                                   
          CALL SD_NNS                                                     
     *       (N,  R, C,                                                 
     *        ISP(IL), ISP(JL), RSP(L),  RSP(D),                        
     *          ISP(IU), ISP(JU), RSP(U),                               
     *        Z,  B,  RSP(TMP))                                         
C                                                                       
4     IF((PATH-4).NE.0) GO TO 5                                   
          CALL SD_NNT(N,R,C,ISP(IL),ISP(JL),RSP(L),RSP(D),                        
     *                ISP(IU),ISP(JU),RSP(U),Z,B,RSP(TMP))                                         
5     RETURN                                                          
C                                                                       
C ** ERROR:  ERROR DETECTED IN NSF, NNF, NNS, OR NNT                    
100   RETURN                                                          
C ** ERROR:  INSUFFICIENT STORAGE                                       
110   FLAG = 10*N + 1                                                 
      RETURN                                                          
C ** ERROR:  ILLEGAL PATH SPECIFICATION                                 
111   FLAG = 11*N + 1                                                 
      RETURN                                                          
      END
