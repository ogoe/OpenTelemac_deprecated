C		      *********************
                      SUBROUTINE SD_SOLVE_1
C		      *********************
C
     &(NPOIN,NSEGB,GLOSEG,MAXSEG,DA,XA,XINC,RHS,INFOGR,TYPEXT)
C
C***********************************************************************
C BIEF VERSION 5.7     20/11/06   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C
C  FONCTION: RESOLUTION DIRECTE D'UN SYSTEME LINEAIRE SYMETRIQUE PAR
C            PERMUTATION MINIMUM DEGRE ET DECOMPOSITION LDLT 
C
C            TRANSFORMATION STOCKAGE SEGMENT EN STOCKAGE COMPACT (MORSE)
C
C-----------------------------------------------------------------------
C                                                                       
C                                                                      
C               YALE SPARSE MATRIX PACKAGE - NONSYMMETRIC CODES         
C                    SOLVING THE SYSTEM OF EQUATIONS MX = B             
C                        (UNCOMPRESSED POINTER STORAGE)                 
C                                                                       
C    I.   CALLING SEQUENCES                                             
C         THE COEFFICIENT MATRIX CAN BE PROCESSED BY AN ORDERING ROUTINE
C    (E.G., TO REDUCE FILLIN OR ENSURE NUMERICAL STABILITY) BEFORE USING
C    THE REMAINING SUBROUTINES.  IF NO REORDERING IS DONE, THEN SET     
C    R(I) = C(I) = IC(I) = I  FOR I=1,...,N.  THE CALLING SEQUENCE IS --
C        (      (MATRIX ORDERING))                                      
C         NSF   (SYMBOLIC FACTORIZATION TO DETERMINE WHERE FILLIN WILL  
C                 OCCUR DURING NUMERIC FACTORIZATION)                   
C         NNF   (NUMERIC FACTORIZATION INTO PRODUCT LDU OF UNIT LOWER   
C                 TRIANGULAR MATRIX L, DIAGONAL MATRIX D, AND UNIT UPPER
C                 TRIANGULAR MATRIX U, AND SOLUTION OF LINEAR SYSTEM)   
C         NNS   (SOLUTION OF LINEAR SYSTEM FOR ADDITIONAL RIGHT-HAND    
C     OR          SIDE USING LDU FACTORIZATION FROM NNF)                
C         NNT   (SOLUTION OF TRANSPOSED LINEAR SYSTEM FOR ADDITIONAL    
C                 RIGHT-HAND SIDE USING LDU FACTORIZATION FROM NNF)     
C                                                                       
C    II.  STORAGE OF SPARSE MATRICES                                    
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
C         THE STRICT TRIANGULAR PORTIONS OF THE MATRICES L AND U ARE    
C    STORED IN THE SAME FASHION USING THE ARRAYS  IL, JL, L  AND        
C    IU, JU, U  RESPECTIVELY.  THE DIAGONAL ENTRIES OF L AND U ARE      
C    ASSUMED TO BE EQUAL TO ONE AND ARE NOT STORED.  THE ARRAY D        
C    CONTAINS THE RECIPROCALS OF THE DIAGONAL ENTRIES OF THE MATRIX D.  
C                                                                       
C    III. ADDITIONAL STORAGE SAVINGS                                    
C         IN NSF, R AND IC CAN BE THE SAME ARRAY IN THE CALLING         
C    SEQUENCE IF NO REORDERING OF THE COEFFICIENT MATRIX HAS BEEN DONE. 
C         IN NNF, R, C AND IC CAN ALL BE THE SAME ARRAY IF NO REORDERING
C    HAS BEEN DONE.  IF ONLY THE ROWS HAVE BEEN REORDERED, THEN C AND IC
C    CAN BE THE SAME ARRAY.  IF THE ROW AND COLUMN ORDERINGS ARE THE    
C    SAME, THEN R AND C CAN BE THE SAME ARRAY.  Z AND ROW CAN BE THE    
C    SAME ARRAY.                                                        
C         IN NNS OR NNT, R AND C CAN BE THE SAME ARRAY IF NO REORDERING 
C    HAS BEEN DONE OR IF THE ROW AND COLUMN ORDERINGS ARE THE SAME.  Z  
C    AND B CAN BE THE SAME ARRAY;  HOWEVER, THEN B WILL BE DESTROYED.   
C                                                                       
C    IV.  PARAMETERS                                                    
C         FOLLOWING IS A LIST OF PARAMETERS TO THE PROGRAMS.  NAMES ARE 
C    UNIFORM AMONG THE VARIOUS SUBROUTINES.  CLASS ABBREVIATIONS ARE -- 
C       N - INTEGER VARIABLE                                            
C       F - REAL VARIABLE                                               
C       V - SUPPLIES A VALUE TO A SUBROUTINE                            
C       R - RETURNS A RESULT FROM A SUBROUTINE                          
C       I - USED INTERNALLY BY A SUBROUTINE                             
C       A - ARRAY                                                       
C                                                                       
C CLASS \ PARAMETER                                                     
C ------+----------                                                     
C FVA   \ A     - NONZERO ENTRIES OF THE COEFFICIENT MATRIX M, STORED   
C       \           BY ROWS.                                            
C       \           SIZE = NUMBER OF NONZERO ENTRIES IN M.              
C FVA   \ B     - RIGHT-HAND SIDE B.                                    
C       \           SIZE = N.                                           
C NVA   \ C     - ORDERING OF THE COLUMNS OF M.                         
C       \           SIZE = N.                                           
C FVRA  \ D     - RECIPROCALS OF THE DIAGONAL ENTRIES OF THE MATRIX D.  
C       \           SIZE = N.                                           
C NR    \ FLAG  - ERROR FLAG;  VALUES AND THEIR MEANINGS ARE --         
C       \            0     NO ERRORS DETECTED                           
C       \            N+K   NULL ROW IN A  --  ROW = K                   
C       \           2N+K   DUPLICATE ENTRY IN A  --  ROW = K            
C       \           3N+K   INSUFFICIENT STORAGE FOR JL  --  ROW = K     
C       \           4N+1   INSUFFICIENT STORAGE FOR L                   
C       \           5N+K   NULL PIVOT  --  ROW = K                      
C       \           6N+K   INSUFFICIENT STORAGE FOR JU  --  ROW = K     
C       \           7N+1   INSUFFICIENT STORAGE FOR U                   
C       \           8N+K   ZERO PIVOT  --  ROW = K                      
C NVA   \ IA    - POINTERS TO DELIMIT THE ROWS IN A.                    
C       \           SIZE = N+1.                                         
C NVA   \ IC    - INVERSE OF THE ORDERING OF THE COLUMNS OF M;  I.E.,   
C       \           IC(C(I) = I  FOR I=1,...N.                          
C       \           SIZE = N.                                           
C NVRA  \ IL    - POINTERS TO DELIMIT THE ROWS IN L.                    
C       \           SIZE = N+1.                                         
C NVRA  \ IU    - POINTERS TO DELIMIT THE ROWS IN U.                    
C       \           SIZE = N+1.                                         
C NVA   \ JA    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF A.    
C       \           SIZE = SIZE OF A.                                   
C NVRA  \ JL    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF L.    
C       \           SIZE = JLMAX.                                       
C NV    \ JLMAX - DECLARED DIMENSION OF JL;  JLMAX MUST BE LARGER THAN  
C       \           THE NUMBER OF NONZERO ENTRIES IN THE STRICT LOWER   
C       \           TRIANGLE OF M PLUS FILLIN  (IL(N+1)-1 AFTER NSF).   
C NVRA  \ JU    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF U.    
C       \           SIZE = JUMAX.                                       
C NV    \ JUMAX - DECLARED DIMENSION OF JU;  JUMAX MUST BE LARGER THAN  
C       \           THE NUMBER OF NONZERO ENTRIES IN THE STRICT UPPER   
C       \           TRIANGLE OF M PLUS FILLIN  (IU(N+1)-1 AFTER NSF).   
C FVRA  \ L     - NONZERO ENTRIES IN THE STRICT LOWER TRIANGULAR PORTION
C       \           OF THE MATRIX L, STORED BY ROWS.                    
C       \           SIZE = LMAX                                         
C NV    \ LMAX  - DECLARED DIMENSION OF L;  LMAX MUST BE LARGER THAN    
C       \           THE NUMBER OF NONZERO ENTRIES IN THE STRICT LOWER   
C       \           TRIANGLE OF M PLUS FILLIN  (IL(N+1)-1 AFTER NSF).   
C NV    \ N     - NUMBER OF VARIABLES/EQUATIONS.                        
C NVA   \ R     - ORDERING OF THE ROWS OF M.                            
C       \           SIZE = N.                                           
C FVRA  \ U     - NONZERO ENTRIES IN THE STRICT UPPER TRIANGULAR PORTION
C       \           OF THE MATRIX U, STORED BY ROWS.                    
C       \           SIZE = UMAX.                                        
C NV    \ UMAX  - DECLARED DIMENSION OF U;  UMAX MUST BE LARGER THAN    
C       \           THE NUMBER OF NONZERO ENTRIES IN THE STRICT UPPER   
C       \           TRIANGLE OF M PLUS FILLIN  (IU(N+1)-1 AFTER NSF).   
C FRA   \ Z     - SOLUTION X.                                           
C       \           SIZE = N.                                           
C                                         
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   NPOIN        | -->| NOMBRE D'INCONNUES
C |   NSEGB        | -->| NOMBRE DE SEGMENTS 
C |   GLOSEG       | -->| NUMEROS GLOBAUX DES POINTS DES SEGMENTS
C |   DA,XA        | -->| DIAGONALE ET TERMES EXTRA-DIAGONAUX DE LA MATRICE
C |   XINC         |<-- | SOLUTION
C |   RHS          | -->| SECOND MEMBRE
C |   INFOGR       | -->| IF, YES INFORMATIONS ON LISTING
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_SD_SOLVE_1 => SD_SOLVE_1
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NPOIN,NSEGB,MAXSEG
      INTEGER, INTENT(IN)             :: GLOSEG(MAXSEG,2)
      LOGICAL, INTENT(IN)             :: INFOGR
      DOUBLE PRECISION, INTENT(IN)    :: XA(*),RHS(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: XINC(NPOIN),DA(NPOIN)
      CHARACTER(LEN=1), INTENT(IN)    :: TYPEXT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C     
      INTEGER MEMFACTOR,IERR,NPBLK,NSEGBLK,NSP,ESP,INDIC,FLAG,I,IERRK
C
      INTEGER, ALLOCATABLE :: INDTRI(:),INX(:),IPX(:),IW1(:)
      INTEGER, ALLOCATABLE :: IN(:),IP(:),ISP(:),ISEGIP(:)
C
      DOUBLE PRECISION, ALLOCATABLE :: AC(:),ACTRI(:),RSP(:)
C
C     GESTION DE LA TAILLE DES TABLEAUX ALLOUABLES :
C
      INTEGER SIZE_IN,SIZE_IP,SIZE_ISEGIP,SIZE_IW1,SIZE_INDTRI
      INTEGER SIZE_INX,SIZE_IPX,SIZE_AC,SIZE_ACTRI,SIZE_ISP,SIZE_RSP
C
      DATA SIZE_IN    /0/
      DATA SIZE_IP    /0/
      DATA SIZE_ISEGIP/0/
      DATA SIZE_IW1   /0/
      DATA SIZE_INDTRI/0/
      DATA SIZE_INX   /0/
      DATA SIZE_IPX   /0/
      DATA SIZE_AC    /0/
      DATA SIZE_ACTRI /0/
      DATA SIZE_ISP   /0/
      DATA SIZE_RSP   /0/
C 
      SAVE          
C
C-----------------------------------------------------------------------
C
C     CORRECTION OF DIAGONALS (TIDAL FLATS WITH MASKING)
C
      DO I=1,NPOIN
        IF(ABS(DA(I)).LT.1.D-15) DA(I)=1.D0
      ENDDO
C
C-----------------------------------------------------------------------
C    
      IF(INFOGR) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) '                       RESOLUTION DIRECTE'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) '                      DIRECT SYSTEM SOLVER'
        ENDIF
      ENDIF
C
      NPBLK=NPOIN
      NSEGBLK=NSEGB
C       
C     1. MEMFACTOR : FACTEUR MEMOIRE POUR TAILLE ISP ET RSP DANS ODRV ET SDRV
C     =======================================================================
C 
      IF(TYPEXT.EQ.'S') THEN
        MEMFACTOR = 5
      ELSE 
        MEMFACTOR = 15
      ENDIF
      NSP=MEMFACTOR*(NPBLK+4*NSEGBLK) 
      ESP=MEMFACTOR*(NPBLK+4*NSEGBLK)      
C       
C     2. ALLOCATION DES TABLEAUX (OU REALLOCATION SI TROP PETITS)
C     ======================================================================= 
C
      IF(SIZE_IN.EQ.0) THEN      
        ALLOCATE(IN(NPBLK+1))
        SIZE_IN=    NPBLK+1
      ELSEIF(       NPBLK+1.GT.SIZE_IN) THEN
        DEALLOCATE(IN)
        ALLOCATE(IN(NPBLK+1))
        SIZE_IN=    NPBLK+1
      ENDIF
C
      IF(SIZE_IP.EQ.0) THEN
        ALLOCATE(IP(NSEGBLK*2))
        SIZE_IP=    NSEGBLK*2
      ELSEIF(       NSEGBLK*2.GT.SIZE_IP) THEN
        DEALLOCATE(IP)
        ALLOCATE(IP(NSEGBLK*2))
        SIZE_IP=    NSEGBLK*2
      ENDIF
C
      IF(SIZE_ISEGIP.EQ.0) THEN
        ALLOCATE(ISEGIP(NSEGBLK*2+1))
        SIZE_ISEGIP=    NSEGBLK*2+1
      ELSEIF(           NSEGBLK*2+1.GT.SIZE_ISEGIP) THEN
        DEALLOCATE(ISEGIP)
        ALLOCATE(ISEGIP(NSEGBLK*2+1))
        SIZE_ISEGIP=    NSEGBLK*2+1
      ENDIF
C
      IF(SIZE_IW1.EQ.0) THEN
        ALLOCATE(IW1(NPBLK))
        SIZE_IW1=    NPBLK
      ELSEIF(        NPBLK.GT.SIZE_IW1) THEN
        DEALLOCATE(IW1)
        ALLOCATE(IW1(NPBLK))
        SIZE_IW1=    NPBLK
      ENDIF
C
      IF(SIZE_INDTRI.EQ.0) THEN
        ALLOCATE(INDTRI(NPBLK))
        SIZE_INDTRI=    NPBLK
      ELSEIF(           NPBLK.GT.SIZE_INDTRI) THEN
        DEALLOCATE(INDTRI)
        ALLOCATE(INDTRI(NPBLK))
        SIZE_INDTRI=    NPBLK
      ENDIF
C
      IF(SIZE_INX.EQ.0) THEN
        ALLOCATE(INX(NPBLK+1))
        SIZE_INX=    NPBLK+1
      ELSEIF(        NPBLK+1.GT.SIZE_INX) THEN
        DEALLOCATE(INX)
        ALLOCATE(INX(NPBLK+1))
        SIZE_INX=    NPBLK+1
      ENDIF
C
      IF(SIZE_IPX.EQ.0) THEN
        ALLOCATE(IPX(NSEGBLK*2+NPBLK+1))
        SIZE_IPX=    NSEGBLK*2+NPBLK+1
      ELSEIF(        NSEGBLK*2+NPBLK+1.GT.SIZE_IPX) THEN
        DEALLOCATE(IPX)
        ALLOCATE(IPX(NSEGBLK*2+NPBLK+1))
        SIZE_IPX=    NSEGBLK*2+NPBLK+1
      ENDIF
C
      IF(SIZE_AC.EQ.0) THEN
        ALLOCATE(AC(NSEGBLK*2+NPBLK+1))
        SIZE_AC=    NSEGBLK*2+NPBLK+1
      ELSEIF(       NSEGBLK*2+NPBLK+1.GT.SIZE_AC) THEN
        DEALLOCATE(AC)
        ALLOCATE(AC(NSEGBLK*2+NPBLK+1))
        SIZE_AC=    NSEGBLK*2+NPBLK+1
      ENDIF
C
      IF(SIZE_ACTRI.EQ.0) THEN
        ALLOCATE(ACTRI(NPBLK))
        SIZE_ACTRI=    NPBLK
      ELSEIF(          NPBLK.GT.SIZE_ACTRI) THEN
        DEALLOCATE(ACTRI)
        ALLOCATE(ACTRI(NPBLK))
        SIZE_ACTRI=    NPBLK
      ENDIF
C
      IF(SIZE_ISP.EQ.0) THEN
        ALLOCATE(ISP(NSP))
        SIZE_ISP=    NSP
      ELSEIF(        NSP.GT.SIZE_ISP) THEN
        DEALLOCATE(ISP)
        ALLOCATE(ISP(NSP))
        SIZE_ISP=    NSP
      ENDIF
C
      IF(SIZE_RSP.EQ.0) THEN      
        ALLOCATE(RSP(ESP))
        SIZE_RSP=    ESP
      ELSEIF(        ESP.GT.SIZE_RSP) THEN
        DEALLOCATE(RSP)
        ALLOCATE(RSP(ESP))
        SIZE_RSP=    ESP
      ENDIF     
C
C     3. CONSTRUCTION STOCKAGE COMPACT NON SYMETRIQUE (IN,IP)
C     SANS LA DIAGONALE ET (INX,IPX) AVEC DIAGONALE
C     =======================================================================        
C
      CALL SD_STRSSD(NPBLK,NSEGBLK,GLOSEG(1,1),GLOSEG(1,2),
     *               IN,IP,ISEGIP,IW1)   
C
      IF(TYPEXT.EQ.'S') THEN 
        CALL SD_FABCAD(NPBLK,NSEGBLK,IN,IP,ISEGIP,
     *                 INDTRI,IW1,INX,IPX,ACTRI,XA,XA,DA,AC)
C                             ISTRI
      ELSE
        CALL SD_FABCAD(NPBLK,NSEGBLK,IN,IP,ISEGIP,
     *                 INDTRI,IW1,INX,IPX,ACTRI,XA,XA(NSEGBLK+1),DA,AC) 
      ENDIF     
C       
C     4. PERMUTATION MINIMUM DEGRE (YSMP PACKAGE)
C     ======================================================================= 
C      
      INDIC=1
C       
      CALL SD_ODRV(NPBLK,INX,IPX,AC,IN  ,IW1 ,NSP,ISP,INDIC,FLAG)
C                                   PERM,INVP
C
      IF(FLAG.NE.0) THEN
        IF(LNG.EQ.1) THEN
           WRITE(LU,*) 'AUGMENTER LE FACTEUR MEMOIRE (MEMFACTOR)', 
     *                 ' DANS LA ROUTINE SD_SOLVE_1'
        ELSEIF(LNG.EQ.2)THEN
           WRITE(LU,*) 'INCREASE THE MEMORY FACTOR (MEMFACTOR)', 
     *                 ' IN THE ROUTINE SD_SOLVE_1'
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
C     
C---> SECOND MEMBRE DU SYSTEME 
C         
C     5. DECOMPOSITION LDLT ET RESOLUTION (YSMP PACKAGE)
C
      IF(TYPEXT.EQ.'S') THEN
C                          PERM,INVP 
        CALL SD_SDRV(NPBLK,IN  ,IW1 ,INX,IPX,AC,RHS,XINC, 
     *               NSP,ISP,RSP,ESP,INDIC,FLAG)
      ELSE
        CALL SD_NDRV(NPBLK,IN  ,IN  ,IW1,INX,IPX,AC,RHS,XINC, 
     *               NSP,ISP,RSP,ESP,INDIC,FLAG)
      ENDIF
C
      IF(TYPEXT.EQ.'S') THEN
      IF(FLAG.NE.0) THEN
        IERR=FLAG-8*NPBLK
        IF(IERR.GT.0) THEN
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'MATRICE AVEC PIVOT NUL A LA LIGNE'
            WRITE(LU,*) IERR 
          ELSEIF(LNG.EQ.2) THEN
            WRITE(LU,*) 'MATRIX WITH ZERO PIVOT AT ROW'
            WRITE(LU,*) IERR 
          ENDIF
        ELSE 
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'AUGMENTER LE FACTEUR MEMOIRE (MEMFACTOR)', 
     *                  ' DE 1 DANS LA ROUTINE SD_SOLVE_1'
          ELSEIF(LNG.EQ.2) THEN
            WRITE(LU,*) 'ADD 1 TO THE MEMORY FACTOR (MEMFACTOR)', 
     *                  ' IN SUBROUTINE SD_SOLVE_1'
          ENDIF
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
      ELSE
C    
C---> COMMENTAIRES DE L'ERREUR : FLAG_SD_NDRV :
C    
C FLAG  - ERROR FLAG;  VALUES AND THEIR MEANINGS ARE --         
C             0     NO ERRORS DETECTED                          
C             N+K   NULL ROW IN A  --  ROW = K                  
C            2N+K   DUPLICATE ENTRY IN A  --  ROW = K           
C            3N+K   INSUFFICIENT STORAGE IN NSF  --  ROW = K    
C            4N+1   INSUFFICIENT STORAGE IN NNF                 
C            5N+K   NULL PIVOT  --  ROW = K                     
C            6N+K   INSUFFICIENT STORAGE IN NSF  --  ROW = K    
C            7N+1   INSUFFICIENT STORAGE IN NNF                 
C            8N+K   ZERO PIVOT  --  ROW = K                     
C           10N+1   INSUFFICIENT STORAGE IN NDRV                
C           11N+1   ILLEGAL PATH SPECIFICATION (INDIC)
      IF(FLAG.NE.0) THEN
        IERR=INT(FLAG/NPBLK)	
        IF(IERR.EQ.3.OR.IERR.EQ.5.OR.IERR.EQ.8) THEN
          IERRK=FLAG-IERR*NPBLK
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'MATRICE AVEC PIVOT NUL A LA LIGNE'
            WRITE(LU,*) IERRK 
          ELSEIF(LNG.EQ.2) THEN
            WRITE(LU,*) 'MATRIX WITH ZERO PIVOT AT ROW'
            WRITE(LU,*) IERRK 
          ENDIF
        ELSE 
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'AUGMENTER LE FACTEUR MEMOIRE (MEMFACTOR)', 
     *                  ' DANS LE SOUS-PROGRAMME SD_SOLVE_1'
          ELSEIF(LNG.EQ.2) THEN
            WRITE(LU,*) 'INCREASE THE MEMORY FACTOR (MEMFACTOR)', 
     *                  ' IN SUBROUTINE SD_SOLVE_1'
          ENDIF
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
