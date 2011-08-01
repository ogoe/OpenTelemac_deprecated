C                         *****************
                          SUBROUTINE SD_MDU
C                         *****************
C
     *(EK,DMIN,V,L,HEAD,LAST,NEXT,MARK)
C
C***********************************************************************
C BIEF VERSION 5.7     20/11/06   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C
C  FONCTION : UPDATE DEGREES OF UNELIMINATED VERTICES IN EK 
C              
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
      USE BIEF, EX_SD_MDU => SD_MDU
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    ::  EK,V(*),L(*)
      INTEGER, INTENT(INOUT) ::  DMIN,HEAD(*),LAST(*),NEXT(*),MARK(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                          
      INTEGER TAG,VI,EVI,DVI,S,VS,ES,B,VB,ILP,ILPMAX,BLP,BLPMAX,I                                                   
      EQUIVALENCE (VS,ES)                                           
C                                                                       
C----INITIALIZE TAG 
C                                                    
      TAG = MARK(EK) - LAST(EK)                                       
C                                                                       
C----FOR EACH VERTEX VI IN EK
C                                           
      I = EK                                                          
      ILPMAX = LAST(EK)                                               
      IF(ILPMAX.LE.0)  GO TO 11                                      
      DO 10 ILP=1,ILPMAX                                              
        I = L(I)                                                      
        VI = V(I)                                                     
        IF (LAST(VI))  1, 10, 8                                       
C                                                                       
C------IF VI NEITHER PROTOTYPE NOR DUPLICATE VERTEX, THEN MERGE ELEMENTS
C------TO COMPUTE DEGREE
C                                                
1       TAG = TAG + 1                                               
        DVI = LAST(EK)                                              
C                                                                       
C--------FOR EACH VERTEX/ELEMENT VS/ES IN ELEMENT LIST OF VI
C            
        S = L(VI)                                                   
2       S = L(S)                                                    
        IF(S.EQ.0) GO TO 9                                        
        VS = V(S)                                                 
        IF(NEXT(VS).LT.0)  GO TO 3                               
C                                                                       
C----------IF VS IS UNELIMINATED VERTEX, THEN TAG AND ADJUST DEGREE
C     
        MARK(VS) = TAG                                          
        DVI = DVI + 1                                           
        GO TO 5                                                 
C                                                                       
C----------IF ES IS ACTIVE ELEMENT, THEN EXPAND                         
C------------CHECK FOR OUTMATCHED VERTEX
C                                
3       IF(MARK(ES).LT.0)  GO TO 6                             
C                                                                       
C------------FOR EACH VERTEX VB IN ES
C                                   
        B = ES                                                  
        BLPMAX = LAST(ES)                                       
        DO 4 BLP=1,BLPMAX                                       
          B = L(B)                                              
          VB = V(B)                                             
C                                                                       
C--------------IF VB IS UNTAGGED, THEN TAG AND ADJUST DEGREE
C            
          IF(MARK(VB).GE.TAG)  GO TO 4                         
          MARK(VB) = TAG                                      
          DVI = DVI + 1                                       
4       CONTINUE                                             
C                                                                       
5       GO TO 2                                                   
C                                                                       
C------ELSE IF VI IS OUTMATCHED VERTEX, THEN ADJUST OVERLAPS BUT DO NOT 
C------COMPUTE DEGREE
C                                                   
6       LAST(VI) = 0                                                
        MARK(ES) = MARK(ES) - 1                                     
7       S = L(S)                                                    
        IF(S.EQ.0)  GO TO 10                                       
        ES = V(S)                                                 
        IF(MARK(ES).LT.0)  MARK(ES) = MARK(ES) - 1               
        GO TO 7                                                   
C                                                                       
C------ELSE IF VI IS PROTOTYPE VERTEX, THEN CALCULATE DEGREE BY         
C------INCLUSION/EXCLUSION AND RESET OVERLAP COUNT
C                      
8       EVI = LAST(VI)                                              
        DVI = LAST(EK) + LAST(EVI) + MARK(EVI)                      
        MARK(EVI) = 0                                               
C                                                                       
C------INSERT VI IN APPROPRIATE DEGREE LIST
C                             
9       NEXT(VI)  = HEAD(DVI)                                          
        HEAD(DVI) = VI                                                
        LAST(VI)  = -DVI                                               
        IF(NEXT(VI).GT.0)  LAST(NEXT(VI)) = VI                       
        IF(DVI.LT.DMIN)  DMIN = DVI                                  
C                                                                       
10    CONTINUE                                                      
C                                                                       
11    CONTINUE
C                                                                       
C-----------------------------------------------------------------------
C 
      RETURN                                                          
      END
