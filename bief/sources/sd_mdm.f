C                        *****************
                         SUBROUTINE SD_MDM
C                        *****************
C
     *(VK,TAIL,V,L,LAST,NEXT,MARK)
C
C***********************************************************************
C BIEF VERSION 5.7     20/11/06   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C
C  FONCTION : FORM ELEMENT FROM UNELIMINATED NEIGHBORS OF VK
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
      USE BIEF, EX_SD_MDM => SD_MDM
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: VK,LAST(*),NEXT(*),V(*)
      INTEGER, INTENT(INOUT) :: TAIL,L(*),MARK(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                               
      INTEGER TAG,S,LS,VS,ES,B,LB,VB,BLP,BLPMAX                         
      EQUIVALENCE (VS,ES)                                           
C                                                                       
C----INITIALIZE TAG AND LIST OF UNELIMINATED NEIGHBORS
C                  
      TAG = MARK(VK)                                                  
      TAIL = VK                                                       
C                                                                       
C----FOR EACH VERTEX/ELEMENT VS/ES IN ELEMENT LIST OF VK
C                
      LS = L(VK)                                                      
1     S = LS                                                          
      IF(S.EQ.0)  GO TO 5                                            
      LS = L(S)                                                     
      VS = V(S)                                                     
      IF(NEXT(VS).LT.0)  GO TO 2                                   
C                                                                       
C------IF VS IS UNELIMINATED VERTEX, THEN TAG AND APPEND TO LIST OF     
C------UNELIMINATED NEIGHBORS
C                                           
      MARK(VS) = TAG                                              
      L(TAIL) = S                                                 
      TAIL = S                                                    
      GO TO 4                                                     
C                                                                       
C------IF ES IS ACTIVE ELEMENT, THEN ...                                
C--------FOR EACH VERTEX VB IN BOUNDARY LIST OF ELEMENT ES
C              
2     LB = L(ES)                                                  
      BLPMAX = LAST(ES)                                           
      DO 3 BLP=1,BLPMAX                                           
        B = LB                                                    
        LB = L(B)                                                 
        VB = V(B)                                                 
C                                                                       
C----------IF VB IS UNTAGGED VERTEX, THEN TAG AND APPEND TO LIST OF     
C----------UNELIMINATED NEIGHBORS
C                                       
        IF(MARK(VB).GE.TAG)  GO TO 3                             
        MARK(VB) = TAG                                          
        L(TAIL) = B                                             
        TAIL = B                                                
3     CONTINUE                                                  
C                                                                       
C--------MARK ES INACTIVE
C                                               
      MARK(ES) = TAG                                              
C                                                                       
4     GO TO 1                                                       
C                                                                       
C----TERMINATE LIST OF UNELIMINATED NEIGHBORS
C                           
5     L(TAIL) = 0                                                     
C                                                                       
C-----------------------------------------------------------------------
C                                                                        
      RETURN                                                          
      END
