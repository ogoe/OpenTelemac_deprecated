C                         *****************
                          SUBROUTINE SD_MDI
C                         *****************
C
     *(N,IA,JA,MAX,V,L,HEAD,LAST,NEXT,MARK,TAG,FLAG)
C
C***********************************************************************
C BIEF VERSION 5.7     20/11/06   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C
C  FONCTION : INITIALIZATION 
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
      USE BIEF, EX_SD_MDI => SD_MDI
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: N,MAX,IA(*),JA(*)
      INTEGER, INTENT(INOUT) :: V(*),L(*),HEAD(*),LAST(*) 
      INTEGER, INTENT(INOUT) :: NEXT(*),MARK(*),TAG,FLAG
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C           
      INTEGER SFS,VI,DVI,VJ,JMIN,JMAX,J                        
C                                                                       
C----INITIALIZE DEGREES, ELEMENT LISTS, AND DEGREE LISTS
C                
      DO 1 VI=1,N                                                     
        MARK(VI) = 1                                                  
        L(VI) = 0                                                     
        HEAD(VI) = 0 
1     CONTINUE                                                 
      SFS = N+1                                                       
C                                                                       
C----CREATE NONZERO STRUCTURE                                           
C----FOR EACH NONZERO ENTRY A(VI,VJ) IN STRICT UPPER TRIANGLE
C           
      DO 3 VI=1,N                                                     
        JMIN = IA(VI)                                                 
        JMAX = IA(VI+1) - 1                                           
        IF(JMIN.GT.JMAX)  GO TO 3                                    
        DO 2 J=JMIN,JMAX                                              
          VJ = JA(J)                                                  
          IF(VI.GE.VJ) GO TO 2                                      
          IF(SFS.GE.MAX) GO TO 101                                
C                                                                       
C------ENTER VJ IN ELEMENT LIST FOR VI
C                                  
          MARK(VI) = MARK(VI) + 1                                   
          V(SFS) = VJ                                               
          L(SFS) = L(VI)                                            
          L(VI) = SFS                                               
          SFS = SFS+1                                               
C                                                                       
C------ENTER VI IN ELEMENT LIST FOR VJ
C                                  
          MARK(VJ) = MARK(VJ) + 1                                   
          V(SFS) = VI                                               
          L(SFS) = L(VJ)                                            
          L(VJ) = SFS                                               
          SFS = SFS+1                                               
2       CONTINUE                                                    
3     CONTINUE                                                        
C                                                                       
C----CREATE DEGREE LISTS AND INITIALIZE MARK VECTOR
C                     
      DO 4 VI=1,N                                                     
        DVI = MARK(VI)                                                
        NEXT(VI) = HEAD(DVI)                                          
        HEAD(DVI) = VI                                                
        LAST(VI) = -DVI                                               
        IF(NEXT(VI).GT.0)  LAST(NEXT(VI)) = VI                       
        MARK(VI) = TAG
4     CONTINUE                                                
C                                                                       
      RETURN                                                          
C                                                                       
C ** ERROR -- INSUFFICIENT STORAGE
C                                      
101   FLAG = 9*N + VI
C                                                                       
C-----------------------------------------------------------------------
C                                                  
      RETURN                                                          
      END
