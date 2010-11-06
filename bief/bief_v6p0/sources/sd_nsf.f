C                          *****************                                                                          
                           SUBROUTINE SD_NSF
C                          *****************
C
     *(N,R,IC,IA,JA,IL,JL,JLMAX,IU,JU,JUMAX,Q,IM,FLAG)      
C
C***********************************************************************
C BIEF VERSION 5.9     18/02/08   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C
C  FUNCTION: SYMBOLIC LDU-FACTORIZATION OF A NONSYMMETRIC SPARSE MATRIX         
C            (UNCOMPRESSED POINTER STORAGE)
C                                                                        
C       INPUT VARIABLES:   N, R,IC, IA,JA, JLMAX, JUMAX.                
C       OUTPUT VARIABLES:  IL,JL, IU,JU, FLAG.                          
C                                                                       
C       PARAMETERS USED INTERNALLY:                                     
C NIA   \ Q     - SUPPOSE M' IS THE RESULT OF REORDERING M;  IF         
C       \           PROCESSING OF THE KTH ROW OF M' (HENCE THE KTH ROWS 
C       \           OF L AND U) IS BEING DONE, THEN Q(J) IS INITIALLY   
C       \           NONZERO IF M'(K,J) IS NONZERO;  SINCE VALUES NEED   
C       \           NOT BE STORED, EACH ENTRY POINTS TO THE NEXT        
C       \           NONZERO;  FOR EXAMPLE, IF  N=9  AND THE 5TH ROW OF  
C       \           M' IS                                               
C       \                   0 X X 0 X 0 0 X 0,                          
C       \           THEN Q WILL INITIALLY BE                            
C       \                   A 3 5 A 8 A A 10 A 2        (A - ARBITRARY);
C       \           Q(N+1) POINTS TO THE FIRST NONZERO IN THE ROW AND   
C       \           THE LAST NONZERO POINTS TO  N+1;  AS THE ALGORITHM  
C       \           PROCEEDS, OTHER ELEMENTS OF Q ARE INSERTED IN THE   
C       \           LIST BECAUSE OF FILLIN.                             
C       \           SIZE = N+1.                                         
C NIA   \ IM    - AT EACH STEP IN THE FACTORIZATION, IM(I) IS THE LAST  
C       \           ELEMENT IN THE ITH ROW OF U WHICH NEEDS TO BE       
C       \           CONSIDERED IN COMPUTING FILLIN.                     
C       \           SIZE = N.                                           
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
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+                 
C 
      INTEGER N                                                                      
      INTEGER R(*),IC(*),IA(*),JA(*),IL(*),JL(*)             
      INTEGER IU(*),JU(*),Q(*),IM(*),FLAG,QM,VJ                   
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                               
C     JLPTR - POINTS TO THE LAST POSITION USED IN  JL.                   
C     JUPTR - POINTS TO THE LAST POSITION USED IN  JU.
C  
      INTEGER JLPTR,JUPTR,JLMAX,JUMAX,I,J,K,M,JMIN,JMAX
C
C-----------------------------------------------------------------------
C                                                                       
C     INITIALIZE POINTERS
C
      JLPTR = 0                                                       
      IL(1) = 1                                                       
      JUPTR = 0                                                       
      IU(1) = 1                                                       
C                                                                       
C     FOR EACH ROW OF L AND U  
C
      DO 10 K=1,N                                                     
C       SET Q TO THE REORDERED ROW OF A 
        Q(N+1) = N+1                                                  
        JMIN = IA(R(K))                                               
        JMAX = IA(R(K)+1) - 1                                         
        IF(JMIN.GT.JMAX) GO TO 101                                  
        DO J=JMIN,JMAX                                              
          VJ = IC(JA(J))                                              
          QM = N+1
1         CONTINUE                                                    
          M = QM                                                      
          QM = Q(M)                                                   
          IF(QM.LT.VJ) GO TO 1                                      
          IF(QM.EQ.VJ) GO TO 102                                    
          Q(M) = VJ                                                 
          Q(VJ) = QM                                                
        ENDDO                                                    
C                                                                       
C       FOR EACH ENTRY IN THE LOWER TRIANGLE 
C
        I = N+1 
3       CONTINUE                                                      
        I = Q(I)                                                      
        IF(I.GE.K) GO TO 7                                          
C       L(K,I) WILL BE NONZERO, SO ADD IT TO JL  
        JLPTR = JLPTR+1                                             
        IF(JLPTR.GT.JLMAX) GO TO 103                              
        JL(JLPTR) = I                                               
        QM = I                                                      
C       INSPECT ITH ROW FOR FILLIN, ADJUST IM IF POSSIBLE 
        JMIN = IU(I)                                                
        JMAX = IM(I)                                                
        IF(JMIN.GT.JMAX)  GO TO 3                                  
        DO 5 J=JMIN,JMAX                                            
          VJ = JU(J)                                                
          IF (VJ.EQ.K)  IM(I) = J 
4         CONTINUE                                  
          M = QM                                                    
          QM = Q(M)                                                 
          IF (QM.LT.VJ)  GO TO 4                                    
          IF (QM.EQ.VJ)  GO TO 5                                    
          Q(M) = VJ                                               
          Q(VJ) = QM                                              
          QM = VJ                                                 
5       CONTINUE                                                  
        GO TO 3                                                     
C                                                                       
C       CHECK FOR NULL PIVOT
C
7       CONTINUE
        IF(I.NE.K) GO TO 105                                        
C       REMAINING ELEMENTS OF Q DEFINE STRUCTURE OF U(K, )
8       CONTINUE
        I = Q(I)                                                      
        IF(I.GT.N) GO TO 9                                          
        JUPTR = JUPTR+1                                             
        IF (JUPTR.GT.JUMAX)  GO TO 106                              
        JU(JUPTR) = I                                               
        GO TO 8                                                     
C       GET READY FOR NEXT ROW 
9       CONTINUE
        IM(K) = JUPTR                                                 
        IL(K+1) = JLPTR+1                                             
        IU(K+1) = JUPTR+1 
10    CONTINUE                                            
C                                                                       
        FLAG = 0                                                        
        RETURN                                                          
C                                                                       
C ** ERROR:  NULL ROW IN A                                              
 101    FLAG = N + R(K)                                                 
        RETURN                                                          
C ** ERROR:  DUPLICATE ENTRY IN A                                       
 102    FLAG = 2*N + R(K)                                               
        RETURN                                                          
C ** ERROR:  INSUFFICIENT STORAGE FOR JL                                
 103    FLAG = 3*N + K                                                  
        RETURN                                                          
C ** ERROR:  NULL PIVOT                                                 
 105    FLAG = 5*N + K                                                  
        RETURN                                                          
C ** ERROR:  INSUFFICIENT STORAGE FOR JU                                
 106    FLAG = 6*N + K 
C
C-----------------------------------------------------------------------
C                                                   
      RETURN                                                          
      END
