C                       *****************                               
                        SUBROUTINE ROTNE0                               
C                       *****************                               
C                                                                       
     *(MESH,M1,A11,A12,A21,A22,SMU,SMV,UN,VN,H0,MSK,MASKEL,S,DT)         
C                                                                       
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    17/08/94    J-M HERVOUET (LNH) 30 87 80 18 
C                                          C MOULIN   (LNH) 30 87 83 81 
C***********************************************************************
C                                                                       
C      FONCTION:                                                        
C      =========                                                        
C                                                                       
C      CALCUL DES MATRICES POUR LA RESOLUTION DES EQUATIONS D'HELMHOLTZ
C      (ETAPES 1 ET 3 DE L'ALGORITHME POUR BOUSSINESQ)
C      
C      A11 EST LA MATRICE QUI MULTIPLIE U DANS L'EQUATION DE U
C
C      A12 EST LA MATRICE QUI MULTIPLIE V DANS L'EQUATION DE U                                                                          
C
C      A21 EST LA MATRICE QUI MULTIPLIE U DANS L'EQUATION DE V
C
C      A22 EST LA MATRICE QUI MULTIPLIE V DANS L'EQUATION DE V 
C
C      SMU EST LE SECOND MEMBRE DE L'EQUATION DE U
C
C      SMV EST LE SECOND MEMBRE DE L'EQUATION DE V
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   U ,V ,H      |<-- |  VALEURS A L' ETAPE N+1.                     |
C |   UN,VN,HN     | -->|  VALEURS A L' ETAPE N.                        
C |   VR,V.,H.     | -->|  VALEURS APRES LA CONVECTION.                 
C |                | -- |                                               
C |   ZF           | -->|  COTE DU FONT AU NOEUD DE MAILLAGE .          
C |   AM1,2,3      |<-->|  MATRICES                                     
C |   TM1          |<-->|  MATRICE                                      
C |   SMU,SMV      |<-- |  SECONDS MEMBRES DU SYSTEME.                                                      
C |   AT,DT        | -->|  TEMPS, PAS DE TEMPS      
C |---------------------------------------------------------------------
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C APPELE PAR : TELMAC                                                   
C                                                                       
C SOUS-PROGRAMMES APPELES :                                             
C                                                                       
C********************************************************************** 
C
      USE BIEF
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      LOGICAL, INTENT(IN) :: MSK  
      DOUBLE PRECISION, INTENT(IN)   :: DT
      TYPE(BIEF_OBJ), INTENT(IN)     :: MASKEL,H0,S,UN,VN
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: SMU,SMV
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: A11,A12,A21,A22,M1
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                       
      INTEGER IELMU,IELMH                                               
C                                                                                                                                             
      DOUBLE PRECISION SL11,C,SURDT                                                    
C                                                                       
      CHARACTER*16 FORMUL                                                   
C                                                                       
C-----------------------------------------------------------------------
C                                                                    
      IELMH=H0%ELM                                              
      IELMU=UN%ELM                                 
C
C------------------------------------------------------------------     
C
      SURDT = 1.D0 / DT
C
C     MATRICE DE U DANS L'EQUATION DE U (STOCKEE D'ABORD DANS M1)
C                                                                                                                                                                    
      FORMUL='FFBT        0XX0'                                                                                                                                                                                     
      SL11 = SURDT / 6.D0                                                                                                                     
      CALL MATRIX(M1,'M=N     ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                                
C                                                                                                                                                    
      FORMUL='FFBT        XX00'                                                                                                                
      SL11 = SURDT / 2.D0                                                                                                                      
      CALL MATRIX(M1,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                                
C                                                                       
C     FORMUL='FFBT        0X0X' 
      FORMUL='FFBT        0XX0'                                                                                                                 
      SL11 = SURDT / 2.D0                                                                                                                     
      CALL MATRIX(M1,'M=M+TN  ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C
      FORMUL='FFBT        00XX'                                                                                                               
      SL11 = SURDT / 3.D0                                                                                                                     
      CALL MATRIX(M1,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C
C     AJOUT A LA MATRICE A11
C
      IF(A11%TYPEXT.EQ.'S') CALL OM( 'M=X(M)  ' ,A11,A11,S,C,MESH)
      CALL OM( 'M=M+N   ' , A11 , M1 , S , C , MESH )  
C                                                                                                                                                         
C     SECOND MEMBRE SMU                                       
C                                                                       
      CALL MATVEC( 'X=X+AY  ',SMU,M1,UN,C,MESH) 
C                          
C------------------------------------------------------------------     
C                                                                                                                                                                              
C     MATRICE DE V DANS L'EQUATION DE U (STOCKEE D'ABORD DANS M1)                                     
C
      FORMUL='FFBT        0Y0X'                                                                                                               
      SL11 = SURDT / 2.D0                                                                                                                      
      CALL MATRIX(M1,'M=N     ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C
      FORMUL='FFBT        XY00'                                                                                                               
      SL11 = SURDT / 2.D0                                                                                                                      
      CALL MATRIX(M1,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C
      FORMUL='FFBT        00XY'                                                                                                               
      SL11 = SURDT / 3.D0                                                                                                                      
      CALL MATRIX(M1,'M=M+TN  ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C                                                                       
      FORMUL='FFBT        0X0Y'                                                                                                               
      SL11 = SURDT / 6.D0                                                                                                                      
      CALL MATRIX(M1,'M=M+TN  ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                              
C
      CALL OM( 'M=N     ' , A12 , M1 , S , C , MESH )  
C                                                                                                                                                         
C     SECOND MEMBRE SMU                                       
C                                                                       
      CALL MATVEC( 'X=X+AY  ',SMU,M1,VN,C,MESH) 
C                                                                                                                               
C------------------------------------------------------------------
C                                                                                                                                              
C     MATRICE DE V DANS L'EQUATION DE V (STOCKEE D'ABORD DANS M1)  
C                                                                                                                                                                                                                                           
      FORMUL='FFBT        0YY0'                                                                                                            
      SL11 = SURDT / 6.D0                                                                                                                  
      CALL MATRIX(M1,'M=N     ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                                
C                                                                                                                                                     
      FORMUL='FFBT        YY00'                                                                                                             
      SL11 = SURDT / 2.D0                                                                                                                     
      CALL MATRIX(M1,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                                
C                                                                       
      FORMUL='FFBT        0YY0'                                                                                                             
      SL11 = SURDT / 2.D0                                                                                                                  
      CALL MATRIX(M1,'M=M+TN  ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                                
C 
      FORMUL='FFBT        00YY'                                                                                                               
      SL11 = SURDT / 3.D0                                                                                                                       
      CALL MATRIX(M1,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C
C     AJOUT A LA MATRICE A22
C
      IF(A22%TYPEXT.EQ.'S') CALL OM( 'M=X(M)  ' ,A22,A22,S,C,MESH)
      CALL OM( 'M=M+N   ' , A22 , M1 , S , C , MESH )                                  
C                                                                                                                                                         
C     SECOND MEMBRE SMV                                       
C                                                                       
      CALL MATVEC( 'X=X+AY  ',SMV,M1,VN,C,MESH) 
C
C------------------------------------------------------------------     
C                                                                                                                                            
C     MATRICE DE U DANS L'EQUATION DE V (STOCKEE D'ABORD DANS M1)                            
C
      FORMUL='FFBT        0X0Y'                                                                                                               
      SL11 = SURDT / 2.D0                                                                                                                      
      CALL MATRIX(M1,'M=N     ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C 
      FORMUL='FFBT        00XY'                                                                                                                
      SL11 = SURDT / 3.D0                                                                                                                      
      CALL MATRIX(M1,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C
      FORMUL='FFBT        XY00'                                                                                                                
      SL11 = SURDT / 2.D0                                                                                                                      
      CALL MATRIX(M1,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C                                                                          
      FORMUL='FFBT        0Y0X'                                                                                                               
      SL11 = SURDT / 6.D0                                                                                                                      
      CALL MATRIX(M1,'M=M+TN  ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                                  
C
      CALL OM( 'M=N     ' , A21 , M1 , S , C , MESH ) 
C                                                                                                                                                        
C     SECOND MEMBRE SMV                                       
C                                                                       
      CALL MATVEC( 'X=X+AY  ',SMV,M1,UN,C,MESH) 
C                                                        
C------------------------------------------------------------------     
C                                                                             
      RETURN                                                            
      END     
