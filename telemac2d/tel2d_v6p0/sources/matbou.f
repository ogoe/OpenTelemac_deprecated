C                       *****************                               
                        SUBROUTINE MATBOU                               
C                       *****************                               
C                                                                       
     *(MESH,M1,M2,A11,A12,A21,A22,SMU,SMV,VR,VS,H0,MSK,MASKEL,S)         
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
C |   VRN,VSN      |<-->|                                               
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
C                                                                       
C  STRUCTURES DE MATRICES                                               
C                                                                       
      TYPE(BIEF_OBJ), INTENT(INOUT) :: A11,A12,A21,A22,M1,M2
C                                                                                                                                              
C  STRUCTURE DE MAILLAGE                                                
C                                                                       
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C                                                                       
C  STRUCTURES DE VECTEURS                                               
C                                                                       
      TYPE(BIEF_OBJ), INTENT(INOUT) :: SMU,SMV
      TYPE(BIEF_OBJ), INTENT(IN)    :: MASKEL,H0,S,VR,VS    
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                      
      INTEGER IELMU,IELMH                                               
C                                                                                                                                             
      DOUBLE PRECISION SL1,SL11,C                                              
C                                                                       
      CHARACTER*16 FORMUL                                               
C                                                                                                                                            
C-----------------------------------------------------------------------                                                                         
C                                                                   
      IELMH=H0%ELM                                              
      IELMU=VR%ELM                                         
C                                                                                                                                                                                                      
C     MATRICE DE MASSE (POUR CALCUL DE A11 ET DU SECOND MEMBRE SMU)                                          
C                                                                       
      FORMUL='MATMAS          '
      SL1 = 1.D0                                         
      CALL MATRIX(M1,'M=N     ',FORMUL,IELMU,IELMU,                    
     *            SL1,S,S,S,S,S,S,MESH,MSK,MASKEL)                             
C
C------------------------------------------------------------------     
C
C     SECOND MEMBRE DE L'EQUATION DE U
C
      CALL MATVEC( 'X=AY    ',SMU,M1,VR,C,MESH)  
C
C------------------------------------------------------------------     
C
C     MATRICE DE U DANS L'EQUATION DE U (MATRICE DE MASSE FAITE
C     PLUS HAUT)
C                                                                                                                                                                    
      FORMUL='FFBT        0XX0'                                                                                                                                                                                     
      SL11 = 1.D0 / 6.D0                                                                                                                     
      CALL MATRIX(A11,'M=N     ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                                
C                                                                               
      FORMUL='FFBT        0YY0'                                                                                                                
      SL11 = 2.D0 / 3.D0                                                                                                                      
      CALL MATRIX(A11,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                                
C                                                                      
      FORMUL='FFBT        XX00'                                                                                                                
      SL11 = 1.D0 / 2.D0                                                                                                                      
      CALL MATRIX(A11,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                                
C                                                                       
C     FORMUL='FFBT        0X0X' 
      FORMUL='FFBT        0XX0'                                                                                                                 
      SL11 = 1.D0 / 2.D0                                                                                                                     
      CALL MATRIX(A11,'M=M+TN  ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C
      FORMUL='FFBT   00XX+00YY'                                                                                                               
      SL11 = 1.D0 / 3.D0                                                                                                                     
      CALL MATRIX(A11,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C
C     AJOUT DE LA MATRICE DE MASSE
C
      CALL OM( 'M=M+N   ' , A11 , M1 , S , C , MESH )  
C                           
C------------------------------------------------------------------     
C                                                                                                                                                                              
C     MATRICE DE V DANS L'EQUATION DE U                                     
C
      FORMUL='FFBT        0Y0X'                                                                                                               
      SL11 = 1.D0 / 2.D0                                                                                                                      
      CALL MATRIX(A12,'M=N     ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C
      FORMUL='FFBT        XY00'                                                                                                               
      SL11 = 1.D0 / 2.D0                                                                                                                      
      CALL MATRIX(A12,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C                                                                       
C     FORMUL='FFBT        0XY0'
      FORMUL='FFBT        0X0Y'                                                                                                               
      SL11 = -1.D0 / 2.D0                                                                                                                      
      CALL MATRIX(A12,'M=M+TN  ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                              
C                                                                                                                               
C------------------------------------------------------------------     
C                                                                                                                                                                                                      
C     MATRICE DE MASSE (POUR CALCUL DE A22 ET DU SECOND MEMBRE SMV)                                          
C                                                                       
      FORMUL='MATMAS          '
      SL1 = 1.D0                                         
      CALL MATRIX(M2,'M=N     ',FORMUL,IELMU,IELMU,                    
     *            SL1,S,S,S,S,S,S,MESH,MSK,MASKEL)
C                                                                                                                                                          
C     SECOND MEMBRE SMV                                       
C                                                                       
      CALL MATVEC( 'X=AY    ',SMV,M2,VS,C,MESH)
C                                 
C------------------------------------------------------------------     
C                                                                                                                                              
C     MATRICE DE V DANS L'EQUATION DE V  
C                                                                                                                                                                                                                                           
      FORMUL='FFBT        0YY0'                                                                                                            
      SL11 = 1.D0 / 6.D0                                                                                                                  
      CALL MATRIX(A22,'M=N     ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                                
C                                                                                 
      FORMUL='FFBT        0XX0'                                                                                                               
      SL11 = 2.D0 / 3.D0                                                                                                                     
      CALL MATRIX(A22,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                                
C                                                                       
      FORMUL='FFBT        YY00'                                                                                                             
      SL11 = 1.D0 / 2.D0                                                                                                                     
      CALL MATRIX(A22,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                                
C                                                                       
C     FORMUL='FFBT        0Y0Y'
      FORMUL='FFBT        0YY0'                                                                                                             
      SL11 = 1.D0 / 2.D0                                                                                                                  
      CALL MATRIX(A22,'M=M+TN  ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                                
C 
      FORMUL='FFBT   00XX+00YY'                                                                                                               
      SL11 = 1.D0 / 3.D0                                                                                                                       
      CALL MATRIX(A22,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C
C     AJOUT DE LA MATRICE DE MASSE
C
      CALL OM( 'M=M+N   ' , A22 , M2 , S , C , MESH )                                  
C
C------------------------------------------------------------------     
C                                                                                                                                            
C     MATRICE DE U DANS L'EQUATION DE V                            
C
      FORMUL='FFBT        0X0Y'                                                                                                               
      SL11 = 1.D0 / 2.D0                                                                                                                      
      CALL MATRIX(A21,'M=N     ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C 
      FORMUL='FFBT        XY00'                                                                                                                
      SL11 = 1.D0 / 2.D0                                                                                                                      
      CALL MATRIX(A21,'M=M+N   ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)
C                                                                       
C     FORMUL='FFBT        0YX0'    
      FORMUL='FFBT        0Y0X'                                                                                                               
      SL11 = -1.D0 / 2.D0                                                                                                                      
      CALL MATRIX(A21,'M=M+TN  ',FORMUL,IELMU,IELMU,                    
     *            SL11,H0,S,S,S,S,S,MESH,MSK,MASKEL)                                  
C                                                         
C------------------------------------------------------------------     
C                                                                             
      RETURN                                                            
      END                       
