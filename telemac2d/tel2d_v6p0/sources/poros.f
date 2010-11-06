C                       ****************                               
                        SUBROUTINE POROS                               
C                       ****************                               
C                                                                       
     *(TETA,ZF,HN,MESH)                                                 
C                                                                       
C***********************************************************************
C TELEMAC 2D VERSION 5.2     01/08/97    J-M HERVOUET (LNH) 30 87 80 18
C                                  PAUL BATES (BRISTOL) 44 117 928 9108 
C***********************************************************************
C                                                                       
C FONCTION : MARQUAGE DES BANCS DECOUVRANTS                             
C                                                                       
C MODIFIED TO IMPLEMENT WETTING/DRYING ALGORITHM
C OF DELFINA ET AL                          
C                                                                       
C            PARTIALLY WET ELEMENT : TETA = 0 < NU < 1                            
C            WET ELEMENT           : TETA = NU = 1.0
C            DRY ELEMENT           : TETA = NU = 0.0    
C           
C            LE CRITERE DE DECOUVREMENT EST CELUI DE J.-M. JANIN :      
C            LE FOND D'UN POINT DE L'ELEMENT EST PLUS HAUT QUE LA       
C            SURFACE LIBRE D'UN AUTRE.                                  
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      TETA      |<-- |  NU (PAR ELEMENT)                   
C |      HN,ZF     | -->|  HAUTEUR ET FOND     
C |      MESH      | -->|  STRUCTURE DE MAILLAGE                        
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C  APPELE PAR : PROPAG                                       
C                                                                       
C  SOUS-PROGRAMME APPELE : PORO11                              
C                                                                       
C********************************************************************** 
C
      USE BIEF
      USE INTERFACE_TELEMAC2D, EX_POROS => POROS
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(IN)    :: ZF,HN
      TYPE(BIEF_OBJ), INTENT(INOUT) :: TETA
      TYPE(BIEF_MESH), INTENT(IN)   :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                      
      INTEGER IELMZ,IELMH                              
C                                                                                                                                             
      IELMZ=ZF%ELM
      IELMH=HN%ELM
C                                                                       
C-----------------------------------------------------------------------
C
C     1) COMPUTATION OF POROSITY ON TIDAL FLATS
C                                                                                                                                    
      IF(IELMZ.EQ.11.AND.IELMH.EQ.11) THEN
C                                                                       
        CALL PORO11(TETA%R,ZF%R,                  
     *              HN%R,MESH%IKLE%I,MESH%NELEM,MESH%NELMAX) 
C                                                                                                                                             
      ELSE                                                              
C                                                                       
        IF(LNG.EQ.1) WRITE(LU,10) IELMH,IELMZ                           
        IF(LNG.EQ.2) WRITE(LU,11) IELMH,IELMZ                           
10      FORMAT(1X,'POROS : DISCRETISATION NON PREVUE :'    ,I6,' ',I6) 
11      FORMAT(1X,'POROS : DISCRETIZATION NOT IMPLEMENTED:',I6,' ',I6) 
        CALL PLANTE(1)                                                  
        STOP                                                            
C                                                                       
      ENDIF
C
C     2) CORRECTION BY A USER SUBROUTINE
C
      CALL CORPOR(TETA)                                                              
C                                                                        
C-----------------------------------------------------------------------
C                                                       
      RETURN                                                            
      END
