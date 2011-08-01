C                       *************************                       
                        SUBROUTINE CORSTR_SISYPHE
C                       *************************                       
C                                                                       
C                                                                       
C***********************************************************************
C  SISYPHE VERSION 5.1    12/11/97    C. LE NORMANT (LNH) 30 87 78 54 
C                                                
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT                                                                          
C***********************************************************************
C                                                                       
C      FONCTION: CORRECTION DU COEFFICIENT DE FROTTEMENT SUR LE FOND    
C                QUAND IL EST VARIABLE EN TEMPS.                        
C                                                                         
C      FUNCTION: CORRECTION OF BOTTOM FRICTION COEFFICIENT    
C                (IF TIME VARIATION).                        
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    CHESTR      |<-- |  BOTTOM FRICTION COEFFICIENT                 |
C |    X,Y         | -->|  2D COORDINATES                              |
C |    NPOIN       | -->|  NUMBER OF GRID POINTS                       |
C |    PRIVE       | -->|  ARRAY FOR USER                               
C |    ZF          | -->|  BOTTOM ELEVATION                            |
C |   QU , QV      |<-- |  FLOW RATE ALONG X AND Y                     |
C |    H           | -->|  WATER DEPTH                                 | 
C |    TIME        | -->|  TIME                                        | 
C |________________|____|______________________________________________|
C MODE : -->(INPUT), <--(RESULT), <-->(MODIFIED DATA) 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C  CALLED BY : SISYPHE                                               
C                                                                                                                                              
C********************************************************************** 
C
      USE BIEF
      USE DECLARATIONS_SISYPHE
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END SUBROUTINE CORSTR_SISYPHE 
