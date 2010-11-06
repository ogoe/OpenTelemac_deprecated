C                       *****************                               
                        SUBROUTINE FWSPEC                           
C                       *****************                               
C                                                                       
     *(FW,FWCOEF,X,Y,NPOIN,PRIVE,ZF)                      
C
C***********************************************************************
C
C  ARTEMIS VERSION 5.1    02/06/99   D. AELBRECHT (LNH) 01 30 87 74 12 
C
C  LINKED TO BIEF VERS. 5.0          J-M HERVOUET (LNH) 01 30 87 80 18
C                                                                       
C***********************************************************************
C                                                                       
C      FONCTION: SPECIFICATION DU COEFFICIENT DE FROTTEMENT SUR LE FOND 
C                SI IL  EST VARIABLE EN ESPACE.                         
C                                                                       
C      CE SOUS-PROGRAMME EST SIMPLEMENT UN MODELE                       
C      IL DOIT ETRE REMPLI PAR L'UTILISATEUR                            
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    FW          |<-- |  COEFFICIENT DE FROTTEMENT                   |
C |    FWCOEF      | -->|  COEFFICIENT DE FROTTEMENT CONSTANT IMPOSE   |
C |    X,Y         | -->|  COORDONNEE DU MAILLAGE .                    |
C |    NPOIN       | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C |    PRIVE       | -->|  TABLEAU DE TRAVAIL DEFINI DANS PRINCI       |
C |    ZF          | -->|  COTE DU FOND                                |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C  APPELE PAR : BERKHO                                                  
C                                                                       
C  SOUS-PROGRAMME APPELE : OS                                       
C                                                                       
C********************************************************************** 
C                                                                       
      USE BIEF
      USE INTERFACE_ARTEMIS, EX_FWSPEC => FWSPEC 
C
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C                                                                       
      INTEGER NPOIN                                  
C                                                                       
      DOUBLE PRECISION FW(NPOIN),X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION ZF(NPOIN),FWCOEF   
C
      TYPE(BIEF_OBJ) :: PRIVE
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C  ICI ON MET UN COEFFICIENT DE FROTTEMENT CONSTANT (EXEMPLE)           
C                                                                       
      CALL OV( 'X=C     ' , FW , X , X , FWCOEF , NPOIN )
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END                                                               
