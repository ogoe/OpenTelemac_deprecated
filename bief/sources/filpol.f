C                       *****************                               
                        SUBROUTINE FILPOL                               
C                       *****************                               
C                                                                       
     *( F , C , XSOM , YSOM , NSOM , MESH )                     
C                                                                       
C***********************************************************************
C BIEF VERSION 5.1           08/12/94    J-M HERVOUET (LNH) 30 87 80 18 
C***********************************************************************
C                                                                       
C  FONCTION : LES POINTS DU VECTEUR F QUI SONT A L'INTERIEUR DE LA      
C             ZONE DEFINIE PAR LE POLYGONE DE SOMMETS XSOM ET YSOM      
C             SONT INITIALISES A LA VALEUR C.                           
C                                                                       
C             LE POLYGONE EST PARCOURU DANS LE SENS TRIGONOMETRIQUE.    
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     F          |<-->|  VECTEUR RESULTAT.                            
C |     C          | -->|  CONSTANTE IMPOSEE A L'INTERIEUR DU POLYGONE. 
C |     XSOM, YSOM | -->|  TABLEAUX DES COORDONNEES DES SOMMETS.        
C |     NSOM       | -->|  NOMBRE DE SOMMETS.                           
C |     MESH       | -->|  BLOC DES TABLEAUX D'ENTIERS DU MAILLAGE.     
C |     XMESH      | -->|  BLOC DES TABLEAUX DE REELS DU MAILLAGE.      
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C-----------------------------------------------------------------------
C                                                                       
C  PROGRAMMES APPELES :                                                 
C                                                                       
C********************************************************************** 
C                                                                       
      USE BIEF, EX_FILPOL => FILPOL
C
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NSOM
      DOUBLE PRECISION, INTENT(IN) :: C,XSOM(*),YSOM(*)
      TYPE(BIEF_OBJ), INTENT(INOUT) :: F
      TYPE(BIEF_MESH), INTENT(IN) :: MESH     
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                       
      INTEGER NPOIN,IELM                    
C                                                                                                                                             
C-----------------------------------------------------------------------
C                                                                                                                                              
C  TYPE D'ELEMENT DE LA VITESSE                                         
C                                                                                                             
      IELM = F%ELM
C                                                                       
C-----------------------------------------------------------------------
C                                                                                                                                              
      IF(IELM.EQ.10) THEN                                               
C                                                                       
        CALL FILP10( F%R , C , XSOM , YSOM , NSOM ,                       
     *       MESH%X%R,MESH%Y%R,MESH%NELEM,MESH%NELMAX,MESH%IKLE%I)
C                                                                       
      ELSEIF(IELM.EQ.11) THEN                                           
C                                                                                                                   
        NPOIN = F%DIM1                                               
        CALL FILP11(F%R,C,XSOM,YSOM,NSOM,MESH%X%R,MESH%Y%R,NPOIN)                    
C                                                                       
      ELSEIF(IELM.EQ.12) THEN                                           
C                                                                                                            
        NPOIN = F%DIM1 - MESH%NELEM                                       
        CALL FILP12( F%R , C , XSOM , YSOM , NSOM , MESH%X%R , MESH%Y%R,
     *               NPOIN , MESH%NELEM , MESH%NELMAX , MESH%IKLE%I)               
C                                                                       
      ELSE                                                              
C                                                                       
C       ELM NON PREVU : ERREUR                                         
C                                                                       
        IF (LNG.EQ.1) WRITE(LU,100) IELM                                
        IF (LNG.EQ.2) WRITE(LU,101) IELM                                
100     FORMAT(1X,'FILPOL (BIEF) : IELM = ',1I6,' ELEMENT NON PREVU')   
101     FORMAT(1X,'FILPOL (BIEF): IELM = ',1I6,' ELEMENT NOT AVAILABLE')
        CALL PLANTE(0)                                                  
        STOP                                                            
C                                                                       
      ENDIF                                                             
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END                                                               
 
 
