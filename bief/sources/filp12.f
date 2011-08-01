C                       *****************                               
                        SUBROUTINE FILP12                               
C                       *****************                               
C                                                                       
     *(F,C,XSOM,YSOM,NSOM,X,Y,NPOIN,NELEM,NELMAX,IKLE)                  
C                                                                       
C***********************************************************************
C  BIEF VERSION 5.1    06/12/94              C MOULIN (LNH) 30 87 83 81
C                                                                       
C***********************************************************************
C                                                                       
C    FONCTION  : INITIALISE UNE FONCTION A UNE VALEUR CONSTANTE         
C                A L'INTERIEUR D'UN POLYGONE.                           
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C :      NOM       :MODE:                   ROLE                       :
C :________________:____:______________________________________________:
C :   F            :<-->: FONCTION A INITIALISER                       :
C :   C            : -->: VALEUR CONSTANTE DONNEE A F DANS LE POLYGONE :
C :   XSOM         : -->: TABLEAU DES ABSCISSES DES SOMMETS DU POLYGONE:
C :   YSOM         : -->: TABLEAU DES ORDONNEES DES SOMMETS DU POLYGONE:
C :   X,Y          : -->: COORDONNEES DES POINTS DU MAILLAGE           :
C :   NSOM         : -->: NOMBRE DE SOMMETS                             
C :   NPOIN        : -->: NOMBRE DE POINTS DU MAILLAGE                  
C :   NELEM        : -->: NOMBRE D'ELEMENTS DU MAILLAGE                 
C :   NELMAX       : -->: NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE         
C :   IKLE         : -->: TABLE DE CONNECTIVITE.                        
C :________________:____:______________________________________________:
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C  APPELE PAR :                                                         
C                                                                       
C  SOUS-PROGRAMME APPELE :                                              
C                                                                       
C***********************************************************************
C
      USE BIEF, EX_FILP12 => FILP12
C                                                                      
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NSOM , NELEM , NELMAX , NPOIN
      DOUBLE PRECISION, INTENT(INOUT) :: F(*)
      DOUBLE PRECISION, INTENT(IN) :: X(*) , Y(*)                           
      DOUBLE PRECISION, INTENT(IN) :: XSOM(NSOM) , YSOM(NSOM) , C
      INTEGER, INTENT(IN) :: IKLE(NELMAX,4)      
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                        
      INTEGER I , I1 , I2 , I3 , IELEM                                             
C                                                                                             
      DOUBLE PRECISION XX , YY                  
C                                                                                                                                             
C-----------------------------------------------------------------------
C                                                                       
      DO 10 I = 1 , NPOIN                                               
C                                                                       
        IF(INPOLY(X(I),Y(I),XSOM,YSOM,NSOM)) F(I) = C                   
C                                                                       
10    CONTINUE                                                          
C                                                                       
      DO 20 IELEM = 1 , NELEM                                           
C                                                                       
        I1 = IKLE(IELEM,1)                                              
        I2 = IKLE(IELEM,2)                                              
        I3 = IKLE(IELEM,3)                                              
        XX = 0.3333333333D0 * ( X(I1)+X(I2)+X(I3) )                     
        YY = 0.3333333333D0 * ( Y(I1)+Y(I2)+Y(I3) )                     
        IF(INPOLY(XX,YY,XSOM,YSOM,NSOM)) F(IELEM+NPOIN) = C             
C                                                                       
20    CONTINUE                                                          
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END                                                               
 
