C                       *****************                               
                        SUBROUTINE FILP10                               
C                       *****************                               
C                                                                       
     *( F , C , XSOM , YSOM , NSOM , X , Y , NELEM , NELMAX , IKLE )    
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
C :   NSOM         : -->: NOMBRE DE SOMMETS                            :
C :   NELEM        : -->: NOMBRE D'ELEMENTS DU MAILLAGE.                
C :   NELMAX       : -->: NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE.        
C :   IKLE         : -->: TABLEAU DES CONNECTIVITES.                    
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
      USE BIEF, EX_FILP10 => FILP10
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NSOM , NELEM , NELMAX
      DOUBLE PRECISION, INTENT(INOUT) :: F(*)
      DOUBLE PRECISION, INTENT(IN) :: X(*) , Y(*)                           
      DOUBLE PRECISION, INTENT(IN) :: XSOM(NSOM) , YSOM(NSOM) , C
      INTEGER, INTENT(IN) :: IKLE(NELMAX,3)      
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                       
      INTEGER IELEM , I1 , I2 , I3                                                        
C                                                                                                                                            
C-----------------------------------------------------------------------
C                                                                       
      DO 10 IELEM = 1 , NELEM                                           
C                                                                       
        I1 = IKLE(IELEM,1)                                              
        I2 = IKLE(IELEM,2)                                              
        I3 = IKLE(IELEM,3)                                              
        IF( INPOLY(X(I1),Y(I1),XSOM,YSOM,NSOM) .AND.                    
     *      INPOLY(X(I2),Y(I2),XSOM,YSOM,NSOM) .AND.                    
     *      INPOLY(X(I3),Y(I3),XSOM,YSOM,NSOM) ) F(IELEM) = C           
C                                                                       
10    CONTINUE                                                          
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END                                                               
 
