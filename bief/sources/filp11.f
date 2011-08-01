C                       *****************                               
                        SUBROUTINE FILP11                               
C                       *****************                               
C                                                                       
     *( F , C , XSOM , YSOM , NSOM , X , Y , NPOIN )                    
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
C :   NPOIN        : -->: NOMBRE DE POINTS DU MAILLAGE                 :
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
      USE BIEF, EX_FILP11 => FILP11
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NSOM , NPOIN
      DOUBLE PRECISION, INTENT(INOUT) :: F(*)
      DOUBLE PRECISION, INTENT(IN) :: X(*) , Y(*)                           
      DOUBLE PRECISION, INTENT(IN) :: XSOM(NSOM) , YSOM(NSOM) , C     
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                         
      INTEGER I                                          
C                                                                                                                        
C-----------------------------------------------------------------------
C                                                                       
      DO 10 I = 1 , NPOIN                                               
C                                                                       
        IF(INPOLY(X(I),Y(I),XSOM,YSOM,NSOM)) F(I) = C                   
C                                                                       
10    CONTINUE                                                          
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END                                                               
 
