C                       *****************                                       
                        SUBROUTINE CALCUW                                       
C                       *****************
C                                      
     * ( UW, H, HW, TW, GRAV ,NPOIN)                                            
C                                                                               
C***********************************************************************        
C SISYPHE VERSION 5.1                                           20/05/96        
C                                                
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT                                                             
C***********************************************************************        
C     FONCTION  : CALCUL DE LA VITESSE ORBITALE DE HOULE                                                               
C                               
C-----------------------------------------------------------------------        
C                             ARGUMENTS                                         
C .________________.____.______________________________________________         
C |      NOM       |MODE|                   ROLE                                
C |________________|____|______________________________________________         
C |   UW           |<-- | VITESSE ORBITALE DE LA HOULE                       
C |   H            | -->| HAUTEUR D'EAU                                         
C |   HW           | -->| HAUTEUR DE HOULE                                      
C |   TW           | -->| PERIODE DE LA HOULE                            
C |   GRAV         | -->| GRAVITE                                               
C |   NPOIN        | -->| NOMBRE DE POINTS                                      
C |________________|____|______________________________________________         
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)         
C-----------------------------------------------------------------------        
C PROGRAMME APPELANT : BIJKER                                                    
C PROGRAMMES APPELES : /                                                        
C***********************************************************************
C        
      IMPLICIT NONE                                                             
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                               
      INTEGER, INTENT(IN) :: NPOIN                                                             
      DOUBLE PRECISION, INTENT(INOUT) :: UW(NPOIN)                         
      DOUBLE PRECISION, INTENT(IN) :: TW(NPOIN),H(NPOIN), HW(NPOIN) 
      DOUBLE PRECISION, INTENT(IN) :: GRAV                                        
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                               
      DOUBLE PRECISION   PI,DPI2                                                
      PARAMETER ( PI = 3.141592653589793D0 , DPI2 = (4.D0*PI*PI) )              
      DOUBLE PRECISION   POL, Y ,X                                              
      INTEGER I                                                                 
      INTRINSIC SQRT, SINH                                                    
C                                                                              
C  RESOLUTION DE Y=X*TH(X) AVEC Y=(2*PI/TW)**2*H/G ET X=(2*PI/L)*H              
C  PAR UNE FONCTION POLYNOMIALE (METHODE DE HUNT -ORDRE 9)                      
C                                                                               
      DO 10 I=1,NPOIN                                                           
       IF ( (TW(I) .GT. 0.D0).AND.(HW(I).GT.0.D0) ) THEN                        
         Y = DPI2 / GRAV * H(I) / (TW(I) * TW(I))                               
         POL = 1.D0 + Y * ( 0.66667D0 +                                         
     *                Y * ( 0.35550D0 +                                         
     *                Y * ( 0.16084D0 +                                         
     *                Y * ( 0.06320D0 +                                         
     *                Y * ( 0.02174D0 +                                         
     *                Y * ( 0.00654D0 +                                         
     *                Y * ( 0.00171D0 +                                         
     *                Y * ( 0.00039D0 +                                         
     *                Y * ( 0.00011D0 ) ))))))))                                
         X = SQRT( Y*Y + Y / POL )                                             
C                                                                               
         IF ( X .GT. 10.D0) THEN                                                
            UW(I) = 0.D0                                                        
         ELSE                                                                   
            UW(I) = PI / TW(I) * HW(I) / (SINH(X))                             
         ENDIF                                                                  
       ELSE                                                                     
         UW(I) = 0.D0                                                           
       ENDIF                                                                    
 10   CONTINUE                                                                  

      RETURN                                                                    
      END SUBROUTINE CALCUW
