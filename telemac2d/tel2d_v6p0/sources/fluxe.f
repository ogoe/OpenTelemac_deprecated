C                       ****************                               
                        SUBROUTINE FLUXE                                 
C                       ****************                               
C                                                                       
     *(HJ,UJ,VJ,HI,UI,VI,XN,YN,RNORM,A1,A2,G,FLULOC)                   
C                                                                       
C***********************************************************************
C TELEMAC 2D VERSION 5.2                  N. GOUTAL 24/11/97
C-----------------------------------------------------------------------
C                                                                       
C  FONCTION  : . INTEGRATION EN TEMPS                                   
C 
C  NOTE JMH : RNORM, A ET A2 NON UTILISES
C                                                                      
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C ! . AIRS         ! -->! TABLEAU DES AIRES DES CELLULES               !
C ! . H            !<-->! ENTHALPIE                                    !
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)   
C        -- (TABLEAU DE TRAVAIL)                                        
C-----------------------------------------------------------------------
C                                                                      
         IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
         DOUBLE PRECISION, INTENT(INOUT) :: FLULOC(3)
         DOUBLE PRECISION, INTENT(IN) :: G,HI,HJ,UI,UJ,VI,VJ,RNORM
         DOUBLE PRECISION, INTENT(IN) :: XN,YN,A1,A2
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
         DOUBLE PRECISION RI,RJ,CT2,UT,VT,RLAMB0
         DOUBLE PRECISION CT,PRII,PRIJ,ALPHA
         DOUBLE PRECISION RLAMBM,PS,SA,RLAMBP,PROD,TW(3)
         DOUBLE PRECISION TR(3),T(3),CI2,CI,CJ,CJ2,RLAMBI,RLAMBJ
C
C-----------------------------------------------------------------------
C                                                                       
C   --->    CALCUL DES MOYENNES DE ROE DE U,V,H,C**2 ET C               
C           --------------------------------------------- 
C              
         IF(HI.LE.0.D0) THEN
           RI = 0.D0
         ELSE
           RI = SQRT ( HI )                                               
         ENDIF
         IF (HJ.LE.0.D0) THEN
           RJ = 0.D0
         ELSE
           RJ = SQRT ( HJ )                                               
         ENDIF
C                                                                       
         UT = ( RI * UI + RJ * UJ ) /(RI + RJ)                          
         VT = ( RI * VI + RJ * VJ ) /(RI + RJ)                          
         IF ( (HI + HJ) .LE.0.D0)  THEN 
           CT2= 0.D0 
         ELSE
         CT2 = G*(HI+HJ)/2.D0                                             
         ENDIF
         CT = SQRT ( CT2 )                                                
C                                                                       
C   --->  TEST SUR LE SIGNE DE LA VALEUR PROPRE LAMB0 = <UT,N>          
C           ----------------------------------------------------------  
C                                                                       
         RLAMB0 = UT * XN + VT * YN                                     
C                                                                       
CTBTB DEBUT : MODIFICATION DE RLAMB0 SI RLAMB0 < D1                     
CC  A FAIRE SI ON VEUT CAR IL FAUT AJOUTER LES FLUX POUR LA V.P. DOUBLE 
CC                                                                      
CC                                                                      
CTBTB FIN                                                               
C                                                                       
C---------------------------------------------------------------------  
         IF  ( RLAMB0 . GE .-0.000001D0 ) THEN                          
C        ---- SEGMENT SORTIE---------                                   
C                                                                       
C   --->    PETITS CALCULS                                              
C                                                                       
            RLAMBM = RLAMB0 - CT                                        
C                                                                       
            PRII = G*(HI**2)/2.D0                                       
            PRIJ = G*(HJ**2)/2.D0                                       
            ALPHA = UI * XN + VI * YN                                   
C                                                                       
CTBTB DEBUT : MODIFICATION DE RLAMBM SI RLAMBM < D1                     
C                                                                       
           IF (HI.LE.0.D0) THEN
           CI2 = 0.D0
           PRII = 0.D0
           ELSE
           CI2 =  2.D0*PRII / HI                                        
           ENDIF
           IF (HJ.LE.0.D0) THEN
           CJ2 = 0.D0
           PRIJ = 0.D0
           ELSE
           CJ2 =  2.D0*PRIJ / HJ                                        
           ENDIF
           CI = SQRT (CI2)                                              
           CJ = SQRT (CJ2)                                              
           RLAMBI = ALPHA - CI                                          
           RLAMBJ = UJ * XN + VJ * YN - CJ                              
           PROD = RLAMBI * RLAMBJ                                       
C                                                                       
CTBTB : MODIF UNIQUEMENT DANS LA DETENTE :                              
           IF ( RLAMBI .LT. 0.D0 .AND. RLAMBJ .GT. 0.D0                 
CC   *                         .AND. ABS(RLAMBM) .LT. D1                
     *                                                   ) THEN         
C                                                                       
C     : MODIF DANS LA DETENTE OU DANS LE CHOC :                         
C          IF ( PROD . LT . 0.D0 .AND. ABS(RLAMBM).LT.D1 ) THEN         
C                                                                       
            RLAMBM = MIN(0.D0,RLAMBM) - ABS(RLAMBI - RLAMBJ) / 4.D0     
           ENDIF                                                        
C     FIN                                                               
C                                                                       
C   --->    CALCUL DU FLUX 1                                            
C                                                                       
            FLULOC(1) = ALPHA * HI                               
            FLULOC(2) = ALPHA * HI*UI                                
            FLULOC(3) = ALPHA * HI*VI                               
C                                                                       
            FLULOC (2) = FLULOC(2) + PRII * XN                          
            FLULOC (3) = FLULOC(3) + PRII * YN                          
C                                                                       
C   --->    TEST SUR LE SIGNE DE LAMBDAM                                
C           ----------------------------                                
C                                                                       
            IF ( RLAMBM . LT . 0.D0 ) THEN                              
C           - - - - - - - - - - - - - -                                 
C                                                                       
               T (1) = 1.D0                                             
               T (2) = UT - CT * XN                                     
               T (3) = VT - CT * YN                                     
C                                                                       
               TR(1) = HJ-HI                            
               TR(2) = HJ*UJ-HI*UI                            
               TR(3) = HJ*VJ-HI*VI                            
C                                                                       
               TW(1) = (UT*XN + VT*YN)*CT + CT2                         
               TW(2) = -XN*CT                                           
               TW(3) = -YN*CT                                           
C                                                                       
                 PS = TR(1)*TW(1)+TR(2)*TW(2)+TR(3)*TW(3)               
C                                                                       
C   --->    CALCUL DU FLUX LOCAL TOTAL                                  
C           --------------------------                                  
C                                                                       
               SA = PS * RLAMBM / (2.D0 * CT2 )                         
        FLULOC(1)= FLULOC(1)+SA*T(1)                                    
        FLULOC(2)= FLULOC(2)+SA*T(2)                                    
        FLULOC(3)= FLULOC(3)+SA*T(3)                                    
C                                                                       
C                                                                       
            ENDIF                                                       
C           -----                                                       
C                                                                       
C      TESTEST                                                          
         ELSE                                                           
C      TESTEST                                                          
C                                                                       
C   --->    PETITS CALCULS                                              
C           --------------                                              
C                                                                       
            RLAMBP = RLAMB0 + CT                                        
C                                                                       
C                                                                       
            ALPHA = UJ * XN + VJ* YN                                    
C                                                                       
CTBTB DEBUT : MODIFICATION DE RLAMBP SI RLAMBP < D1                     
C                                                                       
            IF (HI.LE.0.D0) THEN
              CI2 = 0.D0
            ELSE
              CI2 = G*HI                                                  
            ENDIF
            CI = SQRT (CI2)                                             
            IF (HJ.LE.0.D0) THEN
            CJ2 = 0.D0
            PRIJ = 0.D0                                      
            ELSE
            CJ2 = G*HJ                                                  
            PRIJ = G*(HJ**2)/2.D0                                       
            ENDIF
            CJ = SQRT (CJ2)                                             
            RLAMBI = UI * XN + VI * YN + CI                             
            RLAMBJ = ALPHA + CJ                                         
            PROD = RLAMBI * RLAMBJ                                      
C                                                                       
CTBTB : MODIF UNIQUEMENT DANS LA DETENTE :                              
            IF ( RLAMBI .LT. 0. .AND. RLAMBJ .GT. 0.                    
CC   *                          .AND. ABS(RLAMBP) .LT. D1               
     *                                                    ) THEN        
C                                                                       
CTBTB : MODIF DANS LA DETENTE OU DANS LE CHOC :                         
C           IF ( PROD . LT . 0. .AND. ABS(RLAMBP).LT.D1 ) THEN          
C                                                                       
               RLAMBP = MAX(0.D0,RLAMBP) + ABS(RLAMBI - RLAMBJ) / 4.    
            ENDIF                                                       
CTBTB FIN                                                               
C                                                                       
C   --->    CALCUL DU FLUX 1                                            
C           ----------------                                            
C                                                                       
            FLULOC(1) = ALPHA * HJ                               
            FLULOC(2) = ALPHA * HJ*UJ                              
            FLULOC(3) = ALPHA * HJ*VJ                              
C                                                                       
            FLULOC (2) = FLULOC(2) + PRIJ * XN                          
            FLULOC (3) = FLULOC(3) + PRIJ * YN                          
C                                                                       
C   --->    TEST SUR LE SIGNE DE LAMBDAP                                
C           ----------------------------                                
C                                                                       
            IF ( RLAMBP . GT . 0.D0 ) THEN                              
C           - - - - - - - - - - - - - -                                 
C                                                                       
               T(1) = 1.D0                                              
               T(2) = UT + CT * XN                                      
               T(3) = VT + CT * YN                                      
C                                                                       
               TR(1) = HJ-HI                            
               TR(2) = HJ*UJ-HI*UI                            
               TR(3) = HJ*VJ-HI*VI                            
C                                                                       
               TW(1) = (-UT*XN - VT*YN)*CT +CT2                         
               TW(2) = CT*XN                                            
               TW(3) = CT*YN                                            
C                                                                       
               PS = TR(1)*TW(1)+TR(2)*TW(2)+TR(3)*TW(3)                 
C                                                                       
C   --->    CALCUL DU FLUX LOCAL TOTAL                                  
C           --------------------------                                  
C                                                                       
               SA = - PS * RLAMBP / (2.D0 * CT2 )                       
        FLULOC(1)= FLULOC(1)+SA*T(1)                                    
        FLULOC(2)= FLULOC(2)+SA*T(2)                                    
        FLULOC(3)= FLULOC(3)+SA*T(3)                                    
C                                                                       
            ENDIF                                                       
C           -----                                                       
        ENDIF
C                                                                       
C                                                                       
        RETURN
        END 
