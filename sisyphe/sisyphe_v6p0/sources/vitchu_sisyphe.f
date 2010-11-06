C                       *************************                               
                        SUBROUTINE VITCHU_SISYPHE
C                       *************************                               
C
     * ( WS , DENS , DM , GRAV , VCE )                  
C                                                                               
C***********************************************************************        
C SISYPHE VERSION 5.1                                           20/05/96 
C                                                
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT          
C***********************************************************************        
C                                                                               
C     FONCTION  : CALCULE LA VITESSE DE CHUTE                              
C-----------------------------------------------------------------------        
C                             ARGUMENTS                                         
C .________________.____.______________________________________________         
C |      NOM       |MODE|                   ROLE                                
C |________________|____|______________________________________________                                                                                          
C |   WS           | -->| VITESSE DE CHUTE DES PARTICULES
C |   DENS         | -->| POIDS DEJAUGE                              
C |   DM           | -->| DIAMETRE MOYEN DU SEDIMENT                            
C |   GRAV         | -->| ACCELERATION DE LA PESANTEUR                          
C |   VCE          | -->| VISCOSITE DE L'EAU                                    
C |________________|____|______________________________________________         
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)         
C-----------------------------------------------------------------------        
C PROGRAMME APPELANT : BIJKER                                                   
C PROGRAMMES APPELES :                                        
C***********************************************************************
C        
      IMPLICIT NONE                                                             
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                                
      DOUBLE PRECISION, INTENT(IN)    :: DENS,  DM,  GRAV, VCE
      DOUBLE PRECISION, INTENT(INOUT) :: WS 
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                       
C                                                                                                                                                      
C VITESSE DE CHUTE                                                              
C ================
C
      IF (DM.LT.1.D-4) THEN
        WS = DENS * DM * DM * GRAV / ( 18.D0 * VCE )                            
      ELSEIF (DM.LT.1D-3) THEN                                                 
        WS = 10.D0 * VCE / DM * (SQRT( 1.D0 + 0.01D0* DENS * GRAV *             
     *       DM**3.D0 / (VCE*VCE) ) -1.D0 )                                     
      ELSE                                                                      
        WS = 1.1D0 * SQRT( DENS * GRAV * DM )                                   
      ENDIF                                                                     
C
C-----------------------------------------------------------------------
C 
      RETURN
      END SUBROUTINE VITCHU_SISYPHE
