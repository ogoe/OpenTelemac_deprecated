C                       *********************** 
                        SUBROUTINE TOBW_SISYPHE
C                       ***********************                       
C                                                                    
     *(TOBW ,CF, FW, UW,TW,HN,NPOIN,XMVE) 
C                                                                       
C***********************************************************************
C  SISYPHE VERSION 5.4                   C. VILLARET (LNHE)   01/10/2003
C                                                
C  COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT                                                                          
C***********************************************************************
C                                                                       
C      FONCTION: CALCUL DE LA CONTRAINTE de frottement GENEREE par la houle
C                Le coefficient de frottement est calculée par la formule de Swart (1976)
C                                                       
C---------------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.___________________________________________________
C |      NOM       |MODE|                   ROLE                            |
C |________________|____|________________________________________________   |
C |    TOBW        |<-- |  CONTRAINTE TOTALE AU FOND                        | 
C !    CF          | -->|  COEFFICIENT DE FROTTEMENT QUADRATIQUE (COURAN    |
C |    FW          |<-- |  COEFFICIENT DE FROTTEMENT quadratique (houle)    | 
C |    UW          | -->| VITESSE ORBITALE DE LA HOULE                      |           
C |    TW          | -->| PERIODEE DE LA HOULE                              | 
C |    HN          | -->|  HAUTEUR D'EAU AU TEMPS N                         | 
C |    NPOIN       | -->|  NOMBRE DE POINTS DU MAILLAGE 2D                  |
C |    XMVE        | -->|  MASSE VOLUMIQUE DE L'EAU                         | 
C |________________|____|_________________________________________________ _|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C  APPELE PAR : SISYPHE                                               
C                                                                       
C********************************************************************** 
C                                                                       
      IMPLICIT NONE 
C                                                    
      INTEGER LNG,LU                                              
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                       
      INTEGER, INTENT(IN) :: NPOIN                                          
C                           
      DOUBLE PRECISION, INTENT(IN)    :: CF(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: UW(NPOIN),TW(NPOIN),HN(NPOIN)                                         
      DOUBLE PRECISION, INTENT(IN)    :: XMVE
      DOUBLE PRECISION, INTENT(INOUT) :: TOBW(NPOIN),FW(NPOIN) 
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
      DOUBLE PRECISION KS,AUX 
      DOUBLE PRECISION PI,AW,KARMAN
      PARAMETER (PI=3.141592653589793D0)
      PARAMETER (KARMAN=0.4D0)
C
C-----------------------------------------------------------------------
C        
      DO  I=1,NPOIN
C       KS : coefficient de Nikuradse (frottement total)      
        AUX=1.D0+KARMAN*SQRT(2.D0/MAX(CF(I),1.D-10))
        KS=30.D0*MAX(HN(I),1.D-8)*EXP(-AUX)    
        AW= UW(I)*TW(I) / (2.D0*PI)
        IF(AW/KS.GT.1.59D0) THEN
          FW(I)=EXP( -6.D0 + 5.2D0 * (AW/KS)**(-0.19D0) )
        ELSE
          FW(I)=0.3D0
        ENDIF
        TOBW(I)=0.5D0 * XMVE * FW(I) * UW(I)*UW(I)
      ENDDO
C                                                                       
C-----------------------------------------------------------------------
C                                                                
      RETURN                                                            
      END SUBROUTINE TOBW_SISYPHE
