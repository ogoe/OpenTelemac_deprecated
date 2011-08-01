C                      ****************** 
                       SUBROUTINE RIDE_VR                                     
C                      ****************** 
C 
     * (KSR,KS,UNORM,HN,GRAV,XMVE,XMVS,NPOIN,ACLADM)
C 
C*********************************************************************** 
C  SISYPHE VERSION 6.0               C. VILLARET (LNHE)  AG DAVIES (UCW)
C                                                  
C  COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT                                                                           
C*********************************************************************** 
C                                                                                  
C  FONCTION:  CALCUL DES DIMENSIONS DES RIDES A L'EQUILIBRE 
C  VAN RIJN (2007) (courant seul) 
C ---------------------------------------------------------------------         
C                             ARGUMENTS                                          
C .________________.____._______________________________________________          
C |      NOM       |MODE|                   ROLE                                 
C |________________|____|_______________________________________________          
C |   KSR           |<-- | COEFFICIENT DE RUGOSITE DE PEAU
C |   KS           |<-- | COEFFICIENT DE RUGOSITE TOTALE
C |   HN           | -->| HAUTEUR D'EAU 
C |   UNORM        | -->| INTENSITE DU COURANT 
C |   ACLAD50           | -->| DIAMTRE MOYEN DU SEDIMENT 
C |   GRAV         | -->| ACCELERATION DE LA PESANTEUR                           
C |   XMVS         | -->| MASSE VOLUMIQUE DU SEDIMENT                   
C |   XMVE         | -->| MASSE VOLUMIQUE DE L'EAU 
C |   NPOIN        | -->| NOMBRE DE POINTS                                       
C |________________|____|_______________________________________________          
C  PROGRAMME APPELANT : SISYPHE                                                    
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C 
      IMPLICIT NONE                                                              
C 
      INTEGER I,NPOIN 
C
      DOUBLE PRECISION, INTENT(INOUT)  :: KSR(NPOIN),KS(NPOIN)
      DOUBLE PRECISION, INTENT(IN)     :: GRAV,XMVE,XMVS                 
      DOUBLE PRECISION, INTENT(IN)     :: HN(NPOIN)
      DOUBLE PRECISION, INTENT(IN)     :: ACLADM(NPOIN),UNORM(NPOIN) 
C            
C local variables
C
      DOUBLE PRECISION UC,AI,ZERO,KSCR,KSCD,KSCMR,MOB,FES,FFS
      DOUBLE PRECISION DSAND,DGRAVEL,DSILT
C
C--------------------------------------------------------------------- 
C 
      ZERO=1.D-6 
      DSILT=0.000032D0   
      DGRAVEL=0.002D0
      DSAND=0.000062D0
C
C CALCULATION OF CURRENT-DOMINATED ROUGHNESS USING VAN RIJN (2007)
C
      DO I=1,NPOIN
C
C Mobility number for current only
C
         AI  = ACLADM(I)*GRAV*(XMVS-XMVE)/XMVE             
         MOB = UNORM(I)**2/AI
C
C RIPPLE ROUGHNESS
C
         IF(ACLADM(I).LE.0.25D0*DGRAVEL)THEN        
           FES=1.D0
         ELSE
           FES=(0.25D0*DGRAVEL/ACLADM(I))**1.5D0             
         ENDIF 
C 
         IF(ACLADM(I).LT.DSILT)THEN 
           KSCR=20.D0*DSILT
         ELSE
           AI= TANH(0.015D0*(MOB-150.D0))
           KSCR=FES*ACLADM(I)*(85.D0-65.D0*AI)
         ENDIF
C
C MEGARIPPLE ROUGHNESS
C
         IF(ACLADM(I).GE.(1.5D0*DSAND))THEN
           FFS=1.D0
         ELSE
           FFS=ACLADM(I)/1.5D0/DSAND
         ENDIF
         IF(ACLADM(I).LE.DSILT)THEN
           KSCMR=0.D0            
         ELSE
           KSCMR=0.00002D0*FFS*HN(I)*(1.D0-EXP(-0.05D0*MOB))
           IF(MOB.GT.550.D0.AND.ACLADM(I).GE.1.5D0*DSAND)THEN
             KSCMR=0.02D0
           ELSEIF(MOB.GT.550D0.AND.ACLADM(I).LT.1.5D0*DSAND)THEN
             KSCMR=200.D0*ACLADM(I)
           ENDIF
         ENDIF
C
C DUNE ROUGHNESS
C
         IF(ACLADM(I).LT.DSILT) THEN
           KSCD=0.D0
         ELSE
           AI=(1.D0-EXP(-0.02D0*MOB))*(600.D0-MOB)
           KSCD=0.00008D0*FFS*HN(I)* AI
         ENDIF
         IF(MOB.GT.600.D0) KSCD=0.D0
         IF(KSCD.GT.1.D0) KSCD=1.D0
C
C ***RIPPLE BED ROUGHNESS FOR SEDIMENT COMPUTATIONS IN SISYPHE ***
C
         KSR(I)=KSCR
C 
C *** TOTAL ROUGHNESS FOR COMPUTATIONS IN TELEMAC2D **
C
         KS(I)=SQRT(KSCR**2+KSCMR**2+KSCD**2)
C
      ENDDO
C                          
CAGD****************************************************       
C
C--------------------------------------------------------------------- 
C        
      RETURN  
      END 
