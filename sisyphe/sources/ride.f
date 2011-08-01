C                       ***************                                    
                        SUBROUTINE RIDE                                       
C                       *************** 
C  
     * (KS,TW,UW,UNORM,GRAV,XMVE,XMVS,VCE,NPOIN,KSPRATIO,ACLADM) 
C                                                                        
C*********************************************************************** 
C  SISYPHE VERSION 6.0                   C. VILLARET (LNHE)   01/10/2003 
C                                                 
C  COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT                                                                           
C*********************************************************************** 
C                                                                                  
C  FONCTION:  CALCUL DES DIMENSIONS DES RIDES A L'EQUILIBRE 
C  METHODE DE WIBERG ET HARRIS, JGR 1994  
C ---------------------------------------------------------------------         
C                             ARGUMENTS                                          
C .________________.____._______________________________________________          
C |      NOM       |MODE|                   ROLE                                 
C |________________|____|_______________________________________________          
C |   KS           |<-- | COEFFICIENT DE RUGOSITE 
C |   HN           | -->| HAUTEUR D'EAU 
C |   Q            | -->| DEBIT MOYEN                                          
C |   UW           | -->| COURANT ORBITAL                                      
C |   TW           | -->| PERIODE DE HOULE 
C |   DM           | -->| DIAMTRE MOYEN DU SEDIMENT 
C |   GRAV         | -->| ACCELERATION DE LA PESANTEUR                           
C |   XMVS         | -->| MASSE VOLUMIQUE DU SEDIMENT                   
C |   XMVE         | -->| MASSE VOLUMIQUE DE L'EAU 
C |   VCE          | -->| VISCOSITE DE L'EAU  
C |   NPOIN        | -->| NOMBRE DE POINTS                                       
C |________________|____|_______________________________________________          
C  PROGRAMME APPELANT : SISYPHE                                                    
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C 
      IMPLICIT NONE                                                              
C 
      INTEGER I,NPOIN 
      DOUBLE PRECISION KS(NPOIN) 
C     SI GRANULO ETENDUE DM(NPOIN) 
      DOUBLE PRECISION  GRAV,XMVE,XMVS, VCE                     
      DOUBLE PRECISION UNORM(NPOIN), UW(NPOIN), TW(NPOIN) 
C                                      
      DOUBLE PRECISION PI, ZERO,AI 
C 
      DOUBLE PRECISION ETA, LAMBDA 
C 
      DOUBLE PRECISION AA,BB,CC,DD 
      DOUBLE PRECISION ALPHA,S,M,A0 
      DOUBLE PRECISION WH1,WH2,WH3 
      DOUBLE PRECISION VAR1,DHRA,LRA,HRA,LRO  
      DOUBLE PRECISION UC,KSP  
      DOUBLE PRECISION, INTENT(IN) :: KSPRATIO 
      DOUBLE PRECISION, INTENT(IN) :: ACLADM(NPOIN) 
C 
C--------------------------------------------------------------------- 
C 
      PI=4.D0*DATAN(1.D0) 
      ZERO=1.D-6 
C        
C     COEFFICIENTS 
C 
      WH1=0.095D0 
      WH2=0.442D0 
      WH3=2.28D0 
      AA=(WH2+1.D0)/2.D0/WH1 
      BB=AA*AA-WH3/WH1                 
      CC=1.D0/WH1 
C 
C     BOUCLE SUR LES POINTS 
C    
      DO I=1,NPOIN 
C         
C       FROTTEMENT DE PEAU 
C 
        KSP = KSPRATIO * ACLADM(I) 
        AI  = ACLADM(I)*GRAV*(XMVS-XMVE)/XMVE  
C 
C       MOBILITY NUMBER 
C 
        M=UW(I)**2/AI 
C 
        IF(M.LT.1.69D0) THEN 
C 
          KS(I)=KSP 
C 
        ELSE        
C 
C         WIBERG AND HARRIS 
C 
          A0=UW(I)*TW(I)/(2.D0*PI) 
          S=ACLADM(I)*SQRT(AI)/4.D0/VCE 
          LRA=535.D0*ACLADM(I) 
CJMB*************************************************** 
CJMB Line of code moved so alpha calculated before VAR1 
CJMB  TANNAKA AND DANG (1996) 
          UC=UNORM(I)
          IF(UW(I).GT.ZERO) THEN 
            ALPHA=(TANH(0.3D0*S**(2.D0/3.D0)))**2.5D0 
            ALPHA=1.D0+0.81D0*ALPHA*(UC/UW(I))**1.9D0 
          ELSE 
            ALPHA=1.D0 
          ENDIF 
CJMB******************************************************* 
 
          VAR1=LOG(ALPHA*2.D0*A0/LRA)           
          DD=MAX((BB-CC*VAR1),0.D0)          
          DHRA=EXP(AA-SQRT(DD)) 
          HRA=ALPHA*2.D0*A0/DHRA 
 
C 
          IF(DHRA.LE.20.D0) THEN 
C           ORBITAL RIPPLES DHRA<20 
            LRO=0.62D0*2.D0*A0*ALPHA 
            LAMBDA=LRO 
            ETA=0.17D0*LAMBDA 
          ELSEIF(DHRA.LE.100.D0) THEN  
C           SUB ORBITAL RIPPLES 20<DHRA<100         
            LRO=0.62D0*2.D0*A0*ALPHA 
            VAR1=(LOG(DHRA)-LOG(100.D0))/(LOG(20.D0)-LOG(100.D0)) 
            VAR1=LOG(LRA)+VAR1*(LOG(LRO)-LOG(LRA)) 
            LAMBDA=EXP(VAR1) 
            VAR1=LOG(ALPHA*2.D0*A0/LAMBDA) 
CCV 25/05               ETA=ALPHA*2.D0*A0/EXP(AA-SQRT(BB-CC*VAR1)) 
            DD=MAX((BB-CC*VAR1),0.D0) 
            ETA=ALPHA*2.D0*A0/EXP(AA-SQRT(DD)) 
          ELSE 
C           ANORBITAL RIPPLES DHRA>100 
C           LAMBDA NOT USED HERE BUT KEPT FOR OTHER FORMULAS 
C           LAMBDA=LRA 
            ETA=HRA 
          ENDIF 
C 
          KS(I)=MAX(ETA,KSP)    
C 
        ENDIF 
C       
      ENDDO 
C 
C--------------------------------------------------------------------- 
C        
      RETURN  
      END 
