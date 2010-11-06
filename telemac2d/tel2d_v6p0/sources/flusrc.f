C                       *****************                               
                        SUBROUTINE FLUSRC                               
C                       *****************                               
C                                                                       
     *(IEL1,IEL2,ISEGIN,VNOIN,W,FLUSCE,X,Y,AIRS,NPOIN,NSEG,ZF,EPS,G)    
C                                                                       
C***********************************************************************
C TELEMAC 2D VERSION 5.2   19/08/94             N.GOUTAL                     
C-----------------------------------------------------------------------
C                                                                       
C  FONCTION  : . CALCUL DES FLUX DUS AUX TERMES SOURCES                 
C                 - DE TYPE  DECENTRE                                   
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C ! . FLUX         !<-- ! TABLEAU DES FLUX DE TYPE ROE                 !
C ! . W            ! -->! VARIABLES CONSERVATIVES DU PB A L'INSTANT N  !
C ! . VNOIN        ! -->! NORMALE DU SEGMENT INTERNE                   !
C !                !    ! (2 PREMIERES COMPOSANTES) ET                 !
C !                !    ! LONGUEUR DE CE SEGMENT (3IEME COMPOSANTE)    !
C ! . AIRS         ! -->! AIRES DES CELLULES DU MAILLAGE.              !
C ! . WINF         ! -->! FLUX FRONTIERE ENTREE-SORTIE SI STAT.        !
C !                !    ! WINF INITIAL SI INSTAT.                      !
C ! . ZF           !    ! FOND.                                        !
C !________________!____!______________________________________________!
C !   /COMMON/     !    !                                              !
C ! . 3            ! -->! NOMBRE DE FONCTIONS INCONNUES DU PROBLEME    !
C ! . NDIM         ! -->! DIMENSION DE L'ESPACE.                       !
C ! . NSEG         ! -->! NOMBRE TOTAL DE SEGMENTS DU MAILLAGE         !
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)   
C        -- (TABLEAU DE TRAVAIL)                                        
C-----------------------------------------------------------------------
C     - SOUS PROGRAMME(S) APPELANT : RESOLU                             
C     - SOUS PROGRAMME(S) APPELES  : CLINST FLUSEW                      
C     - PORTABILITE:                                                    
C***********************************************************************
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NPOIN,NSEG,ISEGIN,IEL1,IEL2 
      DOUBLE PRECISION, INTENT(IN)    :: G,EPS,VNOIN(3,NSEG),ZF(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: AIRS(NPOIN),X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: W(3,NPOIN) 
      DOUBLE PRECISION, INTENT(INOUT) :: FLUSCE(3,NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                                                                  
      INTEGER INDIC(2)                                                  
C                                                                                                                             
      DOUBLE PRECISION XGI, YGI, XGJ, YGJ, DIJ, A1, A2                  
      DOUBLE PRECISION D1,HI,UI,VI,HJ,VJ,UJ,XN,YN,RNORM             
      DOUBLE PRECISION CT2,CT,RLAMB0,RLAMBM,ALPHA,CI2                   
      DOUBLE PRECISION CJ,CJ2,RLAMBJ,PROD,RLAMBI,RLAMBP                 
      DOUBLE PRECISION RI,RJ,UT,VT,CI,UN                                                                                              
      DOUBLE PRECISION T11(3),T21(3),T31(3),TS11(3),TS21(3),TS31(3)     
      DOUBLE PRECISION T12(3),T22(3),T32(3),TS12(3),TS22(3),TS32(3)     
      DOUBLE PRECISION GE(3),FE(3),ZF1,ZF2,PSA1,PSA2,PSA,UM,VM                           
C                                                                       
      INTRINSIC MIN,MAX                                                 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C------                                                                 
C 1. INITIALISATIONS                                                    
C------                                                                 
C                                                                       
      D1 = 0.3D0                                                        
C                                                                       
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                                                                       
C------                                                                 
C 1. CALCUL DES TERMES SOURCES A L'INTERFACE IEL1 , IEL2                
C------                                                                 
C                                                                       
C                                                                       
C                                                                       
C   --->    QUELQUES CALCULS INTERMEDIAIRES                             
C           ------------------------------                              
C                                                                       
       HI = W(1,IEL1)                                                   
       IF (HI.GT.EPS) THEN                                              
         UI = W(2,IEL1) / HI                                            
         VI = W(3,IEL1) / HI                                            
         INDIC(1) = 0                                                   
       ELSE                                                             
         UI = 0.D0                                                      
         VI = 0.D0                                                      
         INDIC(1) = 1                                                   
       ENDIF                                                            
C                                                                       
       HJ = W(1,IEL2)                                                   
       IF ( HJ.GT.EPS) THEN                                             
         UJ = W(2,IEL2) / HJ                                            
         VJ = W(3,IEL2) / HJ                                            
         INDIC(2) = 0                                                   
       ELSE                                                             
         UJ = 0.D0                                                      
         VJ = 0.D0                                                      
         INDIC(2) = 1                                                   
       ENDIF                                                            
C                                                                       
       XN = VNOIN (1,ISEGIN)                                            
       YN = VNOIN (2,ISEGIN)                                            
       RNORM = VNOIN (3,ISEGIN)                                         
C                                                                       
C    MODIFICATION DU FOND : AU REPOS AVEC UN ELEMENT SEC                
C                                                                       
       UN = UJ*XN+VJ*YN                                                 
       ZF1 = ZF(IEL1)                                                   
       IF(INDIC(1).EQ.1) THEN                                           
         IF((ZF1+EPS.GT.ZF(IEL2)+HJ).AND.(UN.GE.-EPS)) THEN             
           ZF1 = ZF(IEL2)+HJ-EPS                                        
         ENDIF                                                          
       ENDIF                                                            
C                                                                       
       UN = UI*XN+VI*YN                                                 
       ZF2 = ZF(IEL2)                                                   
       IF(INDIC(2).EQ.1) THEN                                           
         IF((ZF2+EPS.GT.ZF(IEL1)+HI).AND.(UN.LE.EPS)) THEN              
           ZF2 = ZF(IEL1)+HI-EPS                                        
         ENDIF                                                          
       ENDIF                                                            
C                                                                       
C   --->    CALCUL DES MOYENNES DE ROE DE U,V,H,C**2 ET C               
C           ---------------------------------------------               
C                                                                       
       IF(HI.LE.0.D0) THEN
         HI = 0.D0
C LIGNE SUIVANTE AJOUTEE PAR JMH
         RI = 0.D0 
       ELSE
         RI = SQRT ( HI )                                                 
       ENDIF 
       IF(HJ.LE.0.D0) THEN
         HJ = 0.D0
C LIGNE SUIVANTE AJOUTEE PAR JMH
         RJ = 0.D0 
       ELSE
         RJ = SQRT ( HJ )                                                 
       ENDIF
C      MAX DES DEUX LIGNES SUIVANTES AJOUTE PAR JMH                                                                       
       UT = ( RI * UI + RJ * UJ ) / MAX(RI+RJ,1.D-8)                            
       VT = ( RI * VI + RJ * VJ ) / MAX(RI+RJ,1.D-8)                            
       UM = (UI+UJ)/2.D0                                                
       VM = (VI+VJ)/2.D0                                                
       CT2 = G*(HI+HJ)/2.D0                                             
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
C     CALCUL DES MATRICES DES VALEURS PROPRES                           
C--------------------------------------------                           
         T11(1) = 1.D0                                                  
         T11(2) = UT - CT * XN                                          
         T11(3) = VT - CT * YN                                          
         T21(1) = 0.D0                                                  
         T21(2) = CT * YN                                               
         T21(3) = -CT * XN                                              
         T31(1) = 1.D0                                                  
         T31(2) = UT + CT * XN                                          
         T31(3) = VT + CT * YN                                          
C                                                                       
         T12(1) = 1.D0                                                  
         T12(2) = UT + CT * XN                                          
         T12(3) = VT + CT * YN                                          
         T22(1) = 0.D0                                                  
         T22(2) = -CT * YN                                              
         T22(3) = +CT * XN                                              
         T32(1) = 1.D0                                                  
         T32(2) = UT - CT * XN                                          
         T32(3) = VT - CT * YN                                          
C                                                                       
         TS11(1) = (UT * XN + VT * YN) * CT + CT2                       
         TS21(1) = (2.D0 * VT * XN - 2.D0 * UT * YN) * CT               
         TS31(1) = -(UT * XN + VT * YN) * CT + CT2                      
         TS11(2) = -XN * CT                                             
         TS21(2) = 2.D0 * YN * CT                                       
         TS31(2) = XN * CT                                              
         TS11(3) = -YN * CT                                             
         TS21(3) = -2.D0 * XN * CT                                      
         TS31(3) = YN * CT                                              
C                                                                       
         TS12(1) = -(UT * XN + VT * YN) * CT + CT2                      
         TS22(1) = -(2.D0 * VT * XN - 2.D0 * UT * YN) * CT              
         TS32(1) = +(UT * XN + VT * YN) * CT + CT2                      
         TS12(2) = +XN * CT                                             
         TS22(2) = -2.D0 * YN * CT                                      
         TS32(2) = -XN * CT                                             
         TS12(3) = +YN * CT                                             
         TS22(3) = +2.D0* XN * CT                                       
         TS32(3) = -YN * CT                                             
C                                                                       
C----------CALCULS POUR LES TERMES SOURCES--------------------          
C                                                                       
C                                                                       
          XGI = X(IEL1)                                                 
          YGI = Y(IEL1)                                                 
          XGJ = X(IEL2)                                                 
          YGJ = Y(IEL2)                                                 
C                                                                       
          DIJ = SQRT ((XGJ -XGI)**2 + (YGJ - YGI)**2)                   
          A1  = VNOIN(3,ISEGIN)*DIJ/2.D0                                
          A2  = VNOIN(3,ISEGIN)*DIJ/2.D0                                
C                                                                       
C  GRADIENTS DE FOND                                                    
C                                                                       
        GE(1)=0.D0                                                      
        GE(2)=G*((HI+HJ)/2.D0)*(ZF2-ZF1)*XN/DIJ                         
        GE(3)=G*((HI+HJ)/2.D0)*(ZF2-ZF1)*YN/DIJ                         
C                                                                       
C  TERMES DUS AU FROTTEMENT                                             
C                                                                       
C       CH = 900.D0                                                     
C       H = (CT2/G)**(1./3.)                                            
        FE(1)= 0.D0                                                     
C       FE(2)= G*UT*SQRT(UT**2 + VT**2)/((CH**2)*H)                     
C       FE(3)= G*VT*SQRT(UT**2 + VT**2)/((CH**2)*H)
        FE(2) = 0.D0
        FE(3) = 0.D0                     
C                                                                       
        GE(1)=GE(1)+FE(1)                                               
        GE(2)=GE(2)+FE(2)                                               
        GE(3)=GE(3)+FE(3)
C                                               
C---------------------------------------------------------------------  
         IF(RLAMB0.GE.-.000001D0) THEN                                  
C        ---- SEGMENT SORTIE---------                                   
C                                                                       
C   --->   PETITS CALCULS                                              
C                                                                       
           RLAMBM = RLAMB0 - CT                                        
C                                                                       
C                                                                       
           ALPHA = UI * XN + VI * YN                                   
C                                                                       
CTBTB DEBUT : MODIFICATION DE RLAMBM SI RLAMBM < D1                     
C                                                                       
           CI2 = G*HI                                                   
           CI = SQRT (CI2)                                              
           CJ2 =  G*HJ                                                  
           CJ = SQRT (CJ2)                                              
           RLAMBI = ALPHA - CI                                          
           RLAMBJ = UJ * XN + VJ * YN - CJ                                 
           PROD = RLAMBI * RLAMBJ                                       
C                                                                       
           IF ( RLAMBI .LT. 0.D0 .AND. RLAMBJ .GT. 0.D0                 
     *                                                   ) THEN         
            RLAMBM = MIN(0.D0,RLAMBM) - ABS(RLAMBI - RLAMBJ) / 4.D0     
           ENDIF                                                        
C                                                                       
C------------CALCUL DES TERMES SOURCES ------------------------         
C                                                                       
            FLUSCE (1,IEL1) = 0.D0                                      
            FLUSCE (2,IEL1) = 0.D0                                      
            FLUSCE (3,IEL1) = 0.D0                                      
C                                                                       
            FLUSCE (1,IEL2) = 0.D0                                      
            FLUSCE (2,IEL2) = 0.D0                                      
            FLUSCE (3,IEL2) = 0.D0                                      
C                                                                       
C                                                                       
C                                                                       
C   --->    TEST SUR LE SIGNE DE LAMBDAM                                
C           ----------------------------                                
C                                                                       
            IF ( RLAMBM . LT . 0.D0 ) THEN                              
C           - - - - - - - - - - - - - -                                 
C                                                                       
C----------CALCUL DES TERMES SOURCES --------------------------         
C                                                                       
              PSA = TS11(1)*GE(1)+TS11(2)*GE(2)+TS11(3)*GE(3)           
C                                                                       
         FLUSCE(1,IEL1) = PSA*T11(1)                                    
         FLUSCE(2,IEL1) = PSA*T11(2)                                    
         FLUSCE(3,IEL1) = PSA*T11(3)                                    
C                                                                       
C                                                                       
              PSA1= TS12(1)*GE(1)+TS12(2)*GE(2)+TS12(3)*GE(3)           
              PSA2= TS22(1)*GE(1)+TS22(2)*GE(2)+TS22(3)*GE(3)           
C                                                                       
C                                                                       
         FLUSCE(1,IEL2) = (PSA1*T12(1)+PSA2*T22(1))                     
         FLUSCE(2,IEL2) = (PSA1*T12(2)+PSA2*T22(2))                     
         FLUSCE(3,IEL2) = (PSA1*T12(3)+PSA2*T22(3))                     
C                                                                       
            ELSE                                                        
C           -----                                                       
C                                                                       
C                                                                       
         FLUSCE(1,IEL1) = 0.D0                                          
         FLUSCE(2,IEL1) = 0.D0                                          
         FLUSCE(3,IEL1) = 0.D0                                          
C                                                                       
         FLUSCE(1,IEL2) = GE(1)*CT2*2.D0                                
         FLUSCE(2,IEL2) = GE(2)*CT2*2.D0                                
         FLUSCE(3,IEL2) = GE(3)*CT2*2.D0                                
C                                                                       
            ENDIF                                                       
C          -----                                                        
C      TESTEST                                                          
         ELSE                                                           
C      TESTEST                                                          
C                                                                       
C   --->   PETITS CALCULS                                              
C          --------------                                              
C                                                                       
           RLAMBP = RLAMB0 + CT                                        
C
C          PROD = RLAMBI * RLAMBJ
C                                                                       
C                                                                       
           ALPHA = UI * XN + VI * YN                                   
C                                                                       
CTBTB DEBUT : MODIFICATION DE RLAMBP SI RLAMBM < D1                     
C                                                                       
           CI2 = G*HI                                                   
           CI = SQRT (CI2)                                              
           CJ2 =  G*HJ                                                  
           CJ = SQRT (CJ2)                                              
           RLAMBI = ALPHA - CI                                          
           RLAMBJ = UJ * XN + VJ * YN - CJ                                
           PROD = RLAMBI * RLAMBJ                                                                          
C                                                                       
CTBTB : MODIF UNIQUEMENT DANS LA DETENTE :                              
            IF ( RLAMBI .LT. 0.D0 .AND. RLAMBJ .GT. 0.D0) THEN              
C                                                                       
CTBTB : MODIF DANS LA DETENTE OU DANS LE CHOC :                         
C           IF ( PROD . LT . 0.D0 .AND. ABS(RLAMBP).LT.D1 ) THEN          
C                                                                       
               RLAMBP = MAX(0.D0,RLAMBP) + ABS(RLAMBI - RLAMBJ)/4.D0    
            ENDIF                                                       
CTBTB FIN                                                               
C                                                                       
C                                                                       
C-----------CALCUL DES TERMES SOURCE --------------------------         
C                                                                       
         FLUSCE(1,IEL1) = GE(1)*CT2*2.D0                                
         FLUSCE(2,IEL1) = GE(2)*CT2*2.D0                                
         FLUSCE(3,IEL1) = GE(3)*CT2*2.D0                                
C                                                                       
         FLUSCE(1,IEL2) = 0.D0                                          
         FLUSCE(2,IEL2) = 0.D0                                          
         FLUSCE(3,IEL2) = 0.D0                                          
C                                                                                                                                               
C                                                                       
C   --->    TEST SUR LE SIGNE DE LAMBDAP                                
C           ----------------------------                                
C                                                                       
            IF ( RLAMBP . GT . 0.D0 ) THEN                              
C           - - - - - - - - - - - - - -                                 
C                                                                       
C-----------CALCUL DES TERMES SOURCE ------------------                 
C                                                                       
C                                                                       
         FLUSCE(1,IEL1) = 0.D0                                          
         FLUSCE(2,IEL1) = 0.D0                                          
         FLUSCE(3,IEL1) = 0.D0                                          
C                                                                       
         FLUSCE(1,IEL2) = 0.D0                                          
         FLUSCE(2,IEL2) = 0.D0                                          
         FLUSCE(3,IEL2) = 0.D0                                          
              PSA1= TS11(1)*GE(1)+TS11(2)*GE(2)+TS11(3)*GE(3)           
              PSA2= TS21(1)*GE(1)+TS21(2)*GE(2)+TS21(3)*GE(3)           
C                                                                       
C                                                                       
         FLUSCE(1,IEL1) = (PSA1*T11(1)+PSA2*T21(1))                     
         FLUSCE(2,IEL1) = (PSA1*T11(2)+PSA2*T21(2))                     
         FLUSCE(3,IEL1) = (PSA1*T11(3)+PSA2*T21(3))                     
C                                                                       
                                                                        
              PSA = TS12(1)*GE(1)+TS12(2)*GE(2)+TS12(3)*GE(3)           
C                                                                       
         FLUSCE(1,IEL2) = PSA*T12(1)                                    
         FLUSCE(2,IEL2) = PSA*T12(2)                                    
         FLUSCE(3,IEL2) = PSA*T12(3)                                    
C                                                                       
            ENDIF                                                       
C           -----                                                       
C                                                                       
C       TESTEST                                                         
         ENDIF                                                          
C       TESTEST                                                         
      FLUSCE(1,IEL1)=FLUSCE(1,IEL1)*A1/CT2                              
      FLUSCE(2,IEL1)=FLUSCE(2,IEL1)*A1/CT2                              
      FLUSCE(3,IEL1)=FLUSCE(3,IEL1)*A1/CT2                              
      FLUSCE(1,IEL2)=FLUSCE(1,IEL2)*A2/CT2                              
      FLUSCE(2,IEL2)=FLUSCE(2,IEL2)*A2/CT2                              
      FLUSCE(3,IEL2)=FLUSCE(3,IEL2)*A2/CT2                              
C                                                                       
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                                                                       
C                                                                       
       RETURN                                                           
       END 
