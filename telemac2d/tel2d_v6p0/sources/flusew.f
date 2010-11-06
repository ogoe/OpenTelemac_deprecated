C                       *****************                               
                        SUBROUTINE FLUSEW                               
C                       *****************                               
C                                                                       
     *(AMINF,UBOR,VBOR,NPOIN,EPS,G,W,                                          
     * XNEBOR,YNEBOR,NPTFR,LIMPRO,NBOR,KDIR,KNEU,KDDL)                 
C                                                                       
C***********************************************************************
C TELEMAC 2D VERSION 5.2     09/09/94            N. GOUTAL                        
C-----------------------------------------------------------------------
C 
C  NOTE JMH : UBOR ET VBOR PAS UTILISES                                                                      
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C ! . VALBOR       ! -->! VALBOR I.E. (RHO,U,V,P) INITIAL A L'INFINI   !
C ! . AMINF        !<-- ! A-(WINF,NINF).WINF AVEC  /NINF/ = 1          !
C !                !    !                                              !
C ! . NDIM         ! -->! DIMENSION DE L'ESPACE.                       !
C ! . 3            ! -->! NOMBRE DE FONCTIONS DU PROBLEME              !
C !  (3)           !    !                                              !
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)   
C        -- (TABLEAU DE TRAVAIL)                                        
C-----------------------------------------------------------------------
C     - SOUS PROGRAMME(S) APPELANT : INIGEO FLUROE                      
C     - SOUS PROGRAMME(S) APPELE   : AUCUN                              
C     - PORTABILITE: CRAY                                               
C***********************************************************************
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NPOIN,NPTFR,KDIR,KNEU,KDDL
      INTEGER, INTENT(IN)             :: LIMPRO(NPTFR,6),NBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: XNEBOR(NPTFR),YNEBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: W(3,NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: UBOR(NPTFR),VBOR(NPTFR),EPS,G
      DOUBLE PRECISION, INTENT(INOUT) :: AMINF(3,NPTFR) 
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                       
      INTEGER IEL,K                                                                                             
C                                                                       
      DOUBLE PRECISION HI,UI,VI,XN,YN,R1,RLAMB0,HJ,UJ,VJ,R,PI
C                                                                                                                                              
C------                                                                 
C 1. CALCUL DES FLUX AUX BORDS D'ENTREE                                 
C------                                                                 
C                                                                       
      DO 10 K = 1 , NPTFR                                               
C
C     SI H LIBRE OU DEBIT LIBRE                                                                       
      IF(LIMPRO(K,1).EQ.KDDL.OR.LIMPRO(K,2).EQ.KDDL) THEN 
      IEL = NBOR(K)                                                     
      XN = XNEBOR(K)                                                    
      YN = YNEBOR(K)                                                    
      HJ = W(1,IEL)
      IF(HJ.GT.EPS) THEN
        UJ = W(2,IEL)/HJ
        VJ = W(3,IEL)/HJ
        R  =  UJ*XN + VJ*YN-2.D0*SQRT(G*HJ)
        R1 =  UJ*XN + VJ*YN+2.D0*SQRT(G*HJ)
C
C       SI DEBIT IMPOSE
C        
        IF(LIMPRO(K,2).EQ.KDIR) THEN                                      
C
C   Q Donne ; calcul de H pour une entree subsonique
C
          RLAMB0 = UJ * XN + VJ * YN                                    
          IF ( RLAMB0.LE.0.D0) THEN
            PI = -R+RLAMB0
            HI = (PI**2/4.D0)/G
            UI = AMINF(2,K)/HI                                              
            VI = AMINF(3,K)/HI                                              
            AMINF(1,K) = HI                                        
          ENDIF
        ENDIF
C
C       SI H IMPOSE
C
        IF(LIMPRO(K,1).EQ.KDIR) THEN                                                                                                            
C
C   H Donne ; calcul de Q pour une sortie subsonique
C
           HI = AMINF(1,K)
C
           RLAMB0 = (UJ * XN) + (VJ * YN)                                    
           IF ( RLAMB0.GE.-0.0001D0) THEN
             UI = (R1-2.D0*SQRT(G*HI))*XN 
             VI = (R1-2.D0*SQRT(G*HI))*YN 
             AMINF(2,K) = UI*HI
             AMINF(3,K) = VI*HI
           ENDIF
         ENDIF
       ENDIF  
C                                                                       
       ENDIF
C                                                                       
10    CONTINUE                                                          
C                                                                       
1000  FORMAT(' LE SEGMENT DE NO LOCAL ',I2,' EST UN SEGMENT SORTIE !! ')
1010  FORMAT(' LE SEGMENT DE NO LOCAL ',I2,' EST UN SEGMENT ENTREE !! ')
1001  FORMAT(' VALEUR DE < U , N > : ',E12.5)                           
C                                                                       
      RETURN                                                            
      END  
