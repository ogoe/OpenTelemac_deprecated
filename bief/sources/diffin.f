C                       ******************                               
                        SUBROUTINE DIFFIN                               
C                       ******************                               
C                                                                       
     *(MASKTR,LIMTRA,LITBOR,CLT,U,V,XNEBOR,YNEBOR,NBOR,
     * KP1BOR,NPTFR,KENT,KSORT,KLOG,KINC,KNEU,KDIR,KDDL,
     * ICONV,NELBOR,NPOIN,NELMAX,MSK,MASKEL,
     * NFRLIQ,THOMFR,DEBLIQ,FINLIQ,FRTYPE,TN,TBOR,MESH)
C                                                                       
C***********************************************************************
C  BIEF VERSION 5.9     25/06/2008     J-M HERVOUET (LNH) 01 30 87 80 18
C                                          
C  DEPLACE DE TELEMAC-2D POUR APPEL PAR SISYPHE
C***********************************************************************
C                                                                       
C     FONCTION                                                          
C     ========                                                          
C                                                                       
C     CE SOUS-PROGRAMME INITIALISE LES CONDITIONS AUX LIMITES POUR LA   
C     DIFFUSION DU TRACEUR.                                             
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____._______________________________________________
C |      NOM       |MODE|                   ROLE                        
C |________________|____|_______________________________________________
C |   MASKTR(1,1)  |<-- | MASQUE VALANT UN POUR LES SEGMENTS DIRICHLET  
C |   MASKTR(1,2)  |<-- | MASQUE VALANT UN POUR LES SEGMENTS DDL        
C |   MASKTR(1,3)  |<-- | MASQUE VALANT UN POUR LES SEGMENTS NEUMANN    
C |                |    | (ET ZERO SINON)                               
C |                |    |                                               
C |   LIMTRA       |<-- | CONDITIONS AUX LIMITES POUR LA DISPERSION     
C |   LITBOR       | -->| TYPES DE CONDITIONS AUX LIMITES DU TRACEUR.   
C |   CLT          |<-- | CONDITIONS AUX LIMITES DU TRACEUR POUR        
C |                |    | LE CALCUL (LITBOR MODIFIES).                  
C |   U,V          | -->| COMPOSANTES DU COURANT                        
C |   XNEBOR,YNEBOR| -->| COMPOSANTES DE LA NORMALE EXTERIEURE AU       
C |                |    | DOMAINE PAR POINT DE BORD                     
C |   NBOR         | -->| NUMEROS GLOBAUX DES POINTS DE BORD            
C |   KP1BOR       | -->| UN SEGMENT EST SITUE ENTRE LE POINT DE BORD K 
C |                |    | ET KP1BOR(K).                                 
C |   NPTFR        | -->| DIMENSION DES TABLEAUX .                      
C |   KENT         | -->| INDICATEUR DE POINT D'ENTREE FLUIDE .         
C |   KSORT        | -->| INDICATEUR DE POINT DE SORTIE FLUIDE .        
C |   KLOG         | -->| INDICATEUR DE PAROI SOLIDE .                  
C |   KDIR         | -->| INDICATEUR DE POINT DE DIRICHLET              
C |   KNEU         | -->| INDICATEUR DE POINT DE NEUMANN                
C |   KDDL         | -->| INDICATEUR DE DEGRE DE LIBERTE                
C |   ICONV        | -->| OPTION DE CONVECTION : 1) CARACTERISTIQUES    
C |                |    |                        2) SUPG                
C |                |    |                        3) HYBRIDE             
C |   NELBOR       | -->| NUMEROS DES ELEMENTS ADJACENTS AUX BORDS      
C |   NPOIN        | -->| NOMBRE DE POINTS DU MAILLAGE                  
C |   NELMAX       | -->| NOMBRE MAXIMUM D'ELEMENTS                     
C |   MSK          | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.         
C |   MASKEL       | -->|  TABLEAU DE MASQUAGE DES ELEMENTS             
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE          
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C APPELE PAR : SUBIEF                                                   
C                                                                       
C SOUS-PROGRAMME APPELE : OV                                            
C                                                                       
C********************************************************************** 
C                                                                       
      USE BIEF, EX_DIFFIN => DIFFIN
C 
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C  
      TYPE(BIEF_OBJ), INTENT(INOUT) :: MASKTR,TBOR
      TYPE(BIEF_OBJ), INTENT(IN)    :: TN  
      INTEGER, INTENT(IN)           :: NPOIN,NPTFR,NELMAX,ICONV,NFRLIQ
      INTEGER, INTENT(IN)    :: LITBOR(NPTFR),KP1BOR(NPTFR),NBOR(NPTFR)
      INTEGER, INTENT(INOUT) :: LIMTRA(NPTFR),CLT(NPTFR)   
      INTEGER, INTENT(IN)    :: KENT,KSORT,KLOG,KDIR,KDDL,KNEU,KINC  
      INTEGER, INTENT(IN)    :: NELBOR(NPTFR)       
      INTEGER, INTENT(IN)    :: DEBLIQ(NFRLIQ),FINLIQ(NFRLIQ)
      INTEGER, INTENT(IN)    :: FRTYPE(NFRLIQ)
C
      DOUBLE PRECISION, INTENT(IN) :: U(NPOIN), V(NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: XNEBOR(NPTFR), YNEBOR(NPTFR)  
      DOUBLE PRECISION, INTENT(IN) :: MASKEL(NELMAX)                                                              
C                                                                                                                                              
      LOGICAL, INTENT(IN) :: MSK,THOMFR
C
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K,K1,K2,IELEM
      INTEGER DIR,DDL,NEU,OND,IFRLIQ
      DOUBLE PRECISION USCALN
      LOGICAL DEP  
C
C-----------------------------------------------------------------------
C
      DIR=1
      DDL=2
      NEU=3
      OND=4
C  
C CLT CONTIENT LE TABLEAU LITBOR EVENTUELLEMENT MODIFIE EN FONCTION DU  
C SIGNE DE U.N SUR LES FRONTIERES LIQUIDES, N EST LA NORMALE SORTANTE.
C                                                                     
       DO K=1,NPTFR                                               
         CLT(K) = LITBOR(K)                                             
C        REPERAGE DES FRONTIERES LIQUIDES :        
         IF(CLT(K).EQ.KENT) THEN
          USCALN = U(NBOR(K))*XNEBOR(K) + V(NBOR(K))*YNEBOR(K)
C         VITESSE SORTANTE, TRACEUR LIBRE
          IF(USCALN.GT.0.D0) CLT(K) = KSORT
         ELSEIF(CLT(K).EQ.KSORT) THEN
          USCALN = U(NBOR(K))*XNEBOR(K) + V(NBOR(K))*YNEBOR(K)
C         VITESSE ENTRANTE, TRACEUR IMPOSE A LA VALEUR
C         PRECEDENTE
          IF(USCALN.LT.0.D0) THEN
            TBOR%R(K)=TN%R(NBOR(K))
            CLT(K) = KENT
          ENDIF 
         ENDIF 
       ENDDO                             
C                       
C                                                                       
C  CONSTRUCTION DU TABLEAU MASKTR EN FONCTION DE CLT                    
C                                                                       
C  MASKTR VAUT 1 POUR UN SEGMENT DE TYPE NEUMANN ET ZERO SINON          
C                                                                       
C  UN SEGMENT EST DE TYPE NEUMANN SI UN DE SES POINTS AU MOINS EST      
C  DONNE COMME NEUMANN PAR L'UTILISATEUR.                               
C 
C                                                                       
C     INITIALISATION A ZERO DES MASQUES
C
      CALL OS('X=0     ',MASKTR)
      DO K1 = 1 , NPTFR                                              
        K2 = KP1BOR(K1)
C       K2=K1 => EN PARALLELISME, SUIVANT DANS UN AUTRE SOUS-DOMAINE
C                                 DANS CE CAS LE MASQUE PAR SEGMENT
C                                 NE DOIT PAS SERVIR
        IF(K2.NE.K1) THEN                                                 
          IF(CLT(K1).EQ.KLOG.OR.CLT(K2).EQ.KLOG) THEN               
C           SEGMENTS DE TYPE NEUMANN                                      
            MASKTR%ADR(NEU)%P%R(K1)=1.D0                                      
          ELSEIF(CLT(K1).EQ.KENT.AND.CLT(K2).EQ.KSORT) THEN
C           SEGMENTS DE TYPE SORTIE
            MASKTR%ADR(DDL)%P%R(K1)=1.D0
          ELSEIF(CLT(K1).EQ.KSORT.OR.CLT(K2).EQ.KSORT) THEN         
            MASKTR%ADR(DDL)%P%R(K1)=1.D0                                      
          ELSEIF(CLT(K1).EQ.KSORT.AND.CLT(K2).EQ.KENT) THEN
C           SEGMENTS DE TYPE SORTIE
            MASKTR%ADR(DDL)%P%R(K1)=1.D0
          ELSEIF(CLT(K1).EQ.KENT.OR.CLT(K2).EQ.KENT) THEN           
            MASKTR%ADR(DIR)%P%R(K1)=1.D0                                      
          ELSEIF(CLT(K1).EQ.KINC.OR.CLT(K2).EQ.KINC) THEN           
            MASKTR%ADR(OND)%P%R(K1)=1.D0                                      
          ELSE                                                            
            IF(LNG.EQ.1) WRITE(LU,100)                                   
            IF(LNG.EQ.2) WRITE(LU,101)                                   
100         FORMAT(1X,'DIFFIN : CAS NON PREVU')                           
101         FORMAT(1X,'DIFFIN : UNEXPECTED CASE')                         
            CALL PLANTE(1)                                                
            STOP                                                          
          ENDIF 
        ENDIF                                                          
      ENDDO
C
C     EN PARALLELISME, ON RECUPERE LES 1 DONNES PAR UN AUTRE
C     SOUS-DOMAINE (DONC ON PREND LE MAX DE CHAQUE POINT)
C
      IF(NCSIZE.GT.1) THEN
        CALL PARCOM_BORD(MASKTR%ADR(NEU)%P%R,3,MESH)
        CALL PARCOM_BORD(MASKTR%ADR(DDL)%P%R,3,MESH)
        CALL PARCOM_BORD(MASKTR%ADR(DIR)%P%R,3,MESH)
        CALL PARCOM_BORD(MASKTR%ADR(OND)%P%R,3,MESH)
      ENDIF                                               
C                                                                       
C  MASQUAGE EVENTUEL                                                    
C                                                                       
      IF(MSK) THEN                                                      
        DO K1 = 1 , NPTFR
          IELEM=NELBOR(K1)
          IF(IELEM.GT.0) THEN
            MASKTR%ADR(DIR)%P%R(K1) = MASKTR%ADR(DIR)%P%R(K1) *
     *                                                     MASKEL(IELEM)
            MASKTR%ADR(DDL)%P%R(K1) = MASKTR%ADR(DDL)%P%R(K1) * 
     *                                                     MASKEL(IELEM)
            MASKTR%ADR(NEU)%P%R(K1) = MASKTR%ADR(NEU)%P%R(K1) * 
     *                                                     MASKEL(IELEM)
            MASKTR%ADR(OND)%P%R(K1) = MASKTR%ADR(OND)%P%R(K1) * 
     *                                                     MASKEL(IELEM)
          ENDIF
        ENDDO
C
C       EN PARALLELISME, ON RECUPERE LES 0 DONNES PAR UN AUTRE
C       SOUS-DOMAINE (DONC ON PREND LE MIN DE CHAQUE POINT)
C
        IF(NCSIZE.GT.1) THEN
          CALL PARCOM_BORD(MASKTR%ADR(NEU)%P%R,4,MESH)
          CALL PARCOM_BORD(MASKTR%ADR(DDL)%P%R,4,MESH)
          CALL PARCOM_BORD(MASKTR%ADR(DIR)%P%R,4,MESH)
          CALL PARCOM_BORD(MASKTR%ADR(OND)%P%R,4,MESH)
        ENDIF                                                            
      ENDIF                                                         
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C PASSAGE DES CONDITION PHYSIQUES AUX CONDITIONS TECHNIQUES             
C                                                                       
      DO 4 K=1,NPTFR                                                    
C                                                                       
        IF(CLT(K).EQ.KENT ) THEN                                     
C                                                                       
C  ENTREE DE DOMAINE : TRACEUR IMPOSE                                   
C                                                                       
          LIMTRA(K) = KDIR                                              
C                                                                       
        ELSEIF(CLT(K).EQ.KSORT) THEN                                 
C                                                                       
C  SORTIE DE DOMAINE : LIBRE SI SUPG OU SCHEMA PSI,                     
C                      RESULTAT DE LA CONVECTION IMPOSE SINON           
C                                                                       
          IF(ICONV.EQ.1) THEN
C           VOIR DIFFCL OU ON MET TTILD DANS TBOR
            LIMTRA(K) = KDIR                                            
          ELSE                                                          
            LIMTRA(K) = KDDL                                            
          ENDIF                                                         
C                                                                       
        ELSEIF(CLT(K).EQ.KLOG ) THEN                                 
C                                                                       
C  PAROI : CONDITIONS DE NEUMANN (EN FAIT ON NE S'EN SERT PAS)          
C                                                                       
          LIMTRA(K) = KNEU                                              
C                                                                       
        ELSE                                                            
C                                                                       
C  ERREUR, VALEUR DE LITBOR INCONNUE                                    
C                                                                       
          IF(LNG.EQ.1) WRITE(LU,10) K,LITBOR(K)                         
          IF(LNG.EQ.2) WRITE(LU,12) K,LITBOR(K)                         
10        FORMAT(1X,'DIFFIN: POINT ',1I6,' LITBOR= ',1I6,' ?????')      
12        FORMAT(1X,'DIFFIN: POINT ',1I6,' LITBOR= ',1I6,' ?????') 
          CALL PLANTE(1)     
          STOP                                                          
C                                                                       
        ENDIF                                                           
C                                                                       
4     CONTINUE                                                          
C                                                                       
C---------------------------------------------------------------------- 
C POST-TRAITEMENT POUR LES CONDITIONS AUX LIMITES LIQUIDES TRAITEES
C AVEC LA METHODE DE THOMPSON. LA CONDITION A LA LIMITE POUR LE
C TRACEUR EST ALORS DE TYPE DIRICHLET.
C
C NE MARCHERA PAS EN PARALLELE (MAIS THOMSON A BESOIN DES 
C CARACTERISTIQUES)
C
      IF(NFRLIQ.NE.0.AND.THOMFR) THEN
C
      DO 17 IFRLIQ = 1 , NFRLIQ
        DEP = .FALSE.
        K = DEBLIQ(IFRLIQ)
15      CONTINUE
        IF(FRTYPE(IFRLIQ).EQ.2) LIMTRA(K) = KDIR
        IF(K.EQ.FINLIQ(IFRLIQ).AND.DEP) THEN
          GO TO 16
        ELSE
          DEP=.TRUE.
          K = KP1BOR(K)
          GO TO 15
        ENDIF
16      CONTINUE
17    CONTINUE
C
      ENDIF
C
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END
