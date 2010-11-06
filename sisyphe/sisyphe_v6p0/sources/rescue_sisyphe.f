C                       *************************                       
                        SUBROUTINE RESCUE_SISYPHE
C                       *************************                       
C                                                                       
     *(QU,QV,Q,U,V,H,S,ZF,HW,TW,THETAW,NPOIN,TROUVE,ALIRE,PASS,
     * ICF,LISTI)                                                   
C                                                                       
C***********************************************************************
C SISYPHE VERSION 6.0                             C. LENORMANT
C                                                
C
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT 
C***********************************************************************
C                                                                       
C   FONCTION  : RECONSTITUTION DE DONNEES MANQUANTES POUR 
C               L'HYDRODYNAMIQUE ET/OU LES SUITES DE CALCUL 
C               SEDIMENTOLOGIQUES
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   H            |<-- | HAUTEUR.                                     |
C |   S            |<-- | SURFACE LIBRE.                               |
C |   ZF           |<-- | COTE DES POINTS DU FOND.                     |
C |   KX,KY        |<-- | COEFFICIENTS DE DISPERSION                   |
C |   NPOIN        | -->| NOMBRE DE POINTS DU MAILLAGE                 |
C |   NTRAC        | -->| NOMBRE DE TRACEURS DU CAS A TRAITER          |
C |   NPRIV        | -->| NOMBRE DE TABLEAUX UTILISATEURS DU CAS       |
C |   TROUVE       | -->| LOGIQUE INDIQUANT LES VARIABLES TROUVEES     |
C |                |    | DANS LE SOUS-PROGRAMME SUITE                 |
C |   ALIRE        |    | TABLEAU DES VARIABLES A LIRE                 |
C |   PASS         | -->| LOGIQUE VRAI SI ON EST EN DEBUT DE CALCUL    |
C |   LISTI        | -->| LOGIQUE VRAI SI ON IMPRIME DES MESSAGES      |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C APPELE PAR : SUBIEF                                                   
C                                                                       
C SOUS-PROGRAMME APPELE : OV                                            
C                                                                       
C***********************************************************************
C
      USE BIEF
C
      USE INTERFACE_SISYPHE, EX_RESCUE_SISYPHE 
     *           => RESCUE_SISYPHE
C
      USE DECLARATIONS_SISYPHE, ONLY : MAXVAR
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: TROUVE(MAXVAR),ALIRE(MAXVAR),NPOIN,ICF
      LOGICAL, INTENT(IN) :: PASS,LISTI
C                                                                       
      DOUBLE PRECISION, INTENT(INOUT) :: QU(NPOIN), QV(NPOIN), Q(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: U(NPOIN) , V(NPOIN)   
      DOUBLE PRECISION, INTENT(INOUT) :: S(NPOIN) , ZF(NPOIN), H(NPOIN) 
      DOUBLE PRECISION, INTENT(INOUT) :: HW(NPOIN), TW(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: THETAW(NPOIN) 
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                      
      INTEGER K                                          
C                                                                       
C-----------------------------------------------------------------------
C
C IMPRESSIONS :
C -----------
      IF(PASS.AND.LISTI) THEN 
        WRITE(LU,200)
200     FORMAT(80('-'))            
        IF(ALIRE(8).EQ.1) THEN
          IF(LNG.EQ.1) WRITE(LU,300)                                    
          IF(LNG.EQ.2) WRITE(LU,301)                                    
300       FORMAT(1X,'RESCUE : FICHIER HYDRODYNAMIQUE')                  
301       FORMAT(1X,'RESCUE : HYDRODYNAMIC FILE')
        ELSE 
          IF(LNG.EQ.1) WRITE(LU,310)                                    
          IF(LNG.EQ.2) WRITE(LU,311)                                    
310       FORMAT(1X,'RESCUE : FICHIER SEDIMENTOLOGIQUE')                
311       FORMAT(1X,'RESCUE : SEDIMENTOLOGICAL FILE')
        ENDIF 
      ENDIF                       
C
C ------------------------------------------------------------------    
C  HAUTEUR D'EAU :                                                      
C  -------------                                                        
      IF((ALIRE(3).EQ.1).AND.(TROUVE(3).NE.1)) THEN                     
        IF(TROUVE(4).EQ.1.AND.TROUVE(5).EQ.1) THEN 
          IF (LISTI) THEN                          
            IF(LNG.EQ.1) WRITE(LU,400)                                  
            IF(LNG.EQ.2) WRITE(LU,401)                                  
          ENDIF                  
          CALL OV( 'X=Y-Z   ' , H , S , ZF , 0.D0 , NPOIN )             
        ELSE 
          IF (LISTI) THEN                                               
            IF(LNG.EQ.1) WRITE(LU,420)                                  
            IF(LNG.EQ.2) WRITE(LU,421)
          ENDIF 
          CALL PLANTE(1)                                      
          STOP                                                          
        ENDIF                                                           
      ENDIF 
C
400       FORMAT(1X,'HAUTEUR D''EAU CALCULEE AVEC LE FOND',    
     *         /,1X,'ET LA SURFACE LIBRE')                     
401       FORMAT(1X,'WATER DEPTH COMPUTED WITH BATHYMETRY',    
     *         /,1X,' AND SURFACE ELEVATION')
420       FORMAT(1X,'IMPOSSIBLE DE CALCULER LA HAUTEUR D''EAU')
421       FORMAT(1X,'WATER DEPTH UNABLE TO BE COMPUTED')                
C  
C ----------------------------------------------------------------------
C
C CLIPPING DES HAUTEURS D'EAU NEGATIVES :
C -------------------------------------
C
      DO K = 1,NPOIN
        H(K) = MAX(H(K),0.D0)
      ENDDO 
C     
C------------------------------------------------------------------------
C
C  HAUTEUR DE HOULE ET PERIODE DE HOULE                                                         
C 
      IF(ICF==4.OR.ICF==5.OR.ICF==8.OR.ICF==9) THEN
C
        IF(ALIRE(12).EQ.1.AND.TROUVE(12).EQ.0) THEN         
          IF(LNG.EQ.1) WRITE(LU,900)                       
          IF(LNG.EQ.2) WRITE(LU,901) 
          CALL OV( 'X=C     ' , HW , U , V , 0.D0 , NPOIN )
        ENDIF
C         
900     FORMAT(1X,'CALCUL PRECEDENT SANS LA HAUTEUR DE HOULE : ON',       
     *          ' PREND ZERO')                   
901     FORMAT(1X,'PREVIOUS COMPUTATION WITHOUT WAVE HEIGHT : IT IS',     
     *          ' FIXED TO ZERO')     
C
         IF(ALIRE(13).EQ.1.AND.TROUVE(13).EQ.0) THEN         
           IF(LNG.EQ.1) WRITE(LU,902)                        
           IF(LNG.EQ.2) WRITE(LU,903)
           CALL OV( 'X=C     ' , TW , U , V , 0.D0 , NPOIN )
         ENDIF
902     FORMAT(1X,'CALCUL PRECEDENT SANS LA PERIODE DE HOULE : ON',       
     *          ' PREND ZERO')                       
903     FORMAT(1X,'PREVIOUS COMPUTATION WITHOUT WAVE PERIOD : IT IS',     
     *          ' FIXED TO ZERO')     
C                                        
         IF(ALIRE(14).EQ.1.AND.TROUVE(14).EQ.0) THEN         
           IF(LNG.EQ.1) WRITE(LU,902)                        
           IF(LNG.EQ.2) WRITE(LU,903)
           CALL OV( 'X=C     ' , THETAW , U , V , 90.D0  , NPOIN )
         ENDIF
      ENDIF
909   FORMAT(1X,'CALCUL PRECEDENT SANS ANGLE DE HOULE : ON',       
     *          ' PREND ZERO')                       
910   FORMAT(1X,'PREVIOUS COMPUTATION WITHOUT WAVE ANGLE : IT IS',     
     *          ' FIXED TO ZERO')                                 
C 
C-----------------------------------------------------------------------
C  FOND NON ERODABLE                                                       
C                                                                     
      IF(ALIRE(9).EQ.1.AND.TROUVE(9).EQ.0) THEN         
        IF(LNG.EQ.1) WRITE(LU,907)                       
        IF(LNG.EQ.2) WRITE(LU,908)                     
      ENDIF
907   FORMAT(1X,'CALCUL PRECEDENT SANS FOND NON ERODABLE')                  
908   FORMAT(1X,'PREVIOUS CALCULATION WITHOUT NON ERODABLE', 
     *         /,1X,'BOTTOM')                                                                   
C
C-----------------------------------------------------------------------                                                                       
C  COTE DU FOND                                                         
C                                                                       
      IF(ALIRE(5).EQ.1.AND.TROUVE(5).EQ.0) THEN  
C        
        IF(TROUVE(4).EQ.1.AND.TROUVE(3).EQ.1) THEN 
          IF (LISTI) THEN                      
          IF(LNG.EQ.1) WRITE(LU,410)                                    
          IF(LNG.EQ.2) WRITE(LU,411)                                    
410       FORMAT(1X,'FOND CALCULE AVEC LA HAUTEUR D''EAU',    
     *         /,1X,'ET LA SURFACE LIBRE')                     
411       FORMAT(1X,'BATHYMETRY COMPUTED FROM WATER DEPTH',    
     *         /,1X,'AND SURFACE ELEVATION')    
          ENDIF               
          CALL OV( 'X=Y-Z   ' , ZF , S , H , 0.D0 , NPOIN )             
        ELSE                                                                 
          CALL  OV( 'X=C     ' , ZF , ZF, ZF, 0.D0 , NPOIN )             
          IF(LNG.EQ.1) WRITE(LU,960)                                   
          IF(LNG.EQ.2) WRITE(LU,961)  
        ENDIF                                    
960     FORMAT(1X,'COTE DU FOND NON TROUVEE',/,
     *            'LA COTE EST INITIALISEE A ZERO')         
961     FORMAT(1X,'BOTTOM TOPOGRAPHY NOT FOUND',/,
     *            'IT IS SET TO ZERO')            
C       
      ENDIF                                                            
C
      IF (PASS.AND.LISTI) THEN 
        WRITE(LU,970)
970     FORMAT(80('-')) 
      ENDIF
C                                                             
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END SUBROUTINE RESCUE_SISYPHE
