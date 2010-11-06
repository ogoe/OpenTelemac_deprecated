C                    **********************************                       
                     SUBROUTINE RESCUE_SISYPHE_NOTPERMA
C                    **********************************
C                                                                                                
     *(QU,QV,Q,U,V,H,S,ZF,HW,TW,THETAW,NPOIN,TROUVE,ALIRE,ICF,ENTET)                                                   
C                                                                       
C***********************************************************************
C SISYPHE VERSION 5.7                           C.VILLARET
C                                                     
C                                                
C RESCUE_sisyphe doit être modifié en non-permmanent pour tenir compte 
C des évolutions du fond
C --> la hauteur d'eau doit être recalulcée à partir du fond (SISYPHE) 
C    et de la surface libre (FICHIER HYDRO
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
      USE INTERFACE_SISYPHE, EX_RESCUE_SISYPHE_NOTPERMA 
     *           => RESCUE_SISYPHE_NOTPERMA
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
C                                                                       
      DOUBLE PRECISION, INTENT(INOUT) :: QU(NPOIN), QV(NPOIN), Q(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: U(NPOIN) , V(NPOIN)   
      DOUBLE PRECISION, INTENT(INOUT) :: S(NPOIN) , ZF(NPOIN), H(NPOIN) 
      DOUBLE PRECISION, INTENT(INOUT) :: HW(NPOIN), TW(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: THETAW(NPOIN)
      LOGICAL, INTENT(IN)             :: ENTET 
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     LES VARIABLES ESSENTIELLES NON MODIFIES PAR LES EVOLUTIONS 
C     DU FOND SONT LA SURFACVE LIBRE ET LE DEBIT
C
C     reconstitutuion surface libre(4) à partir de la hauteur d'eau 
C     non modifiée et du fond (non modifié)
C                                                
      IF(TROUVE(4).EQ.0.AND.ALIRE(4).EQ.1) THEN 
        IF(TROUVE(3).EQ.1.AND.TROUVE(5).EQ.1) THEN 
          CALL OV( 'X=Y+Z   ' , S , H , ZF , 0.D0 , NPOIN )               
        ELSE 
          WRITE(LU,*) 'UNABLE TO BUILD FREE SURFACE'
          CALL PLANTE(1)
          STOP                                                         
        ENDIF                                                           
      ENDIF 
C 
C reconstitution hauteur d'eau (3) du calcul Telemac
C à partir du fond (non modifié) et de la surface libre 
C
      IF(TROUVE(3).EQ.0.AND.ALIRE(3).EQ.1) THEN 
        IF(TROUVE(5).EQ.1) THEN
          CALL OV( 'X=Y-Z   ' , H , S , ZF , 0.D0 , NPOIN ) 
        ELSE 
          WRITE(LU,*) 'MISSING BOTTOM OR WATER DEPTH'
          CALL PLANTE(1)
          STOP             
        ENDIF 
      ENDIF 
C                                                                       
C  DEBIT VECTORIEL SUIVANT X :  QU
C
      IF (ALIRE(7).EQ.1.AND.TROUVE(7).EQ.0)  THEN
        IF (TROUVE(1).EQ.1) THEN
          IF(ENTET) THEN
          IF (LNG.EQ.1) WRITE(LU,150)
          IF (LNG.EQ.2) WRITE(LU,151)
150       FORMAT(1X,'DEBIT /X CALCULE AVEC LA VITESSE ET LA HAUTEUR')
151       FORMAT(1X,'DISCHARGE /X COMPUTED WITH DEPTH AND VELOCITY')
          ENDIF
          CALL OV( 'X=YZ    ' , QU , U , H  , 0.D0 , NPOIN )
        ELSE
          IF(ENTET) THEN
          IF (LNG.EQ.1) WRITE(LU,190)
          IF (LNG.EQ.2) WRITE(LU,191)
190       FORMAT(1X,'CALCUL PRECEDENT SANS LA VITESSE U: ON PREND ZERO')
191       FORMAT(1X,'PREVIOUS COMPUTATION WITHOUT VELOCITY U: U IS
     *               EQUAL TO ZERO')
          ENDIF
          CALL OV( 'X=C     ' , QU , U , H, 0.D0 , NPOIN )
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
C  DEBIT VECTORIEL SUIVANT Y :  QV
C
      IF (ALIRE(8).EQ.1.AND.TROUVE(8).EQ.0)  THEN
        IF (TROUVE(2).EQ.1) THEN
          IF(ENTET) THEN
          IF (LNG.EQ.1) WRITE(LU,160)
          IF (LNG.EQ.2) WRITE(LU,161)
160       FORMAT(1X,'DEBIT /Y CALCULE AVEC LA VITESSE ET LA HAUTEUR')
161       FORMAT(1X,'DISCHARGE /Y COMPUTED WITH DEPTH AND VELOCITY')
          ENDIF
          CALL OV( 'X=YZ    ' , QV , V , H , 0.D0 , NPOIN )
        ELSE
          IF(ENTET) THEN
          IF (LNG.EQ.1) WRITE(LU,210)
          IF (LNG.EQ.2) WRITE(LU,211)
210       FORMAT(1X,'CALCUL PRECEDENT SANS LA VITESSE V: ON PREND ZERO')
211       FORMAT(1X,'PREVIOUS COMPUTATION WITHOUT VELOCITY V:
     *               FIXED TO ZERO')
          ENDIF
          CALL OV( 'X=C     ' , QV , V , H , 0.D0, NPOIN )
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
C  DEBIT (m2/s)
C
      IF ((ALIRE(6).EQ.1).AND.(TROUVE(6).EQ.0))  THEN
       CALL OV( 'X=N(Y,Z)' , Q  , QU  , QV  , 0.D0 , NPOIN )
      ENDIF
C     
C-----------------------------------------------------------------------
C  HAUTEUR DE HOULE ET PERIODE DE HOULE                                                         
C 
      IF ((ICF==4).OR.(ICF==5).OR.
     &    (ICF==8).OR.(ICF==9)    ) THEN
C
         IF ((ALIRE(12).EQ.1).AND.(TROUVE(12).EQ.0)) THEN         
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
         IF ((ALIRE(13).EQ.1).AND.(TROUVE(13).EQ.0)) THEN         
            IF(LNG.EQ.1) WRITE(LU,902)                        
            IF(LNG.EQ.2) WRITE(LU,903)
            CALL OV( 'X=C     ' , TW , U , V , 0.D0 , NPOIN )
         ENDIF
902     FORMAT(1X,'CALCUL PRECEDENT SANS LA PERIODE DE HOULE : ON',       
     *          ' PREND ZERO')                       
903     FORMAT(1X,'PREVIOUS COMPUTATION WITHOUT WAVE PERIOD : IT IS',     
     *          ' FIXED TO ZERO')     
C                                        
         IF ((ALIRE(14).EQ.1).AND.(TROUVE(14).EQ.0)) THEN         
            IF(LNG.EQ.1) WRITE(LU,902)                        
            IF(LNG.EQ.2) WRITE(LU,903)
            CALL OV( 'X=C     ' , THETAW , U , V , 90.D0  , NPOIN )
         ENDIF
      ENDIF
909     FORMAT(1X,'CALCUL PRECEDENT SANS ANGLE DE HOULE : ON',       
     *          ' PREND ZERO')                       
910     FORMAT(1X,'PREVIOUS COMPUTATION WITHOUT WAVE ANGLE : IT IS',     
     *          ' FIXED TO ZERO')                                 
C 
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END SUBROUTINE RESCUE_SISYPHE_NOTPERMA
