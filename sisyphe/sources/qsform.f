C                       *****************                                       
                        SUBROUTINE QSFORM
C                       *****************                                                       
C             
C                                                                               
C***********************************************************************        
C SISYPHE VERSION 5.4                                            20/05/96
C                                                
C
C MODIFICATIONS F. HUVELIN NOV 2003
C
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT        
C***********************************************************************        
C                                                                               
C     FONCTION: PERMET A L'UTILISATEUR DE PROGRAMMER LA FORMULE
C               DE TRANSPORT SOLIDE ADAPTEE AU CAS ETUDIE               
C
C              ENABLE USERS TO PROGRAM THE BED-LOAD TRANSPORT FORMULA
C              SUITED FOR THEIR APPLICATION.
C
C    FUNCTION: USER'S SAND TRANSPORT FORMULA
C                                 
C-----------------------------------------------------------------------        
C                             ARGUMENTS                                         
C .________________.____.______________________________________________         
C |      NOM       |MODE|                   ROLE                                
C |________________|____|______________________________________________         
C |   QSS,QSC      |<-- | SUSPENDED AND BED LOAD SAND TRANSPORT                                
C |   HN           | -->| WATER DEPTH
C |   Q            | -->| WATER FLOW RATE                                         
C |   NPOIN        | -->| NUMBER OF POINTS                                      
C |   HW           | -->| WAVE HEIGHT                                     
C |   TW           | -->| WAVE PERIOD
C |   XMVS         | -->| SAND DENSITY
C |   XMVE         | -->| FLOW DENSITY
C |   AC           | -->|  CRITICAL SHIELDS PARAMETER 
C |   DM           | -->|  MEAN SAND DIAMETER                           
C |   GRAV         | -->| GRAVITY                         
C |   VCE          | -->| FLOW VISCOSITY
C |   CHESTR       | -->| FRICTION COEFFICIENT (hezy, Nikuradse or Stickler)
C |   T3,..T10     |<-->| WORKING ARRAYS    
C |   TOB          | -->| MEAN BOTTOM FRICTION
C |   TOBW         | -->| WAVE INDUCED BOTTOM FRICTION
C |   CF           | -->| QUADRATIC FRICTION COEFFICIENT (TOTAL FRICTION)
C |    CFP         | -->| QUADRATIC FRICTION COEFFICIENT (SKIN FRICTION)
C |________________|____|______________________________________________         
C MODE : -->(INPUT), <--(RESULT), <-->(MODIFIED DATA)
C-----------------------------------------------------------------------        
C CALLED BY : BED_LOAD                                                   
C***********************************************************************        
C IF COUPLING, ONLY THE BED LOAD COMPONENT WILL BE CONSERVED 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      USE DECLARATIONS_SISYPHE
      IMPLICIT NONE 
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU                                                            
C
C  FOLLOWING LINES NEED TO BE COMMENTED
C
      IF(LNG.EQ.1) WRITE(LU,52)
      IF(LNG.EQ.2) WRITE(LU,53)
C
52    FORMAT(/,1X,' STOP :',/
     *     ,1X,' LE TAUX DE TRANSPORT DOIT ETRE 
     *       CALCULE DANS QSFORM')
53    FORMAT(/,1X,'SISYPHE IS STOPPED : ',/
     *      ,1X,' SAND TRANSPORT MUST BE CALCULATED IN QSFORM')
      CALL PLANTE(1)
      STOP
C
      RETURN
      END
