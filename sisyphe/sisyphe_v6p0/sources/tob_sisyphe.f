CV : 04/05: correction pour mode Sisyphe seul: ne pas modifier le CHESTR 
C           sauf si KFROT = 0 ou 1
C      
C                       **********************                       
                        SUBROUTINE TOB_SISYPHE
C                       **********************                       
C                                                                     
     * (TOB, TOBW, MU, KS,KSP, KSR,CF,FW,CHESTR,UETCAR,CF_TEL,CODE,
     *  KFROT,ICR, KSPRATIO, HOULE,GRAV,XMVE,XMVS, VCE, KARMAN,
     *  ZERO,HMIN,HN, ACLADM, UNORM,UW, TW, NPOIN)
C 
C***********************************************************************
C  SISYPHE VERSION 6.0  29/11/06       C. VILLARET (LNHE) 01 30 87 83 28
C                                                
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT                                                                          
C***********************************************************************
C                                                                       
C      FONCTION: CALCUL DE LA CONTRAINTE TOTALE AU FOND SELON QUE 
C                    SISYPHE EST COUPLE OU NON                                                   
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                        
C |________________|____|______________________________________________ 
C |    TOB         |<-- |  CONTRAINTE DE FROTTEMENT TOTAL EN COURANT SEUL 
C |    TOBW        |<-- |  CONTRAINTE DE FROTTEMENT  EN HOULE SEULE 
C !    MU          |<-- |  RAPPORT ENTRE LA CONTRAINTE DE FROTTEMENT DE PEAU ET 
C                          LA CONTRAINTE TOTALE
C !    KS          |<-- |  RUGOSITE TOTALE
C !    KSP         |<-- |  RUGOSITE DE PEAU
C !    KSR         |<-- |  RUGOSITE DE RIDE 
C
C !    CF          |<-- |  COEFFICIENT DE FROTTEMENT QUADRATIQUE DU COURANT
C !    FW          |<-- |  COEFFICIENT DE FROTTEMENT QUADRATIQUE DE LA HOULE       
C |    CHESTR      | -->|  COEFFICIENT DE FROTTEMENT (mot clé)                   
C |    UETCAR      | -->|  VITESSE DE FROTTEMENT AU CARRE SI COUPL. T3D
C |    CF_TEL      | -->|  COEFFICIENT DE FROTTMT CF      SI COUPL. T2D
C |    NPOIN       | -->|  NOMBRE DE POINTS DU MAILLAGE 2D 
C !    ICR         | -->|  PREDICTEUR DE RIDE POUR LE FROTTEMENT DE PEAU 
C |    KFROT       | -->|  LOI     DE FROTTEMENT                     
C |    GRAV        | -->|  GRAVITE                                   
C |    ACLADM      ! -->|  DIAMETRE MOYEN  DU SEDIMENT                     
C |    XMVE,XMVS   | -->|  MASSE VOLUMIQUE DE L'EAU, DU SEDIMENT
C !    VCE         | -->|  VISCOSITE DE L'EAU                    
C |    Q           | -->|  DEBIT LIQUIDE                              
C |    HN          | -->|  HAUTEUR D'EAU AU TEMPS N                    
C |    CODE        | -->|  CALLING PROGRAM IN COUPLING
C !    HOULE       | -->|  PRISE EN COMPTE DE LA HOULE
C |    UNORM       | -->|  INTENSITE DU COURANT                      
C |    TW,UW       | -->|  PERIODE DE LA HOULE ET VITESSE ORBITALE  
C |________________|____|______________________________________________ 
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C  APPELE PAR : SISYPHE                                               
C                                                                       
C  SOUS-PROGRAMME APPELE : COEFRO_SISYPHE
C
C********************************************************************** 
C
      USE BIEF
      USE INTERFACE_SISYPHE, EX_TOB_SISYPHE=>TOB_SISYPHE
C                                                                       
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                    
      INTEGER,            INTENT(IN)  :: NPOIN,KFROT,ICR
C      LOGICAL,            INTENT(IN) :: LCONDIS
      LOGICAL,            INTENT(IN)  :: HOULE
      CHARACTER(LEN=24),  INTENT(IN)  :: CODE
      DOUBLE PRECISION,   INTENT(IN)  :: XMVE,XMVS, VCE,GRAV,KARMAN
      DOUBLE PRECISION,   INTENT(IN)  :: ZERO,HMIN,KSPRATIO
      TYPE(BIEF_OBJ), INTENT(IN)      :: UETCAR
      TYPE(BIEF_OBJ), INTENT(IN)      :: HN,UNORM
      TYPE(BIEF_OBJ), INTENT(IN)      :: TW,UW
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: KS,KSP,KSR
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: CHESTR,MU
      TYPE(BIEF_OBJ), INTENT(IN)      :: ACLADM
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: CF,TOB
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: FW,TOBW 
      TYPE(BIEF_OBJ), INTENT(IN)      :: CF_TEL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER                     :: I
      DOUBLE PRECISION            :: A,B,C, HCLIP
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BED ROUGHNESS PREDICTOR
C                         SKIN   : KSP   
C                         TOTAL  : KS 
C                         RIPPLES : KSR   
C                         KS mis dans CHESTR si non-couplage, sinon recalculé 
C  Note: il est conseille d'utiliser loi de frottement 3 en cas de couplage
C        pour eviter un calcul inutil
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C SKIN BED ROUGHNESS --> KSP
C
        CALL OS('X=CY    ', X=KSP, Y=ACLADM, C=KSPRATIO)
C
C RIPPLED BED ROUGHNESS --> KSR =KSP
C
        CALL OS('X=CY    ', X=KSR, Y=ACLADM, C=KSPRATIO)
C
C TOTAL BED ROUGHNESS --> KS
C        KFROT= 0: FLAT BED     KS=KSP
C        KFROT = 1: RIPPLED BED KS= KSP + KSR +KSMR
C
         IF(KFROT.EQ.0) THEN
           CALL OS('X=Y     ', X=KS, Y=KSP)  
          ENDIF
C
          IF(KFROT.EQ.1.OR.ICR.EQ.2) THEN
C           
            IF(HOULE) THEN 
C Wiberg et Harris: KSR (rides)
C                   KS (rides + peau)
              CALL RIDE(KSR%R, TW%R, UW%R, UNORM%R, GRAV, XMVE,
     *                XMVS, VCE, NPOIN, KSPRATIO, ACLADM%R)
              CALL OS('X=Y+Z   ', X=KS, Y=KSP, Z=KSR)              
            ELSE 
C predicteur de VR : KSR (rides)
C                    KS (rides+dunes+megarides) +peau       
              CALL RIDE_VR(KSR%R,KS%R,UNORM%R,HN%R,GRAV,XMVE,
     *                     XMVS,NPOIN,ACLADM%R)
              CALL OS('X=X+Y   ', X=KS, Y=KSP)              
            ENDIF 
C
          ENDIF
C          
C mode Sisyphe seul: on change la valeur du CHESTR seulement si KFROT =1 ou 0
C                                                
          IF(KFROT.EQ.1.OR.KFROT.EQ.0) 
     *      CALL OS('X=Y     ', X=CHESTR, Y=KS) 
C            

C
C ----------------------------------------------------------------------------------------------
C TOTAL HYDRODYNAMIC FRICTION :  --> TOB
C  QUADRATIC COEFT            :  ---> CF
C
C-----------------------------------------------------------------------
C
C     INTERNAL COUPLING WITH TELEMAC2D
C     UETCAR IS CF IN TELEMAC-2D
C
      IF(CODE(1:9).EQ.'TELEMAC2D') THEN
         CALL OV('X=Y     ',CF%R,CF_TEL%R,CF_TEL%R,0.D0,CF%DIM1)     
         DO I=1,NPOIN 
           TOB%R(I) = XMVE*0.5D0*CF%R(I)*UNORM%R(I)**2
         ENDDO
C
C     INTERNAL COUPLING WITH TELEMAC3D
C     UETCAR CORRESPONDS TO THE FRICTION VELOCITY SQUARED
C 
      ELSEIF(CODE(1:9).EQ.'TELEMAC3D') THEN
        CALL OS( 'X=CY     ',X=TOB,Y=UETCAR,C=XMVE)        
        CALL OV('X=Y     ',CF%R,CF_TEL%R,CF_TEL%R,0.D0,CF%DIM1)
C
C   NO  COUPLING : USE KFROT AND CHESTR
C
      ELSE
C
        CALL COEFRO_SISYPHE(CF,HN,KFROT,CHESTR,GRAV,NPOIN,HMIN,KARMAN)
        DO I=1,NPOIN 
          TOB%R(I) = XMVE*0.5D0*CF%R(I)*UNORM%R(I)**2
        ENDDO
C
      ENDIF                  
C 
C ---------------------------------------------------------------------
C TOTAL BED ROUGHNESS CALCULATED AS A FUNCTION OF QUADRATIC BED FRICTION
C UNECESSARY IF (KFROT =0, 1 OR 5) AND (NO COUPLING) 
C              ---->   KS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (CODE(1:8).EQ.'TELEMAC'.OR.KFROT.GE.2.OR.KFROT.LE.4) THEN
         DO I=1,NPOIN
           A = KARMAN*SQRT(2.D0/MAX(CF%R(I),ZERO))
           KS%R(I)=12.D0*HN%R(I)/EXP(A)
         ENDDO
      ENDIF   
C ------------------------------------------------------------------------
C SKIN FRICTION CORRECTOR
C                ---> MU = TOP/TOB 
C ICR=0:    MU=1
C ICR=1     : Skin friction correction USE KSP
C ICR= 2    : Ripple roughness USE KSR, KSR
C En couplage avec telemac: MU>1 est acceptable 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(ICR.EQ.0) THEN
          CALL OS('X=C     ', X=MU, C=1.D0)      
      ELSE IF(ICR.EQ.1) THEN            
         DO I= 1, NPOIN
          IF((CF%R(I) > ZERO).AND.(HN%R(I).GT.KSP%R(I))) THEN
            HCLIP=MAX(HN%R(I),KSP%R(I))  
            A = 2.5D0*LOG(12.D0*HCLIP/ KSP%R(I))
            C =2.D0/A**2  
            MU%R(I) = C/CF%R(I)
          ELSE
             MU%R(I) = 0.D0
          ENDIF
        ENDDO
      ELSE IF(ICR.EQ.2) THEN
        DO I= 1, NPOIN
           IF(HN%R(I).GT.MAX(KSR%R(I),KSP%R(I)).AND.
     *        CF%R(I).GT.ZERO)THEN
                A = LOG(12.D0*HN%R(I)/ KSP%R(I))
                B = LOG(12.D0*HN%R(I)/ KSR%R(I))
                C =0.32D0/CF%R(I)
                MU%R(I) = C/SQRT(B*A**3)
           ELSE
             MU%R(I) = 0.D0
           ENDIF
        ENDDO 
      ENDIF
C  
C -----FROTTEMENT DU A LA HOULE -----------------------------
C 					--> TOBW
C
      IF(HOULE) THEN
        CALL TOBW_SISYPHE
     &          (TOBW%R,CF%R,FW%R,UW%R,TW%R,HN%R,NPOIN,XMVE)
      ENDIF
C 
C------------------------------------------------------------
C
      RETURN                                                            
      END
